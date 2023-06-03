#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <omp.h>
#include <random>
#include <chrono>
#include <fstream>
#include <string>
#include <sstream>
#define G_CONST 6.67428e-11
// Define particle structure
struct Particle {
    float posi[2];
    float velocity[2];
    float acceleration[2];
    float acceleration_prev[2];
    float mass;
};

// 自定义reduction运算符：向量相加
#pragma omp declare reduction(vec_float_plus : std::vector<float> : \
    std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<float>())) \
    initializer(omp_priv = decltype(omp_orig)(omp_orig.size(), 0.0))

// 计算系统动量
std::vector<float> calculate_system_momentum(const std::vector<Particle>& particles) {
    std::vector<float> system_momentum(2, 0.0);

    #pragma omp parallel for reduction(vec_float_plus:system_momentum)
    for (int i = 0; i < particles.size(); ++i) {
        const Particle& particle = particles[i];
        system_momentum[0] += particle.mass * particle.velocity[0];
        system_momentum[1] += particle.mass * particle.velocity[1];
    }

    return system_momentum;
}

// 计算系统角动量
std::vector<float> calculate_system_angular_momentum(const std::vector<Particle>& particles) {
    std::vector<float> system_angular_momentum(1, 0.0);

    #pragma omp parallel for reduction(vec_float_plus:system_angular_momentum)
    for (int i = 0; i < particles.size(); ++i) {
        const Particle& particle = particles[i];
        system_angular_momentum[0] += particle.mass * (particle.posi[0] * particle.velocity[1] - particle.posi[1] * particle.velocity[0]);
    }

    return system_angular_momentum;
}


float calculate_system_energy(const std::vector<Particle>& particles, float G) {
    float total_kinetic_energy = 0.0;
    float total_potential_energy = 0.0;

    #pragma omp parallel for reduction(+:total_kinetic_energy, total_potential_energy)
    for (size_t i = 0; i < particles.size(); i++) {
        const auto& p = particles[i];

        // Calculate kinetic energy
        float speed_squared = p.velocity[0] * p.velocity[0] + p.velocity[1] * p.velocity[1];
        float kinetic_energy = 0.5 * p.mass * speed_squared;
        #pragma omp atomic
        total_kinetic_energy += kinetic_energy;

        // Calculate potential energy
        for (size_t j = 0; j < particles.size(); j++) {
            if (i == j) {
                continue;
            }

            const auto& other_p = particles[j];
            float dx = other_p.posi[0] - p.posi[0];
            float dy = other_p.posi[1] - p.posi[1];
            float distance = std::sqrt(dx * dx + dy * dy);
            float potential_energy = -G * p.mass * other_p.mass / distance;
            #pragma omp atomic
            total_potential_energy += potential_energy;
        }
    }

    return total_kinetic_energy + total_potential_energy;
}

__global__ void simulate_particles_cuda_kernel(Particle* particles, const float* particle_masses, float G, float dt, int num_particles) {
    extern __shared__ Particle sharedParticles[];

    int index = blockIdx.x * blockDim.x + threadIdx.x;

    if (index < num_particles) {
        Particle& p = particles[index];
        Particle& shared_p = sharedParticles[threadIdx.x];

        shared_p = p;

        __syncthreads();

        for (int j = 0; j < num_particles; j++) {
            if (j != index) {
                Particle& other = sharedParticles[j];
                float dx = other.posi[0] - shared_p.posi[0];
                float dy = other.posi[1] - shared_p.posi[1];
                float dist_squared = dx * dx + dy * dy;
                float dist_cubed = dist_squared * sqrt(dist_squared);

                float force_magnitude = G * shared_p.mass * particle_masses[j] / dist_cubed;
                float force_x = force_magnitude * dx;
                float force_y = force_magnitude * dy;

                atomicAdd(&shared_p.acceleration[0], force_x / shared_p.mass);
                atomicAdd(&shared_p.acceleration[1], force_y / shared_p.mass);
            }
        }

        p.posi[0] += p.velocity[0] * dt + 0.5 * p.acceleration[0] * dt * dt;
        p.posi[1] += p.velocity[1] * dt + 0.5 * p.acceleration[1] * dt * dt;
        p.velocity[0] += 0.5 * (p.acceleration[0] + p.acceleration_prev[0]) * dt;
        p.velocity[1] += 0.5 * (p.acceleration[1] + p.acceleration_prev[1]) * dt;
        p.acceleration_prev[0] = p.acceleration[0];
        p.acceleration_prev[1] = p.acceleration[1];

        p.acceleration[0] = 0.0f;
        p.acceleration[1] = 0.0f;
    }
}


int main() {
    // Define simulation parameters
    // const float G = 6.674e-11;

    // std::cout << std::fixed << std::setprecision(12);

    // // Input time and time step
    // float t = 0.001;
    // float dt = 0.001;

    // // Perform simulation
    // int num_steps = t / dt;
    // int particles_num = 100000;

    // Random number generator
    // std::random_device rd;
    // std::mt19937 gen(rd());
    // std::uniform_real_distribution<float> dist_pos(-1.0f, 1.0f);
    // std::uniform_real_distribution<float> dist_vel(-10.0f, 10.0f);
    // std::uniform_real_distribution<float> dist_mass(1.0e6f, 1.0e7f);

    // // Generate particles
    // std::vector<Particle> particles;
    // particles.reserve(particles_num);
    // for (int i = 0; i < particles_num; ++i) {
    //     Particle particle;
    //     particle.posi[0] = dist_pos(gen);
    //     particle.posi[1] = dist_pos(gen);
    //     particle.velocity[0] = dist_vel(gen);
    //     particle.velocity[1] = dist_vel(gen);
    //     particle.acceleration[0] = 0.0f;
    //     particle.acceleration[1] = 0.0f;
    //     particle.mass = dist_mass(gen);
    //     particles.push_back(particle);
    // }




    std::ifstream file("one_step_data.csv");
    std::vector<Particle> particles; // store the particles
    std::string line; // each line of the file data
    bool isFirstLine = true; // if line is the first line or not
    int particleNum =0;
    // to get all of the data
    while (std::getline(file, line)) {
        // delete the first line because the first line is
        // "Time,Particle,Mass,PositionX,PositionY,VectorX,VectorY,AccelerationX,AccelerationY"
        if (isFirstLine) {
            isFirstLine = false;
            continue; 
        }

        // set the parameter to store data
        std::istringstream iss(line);
        std::string element;
        std::vector<std::string> elements; 
        // change the line string into a vector
        while (std::getline(iss, element, ',')) {
            elements.push_back(element);
        }

        Particle particle;
        // input data
        double m = std::stod(elements[2]);
        double x = std::stod(elements[3]);
        double y = std::stod(elements[4]);
        double vx = std::stod(elements[5]);
        double vy = std::stod(elements[6]);
        double ax = std::stod(elements[7]);
        double ay = std::stod(elements[8]);
        // add the particle initial condition
        particle.mass = m;
        particle.posi[0]=(x);
        particle.posi[1]=(y);
        particle.velocity[0]=(vx);
        particle.velocity[1]=(vy);
        particle.acceleration[0]=(ax);
        particle.acceleration[1]=(ay);

        // add the particle
        particleNum++;
        particles.push_back(particle);
    }

    // int index = 1;
    // for (const auto& particle : particles) {        
    //     std::cout << "Particle: " << index << ", Mass: " << particle.mass 
    //     << ", PositionX: " << particle.posi[0] << ", PositionY: " << particle.posi[1] 
    //     << ", VectorX: " << particle.velocity[0] << ", VectorY: " << particle.velocity[1] 
    //     << ", AccelerationX: " << particle.acceleration[0] << ", AccelerationY: " << particle.acceleration[1] << std::endl;
    //     index++;
    // }
    double t=0.005 ,dt = 0.005;
    int num_steps = t / dt;
    std::cout << "Physical Time: " << (num_steps)*dt <<"seconds"<< std::endl;
    std::cout << particleNum <<  "particles"<< std::endl;
    std::cout << num_steps <<  "steps"<< std::endl;

    // Start timer
    auto start = std::chrono::high_resolution_clock::now();

    // Transfer particles to device memory
    Particle* d_particles;
    size_t particlesSize = particles.size() * sizeof(Particle);
    cudaMalloc((void**)&d_particles, particlesSize);
    cudaMemcpy(d_particles, particles.data(), particlesSize, cudaMemcpyHostToDevice);

    // Transfer particle masses to constant memory
    float* d_masses;
    size_t particleMassesSize = particles.size() * sizeof(float);
    cudaMalloc((void**)&d_masses, particleMassesSize);
    std::vector<float> particle_masses(particles.size());
    for (int i = 0; i < particles.size(); ++i) {
        particle_masses[i] = particles[i].mass;
    }
    cudaMemcpy(d_masses, particle_masses.data(), particleMassesSize, cudaMemcpyHostToDevice);

    // Launch kernel
    int blockSize = 1024;
    int gridSize = (particles.size() + blockSize - 1) / blockSize;

    int sharedMemSize = blockSize * sizeof(Particle);
    simulate_particles_cuda_kernel<<<gridSize, blockSize, sharedMemSize>>>(d_particles, d_masses, G_CONST, dt, particles.size());

    // Transfer particles back to host memory
    cudaMemcpy(particles.data(), d_particles, particlesSize, cudaMemcpyDeviceToHost);

    // Cleanup
    cudaFree(d_particles);
    cudaFree(d_masses);

    // End timer
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate duration
    std::chrono::duration<float> duration = end - start;
    float seconds = duration.count();

    std::vector<float> system_momentum = calculate_system_momentum(particles);
    float system_energy = calculate_system_energy(particles, G_CONST);

    std::cout << "Verlet_velocity Time: " << (num_steps) * dt << std::endl;
    std::cout << "System Momentum (X, Y): (" << system_momentum[0] << ", " << system_momentum[1] << ")" << std::endl;
    std::cout << "System Energy: " << system_energy << std::endl;
    std::cout << "Particles: " << particles.size() << std::endl;
    // Output the runtime
     printf("Runtime: %f seconds\n", seconds);


    return 0;
}
