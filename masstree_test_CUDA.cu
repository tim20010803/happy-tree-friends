#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <omp.h>
#include <random>
#include <chrono>


// Define particle structure
struct Particle {
    float posi[2];
    float velocity[2];
    float acceleration[2];
    float acceleration_prev[2];
    float mass;
};

// Calculate system momentum
std::vector<float> calculate_system_momentum(const std::vector<Particle>& particles) {
    float total_momentum_x = 0.0;
    float total_momentum_y = 0.0;

    for (const Particle& p : particles) {
        total_momentum_x += p.mass * p.velocity[0];
        total_momentum_y += p.mass * p.velocity[1];
    }

    return { total_momentum_x, total_momentum_y };
}

// Calculate system energy
float calculate_system_energy(const std::vector<Particle>& particles, float G) {
    float total_energy = 0.0;

    for (const Particle& p : particles) {
        float kinetic_energy = 0.5 * p.mass * (p.velocity[0] * p.velocity[0] + p.velocity[1] * p.velocity[1]);
        float potential_energy = 0.0;

        for (const Particle& other : particles) {
            if (&p != &other) {
                float dx = other.posi[0] - p.posi[0];
                float dy = other.posi[1] - p.posi[1];
                float dist = sqrt(dx * dx + dy * dy);
                potential_energy -= G * p.mass * other.mass / dist;
            }
        }

        total_energy += kinetic_energy + potential_energy;
    }

    return total_energy;
}



__global__ void simulate_particles_cuda_kernel(Particle* particles, const float* particle_masses, double G, double dt, int num_particles) {
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
    const float G = 6.674e-11;
    const float M_PI = 3.14159;

    std::cout << std::fixed << std::setprecision(12);

    // Input time and time step
    float t = 0.001;
    float dt = 0.001;

    // Perform simulation
    int num_steps = t / dt;

    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dist_pos(-1.0f, 1.0f);
    std::uniform_real_distribution<float> dist_vel(-10.0f, 10.0f);
    std::uniform_real_distribution<float> dist_mass(1.0e6f, 1.0e7f);

    // Generate particles
    std::vector<Particle> particles;
    particles.reserve(131072);
    for (int i = 0; i < 131072; ++i) {
        Particle particle;
        particle.posi[0] = dist_pos(gen);
        particle.posi[1] = dist_pos(gen);
        particle.velocity[0] = dist_vel(gen);
        particle.velocity[1] = dist_vel(gen);
        particle.acceleration[0] = 0.0f;
        particle.acceleration[1] = 0.0f;
        particle.mass = dist_mass(gen);
        particles.push_back(particle);
    }

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
    simulate_particles_cuda_kernel<<<gridSize, blockSize, sharedMemSize>>>(d_particles, d_masses, G, dt, particles.size());

    // Transfer particles back to host memory
    cudaMemcpy(particles.data(), d_particles, particlesSize, cudaMemcpyDeviceToHost);

    // Cleanup
    cudaFree(d_particles);
    cudaFree(d_masses);

    // End timer
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate duration
    std::chrono::duration<double> duration = end - start;
    double seconds = duration.count();

    std::vector<float> system_momentum = calculate_system_momentum(particles);
    float system_energy = calculate_system_energy(particles, G);

    std::cout << "Verlet_velocity Time: " << (num_steps) * dt << std::endl;
    std::cout << "System Momentum (X, Y): (" << system_momentum[0] << ", " << system_momentum[1] << ")" << std::endl;
    std::cout << "System Energy: " << system_energy << std::endl;
    std::cout << "Particles: " << particles.size() << std::endl;
    // Output the runtime
     printf("Runtime: %f seconds\n", seconds);


    return 0;
}

