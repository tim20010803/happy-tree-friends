#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <omp.h>

// Define particle structure
struct Particle {
    float posi[2];
    float velocity[2];
    float acceleration[2];
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

void calculate_gravity(std::vector<Particle>& particles, float G) {
    for (Particle& p1 : particles) {
        for (Particle& p2 : particles) {
            if (&p1 != &p2) {
                float dx = p2.posi[0] - p1.posi[0];
                float dy = p2.posi[1] - p1.posi[1];
                float dist_squared = dx * dx + dy * dy;
                float dist_cubed = dist_squared * sqrt(dist_squared);

                float force_magnitude = G * p1.mass * p2.mass / dist_cubed;
                float force_x = force_magnitude * dx;
                float force_y = force_magnitude * dy;

                p1.acceleration[0] += force_x / p1.mass;
                p1.acceleration[1] += force_y / p1.mass;
            }
        }
    }
}

void Verlet_velocity(std::vector<Particle>& particles, float G, float dt) {
    for (Particle& p : particles) {
        p.velocity[0] += p.acceleration[0] * dt;
        p.velocity[1] += p.acceleration[1] * dt;

        p.posi[0] += p.velocity[0] * dt;
        p.posi[1] += p.velocity[1] * dt;

        p.acceleration[0] = 0.0;
        p.acceleration[1] = 0.0;
    }
}

__global__ void calculate_gravity_cuda_kernel(Particle* particles, float G, int num_particles) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    if (index < num_particles) {
        Particle& p1 = particles[index];

        for (int j = 0; j < num_particles; j++) {
            if (j != index) {
                Particle& p2 = particles[j];
                float dx = p2.posi[0] - p1.posi[0];
                float dy = p2.posi[1] - p1.posi[1];
                float dist_squared = dx * dx + dy * dy;
                float dist_cubed = dist_squared * sqrt(dist_squared);

                float force_magnitude = G * p1.mass * p2.mass / dist_cubed;
                float force_x = force_magnitude * dx;
                float force_y = force_magnitude * dy;

                atomicAdd(&p1.acceleration[0], force_x / p1.mass);
                atomicAdd(&p1.acceleration[1], force_y / p1.mass);


            }
        }
    }
}

__global__ void Verlet_velocity_cuda_kernel(Particle* particles, float G, float dt, int num_particles) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    if (index < num_particles) {
        Particle& p = particles[index];

        p.velocity[0] += p.acceleration[0] * dt;
        p.velocity[1] += p.acceleration[1] * dt;

        p.posi[0] += p.velocity[0] * dt;
        p.posi[1] += p.velocity[1] * dt;

        p.acceleration[0] = 0.0;
        p.acceleration[1] = 0.0;
    }
}

int main() {
    // Define simulation parameters
    const float G = 6.674e-11;
    const float M_PI = 3.14159;

    std::cout << std::fixed << std::setprecision(12);

    // Input time and time step
    float t = 4.0 * M_PI / sqrt(G * 1000000000);
    float dt = 0.001;

    // Perform simulation
    int num_steps = t / dt;

    std::vector<Particle> particles = {
    { {0.0f, 1.0f}, {sqrtf(G*1000000000.0f)/2.0f, 0.0f}, {0.0f, 0.0f}, 1000000000.0f },
    { {0.0f, -1.0f}, {-sqrtf(G*1000000000.0f)/2.0f, 0.0f}, {0.0f, 0.0f}, 1000000000.0f }
    };


    float start_time = omp_get_wtime();

    // Perform Verlet simulation using CUDA
    Particle* d_particles;
    size_t particlesSize = particles.size() * sizeof(Particle);

    cudaMalloc((void**)&d_particles, particlesSize);
    cudaMemcpy(d_particles, particles.data(), particlesSize, cudaMemcpyHostToDevice);

    int blockSize = 256;
    int gridSize = (particles.size() + blockSize - 1) / blockSize;

    for (int i = 0; i <= num_steps; i++) {
        calculate_gravity_cuda_kernel<<<gridSize, blockSize>>>(d_particles, G, particles.size());
        Verlet_velocity_cuda_kernel<<<gridSize, blockSize>>>(d_particles, G, dt, particles.size());
    }

    cudaMemcpy(particles.data(), d_particles, particlesSize, cudaMemcpyDeviceToHost);
    cudaFree(d_particles);

    float end_time = omp_get_wtime();

    std::vector<float> system_momentum = calculate_system_momentum(particles);
    float system_energy = calculate_system_energy(particles, G);

    std::cout << "Verlet_velocity Time: " << (num_steps) * dt << std::endl;
    std::cout << "System Momentum (X, Y): (" << system_momentum[0] << ", " << system_momentum[1] << ")" << std::endl;
    std::cout << "System Energy: " << system_energy << std::endl;
    std::cout << "Total Execution Time: " << end_time - start_time << " seconds" << std::endl;

    return 0;
}
