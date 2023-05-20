#include <iostream>
#include <queue>
#include <vector>
#include <cmath>
#include <iomanip>
#include <omp.h>
#include <algorithm>

struct Particle{
    std::vector<double> posi;
    std::vector<double> velocity;
    std::vector<double> acceleration;

    double mass;
};

void calculate_gravity(std::vector<Particle>& particles, double G) {
    #pragma omp parallel for
    for (int i = 0; i < particles.size(); i++) {
        Particle& p1 = particles[i];
        for (int j = 0; j < particles.size(); j++) {
            if (i == j) {
                continue; // Skip self-interaction
            }
            Particle& p2 = particles[j];
            
            // Calculate distance between particles
            double dx = p2.posi[0] - p1.posi[0];
            double dy = p2.posi[1] - p1.posi[1];
            double dist_squared = dx * dx + dy * dy;
            double dist_cubed = dist_squared * std::sqrt(dist_squared);

            // Calculate gravitational force
            double force_magnitude = G * p1.mass * p2.mass / dist_cubed;
            double force_x = force_magnitude * dx;
            double force_y = force_magnitude * dy;

            // Update particle accelerations
            #pragma omp atomic
            p1.acceleration[0] += force_x / p1.mass;
            #pragma omp atomic
            p1.acceleration[1] += force_y / p1.mass;
        }
    }
}

// 自定义reduction运算符：向量相加
#pragma omp declare reduction(vec_double_plus : std::vector<double> : \
    std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
    initializer(omp_priv = decltype(omp_orig)(omp_orig.size(), 0.0))

// 计算系统动量
std::vector<double> calculate_system_momentum(const std::vector<Particle>& particles) {
    std::vector<double> system_momentum(2, 0.0);

    #pragma omp parallel for reduction(vec_double_plus:system_momentum)
    for (int i = 0; i < particles.size(); ++i) {
        const Particle& particle = particles[i];
        system_momentum[0] += particle.mass * particle.velocity[0];
        system_momentum[1] += particle.mass * particle.velocity[1];
    }

    return system_momentum;
}

// 计算系统角动量
std::vector<double> calculate_system_angular_momentum(const std::vector<Particle>& particles) {
    std::vector<double> system_angular_momentum(1, 0.0);

    #pragma omp parallel for reduction(vec_double_plus:system_angular_momentum)
    for (int i = 0; i < particles.size(); ++i) {
        const Particle& particle = particles[i];
        system_angular_momentum[0] += particle.mass * (particle.posi[0] * particle.velocity[1] - particle.posi[1] * particle.velocity[0]);
    }

    return system_angular_momentum;
}


double calculate_system_energy(const std::vector<Particle>& particles, double G) {
    double total_kinetic_energy = 0.0;
    double total_potential_energy = 0.0;

    #pragma omp parallel for reduction(+:total_kinetic_energy, total_potential_energy)
    for (size_t i = 0; i < particles.size(); i++) {
        const auto& p = particles[i];

        // Calculate kinetic energy
        double speed_squared = p.velocity[0] * p.velocity[0] + p.velocity[1] * p.velocity[1];
        double kinetic_energy = 0.5 * p.mass * speed_squared;
        #pragma omp atomic
        total_kinetic_energy += kinetic_energy;

        // Calculate potential energy
        for (size_t j = 0; j < particles.size(); j++) {
            if (i == j) {
                continue;
            }

            const auto& other_p = particles[j];
            double dx = other_p.posi[0] - p.posi[0];
            double dy = other_p.posi[1] - p.posi[1];
            double distance = std::sqrt(dx * dx + dy * dy);
            double potential_energy = -G * p.mass * other_p.mass / distance;
            #pragma omp atomic
            total_potential_energy += potential_energy;
        }
    }

    return total_kinetic_energy + total_potential_energy;
}

void RK4(std::vector<Particle>& particles, double G, double dt) {
    #pragma omp parallel for
    for (int i = 0; i < particles.size(); i++) {
        Particle& p = particles[i];
        // Get the current velocity and acceleration of the particle
        std::vector<double> current_velocity = p.velocity;
        std::vector<double> current_acceleration = p.acceleration;

        // Calculate the k1 values for velocity and position
        std::vector<double> k1_velocity = current_acceleration;
        std::vector<double> k1_position = current_velocity;

        // Calculate the k2 values for velocity and position
        std::vector<double> k2_velocity(p.velocity.size());
        std::vector<double> k2_position(p.posi.size());
        for (size_t j = 0; j < p.posi.size(); j++) {
            k2_velocity[j] = current_acceleration[j] + 0.5 * dt * k1_velocity[j];
            k2_position[j] = current_velocity[j] + 0.5 * dt * k1_position[j];
        }

        // Calculate the k3 values for velocity and position
        std::vector<double> k3_velocity(p.velocity.size());
        std::vector<double> k3_position(p.posi.size());
        for (size_t j = 0; j < p.posi.size(); j++) {
            k3_velocity[j] = current_acceleration[j] + 0.5 * dt * k2_velocity[j];
            k3_position[j] = current_velocity[j] + 0.5 * dt * k2_position[j];
        }

        // Calculate the k4 values for velocity and position
        std::vector<double> k4_velocity(p.velocity.size());
        std::vector<double> k4_position(p.posi.size());
        for (size_t j = 0; j < p.posi.size(); j++) {
            k4_velocity[j] = current_acceleration[j] + dt * k3_velocity[j];
            k4_position[j] = current_velocity[j] + dt * k3_position[j];
        }

        // Update the particle's position and velocity using the k values
        for (size_t j = 0; j < p.posi.size(); j++) {
            p.velocity[j] += (1.0 / 6.0) * dt * (k1_velocity[j] + 2.0 * k2_velocity[j] + 2.0 * k3_velocity[j] + k4_velocity[j]);
            p.posi[j] += (1.0 / 6.0) * dt * (k1_position[j] + 2.0 * k2_position[j] + 2.0 * k3_position[j] + k4_position[j]);
        }
        p.acceleration[0] = 0.0;
        p.acceleration[1] = 0.0;
    }
    calculate_gravity(particles, G);
}

void Verlet_velocity(std::vector<Particle>& particles, double G, double dt) {
    std::vector<std::vector<double>> acceleration_prevs(particles.size(), std::vector<double>(2, 0.0));

    acceleration_prevs.reserve(particles.size());

    #pragma omp parallel for
    for (int i = 0; i < particles.size(); i++) {
        Particle& p = particles[i];
        // Update position using current velocity and acceleration
        p.posi[0] += p.velocity[0] * dt + 0.5 * p.acceleration[0] * dt * dt;
        p.posi[1] += p.velocity[1] * dt + 0.5 * p.acceleration[1] * dt * dt;

        // Save the previous acceleration in the vector
        acceleration_prevs[i] = p.acceleration;

        p.acceleration[0] = 0.0;
        p.acceleration[1] = 0.0;
    }

    calculate_gravity(particles, G);

    // Use the saved acceleration_prevs to update velocities
    #pragma omp parallel for
    for (int i = 0; i < particles.size(); i++) {
        Particle& p = particles[i];
        // Update velocity using average of old and new accelerations
        p.velocity[0] += 0.5 * (acceleration_prevs[i][0] + p.acceleration[0]) * dt;
        p.velocity[1] += 0.5 * (acceleration_prevs[i][1] + p.acceleration[1]) * dt;
    }
}



void AB(std::vector<Particle>& particles, double G, double dt) {
    #pragma omp parallel for
    for (int i = 0; i < particles.size(); i++) {
        Particle& p = particles[i];
        // Store original positions and velocities
        std::vector<double> prev_velocity = p.velocity;
        std::vector<double> prev_acceleration = p.acceleration;

        // Update positions using Adams-Bashforth
        p.posi[0] += dt * (1.5 * p.velocity[0] - 0.5 * prev_velocity[0]);
        p.posi[1] += dt * (1.5 * p.velocity[1] - 0.5 * prev_velocity[1]);

        // Update velocities
        p.velocity[0] += dt * (1.5 * p.acceleration[0] - 0.5 * prev_acceleration[0]);
        p.velocity[1] += dt * (1.5 * p.acceleration[1] - 0.5 * prev_acceleration[1]);

        // Reset accelerations for the next iteration
        p.acceleration[0] = 0.0;
        p.acceleration[1] = 0.0;
    }
    calculate_gravity(particles, G);
}



// main function is for testing
int main() {
    // Define simulation parameters
    const double G = 6.674e-11;
    std::cout << std::fixed << std::setprecision(12);

    // Input time and time step
    double t=4.*M_PI/pow(G*1000000000,0.5), dt=0.001;

    // Perform simulation
    int num_steps = t / dt;
    
    //RK4
    std::vector<Particle> particles = {
        {{0.0, 1.0}, {pow(G*1000000000,0.5)/2., 0.0}, {0.0, 0.0}, 1000000000},
        {{0.0, -1.0}, {-pow(G*1000000000,0.5)/2., 0.0}, {0.0, 0.0}, 1000000000}
    };
        std::cout << "theoretical energy: " << calculate_system_energy(particles, G) << std::endl;

    calculate_gravity(particles, G);
    for (int i = 0; i <= num_steps; i++) {
        RK4(particles, G, dt);        
    }
    std::vector<double> system_momentum = calculate_system_momentum(particles);
    double system_energy = calculate_system_energy(particles, G);
    std::cout << "RK4 Time: " << (num_steps)*dt << std::endl;
    std::cout << "Particle 1 mass: " << particles[0].mass << std::endl;
    std::cout << "Particle 2 mass: " << particles[1].mass << std::endl;
    std::cout << "Particle 1 position: " << particles[0].posi[0] << ", " << particles[0].posi[1] << std::endl;
    std::cout << "Particle 2 position: " << particles[1].posi[0] << ", " << particles[1].posi[1] << std::endl;
    std::cout << "Particle 1 velocity: " << particles[0].velocity[0] << ", " << particles[0].velocity[1] << std::endl;
    std::cout << "Particle 2 velocity: " << particles[1].velocity[0] << ", " << particles[1].velocity[1] << std::endl;
    std::cout << "Particle 1 acceleration: " << particles[0].acceleration[0] << ", " << particles[0].acceleration[1] << std::endl;
    std::cout << "Particle 2 acceleration: " << particles[1].acceleration[0] << ", " << particles[1].acceleration[1] << std::endl;
    std::cout << "System momentum: " << system_momentum[0] << ", " << system_momentum[1] << std::endl;
    std::cout << "System energy: " << system_energy << std::endl;

    std::cout<<"//////////////////////////////////"<<std::endl;

    //verlet
    std::vector<Particle> particles = {
        {{0.0, 1.0}, {pow(G*1000000000,0.5)/2., 0.0}, {0.0, 0.0}, 1000000000},
        {{0.0, -1.0}, {-pow(G*1000000000,0.5)/2., 0.0}, {0.0, 0.0}, 1000000000}
    };
    calculate_gravity(particles, G);
    for (int i = 0; i <= num_steps; i++) {
        Verlet_velocity(particles, G, dt);
    }
    std::vector<double> system_momentum = calculate_system_momentum(particles);
    double system_energy = calculate_system_energy(particles, G);
    std::cout << "Verlet_velocity Time: " << (num_steps)*dt << std::endl;
    std::cout << "Particle 1 mass: " << particles[0].mass << std::endl;
    std::cout << "Particle 2 mass: " << particles[1].mass << std::endl;
    std::cout << "Particle 1 position: " << particles[0].posi[0] << ", " << particles[0].posi[1] << std::endl;
    std::cout << "Particle 2 position: " << particles[1].posi[0] << ", " << particles[1].posi[1] << std::endl;
    std::cout << "Particle 1 velocity: " << particles[0].velocity[0] << ", " << particles[0].velocity[1] << std::endl;
    std::cout << "Particle 2 velocity: " << particles[1].velocity[0] << ", " << particles[1].velocity[1] << std::endl;
    std::cout << "Particle 1 acceleration: " << particles[0].acceleration[0] << ", " << particles[0].acceleration[1] << std::endl;
    std::cout << "Particle 2 acceleration: " << particles[1].acceleration[0] << ", " << particles[1].acceleration[1] << std::endl;
    std::cout << "System momentum: " << system_momentum[0] << ", " << system_momentum[1] << std::endl;
    std::cout << "System energy: " << system_energy << std::endl;
    
    
    //AB
    particles = {
        {{0.0, 1.0}, {pow(G*1000000000,0.5)/2., 0.0}, {0.0, 0.0}, 1000000000},
        {{0.0, -1.0}, {-pow(G*1000000000,0.5)/2., 0.0}, {0.0, 0.0}, 1000000000}
    };
    calculate_gravity(particles, G);
    for (int i = 0; i <= num_steps; i++) {
        AB(particles, G, dt);
    }
    system_momentum = calculate_system_momentum(particles);
    system_energy = calculate_system_energy(particles, G);
    std::cout << "AB Time: " << (num_steps)*dt << std::endl;
    std::cout << "Particle 1 mass: " << particles[0].mass << std::endl;
    std::cout << "Particle 2 mass: " << particles[1].mass << std::endl;
    std::cout << "Particle 1 position: " << particles[0].posi[0] << ", " << particles[0].posi[1] << std::endl;
    std::cout << "Particle 2 position: " << particles[1].posi[0] << ", " << particles[1].posi[1] << std::endl;
    std::cout << "Particle 1 velocity: " << particles[0].velocity[0] << ", " << particles[0].velocity[1] << std::endl;
    std::cout << "Particle 2 velocity: " << particles[1].velocity[0] << ", " << particles[1].velocity[1] << std::endl;
    std::cout << "Particle 1 acceleration: " << particles[0].acceleration[0] << ", " << particles[0].acceleration[1] << std::endl;
    std::cout << "Particle 2 acceleration: " << particles[1].acceleration[0] << ", " << particles[1].acceleration[1] << std::endl;
    std::cout << "System momentum: " << system_momentum[0] << ", " << system_momentum[1] << std::endl;
    std::cout << "System energy: " << system_energy << std::endl;
    

    return 0;
}
