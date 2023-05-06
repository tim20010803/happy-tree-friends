#include <iostream>
#include <queue>
#include <vector>
#include <cmath>
struct Particle{
    std::vector<float> posi;
    std::vector<float> velocity;
    std::vector<float> acceleration;
    float mass;
};

void calculate_gravity(std::vector<Particle>& particles, float G) {
    for (auto& p1 : particles) {
        for (auto& p2 : particles) {
            if (&p1 == &p2) {
                continue; // Skip self-interaction
            }
            // Calculate distance between particles
            float dx = p2.posi[0] - p1.posi[0];
            float dy = p2.posi[1] - p1.posi[1];
            float dist_squared = dx*dx + dy*dy;
            float dist_cubed = dist_squared * std::sqrt(dist_squared);

            // Calculate gravitational force
            float force_magnitude = G * p1.mass * p2.mass / dist_cubed;
            float force_x = force_magnitude * dx;
            float force_y = force_magnitude * dy;

            // Update particle accelerations
            p1.acceleration[0] += force_x / p1.mass;
            p1.acceleration[1] += force_y / p1.mass;
        }
    }
}

std::vector<float> calculate_system_momentum(const std::vector<Particle>& particles) {
    std::vector<float> system_momentum(2, 0.0);
    for (const auto& p : particles) {
        system_momentum[0] += p.mass * p.velocity[0];
        system_momentum[1] += p.mass * p.velocity[1];
    }
    return system_momentum;
}

float calculate_system_energy(const std::vector<Particle>& particles, float G) {
    float total_kinetic_energy = 0.0;
    float total_potential_energy = 0.0;

    for (const auto& p : particles) {
        // Calculate kinetic energy
        float speed_squared = p.velocity[0]*p.velocity[0] + p.velocity[1]*p.velocity[1];
        float kinetic_energy = 0.5 * p.mass * speed_squared;
        total_kinetic_energy += kinetic_energy;

        // Calculate potential energy
        for (const auto& other_p : particles) {
            if (&p == &other_p) {
                continue;
            }
            float dx = other_p.posi[0] - p.posi[0];
            float dy = other_p.posi[1] - p.posi[1];
            float distance = std::sqrt(dx*dx + dy*dy);
            float potential_energy = -G * p.mass * other_p.mass / distance;
            total_potential_energy += potential_energy;
        }
    }

    return total_kinetic_energy + total_potential_energy;
}


void RK4(std::vector<Particle>& particles, float G, float dt) {
    for (auto& p : particles) {
        calculate_gravity(particles, G);
        // Get the current velocity and acceleration of the particle
        std::vector<float> current_velocity = p.velocity;
        std::vector<float> current_acceleration = p.acceleration;

        // Calculate the k1 values for velocity and position
        std::vector<float> k1_velocity = current_acceleration;
        std::vector<float> k1_position = current_velocity;

        // Calculate the k2 values for velocity and position
        std::vector<float> k2_velocity(p.velocity.size());
        std::vector<float> k2_position(p.posi.size());
        for (size_t i = 0; i < p.posi.size(); i++) {
            k2_velocity[i] = current_acceleration[i] + 0.5 * dt * k1_velocity[i];
            k2_position[i] = current_velocity[i] + 0.5 * dt * k1_position[i];
        }

        // Calculate the k3 values for velocity and position
        std::vector<float> k3_velocity(p.velocity.size());
        std::vector<float> k3_position(p.posi.size());
        for (size_t i = 0; i < p.posi.size(); i++) {
            k3_velocity[i] = current_acceleration[i] + 0.5 * dt * k2_velocity[i];
            k3_position[i] = current_velocity[i] + 0.5 * dt * k2_position[i];
        }

        // Calculate the k4 values for velocity and position
        std::vector<float> k4_velocity(p.velocity.size());
        std::vector<float> k4_position(p.posi.size());
        for (size_t i = 0; i < p.posi.size(); i++) {
            k4_velocity[i] = current_acceleration[i] + dt * k3_velocity[i];
            k4_position[i] = current_velocity[i] + dt * k3_position[i];
        }

        // Update the particle's position and velocity using the k values
        for (size_t i = 0; i < p.posi.size(); i++) {
            p.velocity[i] += (1.0 / 6.0) * dt * (k1_velocity[i] + 2.0 * k2_velocity[i] + 2.0 * k3_velocity[i] + k4_velocity[i]);
            p.posi[i] += (1.0 / 6.0) * dt * (k1_position[i] + 2.0 * k2_position[i] + 2.0 * k3_position[i] + k4_position[i]);
        }
    }
}


// main function is for testing
int main() {
    // Define simulation parameters
    const float G = 6.674e-11;
    std::vector<Particle> particles = {
        {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, 1.0},
        {{1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, 1.0}
    };

    // Input time and time step
    float t=0.1, dt=0.01;

    // Perform simulation
    int num_steps = t / dt;
    for (int i = 0; i < num_steps; i++) {
        RK4(particles, G, dt);
        std::vector<float> system_momentum = calculate_system_momentum(particles);
        float system_energy = calculate_system_energy(particles, G);

        std::cout << "Time: " << i*dt << std::endl;
        std::cout << "Particle 1 position: " << particles[0].posi[0] << ", " << particles[0].posi[1] << std::endl;
        std::cout << "Particle 2 position: " << particles[1].posi[0] << ", " << particles[1].posi[1] << std::endl;
        std::cout << "Particle 1 velocity: " << particles[0].velocity[0] << ", " << particles[0].velocity[1] << std::endl;
        std::cout << "Particle 2 velocity: " << particles[1].velocity[0] << ", " << particles[1].velocity[1] << std::endl;
        std::cout << "Particle 1 acceleration: " << particles[0].acceleration[0] << ", " << particles[0].acceleration[1] << std::endl;
        std::cout << "Particle 2 acceleration: " << particles[1].acceleration[0] << ", " << particles[1].acceleration[1] << std::endl;
        std::cout << "System momentum: " << system_momentum[0] << ", " << system_momentum[1] << std::endl;
        std::cout << "System energy: " << system_energy << std::endl;
    }
    return 0;
}
