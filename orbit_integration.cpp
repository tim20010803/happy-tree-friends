#include <iostream>
#include <queue>
#include <vector>
#include <cmath>
#include <iomanip>

struct Particle{
    std::vector<double> posi;
    std::vector<double> velocity;
    std::vector<double> acceleration;

    double mass;
};

void calculate_gravity(std::vector<Particle>& particles, double G) {
    for (auto& p1 : particles) {
        for (auto& p2 : particles) {
            if (&p1 == &p2) {
                continue; // Skip self-interaction
            }
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
            p1.acceleration[0] += force_x / p1.mass;
            p1.acceleration[1] += force_y / p1.mass;
        }
    }
}

void RK4(std::vector<Particle>& particles, double G, double dt) {
    for (auto& p : particles) {
        // Get the current velocity and acceleration of the particle
        std::vector<double> current_velocity = p.velocity;
        std::vector<double> current_acceleration = p.acceleration;

        // Calculate the k1 values for velocity and position
        std::vector<double> k1_velocity = current_acceleration;
        std::vector<double> k1_position = current_velocity;

        // Calculate the k2 values for velocity and position
        std::vector<double> k2_velocity(p.velocity.size());
        std::vector<double> k2_position(p.posi.size());
        for (size_t i = 0; i < p.posi.size(); i++) {
            k2_velocity[i] = current_acceleration[i] + 0.5 * dt * k1_velocity[i];
            k2_position[i] = current_velocity[i] + 0.5 * dt * k1_position[i];
        }

        // Calculate the k3 values for velocity and position
        std::vector<double> k3_velocity(p.velocity.size());
        std::vector<double> k3_position(p.posi.size());
        for (size_t i = 0; i < p.posi.size(); i++) {
            k3_velocity[i] = current_acceleration[i] + 0.5 * dt * k2_velocity[i];
            k3_position[i] = current_velocity[i] + 0.5 * dt * k2_position[i];
        }

        // Calculate the k4 values for velocity and position
        std::vector<double> k4_velocity(p.velocity.size());
        std::vector<double> k4_position(p.posi.size());
        for (size_t i = 0; i < p.posi.size(); i++) {
            k4_velocity[i] = current_acceleration[i] + dt * k3_velocity[i];
            k4_position[i] = current_velocity[i] + dt * k3_position[i];
        }

        // Update the particle's position and velocity using the k values
        for (size_t i = 0; i < p.posi.size(); i++) {
            p.velocity[i] += (1.0 / 6.0) * dt * (k1_velocity[i] + 2.0 * k2_velocity[i] + 2.0 * k3_velocity[i] + k4_velocity[i]);
            p.posi[i] += (1.0 / 6.0) * dt * (k1_position[i] + 2.0 * k2_position[i] + 2.0 * k3_position[i] + k4_position[i]);
        }
        p.acceleration[0] = 0.0;
        p.acceleration[1] = 0.0;
    }
    calculate_gravity(particles, G);
}

void Verlet_velocity(std::vector<Particle>& particles, double G, double dt) {
    std::vector<std::vector<double>> acceleration_prevs;
    acceleration_prevs.reserve(particles.size());

    for (auto& p : particles) {
        // Update position using current velocity and acceleration
        p.posi[0] += p.velocity[0] * dt + 0.5 * p.acceleration[0] * dt * dt;
        p.posi[1] += p.velocity[1] * dt + 0.5 * p.acceleration[1] * dt * dt;

        // Save the previous acceleration in the vector
        acceleration_prevs.push_back(p.acceleration);

        p.acceleration[0] = 0.0;
        p.acceleration[1] = 0.0;
    }

    calculate_gravity(particles, G);

    // Use the saved acceleration_prevs to update velocities
    auto prev_iter = acceleration_prevs.begin();
    for (auto& p : particles) {
        // Update velocity using average of old and new accelerations
        p.velocity[0] += 0.5 * ((*prev_iter)[0] + p.acceleration[0]) * dt;
        p.velocity[1] += 0.5 * ((*prev_iter)[1] + p.acceleration[1]) * dt;

        ++prev_iter;
    }
}
void AB(std::vector<Particle>& particles, double G, double dt) {
    for (auto& p : particles) {
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