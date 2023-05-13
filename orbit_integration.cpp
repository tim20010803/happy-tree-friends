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
    for (auto& p : particles) {
        // Update position using current velocity and acceleration
        p.posi[0] += p.velocity[0] * dt + 0.5 * p.acceleration[0] * dt * dt;
        p.posi[1] += p.velocity[1] * dt + 0.5 * p.acceleration[1] * dt * dt;
        
        // Calculate updated acceleration based on new positions
        p.acceleration_prev[0] = p.acceleration[0];
        p.acceleration_prev[1] = p.acceleration[1];

        p.acceleration[0] = 0.0;
        p.acceleration[1] = 0.0;
    }
        calculate_gravity(particles, G);
    for (auto& p : particles) {
        // Update velocity using average of old and new accelerations
        p.velocity[0] += 0.5 * (p.acceleration_prev[0] + p.acceleration[0]) * dt;
        p.velocity[1] += 0.5 * (p.acceleration_prev[1] + p.acceleration[1]) * dt;
    }
}

void AB(std::vector<Particle>& particles, double G, double dt) {
    for (auto& p : particles) {
        // Store original positions and velocities
        double original_vel_x = p.velocity[0];
        double original_vel_y = p.velocity[1];

        // Update positions using Adams-Bashforth
        p.posi[0] += dt * (1.5 * p.velocity[0] - 0.5 * p.velocity_prev[0]);
        p.posi[1] += dt * (1.5 * p.velocity[1] - 0.5 * p.velocity_prev[1]);

        // Update velocities
        p.velocity[0] += dt * (1.5 * p.acceleration[0] - 0.5 * p.acceleration_prev[0]);
        p.velocity[1] += dt * (1.5 * p.acceleration[1] - 0.5 * p.acceleration_prev[1]);

        // Store current accelerations and velocities for the next iteration
        p.acceleration_prev[0] = p.acceleration[0];
        p.acceleration_prev[1] = p.acceleration[1];
        p.velocity_prev[0] = original_vel_x;
        p.velocity_prev[1] = original_vel_y;

        // Reset accelerations for the next iteration
        p.acceleration[0] = 0.0;
        p.acceleration[1] = 0.0;
    }
    calculate_gravity(particles, G);
}