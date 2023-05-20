#include <iostream>
#include <queue>
#include <vector>
#include <cmath>
#include <iomanip>
#include <omp.h>
#include "orbit_integration_OMP.h"
#include "quadrupleTree.h"

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


void RK4_Tree(std::vector<Particle>& particles, double G, double dt,double mX,double mY,double mZ,double MX,double MY, double MZ) {
    QuadrupleTree quadTree(particles,mX,mY,mZ,MX,MY,MZ);
    #pragma omp parallel for
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
    quadTree.TreeForce();
    // std::cout<<"particle 1 position:("<<particles[0].posi[0]<<","<<particles[0].posi[1]<<")\n";
}


void Verlet_velocity_Tree(std::vector<Particle>& particles, double G, double dt,double mX,double mY,double mZ,double MX,double MY, double MZ) {
    QuadrupleTree quadTree(particles,mX,mY,mZ,MX,MY,MZ);
    std::vector<std::vector<double>> acceleration_prevs(particles.size(), std::vector<double>(2, 0.0));
    acceleration_prevs.reserve(particles.size());
    #pragma omp parallel for
    for (auto& p : particles) {
        // Update position using current velocity and acceleration
        p.posi[0] += p.velocity[0] * dt + 0.5 * p.acceleration[0] * dt * dt;
        p.posi[1] += p.velocity[1] * dt + 0.5 * p.acceleration[1] * dt * dt;

        // Save the previous acceleration in the vector
        acceleration_prevs.push_back(p.acceleration);

        p.acceleration[0] = 0.0;
        p.acceleration[1] = 0.0;
    }
    // std::cout<<"particle 1 position:("<<particles[0].posi[0]<<","<<particles[0].posi[1]<<")\n";
    quadTree.TreeForce();
    // Use the saved acceleration_prevs to update velocities
    #pragma omp parallel for
    for (int i = 0; i < particles.size(); i++) {
        Particle& p = particles[i];
        // Update velocity using average of old and new accelerations
        p.velocity[0] += 0.5 * (acceleration_prevs[i][0] + p.acceleration[0]) * dt;
        p.velocity[1] += 0.5 * (acceleration_prevs[i][1] + p.acceleration[1]) * dt;
    }
}
void AB_Tree(std::vector<Particle>& particles, double G, double dt,double mX,double mY,double mZ,double MX,double MY, double MZ) {
    QuadrupleTree quadTree(particles,mX,mY,mZ,MX,MY,MZ);
    #pragma omp parallel for
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
    quadTree.TreeForce();
    // std::cout<<"particle 1 position:("<<particles[0].posi[0]<<","<<particles[0].posi[1]<<")\n";
}