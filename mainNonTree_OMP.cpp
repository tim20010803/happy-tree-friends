
#include <iostream>
#include <queue>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include "orbit_integration_OMP.h"

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

// // main function is for testing
// int main() {
//     Particle a;
//     a.posi = {1.,6.,4.};
//     a.velocity = {4.,3.,2.};
//     a.mass = {12.};
//     a.acceleration = {0., 0., 0.};
//     Particle b;
//     b.posi = {2.,7.,8.};
//     b.velocity = {1.,6.,7.};
//     b.mass = {23.};
//     b.acceleration = {0., 0., 0.};
//     Particle c;
//     c.posi = {3.,5.,8.};
//     c.velocity = {1.,6.,7.};
//     c.mass = {212.};
//     c.acceleration = {0., 0., 0.};
//     Particle d;
//     d.posi = {4.,6.,8.};
//     d.velocity = {8.,6.,7.};
//     d.mass = {62.};
//     d.acceleration = {0., 0., 0.};
    
//     std::vector<Particle> Pvec = {a,b,c,d};
 
//     std::cout << std::fixed << std::setprecision(12);

//     std::vector<Particle> particles = Pvec;
//     // Input time and time step
//     double t=40.*M_PI/pow( G_CONST*1000000000,0.5), dt=0.001;

//     // Perform simulation
//     int num_steps = t / dt;
    
//     calculate_gravity(particles, G_CONST);
//     for (int i = 0; i <= num_steps; i++) {
//         Verlet_velocity(particles,  G_CONST, dt);
//     }

//     std::vector<double> system_momentum = calculate_system_momentum(particles);
//     double system_energy = calculate_system_energy(particles, G_CONST);
//     std::cout << "Verlet_velocity Time: " << (num_steps)*dt << std::endl;
//     std::cout << "Particle 1 mass: " << particles[0].mass << std::endl;
//     std::cout << "Particle 2 mass: " << particles[1].mass << std::endl;
//     std::cout << "Particle 1 position: " << particles[0].posi[0] << ", " << particles[0].posi[1] << std::endl;
//     std::cout << "Particle 2 position: " << particles[1].posi[0] << ", " << particles[1].posi[1] << std::endl;
//     std::cout << "Particle 1 velocity: " << particles[0].velocity[0] << ", " << particles[0].velocity[1] << std::endl;
//     std::cout << "Particle 2 velocity: " << particles[1].velocity[0] << ", " << particles[1].velocity[1] << std::endl;
//     std::cout << "Particle 1 acceleration: " << particles[0].acceleration[0] << ", " << particles[0].acceleration[1] << std::endl;
//     std::cout << "Particle 2 acceleration: " << particles[1].acceleration[0] << ", " << particles[1].acceleration[1] << std::endl;
//     std::cout << "System momentum: " << system_momentum[0] << ", " << system_momentum[1] << std::endl;
//     std::cout << "System energy: " << system_energy << std::endl;

//     return 0;
// }
