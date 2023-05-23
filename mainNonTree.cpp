
#include <iostream>
#include <queue>
#include <vector>
#include <cmath>
#include <iomanip>
#include "orbit_integration.h"

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
std::vector<double> calculate_system_momentum(const std::vector<Particle>& particles) {
    std::vector<double> system_momentum(2, 0.0);
    for (const auto& p : particles) {
        system_momentum[0] += p.mass * p.velocity[0];
        system_momentum[1] += p.mass * p.velocity[1];
    }
    return system_momentum;
}

double calculate_system_energy(const std::vector<Particle>& particles) {
  double total_kinetic_energy = 0.0;
  double total_potential_energy = 0.0;
  
    for (const auto& p : particles) {
        // Calculate kinetic energy
        double speed_squared = p.velocity[0]*p.velocity[0] + p.velocity[1]*p.velocity[1];
        double kinetic_energy = 0.5 * p.mass * speed_squared;
        total_kinetic_energy += kinetic_energy;

        // Calculate potential energy
        for (const auto& other_p : particles) {
            if (&p == &other_p) {
                continue;
            }
            double dx = other_p.posi[0] - p.posi[0];
            double dy = other_p.posi[1] - p.posi[1];
            double distance = std::sqrt(dx*dx + dy*dy);
            double potential_energy = -G_CONST * p.mass * other_p.mass / distance;
            total_potential_energy += potential_energy;
        }
    }

    return total_kinetic_energy + total_potential_energy;
}
// main function is for testing
int main() {
    clock_t start_t, end_t;
    start_t = clock();
    Particle a;
    a.posi = {1.,6.,4.};
    a.velocity = {4.,3.,2.};
    a.mass = {12.};
    a.acceleration = {0., 0., 0.};
    Particle b;
    b.posi = {2.,7.,8.};
    b.velocity = {1.,6.,7.};
    b.mass = {23.};
    b.acceleration = {0., 0., 0.};
    Particle c;
    c.posi = {3.,5.,8.};
    c.velocity = {1.,6.,7.};
    c.mass = {212.};
    c.acceleration = {0., 0., 0.};
    Particle d;
    d.posi = {4.,6.,8.};
    d.velocity = {8.,6.,7.};
    d.mass = {62.};
    d.acceleration = {0., 0., 0.};
    
    std::vector<Particle> Pvec = {a,b,c,d};
    const int particleNum = 100;//74.801941000000secondsif100
    for (int i = 0; i < particleNum-4; i++){
        Particle a;
        double rand3[3]={0};
        for (int j = 0; j < 3; j++){
            rand3[j] = static_cast<double>(rand() % 1000)/100.;
        }
        a.posi = {rand3[0],rand3[1],rand3[2]};
        for (int j = 0; j < 3; j++){
            rand3[j] = static_cast<double>(rand() % 1000)/100.;
        }
        a.velocity = {rand3[0],rand3[1],rand3[2]};
        a.mass = static_cast<double>(rand() % 1000)/10.;
        a.acceleration={0.,0.,0.};
        Pvec.push_back(a);
    }
    
    
    for (int i = 0; i < Pvec.size(); i++)
    {
        std::cout << "Particle mass: " << Pvec[i].mass << std::endl;
        std::cout << "Particle position: " << Pvec[i].posi[0] << ", " << Pvec[i].posi[1] << std::endl;
        std::cout << "Particle velocity: " << Pvec[i].velocity[0] << ", " << Pvec[i].velocity[1] << std::endl;
        std::cout << "Particle acceleration: " << Pvec[i].acceleration[0] << ", " << Pvec[i].acceleration[1] << std::endl;

    }
    std::cout << std::fixed << std::setprecision(12);

    std::vector<Particle> particles = Pvec;
    // Input time and time step
    double t=40.*M_PI/pow( G_CONST*1000000000,0.5), dt=0.001;

    // Perform simulation
    int num_steps = t / dt;
    
    calculate_gravity(particles, G_CONST);
    for (int i = 0; i <= num_steps; i++) {
        Verlet_velocity(particles,  G_CONST, dt);
    }
    end_t = clock();
    double total_t = static_cast<double>(end_t - start_t) / CLOCKS_PER_SEC;
    std::vector<double> system_momentum = calculate_system_momentum(particles);
    double system_energy = calculate_system_energy(particles);
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
    std::cout << total_t <<  "seconds"<< std::endl;
    std::cout << particleNum <<  "particles"<< std::endl;

    return 0;
}
