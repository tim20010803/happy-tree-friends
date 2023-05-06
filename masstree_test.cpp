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


void perform_rk45_step(std::vector<Particle>& particles, float G, float dt) {
    // Define coefficients for RK45
    std::vector<float> a = {0.0, 1.0/4.0, 3.0/8.0, 12.0/13.0, 1.0, 1.0/2.0};
    std::vector<float> b = {25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0, 0.0};
    std::vector<float> c = {1.0/360.0, 0.0, -128.0/4275.0, -2197.0/75240.0, 1.0/50.0, 2.0/55.0};

    // Perform RK45 step
    std::vector<Particle> orig_particles = particles;
    std::vector<std::vector<float>> k(6, std::vector<float>(4*particles.size()));
    for (int i = 0; i < 6; i++) {
        particles = orig_particles;
        for (int j = 0; j < particles.size(); j++) {
            particles[j].posi[0] += a[i]*dt*particles[j].velocity[0];
            particles[j].posi[1] += a[i]*dt*particles[j].velocity[1];
            particles[j].velocity[0] += a[i]*dt*particles[j].acceleration[0];
            particles[j].velocity[1] += a[i]*dt*particles[j].acceleration[1];
        }
        calculate_gravity(particles, G);
        for (int j = 0; j < particles.size(); j++) {
            k[i][j] = particles[j].velocity[0];
            k[i][j+particles.size()] = particles[j].velocity[1];
            k[i][j+2*particles.size()] = particles[j].acceleration[0];
            k[i][j+3*particles.size()] = particles[j].acceleration[1];
        }
    }
    // Calculate new particle positions and velocities
    for (int i = 0; i < particles.size(); i++) {
        particles[i].posi[0] += dt*(b[0]*k[0][i] + b[1]*k[1][i] + b[2]*k[2][i] + b[3]*k[3][i] + b[4]*k[4][i] + b[5]*k[5][i]);
        particles[i].posi[1] += dt*(b[0]*k[0][i+particles.size()] + b[1]*k[1][i+particles.size()] + b[2]*k[2][i+particles.size()] + b[3]*k[3][i+particles.size()] + b[4]*k[4][i+particles.size()] + b[5]*k[5][i+particles.size()]);
        particles[i].velocity[0] += dt*(c[0]*k[0][i+2*particles.size()] + c[1]*k[1][i+2*particles.size()] + c[2]*k[2][i+2*particles.size()] + c[3]*k[3][i+2*particles.size()] + c[4]*k[4][i+2*particles.size()] + c[5]*k[5][i+2*particles.size()]);
        particles[i].velocity[1] += dt*(c[0]*k[0][i+3*particles.size()] + c[1]*k[1][i+3*particles.size()] + c[2]*k[2][i+3*particles.size()] + c[3]*k[3][i+3*particles.size()] + c[4]*k[4][i+3*particles.size()] + c[5]*k[5][i+3*particles.size()]);
        particles[i].acceleration[0] += dt*(c[0]*k[0][i+4*particles.size()] + c[1]*k[1][i+4*particles.size()] + c[2]*k[2][i+4*particles.size()] + c[3]*k[3][i+4*particles.size()] + c[4]*k[4][i+4*particles.size()] + c[5]*k[5][i+4*particles.size()]);
        particles[i].acceleration[1] += dt*(c[0]*k[0][i+5*particles.size()] + c[1]*k[1][i+5*particles.size()] + c[2]*k[2][i+5*particles.size()] + c[3]*k[3][i+5*particles.size()] + c[4]*k[4][i+5*particles.size()] + c[5]*k[5][i+5*particles.size()]);
    }
}


// main function is for testing
int main() {
    // Define simulation parameters
    const float G = 6.674e-11;
    std::vector<Particle> particles = {
        {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, 1.0},
        {{10.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, 1.0}
    };

    // Input time and time step
    float t=1., dt=0.01;

    // Perform simulation
    int num_steps = t / dt;
    for (int i = 0; i < num_steps; i++) {
        perform_rk45_step(particles, G, dt);
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
