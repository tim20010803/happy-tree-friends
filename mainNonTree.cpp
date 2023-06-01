
#include <iostream>
#include <queue>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include "orbit_integration.h"
#include "orbit_integration_OMP.h"
#include <random>
#include <fstream>
// main function is for testing
int main() {
    clock_t start_t, end_t;
    // double starttime =omp_get_wtime();
    start_t = clock();
    
    int particleNum = 100; // number of particle
    double r = 1000; // radius of the initial disk
    std::vector<Particle> particles;
    
    // produce random seed
    std::random_device rd;
    std::mt19937 gen(rd());
    // set random range
    std::uniform_real_distribution<double> massDist(1.0, 10.0);
    std::uniform_real_distribution<double> xDist(-1000.0, 1000.0); 
    std::uniform_real_distribution<double> yDist(-1000.0, 1000.0);
    std::uniform_real_distribution<double> velocityXDist(-1000.0, 1000.0);
    std::uniform_real_distribution<double> velocityYDist(-1000.0, 1000.0);

    while ( particles.size() < particleNum ){ 
        // set parameter  
        double m = massDist(gen);
        double x = xDist(gen);
        double y = yDist(gen);
        double delta_vx = velocityXDist(gen);
        double delta_vy = velocityYDist(gen);
        double d = sqrt(x*x + y*y);
        // let the particles in the range of a disk
        if ( d <=  r ){
            Particle particle;
            
            // compute the velocity and the acceleratiion
            double vx = - d * y / d + delta_vx;
            double vy = d * x / d + delta_vy;
            double a = m / (d * d);
            double ax = - a * x / d ;
            double ay = - a * x / d;

            // add the particle initial condition
            particle.mass = m;
            particle.posi.push_back(x);
            particle.posi.push_back(y);
            particle.velocity.push_back(vx);
            particle.velocity.push_back(vy);
            particle.acceleration.push_back(ax);
            particle.acceleration.push_back(ay);

            // add the particle
            particles.push_back(particle);
        }
    }


    // Input time and time step
    double t=0.04*M_PI/pow( G_CONST*1000000000,0.5), dt=0.001;

    // Perform simulation
    int num_steps = 10000;
    calculate_gravity(particles);
    // make a new file
    std::ofstream file("mainNonTree_data.csv");
    file << "Time,Particle,Mass,PositionX,PositionY,VectorX,VectorY,AccelerationX,AccelerationY,SystemEnergy,SystemMomentumX,SystemMomentumY,SystemAngularMomentumX,SystemAngularMomentumY" << std::endl;
    
    // put the data into the file
    for (int time = 0; time <= num_steps; time += 1000){
        Verlet_velocity(particles, dt);
        for (int i = 0; i < particles.size(); i++)
        {
            double p = i+1;
            double mass = particles[i].mass;
            double posX = particles[i].posi[0];
            double posY = particles[i].posi[1];
            double vecX = particles[i].velocity[0];
            double vecY = particles[i].velocity[1];
            double accX = particles[i].acceleration[0];
            double accY = particles[i].acceleration[1];
            double system_energy = calculate_system_energy(particles);
            std::vector<double> system_momentum = calculate_system_momentum(particles);
            std::vector<double> system_angular_momentum = calculate_system_angular_momentum(particles);


            file << time << "," << p  << "," << mass << "," 
                << posX << "," << posY << "," << vecX << "," << vecY << "," 
                << accX << "," << accY << "," << system_energy << "," 
                << system_momentum[0] << "," << system_momentum[1] << "," 
                << system_angular_momentum[0] << "," << system_angular_momentum[1] << std::endl;
        }
    }
    file.close();
    
    end_t = clock();
    std::vector<double> system_momentum = calculate_system_momentum(particles);
    double system_energy = calculate_system_energy(particles);
    // double endtime =omp_get_wtime();
    double total_t = static_cast<double>(end_t - start_t) / CLOCKS_PER_SEC;
    //std::cout << "Physical Time: " << (num_steps)*dt <<"seconds"<< std::endl;
    //std::cout << "Particle 1 mass: " << particles[0].mass << std::endl;
    //std::cout << "Particle 2 mass: " << particles[1].mass << std::endl;
    //std::cout << "Particle 1 position: " << particles[0].posi[0] << ", " << particles[0].posi[1] << std::endl;
    //std::cout << "Particle 2 position: " << particles[1].posi[0] << ", " << particles[1].posi[1] << std::endl;
    //std::cout << "Particle 1 velocity: " << particles[0].velocity[0] << ", " << particles[0].velocity[1] << std::endl;
    //std::cout << "Particle 2 velocity: " << particles[1].velocity[0] << ", " << particles[1].velocity[1] << std::endl;
    //std::cout << "Particle 1 acceleration: " << particles[0].acceleration[0] << ", " << particles[0].acceleration[1] << std::endl;
    //std::cout << "Particle 2 acceleration: " << particles[1].acceleration[0] << ", " << particles[1].acceleration[1] << std::endl;
    //std::cout << "System momentum: " << system_momentum[0] << ", " << system_momentum[1] << std::endl;
    //std::cout << "System energy: " << system_energy << std::endl;
    //std::cout <<"calculation time:" << total_t <<  "seconds"<< std::endl;
    //std::cout << particleNum <<  "particles"<< std::endl;
    //std::cout << num_steps <<  "steps"<< std::endl;
    return 0;
}