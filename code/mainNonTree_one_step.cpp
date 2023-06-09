
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
#include <string>
#include <sstream>
#include <chrono>
// main function is for testing
int main() {
    clock_t start_t, end_t;
    // input the data from the file
    std::ifstream file("../initial_condition/one_step_data.csv");
    std::vector<Particle> particles; // store the particles
    std::string line; // each line of the file data
    bool isFirstLine = true; // if line is the first line or not
    int particleNum =0;
    // to get all of the data
    while (std::getline(file, line)) {
        // delete the first line because the first line is
        // "Time,Particle,Mass,PositionX,PositionY,VectorX,VectorY,AccelerationX,AccelerationY"
        if (isFirstLine) {
            isFirstLine = false;
            continue; 
        }

        // set the parameter to store data
        std::istringstream iss(line);
        std::string element;
        std::vector<std::string> elements; 
        // change the line string into a vector
        while (std::getline(iss, element, ',')) {
            elements.push_back(element);
        }

        Particle particle;
        // input data
        double m = std::stod(elements[2]);
        double x = std::stod(elements[3]);
        double y = std::stod(elements[4]);
        double vx = std::stod(elements[5]);
        double vy = std::stod(elements[6]);
        double ax = std::stod(elements[7]);
        double ay = std::stod(elements[8]);
        // add the particle initial condition
        particle.mass = m;
        particle.posi.push_back(x);
        particle.posi.push_back(y);
        particle.velocity.push_back(vx);
        particle.velocity.push_back(vy);
        particle.acceleration.push_back(ax);
        particle.acceleration.push_back(ay);

        // add the particle
        particleNum++;
        particles.push_back(particle);
    }

    // int index = 1;
    // for (const auto& particle : particles) {        
    //     std::cout << "Particle: " << index << ", Mass: " << particle.mass 
    //     << ", PositionX: " << particle.posi[0] << ", PositionY: " << particle.posi[1] 
    //     << ", VectorX: " << particle.velocity[0] << ", VectorY: " << particle.velocity[1] 
    //     << ", AccelerationX: " << particle.acceleration[0] << ", AccelerationY: " << particle.acceleration[1] << std::endl;
    //     index++;
    // }
    double t=3. ,dt = 0.005;
    // int num_steps = t / dt;
    int num_steps =t/dt;
    std::cout << "Physical Time: " << (num_steps)*dt <<"seconds"<< std::endl;
    std::cout << particleNum <<  "particles"<< std::endl;
    std::cout << num_steps <<  "steps"<< std::endl;

    // start_t = clock();
    auto start = std::chrono::high_resolution_clock::now();
    calculate_gravity(particles);
    // make a new file
    std::ofstream fileOut("mainNonTree_data.csv");
    fileOut << "Time,Particle,Mass,PositionX,PositionY,VectorX,VectorY,AccelerationX,AccelerationY,SystemEnergy,SystemMomentumX,SystemMomentumY,SystemAngularMomentumX,SystemAngularMomentumY" << std::endl;
    int precision = std::min<int>(std::numeric_limits<double>::digits10 + 1,10);
    // put the data into the file
    for (int step = 0; step <= num_steps; step++){
        for (int i = 0; i < particles.size(); i++)
        {
            int p = i+1;
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


            fileOut  << std::setprecision(precision)
                << std::scientific << step * dt  << "," << p  << "," << mass << "," 
                << posX << "," << posY << "," << vecX << "," << vecY << "," 
                << accX << "," << accY << "," << system_energy << "," 
                << system_momentum[0] << "," << system_momentum[1] << "," 
                << system_angular_momentum[0] << "," << system_angular_momentum[1] << std::endl;
        }
        Verlet_velocity(particles, dt);
    }
    fileOut.close();
    // End timer
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate duration
    std::chrono::duration<float> duration = end - start;
    float seconds = duration.count();
    

    // for (int i = 0; i < particles.size(); i++)
    // {
    //     std::cout << "Particle mass: " << particles[i].mass << std::endl;
    //     std::cout << "Particle position: " << particles[i].posi[0] << ", " << particles[i].posi[1] << std::endl;
    //     std::cout << "Particle velocity: " << particles[i].velocity[0] << ", " << particles[i].velocity[1] << std::endl;
    //     std::cout << "Particle acceleration: " << particles[i].acceleration[0] << ", " << particles[i].acceleration[1] << std::endl;

    // }


    // end_t = clock();
    // std::vector<double> system_momentum = calculate_system_momentum(particles);
    // double system_energy = calculate_system_energy(particles);
    // double total_t = static_cast<double>(end_t - start_t) / CLOCKS_PER_SEC;
    //std::cout << "Physical Time: " << (num_steps)*dt <<"seconds"<< std::endl;
    //std::cout << "System momentum: " << system_momentum[0] << ", " << system_momentum[1] << std::endl;
    //std::cout << "System energy: " << system_energy << std::endl;
    std::cout <<"calculation time:" << seconds <<  "seconds"<< std::endl;
    std::cout << particleNum <<  "particles"<< std::endl;
    std::cout << num_steps <<  "steps"<< std::endl;
    return 0;
}
