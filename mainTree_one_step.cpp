
#include <iostream>
#include <queue>
#include <vector>
#include <cmath>
#include <ctime>
#include <iomanip>
#include "orbit_integration_OMP.h"
#include "orbit_integration.h"
#include "quadrupleTree.h"
#include <omp.h>
#include <random>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
// main function is for testing
int main() {
    clock_t start_t, end_t;

    
    double t=3. ,dt = 0.005;
    // int num_steps = t / dt;
    // std::cout << "Physical Time: " << (num_steps)*dt <<"seconds"<< std::endl;
    // std::cout << particleNum <<  "particles"<< std::endl;
    // std::cout << num_steps <<  "steps"<< std::endl;

    // int index = 1;
    // for (const auto& particle : particles) {    
    //     std::cout << "Particle: " << index << ", Mass: " << particle.mass 
    //     << ", PositionX: " << particle.posi[0] << ", PositionY: " << particle.posi[1] 
    //     << ", VectorX: " << particle.velocity[0] << ", VectorY: " << particle.velocity[1] 
    //     << ", AccelerationX: " << particle.acceleration[0] << ", AccelerationY: " << particle.acceleration[1] << std::endl;
    //     index++;
    // }

    // build the tree
    

    // make a new file
    std::ofstream fileOut("RuntimeTreeComplex_uni2.csv");
    fileOut << "Type,NThread,theta,particleNumber,construct_time(s),force_time(s),total_time(s)" << std::endl;
    int pN=1000;
    for(int j=0;j<4;j++ )
    {
        // put the data into the fileOut
        for (int Thread = 1; Thread <= 8; Thread++){

            // input the data from the file
            std::string str = "one_step_data.csv";
            std::string str1 = "_uni_";
            std::string str2 = std::to_string(pN);
            std::string str3 = ".csv";
            str=str+str1+str2+str3;
            std::ifstream file(str);
            std::vector<Particle> particles; // store the particles
            std::string line; // each line of the file data
            bool isFirstLine = true; // if line is the first line or not
            int particleNum =0;
            int precision = std::min<int>(std::numeric_limits<double>::digits10 + 1,10);
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
                particles.push_back(particle);
                particleNum++;
            }
            file.close();
            double grid_range = 1000.; // the range of the Tree
            auto start = std::chrono::high_resolution_clock::now();

            QuadrupleTree T(1.,Thread,particles, -1 * grid_range, -1 * grid_range, -1 * grid_range, grid_range, grid_range, grid_range);


            auto mid = std::chrono::high_resolution_clock::now();
            T.TreeForce();
            T.~QuadrupleTree();

            auto end = std::chrono::high_resolution_clock::now();

            // Calculate duration
            std::chrono::duration<float> duration = end - start;
            std::chrono::duration<float> ForceDuration = end - mid;
            std::chrono::duration<float> ConsDuration = mid - start;
            float seconds = duration.count();
            float Fseconds = ForceDuration.count();
            float Cseconds = ConsDuration.count();

            fileOut<< std::setprecision(precision)<<"Uni"
            << ","<<T.NThread<< ","<<T.THETA<< ","<< particleNum<< ","<<Cseconds<< ","<<Fseconds<< ","<<seconds << std::endl;
            std::cout<< std::setprecision(precision)<<"Uni"
            << ","<<T.NThread<< ","<<T.THETA<< ","<< particleNum<< ","<<Cseconds<< ","<<Fseconds<< ","<<seconds << std::endl;
        }
        // put the data into the fileOut
        for (float i= 1.; i <= 6.; i++){

            std::string str = "one_step_data.csv";
            std::string str1 = "_uni_";
            std::string str2 = std::to_string(pN);
            std::string str3 = ".csv";
            str=str+str1+str2+str3;
            std::ifstream file(str);
            std::vector<Particle> particles; // store the particles
            std::string line; // each line of the file data
            bool isFirstLine = true; // if line is the first line or not
            int particleNum =0;
            int precision = std::min<int>(std::numeric_limits<double>::digits10 + 1,10);
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
                particles.push_back(particle);
                particleNum++;
            }
            file.close();
            double grid_range = 1000.; // the range of the Tree
            auto start = std::chrono::high_resolution_clock::now();

            QuadrupleTree T(0.2*i,8,particles, -1 * grid_range, -1 * grid_range, -1 * grid_range, grid_range, grid_range, grid_range);


            auto mid = std::chrono::high_resolution_clock::now();
            T.TreeForce();
            T.~QuadrupleTree();

            auto end = std::chrono::high_resolution_clock::now();

            // Calculate duration
            std::chrono::duration<float> duration = end - start;
            std::chrono::duration<float> ForceDuration = end - mid;
            std::chrono::duration<float> ConsDuration = mid - start;
            float seconds = duration.count();
            float Fseconds = ForceDuration.count();
            float Cseconds = ConsDuration.count();

            fileOut<< std::setprecision(precision)<<"Uni"
            << ","<<T.NThread<< ","<<T.THETA<< ","<< particleNum<< ","<<Cseconds<< ","<<Fseconds<< ","<<seconds << std::endl;
            std::cout<< std::setprecision(precision)<<"Uni"
            << ","<<T.NThread<< ","<<T.THETA<< ","<< particleNum<< ","<<Cseconds<< ","<<Fseconds<< ","<<seconds << std::endl;
        }
        pN*=10;
    }
    fileOut.close();

    // make a new file
    // std::ofstream fileOut("mainTree_data.csv");
    // fileOut << "Time,Particle,Mass,PositionX,PositionY,VectorX,VectorY,AccelerationX,AccelerationY,SystemEnergy,SystemMomentumX,SystemMomentumY,SystemAngularMomentumX,SystemAngularMomentumY" << std::endl;

    // // put the data into the fileOut
    // for (int step = 0; step <= num_steps; step++){
    //     for (int i = 0; i < particles.size(); i++)
    //     {
    //         int p = i+1;
    //         double mass = particles[i].mass;
    //         double posX = particles[i].posi[0];
    //         double posY = particles[i].posi[1];
    //         double vecX = particles[i].velocity[0];
    //         double vecY = particles[i].velocity[1];
    //         double accX = particles[i].acceleration[0];
    //         double accY = particles[i].acceleration[1];
    //         double system_energy = calculate_system_energy(particles);
    //         std::vector<double> system_momentum = calculate_system_momentum(particles);
    //         std::vector<double> system_angular_momentum = calculate_system_angular_momentum(particles);


    //         fileOut<< std::setprecision(precision)
    //         << std::scientific  << step * dt << "," << p  << "," << mass << "," 
    //         << posX << "," << posY << "," << vecX << "," << vecY << "," 
    //         << accX << "," << accY << "," << system_energy << "," 
    //         << system_momentum[0] << "," << system_momentum[1] << "," 
    //         << system_angular_momentum[0] << "," << system_angular_momentum[1] << std::endl;
    //     }
    //     Verlet_velocity_Tree(particles, dt, -1 * grid_range, -1 * grid_range, -1 * grid_range, grid_range, grid_range, grid_range);
    //     grid_range *= 1.0023; //graduate expand region
    // }
    // file.close();

    // std::vector<double> system_momentum = calculate_system_momentum(particles);
    // double system_energy = calculate_system_energy(particles);
    //std::cout << "Physical Time: " << (num_steps)*dt <<"seconds"<< std::endl;
    //std::cout << "System momentum: " << system_momentum[0] << ", " << system_momentum[1] << std::endl;
    // std::cout << "System energy: " << system_energy << std::endl;
    // std::cout <<"calculation time:" << seconds <<  "seconds"<< std::endl;
    // std::cout <<"calculationF time:" << Fseconds <<  "seconds"<< std::endl;
    // std::cout <<"calculationC time:" << Cseconds <<  "seconds"<< std::endl;
    // std::cout << particleNum <<  "particles"<< std::endl;
    
    // std::cout << num_steps <<  "steps"<< std::endl;
    return 0;
}
