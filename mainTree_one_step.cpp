
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
// main function is for testing
int main() {
    clock_t start_t, end_t;
    // double starttime =omp_get_wtime();

    // Particle a;
    // a.posi = {1.,6.,4.};
    // a.velocity = {4.,3.,2.};
    // a.mass = {12.};
    // a.acceleration = {0., 0., 0.};
    // Particle b;
    // b.posi = {2.,7.,8.};
    // b.velocity = {1.,6.,7.};
    // b.mass = {23.};
    // b.acceleration = {0., 0., 0.};
    // Particle c;
    // c.posi = {3.,5.,8.};
    // c.velocity = {1.,6.,7.};
    // c.mass = {212.};
    // c.acceleration = {0., 0., 0.};
    // Particle d;
    // d.posi = {4.,6.,8.};
    // d.velocity = {8.,6.,7.};
    // d.mass = {62.};
    // d.acceleration = {0., 0., 0.};
    // std::vector<Particle> particles = {a,b,c,d};
    // std::random_device rd;
    // std::mt19937_64 gen(rd());
    // std::uniform_real_distribution<> dist(-1000,1000);
    // const int particleNum = 100000;//11.653364000000secondsif100
    // for (int i = 0; i < particleNum-4; i++){
    //     Particle a;
    //     double rand3[3]={0};
    //     
    //     for (int j = 0; j < 3; j++){
    //         // rand3[j] = static_cast<double>(rand() % (int)(1e10))/1e10;
    //         rand3[j] = dist(gen);
    //     }
    //     a.posi = {rand3[0],rand3[1],rand3[2]};
    //     for (int j = 0; j < 3; j++){
    //         // rand3[j] = static_cast<double>(rand() % (int)(1e9))/1e9;
    //         rand3[j] = dist(gen);
    //     }
    //     a.velocity = {rand3[0],rand3[1],rand3[2]};
    //     a.mass = static_cast<double>(rand() % 1000)/10.;
    //     a.acceleration={0.,0.,0.};
    //     particles.push_back(a);
    // }
    // std::cout << std::fixed << std::setprecision(12);

    // set parameter
    int particleNum = 100; // number of particle
    double grid_range = 1000.; // the range of the Tree
    double r = grid_range/10.; // radius of the initial disk
    double v_norm = 1e-3; // the coefficient of velocity
    double pos_range = 1000.; // the position range of the random number
    double vec_range = 100.; // the velocity range of the random number
    double t=0.04*M_PI/pow( G_CONST*1000000000,0.5), dt=0.001; // Input time and time step
    int num_steps = t / dt;
    

    // double mX = -10000.;
    // double mY = -10000.;
    // double mZ = -10000.;
    // double MX = 10000.;
    // double MY = 10000.;
    // double MZ = 10000.;

    // std::vector<Particle> particles;
    
    // produce random seed
    // std::random_device rd;
    // std::mt19937 gen(rd());
    // set random range
    // std::uniform_real_distribution<double> massDist(1.0, 10.0);
    // std::uniform_real_distribution<double> xDist(-1 * pos_range, pos_range); 
    // std::uniform_real_distribution<double> yDist(-1 * pos_range, pos_range);
    // std::uniform_real_distribution<double> wDist(0, 1);
    // std::uniform_real_distribution<double> velocityXDist(-1 * vec_range, vec_range);
    // std::uniform_real_distribution<double> velocityYDist(-1 * vec_range, vec_range);

    // while ( particles.size() < particleNum ){   
    //     // set parameter 
    //     double m = massDist(gen);
    //     double x = xDist(gen);
    //     double y = yDist(gen);
    //     double w = wDist(gen);
    //     double delta_vx = velocityXDist(gen);
    //     double delta_vy = velocityYDist(gen);
    //     double d = sqrt(x*x + y*y);
    //     // let the particles in the range of a disk
    //     if ( d <=  r and  w*pos_range*pos_range*pos_range >= d*d*d ){
    //         Particle particle;
            
    //         // compute the velocity and the acceleratiion
    //         double vx = - v_norm * d * d * y / d + delta_vx;
    //         double vy = v_norm * d * d * x / d + delta_vy;
    //         double ax = 0;
    //         double ay = 0;

    //         // add the particle initial condition
    //         particle.mass = m;
    //         particle.posi.push_back(x);
    //         particle.posi.push_back(y);
    //         particle.velocity.push_back(vx);
    //         particle.velocity.push_back(vy);
    //         particle.acceleration.push_back(ax);
    //         particle.acceleration.push_back(ay);

    //         // add the particle
    //         particles.push_back(particle);
    //     }
    // }
    start_t = clock();

    // input the data from the file
    std::ifstream file("one_step_data.csv");
    std::vector<Particle> particles; // store the particles
    std::string line; // each line of the file data
    bool isFirstLine = true; // if line is the first line or not

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
    }

    int index = 1;
    for (const auto& particle : particles) {       
        std::cout << "Particle: " << index << ", Mass: " << particle.mass 
        << ", PositionX: " << particle.posi[0] << ", PositionY: " << particle.posi[1] 
        << ", VectorX: " << particle.velocity[0] << ", VectorY: " << particle.velocity[1] 
        << ", AccelerationX: " << particle.acceleration[0] << ", AccelerationY: " << particle.acceleration[1] << std::endl;
        index++;
    }

    // build the tree
    QuadrupleTree T(particles, -1 * grid_range, -1 * grid_range, -1 * grid_range, grid_range, grid_range, grid_range); 
    std::cout << "theta: " << T.THETA<<"\n";
    // for (int i = 0; i < particles.size(); i++)
    // {
    //     std::cout << "Particle mass: " << particles[i].mass << std::endl;
    //     std::cout << "Particle position: " << particles[i].posi[0] << ", " << particles[i].posi[1] << std::endl;
    //     std::cout << "Particle velocity: " << particles[i].velocity[0] << ", " << particles[i].velocity[1] << std::endl;
    //     std::cout << "Particle acceleration: " << particles[i].acceleration[0] << ", " << particles[i].acceleration[1] << std::endl;

    // }

    T.TreeForce();
    T.~QuadrupleTree();

    // for (int i = 0; i < particles.size(); i++)
    // {
    //     std::cout<<i<<"th" << "Particle mass: " << particles[i].mass << std::endl;
    //     std::cout<<i<<"th" << "Particle position: " << particles[i].posi[0] << ", " << particles[i].posi[1] << std::endl;
    //     std::cout<<i<<"th" << "Particle velocity: " << particles[i].velocity[0] << ", " << particles[i].velocity[1] << std::endl;
    //     std::cout<<i<<"th" << "Particle acceleration: " << particles[i].acceleration[0] << ", " << particles[i].acceleration[1] << std::endl;

    // }

    // if(num_steps != 0){
    //     for (int i = 0; i <= num_steps; i++) {
    //         Verlet_velocity_Tree(particles, dt,-10000.,-10000.,-10000.,10000.,10000.,10000.);
    //         std::cout<<"running at "<<i<<"-th step\n";
    //         // AB_Tree(particles, dt,-10000.,-10000.,-10000.,10000.,10000.,10000.);
    //         // RK4_Tree(particles, dt,-10000.,-10000.,-10000.,10000.,10000.,10000.);
    //     }
    // }


    // make a new file
    // std::ofstream file("mainTree_data.csv");
    // file << "Time,Particle,Mass,PositionX,PositionY,VectorX,VectorY,AccelerationX,AccelerationY,SystemEnergy,SystemMomentumX,SystemMomentumY,SystemAngularMomentumX,SystemAngularMomentumY" << std::endl;

    // put the data into the file
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


    //         file << time << "," << p  << "," << mass << "," 
    //         << posX << "," << posY << "," << vecX << "," << vecY << "," 
    //         << accX << "," << accY << "," << system_energy << "," 
    //         << system_momentum[0] << "," << system_momentum[1] << "," 
    //         << system_angular_momentum[0] << "," << system_angular_momentum[1] << std::endl;
    //     }
    //     Verlet_velocity_Tree(particles, dt, -1 * grid_range, -1 * grid_range, -1 * grid_range, grid_range, grid_range, grid_range);
    // }
    // file.close();
    
    end_t = clock();
    // std::vector<double> system_momentum = calculate_system_momentum(particles);
    // double system_energy = calculate_system_energy(particles);
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
