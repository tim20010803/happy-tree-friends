
#include <iostream>
#include <queue>
#include <vector>
#include <cmath>
#include <ctime>
#include <random>
#include <iomanip>
#include "orbit_integration_OMP.h"
#include "orbit_integration.h"
#include "quadrupleTree.h"
#include <omp.h>
// main function is for testing
int main() {
    clock_t start_t, end_t;
    // double starttime =omp_get_wtime();
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
    std::vector<Particle> particles = {a,b,c,d};
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<> dist(-1000,1000);
    const int particleNum = 100000;//11.653364000000secondsif100
    for (int i = 0; i < particleNum-4; i++){
        Particle a;
        double rand3[3]={0};
        
        for (int j = 0; j < 3; j++){
            // rand3[j] = static_cast<double>(rand() % (int)(1e10))/1e10;
            rand3[j] = dist(gen);
        }
        a.posi = {rand3[0],rand3[1],rand3[2]};
        for (int j = 0; j < 3; j++){
            // rand3[j] = static_cast<double>(rand() % (int)(1e9))/1e9;
            rand3[j] = dist(gen);
        }
        a.velocity = {rand3[0],rand3[1],rand3[2]};
        a.mass = static_cast<double>(rand() % 1000)/10.;
        a.acceleration={0.,0.,0.};
        particles.push_back(a);
    }
    std::cout << std::fixed << std::setprecision(12);
    start_t = clock();
    QuadrupleTree T(0.8,particles,-10000.,-10000.,-10000.,10000.,10000.,10000.); 
    std::cout << "theta: " << T.THETA<<"\n";
    // for (int i = 0; i < particles.size(); i++)
    // {
    //     std::cout << "Particle mass: " << particles[i].mass << std::endl;
    //     std::cout << "Particle position: " << particles[i].posi[0] << ", " << particles[i].posi[1] << std::endl;
    //     std::cout << "Particle velocity: " << particles[i].velocity[0] << ", " << particles[i].velocity[1] << std::endl;
    //     std::cout << "Particle acceleration: " << particles[i].acceleration[0] << ", " << particles[i].acceleration[1] << std::endl;

    // }
    
    // Input time and time step
    double t=0.04*M_PI/pow( G_CONST*1000000000,0.5), dt=0.001;

    // Perform simulation
    // int num_steps = t / dt;
    int num_steps =0;

    T.TreeForce();

    T.~QuadrupleTree();
    

    // for (int i = 0; i < particles.size(); i++)
    // {
    //     std::cout<<i<<"th" << "Particle mass: " << particles[i].mass << std::endl;
    //     std::cout<<i<<"th" << "Particle position: " << particles[i].posi[0] << ", " << particles[i].posi[1] << std::endl;
    //     std::cout<<i<<"th" << "Particle velocity: " << particles[i].velocity[0] << ", " << particles[i].velocity[1] << std::endl;
    //     std::cout<<i<<"th" << "Particle acceleration: " << particles[i].acceleration[0] << ", " << particles[i].acceleration[1] << std::endl;

    // }
    if(num_steps != 0){
        for (int i = 0; i <= num_steps; i++) {
            Verlet_velocity_Tree(particles, dt,-10000.,-10000.,-10000.,10000.,10000.,10000.);
            std::cout<<"running at "<<i<<"-th step\n";
            // AB_Tree(particles, dt,-10000.,-10000.,-10000.,10000.,10000.,10000.);
            // RK4_Tree(particles, dt,-10000.,-10000.,-10000.,10000.,10000.,10000.);
        }
    }
    
    end_t = clock();

    // std::vector<double> system_momentum = calculate_system_momentum(particles);
    // double system_energy = calculate_system_energy(particles);
    // double endtime =omp_get_wtime();
    double total_t = static_cast<double>(end_t - start_t) / CLOCKS_PER_SEC;
    std::cout << "Physical Time: " << (num_steps)*dt <<"seconds"<< std::endl;
    std::cout << "Particle 1 mass: " << particles[0].mass << std::endl;
    std::cout << "Particle 2 mass: " << particles[1].mass << std::endl;
    std::cout << "Particle 1 position: " << particles[0].posi[0] << ", " << particles[0].posi[1] << std::endl;
    std::cout << "Particle 2 position: " << particles[1].posi[0] << ", " << particles[1].posi[1] << std::endl;
    std::cout << "Particle 1 velocity: " << particles[0].velocity[0] << ", " << particles[0].velocity[1] << std::endl;
    std::cout << "Particle 2 velocity: " << particles[1].velocity[0] << ", " << particles[1].velocity[1] << std::endl;
    std::cout << "Particle 1 acceleration: " << particles[0].acceleration[0] << ", " << particles[0].acceleration[1] << std::endl;
    std::cout << "Particle 2 acceleration: " << particles[1].acceleration[0] << ", " << particles[1].acceleration[1] << std::endl;
    // std::cout << "System momentum: " << system_momentum[0] << ", " << system_momentum[1] << std::endl;
    // std::cout << "System energy: " << system_energy << std::endl;
    std::cout <<"calculation time:" << total_t <<  "seconds"<< std::endl;
    std::cout << particleNum <<  "particles"<< std::endl;
    std::cout << num_steps <<  "steps"<< std::endl;
    return 0;
}
