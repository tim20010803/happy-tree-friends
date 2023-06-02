#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>


void random_data( int particleNum, double r, double v_norm, double pos_range, double vec_range, std::string& str)
{
    std::vector<double> Mass;
    std::vector<double> PositionX;
    std::vector<double> PositionY;
    std::vector<double> VectorX;
    std::vector<double> VectorY;
    std::vector<double> AccelerationX;
    std::vector<double> AccelerationY;

    // produce random seed
    std::random_device rd;
    std::mt19937 gen(rd());
    // set random range
    std::uniform_real_distribution<double> massDist(1.0, 10.0);
    std::uniform_real_distribution<double> xDist(-1 * pos_range, pos_range); 
    std::uniform_real_distribution<double> yDist(-1 * pos_range, pos_range);
    std::uniform_real_distribution<double> wDist(0, 1);
    std::uniform_real_distribution<double> velocityXDist(-1 * vec_range, vec_range);
    std::uniform_real_distribution<double> velocityYDist(-1 * vec_range, vec_range);


    while ( PositionX.size() < particleNum )
    {   
        // set parameter  
        double m = massDist(gen);
        double x = xDist(gen);
        double y = yDist(gen);
        double w = wDist(gen);
        double delta_vx = velocityXDist(gen);
        double delta_vy = velocityYDist(gen);
        double d = sqrt(x*x + y*y);
        // let the particles in the range of a disk
        if ( d <=  r and w*pos_range*pos_range*pos_range >= d*d*d )
        {
            // compute the velocity and the acceleratiion
            double vx = - v_norm * d * d * y / d + delta_vx;
            double vy = v_norm * d * d * x / d + delta_vy;
            double ax = 0 ;
            double ay = 0;

            // add the particle initial condition
            Mass.push_back(m);
            PositionX.push_back(x);
            PositionY.push_back(y);
            VectorX.push_back(vx);
            VectorY.push_back(vy);
            AccelerationX.push_back(ax);
            AccelerationY.push_back(ay);
        }
    }

    std::ofstream file(str);
    file << "Time,Particle,Mass,PositionX,PositionY,VectorX,VectorY,AccelerationX,AccelerationY" << std::endl;
    for (int i = 0; i < particleNum; i++)
    {
        double t = 0;
        int p = i + 1;
        double mass = Mass[i];
        double posX = PositionX[i];
        double posY = PositionY[i];
        double vecX = VectorX[i];
        double vecY = VectorY[i];
        double accX = AccelerationX[i];
        double accY = AccelerationY[i];
        //std::cout << i << " * " << n << " + " << j << " = " << i*n + j << std::endl;
        file << t << "," << p  << "," << mass << "," 
                << posX << "," << posY << "," << vecX << "," << vecY << "," 
                << accX << "," << accY << std::endl;
    }
    file.close();
}


int main()
{
    int particleNum = 10; // number of particle
    double r = 1000; // radius of the initial disk
    double v_norm = 1e-3; // the coefficient of velocity
    double pos_range = 1000.; // the position range of the random number
    double vec_range = 100.; // the velocity range of the random number

    std::string str = "one_step_data.csv";
    
    random_data( particleNum, r, v_norm, pos_range, vec_range, str );
    
    return 0;
}