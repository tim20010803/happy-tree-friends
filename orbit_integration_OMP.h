#ifndef G_CONST
#define G_CONST 6.67428e-11
#endif

#ifndef PARTICLE
#define PARTICLE
struct Particle{
    std::vector<double> posi;
    std::vector<double> velocity;
    std::vector<double> acceleration;

    double mass;
};
#endif
void calculate_gravity(std::vector<Particle>& particles);
void Verlet_velocity(std::vector<Particle>& particles, double dt);
void RK4(std::vector<Particle>& particles, double dt);
void AB(std::vector<Particle>& particles, double dt);
void Verlet_velocity_Tree(std::vector<Particle>& particles, double dt,double mX,double mY,double mZ,double MX,double MY, double MZ);
void RK4_Tree(std::vector<Particle>& particles, double dt,double mX,double mY,double mZ,double MX,double MY, double MZ);
void AB_Tree(std::vector<Particle>& particles, double dt,double mX,double mY,double mZ,double MX,double MY, double MZ);
std::vector<double> calculate_system_momentum(const std::vector<Particle>& particles);
std::vector<double> calculate_system_angular_momentum(const std::vector<Particle>& particles);
double calculate_system_energy(const std::vector<Particle>& particles);

