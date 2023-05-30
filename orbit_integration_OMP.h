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
void calculate_gravity(std::vector<Particle>& particles, double G);
void Verlet_velocity(std::vector<Particle>& particles, double G, double dt);
void RK4(std::vector<Particle>& particles, double G, double dt);
void AB(std::vector<Particle>& particles, double G, double dt);
void Verlet_velocity_Tree(std::vector<Particle>& particles, double G, double dt,double mX,double mY,double mZ,double MX,double MY, double MZ);
void RK4_Tree(std::vector<Particle>& particles, double G, double dt,double mX,double mY,double mZ,double MX,double MY, double MZ);
void AB_Tree(std::vector<Particle>& particles, double G, double dt,double mX,double mY,double mZ,double MX,double MY, double MZ);

