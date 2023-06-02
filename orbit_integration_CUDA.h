#ifndef G_CONST
#define G_CONST 6.67428e-11
#endif



#ifndef PARTICLE
#define PARTICLE
struct Particle{
    std::vector<float> posi;
    std::vector<float> velocity;
    std::vector<float> acceleration;

    float mass;
};
#endif
void calculate_gravity(std::vector<Particle>& particles, float G);
__global__ void Verlet_velocity_cuda_kernel(Particle* particles, float G, float dt, int num_particles);
void RK4(std::vector<Particle>& particles, float G, float dt);
void AB(std::vector<Particle>& particles, float G, float dt);


