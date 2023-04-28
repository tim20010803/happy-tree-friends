#include <iostream>
#include <vector>
#include <cmath>

struct Particle {
    double mass;
    double x, y;
    double vx, vy;
    double ax, ay;  // add acceleration components
    double ax_new, ay_new;
};

struct Quadtree {
    double x, y, width, height;
    Particle* particle;
    Quadtree* nw;
    Quadtree* ne;
    Quadtree* sw;
    Quadtree* se;

    Quadtree(double x_, double y_, double width_, double height_)
        : x(x_), y(y_), width(width_), height(height_), particle(NULL), nw(NULL), ne(NULL), sw(NULL), se(NULL) {}

    ~Quadtree() {
        if (nw != NULL) delete nw;
        if (ne != NULL) delete ne;
        if (sw != NULL) delete sw;
        if (se != NULL) delete se;
    }

    bool contains(Particle* p) {
        return (p->x >= x && p->x <= x+width && p->y >= y && p->y <= y+height);
    }

    void insert(Particle* p) {
        if (!contains(p)) return;

        if (particle == NULL) {
            particle = p;
        } else {
            if (nw == NULL) subdivide();
            nw->insert(p);
            ne->insert(p);
            sw->insert(p);
            se->insert(p);
        }
    }

    void subdivide() {
        double sub_width = width/2.0;
        double sub_height = height/2.0;
        nw = new Quadtree(x, y, sub_width, sub_height);
        ne = new Quadtree(x+sub_width, y, sub_width, sub_height);
        sw = new Quadtree(x, y+sub_height, sub_width, sub_height);
        se = new Quadtree(x+sub_width, y+sub_height, sub_width, sub_height);
    }

    void calculate_gravity(Particle* p) {
        if (particle != NULL && particle != p) {
            double dx = particle->x - p->x;
            double dy = particle->y - p->y;
            double r_squared = dx*dx + dy*dy;
            double r = sqrt(r_squared);
            double f = p->mass * particle->mass / r_squared;
            p->ax += f*dx/p->mass;
            p->ay += f*dy/p->mass;
        } else {
            if (nw != NULL) {
                nw->calculate_gravity(p);
                ne->calculate_gravity(p);
                sw->calculate_gravity(p);
                se->calculate_gravity(p);
            }
        }
    }

    void remove(Particle* p) {
        if (!contains(p)) return;
        // If the particle is in a leaf node, remove it
        if (particle == p) {
            particle = NULL;
            return;
        }

        // If the particle is not in a leaf node, recurse into the child nodes
        if (nw != NULL) {
            nw->remove(p);
            ne->remove(p);
            sw->remove(p);
            se->remove(p);

            // If all the child nodes are empty, delete them and set pointers to NULL
            if (nw->particle == NULL && ne->particle == NULL && sw->particle == NULL && se->particle == NULL) {
                delete nw;
                delete ne;
                delete sw;
                delete se;
                nw = NULL;
                ne = NULL;
                sw = NULL;
                se = NULL;
            }
        }
    }
};

int main() {
    std::vector<Particle> particles = {{100, 100, 100, 0, 0, 0, 0, 0, 0}, {200, 900, 900, 0, 0, 0, 0, 0, 0}, {300, 300, 700, 0, 0, 0, 0, 0, 0}};

    // Add particles to vector
    Quadtree tree(0, 0, 1000, 1000);
    for (Particle& p : particles) {
        tree.insert(&p);
    }

    double dt = 0.1; // time step
    int n_steps = 100; // number of time steps to simulate

    // Simulate motion of particles
    for (int i = 0; i < n_steps; i++) {
        // Calculate gravitational forces
        for (Particle& p : particles) {
            tree.calculate_gravity(&p);
        }

        // Integrate motion using Verlet integration
        for (Particle& p : particles) {
            double x_new = p.x + p.vx*dt + 0.5*p.ax*dt*dt;
            double y_new = p.y + p.vy*dt + 0.5*p.ay*dt*dt;

            p.ax_new = 0; // reset new acceleration
            p.ay_new = 0;

            // update velocity using new acceleration
            p.vx += 0.5*dt*(p.ax + p.ax_new);
            p.vy += 0.5*dt*(p.ay + p.ay_new);

            p.x = x_new;
            p.y = y_new;
        }

        // Update quadtree
        for (Particle& p : particles) {
            tree.remove(&p);
            tree.insert(&p);
        }
    }

    // Output final positions and velocities
    for (Particle& p : particles) {
        std::cout << "Particle at (" << p.x << ", " << p.y << ") has velocity (" << p.vx << ", " << p.vy << ")" << std::endl;
    }

    return 0;
}
