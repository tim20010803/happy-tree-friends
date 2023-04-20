#include <iostream>
#include <vector>
#include <cmath>

struct Particle {
    double mass;
    double x, y;
    double vx, vy;
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
            p->vx += f*dx/r;
            p->vy += f*dy/r;
        } else {
            if (nw != NULL) {
                nw->calculate_gravity(p);
                ne->calculate_gravity(p);
                sw->calculate_gravity(p);
                se->calculate_gravity(p);
            }
        }
    }
};

int main() {
    std::vector<Particle> particles;
    // add particles to vector

    Quadtree tree(0, 0, 1000, 1000);
    for (Particle& p : particles) {
        tree.insert(&p);
    }

    // calculate gravity for each particle
    for (Particle& p : particles) {
        tree.calculate_gravity(&p);
    }

    // output results
    for (Particle& p : particles) {
        std::cout << "Particle at (" << p.x << ", " << p.y << ") has velocity (" << p.vx << ", " << p.vy << ")" << std::endl;
    }

    return 0;
}
