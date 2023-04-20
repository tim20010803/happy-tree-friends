#include <iostream>
#include <vector>
#include <cmath>
#include "matplotlibcpp.h"
#include <chrono>
#include <omp.h>

namespace plt = matplotlibcpp;

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

void plot_particles(std::vector<Particle>& particles) {
    std::vector<double> x, y;
    for (const Particle& p : particles) {
        x.push_back(p.x);
        y.push_back(p.y);
    }

    plt::plot(x, y, "ro");
    plt::title("Initial particle positions");
    plt::show();
}

void plot_final_positions(std::vector<Particle>& particles) {
    std::vector<double> x, y;
    for (const Particle& p : particles) {
        x.push_back(p.x);
        y.push_back(p.y);
    }

    plt::plot(x, y, "bo");
    plt::title("Final particle positions");
    plt::show();
}

int main() {
    const int num_particles = 1000;
    const int num_steps = 100;
    const double dt=0.1;
    std::vector<Particle> particles(num_particles);
for (int i = 0; i < num_particles; i++) {
    particles[i].mass = 1.0;
    particles[i].x = rand() % 100 + 1;
    particles[i].y = rand() % 100 + 1;
    particles[i].vx = 0;
    particles[i].vy = 0;
}

plot_particles(particles);

auto start_time = std::chrono::high_resolution_clock::now();

#pragma omp parallel for
for (int i = 0; i < num_steps; i++) {
    Quadtree quadtree(0, 0, 100, 100);

    for (Particle& p : particles) {
        quadtree.insert(&p);
    }

    for (Particle& p : particles) {
        p.vx = 0;
        p.vy = 0;

        quadtree.calculate_gravity(&p);

        p.x += p.vx*dt;
        p.y += p.vy*dt;
    }
}

auto end_time = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> duration = end_time - start_time;

std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;

plot_final_positions(particles);

return 0;
}
