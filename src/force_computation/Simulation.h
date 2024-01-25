//
// Created by noahe on 25/01/2024.
//

#ifndef NBODYPROBLEM_SIMULATION_H
#define NBODYPROBLEM_SIMULATION_H

#include <vector>
#include "../shared/particle.h"
#include "force_helper.h"
#include "Octree.h"
#include <iostream>


class Simulation {
public:
    enum force_type{
        Direct,
        Octree,
    };
    enum integrator_type{
        Euler,
        Leapfrog,
    };

    Simulation(std::vector<Particle> particles, double dt, double softening, force_type force, integrator_type integrator);
    void run(unsigned n_steps);
    void leapfrog_step();
    void euler_step();
    void forces_step();

private:
    force_type force;
    integrator_type integrator;
    std::vector<Particle> particles;
    double dt;
    double softening;
};


#endif //NBODYPROBLEM_SIMULATION_H
