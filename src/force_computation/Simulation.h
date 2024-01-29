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

    Simulation(std::string name, std::vector<Particle> particles, double dt, double softening,
               force_type force, integrator_type integrator);
    void run(unsigned n_steps, unsigned steps_per_frame = 1);
    void leapfrog_step();
    void euler_step();
    void forces_step();
    void frame_to_file(int frame_number);

    std::vector<Particle> particles;
private:
    std::string name;
    force_type force;
    integrator_type integrator;
    double dt;
    double softening;
};


#endif //NBODYPROBLEM_SIMULATION_H
