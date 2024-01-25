//
// Created by noahe on 25/01/2024.
//

#include "Simulation.h"

Simulation::Simulation(std::vector<Particle> particles, double dt, double softening, force_type force, integrator_type integrator) {
    this->particles = particles;
    this->dt = dt;
    this->softening = softening;
    this->force = force;
    this->integrator = integrator;
}

void Simulation::run(unsigned n_steps) {
    for(int i = 0; i < n_steps; ++i) {
        if(integrator == Euler) {
            euler_step();
        } else if(integrator == Leapfrog) {
            leapfrog_step();
        }
        std::cout << "step: " << i << std::endl;
    }
}

void Simulation::leapfrog_step() {
    //x_i+1 = x_i + v_i+1/2 * dt
    std::vector<Eigen::Vector3d> v_half;
    for(int i = 0; i < particles.size(); ++i) {
        v_half[i] = particles[i].vel + 0.5 * particles[i].force / particles[i].mass * dt;
        particles[i].pos += v_half[i] * dt;
    }
    forces_step();
    //v_i+1 = v_i+1/2 + 1/2 * a_i+1 * dt
    for(int i = 0; i < particles.size(); ++i) {
        particles[i].vel = v_half[i] + 0.5 * particles[i].force / particles[i].mass * dt;
    }
}

void Simulation::euler_step() {
    forces_step();
    for(auto &p: particles) {
        p.pos += p.vel * dt;
        p.vel += p.force / p.mass * dt;
    }
}

void Simulation::forces_step() {
    if(force == Direct) {
        compute_forces(particles, softening);
    } else if(force == Octree) {
        Octree_force_approximation(particles, softening);
    }
}

