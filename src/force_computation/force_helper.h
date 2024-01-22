//
// Created by noahe on 10/01/2024.
//

#ifndef NBODYPROBLEM_FORCE_HELPER_H
#define NBODYPROBLEM_FORCE_HELPER_H

#include "../shared/particle.h"
#include "../shared/data_reader.h"
//#include "Octree.h"
#include <cmath>
#include <iostream>

Eigen::Vector3d gravitational_force(Particle p1, Particle p2, double softening);

void compute_forces(std::vector<Particle> &particles, double softening);

double mean_distance(std::vector<Particle> &particles);

std::vector<double> forces_to_bins(std::vector<Particle> &particles, unsigned bins, double max_radius);

std::vector<double> approximate_forces_in_bins(std::vector<Particle> &particles, int n_bins, double max_radius);

std::vector<double> Octree_force_approximation(std::vector<Particle> &particles, double softening, double theta = 0.5);

#endif //NBODYPROBLEM_FORCE_HELPER_H
