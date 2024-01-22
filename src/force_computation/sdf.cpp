//
// Created by noahe on 16/01/2024.
//
#include <vector>
#include "../shared/particle.h"
#include "../shared/data_reader.h"

std::vector<double> approximate_forces_in_bins(std::vector<Particle> particles, int bins, double max_radius){
    std::vector<double> radii = getRadii(particles);
    std::vector<double> mass(bins ,0);
    double bin_size = max_radius/bins;
    for(int i = 0; i < bins; ++i){
        if(radii[i] > max_radius) continue;
        int index = floor(radii[i]/bin_size);
        mass[index] += particles[i].mass;
    }
    for(int i = 1; i < mass.size(); ++i){
        mass[i] += mass[i-1];
    }
    std::vector<double> force(bins);
    for(int i = 0; i < force.size(); ++i){
        radius = bin_size
        force[i] = mass[i]/
    }
}
