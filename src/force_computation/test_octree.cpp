//
// Created by noahe on 15/01/2024.
//
#include "../shared/particle.h"
#include "../shared/data_reader.h"
#include "Octree.h"
#include "force_helper.h"
#include <iostream>
#include <chrono>

void test_theta(std::vector<Particle> &particles, unsigned bins, double max_radius, double softening, double N){
    double time_direct = 0.0;
    for(int i = 0; i < N; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        compute_forces(particles, softening);
        auto stop = std::chrono::high_resolution_clock::now();
        time_direct += std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    }
    time_direct /= N;
    std::vector<double> observed_forces = forces_to_bins(particles, bins, max_radius);

    std::cout << "direct: " << time_direct << std::endl;

    for(auto i: {0.5, 1.0, 1.5}) {
        double time_octree = 0.0;
        for(int n = 0; n < N; ++n) {
            auto start = std::chrono::high_resolution_clock::now();
            Octree_force_approximation(particles, softening, i);
            auto stop = std::chrono::high_resolution_clock::now();
            time_octree = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
        }
        time_octree /= N;
        double mean_relative_error = 0.0;
        std::vector<double> octree_forces = forces_to_bins(particles, bins, max_radius);
        for (int j = 0; j < bins; ++j) {
            mean_relative_error += std::abs(observed_forces[j] - octree_forces[j])/observed_forces[j];
        }
        mean_relative_error /= bins;
        std::cout << "theta: " << i << " time: " << time_octree << " error: " << mean_relative_error << std::endl;
    }

}

int main() {
    std::vector<Particle> particles = buildFromData();
    unsigned bins = 200;
    double max_radius = 2.5;
    Octree tree = Octree(particles);
    double softening = mean_distance(particles);

    Octree_force_approximation(particles, softening, 0);

//    test_theta(particles, bins, max_radius, 1 * softening, 1);


}