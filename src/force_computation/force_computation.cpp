//
// Created by noahe on 10/01/2024.
//

#include <vector>
#include <fstream>
#include <iomanip>
#include "../shared/particle.h"
#include "../shared/data_reader.h"
#include "force_helper.h"

int main () {
    std::vector<Particle> particles = buildFromData();
//    particles = std::vector<Particle>(particles.begin(), particles.begin() + 1000);

    double max_radius = 2.5;
    int bins = 200;

    //calculate the observed force per bin directly
    double softening = mean_distance(particles);
    compute_forces(particles, softening);
    std::vector<double> observed_forces = forces_to_bins(particles, bins, max_radius);

    //calculate the approximated force per bin
    std::vector<double> approximated_forces = approximate_forces_in_bins(particles, bins, max_radius);

    //octree approximation
    Octree_force_approximation(particles, softening, 1.5);
    std::vector<double> octree_forces = forces_to_bins(particles, bins, max_radius);


    //write the data to a file
    double bin_size = max_radius / bins;
    std::ofstream myfile;
    myfile.open ("../output/direct_vs_approx.txt");
    for(int i = 0; i < bins; ++i){
            myfile << std::setw(6) << bin_size*i << "  "
                    << std::setw(11) << observed_forces[i] << "  "
                    << approximated_forces[i] << std::endl;
    }
    myfile.close();

    std::ofstream myfile2;
    myfile2.open ("../output/direct_vs_octree.txt");
    for(int i = 0; i < bins; ++i){
        myfile2 << std::setw(6) << bin_size*i << "  "
               << std::setw(11) << observed_forces[i] << "  "
               << octree_forces[i] << std::endl;
    }
    myfile2.close();

    double mean_relative_error;
    double mean_octree_error;
    for (int i = 0; i < bins; ++i) {
        mean_relative_error += std::abs(observed_forces[i] - approximated_forces[i])/observed_forces[i];
        mean_octree_error += std::abs(observed_forces[i] - octree_forces[i])/observed_forces[i];
    }
    mean_relative_error /= bins;
    std::cout << "mean relative error: " << mean_relative_error << std::endl;
    mean_octree_error /= bins;
    std::cout << "mean octree error: " << mean_octree_error << std::endl;


    return 0;
}