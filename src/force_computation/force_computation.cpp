//
// Created by noahe on 10/01/2024.
//

#include <vector>
#include <fstream>
#include <numeric>
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

    //write the data to a file
    double bin_size = max_radius / bins;
    std::ofstream myfile;
    myfile.open ("../data/force.txt");
    for(int i = 0; i < bins; ++i){
        myfile << bin_size*i << " " << observed_forces[i] << " " << approximated_forces[i] << std::endl;
    }
    myfile.close();

    return 0;
}