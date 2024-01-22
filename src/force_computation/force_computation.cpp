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

//    //TODO: remove this line
//    particles = std::vector<Particle>(particles.begin(), particles.begin() + 20000);


    std::vector<double> radii = getR::adii(particles);
    double max_radius = 2.5;
    int bins = 200;
    double bin_size = max_radius / bins;



    //calculate the observed forces per bin directly
    double softening = mean_distance(particles);
    std::vector<double> forces = compute_forces(particles, softening);
    std::vector<double> observed_forces = forces_to_bins(forces, radii, bins, max_radius);

    //calculate the approximated forces per bin
    std::vector<double> approximated_forces = approximate_forces_in_bins(particles, bins, max_radius);

    //calculate the approximated forces per bin using the Octree
    //TODO: clear the particle vector
    //WARNING: this function changes the particles vector
//    std::vector<double> Octree_forces = Octree_force_approximation(particles, softening, 0.1);
//    std::vector<double> Octree_forces_in_bins = forces_to_bins(Octree_forces, radii, bins, max_radius);

    //print the forces (testing purposes)
    double observed_force = std::accumulate(approximated_forces.begin() + 1, approximated_forces.end(), 0.0);
    std::cout << "observed_forces: " <<  observed_force << std::endl;
    std::cout << "approximated_forces: " << approximated_forces[1] << std::endl;
//    std::cout << "Octree_forces: " << Octree_forces_in_bins[1] << std::endl;
//    std::cout << "Octree_forces_vector: " << particles[0].get_forces()[0] << " " << particles[0].get_forces()[1] << " " << particles[0].get_forces()[2] << std::endl;

    //write the data to a file
    std::ofstream myfile;
    myfile.open ("../data/forces.txt");
    for(int i = 0; i < bins; ++i){
        myfile << bin_size*i << " " << observed_forces[i] << " " << approximated_forces[i] << std::endl;
    }
    myfile.close();

//    //write Octree forces vector to a file
//    std::ofstream myfile2;
//    myfile2.open ("../data/Octree_forces.txt");
//    for(int i = 0; i < particles.size(); ++i){
//        myfile2 << particles[i].get_forces()[0] << " " << particles[i].get_forces()[1] << " " << particles[i].get_forces()[2] << std::endl;
//    }
//
//    //relative error between Octree and observed forces
//    std::ofstream myfile3;
//    myfile3.open ("../data/differences.txt");
//    for(int i = 0; i < particles.size(); ++i){
//        myfile3 << Octree_forces[i] - forces[i] << std::endl;
//    }

    return 0;
}