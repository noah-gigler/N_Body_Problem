//
// Created by noahe on 12/01/2024.
//
#include "force_helper.h"

Eigen::Vector3d newtonian_force(Particle p1, Particle p2, double softening){
    Eigen::Vector3d direction = p1.get_position() - p2.get_position();
    double distance = direction.norm();
    double force = - p1.get_mass() * p2.get_mass() / (distance * distance + softening * softening);
    return force * direction.normalized();
}

std::vector<double> compute_forces(std::vector<Particle> &particles, double softening){
    //make sure particle forces are 0 before adding new forces
    for(int i = 0; i < particles.size(); ++i){
        particles[i].set_forces({0,0,0});
    }
    std::vector<double> forces(particles.size());
    #pragma omp parallel for default(none) shared(particles) firstprivate(softening) collapse(2) schedule(guided)
    for (int i = 0; i < particles.size(); ++i) {
        for(int j = i + 1; j < particles.size(); ++j) {
            Eigen::Vector3d force = newtonian_force(particles[i], particles[j], softening);
            particles[i].add_forces(force);
            particles[j].add_forces(-force);
        }
    }
    for(int i = 0; i < particles.size(); ++i){
        forces[i] = particles[i].get_forces().norm();
    }
    return forces;
}

double mean_distance(std::vector<Particle> particles){
    std::vector<double> radii = getRadii(particles);
    std::sort(radii.begin(), radii.end());
    //PRE: all particles have the same mass
    double half_mass_radius = radii[radii.size()/2];
//    std::cout << "half mass radius: " << half_mass_radius << std::endl;
    std::vector<Particle> particles_in_hmr;
    for(const Particle& p : particles){
        if(p.get_radius() < half_mass_radius){
            particles_in_hmr.push_back(p);
        }
    }
//    std::cout << "particles in hmr: " << particles_in_hmr.size() << std::endl;
    double mean_distance = 0;
    #pragma omp parallel for default(none) shared(mean_distance, std::cout) firstprivate(particles_in_hmr) collapse(2) schedule(guided)
    for(int i = 0; i < particles_in_hmr.size(); ++i){
        for(int j = i + 1; j < particles_in_hmr.size(); ++j){
            #pragma omp atomic
            mean_distance += particles_in_hmr[i].get_distance(particles_in_hmr[j]);
        }
//        if(i % 1000 == 0) std::cout << i << std::endl;
    }
    mean_distance /= particles_in_hmr.size()*(particles_in_hmr.size()-1)/2;
    std::cout << "mean distance: " << mean_distance << std::endl;
    return mean_distance;
}

std::vector<double> forces_to_bins(std::vector<double> forces, std::vector<double> radii, unsigned bins, double max_radius){
    std::vector<double> forces_in_bins(bins);
    double bin_size = max_radius/bins;
    for(int i = 0; i < forces.size(); ++i){
        if(radii[i] > max_radius) continue;
        int b = floor(radii[i]/bin_size);
        forces_in_bins[b] += abs(forces[i]);
    }
    return forces_in_bins;
}

std::vector<double> approximate_forces_in_bins(std::vector<Particle> particles, int bins, double max_radius){
    double bin_size = max_radius/bins;
    std::vector<double> radii = getRadii(particles);

    //cumulative sum of mass vector
    std::vector<double> mass(bins, 0);
    for(int i = 0; i < radii.size(); ++i){
        if(radii[i] > max_radius) continue;
        int b = floor(radii[i]/bin_size);
        mass[b] += particles[i].get_mass();
    }
    for(int i = 1; i < bins; ++i){
        mass[i] += mass[i-1];
    }

    //calculate the approximated forces per bin
    std::vector<double> approximated_forces(bins);
    for(int i = 0; i < bins; ++i){
        double radius = bin_size*i;
        //F(r) = M(r)/r^2
        approximated_forces[i] = mass[i] / (radius * radius);
    }
    return approximated_forces;
};

std::vector<double> Octree_force_approximation(std::vector<Particle>& particles, double softening, double theta){
    Octree tree(particles);
    std::vector<double> forces(particles.size());
    for(int i = 0; i < particles.size(); ++i){
        Eigen::Vector3d force = tree.force_on_particle(particles[i], softening, theta);
        particles[i].set_forces(force);
        forces[i] = force.norm();
    }
    return forces;
}