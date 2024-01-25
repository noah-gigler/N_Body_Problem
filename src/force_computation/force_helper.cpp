//
// Created by noahe on 12/01/2024.
//
#include "force_helper.h"

Eigen::Vector3d gravitational_force(const Particle &p1, const Particle &p2, double softening){
    Eigen::Vector3d direction = p1.pos - p2.pos;
    assert(direction.norm() != 0);
    double distance = direction.norm();
    double force = - p1.mass * p2.mass / (distance * distance + softening * softening);
    return force * direction.normalized();
}

Eigen::Vector3d atomic_vec_add(Eigen::Vector3d &vec, Eigen::Vector3d add){
    #pragma omp atomic
    vec[0] += add[0];
    #pragma omp atomic
    vec[1] += add[1];
    #pragma omp atomic
    vec[2] += add[2];
    return vec;
}

void compute_forces(std::vector<Particle> &particles, double softening){
    //make sure particle force are 0 before adding new force
    for(int i = 0; i < particles.size(); ++i){
        particles[i].force = {0, 0, 0};
    }
    std::vector<double> forces(particles.size());
    #pragma omp parallel for default(none) shared(particles) firstprivate(softening) collapse(2) schedule(guided)
    for (int i = 0; i < particles.size(); ++i) {
        for(int j = i + 1; j < particles.size(); ++j) {
            Eigen::Vector3d force = gravitational_force(particles[i], particles[j], softening);
            atomic_vec_add(particles[i].force, force);
            atomic_vec_add(particles[j].force, -force);
        }
    }
}


double mean_distance(std::vector<Particle> &particles){
    std::sort(particles.begin(), particles.end(), [](Particle p1, Particle p2){return p1.get_radius() < p2.get_radius();});
    //PRE: all particles have the same mass
    unsigned hmr_index = particles.size()/2;
    double mean_distance = 0;
    for(int i = 0; i < hmr_index; ++i){
        for(int j = i + 1; j < hmr_index; ++j){
            mean_distance += particles[i].get_distance(particles[j]);
        }
    }
    mean_distance /= hmr_index*(hmr_index-1)/2;
    std::cout << "mean distance: " << mean_distance << std::endl;
    return mean_distance;
}

std::vector<double> forces_to_bins(std::vector<Particle> &particles, unsigned n_bins, double max_radius){
    std::vector<double> forces_in_bins(n_bins, 0);
    std::vector<double> bin_pop(n_bins, 0);
    double bin_size = max_radius / n_bins;
    for(const auto& p: particles){
        if(p.get_radius() > max_radius) continue;
        int b = floor(p.get_radius()/bin_size);
        forces_in_bins[b] += p.force.norm();
        bin_pop[b] += 1;
    }
    for (int i = 0; i < n_bins; ++i) {
        if(bin_pop[i] == 0) continue;
        forces_in_bins[i] /= bin_pop[i];
    }
    return forces_in_bins;
}

std::vector<double> approximate_forces_in_bins(std::vector<Particle> &particles, int n_bins, double max_radius){
    double bin_size = max_radius / n_bins;

    //cumulative sum of mass vector
    std::vector<double> mass(n_bins, 0);
    for(int i = 0; i < particles.size(); ++i){
        if(particles[i].get_radius() > max_radius) continue;
        int b = floor(particles[i].get_radius()/bin_size);
        mass[b] += particles[i].mass;
    }
    //normalize the mass vector
    for(int i = 1; i < n_bins; ++i){
        mass[i] += mass[i-1];
    }

    //calculate the approximated force per bin
    std::vector<double> approximated_forces(n_bins);
    for(int i = 0; i < n_bins; ++i){
        double radius = bin_size * (i+1);
        //PRE: all particles have the same mass
        //F(r) = M(r)/r^2 * particle_mass
        approximated_forces[i] = mass[i] / (radius * radius) * particles[0].mass;
    }
    return approximated_forces;
};

void Octree_force_approximation(std::vector<Particle> &particles, double softening, double theta){
    //make sure particle force are 0 before adding new force
    for(int i = 0; i < particles.size(); ++i){
        particles[i].force = {0, 0, 0};
    }
    Octree tree(particles);
    std::vector<double> forces(particles.size());
    #pragma omp parallel for default(none) shared(particles, tree, std::cout) firstprivate(softening, theta) schedule(guided)
    for(int i = 0; i < particles.size(); ++i){
        Eigen::Vector3d force = tree.force_on_particle(particles[i], softening, theta);
        atomic_vec_add(particles[i].force, force);
        if(i % 1000 == 0){
            std::cout << "particle " << i << " done" << std::endl;
        }
    }
}