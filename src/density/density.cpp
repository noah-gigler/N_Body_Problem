#include <Eigen/Core>
#include <algorithm>
#include "../shared/particle.h"
#include "../shared/data_io.h"
#include "non_linear_solver.h"

int main() {
    int bins = 200;
    double max_radius = 1;
    double bin_size = max_radius/bins;

    //PRE: all particles have the same mass
    std::vector<Particle> particles = buildFromData();
    double particle_mass = particles[0].mass;
    double total_mass = 0;

    //find the number of particles in every bin
    std::vector<double> num_particles(bins, 0);
    for(auto i:particles){
        if(i.get_radius() > max_radius) continue;
        total_mass += particle_mass;
        int b = floor(i.get_radius()/bin_size);
        num_particles[b] += 1;
    }

    //multiply by mass and divide by the volume of the bin to obtain the density
    std::vector<double> density(bins, 0);
    for(int i = 0; i < bins; ++i){
        double inner_radius = bin_size*i;
        double outer_radius = bin_size*(i+1);
        double bin_volume = 4.0/3.0*M_PI*(pow(outer_radius,3)-pow(inner_radius,3));
        //measured density = num_particles * m/V
        density[i] = num_particles[i] * particle_mass/bin_volume;
    }

    //Fit the Hernquist profile to the measured data
    Eigen::VectorXd a(1);
    a(0) = 1; // initial guess for scale_length

    my_functor functor(density, bin_size, total_mass);
    Eigen::NumericalDiff<my_functor> numDiff(functor);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<my_functor>,double> lm(numDiff);
    lm.parameters.maxfev = 2000;
    lm.parameters.xtol = 1.0e-15;

    lm.minimize(a);
    std::cout << "number of iterations: " << lm.iter << std::endl;
    std::cout << "Optimal scale_length: " << a(0) << std::endl;

    double scale_length = a(0);
    std::vector<double> density_approximation(bins, 0);
    std::vector<double> num_expected_particles(bins, 0);
    for(int i = 0; i < bins; ++i){
        double radius = bin_size*i;
        //d(r) = M/2pi * a/r * 1/(a + r)^3
        density_approximation[i] = total_mass/(2*M_PI) * scale_length/radius *(1/pow(scale_length+radius,3));
        double bin_volume = 4.0/3.0*M_PI*(pow((i+1)*bin_size,3)-pow(i*bin_size,3));
        num_expected_particles[i] = density_approximation[i] * bin_volume/particle_mass;
    }

    std::vector<double> radius(bins);
    for(int i = 0; i < bins; ++i){
        radius[i] = bin_size*i;
    }

    write_to_file_2(radius, density, density_approximation, "density");
    write_to_file_2(radius, num_particles , num_expected_particles, "num_particles");

    return 0;
}