#include <Eigen/Core>
#include <algorithm>
#include <gnuplot-iostream.h>
#include "../shared/particle.h"
#include "../shared/data_reader.h"
#include "non_linear_solver.h"

int main() {
    int bins = 200;
    double max_radius = 1;
    double bin_size = max_radius/bins;
    //double max_radius = *std::max_element(radii.begin(), radii.end());

    //PRE: all particles have the same mass
    std::vector<Particle> particles = buildFromData();
    double particle_mass = particles[0].get_mass();
    double total_mass = 0;

    //find the number of particles in every bin
    std::vector<double> radii = getRadii(particles);
    std::vector<int> num_particles(bins, 0);
    for(auto i:radii){
        if(i > max_radius) continue;
        total_mass += particle_mass;
        int b = floor(i/bin_size);
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
    for(int i = 0; i < bins; ++i){
        double radius = bin_size*i;
        //d(r) = M/2pi * a/r * 1/(a + r)^3
        density_approximation[i] = total_mass/(2*M_PI) * scale_length/radius *(1/pow(scale_length+radius,3));
    }


    // Plotting
    std::vector<std::pair<double, double>> data;
    std::vector<std::pair<double, double>> data_approximation;
    for (size_t i = 0; i < bins; ++i) {

        data.emplace_back(i, density[i]);
        data_approximation.emplace_back(i, density_approximation[i]);

    }

    Gnuplot gp;
    gp << "set terminal png\n"; // Set the output format to PNG
    gp << "set output 'output.png'\n"; // Set the output file
    //approximation is 0 at r = 0, so we plot starting at 1
    gp << "set logscale y\n"; // Set y axis to logarithmic scale
    gp << "set xrange [1:*]\nset yrange [0:*]\n";
    gp << "plot '-' with lines title 'density', '-' with lines title 'density approximation'\n";
    gp.send1d(data);
    gp.send1d(data_approximation);



    return 0;
}