//
// Created by noahe on 25/01/2024.
//

#include "Simulation.h"
#include "../shared/data_io.h"
#include <chrono>

void test_conversion() {
    std::vector<Particle> particles = buildFromData();
    double steps = 2000;

}

void error_test(){
    std::vector<Particle> particles = buildFromData();
    // particles = std::vector<Particle>(particles.begin(), particles.begin() + 100);
    int steps = 3;
    double step_size = 1e-5;
    double softening = 0.2 * mean_distance(particles);
    Simulation sim_direct = Simulation("sim_basic", particles, step_size, softening,Simulation::Direct, Simulation::Leapfrog);
    sim_direct.run(steps);

//    Simulation sim_octree = Simulation("sim_octree", particles, step_size, softening,Simulation::Octree, Simulation::Leapfrog);
//    sim_octree.run(steps);
//
//    double mean_absolute_error = 0.0;
//    for(int i = 0; i < particles.size(); ++i) {
//        if((sim_direct.particles[i].pos - sim_octree.particles[i].pos).norm() > 0.001){
//            std::cout << sim_direct.particles[i].pos.transpose() << "\n " << sim_octree.particles[i].pos.transpose() << std::endl;
//        }
//        mean_absolute_error += (sim_direct.particles[i].pos - sim_octree.particles[i].pos).norm();
//    }
//    mean_absolute_error /= particles.size();
//    std::cout << "mean absolute error: " << mean_absolute_error << std::endl;
}

int main() {
    auto start = std::chrono::high_resolution_clock::now();
    error_test();
    auto stop = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    std::cout << time << std::endl;
}