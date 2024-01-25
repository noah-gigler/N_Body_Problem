//
// Created by noahe on 25/01/2024.
//

#include "Simulation.h"
#include "../shared/data_reader.h"

int main() {
    std::vector<Particle> particles = buildFromData();
    Simulation sim = Simulation(particles, 0.01, 0.13);
    sim.run(100);
}