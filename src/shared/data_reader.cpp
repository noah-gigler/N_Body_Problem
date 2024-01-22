//
// Created by noahe on 12/01/2024.
//
#include "data_reader.h"

std::vector<Particle> buildFromData() {
    std::string path = "../data/data.txt";
    std::ifstream file(path);

    std::vector<Particle> particles;

    if (file.is_open()) {
        double index, m, softening, potential;
        Eigen::Vector3d x, v;
        //index, softening, and potential are just read but not used
        while (file >> index >> m >> x[0] >> x[1] >> x[2] >> v[0] >> v[1] >> v[2] >> softening >> potential) {
            Particle p(m, x, v);
            particles.emplace_back(p);
        }
    }
    else {
        std::cerr << "Unable to open file";
        return {};
    }

    return particles;
}

std::vector<double> getRadii(std::vector<Particle> particles) {
    std::vector<double> radii;
    for (const auto& particle : particles) {
        radii.push_back(particle.get_radius());
    }
    return radii;
}
