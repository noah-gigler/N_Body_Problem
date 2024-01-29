//
// Created by noahe on 12/01/2024.
//
#include "data_io.h"

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


