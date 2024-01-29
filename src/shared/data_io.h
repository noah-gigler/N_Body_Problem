#ifndef NBODYPROBLEM_DATA_IO_H
#define NBODYPROBLEM_DATA_IO_H

#include "particle.h"
#include <vector>

std::vector<Particle> buildFromData();

template<typename T, typename U>
void write_to_file(std::vector<T> index, std::vector<U> values, std::string name) {
    std::ofstream file;
    file.open("../output/" + name + ".txt");
    for(int i = 0; i < values.size(); ++i) {
        file << index[i] << " " << values[i] << std::endl;
    }
    file.close();
}

template<typename T, typename U>
void write_to_file_2(std::vector<T> index, std::vector<U> values, std::vector<U> values2, std::string name) {
    std::ofstream file;
    file.open("../output/" + name + ".txt");
    for(int i = 0; i < values.size(); ++i) {
        file << index[i] << " " << values[i] << " " << values2[i] << std::endl;
    }
    file.close();
}

#endif //NBODYPROBLEM_DATA_IO_H
