#ifndef NBODYPROBLEM_DATA_READER_H
#define NBODYPROBLEM_DATA_READER_H

#include "particle.h"
#include <vector>

std::vector<Particle> buildFromData();

std::vector<double> getRadii(std::vector<Particle> particles);

#endif //NBODYPROBLEM_DATA_READER_H
