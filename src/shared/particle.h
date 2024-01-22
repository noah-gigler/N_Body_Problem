#ifndef NBODYPROBLEM_PARTICLE_H
#define NBODYPROBLEM_PARTICLE_H

#include <Eigen/Dense>
#include <string>
#include <iostream>
#include <fstream>

class Particle {
public:
    Particle();
    Particle(double m, Eigen::Vector3d x, Eigen::Vector3d v = {}, Eigen::Vector3d forces = {});
    double get_distance(Particle p) const;
    double get_distance(Eigen::Vector3d p) const;
    double get_radius() const;

public:
    double mass;
    Eigen::Vector3d force;
    Eigen::Vector3d pos;
    Eigen::Vector3d vel;
};

#endif // NBODYPROBLEM_PARTICLE_H
