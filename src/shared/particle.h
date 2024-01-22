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
    Particle(double x, double y, double z);
    bool operator==(const Particle& other) const;
    double get_radius() const;
    double get_mass() const;
    double get_x() const;
    double get_y() const;
    double get_z() const;
    Eigen::Vector3d get_position() const;
    double get_distance(Particle p) const;
    double get_distance(Eigen::Vector3d p) const;
    void set_forces(Eigen::Vector3d f);
    Eigen::Vector3d get_forces() const;
    Eigen::Vector3d add_forces(Eigen::Vector3d f);

private:
    double mass;
    double radius;
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
    Eigen::Vector3d forces;
};

#endif // NBODYPROBLEM_PARTICLE_H
