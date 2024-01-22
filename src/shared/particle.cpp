#include "particle.h"

Particle::Particle() {
    mass = 0;
    pos = {0, 0, 0};
    vel = {0, 0, 0};
    force = {0, 0, 0};
}

Particle::Particle(double m, Eigen::Vector3d x, Eigen::Vector3d v, Eigen::Vector3d f) {
    mass = m;
    pos = x;
    vel = v;
    force = f;
}

double Particle::get_distance(Particle p) const {
    return (this->pos - p.pos).norm();
}

double Particle::get_distance(Eigen::Vector3d p) const {
    return (this->pos - p).norm();
}

double Particle::get_radius() const {
    return pos.norm();
}








