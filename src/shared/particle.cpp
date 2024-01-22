#include "particle.h"

Particle::Particle() {
    mass = 0;
    position = {0, 0, 0};
    velocity = {0, 0, 0};
    radius = 0;
    forces = {0, 0, 0};
}

Particle::Particle(double m, Eigen::Vector3d x, Eigen::Vector3d v, Eigen::Vector3d f) {
    mass = m;
    position = x;
    velocity = v;
    forces = f;
    radius = x.norm();
}

Particle::Particle(double x, double y, double z) {
    mass = 0;
    position = {x, y, z};
    velocity = {0, 0, 0};
    forces = {0, 0, 0};
    radius = position.norm();
}

bool Particle::operator==(const Particle& other) const {
    return (this->mass == other.mass && this->position == other.position&& this->velocity == other.velocity
    && this->forces == other.forces);
}

double Particle::get_radius() const {
    return radius;
}

double Particle::get_mass() const {
    return mass;
}

double Particle::get_x() const {
    return position[0];
}

double Particle::get_y() const {
    return position[1];
}

double Particle::get_z() const {
    return position[2];
}

double Particle::get_distance(Particle p) const {
    return (this->position - p.position).norm();
}

double Particle::get_distance(Eigen::Vector3d p) const {
    return (this->position - p).norm();
}

Eigen::Vector3d Particle::get_position() const {
    return position;
}

void Particle::set_forces(Eigen::Vector3d f) {
    forces = f;
}

Eigen::Vector3d Particle::get_forces() const {
    return forces;
}

Eigen::Vector3d Particle::add_forces(Eigen::Vector3d f) {
    #pragma omp atomic
    forces[0] += f[0];
    #pragma omp atomic
    forces[1] += f[1];
    #pragma omp atomic
    forces[2] += f[2];
    return forces;
}




