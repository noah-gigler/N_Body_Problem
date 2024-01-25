//
// Created by noahe on 08/01/2024.
//

#include "Octree.h"

Octree::Octree(std::vector<Particle> &particles) {
    set_limits(particles);
    init_tree(particles);
}

Octree::Octree(std::vector<Particle> &particles, Eigen::Vector3d minimum, Eigen::Vector3d maximum) {
    min = minimum;
    max = maximum;
    init_tree(particles);
}

void Octree::init_tree(std::vector<Particle> &particles) {
    //check if particles are in box
    std::vector<Particle> particles_in_box;
    for (const Particle& p: particles) {
        if (is_in_box(p)) {
            particles_in_box.push_back(p);
        }
    }
    if (particles_in_box.size() > 1) {
        type = Internal;
        init_children(particles_in_box, min, max);
        mass = 0;
        com = {0, 0, 0};
        for(auto &p : particles_in_box){
            mass += p.mass;
            com += p.pos * p.mass;
        }
        com /= mass;
    }
    else{
        if(particles_in_box.size() == 1) {
            particle_in_box = particles_in_box[0];
            type = Leaf;
        }
        else{
            type = Empty;
        }
    }
}

void Octree::set_limits(std::vector<Particle> particles) {
    std::sort(particles.begin(), particles.end(), [](Particle a, Particle b) {
        return a.pos.x() < b.pos.x();
    });
    double x_min = particles[0].pos.x();
    double x_max = particles[particles.size() - 1].pos.x();

    std::sort(particles.begin(), particles.end(), [](Particle a, Particle b) {
        return a.pos.y() < b.pos.y();
    });
    double y_min = particles[0].pos.y();
    double y_max = particles[particles.size() - 1].pos.y();

    std::sort(particles.begin(), particles.end(), [](Particle a, Particle b) {
        return a.pos.z() < b.pos.z();
    });
    double z_min = particles[0].pos.z();
    double z_max = particles[particles.size() - 1].pos.z();

    max = Eigen::Vector3d(x_max, y_max, z_max);
    min = Eigen::Vector3d(x_min, y_min, z_min);
}

void Octree::init_children(std::vector<Particle> &particles, Eigen::Vector3d minimum, Eigen::Vector3d maximum) {
    Eigen::Vector3d center = this->get_center();
    children[0] = std::make_unique<Octree>(particles, Eigen::Vector3d(minimum.x(), minimum.y(), minimum.z()), Eigen::Vector3d(center.x(), center.y(), center.z()));
    children[1] = std::make_unique<Octree>(particles, Eigen::Vector3d(center.x(), minimum.y(), minimum.z()), Eigen::Vector3d(maximum.x(), center.y(), center.z()));
    children[2] = std::make_unique<Octree>(particles, Eigen::Vector3d(minimum.x(), center.y(), minimum.z()), Eigen::Vector3d(center.x(), maximum.y(), center.z()));
    children[3] = std::make_unique<Octree>(particles, Eigen::Vector3d(center.x(), center.y(), minimum.z()), Eigen::Vector3d(maximum.x(), maximum.y(), center.z()));
    children[4] = std::make_unique<Octree>(particles, Eigen::Vector3d(minimum.x(), minimum.y(), center.z()), Eigen::Vector3d(center.x(), center.y(), maximum.z()));
    children[5] = std::make_unique<Octree>(particles, Eigen::Vector3d(center.x(), minimum.y(), center.z()), Eigen::Vector3d(maximum.x(), center.y(), maximum.z()));
    children[6] = std::make_unique<Octree>(particles, Eigen::Vector3d(minimum.x(), center.y(), center.z()), Eigen::Vector3d(center.x(), maximum.y(), maximum.z()));
    children[7] = std::make_unique<Octree>(particles, Eigen::Vector3d(center.x(), center.y(), center.z()), Eigen::Vector3d(maximum.x(), maximum.y(), maximum.z()));
}

Eigen::Vector3d Octree::get_center() {
    return {(max.x() + min.x()) / 2,
            (max.y() + min.y()) / 2,
            (max.z() + min.z()) / 2};
}

bool Octree::is_in_box(const Particle &p) {
    return (p.pos.x() >= min.x() && p.pos.x() <= max.x()
            && p.pos.y() >= min.y() && p.pos.y() <= max.y()
            && p.pos.z() >= min.z()&& p.pos.z() <= max.z());
}

double Octree::side_length() {
    Eigen::Vector3d vec = max - min;
    return vec.lpNorm<Eigen::Infinity>();
}

Eigen::Vector3d Octree::force_on_particle(Particle &p, double softening, double theta) {
    Eigen::Vector3d force = {0, 0, 0};
    if(type == Empty) {
        return force;
    }
    if (type == Leaf) {
        //TODO make sure this works
        if(particle_in_box.pos == p.pos){
            return force;
        }
        else{
            return gravitational_force(p, particle_in_box, softening);
        }
    }
    //check if particle is far enough away (D/r < theta)
    double a = side_length()/p.get_distance(com);
    if (a < theta) {
        Particle p2 = Particle(mass, com);
        force += gravitational_force(p, p2, softening);
    }
    else {
        for (const auto & i : children) {
            force += i->force_on_particle(p, softening, theta);
        }
    }
    return force;
}


