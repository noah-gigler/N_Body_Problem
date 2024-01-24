//
// Created by noahe on 08/01/2024.
//

#include "Octree.h"

Octree::Octree(std::vector<Particle> &particles) {
    set_limits(particles);
    if (particles.size() > 1) {
        type = Internal;
        init_children(particles, min, max);
        mass = 0;
        com = {0, 0, 0};
        for(auto &p : particles){
            mass += p.mass;
            com += p.pos * p.mass;
        }
        com /= mass;
    }
    else{
        if(particles.size() == 1) {
            particle_in_box = particles[0];
            type = Leaf;
        }
        else{
            type = Empty;
        }
    }
}

Octree::Octree(std::vector<Particle> &particles, Eigen::Vector3d minimum, Eigen::Vector3d maximum) {
    min = minimum;
    max = maximum;
    //check if particles are in box
    std::vector<Particle> particles_in_box;
    for (Particle p: particles) {
        if (is_in_box(p)) {
            particles_in_box.push_back(p);
        }
    }
    if (particles_in_box.size() > 1) {
        type = Internal;
        init_children(particles_in_box, min, max);
        mass = 0;
        com = {0, 0, 0};
        for(auto &p : particles){
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

bool Octree::is_leaf() {
    return type == Leaf;
}

bool Octree::is_empty() {
    return type == Empty;
}

double Octree::side_length() {
    Eigen::Vector3d vec = max - min;
    return vec.lpNorm<Eigen::Infinity>();
}

Eigen::Vector3d Octree::force_on_particle(Particle &p, double softening, double theta) {
    Eigen::Vector3d force = {0, 0, 0};
    if(is_empty()) {
        return force;
    }
    if (is_leaf()) {
        //TODO
        if(&particle_in_box == &p){
            return force;
        }
        else{
            return gravitational_force(p, particle_in_box, softening);
        }
    }
    double d = side_length();
    double a = d/p.get_distance(com);
    if (a < theta) {
        Particle p2 = Particle(mass, com);
        force += gravitational_force(p, p2, softening);
    }
    else {
        for (int i = 0; i < 8; ++i) {
            force += children[i]->force_on_particle(p, softening, theta);
        }
    }
    return force;
}

void Octree::tree_boxes(Octree* node, std::vector<Eigen::Vector3d>& boxes, std::vector<Eigen::Vector3d>& minima,
                        std::vector<Eigen::Vector3d>& maxima) {
    if (node->is_leaf() || node->is_empty()) {
        return;
    }
    minima.push_back(node->get_min());
    maxima.push_back(node->get_max());
    boxes.push_back(node->get_center());
    for (int i = 0; i < 8; ++i) {
        tree_boxes(node->children[i].get(), boxes, minima, maxima);
    }
}

void Octree::centers_to_file(std::vector<Eigen::Vector3d>& boxes, const std::string& filename) {
    std::ofstream file;
    file.open(filename);
    for (auto & box : boxes) {
        file << box.x() << " " << box.y() << " " << box.z() << std::endl;
    }
    file.close();
}

void Octree::limits_to_file(std::vector<Eigen::Vector3d>& boxes, std::vector<Eigen::Vector3d>& minima, std::vector<Eigen::Vector3d>& maxima, const std::string& filename) {
    std::ofstream file;
    file.open(filename);
    for (int i = 0; i < boxes.size(); ++i) {
        file << minima[i].x() << " " << minima[i].y() << " " << minima[i].z() << " " << maxima[i].x() << " " << maxima[i].y() << " " << maxima[i].z() << std::endl;
    }
    file.close();
}


double Octree::count_boxes() {
    //doesnt count empty boxes
    if(is_empty()){
        return 0;
    }
    if (is_leaf()) {
        return 1;
    }
    else {
        double nodes = 1;
        for (int i = 0; i < 8; ++i) {
            nodes += children[i]->count_boxes();
        }
        return nodes;
    }
}

int Octree::count_empty_boxes(){
    int empty_boxes = 0;
    if (is_leaf()){
        return 0;
    }
    if(is_empty()){
        return 1;
    }
    else {
        for (int i = 0; i < 8; ++i) {
            empty_boxes += children[i]->count_empty_boxes();
        }
    }
    return empty_boxes;
}

Eigen::Vector3d Octree::get_min() {
    return min;
}

Eigen::Vector3d Octree::get_max() {
    return max;
}


