//
// Created by noahe on 08/01/2024.
//

#include "Octree.h"

//PRE:: particles
Octree::Octree(std::vector<Particle> particles) {
    set_limits(particles);
    if (particles.size() > 1) {
        type = Internal;
        init_children(particles, min, max);
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

Octree::Octree(std::vector<Particle> particles, Particle minimum, Particle maximum) {
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

    max = Particle(x_max, y_max, z_max);
    min = Particle(x_min, y_min, z_min);
}

void Octree::init_children(std::vector<Particle> particles, Particle min, Particle max) {
    Particle center = this->get_center();
    children[0] = new Octree(particles, Particle(min.pos.x(), min.pos.y(), min.pos.z()), Particle(center.pos.x(), center.pos.y(), center.pos.z()));
    children[1] = new Octree(particles, Particle(center.pos.x(), min.pos.y(), min.pos.z()), Particle(max.pos.x(), center.pos.y(), center.pos.z()));
    children[2] = new Octree(particles, Particle(min.pos.x(), center.pos.y(), min.pos.z()), Particle(center.pos.x(), max.pos.y(), center.pos.z()));
    children[3] = new Octree(particles, Particle(center.pos.x(), center.pos.y(), min.pos.z()), Particle(max.pos.x(), max.pos.y(), center.pos.z()));
    children[4] = new Octree(particles, Particle(min.pos.x(), min.pos.y(), center.pos.z()), Particle(center.pos.x(), center.pos.y(), max.pos.z()));
    children[5] = new Octree(particles, Particle(center.pos.x(), min.pos.y(), center.pos.z()), Particle(max.pos.x(), center.pos.y(), max.pos.z()));
    children[6] = new Octree(particles, Particle(min.pos.x(), center.pos.y(), center.pos.z()), Particle(center.pos.x(), max.pos.y(), max.pos.z()));
    children[7] = new Octree(particles, Particle(center.pos.x(), center.pos.y(), center.pos.z()), Particle(max.pos.x(), max.pos.y(), max.pos.z()));
}

//TODO
//Octree::~Octree() {
//    cleanupChildren(this);
//}
//
//void Octree::cleanupChildren(Octree* nodes) {
//    if (nodes->is_leaf() || nodes->is_empty()) {
//        return;
//    }
//    for (int i = 0; i < 8; ++i) {
//        cleanupChildren(nodes->children[i]);
//    }
//    delete nodes;
//}

Particle Octree::get_center() {
    return Particle(0, Eigen::Vector3d((max.pos.x() + min.pos.x()) / 2, (max.pos.y() + min.pos.y()) / 2,
                                       (max.pos.z() + min.pos.z()) / 2), Eigen::Vector3d(0, 0, 0));
}

bool Octree::is_in_box(Particle p) {
    return (p.pos.x() >= min.pos.x() && p.pos.x() <= max.pos.x() && p.pos.y() >= min.pos.y() && p.pos.y() <= max.pos.y()
            && p.pos.z() >= min.pos.z() && p.pos.z() <= max.pos.z());
}

bool Octree::is_leaf() {
    return type == Leaf;
}

bool Octree::is_empty() {
    return type == Empty;
}

double Octree::mass {
    if(is_empty()){
        return 0;
    }
    if(is_leaf()){
        return particle_in_box.mass;
    }
    double mass = 0;
    for (auto & i : children) {
        mass += i->mass;
    }
    return mass;
}

double Octree::side_length() {
    Eigen::Vector3d vec = max.pos - min.pos;
    return vec.lpNorm<Eigen::Infinity>();
}

Eigen::Vector3d Octree::get_center_of_mass() {
    Eigen::Vector3d center_of_mass = {0, 0, 0};
    if(is_empty()){
        return center_of_mass;
    }
    if(is_leaf()){
        return particle_in_box.pos * particle_in_box.mass;
    }
    for(auto & i : children){
        center_of_mass += i->get_center_of_mass();
    }
    return center_of_mass/this->mass;
}

Eigen::Vector3d Octree::force_on_particle(Particle p, double softening, double theta) {
    Eigen::Vector3d force = {0, 0, 0};
    if(is_empty()) {
        return force;
    }
    if (is_leaf()) {
        if(particle_in_box == p){
            return force;
        }
        else{
            return gravitational_force(p, particle_in_box, softening);
        }
    }
    else {
        double d = side_length();
        for (int i = 0; i < 8; ++i) {
            double a = d/p.get_distance(children[i]->get_center_of_mass());
            if (a < theta) {
                Particle p2 = Particle(children[i]->mass, children[i]->get_center_of_mass());
                force += gravitational_force(p, p2, softening);
            }
            else{
                force += children[i]->force_on_particle(p, softening, theta);
            }
        }
    }
    return force;
}

void Octree::tree_boxes(Octree* node, std::vector<Particle>& boxes, std::vector<Particle>& minima, std::vector<Particle>& maxima) {
    if (node->is_leaf() || node->is_empty()) {
        return;
    }
    minima.push_back(node->get_min());
    maxima.push_back(node->get_max());
    boxes.push_back(node->get_center());
    for (int i = 0; i < 8; ++i) {
        tree_boxes(node->children[i], boxes, minima, maxima);
    }
}

void Octree::centers_to_file(std::vector<Particle>& boxes, const std::string& filename) {
    std::ofstream file;
    file.open(filename);
    for (auto & box : boxes) {
        file << box.pos.x() << " " << box.pos.y() << " " << box.pos.z() << std::endl;
    }
    file.close();
}

void Octree::limits_to_file(std::vector<Particle>& boxes, std::vector<Particle>& minima, std::vector<Particle>& maxima, const std::string& filename) {
    std::ofstream file;
    file.open(filename);
    for (int i = 0; i < boxes.size(); ++i) {
        file << minima[i].pos.x() << " " << minima[i].pos.y() << " " << minima[i].pos.z() << " " << maxima[i].pos.x() << " " << maxima[i].pos.y() << " " << maxima[i].pos.z() << std::endl;
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

Particle Octree::get_min() {
    return min;
}

Particle Octree::get_max() {
    return max;
}


