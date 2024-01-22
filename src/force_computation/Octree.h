//
// Created by noahe on 08/01/2024.
//

#ifndef NBODYPROBLEM_OCTREE_H
#define NBODYPROBLEM_OCTREE_H

#include <Eigen/Core>
#include <vector>
#include "../shared/particle.h"
#include "force_helper.h"


class Octree {
public:
    enum node_type{
        Leaf,
        Internal,
        Empty
    };

    Octree() = default;
//    ~Octree();

    explicit Octree(std::vector<Particle> particles);
    Octree(std::vector<Particle> particles, Particle min, Particle max);

    //getters
    Particle get_center();
    double mass;
    Eigen::Vector3d get_center_of_mass();
    bool is_leaf();
    bool is_empty();

    double side_length();
    bool is_in_box(Particle p);
    Eigen::Vector3d force_on_particle(Particle p, double softening, double theta);

    //TESTING
    double count_boxes();
    int count_empty_boxes();
    void centers_to_file(std::vector<Particle> &boxes, const std::string &filename);
    void limits_to_file(std::vector<Particle> &boxes, std::vector<Particle> &minima,
                        std::vector<Particle> &maxima, const std::string &filename);
    void tree_boxes(Octree *node, std::vector<Particle> &boxes, std::vector<Particle>& minima, std::vector<Particle>& maxima);
    Particle get_min();
    Particle get_max();

private:
    node_type type;
    Particle min;
    Particle max;
    Particle particle_in_box;
    Octree* children[8];

    void cleanupChildren(Octree* nodes);
    void set_limits(std::vector<Particle> particles);
    void init_children(std::vector<Particle> particles, Particle min, Particle max);
};

#endif //NBODYPROBLEM_OCTREE_H
