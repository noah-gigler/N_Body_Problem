//
// Created by noahe on 08/01/2024.
//

#ifndef NBODYPROBLEM_OCTREE_H
#define NBODYPROBLEM_OCTREE_H

#include <Eigen/Core>
#include <vector>
#include <memory>
#include "../shared/particle.h"
#include "force_helper.h"


class Octree {
public:
    enum node_type{
        Empty,
        Leaf,
        Internal,
    };

    Octree() = default;
    explicit Octree(std::vector<Particle> &particles);
    Octree(std::vector<Particle> &particles, Eigen::Vector3d minimum, Eigen::Vector3d maximum);

    //getters
    Eigen::Vector3d get_center();
    bool is_leaf();
    bool is_empty();
    double side_length();
    bool is_in_box(const Particle &p);
    Eigen::Vector3d force_on_particle(Particle &p, double softening, double theta);

    //TESTING
//    double count_boxes();
//    int count_empty_boxes();
//    void centers_to_file(std::vector<Eigen::Vector3d> &boxes, const std::string &filename);
//    void limits_to_file(std::vector<Eigen::Vector3d> &boxes, std::vector<Eigen::Vector3d> &minima,
//                        std::vector<Eigen::Vector3d> &maxima, const std::string &filename);
//    void tree_boxes(Octree *node, std::vector<Eigen::Vector3d> &boxes,
//                    std::vector<Eigen::Vector3d>& minima, std::vector<Eigen::Vector3d>& maxima);
//    Eigen::Vector3d get_min();
//    Eigen::Vector3d get_max();

public:
    Eigen::Vector3d com;
    double mass;

private:
    node_type type;
    Eigen::Vector3d min;
    Eigen::Vector3d max;
    Particle particle_in_box;
    std::unique_ptr<Octree> children[8];

    void set_limits(std::vector<Particle> particles);
    void init_tree(std::vector<Particle> &particles);
    void init_children(std::vector<Particle> &particles, Eigen::Vector3d minimum, Eigen::Vector3d maximum);
};

#endif //NBODYPROBLEM_OCTREE_H
