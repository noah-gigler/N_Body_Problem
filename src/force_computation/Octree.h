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
    Eigen::Vector3d force_on_particle(Particle &p, double softening, double theta);

private:
    bool is_in_box(const Particle &p);
    void set_limits(std::vector<Particle> particles);
    void init_tree(std::vector<Particle> &particles);
    void init_children(std::vector<Particle> &particles, Eigen::Vector3d minimum, Eigen::Vector3d maximum);
    Eigen::Vector3d get_center();
    double side_length();

    node_type type;
    Eigen::Vector3d min;
    Eigen::Vector3d max;
    Eigen::Vector3d com;
    double mass;
    Particle particle_in_box;
    std::unique_ptr<Octree> children[8];
};

#endif //NBODYPROBLEM_OCTREE_H
