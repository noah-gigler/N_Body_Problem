//
// Created by noahe on 15/01/2024.
//
#include "../shared/particle.h"
#include "Octree.h"
#include "vector"

int main() {
    std::vector<Particle> particles = buildFromData();
    Octree tree = Octree(particles);
    //sanity check
    std::cout << "boxes: " << tree.count_boxes() << std::endl;
    std::cout << "empty boxes: " << tree.count_empty_boxes() << std::endl;

    std::vector<Particle> boxes;
    std::vector<Particle> minima;
    std::vector<Particle> maxima;
    tree.tree_boxes(&tree, boxes, minima, maxima);
    //sort boxes by x_value
    std::sort(boxes.begin(), boxes.end(), [](Particle a, Particle b) {
        return a.pos.x() < b.pos.x();
    });
    tree.centers_to_file(boxes, "../data/boxes.txt");
    tree.limits_to_file(boxes, minima, maxima, "../data/limits.txt");
}