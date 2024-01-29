cmake_minimum_required(VERSION 3.22.1)

project(
  NBodyProblem
  VERSION 1.0
  LANGUAGES CXX)

set (CMAKE_CXX_STANDARD 20)

#-fsanitize=address,undefined,leak -D_GLIBCXX_DEBUG ") -Wall -Wextra -Wpedantic -Werror -std=c++17 -O3 -march=native -mtune=native -fopenmp -fno-omit-frame-pointer -fsanitize=address,undefined,leak -D_GLIBCXX_DEBUG
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ")
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(Boost REQUIRED COMPONENTS iostreams)
find_package(OpenMP REQUIRED)

include_directories(/usr/bin/gnuplot-iostream)

add_executable(density src/density/density.cpp
        src/shared/particle.cpp
        src/shared/data_io.cpp)

add_executable(force_computation
        src/force_computation/force_computation.cpp
        src/force_computation/Octree.cpp
        src/force_computation/force_helper.cpp
        src/shared/particle.cpp
        src/shared/data_io.cpp)

add_executable(test_sim
        src/force_computation/test_sim.cpp
        src/force_computation/Simulation.cpp
        src/force_computation/force_helper.cpp
        src/shared/particle.cpp
        src/shared/data_io.cpp
        src/force_computation/Octree.cpp)

add_executable(test
        src/force_computation/test_octree.cpp
        src/force_computation/Octree.cpp
        src/force_computation/force_helper.cpp
        src/shared/particle.cpp
        src/shared/data_io.cpp)

target_link_libraries(density Eigen3::Eigen Boost::iostreams)
target_link_libraries(force_computation Eigen3::Eigen OpenMP::OpenMP_CXX)
target_link_libraries(test_sim Eigen3::Eigen)
target_link_libraries(test Eigen3::Eigen OpenMP::OpenMP_CXX)
