//
// Created by noahe on 10/01/2024.
//

#ifndef NBODYPROBLEM_NON_LINEAR_SOLVER_H
#define NBODYPROBLEM_NON_LINEAR_SOLVER_H

#include <iostream>
#include <Eigen/Dense>

#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

// Generic functor
template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor
{
    typedef _Scalar Scalar;
    enum {
        InputsAtCompileTime = NX,
        ValuesAtCompileTime = NY
    };
    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

    int m_inputs, m_values;

    Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
    Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

    int inputs() const { return m_inputs; }
    int values() const { return m_values; }

};

struct my_functor : Functor<double>
{
    std::vector<double> measured_density;
    double bin_size;
    double total_mass;

    my_functor(const std::vector<double>& measured_density, double bin_size, double total_mass)
        : Functor<double>(1, measured_density.size() - 1), measured_density(measured_density), bin_size(bin_size), total_mass(total_mass) {}

    int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
    {
        double scale_length = x(0);
        for (size_t i = 1; i < measured_density.size(); ++i) {
            double radius = bin_size*i;
            double approximate_density = total_mass/(2*M_PI) * scale_length/radius *(1/pow(scale_length+radius,3));
            fvec(i - 1) = measured_density[i] - approximate_density;
        }
        return 0;
    }
};



//
//int main(int argc, char *argv[])
//{
//    // ... code to calculate measured_density, radii, and total_mass ...
//
//    Eigen::VectorXd a(1);
//    a(0) = 1; // initial guess for scale_length
//
//    my_functor functor(measured_density, radii, total_mass);
//    Eigen::NumericalDiff<my_functor> numDiff(functor);
//    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<my_functor>,double> lm(numDiff);
//    lm.parameters.maxfev = 2000;
//    lm.parameters.xtol = 1.0e-10;
//
//    int ret = lm.minimize(a);
//
//    std::cout << "Optimal scale_length: " << a(0) << std::endl;
//
//    return 0;
//}

#endif //NBODYPROBLEM_NON_LINEAR_SOLVER_H
