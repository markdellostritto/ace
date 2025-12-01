#pragma once
#ifndef THOLE_HPP
#define THOLE_HPP

// c++
#include <iostream>
#include <vector>
// eigen
#include <Eigen/Dense>
// sim
#include "sim/calc_thole_long.hpp"

#ifndef DEBUG_THOLE
#define DEBUG_THOLE 0
#endif

class Thole{
private:
    //calc flags
    bool inter_{true};//whether or not to include interaction matrices
    bool alphar_{true};//whether or not to use volume-dependent polarizabilities
    //thole dipole interaction
    CalcTholeLong calc_;
    //matrix utilities
    Eigen::MatrixXd alphaC_;
    Eigen::MatrixXd identityC_;
    Eigen::MatrixXd A_;
    //atomic radii
    std::vector<double> r_,r0_;
    std::vector<Eigen::Vector3d> drAtom_;
    //polarizability
    std::vector<double> alpha_;
public:
    //==== constructors/destructors ====
    Thole(){}
    ~Thole(){}

    //==== operators ====
    friend std::ostream& operator<<(std::ostream& out, Thole& thole);

    //==== access ===
    //flags
    bool& inter(){return inter_;}
    const bool& inter()const{return inter_;}
    bool& alphar(){return alphar_;}
    const bool& alphar()const{return alphar_;}
    //dipole interaction
    CalcTholeLong& calc(){return calc_;}
    const CalcTholeLong& calc()const{return calc_;}
    //polarizability
    std::vector<double>& alpha(){return alpha_;}
    const std::vector<double>& alpha()const{return alpha_;}
    
    //==== member functions ====
    void init(){calc_.init();}
    void init(const Structure& struc){calc_.init(struc);}
    Eigen::Matrix3d& compute(const Structure& struc, Eigen::Matrix3d& alpha);
};

#endif