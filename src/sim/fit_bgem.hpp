#pragma once
#ifndef FIT_CGEM2_HPP
#define FIT_CGEM2_HPP

// c++
#include <iostream>
#include <string>
// str
#include "str/token.hpp"
// struc
#include "struc/structure.hpp"
#include "struc/grid.hpp"
// sim
#include "sim/engine.hpp"
#include "sim/integrator.hpp"

#ifndef PRINT_STRUC
#define PRINT_STRUC 0
#endif

//*****************************************************
// BGemType
//*****************************************************

class BGemType{
private:
    std::string name_;
    double alpha_;
    double amp_;
    double rep_;
    double mass_;
    int index_;
public:
    //==== constructors/destructors ====
    BGemType();
    ~BGemType(){}

    //==== operators ====
    friend std::ostream& operator<<(std::ostream& out, const BGemType& type);

    //==== access ====
    std::string& name(){return name_;}
    const std::string& name()const{return name_;}
    double& alpha(){return alpha_;}
    const double& alpha()const{return alpha_;}
    double& amp(){return amp_;}
    const double& amp()const{return amp_;}
    double& rep(){return rep_;}
    const double& rep()const{return rep_;}
    double& mass(){return mass_;}
    const double& mass()const{return mass_;}
    int& index(){return index_;}
    const int& index()const{return index_;}

    //==== member functions ====
    void read(Token& token);
};

//*****************************************************
// Data
//*****************************************************

class CubeData{
private:
    Structure strucA_;//structure
    Structure strucN_;//structure
    Grid espA_;//electrostatic potential - analytical 
    Grid espN_;//electrostatic potential - numerical 
public:
    //==== constructors/destructors ====
    CubeData(){}
    ~CubeData(){}

    //==== access ====
    Structure& strucA(){return strucA_;}
    const Structure& strucA()const{return strucA_;}
    Structure& strucN(){return strucN_;}
    const Structure& strucN()const{return strucN_;}
    Grid& espA(){return espA_;}
    const Grid& espA()const{return espA_;}
    Grid& espN(){return espN_;}
    const Grid& espN()const{return espN_;}
};

//*****************************************************
// Function Data
//*****************************************************

struct FunctionData{
    int nsteps;
    double ftol;
    std::vector<BGemType> types;
    std::vector<CubeData> data;
    Engine engine;
    std::shared_ptr<Integrator> intg;
};

#endif