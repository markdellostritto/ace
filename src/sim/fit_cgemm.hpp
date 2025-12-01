#pragma once
#ifndef FIT_CGEMM_HPP
#define FIT_CGEMM_HPP

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
// CGemmType
//*****************************************************

class CGemmType{
private:
    std::string name_{"NONE"};
    double radius_{0.0};
    double rvdw_{0.0};
    double aOver_{0.0};
    double aRep_{0.0};
    double mass_{0.0};
    int index_{-1};
public:
    //==== constructors/destructors ====
    CGemmType(){}
    ~CGemmType(){}

    //==== operators ====
    friend std::ostream& operator<<(std::ostream& out, const CGemmType& type);

    //==== access ====
    std::string& name(){return name_;}
    const std::string& name()const{return name_;}
    double& radius(){return radius_;}
    const double& radius()const{return radius_;}
    double& rvdw(){return rvdw_;}
    const double& rvdw()const{return rvdw_;}
    double& aOver(){return aOver_;}
    const double& aOver()const{return aOver_;}
    double& aRep(){return aRep_;}
    const double& aRep()const{return aRep_;}
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
    Structure strucA_;//structure - analytical (cores only) 
    Structure strucN_;//structure - numerical (cores + electrons)
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
    std::vector<CGemmType> types;
    std::vector<CubeData> data;
    Engine engine;
    std::shared_ptr<Integrator> intg;
};

#endif