#pragma once
#ifndef CALC_COUL_CUT_HPP
#define CALC_COUL_CUT_HPP

// c++
#include <iostream>
// eigen
#include <Eigen/Dense>
// mem
#include "mem/serialize.hpp"
// sim
#include "sim/calc.hpp"

#ifndef CALC_COUL_CUT_PRINT_FUNC
#define CALC_COUL_CUT_PRINT_FUNC 0
#endif

class CalcCoulCut: public Calculator{
private:
	double eps_{1.0};
public:
    //==== contructors/destructors ====
	CalcCoulCut():Calculator(Calculator::Name::COUL_CUT){}
    CalcCoulCut(double rc):Calculator(Calculator::Name::COUL_CUT,rc){}
    ~CalcCoulCut(){}

    //==== operators ====
    friend std::ostream& operator<<(std::ostream& out, const CalcCoulCut& calc);
    
	//==== access ====
	double& eps(){return eps_;}
    const double& eps()const{return eps_;}
	
    //==== member functions ====
    double energy(Structure& struc, const NeighborList& nlist)const;
	double energy(Structure& struc)const;
	double compute(Structure& struc, const NeighborList& nlist)const;
	double compute(Structure& struc)const;
};

//**********************************************
// serialization
//**********************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const CalcCoulCut& obj);
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcCoulCut& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcCoulCut& obj, const char* arr);
	
}

#endif