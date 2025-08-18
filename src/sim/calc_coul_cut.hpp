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
public:
    //==== contructors/destructors ====
	CalcCoulCut():Calculator(Calculator::Name::LJ_CUT){}
    CalcCoulCut(double rc):Calculator(Calculator::Name::LJ_CUT,rc){}
    ~CalcCoulCut(){}

    //==== operators ====
    friend std::ostream& operator<<(std::ostream& out, const CalcCoulCut& calc);
    
    //==== member functions ====
    void resize(int ntypes);
	void init();
	void read(Token& token);
	void coeff(Token& token);
	double energy(Structure& struc, const NeighborList& nlist);
	double energy(Structure& struc);
	double compute(Structure& struc, const NeighborList& nlist);
	double compute(Structure& struc);
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