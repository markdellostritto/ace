#pragma once
#ifndef CALC_DIPOLE_CUT_HPP
#define CALC_DIPOLE_CUT_HPP

// c++
#include <iostream>
// eigen
#include <Eigen/Dense>
// mem
#include "mem/serialize.hpp"
// sim
#include "sim/calc.hpp"

#ifndef CALC_DIPOLE_CUT_PRINT_FUNC
#define CALC_DIPOLE_CUT_PRINT_FUNC 0
#endif

class CalcDipoleCut: public Calculator{
public:
    //==== contructors/destructors ====
	CalcDipoleCut():Calculator(Calculator::Name::DIPOLE_CUT){}
    CalcDipoleCut(double rc):Calculator(Calculator::Name::DIPOLE_CUT,rc){}
    ~CalcDipoleCut(){}

    //==== operators ====
    friend std::ostream& operator<<(std::ostream& out, const CalcDipoleCut& calc);
    
    //==== member functions ====
    double energy(Structure& struc, const NeighborList& nlist)const;
	double compute(Structure& struc, const NeighborList& nlist)const;
};

//**********************************************
// serialization
//**********************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const CalcDipoleCut& obj);
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcDipoleCut& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcDipoleCut& obj, const char* arr);
	
}

#endif