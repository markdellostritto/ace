#pragma once
#ifndef CALC_DIPOLE_LONG_HPP
#define CALC_DIPOLE_LONG_HPP

// c++
#include <iostream>
// eigen
#include <Eigen/Dense>
// mem
#include "mem/serialize.hpp"
// sim
#include "sim/calc.hpp"
#include "sim/kspace_dipole.hpp"

#ifndef CALC_DIPOLE_LONG_PRINT_FUNC
#define CALC_DIPOLE_LONG_PRINT_FUNC 0
#endif

class CalcDipoleLong: public Calculator{
private:
	double eps_{1.0};
	double prec_{1.0e-6};
	KSpace::Dipole kDipole_;
public:
	//==== constructors/destructors ====
	CalcDipoleLong():Calculator(Calculator::Name::DIPOLE_LONG){}
    CalcDipoleLong(double rc):Calculator(Calculator::Name::DIPOLE_LONG,rc){}
	~CalcDipoleLong(){}
	
	//==== operator ====
	friend std::ostream& operator<<(std::ostream& out, const CalcDipoleLong& calc);
	
	//==== access ====
	double& prec(){return prec_;}
	const double& prec()const{return prec_;}
	double& eps(){return eps_;}
	const double& eps()const{return eps_;}
	KSpace::Dipole& kDipole(){return kDipole_;}
	const KSpace::Dipole& kDipole()const{return kDipole_;}
	
	//==== member functions ====
    void init();
	void init(const Structure& struc);
	void read(Token& token);
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
	
	template <> int nbytes(const CalcDipoleLong& obj);
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcDipoleLong& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcDipoleLong& obj, const char* arr);
	
}

#endif