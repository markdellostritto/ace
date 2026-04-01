#pragma once
#ifndef CALC_COUL_LONG_HPP
#define CALC_COUL_LONG_HPP

// c++
#include <iostream>
// eigen
#include <Eigen/Dense>
// mem
#include "mem/serialize.hpp"
// sim
#include "sim/calc.hpp"
#include "sim/kspace_coul.hpp"

#ifndef CALC_COUL_LONG_PRINT_FUNC
#define CALC_COUL_LONG_PRINT_FUNC 0
#endif

class CalcCoulLong: public Calculator{
private:
	double eps_{1.0};
	double prec_{1.0e-6};
	KSpace::Coul kCoul_;
public:
	//==== constructors/destructors ====
	CalcCoulLong():Calculator(Calculator::Name::COUL_LONG){}
	CalcCoulLong(double rc):Calculator(Calculator::Name::COUL_LONG,rc){};
	~CalcCoulLong(){}
	
	//==== operator ====
	friend std::ostream& operator<<(std::ostream& out, const CalcCoulLong& calc);
	
	//==== access ====
	double& prec(){return prec_;}
	const double& prec()const{return prec_;}
	double& eps(){return eps_;}
	const double& eps()const{return eps_;}
	KSpace::Coul& kCoul(){return kCoul_;}
	const KSpace::Coul& kCoul()const{return kCoul_;}
	
	//==== member functions ====
    void init();
	void init(const Structure& struc);
	void read(Token& token);
	void coeff(Token& token){}
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
	
	template <> int nbytes(const CalcCoulLong& obj);
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcCoulLong& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcCoulLong& obj, const char* arr);
	
}

#endif