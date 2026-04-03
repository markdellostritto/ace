#pragma once
#ifndef CALC_LDAMP_LONG_HPP
#define CALC_LDAMP_LONG_HPP

// c++
#include <iostream>
// eigen
#include <Eigen/Dense>
// mem
#include "mem/serialize.hpp"
// sim
#include "sim/calc.hpp"
#include "sim/kspace_london.hpp"

#ifndef CALC_LDAMP_LONG_PRINT_FUNC
#define CALC_LDAMP_LONG_PRINT_FUNC 0
#endif

class CalcLDampLong: public Calculator{
private:
	double prec_{1.0e-6};
	KSpace::London kLondon_;
	Eigen::MatrixXi f_;
	Eigen::MatrixXd rvdw_;
	Eigen::MatrixXd rvdw6_;
	Eigen::MatrixXd c6_;
public:
    //==== contructors/destructors ====
	CalcLDampLong():Calculator(Calculator::Name::LDAMP_LONG){}
    CalcLDampLong(double rc):Calculator(Calculator::Name::LDAMP_LONG,rc){}
    ~CalcLDampLong(){}

    //==== operators ====
    friend std::ostream& operator<<(std::ostream& out, const CalcLDampLong& calc);

	//==== access ====
	double& prec(){return prec_;}
	const double& prec()const{return prec_;}
	Eigen::MatrixXi& f(){return f_;}
	const Eigen::MatrixXi& f()const{return f_;}
	Eigen::MatrixXd& rvdw(){return rvdw_;}
	const Eigen::MatrixXd& rvdw()const{return rvdw_;}
	Eigen::MatrixXd& c6(){return c6_;}
	const Eigen::MatrixXd& c6()const{return c6_;}
	
    //==== member functions ====
    void resize(int ntypes);
	void init();
	void init(const Structure& struc);
	void read(Token& token);
	void coeff(Token& token);
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
	
	template <> int nbytes(const CalcLDampLong& obj);
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcLDampLong& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcLDampLong& obj, const char* arr);
	
}

#endif