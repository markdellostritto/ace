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
	Calculator::Mix mix_{Mix::NONE};//mixing mode
	double prec_{1.0e-6};
	KSpace::London kLondon_;
	Eigen::VectorXd ie_;//ionization energy
	Eigen::VectorXd alpha_;//polarizability
    Eigen::MatrixXd rvdw_;//radius - vdw
	Eigen::MatrixXd rvdw6_;//rvdw^6
	Eigen::MatrixXd c6_;//c6 coefficient
public:
    //==== contructors/destructors ====
	CalcLDampLong():Calculator(Calculator::Name::LDAMP_LONG){}
    CalcLDampLong(double rc):Calculator(Calculator::Name::LDAMP_LONG,rc){}
    ~CalcLDampLong(){}

    //==== operators ====
    friend std::ostream& operator<<(std::ostream& out, const CalcLDampLong& calc);

	//==== access ====
	Calculator::Mix& mix(){return mix_;}
	const Calculator::Mix& mix()const{return mix_;}
    double& prec(){return prec_;}
	const double& prec()const{return prec_;}
	Eigen::VectorXd& ie(){return ie_;}
	const Eigen::VectorXd& ie()const{return ie_;}
	Eigen::VectorXd& alpha(){return alpha_;}
	const Eigen::VectorXd& alpha()const{return alpha_;}
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