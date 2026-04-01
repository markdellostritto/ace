#pragma once
#ifndef CALC_LJ_LONG_HPP
#define CALC_LJ_LONG_HPP

// c++
#include <iostream>
// eigen
#include <Eigen/Dense>
// mem
#include "mem/serialize.hpp"
// sim
#include "sim/calc.hpp"
#include "sim/kspace_london.hpp"

#ifndef CALC_LJ_LONG_PRINT_FUNC
#define CALC_LJ_LONG_PRINT_FUNC 0
#endif

class CalcLJLong: public Calculator{
private:
	double prec_{1.0e-6};
	KSpace::London kLondon_;
	Eigen::MatrixXi f_;//flag (t/f)
    Eigen::MatrixXd e_;//epsilon
    Eigen::MatrixXd s_;//sigma
public:
    //==== contructors/destructors ====
	CalcLJLong():Calculator(Calculator::Name::LJ_LONG){}
    CalcLJLong(double rc):Calculator(Calculator::Name::LJ_LONG,rc){}
    ~CalcLJLong(){}

    //==== operators ====
    friend std::ostream& operator<<(std::ostream& out, const CalcLJLong& calc);

	//==== access ====
	double& prec(){return prec_;}
	const double& prec()const{return prec_;}
	Eigen::MatrixXi& f(){return f_;}
	const Eigen::MatrixXi& f()const{return f_;}
	Eigen::MatrixXd& s(){return s_;}
	const Eigen::MatrixXd& s()const{return s_;}
	Eigen::MatrixXd& e(){return e_;}
	const Eigen::MatrixXd& e()const{return e_;}
	
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
	
	template <> int nbytes(const CalcLJLong& obj);
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcLJLong& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcLJLong& obj, const char* arr);
	
}

#endif