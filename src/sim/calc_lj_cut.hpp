#pragma once
#ifndef CALC_LJ_CUT_HPP
#define CALC_LJ_CUT_HPP

// c++
#include <iostream>
// eigen
#include <Eigen/Dense>
// mem
#include "mem/serialize.hpp"
// sim
#include "sim/calc.hpp"

#ifndef CALC_LJ_CUT_PRINT_FUNC
#define CALC_LJ_CUT_PRINT_FUNC 0
#endif

class CalcLJCut: public Calculator{
private:
    Eigen::MatrixXi f_;//flag (t/f)
    Eigen::MatrixXd e_;//epsilon
    Eigen::MatrixXd s_;//sigma
public:
    //==== contructors/destructors ====
	CalcLJCut():Calculator(Calculator::Name::LJ_CUT){}
    CalcLJCut(double rc):Calculator(Calculator::Name::LJ_CUT,rc){}
    ~CalcLJCut(){}

    //==== operators ====
    friend std::ostream& operator<<(std::ostream& out, const CalcLJCut& calc);

	//==== access ====
	Eigen::MatrixXi& f(){return f_;}
	const Eigen::MatrixXi& f()const{return f_;}
	Eigen::MatrixXd& s(){return s_;}
	const Eigen::MatrixXd& s()const{return s_;}
	Eigen::MatrixXd& e(){return e_;}
	const Eigen::MatrixXd& e()const{return e_;}
	
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
	
	template <> int nbytes(const CalcLJCut& obj);
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcLJCut& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcLJCut& obj, const char* arr);
	
}

#endif