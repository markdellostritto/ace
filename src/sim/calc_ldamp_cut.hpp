#pragma once
#ifndef CALC_LDAMP_CUT_HPP
#define CALC_LDAMP_CUT_HPP

// c++
#include <iostream>
// eigen
#include <Eigen/Dense>
// mem
#include "mem/serialize.hpp"
// sim
#include "sim/calc.hpp"

#ifndef CALC_LDAMP_CUT_PRINT_FUNC
#define CALC_LDAMP_CUT_PRINT_FUNC 0
#endif

class CalcLDampCut: public Calculator{
private:
    Eigen::MatrixXi f_;//flag (t/f)
    Eigen::MatrixXd rvdw_;//vdw radius
	Eigen::MatrixXd rvdw6_;//(vdw radius)^6
	Eigen::MatrixXd c6_;//london c6 coefficient
public:
    //==== contructors/destructors ====
	CalcLDampCut():Calculator(Calculator::Name::LDAMP_CUT){}
    CalcLDampCut(double rc):Calculator(Calculator::Name::LDAMP_CUT,rc){}
    ~CalcLDampCut(){}

    //==== operators ====
    friend std::ostream& operator<<(std::ostream& out, const CalcLDampCut& calc);

	//==== access ====
	Eigen::MatrixXi& f(){return f_;}
	const Eigen::MatrixXi& f()const{return f_;}
	Eigen::MatrixXd& rvdw(){return rvdw_;}
	const Eigen::MatrixXd& rvdw()const{return rvdw_;}
	Eigen::MatrixXd& c6(){return c6_;}
	const Eigen::MatrixXd& c6()const{return c6_;}
	
    //==== member functions ====
    void resize(int ntypes);
	void init();
	void coeff(Token& token);
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
	
	template <> int nbytes(const CalcLDampCut& obj);
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcLDampCut& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcLDampCut& obj, const char* arr);
	
}

#endif