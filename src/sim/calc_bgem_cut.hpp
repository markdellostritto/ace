#pragma once
#ifndef CALC_BGEM_CUT_HPP
#define CALC_BGEM_CUT_HPP

// c++
#include <iostream>
// eigen
#include <Eigen/Dense>
// mem
#include "mem/serialize.hpp"
// sim
#include "sim/calc.hpp"

#ifndef CALC_BGEM_CUT_PRINT_FUNC
#define CALC_BGEM_CUT_PRINT_FUNC 0
#endif

class CalcBGemCut: public Calculator{
private:
	Eigen::VectorXd alpha_;
	Eigen::MatrixXd amp_;//overlap amplitude
	Eigen::MatrixXd rep_;//repulsive amplitude
    Eigen::MatrixXd mu_;//reduced expontent - Overlap
    Eigen::MatrixXd rmu_;//reduced expontent - Coulomb
public:
    //==== contructors/destructors ====
	CalcBGemCut():Calculator(Calculator::Name::BGEM_CUT){}
    CalcBGemCut(double rc):Calculator(Calculator::Name::BGEM_CUT,rc){}
    ~CalcBGemCut(){}

    //==== operators ====
    friend std::ostream& operator<<(std::ostream& out, const CalcBGemCut& calc);

	//==== access ====
    Eigen::VectorXd& alpha(){return alpha_;}
	const Eigen::VectorXd& alpha()const{return alpha_;}
	Eigen::MatrixXd& amp(){return amp_;}
	const Eigen::MatrixXd& amp()const{return amp_;}
	Eigen::MatrixXd& rep(){return rep_;}
	const Eigen::MatrixXd& rep()const{return rep_;}
	Eigen::MatrixXd& mu(){return mu_;}
	const Eigen::MatrixXd& mu()const{return mu_;}
	
    //==== member functions ====
    void resize(int ntypes);
	void init();
	void read(Token& token);
	void coeff(Token& token);
	double energy(Structure& struc, const NeighborList& nlist);
	double compute(Structure& struc, const NeighborList& nlist);
};

//**********************************************
// serialization
//**********************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const CalcBGemCut& obj);
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcBGemCut& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcBGemCut& obj, const char* arr);
	
}

#endif