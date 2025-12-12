#pragma once
#ifndef CALC_CGEMM_CUT_HPP
#define CALC_CGEMM_CUT_HPP

// c++
#include <iostream>
// eigen
#include <Eigen/Dense>
// mem
#include "mem/serialize.hpp"
// sim
#include "sim/calc.hpp"

#ifndef CALC_CGEMM_CUT_PRINT_FUNC
#define CALC_CGEMM_CUT_PRINT_FUNC 0
#endif

class CalcCGemmCut: public Calculator{
private:
	//global parameters
	Calculator::Mix mix_{Mix::NONE};
    double lambdaC_{1.0};//radial scaling factor - Coulomb
    double lambdaS_{1.0};//radial scaling factor - Overlap
	double rRep_{0.0};//repulsive radius
	//type parameters
    Eigen::VectorXd radius_;//radius
    Eigen::MatrixXd aOver_;//overlap amplitude
	Eigen::MatrixXd aRep_;//repulsive amplitude
    Eigen::MatrixXd gammaC_;//reduced expontent - Coulomb
    Eigen::MatrixXd gammaS_;//reduced expontent - Overlap
	Eigen::MatrixXd rgammaC_;//reduced expontent - Coulomb
public:
    //==== contructors/destructors ====
	CalcCGemmCut():Calculator(Calculator::Name::CGEMM_CUT){}
    CalcCGemmCut(double rc):Calculator(Calculator::Name::CGEMM_CUT,rc){}
	CalcCGemmCut(double rc, double lambdaC, double lambdaS, double rRep, Calculator::Mix mix);
    ~CalcCGemmCut(){}

    //==== operators ====
    friend std::ostream& operator<<(std::ostream& out, const CalcCGemmCut& calc);

	//==== access ====
	//global parameters
	Calculator::Mix& mix(){return mix_;}
	const Calculator::Mix& mix()const{return mix_;}
    double& lambdaC(){return lambdaC_;}
    const double& lambdaC()const{return lambdaC_;}
	double& lambdaS(){return lambdaS_;}
    const double& lambdaS()const{return lambdaS_;}
	double& rRep(){return rRep_;}
    const double& rRep()const{return rRep_;}
	//type parameters
	Eigen::VectorXd& radius(){return radius_;}
	const Eigen::VectorXd& radius()const{return radius_;}
	Eigen::MatrixXd& aOver(){return aOver_;}
	const Eigen::MatrixXd& aOver()const{return aOver_;}
	Eigen::MatrixXd& aRep(){return aRep_;}
	const Eigen::MatrixXd& aRep()const{return aRep_;}
	Eigen::MatrixXd& gammaC(){return gammaC_;}
	const Eigen::MatrixXd& gammaC()const{return gammaC_;}
	Eigen::MatrixXd& gammaS(){return gammaS_;}
	const Eigen::MatrixXd& gammaS()const{return gammaS_;}
	
    //==== member functions ====
    void resize(int ntypes);
	void init();
	void read(Token& token);
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
	
	template <> int nbytes(const CalcCGemmCut& obj);
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcCGemmCut& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcCGemmCut& obj, const char* arr);
	
}

#endif