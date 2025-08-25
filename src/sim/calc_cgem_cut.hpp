#pragma once
#ifndef CALC_CGEM_CUT_HPP
#define CALC_CGEM_CUT_HPP

// c++
#include <iostream>
// eigen
#include <Eigen/Dense>
// mem
#include "mem/serialize.hpp"
// sim
#include "sim/calc.hpp"

#ifndef CALC_CGEM_CUT_PRINT_FUNC
#define CALC_CGEM_CUT_PRINT_FUNC 0
#endif

class CalcCGemCut: public Calculator{
private:
    double lambdaC_;//radial scaling factor - Coulomb
    double lambdaS_;//radial scaling factor - Overlap
    Eigen::VectorXd radius_;//radius
    Eigen::MatrixXd aOver_;//overlap amplitude
	Eigen::MatrixXd aRep_;//repulsive amplitude
    Eigen::MatrixXd muC_;//reduced expontent - Coulomb
    Eigen::MatrixXd muS_;//reduced expontent - Overlap
    Eigen::MatrixXd rmuC_;//reduced expontent - Coulomb
public:
    //==== contructors/destructors ====
	CalcCGemCut();
    CalcCGemCut(double rc);
	CalcCGemCut(double rc, double lambdaC, double lambdaS);
    ~CalcCGemCut(){}

    //==== operators ====
    friend std::ostream& operator<<(std::ostream& out, const CalcCGemCut& calc);

	//==== access ====
    double& lambdaC(){return lambdaC_;}
    const double& lambdaC()const{return lambdaC_;}
	double& lambdaS(){return lambdaS_;}
    const double& lambdaS()const{return lambdaS_;}
	Eigen::VectorXd& radius(){return radius_;}
	const Eigen::VectorXd& radius()const{return radius_;}
	Eigen::MatrixXd& aOver(){return aOver_;}
	const Eigen::MatrixXd& aOver()const{return aOver_;}
	Eigen::MatrixXd& aRep(){return aRep_;}
	const Eigen::MatrixXd& aRep()const{return aRep_;}
	Eigen::MatrixXd& muC(){return muC_;}
	const Eigen::MatrixXd& muC()const{return muC_;}
	Eigen::MatrixXd& muS(){return muS_;}
	const Eigen::MatrixXd& muS()const{return muS_;}
	
    //==== member functions ====
    void resize(int ntypes);
	void init();
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
	
	template <> int nbytes(const CalcCGemCut& obj);
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcCGemCut& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcCGemCut& obj, const char* arr);
	
}

#endif