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

/**
 * class CalcCGemmCut
 * Calculation of the energy of a set of particles and shells.
 * This is calculated as both Coulomb and Overlap interactions between
 * Gaussian charge distributions, as well as a short-range repulsive
 * interaction which prevents shells from overlapping in position.
 * The Coulomb interaction is computed using a finite cutoff.
 * @param mix_ mixing mode for computing interaction strengths
 * @param lambdaC_ radial scaling parameter for Coulomb Gaussian width
 * @param lambdaS_ radial scaling parameter for Overlap Gaussian width
 * @param radius_ radii for each type used to determine Gaussian widths
 * @param aOver_ Overlap interaction strength
 * @param aRep_ amplitude of the close-range repulsive interaction
 * @param gammaC_ reduced exponent from Gaussian width for Coulomb interaction
 * @param gammaS_ reduced exponent from Gaussian width for Overalap interaction
 * @param rgammaC_ radial gammaC_
 */
class CalcCGemmCut: public Calculator{
private:
	//global parameters
	double eps_{1.0};//relative permittivity
	Calculator::Mix mix_{Mix::NONE};
    double lambdaC_{1.0};//radial scaling factor - Coulomb
    double lambdaS_{1.0};//radial scaling factor - Overlap
	double rRep_{0.0};//repulsive radius
	//type parameters
    Eigen::VectorXd radius_;//Gaussian radius
    Eigen::MatrixXd aOver_;//overlap amplitude
	Eigen::MatrixXd aRep_;//repulsive amplitude
    Eigen::MatrixXd gammaC_;//reduced expontent - Coulomb
    Eigen::MatrixXd gammaS_;//reduced expontent - Overlap
	Eigen::MatrixXd rgammaC_;//reduced expontent - Coulomb - radical
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
	double& eps(){return eps_;}
    const double& eps()const{return eps_;}
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