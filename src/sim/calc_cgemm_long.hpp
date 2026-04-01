#pragma once
#ifndef CALC_CGEMM_LONG_HPP
#define CALC_CGEMM_LONG_HPP

// c++
#include <iostream>
// eigen
#include <Eigen/Dense>
// mem
#include "mem/serialize.hpp"
// sim
#include "sim/calc.hpp"
#include "sim/kspace_coul.hpp"

#ifndef CALC_CGEMM_LONG_PRINT_FUNC
#define CALC_CGEMM_LONG_PRINT_FUNC 0
#endif

/**
 * class CalcCGemmLong
 * Calculation of the energy of a set of particles and shells.
 * This is calculated as both Coulomb and Overlap interactions between
 * Gaussian charge distributions, as well as a short-range repulsive
 * interaction which prevents shells from overlapping in position.
 * The Coulomb interaction is computed using Ewald summation to
 * take into account long-range interactions.
 * @param eps_ relative permittivity
 * @param prec_ ewald precision
 * @param coul_ Ewald Coulomb calculator
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
class CalcCGemmLong: public Calculator{
private:
	//k-space
	double prec_{1.0e-12};//ewald precision
	KSpace::Coul coul_;//Ewald Coulomb calculator
	//global parameters
	double eps_{1.0};//relative permittivity
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
	CalcCGemmLong():Calculator(Calculator::Name::CGEMM_LONG){}
    CalcCGemmLong(double rc):Calculator(Calculator::Name::CGEMM_LONG,rc){}
	CalcCGemmLong(double rc, double lambdaC, double lambdaS, double rRep, Calculator::Mix mix);
    ~CalcCGemmLong(){}

    //==== operators ====
    friend std::ostream& operator<<(std::ostream& out, const CalcCGemmLong& calc);

	//==== access ====
	//k-space
	double& prec(){return prec_;}
    const double& prec()const{return prec_;}
	KSpace::Coul& coul(){return coul_;}
	const KSpace::Coul& coul()const{return coul_;}
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
	Eigen::MatrixXd& rgammaC(){return rgammaC_;}
	const Eigen::MatrixXd& rgammaC()const{return rgammaC_;}
	
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
	
	template <> int nbytes(const CalcCGemmLong& obj);
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcCGemmLong& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcCGemmLong& obj, const char* arr);
	
}

#endif