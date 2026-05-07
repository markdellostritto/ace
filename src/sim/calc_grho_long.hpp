#pragma once
#ifndef CALC_GRHO_LONG_HPP
#define CALC_GRHO_LONG_HPP

// c++
#include <iostream>
// eigen
#include <Eigen/Dense>
// mem
#include "mem/serialize.hpp"
// sim
#include "sim/calc.hpp"
#include "sim/kspace_coul.hpp"

#ifndef CALC_GRHO_LONG_PRINT_FUNC
#define CALC_GRHO_LONG_PRINT_FUNC 0
#endif

/**
 * class CalcGRhoLong
 * @param eps_ relative permittivity
 * @param prec_ ewald precision
 * @param coul_ Ewald Coulomb calculator
 * @param radius_ radii for each type used to determine Gaussian widths
 * @param gamma_ reduced exponent from Gaussian width for Coulomb interaction
 * @param rgamma_ radial gamma_
 */
class CalcGRhoLong: public Calculator{
private:
	//k-space
	double prec_{1.0e-12};//ewald precision
	KSpace::Coul coul_;//Ewald Coulomb calculator
	//global parameters
	double eps_{1.0};//relative permittivity
	//type parameters
    Eigen::VectorXd radius_;//radius
    Eigen::MatrixXd gamma_;//reduced expontent
    Eigen::MatrixXd rgamma_;//reduced expontent - radical
public:
    //==== contructors/destructors ====
	CalcGRhoLong():Calculator(Calculator::Name::GRHO_LONG){}
    CalcGRhoLong(double rc):Calculator(Calculator::Name::GRHO_LONG,rc){}
	~CalcGRhoLong(){}

    //==== operators ====
    friend std::ostream& operator<<(std::ostream& out, const CalcGRhoLong& calc);

	//==== access ====
	//k-space
	double& prec(){return prec_;}
    const double& prec()const{return prec_;}
	KSpace::Coul& coul(){return coul_;}
	const KSpace::Coul& coul()const{return coul_;}
	//global parameters
	double& eps(){return eps_;}
    const double& eps()const{return eps_;}
	//type parameters
	Eigen::VectorXd& radius(){return radius_;}
	const Eigen::VectorXd& radius()const{return radius_;}
	Eigen::MatrixXd& gammaC(){return gamma_;}
	const Eigen::MatrixXd& gammaC()const{return gamma_;}
	Eigen::MatrixXd& rgammaC(){return rgamma_;}
	const Eigen::MatrixXd& rgammaC()const{return rgamma_;}
	
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
	
	template <> int nbytes(const CalcGRhoLong& obj);
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcGRhoLong& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcGRhoLong& obj, const char* arr);
	
}

#endif