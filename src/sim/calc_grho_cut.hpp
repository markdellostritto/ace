#pragma once
#ifndef CALC_GRHO_CUT_HPP
#define CALC_GRHO_CUT_HPP

// c++
#include <iostream>
// eigen
#include <Eigen/Dense>
// mem
#include "mem/serialize.hpp"
// sim
#include "sim/calc.hpp"

#ifndef CALC_GRHO_CUT_PRINT_FUNC
#define CALC_GRHO_CUT_PRINT_FUNC 0
#endif

/**
 * class CalcGRhoCut
 * Calculation of the energy of a set of particles and shells.
 * This is calculated as both Coulomb and Overlap interactions between
 * Gaussian charge distributions, as well as a short-range repulsive
 * interaction which prevents shells from overlapping in position.
 * The Coulomb interaction is computed using a finite cutoff.
 * @param mix_ mixing mode for computing interaction strengths
 * @param radius_ radii for each type used to determine Gaussian widths
 * @param gamma_ reduced exponent from Gaussian width for Coulomb interaction
 * @param rgamma_ radial gamma_
 */
class CalcGRhoCut: public Calculator{
private:
	//global parameters
	double eps_{1.0};//relative permittivity
	Calculator::Mix mix_{Mix::NONE};
    //type parameters
    Eigen::VectorXd radius_;//Gaussian radius
    Eigen::MatrixXd gamma_;//reduced expontent 
    Eigen::MatrixXd rgamma_;//reduced expontent - radial
public:
    //==== contructors/destructors ====
	CalcGRhoCut():Calculator(Calculator::Name::GRHO_CUT){}
    CalcGRhoCut(double rc):Calculator(Calculator::Name::GRHO_CUT,rc){}
	CalcGRhoCut(double rc, Calculator::Mix mix);
    ~CalcGRhoCut(){}

    //==== operators ====
    friend std::ostream& operator<<(std::ostream& out, const CalcGRhoCut& calc);

	//==== access ====
	//global parameters
	double& eps(){return eps_;}
    const double& eps()const{return eps_;}
	Calculator::Mix& mix(){return mix_;}
	const Calculator::Mix& mix()const{return mix_;}
    //type parameters
	Eigen::VectorXd& radius(){return radius_;}
	const Eigen::VectorXd& radius()const{return radius_;}
	Eigen::MatrixXd& gamma(){return gamma_;}
	const Eigen::MatrixXd& gamma()const{return gamma_;}
	
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
	
	template <> int nbytes(const CalcGRhoCut& obj);
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcGRhoCut& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcGRhoCut& obj, const char* arr);
	
}

#endif