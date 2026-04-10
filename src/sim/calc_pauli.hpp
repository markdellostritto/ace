#pragma once
#ifndef CALC_PAULI_HPP
#define CALC_PAULI_HPP

// c++
#include <iostream>
// eigen
#include <Eigen/Dense>
// mem
#include "mem/serialize.hpp"
// sim
#include "sim/calc.hpp"

#ifndef CALC_PAULI_PRINT_FUNC
#define CALC_PAULI_PRINT_FUNC 0
#endif

/**
 * class CalcPauli
 * @param radius_ radii for each type used to determine Gaussian widths
 * @param gamma_ harmonic average of alpha
 * @param amp_ amplitude of the repulsion
 */
class CalcPauli: public Calculator{
private:
	//global parameters
	double eps_{1.0};//relative permittivity
	//type parameters
    Eigen::VectorXd z_;//Pauli charge
    Eigen::VectorXd radius_;//Gaussian radius
    Eigen::MatrixXd gamma_;//reduced expontent 
    Eigen::MatrixXd amp_;//mixed amplitude
public:
    //==== contructors/destructors ====
	CalcPauli():Calculator(Calculator::Name::PAULI){}
    CalcPauli(double rc):Calculator(Calculator::Name::PAULI,rc){}
	~CalcPauli(){}

    //==== operators ====
    friend std::ostream& operator<<(std::ostream& out, const CalcPauli& calc);

	//==== access ====
	//global parameters
	double& eps(){return eps_;}
    const double& eps()const{return eps_;}
	//type parameters
	Eigen::VectorXd& z(){return z_;}
	const Eigen::VectorXd& z()const{return z_;}
	Eigen::VectorXd& radius(){return radius_;}
	const Eigen::VectorXd& radius()const{return radius_;}
	Eigen::MatrixXd& gamma(){return gamma_;}
	const Eigen::MatrixXd& gamma()const{return gamma_;}
	
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
	
	template <> int nbytes(const CalcPauli& obj);
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcPauli& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcPauli& obj, const char* arr);
	
}

#endif