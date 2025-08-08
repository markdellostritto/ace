#pragma once
#ifndef CALC_CGEM_LONG_HPP
#define CALC_CGEM_LONG_HPP

// c++
#include <iostream>
// eigen
#include <Eigen/Dense>
// mem
#include "mem/serialize.hpp"
// sim
#include "sim/calc.hpp"
#include "sim/kspace_coul.hpp"

#ifndef CALC_CGEM_LONG_PRINT_FUNC
#define CALC_CGEM_LONG_PRINT_FUNC 0
#endif

class CalcCGemLong: public Calculator{
private:
    double eps_;
	double prec_;
    double lambdaC_;//radial scaling factor - Coulomb
    double lambdaS_;//radial scaling factor - Overlap
    Eigen::VectorXd radius_;//radius
    Eigen::MatrixXd aOver_;//overlap amplitude
	Eigen::MatrixXd aRep_;//repulsive amplitude
    Eigen::MatrixXd muC_;//reduced expontent - Coulomb
    Eigen::MatrixXd muS_;//reduced expontent - Overlap
    Eigen::MatrixXd rmuC_;//reduced expontent - Coulomb
    KSpace::Coul coul_;
public:
    //==== contructors/destructors ====
	CalcCGemLong();
    CalcCGemLong(double rc);
	CalcCGemLong(double rc, double lambdaC, double lambdaS);
    ~CalcCGemLong(){}

    //==== operators ====
    friend std::ostream& operator<<(std::ostream& out, const CalcCGemLong& calc);

	//==== access ====
    double& eps(){return eps_;}
    const double& eps()const{return eps_;}
	double& prec(){return prec_;}
    const double& prec()const{return prec_;}
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
	Eigen::MatrixXd& rmuC(){return rmuC_;}
	const Eigen::MatrixXd& rmuC()const{return rmuC_;}
	Eigen::MatrixXd& muS(){return muS_;}
	const Eigen::MatrixXd& muS()const{return muS_;}
	KSpace::Coul& coul(){return coul_;}
	const KSpace::Coul& coul()const{return coul_;}
	
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
	
	template <> int nbytes(const CalcCGemLong& obj);
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcCGemLong& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcCGemLong& obj, const char* arr);
	
}

#endif