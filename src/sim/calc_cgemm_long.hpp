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

class CalcCGemmLong: public Calculator{
public:
	static const double rRep;//repulsive radius
private:
	//k-space
	double eps_;
	double prec_;
	KSpace::Coul coul_;
	//global parameters
	Calculator::Mix mix_;
    double lambdaC_;//radial scaling factor - Coulomb
    double lambdaS_;//radial scaling factor - Overlap
	//type parameters
    Eigen::VectorXd radius_;//radius
    Eigen::MatrixXd aOver_;//overlap amplitude
	Eigen::MatrixXd aRep_;//repulsive amplitude
    Eigen::MatrixXd gammaC_;//reduced expontent - Coulomb
    Eigen::MatrixXd gammaS_;//reduced expontent - Overlap
    Eigen::MatrixXd rgammaC_;//reduced expontent - Coulomb
public:
    //==== contructors/destructors ====
	CalcCGemmLong();
    CalcCGemmLong(double rc);
	CalcCGemmLong(double rc, double lambdaC, double lambdaS);
    ~CalcCGemmLong(){}

    //==== operators ====
    friend std::ostream& operator<<(std::ostream& out, const CalcCGemmLong& calc);

	//==== access ====
	//k-space
	double& eps(){return eps_;}
    const double& eps()const{return eps_;}
	double& prec(){return prec_;}
    const double& prec()const{return prec_;}
	KSpace::Coul& coul(){return coul_;}
	const KSpace::Coul& coul()const{return coul_;}
	//global parameters
	Calculator::Mix& mix(){return mix_;}
	const Calculator::Mix& mix()const{return mix_;}
    double& lambdaC(){return lambdaC_;}
    const double& lambdaC()const{return lambdaC_;}
	double& lambdaS(){return lambdaS_;}
    const double& lambdaS()const{return lambdaS_;}
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