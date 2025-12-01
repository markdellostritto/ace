#pragma once
#ifndef CALC_THOLE_LONG_HPP
#define CALC_THOLE_LONG_HPP

// c++
#include <iostream>
// eigen
#include <Eigen/Dense>
// mem
#include "mem/serialize.hpp"
// sim
#include "sim/calc.hpp"
#include "sim/kspace_dipole.hpp"

#ifndef CALC_THOLE_LONG_PRINT_FUNC
#define CALC_THOLE_LONG_PRINT_FUNC 0
#endif

class CalcTholeLong: public Calculator{
public:
    class IDD{
	public:
		enum Type{
            IDEAL,
            EXP,
            ERF,
            NONE
		};
		//constructor
		IDD():t_(Type::NONE){}
		IDD(Type t):t_(t){}
		//operators
		operator Type()const{return t_;}
		//member functions
		static IDD read(const char* str);
		static const char* name(const IDD& name);
	private:
		Type t_;
		//prevent automatic conversion for other built-in types
		//template<typename T> operator T() const;
	};
private:
    IDD idd_{IDD::NONE};
	double a_{1.0};
	double eps_{1.0};
	double prec_{1.0e-6};
    KSpace::Dipole kDipole_;
    Eigen::VectorXd alpha_;//polarizability
public:
	//==== constructors/destructors ====
	CalcTholeLong():Calculator(Calculator::Name::THOLE_LONG){}
    CalcTholeLong(double rc):Calculator(Calculator::Name::THOLE_LONG,rc){}
	~CalcTholeLong(){}
	
	//==== operator ====
	friend std::ostream& operator<<(std::ostream& out, const CalcTholeLong& calc);
	
	//==== access ====
    IDD& idd(){return idd_;}
    const IDD& idd()const{return idd_;}
	double& a(){return a_;}
    const double& a()const{return a_;}
	double& prec(){return prec_;}
	const double& prec()const{return prec_;}
	double& eps(){return eps_;}
	const double& eps()const{return eps_;}
    KSpace::Dipole& kDipole(){return kDipole_;}
	const KSpace::Dipole& kDipole()const{return kDipole_;}
	Eigen::VectorXd& alpha(){return alpha_;}
	const Eigen::VectorXd& alpha()const{return alpha_;}
	
	//==== member functions ====
    void resize(int ntypes);
	void init();
	void init(const Structure& struc);
	void read(Token& token);
	void coeff(Token& token);
	double energy(Structure& struc, const NeighborList& nlist)const;
	double compute(Structure& struc, const NeighborList& nlist)const;
	Eigen::Matrix3d& interMat(const Eigen::Vector3d& drv, Eigen::Matrix3d& mat, double a);
};

//**********************************************
// serialization
//**********************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const CalcTholeLong& obj);
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcTholeLong& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcTholeLong& obj, const char* arr);
	
}

#endif