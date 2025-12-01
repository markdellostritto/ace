#pragma once
#ifndef CALC_THOLE_CUT_HPP
#define CALC_THOLE_CUT_HPP

// c++
#include <iostream>
// eigen
#include <Eigen/Dense>
// mem
#include "mem/serialize.hpp"
// sim
#include "sim/calc.hpp"

#ifndef CALC_THOLE_CUT_PRINT_FUNC
#define CALC_THOLE_CUT_PRINT_FUNC 0
#endif

class CalcTholeCut: public Calculator{
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
    double a_{1.0};//scaling constant
    IDD idd_{IDD::NONE};//interaction - dipole-dipole
    Eigen::VectorXd alpha_;//polarizability
public:
    //==== contructors/destructors ====
	CalcTholeCut():Calculator(Calculator::Name::THOLE_CUT){}
    CalcTholeCut(double rc):Calculator(Calculator::Name::THOLE_CUT,rc){}
    ~CalcTholeCut(){}

    //==== operators ====
    friend std::ostream& operator<<(std::ostream& out, const CalcTholeCut& calc);
    
    //==== access ====
    double& a(){return a_;}
    const double& a()const{return a_;}
    IDD& idd(){return idd_;}
    const IDD& idd()const{return idd_;}

    //==== member functions ====
    void resize(int ntypes);
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
	
	template <> int nbytes(const CalcTholeCut& obj);
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcTholeCut& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcTholeCut& obj, const char* arr);
	
}

#endif