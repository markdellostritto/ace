#pragma once
#ifndef CALC_HPP
#define CALC_HPP

// c++
#include <iostream>
// string
#include "str/token.hpp"
// structure
#include "struc/structure.hpp"
#include "struc/neighbor.hpp"

#ifndef CALC_PRINT_FUNC
#define CALC_PRINT_FUNC 0
#endif

class Calculator{
public:
	//name
	class Name{
	public:
		enum Type{
			LJ_CUT,
			COUL_CUT,
			COUL_LONG,
			CGEM_CUT,
			CGEM_LONG,
			CGEMM_CUT,
			CGEMM_LONG,
			NONE
		};
		//constructor
		Name():t_(Type::NONE){}
		Name(Type t):t_(t){}
		//operators
		operator Type()const{return t_;}
		friend std::ostream& operator<<(std::ostream& out, const Name& t);
		//member functions
		static Name read(const char* str);
		static const char* name(const Name& t);
	private:
		Type t_;
		//prevent automatic conversion for other built-in types
		//template<typename T> operator T() const;
	};
	//mix
	class Mix{
	public:
		enum Type{
			ARITHMETIC,
			GEOMETRIC,
			HARMONIC,
			NONE
		};
		//constructor
		Mix():t_(Type::NONE){}
		Mix(Type t):t_(t){}
		//operators
		operator Type()const{return t_;}
		friend std::ostream& operator<<(std::ostream& out, const Mix& t);
		//member functions
		static Mix read(const char* str);
		static const char* name(const Mix& t);
	private:
		Type t_;
		//prevent automatic conversion for other built-in types
		//template<typename T> operator T() const;
	};
protected:
	Name name_;//name
	Mix mix_;//name
    int ntypes_;//number of types
    double rc_;//cutoff radius
	double rc2_;//cutoff radius squared
public:
    //==== contructors/destructors ====
	Calculator():name_(Name::NONE),rc_(0.0),rc2_(0.0),ntypes_(0){}
	Calculator(Name name):name_(name),rc_(0.0),rc2_(0.0),ntypes_(0){}
	Calculator(Name name, double rc);
	virtual ~Calculator(){}

	//==== operators ====
	friend std::ostream& operator<<(std::ostream& out, const Calculator& pot);
	
	//==== access ====
	const Name& name()const{return name_;}
	const double& rc()const{return rc_;}
	const double& rc2()const{return rc2_;}
	const int& ntypes()const{return ntypes_;}
	
	//==== member functions ====
    virtual void resize(int ntypes){ntypes_=ntypes;}
	virtual void init(){}
	virtual void init(const Structure& struc){}
	void read(Token& token);
	virtual void coeff(Token& token){}
	virtual double energy(Structure& struc, const NeighborList& nlist)const{}
	virtual double energy(Structure& struc)const{}
	virtual double compute(Structure& struc, const NeighborList& nlist)const{}
	virtual double compute(Structure& struc)const{}
};

//**********************************************
// serialization
//**********************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const Calculator& obj);
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const Calculator& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(Calculator& obj, const char* arr);
	
}

#endif
