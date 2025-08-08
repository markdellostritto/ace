#pragma once
#ifndef TYPE_HPP
#define TYPE_HPP

//c++ libraries
#include <iosfwd>
#include <string>
// str
#include "str/token.hpp"
#include "str/string.hpp"
// mem
#include "mem/serialize.hpp"

#ifndef TYPE_PRINT_FUNC
#define TYPE_PRINT_FUNC 0
#endif

//************************************************************
// TYPE
//************************************************************

class Type{
private:
	//name
	std::string name_;
	int id_;
	//standard
	double energy_;
	double mass_;
	double charge_;
	double z_;
	//radius
	double rvdw_;
	double rcov_;
public:
	//==== constructors/destructors ====
	Type(){clear();}
	Type(const std::string& name){clear();name_=name;}
	~Type(){}
	
	//==== operators ====
	friend std::ostream& operator<<(std::ostream& out, const Type& type);
	
	//==== access ====
	//name
	std::string& name(){return name_;}
	const std::string& name()const{return name_;}
	int& id(){return id_;}
	const int& id()const{return id_;}
	//standard
	double& energy(){return energy_;}
	const double& energy()const{return energy_;}
	double& mass(){return mass_;}
	const double& mass()const{return mass_;}
	double& charge(){return charge_;}
	const double& charge()const{return charge_;}
	double& z(){return z_;}
	const double& z()const{return z_;}
	//radius
	double& rvdw(){return rvdw_;}
	const double& rvdw()const{return rvdw_;}
	double& rcov(){return rcov_;}
	const double& rcov()const{return rcov_;}
	
	//==== member functions ====
	void clear();
	
	//==== static functions ====
	static Type& read(const char* str, Type& type);
	static Type& read(Type& type, Token& token);
	static void write(FILE* out, const Type& type);
};

namespace serialize{

//**********************************************
// byte measures
//**********************************************

template <> int nbytes(const Type& obj);

//**********************************************
// packing
//**********************************************

template <> int pack(const Type& obj, char* arr);

//**********************************************
// unpacking
//**********************************************

template <> int unpack(Type& obj, const char* arr);

}

#endif
