#pragma once
#ifndef ATOM_HPP
#define ATOM_HPP

// c++ libraries
#include <iosfwd>
#include <array>
// mem
#include "mem/serialize.hpp"
// str
#include "str/token.hpp"

//**********************************************************************************************
//Atom
//**********************************************************************************************

class Atom{
public:
	static const int SIZE=14;
private:
	std::array<bool,SIZE> data_;
public:
	
	//==== constructors/destructors ====
	Atom(){clear();}
	~Atom(){}

	//==== operators ====
	friend std::ostream& operator<<(std::ostream& out, const Atom& atom);

	//==== access - mutable ====
	//data
	std::array<bool,14>& data(){return data_;}
	//coordinates
	bool& frac(){return data_[0];} 
	//basic properties
	bool& name(){return data_[1];} 
	bool& an(){return data_[2];} 
	bool& type(){return data_[3];} 
	//serial properties
	bool& mass(){return data_[4];} 
	bool& charge(){return data_[5];} 
	bool& radius(){return data_[6];} 
	bool& eta(){return data_[7];} 
	//vector properties
	bool& image(){return data_[8];} 
	bool& posn(){return data_[9];} 
	bool& vel(){return data_[10];} 
	bool& force(){return data_[11];} 
	bool& dipole(){return data_[12];} 
	//nnp
	bool& symm(){return data_[13];} 

	//==== access - const ====
	//data
	const std::array<bool,14>& data()const{return data_;}
	//coordinates
	const bool& frac()const{return data_[0];} 
	//basic properties
	const bool& name()const{return data_[1];} 
	const bool& an()const{return data_[2];} 
	const bool& type()const{return data_[3];} 
	//serial properties
	const bool& mass()const{return data_[4];} 
	const bool& charge()const{return data_[5];} 
	const bool& radius()const{return data_[6];} 
	const bool& eta()const{return data_[7];} 
	//vector properties
	const bool& image()const{return data_[8];} 
	const bool& posn()const{return data_[9];} 
	const bool& vel()const{return data_[10];} 
	const bool& force()const{return data_[11];} 
	const bool& dipole()const{return data_[12];} 
	//nnp
	const bool& symm()const{return data_[13];} 

	//==== member functions ====
	void clear(){data_.fill(false);}
	Atom& read(Token& token);
};

//**********************************************************************************************
// serialization
//**********************************************************************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const Atom& obj);
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const Atom& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(Atom& obj, const char* arr);
	
}

#endif