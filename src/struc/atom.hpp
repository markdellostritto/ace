#pragma once
#ifndef ATOM_HPP
#define ATOM_HPP

// c++ libraries
#include <iosfwd>
// mem
#include "mem/serialize.hpp"
// str
#include "str/token.hpp"

//**********************************************************************************************
//Atom
//**********************************************************************************************

struct Atom{
	//==== data ====
	//coordinates
	bool frac;
	//basic properties
	bool name;
	bool an;
	bool type;
	//serial properties
	bool mass;
	bool charge;
	bool radius;
	bool eta;
	//vector properties
	bool image;
	bool posn;
	bool vel;
	bool force;
	//nnp
	bool symm;
	//==== constructors/destructors ====
	Atom(){defaults();}
	~Atom(){}
	//==== operators ====
	friend std::ostream& operator<<(std::ostream& out, const Atom& atom);
	//==== member functions ====
	void defaults();
	void clear(){defaults();}
	static Atom read(Token& token);
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