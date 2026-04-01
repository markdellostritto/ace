#pragma once
#ifndef CONSTRAINT_FREEZE_HPP
#define CONSTRAINT_FREEZE_HPP

// ace
#include "sim/constraint.hpp"

class ConstraintFreeze: public Constraint{
private:

public:
    //==== contructors/destructors ====
    ConstraintFreeze():Constraint(Constraint::Name::FREEZE){}

    //==== operators ====
	friend std::ostream& operator<<(std::ostream& out, const ConstraintFreeze& constraint);

	//==== member functions ====
    void read(Token& token);
    double compute(Structure& struc, const NeighborList& nlist);
    double compute(Structure& struc);
};

//**********************************************
// serialization
//**********************************************
	
namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const ConstraintFreeze& obj);
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const ConstraintFreeze& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(ConstraintFreeze& obj, const char* arr);
	
}

#endif