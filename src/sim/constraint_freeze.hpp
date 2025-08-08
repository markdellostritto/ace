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
};

#endif