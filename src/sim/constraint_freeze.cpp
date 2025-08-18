// ace
#include "sim/constraint_freeze.hpp"

//==== operator ====

std::ostream& operator<<(std::ostream& out, const ConstraintFreeze& constraint){
	return out<<static_cast<const Constraint&>(constraint);
}

//==== member functions ====

void ConstraintFreeze::read(Token& token){
    static_cast<Constraint&>(*this).read(token);
}

double ConstraintFreeze::compute(Structure& struc, const NeighborList& nlist){
    for(int i=0; i<indices_.size(); ++i){
        struc.force(indices_[i]).setZero();
    }
    return 0.0;
}

double ConstraintFreeze::compute(Structure& struc){
    for(int i=0; i<indices_.size(); ++i){
        struc.force(indices_[i]).setZero();
    }
    return 0.0;
}