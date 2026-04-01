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

//**********************************************
// serialization
//**********************************************
	
namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const ConstraintFreeze& obj){
		int size=0;
        size+=nbytes(static_cast<const Constraint&>(obj));
		return size;
	}
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const ConstraintFreeze& obj, char* arr){
		int pos=0;
		pos+=pack(static_cast<const Constraint&>(obj),arr+pos);
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(ConstraintFreeze& obj, const char* arr){
		int pos=0;
		pos+=unpack(static_cast<Constraint&>(obj),arr+pos);
		return pos;
	}
	
}