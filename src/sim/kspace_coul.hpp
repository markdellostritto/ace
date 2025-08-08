#pragma once
#ifndef KSPACE_COUL_HPP
#define KSPACE_COUL_HPP

//mem
#include "mem/serialize.hpp"
// structure
#include "struc/structure.hpp"
// sim
#include "sim/kspace.hpp"

#ifndef KSPACEC_PRINT_FUNC
#define KSPACEC_PRINT_FUNC 0
#endif

#ifndef KSPACEC_PRINT_STATUS
#define KSPACEC_PRINT_STATUS 0
#endif

#ifndef KSPACEC_PRINT_DATA
#define KSPACEC_PRINT_DATA 0
#endif

namespace KSpace{

class Coul: public Base{
private:
	bool econst_;
	double eps_;
	double q2_;//sum of squares of charges
	double vc_;//constant term
	mutable std::vector<double> c_;
	mutable std::vector<double> s_;
public:
	//==== constructors/destructors ====
	Coul():eps_(1.0),econst_(true){}
	virtual ~Coul(){}
	
	//==== operators ====
	friend std::ostream& operator<<(std::ostream& out, const Coul& c);
	
	//==== access ====
	bool& econst(){return econst_;}
	const bool& econst()const{return econst_;}
	const double& vc()const{return vc_;}
	double& eps(){return eps_;}
	const double& eps()const{return eps_;}
	
	//==== member functions ====
	void init(const Structure& struc);
	double energy(Structure& struc)const;
	double compute(Structure& struc)const;
	double compute_explicit(Structure& struc)const;
	double compute_brute(Structure& struc)const;
	
};

}

//**********************************************
// serialization
//**********************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const KSpace::Coul& obj);
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const KSpace::Coul& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(KSpace::Coul& obj, const char* arr);
	
}

#endif