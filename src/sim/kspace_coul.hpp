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
	double eps_{1.0};
	double qs_{0.0};//sum of charges
	double q2s_{0.0};//sum of squares of charges
	double vc_{0.0};//constant term
	double vq_{0.0};//net charge term
	double ec_{0.0};//constant energy term
	mutable std::vector<double> c_;//store cos(k*r)
	mutable std::vector<double> s_;//store sin(k*r)
public:
	//==== constructors/destructors ====
	Coul(){}
	virtual ~Coul(){}
	
	//==== operators ====
	friend std::ostream& operator<<(std::ostream& out, const Coul& c);
	
	//==== access ====
	double& eps(){return eps_;}
	const double& eps()const{return eps_;}
	const double& qs()const{return qs_;}
	const double& q2s()const{return q2s_;}
	const double& vc()const{return vc_;}
	const double& vq()const{return vq_;}
	const double& ec()const{return ec_;}
	
	//==== member functions ====
	void init(const Structure& struc);
	double energy(Structure& struc)const;
	double compute(Structure& struc)const;
	double compute_explicit(Structure& struc)const;
	
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