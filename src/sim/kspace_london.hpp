#pragma once
#ifndef KSPACE_LONDON_HPP
#define KSPACE_LONDON_HPP

//mem
#include "mem/serialize.hpp"
// structure
#include "struc/structure.hpp"
// sim
#include "sim/kspace.hpp"

#ifndef KSPACEL_PRINT_FUNC
#define KSPACEL_PRINT_FUNC 0
#endif

#ifndef KSPACEL_PRINT_STATUS
#define KSPACEL_PRINT_STATUS 0
#endif

#ifndef KSPACEL_PRINT_DATA
#define KSPACEL_PRINT_DATA 0
#endif

namespace KSpace{

class London: public Base{
private:
	double a3_{0.0};//(Ewald alpha)^3
	double biis_{0.0};//sum of b(i,i) for all atoms i,i
	double bijs_{0.0};//sum of b(i,j) for all atoms i,j
	double ec_{0.0};//constant energy term
	Eigen::MatrixXd b_;//London c6 coefficients (ntypes,ntypes)
	Eigen::VectorXd br_;//radial diagonal of b_
	mutable std::vector<double> c_;//store cos(k*r)
	mutable std::vector<double> s_;//store sin(k*r)
public:
	//==== constructors/destructors ====
	London(){}
	virtual ~London(){}
	
	//==== operators ====
	friend std::ostream& operator<<(std::ostream& out, const London& l);
	
	//==== access ====
	const double& a3()const{return a3_;}
	Eigen::MatrixXd& b(){return b_;}
	const Eigen::MatrixXd& b()const{return b_;}
	
	//==== static functions ====
	double fv(double pre, double a, double rc, double prec);
	double fd(double pre, double a, double rc, double prec);
	
	//==== member functions ====
	void init(const Structure& struc);
	double energy(Structure& struc)const;
	double compute(Structure& struc)const;
};

}

//**********************************************
// serialization
//**********************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const KSpace::London& obj);
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const KSpace::London& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(KSpace::London& obj, const char* arr);
	
}

#endif