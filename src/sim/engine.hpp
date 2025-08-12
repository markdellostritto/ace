#pragma once
#ifndef ENGINE_HPP
#define ENGINE_HPP

// c++
#include <iostream>
#include <memory>
// sim
#include "sim/calc.hpp"
#include "sim/constraint.hpp"
// struc
#include "struc/neighbor.hpp"

#ifndef ENGINE_PRINT_FUNC
#define ENGINE_PRINT_FUNC 0
#endif

//****************************************************************************
// Engine
//****************************************************************************

class Engine{
private:
	int stride_;
	int ntypes_;
	double rcmax_;
    NeighborList nlist_;
	std::vector<std::shared_ptr<Calculator> > calcs_;//calculators
	std::vector<std::shared_ptr<Constraint> > constraints_;//constraints
public:
	//==== constructors/destructors ====
	Engine():ntypes_(-1){}
	~Engine(){}
	
	//==== operators ====
	friend std::ostream& operator<<(std::ostream& out, const Engine& engine);
	
	//==== access ====
	//types
	int& stride(){return stride_;}
	const int& stride()const{return stride_;}
	const int ntypes()const{return ntypes_;}
	const double& rcmax()const{return rcmax_;}
	NeighborList& nlist(){return nlist_;}
	const NeighborList& nlist()const{return nlist_;}
	//caclulators
	std::vector<std::shared_ptr<Calculator> >& calcs(){return calcs_;}
	const std::vector<std::shared_ptr<Calculator> >& calcs()const{return calcs_;}
	std::shared_ptr<Calculator>& calc(int i){return calcs_[i];}
	const std::shared_ptr<Calculator>& calc(int i)const{return calcs_[i];}
	//constraints
	std::vector<std::shared_ptr<Constraint> >& constraints(){return constraints_;}
	const std::vector<std::shared_ptr<Constraint> >& constraints()const{return constraints_;}
	std::shared_ptr<Constraint>& constraints(int i){return constraints_[i];}
	const std::shared_ptr<Constraint>& constraints(int i)const{return constraints_[i];}
	
	//==== member functions ====
	//reading/writing
	void read(Token& token);
	//setup/initialization
	void clear();
	void resize(int ntypes);
	void init();
	void init(const Structure& struc);
	//energy/forces
	double energy(Structure& struc);
	double compute(Structure& struc);
	
	//==== static functions ====
	static double ke(const Structure& struc);
};

//**********************************************
// serialization
//**********************************************
	
namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const Engine& obj);
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const Engine& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(Engine& obj, const char* arr);
	
}

#endif