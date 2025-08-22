// c
#include <cstdio>
// c++
#include <iostream>
#include <stdexcept>
// str
#include "str/string.hpp"
#include "str/print.hpp"
#include "str/token.hpp"
// math
#include "math/const.hpp"
// sim
#include "sim/engine.hpp"
#include "sim/calc.hpp"
#include "sim/calc_factory.hpp"
#include "sim/constraint_factory.hpp"

//****************************************************************************
// Engine
//****************************************************************************

//==== operators ====

std::ostream& operator<<(std::ostream& out, const Engine& engine){
	char* str=new char[print::len_buf];
	out<<print::buf(str)<<"\n";
	out<<print::title("ENGINE",str)<<"\n";
	out<<"ntypes = "<<engine.ntypes()<<"\n";
	out<<"rcmax  = "<<engine.rcmax()<<"\n";
	out<<"calculators = \n";
	for(int i=0; i<engine.calcs().size(); ++i){
		out<<engine.calc(i)<<"\n";
	}
	out<<"constraints = \n";
	for(int i=0; i<engine.constraints().size(); ++i){
		out<<engine.constraint(i)<<"\n";
	}
	out<<print::buf(str);
	delete[] str;
	return out;
}

//==== member functions ====

//** setup/initialization **

void Engine::clear(){
	if(ENGINE_PRINT_FUNC>0) std::cout<<"Engine::clear():\n";
	ntypes_=-1;
	rcmax_=0;
	calcs_.clear();
	constraints_.clear();
}

void Engine::resize(int ntypes){
	if(ENGINE_PRINT_FUNC>0) std::cout<<"Engine::resize(int):\n";
	if(ntypes<=0) throw std::invalid_argument("Invalid number of types.");
	ntypes_=ntypes;
	for(int i=0; i<calcs_.size(); ++i){
		calcs_[i]->resize(ntypes_);
	}
}

void Engine::init(){
	if(ENGINE_PRINT_FUNC>0) std::cout<<"Engine::init():\n";
	rcmax_=0;
	for(int i=0; i<calcs_.size(); ++i){
		if(calcs_[i]->rc()>rcmax_) rcmax_=calcs_[i]->rc();
		calcs_[i]->init();
	}
	if(rcmax_==0) throw std::invalid_argument("Engine::init(): Invalid max cutoff.\n");
}

void Engine::init(const Structure& struc){
	for(int i=0; i<calcs_.size(); ++i){
		calcs_[i]->init(struc);
	}
}

//** energy/forces **

double Engine::energy(Structure& struc)const{
	if(ENGINE_PRINT_FUNC>0) std::cout<<"Engine::energy(Structure&):\n";
	double energy=0;
	//reset energy
	struc.pe()=0;
	//compute energy
	for(int i=0; i<calcs_.size(); ++i){
		energy+=calcs_[i]->energy(struc);
	}
	//return energy
	return energy;
}

double Engine::energy(Structure& struc, const NeighborList& nlist)const{
	if(ENGINE_PRINT_FUNC>0) std::cout<<"Engine::energy(Structure&,const NeighborList&):\n";
	double energy=0;
	//reset energy
	struc.pe()=0;
	//compute energy
	for(int i=0; i<calcs_.size(); ++i){
		energy+=calcs_[i]->energy(struc,nlist);
	}
	//return energy
	return energy;
}

double Engine::compute(Structure& struc)const{
	if(ENGINE_PRINT_FUNC>0) std::cout<<"Engine::compute(Structure&):\n";
	double energy=0;
	//reset energy/forces
	struc.pe()=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.force(i).setZero();
	}
	//compute energy/forces
	for(int i=0; i<calcs_.size(); ++i){
		energy+=calcs_[i]->compute(struc);
	}
	//compute constraints
	for(int i=0; i<constraints_.size(); ++i){
		constraints_[i]->compute(struc);
	}
	//return energy
	return energy;
}

double Engine::compute(Structure& struc, const NeighborList& nlist)const{
	if(ENGINE_PRINT_FUNC>0) std::cout<<"Engine::compute(Structure&,const NeighborList&):\n";
	double energy=0;
	//reset energy/forces
	struc.pe()=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.force(i).setZero();
	}
	//compute energy/forces
	for(int i=0; i<calcs_.size(); ++i){
		energy+=calcs_[i]->compute(struc,nlist);
	}
	//compute constraints
	for(int i=0; i<constraints_.size(); ++i){
		constraints_[i]->compute(struc,nlist);
	}
	//return energy
	return energy;
}

void Engine::read(Token& token){
	if(ENGINE_PRINT_FUNC>0) std::cout<<"Engine::read(Token&):\n";
	const int stride=std::atoi(token.next().c_str());
	if(stride<=0) throw std::invalid_argument("Engine::read(Token&): invalid neighbor stride.");
	stride_=stride;
}

double Engine::ke(const Structure& struc){
	double energy=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		energy+=struc.vel(i).squaredNorm()*struc.mass(i);
	}
	return 0.5*energy;
}

//**********************************************
// serialization
//**********************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const Engine& obj){
		if(ENGINE_PRINT_FUNC>0) std::cout<<"nbytes(const Engine&):\n";
		int size=0;
		size+=sizeof(int);//stride_
		size+=sizeof(int);//ntypes_
		size+=sizeof(int);//ncalcs
		size+=sizeof(int);//nconstraints
		for(int i=0; i<obj.calcs().size(); ++i){
			size+=nbytes(obj.calc(i));
		}
		for(int i=0; i<obj.constraints().size(); ++i){
			size+=nbytes(obj.constraint(i));
		}
		return size;
	}
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const Engine& obj, char* arr){
		if(ENGINE_PRINT_FUNC>0) std::cout<<"pack(const Engine&,char*):\n";
		int pos=0;
		const int stride=obj.stride();
		const int ntypes=obj.ntypes();
		const int ncalcs=obj.calcs().size();
		const int nconstraints=obj.constraints().size();
		std::memcpy(arr+pos,&stride,sizeof(int)); pos+=sizeof(int);//nstride_
		std::memcpy(arr+pos,&ntypes,sizeof(int)); pos+=sizeof(int);//ntypes_
		std::memcpy(arr+pos,&ncalcs,sizeof(int)); pos+=sizeof(int);//ncalcs
		std::memcpy(arr+pos,&nconstraints,sizeof(int)); pos+=sizeof(int);//nconstraints
		for(int i=0; i<obj.calcs().size(); ++i){
			pos+=pack(obj.calc(i),arr+pos);
		}
		for(int i=0; i<obj.constraints().size(); ++i){
			pos+=pack(obj.constraint(i),arr+pos);
		}
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(Engine& obj, const char* arr){
		if(ENGINE_PRINT_FUNC>0) std::cout<<"unpack(Engine&,const char*):\n";
		int pos=0;
		int stride=0;
		int ntypes=0;
		int ncalcs=0;
		int nconstraints=0;
		std::memcpy(&stride,arr+pos,sizeof(int)); pos+=sizeof(int);//nstride_
		std::memcpy(&ntypes,arr+pos,sizeof(int)); pos+=sizeof(int);//ntypes_
		std::memcpy(&ncalcs,arr+pos,sizeof(int)); pos+=sizeof(int);//ncalcs
		std::memcpy(&nconstraints,arr+pos,sizeof(int)); pos+=sizeof(int);//ncalcs
		obj.resize(ntypes);
		obj.stride()=stride;
		obj.calcs().resize(ncalcs);
		obj.constraints().resize(nconstraints);
		for(int i=0; i<obj.calcs().size(); ++i){
			pos+=unpack(obj.calc(i),arr+pos);
		}
		for(int i=0; i<obj.constraints().size(); ++i){
			pos+=unpack(obj.constraint(i),arr+pos);
		}
		return pos;
	}
	
}
