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
	out<<"nlist  = "<<engine.nlist()<<"\n";
	out<<"calculators = \n";
	for(int i=0; i<engine.calcs().size(); ++i){
		out<<engine.calc(i)<<"\n";
	}
	out<<"constraints = \n";
	for(int i=0; i<engine.constraints().size(); ++i){
		out<<engine.constraints(i)<<"\n";
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
	nlist_.clear();
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
	nlist_=NeighborList(rcmax_);
}

//** energy/forces **

double Engine::compute(Structure& struc){
	if(ENGINE_PRINT_FUNC>0) std::cout<<"Engine::compute(const Structure&):\n";
	double energy=0;
	//reset forces
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.force(i).setZero();
	}
	//compute energy/forces
	for(int i=0; i<calcs_.size(); ++i){
		energy+=calcs_[i]->compute(struc,nlist_);
	}
	//compute constraints
	for(int i=0; i<constraints_.size(); ++i){
		constraints_[i]->compute(struc,nlist_);
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
		size+=nbytes(obj.nlist());
		size+=sizeof(int);//ncalcs
		/*for(int i=0; i<obj.calcs().size(); ++i){
			size+=sizeof(int);//name_
			switch(obj.calc(i)->name()){
				case Calculator::Name::LJ_CUT:{
					size+=nbytes(static_cast<const CalcLJCut&>(*obj.calc(i)));
				}break;
				default: break;
			}
		}*/
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
		std::memcpy(arr+pos,&stride,sizeof(int)); pos+=sizeof(int);//nstride_
		std::memcpy(arr+pos,&ntypes,sizeof(int)); pos+=sizeof(int);//ntypes_
		pos+=pack(obj.nlist(),arr+pos);
		int ncalcs=obj.calcs().size();
		std::memcpy(arr+pos,&ncalcs,sizeof(int)); pos+=sizeof(int);//job_
		/*for(int i=0; i<obj.calcs().size(); ++i){
			Calculator::Name name=obj.calc(i)->name();
			std::memcpy(arr+pos,&name,sizeof(int)); pos+=sizeof(int);//name_
			switch(obj.calc(i)->name()){
				case Calculator::Name::LJ_CUT:{
					pos+=pack(static_cast<const CalcLJCut&>(*obj.calc(i)),arr+pos);
				}break;
				default: break;
			}
		}*/
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
		std::memcpy(&stride,arr+pos,sizeof(int)); pos+=sizeof(int);//nstride_
		std::memcpy(&ntypes,arr+pos,sizeof(int)); pos+=sizeof(int);//ntypes_
		obj.resize(ntypes);
		obj.stride()=stride;
		pos+=unpack(obj.nlist(),arr+pos);
		int ncalcs=0;
		std::memcpy(&ncalcs,arr+pos,sizeof(int)); pos+=sizeof(int);//ntypes_
		/*obj.calcs().resize(ncalcs);
		for(int i=0; i<obj.calcs().size(); ++i){
			Calculator::Name name;
			std::memcpy(&name,arr+pos,sizeof(int)); pos+=sizeof(int);//name_
			switch(name){
				case Calculator::Name::LJ_CUT:{
					CalcLJCut calc; pos+=unpack(calc,arr+pos); obj.calc(i).reset(calc);
				}break;
				default: break;
			}
		}*/
		return pos;
	}
	
}
