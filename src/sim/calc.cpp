// c
#include <cstring>
// c++
#include <iostream>
// ace
#include "sim/calc.hpp"

//==== Name ====

Calculator::Name Calculator::Name::read(const char* str){
	if(std::strcmp(str,"LJ_CUT")==0) return Calculator::Name::LJ_CUT;
	else if(std::strcmp(str,"COUL_CUT")==0) return Calculator::Name::COUL_CUT;
	else if(std::strcmp(str,"COUL_LONG")==0) return Calculator::Name::COUL_LONG;
	else if(std::strcmp(str,"DIPOLE_CUT")==0) return Calculator::Name::DIPOLE_CUT;
	else if(std::strcmp(str,"DIPOLE_LONG")==0) return Calculator::Name::DIPOLE_LONG;
	else if(std::strcmp(str,"THOLE_CUT")==0) return Calculator::Name::THOLE_CUT;
	else if(std::strcmp(str,"THOLE_LONG")==0) return Calculator::Name::THOLE_LONG;
	else if(std::strcmp(str,"CGEMM_CUT")==0) return Calculator::Name::CGEMM_CUT;
	else if(std::strcmp(str,"CGEMM_LONG")==0) return Calculator::Name::CGEMM_LONG;
	else return Calculator::Name::NONE;
}

const char* Calculator::Name::name(const Calculator::Name& t){
	switch(t){
		case Calculator::Name::LJ_CUT: return "LJ_CUT";
		case Calculator::Name::COUL_CUT: return "COUL_CUT";
		case Calculator::Name::COUL_LONG: return "COUL_LONG";
		case Calculator::Name::DIPOLE_CUT: return "DIPOLE_CUT";
		case Calculator::Name::DIPOLE_LONG: return "DIPOLE_LONG";
		case Calculator::Name::THOLE_CUT: return "THOLE_CUT";
		case Calculator::Name::THOLE_LONG: return "THOLE_LONG";
		case Calculator::Name::CGEMM_CUT: return "CGEMM_CUT";
		case Calculator::Name::CGEMM_LONG: return "CGEMM_LONG";
		default: return "NONE";
	}
}

std::ostream& operator<<(std::ostream& out, const Calculator::Name& t){
	switch(t){
		case Calculator::Name::LJ_CUT: out<<"LJ_CUT"; break;
		case Calculator::Name::COUL_CUT: out<<"COUL_CUT"; break;
		case Calculator::Name::COUL_LONG: out<<"COUL_LONG"; break;
		case Calculator::Name::DIPOLE_CUT: out<<"DIPOLE_CUT"; break;
		case Calculator::Name::DIPOLE_LONG: out<<"DIPOLE_LONG"; break;
		case Calculator::Name::THOLE_CUT: out<<"THOLE_CUT"; break;
		case Calculator::Name::THOLE_LONG: out<<"THOLE_LONG"; break;
		case Calculator::Name::CGEMM_CUT: out<<"CGEMM_CUT"; break;
		case Calculator::Name::CGEMM_LONG: out<<"CGEMM_LONG"; break;
		default: out<<"NONE"; break;
	}
	return out;
}

//==== Mix ====

Calculator::Mix Calculator::Mix::read(const char* str){
	if(std::strcmp(str,"ARITHMETIC")==0) return Calculator::Mix::ARITHMETIC;
	else if(std::strcmp(str,"GEOMETRIC")==0) return Calculator::Mix::GEOMETRIC;
	else if(std::strcmp(str,"HARMONIC")==0) return Calculator::Mix::HARMONIC;
	else return Calculator::Mix::NONE;
}

const char* Calculator::Mix::name(const Calculator::Mix& t){
	switch(t){
		case Calculator::Mix::ARITHMETIC: return "ARITHMETIC";
		case Calculator::Mix::GEOMETRIC: return "GEOMETRIC";
		case Calculator::Mix::HARMONIC: return "HARMONIC";
		default: return "NONE";
	}
}

std::ostream& operator<<(std::ostream& out, const Calculator::Mix& t){
	switch(t){
		case Calculator::Mix::ARITHMETIC: out<<"ARITHMETIC"; break;
		case Calculator::Mix::GEOMETRIC: out<<"GEOMETRIC"; break;
		case Calculator::Mix::HARMONIC: out<<"HARMONIC"; break;
		default: out<<"NONE"; break;
	}
	return out;
}

//==== contructors/destructors ====

Calculator::Calculator(Name name, double rc):name_(name){
    if(rc<0) throw std::invalid_argument("Calculator::Calculator(double): Invalid cutoff radius.");
    rc_=rc;
    rc2_=rc_*rc_;
}

//==== operator ====

std::ostream& operator<<(std::ostream& out, const Calculator& calc){
	return out<<"calc name "<<calc.name()<<" rc "<<calc.rc();
}

//==== member functions ====

void Calculator::resize(int ntypes){
	if(ntypes<0) throw std::invalid_argument("Calculator::resize(int): Invalid number of types.");
	ntypes_=ntypes;
}

void Calculator::read(Token& token){
	//calculator name rc 6.0 ...
	rc_=std::atof(token.next().c_str());
	if(rc_<=0) throw std::invalid_argument("Calculator::read(Token&): invalid cutoff.");
	rc2_=rc_*rc_;
}

//**********************************************
// serialization
//**********************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const Calculator& obj){
		if(CALC_PRINT_FUNC>0) std::cout<<"nbytes(const Calculator&):\n";
		int size=0;
		size+=sizeof(Calculator::Name);//name_
		size+=sizeof(double);//rcut_
		return size;
	}
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const Calculator& obj, char* arr){
		if(CALC_PRINT_FUNC>0) std::cout<<"pack(const Calculator&,char*):\n";
		int pos=0;
		std::memcpy(arr+pos,&obj.name(),sizeof(Calculator::Name)); pos+=sizeof(Calculator::Name);//name_
		std::memcpy(arr+pos,&obj.rc(),sizeof(double)); pos+=sizeof(double);//rcut_
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(Calculator& obj, const char* arr){
		if(CALC_PRINT_FUNC>0) std::cout<<"unpack(Calculator&,const char*):\n";
		int pos=0;
		Calculator::Name name;
		double rc=0;
		std::memcpy(&name,arr+pos,sizeof(Calculator::Name)); pos+=sizeof(Calculator::Name);//name_
		std::memcpy(&rc,arr+pos,sizeof(double)); pos+=sizeof(double);//rcut_
		obj=Calculator(name,rc);
		return pos;
	}
	
}