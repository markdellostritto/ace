// c libraries
#include <cstdlib>
#include <cstring>
// c++ libraries
#include <iostream>
// str
#include "str/string.hpp"
// ann - serialize
#include "mem/serialize.hpp"
// type
#include "struc/type.hpp"
	
//************************************************************
// TYPE
//************************************************************

void Type::clear(){
	if(TYPE_PRINT_FUNC>0) std::cout<<"Type::clear():\n";
	//name
	name_=std::string("NULL");
	id_=string::hash(name_);
	//standard
	energy_=0;
	mass_=0;
	charge_=0;
	z_=0;
	//radius
	rvdw_=0;
	rcov_=0;
}

Type& Type::read(const char* str, Type& type){
	if(TYPE_PRINT_FUNC>0) std::cout<<"Type::read(const char*,Type&):\n";
	Token token(str,string::WSC);
	while(!token.end()){
		const std::string tag=string::to_upper(token.next());
		//name
		if(tag=="NAME"){
			type.name()=token.next();
		} 
		//standard
		if(tag=="ENERGY"){
			type.energy()=std::atof(token.next().c_str());
		} else if(tag=="MASS"){
			type.mass()=std::atof(token.next().c_str());
		} else if(tag=="CHARGE"){
			type.charge()=std::atof(token.next().c_str());
		} else if(tag=="Z"){
			type.z()=std::atof(token.next().c_str());
		}
		//radius
		if(tag=="RVDW"){
			type.rvdw()=std::atof(token.next().c_str());
		} else if(tag=="RCOV"){
			type.rcov()=std::atof(token.next().c_str());
		}
	}
	type.id()=string::hash(type.name());
	return type;
}

Type& Type::read(Type& type, Token& token){
	if(TYPE_PRINT_FUNC>0) std::cout<<"Type::read(const char*,Type&):\n";
	while(!token.end()){
		const std::string tag=string::to_upper(token.next());
		//name
		if(tag=="NAME"){
			type.name()=token.next();
		} 
		//standard
		if(tag=="ENERGY"){
			type.energy()=std::atof(token.next().c_str());
		} else if(tag=="MASS"){
			type.mass()=std::atof(token.next().c_str());
		} else if(tag=="CHARGE"){
			type.charge()=std::atof(token.next().c_str());
		} else if(tag=="Z"){
			type.z()=std::atof(token.next().c_str());
		} 
		//radius
		if(tag=="RVDW"){
			type.rvdw()=std::atof(token.next().c_str());
		} else if(tag=="RCOV"){
			type.rcov()=std::atof(token.next().c_str());
		} 
	}
	type.id()=string::hash(type.name());
	return type;
}

std::ostream& operator<<(std::ostream& out, const Type& type){
	//name
	out<<"name "<<type.name_<<" ";
	//standard
	out<<"energy "<<type.energy()<<" ";
	out<<"mass "<<type.mass()<<" ";
	out<<"chg "<<type.charge()<<" ";
	out<<"z "<<type.z()<<" ";
	//radius
	out<<"rvdw "<<type.rvdw()<<" ";
	out<<"rcov "<<type.rcov()<<" ";
	//return
	return out;
}

void Type::write(FILE* out, const Type& type){
	//name
	fprintf(out,"name %s ",type.name().c_str());
	//standard
	fprintf(out,"energy %f ",type.energy());
	fprintf(out,"mass %f ",type.mass());
	fprintf(out,"charge %f ",type.charge());
	fprintf(out,"z %f ",type.z());
	//radius
	fprintf(out,"rvdw %f ",type.rvdw());
	fprintf(out,"rcov %f ",type.rcov());
	//newline
	fprintf(out,"\n");
}

namespace serialize{

//**********************************************
// byte measures
//**********************************************

template <> int nbytes(const Type& obj){
	if(TYPE_PRINT_FUNC>0) std::cout<<"nbytes(const Type&):\n";
	int size=0;
	//name
	size+=nbytes(obj.name());//name
	size+=sizeof(int);//id
	//standard
	size+=sizeof(double)*4;
	//radius
	size+=sizeof(double)*2;
	//return
	return size;
}

//**********************************************
// packing
//**********************************************

template <> int pack(const Type& obj, char* arr){
	if(TYPE_PRINT_FUNC>0) std::cout<<"pack(const Type&,char*):\n";
	int pos=0;
	//name
	pos+=pack(obj.name(),arr+pos);
	std::memcpy(arr+pos,&obj.id(),sizeof(int)); pos+=sizeof(int);
	//standard
	std::memcpy(arr+pos,&obj.energy(),sizeof(double)); pos+=sizeof(double);
	std::memcpy(arr+pos,&obj.mass(),sizeof(double)); pos+=sizeof(double);
	std::memcpy(arr+pos,&obj.charge(),sizeof(double)); pos+=sizeof(double);
	std::memcpy(arr+pos,&obj.z(),sizeof(double)); pos+=sizeof(double);
	//radius
	std::memcpy(arr+pos,&obj.rvdw(),sizeof(double)); pos+=sizeof(double);
	std::memcpy(arr+pos,&obj.rcov(),sizeof(double)); pos+=sizeof(double);
	//return
	return pos;
}

//**********************************************
// unpacking
//**********************************************

template <> int unpack(Type& obj, const char* arr){
	if(TYPE_PRINT_FUNC>0) std::cout<<"unpack(Type&,const char*):\n";
	int pos=0;
	//name
	pos+=unpack(obj.name(),arr+pos);
	std::memcpy(&obj.id(),arr+pos,sizeof(int)); pos+=sizeof(int);//id
	//standard
	std::memcpy(&obj.energy(),arr+pos,sizeof(double)); pos+=sizeof(double);
	std::memcpy(&obj.mass(),arr+pos,sizeof(double)); pos+=sizeof(double);
	std::memcpy(&obj.charge(),arr+pos,sizeof(double)); pos+=sizeof(double);
	std::memcpy(&obj.z(),arr+pos,sizeof(double)); pos+=sizeof(double);
	//radius
	std::memcpy(&obj.rvdw(),arr+pos,sizeof(double)); pos+=sizeof(double);
	std::memcpy(&obj.rcov(),arr+pos,sizeof(double)); pos+=sizeof(double);
	//return
	return pos;
}
	
}
