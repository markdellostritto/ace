//c libraries
#include <cstring>
//c++ libraries
#include <iostream>
// structure
#include "struc/atom.hpp"
// str
#include "str/string.hpp"

//**********************************************************************************************
//Atom
//**********************************************************************************************

std::ostream& operator<<(std::ostream& out, const Atom& atom){
	//basic properties
	if(atom.name)   out<<"name ";
	if(atom.an)     out<<"an ";
	if(atom.type)   out<<"type ";
	//serial properties
	if(atom.mass)   out<<"mass ";
	if(atom.charge) out<<"charge ";
	if(atom.radius) out<<"radius ";
	if(atom.eta)    out<<"eta ";
	//vector properties
	if(atom.image)  out<<"image ";
	if(atom.posn)   out<<"posn ";
	if(atom.vel)    out<<"vel ";
	if(atom.force)  out<<"force ";
	//symm
	if(atom.symm)	out<<"symm ";
	//return
	return out;
}

void Atom::defaults(){
	//basic properties
	name   = false;
	an     = false;
	type   = false;
	//serial properties
	mass   = false;
	charge = false;
	radius = false;
	eta    = false;
	//vector properties
	image  = false;
	posn   = false;
	vel    = false;
	force  = false;
	//nnp
	symm   = false;
}

Atom Atom::read(Token& token){
	Atom atom;
	while(!token.end()){
		const std::string tag=string::to_upper(token.next());
		//basic properties
		if(tag=="NAME") atom.name=true;
		else if(tag=="AN") atom.an=true;
		else if(tag=="TYPE") atom.type=true;
		//serial properties
		else if(tag=="MASS") atom.mass=true;
		else if(tag=="CHARGE") atom.charge=true;
		else if(tag=="RADIUS") atom.radius=true;
		else if(tag=="ETA") atom.eta=true;
		//vector properties
		else if(tag=="IMAGE") atom.image=true;
		else if(tag=="POSN") atom.posn=true;
		else if(tag=="VEL") atom.vel=true;
		else if(tag=="FORCE") atom.force=true;
		//nnp
		else if(tag=="SYMM") atom.symm=true;
	}
	return atom;
}

//**********************************************************************************************
// serialization
//**********************************************************************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const Atom& obj){
		int size=0;
		//basic properties
		size+=sizeof(obj.name);
		size+=sizeof(obj.an);
		size+=sizeof(obj.type);
		//serial properties
		size+=sizeof(obj.mass);
		size+=sizeof(obj.charge);
		size+=sizeof(obj.radius);
		size+=sizeof(obj.eta);
		//vector properties
		size+=sizeof(obj.image);
		size+=sizeof(obj.posn);
		size+=sizeof(obj.vel);
		size+=sizeof(obj.force);
		//nnp
		size+=sizeof(obj.symm);
		//return
		return size;
	}
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const Atom& obj, char* arr){
		int pos=0;
		//basic properties
		std::memcpy(arr+pos,&obj.name,sizeof(bool)); pos+=sizeof(bool);
		std::memcpy(arr+pos,&obj.an,sizeof(bool)); pos+=sizeof(bool);
		std::memcpy(arr+pos,&obj.type,sizeof(bool)); pos+=sizeof(bool);
		//serial properties
		std::memcpy(arr+pos,&obj.mass,sizeof(bool)); pos+=sizeof(bool);
		std::memcpy(arr+pos,&obj.charge,sizeof(bool)); pos+=sizeof(bool);
		std::memcpy(arr+pos,&obj.radius,sizeof(bool)); pos+=sizeof(bool);
		std::memcpy(arr+pos,&obj.eta,sizeof(bool)); pos+=sizeof(bool);
		//vector properties
		std::memcpy(arr+pos,&obj.image,sizeof(bool)); pos+=sizeof(bool);
		std::memcpy(arr+pos,&obj.posn,sizeof(bool)); pos+=sizeof(bool);
		std::memcpy(arr+pos,&obj.vel,sizeof(bool)); pos+=sizeof(bool);
		std::memcpy(arr+pos,&obj.force,sizeof(bool)); pos+=sizeof(bool);
		//nnp
		std::memcpy(arr+pos,&obj.symm,sizeof(bool)); pos+=sizeof(bool);
		//return
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(Atom& obj, const char* arr){
		int pos=0;
		//basic properties
		std::memcpy(&obj.name,arr+pos,sizeof(bool)); pos+=sizeof(bool);
		std::memcpy(&obj.an,arr+pos,sizeof(bool)); pos+=sizeof(bool);
		std::memcpy(&obj.type,arr+pos,sizeof(bool)); pos+=sizeof(bool);
		//serial properties
		std::memcpy(&obj.mass,arr+pos,sizeof(bool)); pos+=sizeof(bool);
		std::memcpy(&obj.charge,arr+pos,sizeof(bool)); pos+=sizeof(bool);
		std::memcpy(&obj.radius,arr+pos,sizeof(bool)); pos+=sizeof(bool);
		std::memcpy(&obj.eta,arr+pos,sizeof(bool)); pos+=sizeof(bool);
		//vector properties
		std::memcpy(&obj.image,arr+pos,sizeof(bool)); pos+=sizeof(bool);
		std::memcpy(&obj.posn,arr+pos,sizeof(bool)); pos+=sizeof(bool);
		std::memcpy(&obj.vel,arr+pos,sizeof(bool)); pos+=sizeof(bool);
		std::memcpy(&obj.force,arr+pos,sizeof(bool)); pos+=sizeof(bool);
		//nnp
		std::memcpy(&obj.symm,arr+pos,sizeof(bool)); pos+=sizeof(bool);
		//return
		return pos;
	}
	
}
