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
	if(atom.name())   out<<"name ";
	if(atom.an())     out<<"an ";
	if(atom.type())   out<<"type ";
	//serial properties
	if(atom.mass())   out<<"mass ";
	if(atom.charge()) out<<"charge ";
	if(atom.radius()) out<<"radius ";
	if(atom.eta())    out<<"eta ";
	//vector properties
	if(atom.image())  out<<"image ";
	if(atom.posn())   out<<"posn ";
	if(atom.vel())    out<<"vel ";
	if(atom.force())  out<<"force ";
	if(atom.dipole()) out<<"dipole ";
	//symm
	if(atom.symm())	out<<"symm ";
	//return
	return out;
}

Atom& Atom::read(Token& token){
	while(!token.end()){
		const std::string tag=string::to_upper(token.next());
		//basic properties
		if(tag=="NAME") name()=true;
		else if(tag=="AN") an()=true;
		else if(tag=="TYPE") type()=true;
		//serial properties
		else if(tag=="MASS") mass()=true;
		else if(tag=="CHARGE") charge()=true;
		else if(tag=="RADIUS") radius()=true;
		else if(tag=="ETA") eta()=true;
		//vector properties
		else if(tag=="IMAGE") image()=true;
		else if(tag=="POSN") posn()=true;
		else if(tag=="VEL") vel()=true;
		else if(tag=="FORCE") force()=true;
		else if(tag=="DIPOLE") dipole()=true;
		//nnp
		else if(tag=="SYMM") symm()=true;
	}
	return *this;
}

//**********************************************************************************************
// serialization
//**********************************************************************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const Atom& obj){
		return Atom::SIZE*sizeof(bool);
	}
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const Atom& obj, char* arr){
		int pos=0;
		std::memcpy(arr+pos,obj.data().data(),Atom::SIZE*sizeof(bool)); pos+=Atom::SIZE*sizeof(bool);
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(Atom& obj, const char* arr){
		int pos=0;
		std::memcpy(obj.data().data(),arr+pos,Atom::SIZE*sizeof(bool)); pos+=Atom::SIZE*sizeof(bool);
		return pos;
	}
	
}
