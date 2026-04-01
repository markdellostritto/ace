#include "sim/constraint.hpp"

//==== Name ====

Constraint::Name Constraint::Name::read(const char* str){
	if(std::strcmp(str,"FREEZE")==0) return Constraint::Name::FREEZE;
	else return Constraint::Name::NONE;
}

const char* Constraint::Name::name(const Constraint::Name& t){
	switch(t){
		case Constraint::Name::FREEZE: return "FREEZE";
		default: return "NONE";
	}
}

std::ostream& operator<<(std::ostream& out, const Constraint::Name& t){
	switch(t){
		case Constraint::Name::FREEZE: out<<"FREEZE"; break;
		default: out<<"NONE"; break;
	}
	return out;
}

//==== operator ====

std::ostream& operator<<(std::ostream& out, const Constraint& constraint){
	return out<<"constraint name "<<constraint.name();
}

//==== member functions ====

void Constraint::read(Token& token){
	//constraint name indices
	indices_.clear();
	while(!token.end()){
		const std::string str=token.next();
		Token atoken(str,":");
		int beg=0,end=0,stride=1;
		beg=std::atoi(atoken.next().c_str());
		if(!atoken.end()){
				end=std::atoi(atoken.next().c_str());
				if(!atoken.end()){
					stride=std::atoi(atoken.next().c_str());
				}
		} else end=beg;
		beg--; end--;
		for(int i=beg; i<=end; i+=stride) indices_.push_back(i);
	}
}

//**********************************************
// serialization
//**********************************************
	
namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const Constraint& obj){
		int size=0;
		size+=sizeof(int);//indices
		size+=sizeof(int)*obj.indices().size();
		return size;
	}
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const Constraint& obj, char* arr){
		int pos=0;
		const int size=obj.indices().size();
		std::memcpy(arr+pos,&size,sizeof(int)); pos+=sizeof(int);//size
		std::memcpy(arr+pos,obj.indices().data(),sizeof(int)*size); pos+=sizeof(int)*size;//indices
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(Constraint& obj, const char* arr){
		int pos=0;
		int size=0;
		std::memcpy(&size,arr+pos,sizeof(int)); pos+=sizeof(int);//size
		obj.indices().resize(size);
		std::memcpy(obj.indices().data(),arr+pos,sizeof(int)*size); pos+=sizeof(int)*size;//indices
		return pos;
	}
	
}