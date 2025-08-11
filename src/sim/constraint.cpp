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
