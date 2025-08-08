#pragma once
#ifndef ALIAS_HPP
#define ALIAS_HPP

// c++ libraries
#include <ostream>
#include <vector>
#include <string>
// string
#include "str/string.hpp"
#include "str/token.hpp"
// mem
#include "mem/serialize.hpp"

class Alias{
private:
	std::string alias_;
	std::vector<std::string> labels_;
public:
	//==== constructors/destructors ====
	Alias(){}
	~Alias(){}
	
	//==== operators ====
	friend std::ostream& operator<<(std::ostream& out, const Alias& alias);
	
	//==== access functions ====
	std::string& alias(){return alias_;}
	const std::string& alias()const{return alias_;}
	std::vector<std::string>& labels(){return labels_;}
	const std::vector<std::string>& labels()const{return labels_;}
	
	//==== member functions ====
	void clear(){alias_.clear(); labels_.clear();}
	
	//==== static functions ====
	static Alias& read(Token& token, Alias& alias);
};

namespace serialize{
		
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const Alias& obj);
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const Alias& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(Alias& obj, const char* arr);
	
}

#endif