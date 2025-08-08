#pragma once
#ifndef CONSTRAINT_HPP
#define CONSTRAINT_HPP

// c++
#include <iostream>
#include <vector>
// ace
#include "struc/structure.hpp"
#include "struc/neighbor.hpp"

class Constraint{
public:
	class Name{
	public:
		enum Type{
			FREEZE,
			NONE
		};
		//constructor
		Name():t_(Type::NONE){}
		Name(Type t):t_(t){}
		//operators
		operator Type()const{return t_;}
		friend std::ostream& operator<<(std::ostream& out, const Name& sys);
		//member functions
		static Name read(const char* str);
		static const char* name(const Name& name);
	private:
		Type t_;
		//prevent automatic conversion for other built-in types
		//template<typename T> operator T() const;
	};
protected:
    Name name_;//name
    std::vector<int> indices_;
public:
    //==== contructors/destructors ====
	Constraint():name_(Name::NONE){}
    Constraint(Name name):name_(name){}
    ~Constraint(){}
    
    //==== operators ====
	friend std::ostream& operator<<(std::ostream& out, const Constraint& constraint);
	
    //==== access ====
	const Name& name()const{return name_;}
    std::vector<int>& indices(){return indices_;}
    const std::vector<int>& indices()const{return indices_;}
    int& index(int i){return indices_[i];}
    const int& index(int i)const{return indices_[i];}

    //==== member functions ====
	void read(Token& token);
    virtual double compute(Structure& struc, const NeighborList& nlist){return 0.0;}
};

#endif