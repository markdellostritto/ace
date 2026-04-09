#pragma once
#ifndef TIME_HPP
#define TIME_HPP

// c++ libraries
#include <iostream>
#include <chrono>

//*************************************************************************
// CLOCK 
//*************************************************************************

struct Clock{
public:
	class Unit{
	public:
		//enum
		enum Type{
			NS,//nano-seconds
			US,//micro-seconds
			MS,//milli-seconds
			S,//seconds
			NONE
		};
		//constructor
		Unit():t_(Type::NONE){}
		Unit(Type t):t_(t){}
		//operators
		operator Type()const{return t_;}
		//member functions
		static Unit read(const char* str);
		static const char* name(const Unit& unit);
	private:
		Type t_;
	};
private:
	std::chrono::high_resolution_clock::time_point tstart_;
	std::chrono::high_resolution_clock::time_point tstop_;
public:
	//==== constructors/destructors ====
	Clock(){}
	~Clock(){}
	
	//==== access ====
	std::chrono::high_resolution_clock::time_point& tstart(){return tstart_;}
	const std::chrono::high_resolution_clock::time_point& tstart()const{return tstart_;}
	std::chrono::high_resolution_clock::time_point& tstop(){return tstop_;}
	const std::chrono::high_resolution_clock::time_point& tstop()const{return tstop_;}
	
	//==== member functions ====
	void start();
	void stop();
	double duration(const Clock::Unit& unit);	
};
std::ostream& operator<<(std::ostream& out, const Clock::Unit& unit);

#endif
