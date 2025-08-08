#pragma once
#ifndef TIME_HPP
#define TIME_HPP

// c++ libraries
#include <iostream>
#include <chrono>

struct Clock{
private:
	std::chrono::high_resolution_clock::time_point tstart_;
	std::chrono::high_resolution_clock::time_point tstop_;
	double time_;
public:
	//==== constructors/destructors ====
	Clock():time_(0){}
	~Clock(){}
	
	//==== access ====
	std::chrono::high_resolution_clock::time_point& tstart(){return tstart_;}
	const std::chrono::high_resolution_clock::time_point& tstart()const{return tstart_;}
	std::chrono::high_resolution_clock::time_point& tstop(){return tstop_;}
	const std::chrono::high_resolution_clock::time_point& tstop()const{return tstop_;}
	double& time(){return time_;}
	const double& time()const{return time_;}
	
	//==== member functions ====
	void start();
	void stop();
	double duration();
	
};

#endif
