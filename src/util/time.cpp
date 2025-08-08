//c++ libraries
#include "time.hpp"

//==== member functions ====

void Clock::start(){
	tstart_ = std::chrono::high_resolution_clock::now();
}

void Clock::stop(){
	tstop_ = std::chrono::high_resolution_clock::now();
}

double Clock::duration(){
	time_ = (std::chrono::duration_cast<std::chrono::duration<double>>(tstop_-tstart_)).count();
	return time_;
}