// c libraries
#include <cstring>
// time
#include "util/time.hpp"

//*************************************************************************
// UNIT 
//*************************************************************************

Clock::Unit Clock::Unit::read(const char* str){
	if(std::strcmp(str,"NS")==0) return Clock::Unit::NS;
	else if(std::strcmp(str,"US")==0) return Clock::Unit::US;
	else if(std::strcmp(str,"MS")==0) return Clock::Unit::MS;
	else if(std::strcmp(str,"S")==0) return Clock::Unit::S;
	else return Clock::Unit::NONE;
}

const char* Clock::Unit::name(const Clock::Unit& t){
	switch(t){
		case Clock::Unit::NS: return "NS";
		case Clock::Unit::US: return "US";
		case Clock::Unit::MS: return "MS";
		case Clock::Unit::S: return "S";
		default: return "NONE";
	}
}

std::ostream& operator<<(std::ostream& out, const Clock::Unit& t){
	switch(t){
		case Clock::Unit::NS: out<<"NS"; break;
		case Clock::Unit::US: out<<"US"; break;
		case Clock::Unit::MS: out<<"MS"; break;
		case Clock::Unit::S: out<<"S"; break;
		default: out<<"NONE"; break;
	}
	return out;
}

//*************************************************************************
// CLOCK 
//*************************************************************************

//==== member functions ====

void Clock::start(){
	tstart_ = std::chrono::high_resolution_clock::now();
}

void Clock::stop(){
	tstop_ = std::chrono::high_resolution_clock::now();
}

double Clock::duration(const Clock::Unit& unit){
	switch(unit){
		case Clock::Unit::NS: 
			return std::chrono::duration_cast<std::chrono::nanoseconds>(tstop_-tstart_).count();
		case Clock::Unit::US: 
			return std::chrono::duration_cast<std::chrono::microseconds>(tstop_-tstart_).count();
		case Clock::Unit::MS: 
			return std::chrono::duration_cast<std::chrono::milliseconds>(tstop_-tstart_).count();
		case Clock::Unit::S: 
			return std::chrono::duration_cast<std::chrono::seconds>(tstop_-tstart_).count();
		default: 
			return 0;
	}
}