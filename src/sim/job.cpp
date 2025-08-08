// c
#include <cstring>
// ace
#include "sim/job.hpp"

//****************************************************************************
// Job
//****************************************************************************

Job Job::read(const char* str){
	if(std::strcmp(str,"SP")==0) return Job::SP;
	else if(std::strcmp(str,"MD")==0) return Job::MD;
	else if(std::strcmp(str,"MIN")==0) return Job::MIN;
	else return Job::UNKNOWN;
}

const char* Job::name(const Job& t){
	switch(t){
		case Job::SP: return "SP";
		case Job::MD: return "MD";
		case Job::MIN: return "MIN";
		default: return "UNKNOWN";
	}
}

std::ostream& operator<<(std::ostream& out, const Job& t){
	switch(t){
		case Job::SP: out<<"SP"; break;
		case Job::MD: out<<"MD"; break;
		case Job::MIN: out<<"MIN"; break;
		default: out<<"UNKNOWN"; break;
	}
	return out;
}
