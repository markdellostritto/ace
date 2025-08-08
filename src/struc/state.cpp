//c libraries
#include <stdexcept>
#include <cstring>
//c++ libraries
#include <iostream>
// structure
#include "struc/state.hpp"

//**********************************************************************************************
//State
//**********************************************************************************************

//==== operators ====

std::ostream& operator<<(std::ostream& out, const State& obj){
	//energy
	out<<"etot  = "<<obj.etot_<<"\n";
	out<<"ecoul = "<<obj.ecoul_<<"\n";
	out<<"evdw  = "<<obj.evdw_<<"\n";
	out<<"erep  = "<<obj.erep_<<"\n";
	out<<"pe    = "<<obj.pe_<<"\n";
	out<<"ke    = "<<obj.ke_<<"\n";
	//state
	out<<"temp  = "<<obj.temp_<<"\n";
	out<<"press = "<<obj.press_<<"\n";
	//charge
	out<<"qtot  = "<<obj.qtot_<<"\n";
	//time
	out<<"dt    = "<<obj.dt_<<"\n";
	out<<"t     = "<<obj.t_;
	return out;
}

//==== member functions ====

void State::clear(){
	//energy
	etot_=0;
	ecoul_=0;
	evdw_=0;
	erep_=0;
	pe_=0;
	ke_=0;
	//state
	temp_=0;
	press_=0;
	//charge
	qtot_=0;
	//time
	dt_=0;
	t_=0;
}

//**********************************************************************************************
// serialization
//**********************************************************************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const State& obj){
		int size=0;
		//energy
		size+=sizeof(obj.etot());
		size+=sizeof(obj.ecoul());
		size+=sizeof(obj.evdw());
		size+=sizeof(obj.erep());
		size+=sizeof(obj.pe());
		size+=sizeof(obj.ke());
		//state
		size+=sizeof(obj.temp());
		size+=sizeof(obj.press());
		//charge
		size+=sizeof(obj.qtot());
		//time
		size+=sizeof(obj.dt());
		size+=sizeof(obj.t());
		//return
		return size;
	}
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const State& obj, char* arr){
		int pos=0;
		//energy
		std::memcpy(arr+pos,&obj.etot(),sizeof(double)); pos+=sizeof(double);
		std::memcpy(arr+pos,&obj.ecoul(),sizeof(double)); pos+=sizeof(double);
		std::memcpy(arr+pos,&obj.evdw(),sizeof(double)); pos+=sizeof(double);
		std::memcpy(arr+pos,&obj.erep(),sizeof(double)); pos+=sizeof(double);
		std::memcpy(arr+pos,&obj.pe(),sizeof(double)); pos+=sizeof(double);
		std::memcpy(arr+pos,&obj.ke(),sizeof(double)); pos+=sizeof(double);
		//state
		std::memcpy(arr+pos,&obj.temp(),sizeof(double)); pos+=sizeof(double);
		std::memcpy(arr+pos,&obj.press(),sizeof(double)); pos+=sizeof(double);
		//charge
		std::memcpy(arr+pos,&obj.qtot(),sizeof(double)); pos+=sizeof(double);
		//time
		std::memcpy(arr+pos,&obj.dt(),sizeof(double)); pos+=sizeof(double);
		std::memcpy(arr+pos,&obj.t(),sizeof(int)); pos+=sizeof(int);
		//return
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(State& obj, const char* arr){
		int pos=0;
		//energy
		std::memcpy(&obj.etot(),arr+pos,sizeof(double)); pos+=sizeof(double);
		std::memcpy(&obj.ecoul(),arr+pos,sizeof(double)); pos+=sizeof(double);
		std::memcpy(&obj.evdw(),arr+pos,sizeof(double)); pos+=sizeof(double);
		std::memcpy(&obj.erep(),arr+pos,sizeof(double)); pos+=sizeof(double);
		std::memcpy(&obj.pe(),arr+pos,sizeof(double)); pos+=sizeof(double);
		std::memcpy(&obj.ke(),arr+pos,sizeof(double)); pos+=sizeof(double);
		//state
		std::memcpy(&obj.temp(),arr+pos,sizeof(double)); pos+=sizeof(double);
		std::memcpy(&obj.press(),arr+pos,sizeof(double)); pos+=sizeof(double);
		//charge
		std::memcpy(&obj.qtot(),arr+pos,sizeof(double)); pos+=sizeof(double);
		//time
		std::memcpy(&obj.dt(),arr+pos,sizeof(double)); pos+=sizeof(double);
		std::memcpy(&obj.t(),arr+pos,sizeof(int)); pos+=sizeof(int);
		//return
		return pos;
	}
	
}
