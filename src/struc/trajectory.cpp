//c++ libraries
#include <iostream>
// str
#include "str/print.hpp"
#include "str/token.hpp"
// traj
#include "struc/trajectory.hpp"

//==== operators ====

std::ostream& operator<<(std::ostream& out, const Trajectory& traj){
	char* str=new char[print::len_buf];
	out<<print::buf(str)<<"\n";
	out<<print::title("TRAJECTORY",str)<<"\n";
	out<<"SIM  = "<<traj.name_<<"\n";
	out<<"TS   = "<<traj.timestep_<<"\n";
	out<<"T    = "<<traj.timesteps_<<"\n";
	out<<"Atom = "<<traj.atom_<<"\n";
	out<<print::buf(str);
	delete[] str;
	return out;
}

//==== member functions ====

void Trajectory::defaults(){
	if(TRAJ_PRINT_FUNC>0) std::cout<<"Trajectory::defaults():\n";
	name_=std::string("SYSTEM");
	timestep_=0;
	timesteps_=0;
}

void Trajectory::clear(){
	if(TRAJ_PRINT_FUNC>0) std::cout<<"Trajectory::clear():\n";
	frames_.clear();
	defaults();
}

void Trajectory::resize(int ts, int nAtoms, const Atom& atom){
	if(TRAJ_PRINT_FUNC>0) std::cout<<"Trajectory::resize(int,int,const Atom&):\n";
	timesteps_=ts;
	atom_=atom;
	frames_.resize(timesteps_,Structure(nAtoms,atom));
}

void Trajectory::resize(int ts){
	if(TRAJ_PRINT_FUNC>0) std::cout<<"Trajectory::resize(int):\n";
	timesteps_=ts;
	frames_.resize(timesteps_);
}

//==== static functions ====

void Trajectory::set_image(Trajectory & traj){
	for(int t=1; t<traj.timesteps(); ++t){
		for(int n=0; n<traj.frame(t).nAtoms(); ++n){
			const Eigen::Vector3d diff=traj.frame(t).RInv()*(traj.frame(t).posn(n)-traj.frame(t-1).posn(n));
			traj.frame(t).image(n)=traj.frame(t-1).image(n);
			if(diff[0]>0.5) traj.frame(t).image(n)[0]--;
			else if(diff[0]<-0.5) traj.frame(t).image(n)[0]++;
			if(diff[1]>0.5) traj.frame(t).image(n)[1]--;
			else if(diff[1]<-0.5) traj.frame(t).image(n)[1]++;
			if(diff[2]>0.5) traj.frame(t).image(n)[2]--;
			else if(diff[2]<-0.5) traj.frame(t).image(n)[2]++;
		}
	}
}

void Trajectory::unwrap(Trajectory & traj){
	for(int t=0; t<traj.timesteps(); ++t){
		for(int n=0; n<traj.frame(t).nAtoms(); ++n){
			traj.frame(t).posn(n).noalias()+=
				traj.frame(t).R().col(0)*traj.frame(t).image(n)[0]+
				traj.frame(t).R().col(1)*traj.frame(t).image(n)[1]+
				traj.frame(t).R().col(2)*traj.frame(t).image(n)[2];
			const Eigen::Vector3d offset=traj.frame(t).R().col(0)*traj.frame(t).image(n)[0]+
				traj.frame(t).R().col(1)*traj.frame(t).image(n)[1]+
				traj.frame(t).R().col(2)*traj.frame(t).image(n)[2];
			traj.frame(t).image(n).setZero();
		}
	}
}

//**********************************************************************************************
// serialization
//**********************************************************************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const Trajectory& traj){
		if(TRAJ_PRINT_FUNC>0) std::cout<<"nbytes(const Trajectory&):\n";
		int size=0;
		size+=sizeof(int);//timesteps_
		size+=sizeof(int);//natoms
		size+=sizeof(double);//timestep
		size+=nbytes(traj.atom());//atom
		size+=nbytes(traj.name());//name
		for(int t=0; t<traj.timesteps(); ++t){
			size+=nbytes(traj.frame(t));
		}
		return size;
	}
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const Trajectory& traj, char* arr){
		if(TRAJ_PRINT_FUNC>0) std::cout<<"pack(const Trajectory&,char*):\n";
		int pos=0,tmpInt=0;
		std::memcpy(arr+pos,&(tmpInt=traj.timesteps()),sizeof(int)); pos+=sizeof(int);
		std::memcpy(arr+pos,&(tmpInt=traj.frame(0).nAtoms()),sizeof(int)); pos+=sizeof(int);
		std::memcpy(arr+pos,&traj.timestep(),sizeof(double)); pos+=sizeof(double);
		pos+=pack(traj.atom(),arr);
		pos+=pack(traj.name(),arr);
		for(int t=0; t<traj.timesteps(); ++t){
			pos+=pack(traj.frame(t),arr);
		}
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(Trajectory& traj, const char* arr){
		if(TRAJ_PRINT_FUNC>0) std::cout<<"unpack(Trajectory&,const char*):\n";
		int pos=0,nAtoms=0,ts=0;
		Atom atom;
		std::memcpy(&ts,arr+pos,sizeof(int)); pos+=sizeof(int);
		std::memcpy(&nAtoms,arr+pos,sizeof(int)); pos+=sizeof(int);
		std::memcpy(&traj.timestep(),arr+pos,sizeof(double)); pos+=sizeof(double);
		pos+=unpack(atom,arr);
		pos+=unpack(traj.name(),arr);
		traj.resize(ts,nAtoms,atom);
		for(int t=0; t<traj.timesteps(); ++t){
			pos+=unpack(traj.frame(t),arr);
		}
		return pos;
	}
	
}