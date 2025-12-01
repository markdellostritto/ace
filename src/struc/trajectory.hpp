#pragma once
#ifndef TRAJECTORY_HPP
#define TRAJECTORY_HPP

//no bounds checking in Eigen
#define EIGEN_NO_DEBUG

//c++ libraries
#include <iosfwd>
//Eigen
#include <Eigen/Dense>
// struc
#include "struc/structure.hpp"
#include "struc/interval.hpp"
// mem 
#include "mem/serialize.hpp"
// string
#include "str/string.hpp"

#ifndef TRAJ_PRINT_FUNC
#define TRAJ_PRINT_FUNC 0
#endif

//**********************************************
// Trajectory
//**********************************************

class Trajectory{
private:
	std::string name_;
	double timestep_;
	int timesteps_;
	Atom atom_;
	std::vector<Structure> frames_;
public:
	//==== constructors/destructors ====
	Trajectory(){defaults();}
	~Trajectory(){}
	
	//==== operators ====
	friend std::ostream& operator<<(std::ostream& out, const Trajectory& traj);
	
	//==== access ====
	std::string& name(){return name_;}
	const std::string& name()const{return name_;}
	double& timestep(){return timestep_;}
	const double& timestep()const{return timestep_;}
	int& timesteps(){return timesteps_;}
	const int& timesteps()const{return timesteps_;}
	Atom& atom(){return atom_;}
	const Atom& atom()const{return atom_;}
	Structure& frame(int i){return frames_[i];}
	const Structure& frame(int i)const{return frames_[i];}
	
	//==== member functions ====
	void defaults();
	void clear();
	void resize(int ts, int nAtoms, const Atom& atom);
	void resize(int ts);
	
	//==== static functions ====
	static void set_image(Trajectory & traj);
	static void unwrap(Trajectory & traj);
};


//**********************************************************************************************
// serialization
//**********************************************************************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const Trajectory& traj);
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const Trajectory& traj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(Trajectory& traj, const char* arr);
	
}

#endif