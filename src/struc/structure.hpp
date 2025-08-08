#pragma once
#ifndef STRUCTURE_HPP
#define STRUCTURE_HPP

//no bounds checking in Eigen
#define EIGEN_NO_DEBUG

// c++ libraries
#include <iosfwd>
// Eigen
#include <Eigen/Dense>
// ann - structure
#include "struc/cell.hpp"
#include "struc/state.hpp"
#include "struc/atom.hpp"
// ann - serialize
#include "mem/serialize.hpp"

#ifndef STRUC_PRINT_FUNC
#define STRUC_PRINT_FUNC 0
#endif

#ifndef STRUC_PRINT_STATUS
#define STRUC_PRINT_STATUS 0
#endif

#ifndef STRUC_PRINT_DATA
#define STRUC_PRINT_DATA 0
#endif

typedef Eigen::Matrix<int,3,1> Vec3i;
typedef Eigen::Matrix<double,3,1> Vec3d;
typedef Eigen::Matrix<double,Eigen::Dynamic,1> VecXd;

//**********************************************************************************************
//AtomData
//**********************************************************************************************

class AtomData{
protected:
	//atom type
	Atom atom_;
	//number of atoms
	int nAtoms_;
	//basic properties
	std::vector<std::string> name_;//name
	std::vector<int>	an_;//atomic_number
	std::vector<int>	type_;//type
	//serial properties
	std::vector<double>	mass_;//mass
	std::vector<double>	charge_;//charge
	std::vector<double> radius_;//radius
	std::vector<double>	eta_;//eta
	//vector properties
	std::vector<Vec3i>	image_;//image
	std::vector<Vec3d>	posn_;//position
	std::vector<Vec3d>	vel_;//velocity
	std::vector<Vec3d>	force_;//force
	//nnp
	std::vector<VecXd>	symm_;
public:
	//==== constructors/destructors ====
	AtomData():nAtoms_(0){}
	~AtomData(){}
	
	//==== operators ====
	friend std::ostream& operator<<(std::ostream& out, const AtomData& ad);
	
	//==== access - global ====
	const Atom& atom()const{return atom_;}
	const int& nAtoms()const{return nAtoms_;}
	
	//==== access - vectors ====
	//basic properties
	std::vector<std::string>& name(){return name_;}
	const std::vector<std::string>& name()const{return name_;}
	std::vector<int>& an(){return an_;}
	const std::vector<int>& an()const{return an_;}
	std::vector<int>& type(){return type_;}
	const std::vector<int>& type()const{return type_;}
	//serial properties
	std::vector<double>& mass(){return mass_;}
	const std::vector<double>& mass()const{return mass_;}
	std::vector<double>& charge(){return charge_;}
	const std::vector<double>& charge()const{return charge_;}
	std::vector<double>& radius(){return radius_;}
	const std::vector<double>& radius()const{return radius_;}
	std::vector<double>& eta(){return eta_;}
	const std::vector<double>& eta()const{return eta_;}
	//vector properties
	std::vector<Vec3i>& image(){return image_;}
	const std::vector<Vec3i>& image()const{return image_;}
	std::vector<Vec3d>& posn(){return posn_;}
	const std::vector<Vec3d>& posn()const{return posn_;}
	std::vector<Vec3d>& vel(){return vel_;}
	const std::vector<Vec3d>& vel()const{return vel_;}
	std::vector<Vec3d>& force(){return force_;}
	const std::vector<Vec3d>& force()const{return force_;}
	//nnp
	std::vector<VecXd>& symm(){return symm_;}
	const std::vector<VecXd>& symm()const{return symm_;}
	
	//==== access - atoms ====
	//basic properties
	std::string& name(int i){return name_[i];}
	const std::string& name(int i)const{return name_[i];}
	int& an(int i){return an_[i];}
	const int& an(int i)const{return an_[i];}
	int& type(int i){return type_[i];}
	const int& type(int i)const{return type_[i];}
	//serial properties
	double& mass(int i){return mass_[i];}
	const double& mass(int i)const{return mass_[i];}
	double& charge(int i){return charge_[i];}
	const double& charge(int i)const{return charge_[i];}
	double& radius(int i){return radius_[i];}
	const double& radius(int i)const{return radius_[i];}
	double& eta(int i){return eta_[i];}
	const double& eta(int i)const{return eta_[i];}
	//vector properties
	Vec3i& image(int i){return image_[i];}
	const Vec3i& image(int i)const{return image_[i];}
	Vec3d& posn(int i){return posn_[i];}
	const Vec3d& posn(int i)const{return posn_[i];}
	Vec3d& vel(int i){return vel_[i];}
	const Vec3d& vel(int i)const{return vel_[i];}
	Vec3d& force(int i){return force_[i];}
	const Vec3d& force(int i)const{return force_[i];}
	//nnp
	VecXd& symm(int i){return symm_[i];}
	const VecXd& symm(int i)const{return symm_[i];}
	
	//==== member functions ====
	void clear();
	void resize(int nAtoms, const Atom& atom);
};

//**********************************************************************************************
//Structure
//**********************************************************************************************

class Structure: public Cell, public State, public AtomData{
public:
	//==== constructors/destructors ====
	Structure(){}
	Structure(int nAtoms, const Atom& atom){resize(nAtoms,atom);}
	~Structure(){}
	
	//==== operators ====
	friend std::ostream& operator<<(std::ostream& out, const Structure& sim);
	
	//==== member functions ====
	void clear();
	
	//==== static functions ====
	static Structure& super(const Structure& struc, Structure& superc, const Eigen::Vector3i nlat);
	static void write_binary(const Structure& struc, const char* file);
	static void read_binary(Structure& struc, const char* file);
};

//**********************************************************************************************
// serialization
//**********************************************************************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const AtomData& obj);
	template <> int nbytes(const Structure& obj);
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const AtomData& obj, char* arr);
	template <> int pack(const Structure& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(AtomData& obj, const char* arr);
	template <> int unpack(Structure& obj, const char* arr);
	
}

#endif
