#pragma once
#ifndef CELL_HPP
#define CELL_HPP

// c++ libraries
#include <iosfwd>
// eigen libraries
#include <Eigen/Dense>
//mem
#include "mem/serialize.hpp"

//****************************************************************
//Cell class
//****************************************************************

class Cell{
private:
	//==== members ====
	double vol_;//the volume of the simulation cell
	double dMax_,dMax2_;//largest distances
	std::vector<Eigen::Vector3d> shifts_;
	Eigen::Vector3i pbc_;
	Eigen::Matrix3d R_;//the lattice vector matrix (lattice vectors are columns of the matrix)
	Eigen::Matrix3d RInv_;//the inverse of the lattice vector matrix
	Eigen::Matrix3d K_;//the repiprocal lattice vector matrix (lattice vectors are columns of the matrix
	Eigen::Matrix3d KInv_;//the inverse of the reciprocal lattice vector matrix
public:	
	//==== constructors/destructors ====
	Cell(){defaults();}
	Cell(const Eigen::Matrix3d& R){init(R);}
	~Cell(){}
	
	//==== operators ====
	friend std::ostream& operator<<(std::ostream& out, const Cell& cell);
	
	//==== access ====
	const double& dMax()const{return dMax_;}
	const double& vol()const{return vol_;}
	Eigen::Vector3i& pbc(){return pbc_;}
	const Eigen::Vector3i& pbc()const{return pbc_;}
	const Eigen::Matrix3d& R()const{return R_;}
	const Eigen::Matrix3d& RInv()const{return RInv_;}
	const Eigen::Matrix3d& K()const{return K_;}
	const Eigen::Matrix3d& KInv()const{return KInv_;}
	
	//==== static functions - vector operations ====
	static Eigen::Vector3d& fracToCart(const Eigen::Vector3d& vFrac, Eigen::Vector3d& vCart, const Eigen::Matrix3d& R);
	static Eigen::Vector3d& cartToFrac(const Eigen::Vector3d& vCart, Eigen::Vector3d& vFrac, const Eigen::Matrix3d& RInv);
	static Eigen::Vector3d& returnToCell(const Eigen::Vector3d& v1, Eigen::Vector3d& v2, const Eigen::Matrix3d& R, const Eigen::Matrix3d& RInv);
	
	//==== vector operations ====
	Eigen::Vector3d& sum(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, Eigen::Vector3d& sum)const;
	Eigen::Vector3d& diff(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, Eigen::Vector3d& diff)const;
	double dist(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2)const;
	double dist(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, Eigen::Vector3d& tmp)const;
	double dist2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, Eigen::Vector3d& tmp)const;
	Eigen::Vector3d& modv(const Eigen::Vector3d& v1, Eigen::Vector3d& v2);
	
	//==== member functions ====
	void defaults();
	void clear(){defaults();}
	void init(const Eigen::Matrix3d& R);
};

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************

	template <> int nbytes(const Cell& obj);
	
	//**********************************************
	// packing
	//**********************************************

	template <> int pack(const Cell& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************

	template <> int unpack(Cell& obj, const char* arr);
	
}

#endif
