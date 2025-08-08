// c++ libraries
#include <iostream>
// ann - math
#include "math/const.hpp"
#include "math/special.hpp"
// ann - cell
#include "struc/cell.hpp"

//==== operators ====

std::ostream& operator<<(std::ostream& out, const Cell& cell){
	const Eigen::Matrix3d& R=cell.R_;
	const Eigen::Matrix3d& K=cell.K_;
	out<<"V = "<<cell.vol_<<"\n";
	out<<"R = "; for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) std::cout<<R(i,j)<<" "; std::cout<<"\n";
	out<<"K = "; for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) std::cout<<K(i,j)<<" ";
	return out;
}

//==== member functions ====

void Cell::defaults(){
	vol_=0;
	pbc_=Eigen::Vector3i::Zero();
	R_=Eigen::Matrix3d::Zero();
	RInv_=Eigen::Matrix3d::Zero();
	K_=Eigen::Matrix3d::Zero();
	KInv_=Eigen::Matrix3d::Zero();
	dMax_=0;
	dMax2_=0;
	shifts_.resize(6,Eigen::Vector3d::Zero());
}

void Cell::init(const Eigen::Matrix3d& R){
	//compute matrix
	R_=R;
	vol_=std::fabs(R_.determinant());
	RInv_.noalias()=R_.inverse();
	//compute reciprocal matrix
	K_.col(0)=2.0*math::constants::PI*R_.col(1).cross(R_.col(2))/(R_.col(0).dot(R_.col(1).cross(R_.col(2))));
	K_.col(1)=2.0*math::constants::PI*R_.col(2).cross(R_.col(0))/(R_.col(1).dot(R_.col(2).cross(R_.col(0))));
	K_.col(2)=2.0*math::constants::PI*R_.col(0).cross(R_.col(1))/(R_.col(2).dot(R_.col(0).cross(R_.col(1))));
	KInv_.noalias()=K_.inverse();
	//compute large distance (half of min lattice vector)
	dMax_=R_.col(0).norm();
	for(int i=1; i<3; ++i){
		const double norm=R.col(i).norm();
		if(norm<dMax_) dMax_=norm;
	}
	dMax_*=0.5;
	dMax2_=dMax_*dMax_;
	//compute shifts
	shifts_[0]=R_.col(0);
	shifts_[1]=-R_.col(0);
	shifts_[2]=R_.col(1);
	shifts_[3]=-R_.col(1);
	shifts_[4]=R_.col(2);
	shifts_[5]=-R_.col(2);
}

//==== static functions - vector operations ====

Eigen::Vector3d& Cell::fracToCart(const Eigen::Vector3d& vFrac, Eigen::Vector3d& vCart, const Eigen::Matrix3d& R){
	vCart=R*vFrac;//not assuming vCart!=vFrac
	return vCart;
}

Eigen::Vector3d& Cell::cartToFrac(const Eigen::Vector3d& vCart, Eigen::Vector3d& vFrac, const Eigen::Matrix3d& RInv){
	vFrac=RInv*vCart;//not assuming vCart!=vFrac
	return vFrac;
}

Eigen::Vector3d& Cell::returnToCell(const Eigen::Vector3d& v1, Eigen::Vector3d& v2, const Eigen::Matrix3d& R, const Eigen::Matrix3d& RInv){
	//convert the vector to fractional coordinates
	v2=RInv*v1;
	//use the method of images to return the vector to the cell
	v2[0]=math::special::mod(v2[0],0.0,1.0);
	v2[1]=math::special::mod(v2[1],0.0,1.0);
	v2[2]=math::special::mod(v2[2],0.0,1.0);
	//return the vector in cartesian coordinates
	v2=R*v2;
	return v2;
}

//==== vector operations ====

Eigen::Vector3d& Cell::sum(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, Eigen::Vector3d& sum)const{
	Eigen::Vector3d v2n=-1*v2;
	return diff(v1,v2n,sum);
}

Eigen::Vector3d& Cell::diff(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, Eigen::Vector3d& diff)const{
	//find difference (in fractional coordinates)
	diff.noalias()=RInv_*(v1-v2);
	//return to cell
	diff[0]=math::special::mod(diff[0],-0.5,0.5);
	diff[1]=math::special::mod(diff[1],-0.5,0.5);
	diff[2]=math::special::mod(diff[2],-0.5,0.5);
	//switch back to Cartesian coordinates
	diff=R_*diff;
	//check distance
	if(diff.squaredNorm()>dMax2_){
		//large shear - manual search required
		Eigen::Vector3d rmin=diff;
		double dmin2=rmin.squaredNorm();
		for(int i=0; i<6; ++i){
			const Eigen::Vector3d tmp=diff+shifts_[i];
			const double dist2=tmp.squaredNorm();
			if(dist2<dmin2){
				dmin2=dist2;
				rmin=tmp;
			}
		}
		diff=rmin;
	}
	return diff;
}

double Cell::dist(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2)const{
	Eigen::Vector3d tmp;
	diff(v1,v2,tmp);
	return tmp.norm();
	return (R_*tmp).norm();
}

double Cell::dist(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, Eigen::Vector3d& tmp)const{
	diff(v1,v2,tmp);
	return tmp.norm();
}

double Cell::dist2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, Eigen::Vector3d& tmp)const{
	diff(v1,v2,tmp);
	return tmp.squaredNorm();
}

Eigen::Vector3d& Cell::modv(const Eigen::Vector3d& v1, Eigen::Vector3d& v2){
	//convert the vector to fractional coordinates
	v2=RInv_*v1;
	//use the method of images to return the vector to the cell
	v2[0]=math::special::mod(v2[0],0.0,1.0);
	v2[1]=math::special::mod(v2[1],0.0,1.0);
	v2[2]=math::special::mod(v2[2],0.0,1.0);
	//return the vector in cartesian coordinates
	v2=R_*v2;
	return v2;
}

//**********************************************************************************************
// serialization
//**********************************************************************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************

	template <> int nbytes(const Cell& cell){
		int size=0;
		size+=3*sizeof(int);
		size+=9*sizeof(double);
		return size;
	};
	
	//**********************************************
	// packing
	//**********************************************

	template <> int pack(const Cell& cell, char* arr){
		int pos=0;
		std::memcpy(arr+pos,cell.pbc().data(),3*sizeof(int)); pos+=3*sizeof(int);
		std::memcpy(arr+pos,cell.R().data(),9*sizeof(double)); pos+=9*sizeof(double);
		return pos;
	};
	
	//**********************************************
	// unpacking
	//**********************************************

	template <> int unpack(Cell& obj, const char* arr){
		Eigen::Vector3i pbc;
		Eigen::Matrix3d R;
		int pos=0;
		std::memcpy(pbc.data(),arr+pos,3*sizeof(int)); pos+=3*sizeof(int);
		std::memcpy(R.data(),arr+pos,9*sizeof(double)); pos+=9*sizeof(double);
		obj.init(R);
		obj.pbc()=pbc;
		return pos;
	}
	
}
