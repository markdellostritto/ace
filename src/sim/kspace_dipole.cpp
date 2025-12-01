// c libraries
#if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
#include <cmath>
#elif defined __ICC || defined __INTEL_COMPILER
#include <mathimf.h> //intel math library
#endif
// c++ libraries
#include <iostream>
// math
#include "math/const.hpp"
// chem
#include "chem/units.hpp"
// sim
#include "sim/kspace_dipole.hpp"

using math::constants::PI;
using math::constants::RadPI;
using math::constants::Rad2;
using math::constants::ZERO;

namespace KSpace{

//==== operators ====

std::ostream& operator<<(std::ostream& out, const Dipole& d){
	return out<<static_cast<const Base&>(d);
}

//==== member functions ====

void Dipole::init(const Structure& struc){
	if(KSPACED_PRINT_FUNC>0) std::cout<<"KSpace::Dipole::init(const Structure&):\n";
	if(prec_<=0) throw std::invalid_argument("KSpace::Dipole::init(const Structure&): invalid precision.");
	if(rc_<=0) throw std::invalid_argument("KSpace::Dipole::init(const Structure&): invalid rcut.");

    //set structural data
	std::cout<<"setting structural data\n";
	const int& N=struc.nAtoms();
	const Eigen::Matrix3d& R=struc.R();
	const Eigen::Matrix3d& K=struc.K();
	const Eigen::Vector3d L=(Eigen::Vector3d()<<R.col(0).norm(),R.col(1).norm(),R.col(2).norm()).finished();
    const double V=struc.vol();
	m2_=0;
	for(int i=0; i<N; ++i) m2_+=struc.dipole(i).squaredNorm();
	if(m2_==0) throw std::invalid_argument("KSpace::Dipole::init(const Structure&): zero abs dipole.");

    //compute alpha
	std::cout<<"computing alpha\n";
    const double rc2=rc_*rc_;
    const double rc4=rc2*rc2;
    const double rc8=rc4*rc4;
    const double rc9=rc8*rc_;
    //alpha_=std::pow(prec_/(32.0*m2_*std::sqrt(N*V*rc9)*rc4*rc2),1.0/8.0);
	alpha_=(1.35 - 0.15*std::log(prec_))/rc_;
    if(KSPACED_PRINT_DATA>0) std::cout<<"alpha_ = "<<alpha_<<"\n";
    const double a2i=1.0/(alpha_*alpha_);

    //set lattice vectors - reciprocal space
	nk_.setZero();
	for(int i=0; i<3; ++i){
		double err=0;
        const double pre=8.0*Rad2*m2_/V*alpha_*sqrt(2.0/15.0*PI/N);
		do{
			++nk_[i];
            const double kv=(K.col(i)*nk_[i]).norm();
			err=pre*pow(kv,1.5)*exp(-PI*kv/(alpha_*L[i]));
		}while(err>prec_);
	}
	//nk_*=5;//needed to match LAMMPS
	if(KSPACED_PRINT_DATA>0) std::cout<<"nk_ = "<<nk_.transpose()<<"\n";

    //compute k-vecs and k-amps
	const int NK=(2*nk_+Eigen::Vector3i::Constant(1)).prod()-1;
	k_.resize(NK);
	int count=0;
	for(int ix=-nk_[0]; ix<=nk_[0]; ++ix){
		for(int iy=-nk_[1]; iy<=nk_[1]; ++iy){
			for(int iz=-nk_[2]; iz<=nk_[2]; ++iz){
				const Eigen::Vector3d Ktmp=ix*K.col(0)+iy*K.col(1)+iz*K.col(2);
				if(Ktmp.norm()>ZERO) k_[count++]=Ktmp;
			}
		}
	}
	const double kc=PI/struc.vol();
	ka_.resize(NK);
	for(int i=0; i<NK; ++i){
		const double k2=0.25*k_[i].squaredNorm();
		ka_[i]=kc*exp(-k2*a2i)/k2;
	}
	
    vc_=-4.0/3.0*alpha_*alpha_*alpha_/RadPI;
}

double Dipole::energy(Structure& struc)const{
	if(KSPACED_PRINT_FUNC>0) std::cout<<"KSpace::Dipole::energy(const Structure&)const:\n";
	const double ke=units::Consts::ke()*eps_;
	double energy=0;
	const int natoms=struc.nAtoms();
	//kspace - energy
	for(int n=0; n<k_.size(); ++n){
		//compute structure factor
		double sfr=0,sfi=0;
        for(int i=0; i<natoms; ++i){
			const double prod=k_[n].dot(struc.posn(i));
            const double mudotk=struc.dipole(i).dot(k_[n]);
			sfr+=mudotk*std::cos(prod);
			sfi+=mudotk*std::sin(prod);
		}
		//add to energy
		energy+=ka_[n]*(sfr*sfr+sfi*sfi);
	}
	//constant
	energy+=m2_*vc_;
	//return total
	const double pe=ke*0.5*energy;
	struc.pe()+=pe;
	return pe;
}

double Dipole::compute(Structure& struc)const{
	if(KSPACED_PRINT_FUNC>0) std::cout<<"KSpace::Dipole::compute(const Structure&)const:\n";
	const double ke=units::Consts::ke()*eps_;
	double energy=0;
	const int natoms=struc.nAtoms();
	c_.resize(natoms);
	s_.resize(natoms);
	//kspace - energy
	for(int n=0; n<k_.size(); ++n){
		//compute structure factor
		double sfr=0,sfi=0;
		for(int i=0; i<natoms; ++i){
			const double prod=k_[n].dot(struc.posn(i));
            const double mudotk=struc.dipole(i).dot(k_[n]);
			c_[i]=std::cos(prod);
			s_[i]=std::sin(prod);
			sfr+=mudotk*c_[i];
			sfi+=mudotk*s_[i];
		}
		//add to energy
		energy+=ka_[n]*(sfr*sfr+sfi*sfi);
		//add to force
		const Eigen::Vector3d kvec=ke*ka_[n]*k_[n];
		for(int i=0; i<natoms; ++i){
			struc.force(i).noalias()-=kvec*struc.dipole(i).dot(k_[n])*(c_[i]*sfi-s_[i]*sfr);
		}
	}
	//constant
	energy+=m2_*vc_;
	//return total
	const double pe=ke*0.5*energy;
	struc.pe()+=pe;
	return pe;
}

Eigen::Matrix3d& Dipole::interMat(const Eigen::Vector3d& drv, Eigen::Matrix3d& mat){
	const double ke=units::Consts::ke()*eps_;
	mat.setZero();
	for(int n=0; n<k_.size(); ++n){
		mat.noalias()-=std::cos(k_[n].dot(drv))*ka_[n]*k_[n]*k_[n].transpose();
	}
	mat*=ke;
	return mat;
}

}

//**********************************************
// serialization
//**********************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const KSpace::Dipole& obj){
		if(KSPACED_PRINT_FUNC>0) std::cout<<"nbytes(const KSpace::Dipole&):\n";
		int size=0;
		size+=nbytes(static_cast<const KSpace::Base&>(obj));
		return size;
	}
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const KSpace::Dipole& obj, char* arr){
		if(KSPACED_PRINT_FUNC>0) std::cout<<"pack(const KSpace::Dipole&,char*):\n";
		int pos=0;
		pos+=pack(static_cast<const KSpace::Base&>(obj),arr+pos);
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(KSpace::Dipole& obj, const char* arr){
		if(KSPACED_PRINT_FUNC>0) std::cout<<"unpack(KSpace::Dipole&,const char*):\n";
		int pos=0;
		pos+=unpack(static_cast<KSpace::Base&>(obj),arr+pos);
		return pos;
	}
	
}