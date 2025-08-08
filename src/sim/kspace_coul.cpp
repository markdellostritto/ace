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
#include "sim/kspace_coul.hpp"

using math::constants::PI;
using math::constants::RadPI;
using math::constants::ZERO;

namespace KSpace{

//==== operators ====

std::ostream& operator<<(std::ostream& out, const Coul& c){
	return out<<static_cast<const Base&>(c);
}

//==== member functions ====

void Coul::init(const Structure& struc){
	if(KSPACEC_PRINT_FUNC>0) std::cout<<"KSpace::Coul::init(const Structure&):\n";
	if(prec_<=0) throw std::invalid_argument("KSpace::Coul::init(const Structure&): invalid precision.");
	if(rc_<=0) throw std::invalid_argument("KSpace::Coul::init(const Structure&): invalid rcut.");
	
	//set structural data
	const int& N=struc.nAtoms();
	const Eigen::Matrix3d& R=struc.R();
	const Eigen::Matrix3d& K=struc.K();
	const Eigen::Vector3d L=(Eigen::Vector3d()<<R.col(0).norm(),R.col(1).norm(),R.col(2).norm()).finished();
	//Eigen::Vector3d L; L<<R.col(0).norm(),R.col(1).norm(),R.col(2).norm();
	q2_=0;
	for(int i=0; i<N; ++i) q2_+=struc.charge(i)*struc.charge(i);
	if(q2_==0) throw std::invalid_argument("KSpace::Coul::init(const Structure&): zero abs charge.");
	
	//set constants
	alpha_=prec_*sqrt(N*rc_*L.prod())/(2.0*q2_);
	if(alpha_!=alpha_) throw std::invalid_argument("Invalid alpha\n");
	if(alpha_>=1.0) alpha_=(1.35-0.15*log(prec_))/rc_;
	else alpha_=sqrt(-log(alpha_))/rc_;
	const double a2i=1.0/(alpha_*alpha_);
	if(KSPACEC_PRINT_DATA>0) std::cout<<"alpha_ = "<<alpha_<<"\n";
	
	//set lattice vectors - reciprocal space
	nk_.setZero();
	for(int i=0; i<3; ++i){
		double err=0;
		const double pre=2.0*q2_*alpha_/(L[i]*L[i]*RadPI*sqrt(N));
		const double expf=exp(-PI*PI*a2i);
		do{
			++nk_[i];
			const double kv=(K.col(i)*nk_[i]).norm();
			err=pre*pow(expf,kv*kv)/sqrt(kv);
		}while(err>prec_);
	}
	nk_*=5;//needed to match LAMMPS
	if(KSPACEC_PRINT_DATA>0) std::cout<<"nk_ = "<<nk_.transpose()<<"\n";
	
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
	
	vc_=-2.0*alpha_/RadPI;
}

double Coul::energy(Structure& struc)const{
	if(KSPACEC_PRINT_FUNC>0) std::cout<<"KSpace::Coul::energy(const Structure&)const:\n";
	const double ke=units::Consts::ke()*eps_;
	double energy=0;
	const int natoms=struc.nAtoms();
	//kspace - energy
	for(int n=0; n<k_.size(); ++n){
		//compute structure factor
		double sfr=0,sfi=0;
		for(int i=0; i<natoms; ++i){
			const double prod=k_[n].dot(struc.posn(i));
			sfr+=struc.charge(i)*std::cos(prod);
			sfi+=struc.charge(i)*std::sin(prod);
		}
		//add to energy
		energy+=ka_[n]*(sfr*sfr+sfi*sfi);
	}
	//constant
	energy+=econst_?q2_*vc_:0.0;
	//net charge
	double qtot=0;
	for(int i=0; i<natoms; ++i){
		qtot+=struc.charge(i);
	}
	energy-=qtot*qtot*PI/(alpha_*alpha_*struc.vol());
	//return total
	const double pe=ke*0.5*energy;
	struc.pe()+=pe;
	return pe;
}

double Coul::compute(Structure& struc)const{
	if(KSPACEC_PRINT_FUNC>0) std::cout<<"KSpace::Coul::compute(const Structure&)const:\n";
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
			c_[i]=std::cos(prod);
			s_[i]=std::sin(prod);
			sfr+=struc.charge(i)*c_[i];
			sfi+=struc.charge(i)*s_[i];
		}
		//add to energy
		energy+=ka_[n]*(sfr*sfr+sfi*sfi);
		//add to force
		const Eigen::Vector3d kvec=ke*ka_[n]*k_[n];
		for(int i=0; i<natoms; ++i){
			struc.force(i).noalias()-=kvec*struc.charge(i)*(c_[i]*sfi-s_[i]*sfr);
		}
	}
	//constant term
	energy+=econst_?q2_*vc_:0.0;
	//net charge
	double qtot=0;
	for(int i=0; i<natoms; ++i){
		qtot+=struc.charge(i);
	}
	energy-=qtot*qtot*PI/(alpha_*alpha_*struc.vol());
	//return total
	const double pe=ke*0.5*energy;
	struc.pe()+=pe;
	return pe;
}

double Coul::compute_explicit(Structure& struc)const{
	if(KSPACEC_PRINT_FUNC>0) std::cout<<"KSpace::Coul::compute_explicit(const Structure&)const:\n";
	const double ke=units::Consts::ke()*eps_;
	double energy=0;
	const int natoms=struc.nAtoms();
	//kspace - energy
	for(int n=0; n<k_.size(); ++n){
		//compute structure factor
		std::complex<double> sf(0,0);
		for(int i=0; i<natoms; ++i){
			sf+=struc.charge(i)*exp(I_*k_[n].dot(struc.posn(i)));
		}
		//add to energy
		energy+=ka_[n]*std::norm(sf);
		//add to force
		for(int i=0; i<natoms; ++i){
			std::complex<double> y=I_*exp(I_*k_[n].dot(struc.posn(i)))*std::conj(sf);
			struc.force(i).noalias()-=0.5*ke*struc.charge(i)*ka_[n]*k_[n]*(y+std::conj(y)).real();
		}
	}
	//constant term
	energy+=q2_*vc_;
	//net charge
	double qtot=0;
	for(int i=0; i<natoms; ++i){
		qtot+=struc.charge(i);
	}
	energy-=qtot*qtot*PI/(alpha_*alpha_*struc.vol());
	//return total
	struc.pe()+=ke*0.5*energy;
	return ke*0.5*energy;
}

double Coul::compute_brute(Structure& struc)const{
	if(KSPACEC_PRINT_FUNC>0) std::cout<<"KSpace::Coul::compute_brute(const Structure&)const:\n";
	const double ke=units::Consts::ke()*eps_;
	double energy=0;
	const int natoms=struc.nAtoms();
	//kspace - energy
	for(int n=0; n<k_.size(); ++n){
		std::complex<double> sf(0,0);
		for(int i=0; i<natoms; ++i){
			sf+=struc.charge(i)*exp(-I_*k_[n].dot(struc.posn(i)));
		}
		energy+=ka_[n]*std::norm(sf);
	}
	//kspace - force
	Eigen::Vector3d rij;
	for(int i=0; i<natoms; ++i){
		const double qi=struc.charge(i);
		for(int j=i+1; j<natoms; ++j){
			const double qj=struc.charge(j);
			struc.diff(struc.posn(i),struc.posn(j),rij);
			Eigen::Vector3d ssin=Eigen::Vector3d::Zero();
			for(int n=0; n<k_.size(); ++n){
				ssin.noalias()+=ka_[n]*k_[n]*sin(k_[n].dot(rij));
			}
			ssin*=ke*qi*qj;
			struc.force(i).noalias()+=ssin;
			struc.force(j).noalias()-=ssin;
		}
	}
	//constant
	energy+=q2_*vc_;
	//net charge
	double qtot=0;
	for(int i=0; i<natoms; ++i){
		qtot+=struc.charge(i);
	}
	energy-=qtot*qtot*PI/(alpha_*alpha_*struc.vol());
	//return total
	struc.pe()+=ke*0.5*energy;
	return ke*0.5*energy;
}

}

//**********************************************
// serialization
//**********************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const KSpace::Coul& obj){
		if(KSPACEC_PRINT_FUNC>0) std::cout<<"nbytes(const KSpace::Coul&):\n";
		int size=0;
		size+=nbytes(static_cast<const KSpace::Base&>(obj));
		return size;
	}
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const KSpace::Coul& obj, char* arr){
		if(KSPACEC_PRINT_FUNC>0) std::cout<<"pack(const KSpace::Coul&,char*):\n";
		int pos=0;
		pos+=pack(static_cast<const KSpace::Base&>(obj),arr+pos);
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(KSpace::Coul& obj, const char* arr){
		if(KSPACEC_PRINT_FUNC>0) std::cout<<"unpack(KSpace::Coul&,const char*):\n";
		int pos=0;
		pos+=unpack(static_cast<KSpace::Base&>(obj),arr+pos);
		return pos;
	}
	
}