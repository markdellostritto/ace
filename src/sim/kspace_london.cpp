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
#include "sim/kspace_london.hpp"

using math::constants::PI;
using math::constants::RadPI;
using math::constants::ZERO;

namespace KSpace{


//==== operators ====

std::ostream& operator<<(std::ostream& out, const London& l){
	return out<<static_cast<const Base&>(l);
}

//==== static functions ====

double London::fv(double pre, double a, double rc, double prec){
	const double b=a*rc;
	const double b2i=1.0/(b*b);
	const double a2=a*a;
	const double a4=a2*a2;
	const double a5=a4*a;
	return pre*a5*(b2i*(b2i+1.0)+0.5)*erfc(b)-prec;
}

double London::fd(double pre, double a, double rc, double prec){
	const double da=a*0.0001;
	const double fp=fv(pre,a+da,rc,prec);
	const double fm=fv(pre,a-da,rc,prec);
	return 0.5*(fp-fm)/da;
}

//==== member functions ====

void London::init(const Structure& struc){
	if(KSPACEL_PRINT_FUNC>0) std::cout<<"KSpace::London::init(const Structure&,const Eigen::MatrixXd&):\n";
	if(prec_<=0) throw std::invalid_argument("KSpace::London::init(const Structure&): invalid precision.");
	if(rc_<=0) throw std::invalid_argument("KSpace::London::init(const Structure&): invalid rcut.");
	
	//==== summations over c6 ====
	bs_.resize(b_.rows());
	for(int i=0; i<b_.rows(); ++i){
		bs_[i]=sqrt(b_(i,i));
	}
	/*double bsum=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const int ti=struc.type(i);
		bsum+=b_(ti,ti);
		for(int j=i+1; j<struc.nAtoms(); ++j){
			const int tj=struc.type(j);
			bsum+=2.0*b_(ti,tj);
		}
	}*/
	bijs_=b_.sum();
	biis_=b_.diagonal().sum();
	
	//==== set structural data ====
	const Eigen::Matrix3d& R=struc.R();
	const Eigen::Matrix3d& K=struc.K();
	const Eigen::Vector3d L=(
		Eigen::Vector3d()<<R.col(0).norm(),R.col(1).norm(),R.col(2).norm()
	).finished();
	
	//==== set real space error and alpha ====
	alpha_=(1.35-0.15*log(prec_))/rc_;
	errEr_=0;
	{
		bool error=false;
		const double pre=RadPI*RadPI*RadPI*bijs_/(struc.vol()*struc.nAtoms());
		const int max=10000;
		const double tol=1.0e-5;
		for(int i=0; i<max; ++i){
			const double da=fv(pre,alpha_,rc_,prec_)/fd(pre,alpha_,rc_,prec_);
			alpha_-=da;
			if(fabs(da) < tol) break;
			if(alpha_<0 || alpha_!=alpha_){
				error=true; break;
			}
		}
		if(error) alpha_=(1.35-0.15*log(prec_))/rc_;
		errEr_=fv(pre,alpha_,rc_,0.0);
	}
	a3_=alpha_*alpha_*alpha_;
	const double a6_=a3_*a3_;
	
	//==== set reciprocal space error and nk ====
	nk_.setZero();
	errEk_=0;
	for(int i=0; i<3; ++i){
		double errK=0;
		const double Kn=K.col(i).norm();
		const double pre=bijs_/(6.0*RadPI*struc.nAtoms())*a6_;
		do{
			nk_[i]++;
			const double kv=nk_[i]*Kn;
			const double arg=0.5*kv/alpha_;
			errK=pre*(2.0*arg*exp(-arg*arg)+RadPI*erfc(arg));
		}while(errK>prec_);
		if(errK>errEk_) errEk_=errK;
	}
	
	//==== compute k-vecs and k-amps ====
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
	ka_.resize(NK);
	const double c=RadPI*RadPI*RadPI/(12.0*struc.vol());
	for(int i=0; i<NK; ++i){
		const double kn=k_[i].norm();
		const double b=0.5*kn/alpha_;
		ka_[i]=c*kn*kn*kn*(RadPI*erfc(b)+(1.0/(2.0*b*b*b)-1.0/b)*exp(-b*b));
	}
	
	if(KSPACEL_PRINT_STATUS>0){
		std::cout<<"alpha = "<<alpha_<<"\n";
		std::cout<<"errEr = "<<errEr_<<"\n";
		std::cout<<"errEk = "<<errEk_<<"\n";
		std::cout<<"nk    = "<<nk_.transpose()<<"\n";
		std::cout<<"NK    = "<<NK<<"\n";
	}
	
}

double London::energy(Structure& struc)const{
	if(KSPACEL_PRINT_FUNC>0) std::cout<<"KSpace::London::energy(const Structure&)const:\n";
	double energy=0;
	const int natoms=struc.nAtoms();
	//kspace
	for(int n=0; n<k_.size(); ++n){
		double sfr=0,sfi=0;
		for(int i=0; i<natoms; ++i){
			const double prod=k_[n].dot(struc.posn(i));
			sfr+=bs_[struc.type(i)]*std::cos(prod);
			sfi+=bs_[struc.type(i)]*std::sin(prod);
		}
		energy+=ka_[n]*(sfr*sfr+sfi*sfi);
	}
	//constant
	energy+=a3_/3.0*(RadPI*RadPI*RadPI/struc.vol()*bijs_-a3_/2.0*biis_);
	//return total
	const double pe=-0.5*energy;
	struc.pe()+=pe;
	return pe;
}

double London::compute(Structure& struc)const{
	if(KSPACEL_PRINT_FUNC>0) std::cout<<"KSpace::London::compute(const Structure&)const:\n";
	double energy=0;
	const int natoms=struc.nAtoms();
	c_.resize(natoms);
	s_.resize(natoms);
	//kspace
	for(int n=0; n<k_.size(); ++n){
		//compute structure factor
		double sfr=0,sfi=0;
		for(int i=0; i<natoms; ++i){
			const double prod=k_[n].dot(struc.posn(i));
			c_[i]=std::cos(prod);
			s_[i]=std::sin(prod);
			sfr+=bs_[struc.type(i)]*c_[i];
			sfi+=bs_[struc.type(i)]*s_[i];
		}
		//add to energy
		energy+=ka_[n]*(sfr*sfr+sfi*sfi);
		//add to force
		for(int i=0; i<natoms; ++i){
			struc.force(i).noalias()-=bs_[struc.type(i)]*ka_[n]*k_[n]*(sfr*s_[i]-sfi*c_[i]);
		}
	}
	//constant
	energy+=a3_/3.0*(RadPI*RadPI*RadPI/struc.vol()*bijs_-a3_/2.0*biis_);
	//return total
	return -0.5*energy;
}

}

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const KSpace::London& obj){
		if(KSPACEL_PRINT_FUNC>0) std::cout<<"nbytes(const KSpace::London&):\n";
		int size=0;
		size+=nbytes(static_cast<const KSpace::Base&>(obj));
		size+=nbytes(obj.b());
		return size;
	}
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const KSpace::London& obj, char* arr){
		if(KSPACEL_PRINT_FUNC>0) std::cout<<"pack(const KSpace::London&,char*):\n";
		int pos=0;
		pos+=pack(static_cast<const KSpace::Base&>(obj),arr+pos);
		pos+=pack(obj.b(),arr+pos);
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(KSpace::London& obj, const char* arr){
		if(KSPACEL_PRINT_FUNC>0) std::cout<<"unpack(KSpace::London&,const char*):\n";
		int pos=0;
		Eigen::MatrixXd b;
		pos+=unpack(static_cast<KSpace::Base&>(obj),arr+pos);
		pos+=unpack(b,arr+pos);
		obj.b()=b;
		return pos;
	}
	
}