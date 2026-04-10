// math
#include "chem/units.hpp"
// math
#include "math/const.hpp"
#include "math/special.hpp"
// string
#include "str/string.hpp"
// sim
#include "sim/calc_pauli.hpp"

using math::constants::RadPI;
using math::constants::ZERO;
using math::special::fmexp;

//==== operator ====

std::ostream& operator<<(std::ostream& out, const CalcPauli& calc){
	return out<<static_cast<const Calculator&>(calc)<<" eps "<<calc.eps_;
}

//==== member functions ====

void CalcPauli::resize(int ntypes){
    if(CALC_PAULI_PRINT_FUNC>0) std::cout<<"CalcPauli::resize(int):\n";
	Calculator::resize(ntypes);
	if(ntypes_>0){
        z_=Eigen::VectorXd::Zero(ntypes_);
        radius_=Eigen::VectorXd::Zero(ntypes_);
		gamma_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
		amp_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
	}
}

void CalcPauli::init(){
    if(CALC_PAULI_PRINT_FUNC>0) std::cout<<"CalcPauli::init():\n";
    Eigen::VectorXd alpha=Eigen::VectorXd::Zero(ntypes_);
    for(int i=0; i<ntypes_; ++i){
        alpha[i]=1.0/(2.0*radius_[i]*radius_[i]);
    }
    for(int i=0; i<ntypes_; ++i){
		for(int j=0; j<ntypes_; ++j){
            gamma_(i,j)=2.0*alpha[i]*alpha[j]/(alpha[i]+alpha[j]);
            amp_(i,j)=std::pow(2.0*sqrt(alpha[i]*alpha[j])/(alpha[i]+alpha[j]),3.0);
		}
	}
}

void CalcPauli::coeff(Token& token){
    if(CALC_PAULI_PRINT_FUNC>0) std::cout<<"CalcPauli::coeff(Token&):\n";
	//coeff type radius 
	const int type=std::atoi(token.next().c_str())-1;
	const double radius=std::atof(token.next().c_str());
    const double z=std::atof(token.next().c_str());
    //check type
	if(type>=ntypes_) throw std::invalid_argument("Invalid type.");
    if(radius<=0) throw std::invalid_argument("Invalid radius.");
    //assign
    radius_[type]=radius;
    z_[type]=z;
}

double CalcPauli::energy(Structure& struc, const NeighborList& nlist)const{
    if(CALC_PAULI_PRINT_FUNC>0) std::cout<<"CalcPauli::energy(const Structure&,const NeighborList&):\n";
    const double ke=units::Consts::ke()*eps_;
	double energy=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const int ti=struc.type(i);
        const double qi=struc.charge(i);
        for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
			const int tj=struc.type(jj);
            const double qj=struc.charge(jj);
            const double dr2=(struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j))).squaredNorm();
			if(dr2<rc2_){
				energy+=z_[ti]*z_[tj]*amp_(ti,tj)*fmexp(-0.5*gamma_(ti,tj)*dr2)/sqrt(dr2);
			}
		}
	}
	energy*=0.5;
    struc.pe()+=energy;
	return energy;
}

double CalcPauli::energy(Structure& struc)const{
    if(CALC_PAULI_PRINT_FUNC>0) std::cout<<"CalcPauli::energy(const Structure&):\n";
    const double ke=units::Consts::ke()*eps_;
	double energy=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const int ti=struc.type(i);
        const double qi=struc.charge(i);
        for(int j=i+1; j<struc.nAtoms(); ++j){
			const int tj=struc.type(j);
            const double qj=struc.charge(j);
            const double dr2=(struc.posn(i)-struc.posn(j)).squaredNorm();
			if(dr2<rc2_){
				energy+=z_[ti]*z_[tj]*amp_(ti,tj)*fmexp(-0.5*gamma_(ti,tj)*dr2)/sqrt(dr2);
			}
		}
	}
	struc.pe()+=energy;
	return energy;
}

double CalcPauli::compute(Structure& struc, const NeighborList& nlist)const{
    if(CALC_PAULI_PRINT_FUNC>0) std::cout<<"CalcPauli::compute(const Structure&,const NeighborList&):\n";
    const double ke=units::Consts::ke()*eps_;
	double energy=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		//store atom data
		const int ti=struc.type(i);
        const double qi=struc.charge(i);
        for(int j=0; j<nlist.size(i); ++j){
			//store atom data
			const int jj=nlist.index(i,j);
			const int tj=struc.type(jj);
            const double qj=struc.charge(jj);
            //compute distance
			const Eigen::Vector3d drv=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
            const double dr2=drv.squaredNorm();
			//compute interaction
			if(dr2<rc2_){
				const double dr=sqrt(dr2);
                const double amp=z_[ti]*z_[tj]*amp_(ti,tj);
                const double gdr2=gamma_(ti,tj)*dr2;
                const double fexp=fmexp(-0.5*gdr2);
                energy+=amp*fexp/dr;
                struc.force(i).noalias()+=amp*fexp*(gdr2+1.0)/(dr2*dr)*drv;
			}
		}
	}
	energy*=0.5;
    struc.pe()+=energy;
	return energy;
}

double CalcPauli::compute(Structure& struc)const{
    if(CALC_PAULI_PRINT_FUNC>0) std::cout<<"CalcPauli::compute(const Structure&):\n";
    const double ke=units::Consts::ke()*eps_;
	double energy=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		//store atom data
		const int ti=struc.type(i);
        const double qi=struc.charge(i);
        for(int j=i+1; j<struc.nAtoms(); ++j){
			//store atom data
			const int tj=struc.type(j);
            const double qj=struc.charge(j);
            //compute distance
			const Eigen::Vector3d drv=struc.posn(i)-struc.posn(j);
            const double dr2=drv.squaredNorm();
			//compute interaction
			if(dr2<rc2_){
				const double dr=sqrt(dr2);
                const double amp=z_[ti]*z_[tj]*amp_(ti,tj);
                const double gdr2=gamma_(ti,tj)*dr2;
                const double fexp=fmexp(-0.5*gdr2);
                energy+=amp*fexp/dr;
                const Eigen::Vector3d fij=amp*fexp*(gdr2+1.0)/(dr2*dr)*drv;
				struc.force(i).noalias()+=fij;
				struc.force(j).noalias()-=fij;
			}
		}
	}
	struc.pe()+=energy;
	return energy;
}

//**********************************************
// serialization
//**********************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const CalcPauli& obj){
		if(CALC_PAULI_PRINT_FUNC>0) std::cout<<"nbytes(const CalcPauli&):\n";
		int size=0;
		const int nt=obj.ntypes();
		size+=nbytes(static_cast<const Calculator&>(obj));
		size+=sizeof(int);//ntypes_
        size+=sizeof(double);//eps_
		size+=nt*sizeof(double);//z_
        size+=nt*sizeof(double);//radius_
		return size;
	}
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcPauli& obj, char* arr){
		if(CALC_PAULI_PRINT_FUNC>0) std::cout<<"pack(const CalcPauli&,char*):\n";
		int pos=0;
		pos+=pack(static_cast<const Calculator&>(obj),arr+pos);
		std::memcpy(arr+pos,&obj.ntypes(),sizeof(int)); pos+=sizeof(int);//ntypes_
        std::memcpy(arr+pos,&obj.eps(),sizeof(double)); pos+=sizeof(double);//eps_
		const int nt=obj.ntypes();
		if(nt>0){
            std::memcpy(arr+pos,obj.z().data(),nt*sizeof(double)); pos+=nt*sizeof(double);//z_
            std::memcpy(arr+pos,obj.radius().data(),nt*sizeof(double)); pos+=nt*sizeof(double);//radius_
		}
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcPauli& obj, const char* arr){
		if(CALC_PAULI_PRINT_FUNC>0) std::cout<<"unpack(CalcPauli&,const char*):\n";
		int pos=0,nt=0;
		pos+=unpack(static_cast<Calculator&>(obj),arr+pos);
		if(obj.name()!=Calculator::Name::GRHO_CUT) throw std::invalid_argument("serialize::unpack(CalcPauli&,const char*): Invalid name.");
		std::memcpy(&nt,arr+pos,sizeof(int)); pos+=sizeof(int);//ntypes_
		std::memcpy(&obj.eps(),arr+pos,sizeof(double)); pos+=sizeof(double);//eps_
        obj.resize(nt);
		if(nt>0){
            std::memcpy(obj.z().data(),arr+pos,nt*sizeof(double)); pos+=nt*sizeof(double);//z_
            std::memcpy(obj.radius().data(),arr+pos,nt*sizeof(double)); pos+=nt*sizeof(double);//radius_
		}
		obj.init();
		return pos;
	}
	
}