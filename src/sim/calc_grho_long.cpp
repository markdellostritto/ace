// math
#include "chem/units.hpp"
// math
#include "math/const.hpp"
#include "math/special.hpp"
// string
#include "str/string.hpp"
// sim
#include "sim/calc_grho_long.hpp"

using math::constants::RadPI;
using math::constants::ZERO;
using math::special::fmexp;

//==== contructors/destructors ====

CalcGRhoLong::CalcGRhoLong(double rc, Calculator::Mix mix):Calculator(Calculator::Name::GRHO_LONG,rc){
    mix_=Calculator::Mix::HARMONIC;
	mix_=mix;
    if(mix_==Calculator::Mix::NONE) throw std::invalid_argument("CalcGRhoLong::CalcGRhoLong(double,Calculator::Mix): Invalid mixing.\n");
}

//==== operator ====

std::ostream& operator<<(std::ostream& out, const CalcGRhoLong& calc){
	return out<<static_cast<const Calculator&>(calc)<<
        " mix "<<calc.mix_<<" eps "<<calc.eps_<<" prec "<<calc.prec_;
}

//==== member functions ====

void CalcGRhoLong::resize(int ntypes){
    if(CALC_GRHO_LONG_PRINT_FUNC>0) std::cout<<"CalcGRhoLong::resize(int):\n";
    Calculator::resize(ntypes);
	if(ntypes_>0){
        radius_=Eigen::VectorXd::Zero(ntypes_);
		gamma_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
        rgamma_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
	}
}

void CalcGRhoLong::init(){
    if(CALC_GRHO_LONG_PRINT_FUNC>0) std::cout<<"CalcGRhoLong::init():\n";
    Eigen::VectorXd alpha=Eigen::VectorXd::Zero(ntypes_);
    for(int i=0; i<ntypes_; ++i){
        alpha[i]=1.0/(2.0*radius_[i]*radius_[i]);
    }
    for(int i=0; i<ntypes_; ++i){
		for(int j=0; j<ntypes_; ++j){
            gamma_(i,j) = 2.0*alpha[i]*alpha[j]/(alpha[i]+alpha[j]);
            rgamma_(i,j) = std::sqrt(0.5*gamma_(i,j));
		}
	}
    coul_.prec()=prec_;
	coul_.rc()=rc_;
}

void CalcGRhoLong::read(Token& token){
    static_cast<Calculator&>(*this).read(token);
    prec_=std::atof(token.next().c_str());
    mix_=Calculator::Mix::read(string::to_upper(token.next()).c_str());
    if(prec_<=0) throw std::invalid_argument("CalcGRhoLong::read(Token&): Invalid prec.");
    if(eps_<=0) throw std::invalid_argument("CalcGRhoLong::read(Token&): Invalid eps.");
    if(mix_==Calculator::Mix::NONE) throw std::invalid_argument("CalcGRhoLong::read(Token&): Invalid mixing type.");
}

void CalcGRhoLong::coeff(Token& token){
    if(CALC_GRHO_LONG_PRINT_FUNC>0) std::cout<<"CalcGRhoLong::coeff(Token&):\n";
	//coeff type radius 
	const int type=std::atoi(token.next().c_str())-1;
	const double radius=std::atof(token.next().c_str());
    //check type
	if(type>=ntypes_) throw std::invalid_argument("Invalid type.");
    if(radius<=0) throw std::invalid_argument("Invalid radius.");
    //assign
    radius_[type]=radius;
}

void CalcGRhoLong::init(const Structure& struc){
	coul_.init(struc);
}

double CalcGRhoLong::energy(Structure& struc, const NeighborList& nlist)const{
    if(CALC_GRHO_LONG_PRINT_FUNC>0) std::cout<<"CalcGRhoLong::energy(const Structure&,const NeighborList&):\n";
    const double ke=units::Consts::ke()*eps_;
    // k-space
	//coul_.init(struc);
	const double energyK=coul_.energy(struc);
	const double a=coul_.alpha();
    // r-space
    double energyR=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const int ti=struc.type(i);
        const double qi=struc.charge(i);
        for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
			const int tj=struc.type(jj);
            const double qj=struc.charge(jj);
            const double dr2=(struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j))).squaredNorm();
			if(dr2<rc2_){
				const double dr=sqrt(dr2);
                double eCoul=0.0;
                if(dr>ZERO){
					const double pf=ke*qi*qj/dr;
                	const double ferfg=std::erf(rgamma_(ti,tj)*dr);
                    const double ferfp=std::erf(a*dr);
                    eCoul=pf*(ferfg-ferfp);
				} else eCoul=ke*qi*qj*2.0/RadPI*(rgamma_(ti,tj)-a);
                //compute energy
				energyR+=eCoul;
			}
		}
	}
	energyR*=0.5;
    struc.pe()+=energyR;
	//return total energy
	return energyR+energyK;
}

double CalcGRhoLong::compute(Structure& struc, const NeighborList& nlist)const{
    if(CALC_GRHO_LONG_PRINT_FUNC>0) std::cout<<"CalcGRhoLong::compute(const Structure&,const NeighborList&):\n";
	const double ke=units::Consts::ke()*eps_;
    // k-space
    //coul_.init(struc);
    const double energyK=coul_.compute(struc);
    const double a=coul_.alpha();
    // r-space
	double energyR=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const int ti=struc.type(i);
        const double qi=struc.charge(i);
        for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
			const int tj=struc.type(jj);
            const double qj=struc.charge(jj);
            const Eigen::Vector3d drv=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
            const double dr2=drv.squaredNorm();
			if(dr2<rc2_){
				const double dr=sqrt(dr2);
                double eCoul=0.0;
                double fCoul=0.0;
				if(dr>ZERO){
					const double pf=ke*qi*qj/dr;
                	const double ferfg=std::erf(rgamma_(ti,tj)*dr);
                    const double ferfp=std::erf(a*dr);
                    eCoul=pf*(ferfg-ferfp);
                    fCoul=pf/dr2*(
                        (ferfg-ferfp)
                        +2.0/RadPI*dr*(
                            -rgamma_(ti,tj)*fmexp(-0.5*gamma_(ti,tj)*dr2)
                            +a*fmexp(-a*a*dr2)
                        )
                    );
				} else eCoul=ke*qi*qj*2.0/RadPI*(rgamma_(ti,tj)-a);
                //compute energy
				energyR+=eCoul;
				//compute force
				struc.force(i).noalias()+=fCoul*drv;
			}
		}
	}
    energyR*=0.5;
	struc.pe()+=energyR;
	//return total energy
	return energyR+energyK;
}

//**********************************************
// serialization
//**********************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const CalcGRhoLong& obj){
		if(CALC_GRHO_LONG_PRINT_FUNC>0) std::cout<<"nbytes(const CalcGRhoLong&):\n";
		int size=0;
		const int nt=obj.ntypes();
		size+=nbytes(static_cast<const Calculator&>(obj));
        size+=sizeof(obj.mix());//mix_
		size+=sizeof(int);//ntypes_
        size+=sizeof(double);//eps_
        size+=sizeof(double);//prec_
        size+=nt*sizeof(double);//radius_
		return size;
	}
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcGRhoLong& obj, char* arr){
		if(CALC_GRHO_LONG_PRINT_FUNC>0) std::cout<<"pack(const CalcGRhoLong&,char*):\n";
		int pos=0;
		pos+=pack(static_cast<const Calculator&>(obj),arr+pos);
        std::memcpy(arr+pos,&obj.mix(),sizeof(obj.mix())); pos+=sizeof(obj.mix());//mix_
		std::memcpy(arr+pos,&obj.ntypes(),sizeof(int)); pos+=sizeof(int);//ntypes_
        std::memcpy(arr+pos,&obj.eps(),sizeof(double)); pos+=sizeof(double);//eps_
        std::memcpy(arr+pos,&obj.prec(),sizeof(double)); pos+=sizeof(double);//prec_
		const int nt=obj.ntypes();
		if(nt>0){
            std::memcpy(arr+pos,obj.radius().data(),nt*sizeof(double)); pos+=nt*sizeof(double);//radius_
		}
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcGRhoLong& obj, const char* arr){
		if(CALC_GRHO_LONG_PRINT_FUNC>0) std::cout<<"unpack(CalcGRhoLong&,const char*):\n";
		int pos=0,nt=0;
		pos+=unpack(static_cast<Calculator&>(obj),arr+pos);
		if(obj.name()!=Calculator::Name::GRHO_LONG) throw std::invalid_argument("serialize::unpack(CalcGRhoLong&,const char*): Invalid name.");
        std::memcpy(&obj.mix(),arr+pos,sizeof(obj.mix())); pos+=sizeof(obj.mix());//mix_
		std::memcpy(&nt,arr+pos,sizeof(int)); pos+=sizeof(int);//ntypes_
		std::memcpy(&obj.eps(),arr+pos,sizeof(double)); pos+=sizeof(double);//eps_
        std::memcpy(&obj.prec(),arr+pos,sizeof(double)); pos+=sizeof(double);//prec_
		obj.resize(nt);
		if(nt>0){
            std::memcpy(obj.radius().data(),arr+pos,nt*sizeof(double)); pos+=nt*sizeof(double);//radius_
		}
		obj.init();
		return pos;
	}
	
}