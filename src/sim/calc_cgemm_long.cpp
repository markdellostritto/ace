// math
#include "chem/units.hpp"
// math
#include "math/const.hpp"
#include "math/special.hpp"
// string
#include "str/string.hpp"
// sim
#include "sim/calc_cgemm_long.hpp"

using math::constants::PI;
using math::constants::RadPI;
using math::constants::ZERO;
using math::special::fmexp;

//==== contructors/destructors ====

CalcCGemmLong::CalcCGemmLong(double rc, double lambdaC, double lambdaS):Calculator(Calculator::Name::CGEMM_LONG,rc){
    mix_=Calculator::Mix::HARMONIC;
	lambdaC_=lambdaC;
    lambdaS_=lambdaS;
    if(lambdaC_<=0) throw std::invalid_argument("CalcCGemmLong::CalcCGemmLong(double,double,double): Invalid lambdaC\n");
    if(lambdaS_<=0) throw std::invalid_argument("CalcCGemmLong::CalcCGemmLong(double,double,double): Invalid lambdaS\n");
}

//==== operator ====

std::ostream& operator<<(std::ostream& out, const CalcCGemmLong& calc){
	return out<<static_cast<const Calculator&>(calc)
        <<" mix "<<calc.mix_<<" lambdaC "<<calc.lambdaC_<<" lambdaS "
        <<calc.lambdaS_<<" rRep "<<calc.rRep_
        <<" eps "<<calc.eps_<<" prec "<<calc.prec_;
}

//==== member functions ====

void CalcCGemmLong::resize(int ntypes){
    if(CALC_CGEMM_LONG_PRINT_FUNC>0) std::cout<<"CalcCGemmLong::resize(int):\n";
    Calculator::resize(ntypes);
	if(ntypes_>0){
        radius_=Eigen::VectorXd::Zero(ntypes_);
		aOver_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
		aRep_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
		gammaS_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
		gammaC_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
        rgammaC_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
	}
}

void CalcCGemmLong::init(){
    if(CALC_CGEMM_LONG_PRINT_FUNC>0) std::cout<<"CalcCGemmLong::init():\n";
    Eigen::VectorXd alphaS=Eigen::VectorXd::Zero(ntypes_);
    Eigen::VectorXd alphaC=Eigen::VectorXd::Zero(ntypes_);
    for(int i=0; i<ntypes_; ++i){
        alphaS[i]=lambdaS_/(2.0*radius_[i]*radius_[i]);
        alphaC[i]=lambdaC_/(2.0*radius_[i]*radius_[i]);
    }
    for(int i=0; i<ntypes_; ++i){
		for(int j=0; j<ntypes_; ++j){
            switch(mix_){
                case Calculator::Mix::ARITHMETIC:
                    aOver_(i,j) = 0.5*(aOver_(i,i)+aOver_(j,j));
                    aRep_(i,j) = 0.5*(aRep_(i,i)+aRep_(j,j));
                break;
                case Calculator::Mix::GEOMETRIC:
			        aOver_(i,j) = std::sqrt(aOver_(i,i)*aOver_(j,j));
                    aRep_(i,j) = std::sqrt(aRep_(i,i)*aRep_(j,j));
                break;
                case Calculator::Mix::HARMONIC:
                    if(aOver_(i,i)<ZERO && aOver_(j,j)<ZERO) aOver_(i,j)=0.0;
                    else aOver_(i,j) = 2.0*(aOver_(i,i)*aOver_(j,j))/(aOver_(i,i)+aOver_(j,j));
                    if(aRep_(i,i)<ZERO && aRep_(j,j)<ZERO) aRep_(i,j)=0.0;
                    else aRep_(i,j) = 2.0*(aRep_(i,i)*aRep_(j,j))/(aRep_(i,i)+aRep_(j,j));
                break;
            }
			gammaS_(i,j) = 2.0*alphaS[i]*alphaS[j]/(alphaS[i]+alphaS[j]);
            gammaC_(i,j) = 2.0*alphaC[i]*alphaC[j]/(alphaC[i]+alphaC[j]);
            rgammaC_(i,j) = std::sqrt(0.5*gammaC_(i,j));
		}
	}
    coul_.prec()=prec_;
	coul_.rc()=rc_;
}

void CalcCGemmLong::read(Token& token){
    static_cast<Calculator&>(*this).read(token);
    lambdaC_=std::atof(token.next().c_str());
    lambdaS_=std::atof(token.next().c_str());
    eps_=std::atof(token.next().c_str());
    prec_=std::atof(token.next().c_str());
    if(!token.end()){
        rRep_=std::atof(token.next().c_str());
    }
    if(!token.end()){
        mix_=Calculator::Mix::read(string::to_upper(token.next()).c_str());
    }
    if(lambdaC_<=0) throw std::invalid_argument("CalcCGemmLong::read(Token&): Invalid lambdaC.");
    if(lambdaS_<=0) throw std::invalid_argument("CalcCGemmLong::read(Token&): Invalid lambdaS.");
    if(prec_<=0) throw std::invalid_argument("CalcCGemmLong::read(Token&): Invalid prec.");
    if(eps_<=0) throw std::invalid_argument("CalcCGemmLong::read(Token&): Invalid eps.");
    if(rRep_<=0) throw std::invalid_argument("CalcCGemmLong::read(Token&): Invalid rRep.");
    if(mix_==Calculator::Mix::NONE) throw std::invalid_argument("CalcCGemmLong::read(Token&): Invalid mixing type.");
}

void CalcCGemmLong::coeff(Token& token){
    if(CALC_CGEMM_LONG_PRINT_FUNC>0) std::cout<<"CalcCGemmLong::coeff(Token&):\n";
	//coeff lj_cut type radius amplitude
	const int type=std::atoi(token.next().c_str())-1;
	const double radius=std::atof(token.next().c_str());
    const double aOver=std::atof(token.next().c_str());
	const double aRep=std::atof(token.next().c_str());
	//check type
	if(type>=ntypes_) throw std::invalid_argument("Invalid type.");
    if(radius<=0) throw std::invalid_argument("Invalid radius.");
    if(aOver<0) throw std::invalid_argument("Invalid overlap amplitude.");
	if(aRep<0) throw std::invalid_argument("Invalid repulsive amplitude.");
    //assign
    radius_[type]=radius;
    aOver_(type,type)=aOver;
	aRep_(type,type)=aRep;
}

void CalcCGemmLong::init(const Structure& struc){
	coul_.init(struc);
}

double CalcCGemmLong::energy(Structure& struc, const NeighborList& nlist)const{
    if(CALC_CGEMM_LONG_PRINT_FUNC>0) std::cout<<"CalcCGemmLong::energy(const Structure&,const NeighborList&):\n";
    const double ke=units::Consts::ke()*eps_;
    // k-space
	//coul_.init(struc);
	const double energyK=coul_.energy(struc);
	const double a=coul_.alpha();
    // r-space
    const double cRep_=1.0/rRep_;
    double energyR=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const int ti=struc.type(i);
        const double qi=struc.charge(i);
        const double zi=std::fabs(qi);
        for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
			const int tj=struc.type(jj);
            const double qj=struc.charge(jj);
            const double zj=std::fabs(qj);
			const double dr2=(struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j))).squaredNorm();
			if(dr2<rc2_){
				const double dr=sqrt(dr2);
                // Coulomb
				double eCoul=0.0;
                if(dr>ZERO){
					const double pf=ke*qi*qj/dr;
                	const double ferfg=std::erf(rgammaC_(ti,tj)*dr);
                    const double ferfp=std::erf(a*dr);
                    eCoul=pf*(ferfg-ferfp);
				} else eCoul=ke*qi*qj*2.0/RadPI*(rgammaC_(ti,tj)-a);
                // overlap
                const double eOver=aOver_(ti,tj)*zi*zj*fmexp(-0.5*gammaS_(ti,tj)*dr2);
				// repulsion
				const double eRep=aRep_(ti,tj)*fmexp(-cRep_*dr);
                //const double eRep=aRep_(ti,tj)*cRep_*fmexp(-cRep_*dr);
				//compute energy
				energyR+=eCoul+eOver+eRep;
			}
		}
	}
	energyR*=0.5;
    struc.pe()+=energyR;
	//return total energy
	return energyR+energyK;
}

double CalcCGemmLong::compute(Structure& struc, const NeighborList& nlist)const{
    if(CALC_CGEMM_LONG_PRINT_FUNC>0) std::cout<<"CalcCGemmLong::compute(const Structure&,const NeighborList&):\n";
	const double ke=units::Consts::ke()*eps_;
    // k-space
    //coul_.init(struc);
    const double energyK=coul_.compute(struc);
    const double a=coul_.alpha();
    // r-space
    const double cRep_=1.0/rRep_;
    double energyR=0;
	//double energyCoul=0;
	//double energyOver=0;
	//double energyRep=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const int ti=struc.type(i);
        const double qi=struc.charge(i);
        const double zi=std::fabs(qi);
		for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
			const int tj=struc.type(jj);
            const double qj=struc.charge(jj);
            const double zj=std::fabs(qj);
			const Eigen::Vector3d drv=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
            const double dr2=drv.squaredNorm();
			if(dr2<rc2_){
				const double dr=sqrt(dr2);
                // Coulomb
				double eCoul=0.0;
                double fCoul=0.0;
				if(dr>ZERO){
					const double pf=ke*qi*qj/dr;
                	const double ferfg=std::erf(rgammaC_(ti,tj)*dr);
                    const double ferfp=std::erf(a*dr);
                    eCoul=pf*(ferfg-ferfp);
                    fCoul=pf/dr2*(
                        (ferfg-ferfp)
                        +2.0/RadPI*dr*(
                            -rgammaC_(ti,tj)*fmexp(-0.5*gammaC_(ti,tj)*dr2)
                            +a*fmexp(-a*a*dr2)
                        )
                    );
				} else eCoul=ke*qi*qj*2.0/RadPI*(rgammaC_(ti,tj)-a);
                // overlap
                const double eOver=aOver_(ti,tj)*zi*zj*fmexp(-0.5*gammaS_(ti,tj)*dr2);
				const double fOver=gammaS_(ti,tj)*eOver;
				// repulsion
				const double eRep=aRep_(ti,tj)*fmexp(-cRep_*dr);
				//const double eRep=aRep_(ti,tj)*cRep_*fmexp(-cRep_*dr);
				const double fRep=cRep_*eRep;
				//compute energy
				energyR+=eCoul+eOver+eRep;
				//energyCoul+=eCoul;
				//energyOver+=eOver;
				//energyRep+=eRep;
                //compute force
				struc.force(i).noalias()+=(fCoul+fOver+fRep)*drv;
			}
		}
	}
    energyR*=0.5;
	//energyCoul*=0.5;
	//energyOver*=0.5;
	//energyRep*=0.5;
    struc.pe()+=energyR;
	//return total energy
	//std::cout<<"eKspace "<<energyK<<" eCoul "<<energyCoul<<" eOver "<<energyOver<<" eRep "<<energyRep<<" alpha "<<a<<"\n";
	return energyR+energyK;
}

//**********************************************
// serialization
//**********************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const CalcCGemmLong& obj){
		if(CALC_CGEMM_LONG_PRINT_FUNC>0) std::cout<<"nbytes(const CalcCGemmLong&):\n";
		int size=0;
		const int nt=obj.ntypes();
		size+=nbytes(static_cast<const Calculator&>(obj));
        size+=sizeof(obj.mix());//mix_
		size+=sizeof(int);//ntypes_
        size+=sizeof(double);//lambdaC_
        size+=sizeof(double);//lambdaS_
        size+=sizeof(double);//eps_
        size+=sizeof(double);//prec_
        size+=nt*sizeof(double);//radius_
		size+=nt*nt*sizeof(double);//aOver_
        size+=nt*nt*sizeof(double);//aRep_
        return size;
	}
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcCGemmLong& obj, char* arr){
		if(CALC_CGEMM_LONG_PRINT_FUNC>0) std::cout<<"pack(const CalcCGemmLong&,char*):\n";
		int pos=0;
		pos+=pack(static_cast<const Calculator&>(obj),arr+pos);
        std::memcpy(arr+pos,&obj.mix(),sizeof(obj.mix())); pos+=sizeof(obj.mix());//mix_
		std::memcpy(arr+pos,&obj.ntypes(),sizeof(int)); pos+=sizeof(int);//ntypes_
        std::memcpy(arr+pos,&obj.lambdaC(),sizeof(double)); pos+=sizeof(double);//lambdaC_
        std::memcpy(arr+pos,&obj.lambdaS(),sizeof(double)); pos+=sizeof(double);//lambdaS_
        std::memcpy(arr+pos,&obj.eps(),sizeof(double)); pos+=sizeof(double);//lambdaC_
        std::memcpy(arr+pos,&obj.prec(),sizeof(double)); pos+=sizeof(double);//lambdaS_
		const int nt=obj.ntypes();
		if(nt>0){
            std::memcpy(arr+pos,obj.radius().data(),nt*sizeof(double)); pos+=nt*sizeof(double);//radius_
			std::memcpy(arr+pos,obj.aOver().data(),nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//aOver_
            std::memcpy(arr+pos,obj.aRep().data(),nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//aRep_
		}
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcCGemmLong& obj, const char* arr){
		if(CALC_CGEMM_LONG_PRINT_FUNC>0) std::cout<<"unpack(CalcCGemmLong&,const char*):\n";
		int pos=0,nt=0;
		pos+=unpack(static_cast<Calculator&>(obj),arr+pos);
		if(obj.name()!=Calculator::Name::CGEMM_LONG) throw std::invalid_argument("serialize::unpack(CalcCGemmLong&,const char*): Invalid name.");
        std::memcpy(&obj.mix(),arr+pos,sizeof(obj.mix())); pos+=sizeof(obj.mix());//mix_
		std::memcpy(&nt,arr+pos,sizeof(int)); pos+=sizeof(int);//ntypes_
		std::memcpy(&obj.lambdaC(),arr+pos,sizeof(double)); pos+=sizeof(double);//lambdaC_
        std::memcpy(&obj.lambdaS(),arr+pos,sizeof(double)); pos+=sizeof(double);//lambdaS_
        std::memcpy(&obj.eps(),arr+pos,sizeof(double)); pos+=sizeof(double);//lambdaC_
        std::memcpy(&obj.prec(),arr+pos,sizeof(double)); pos+=sizeof(double);//lambdaS_
		obj.resize(nt);
		if(nt>0){
            std::memcpy(obj.radius().data(),arr+pos,nt*sizeof(double)); pos+=nt*sizeof(double);//radius_
			std::memcpy(obj.aOver().data(),arr+pos,nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//aOver_
            std::memcpy(obj.aRep().data(),arr+pos,nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//aRep_
		}
		obj.init();
		return pos;
	}
	
}