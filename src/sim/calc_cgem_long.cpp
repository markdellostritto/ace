// math
#include "chem/units.hpp"
// math
#include "math/const.hpp"
#include "math/special.hpp"
// sim
#include "sim/calc_cgem_long.hpp"

using math::constants::PI;
using math::constants::RadPI;
using math::constants::ZERO;
using math::special::fmexp;

//==== contructors/destructors ====

CalcCGemLong::CalcCGemLong():Calculator(Calculator::Name::CGEM_CUT){
    lambdaC_=1.0;
    lambdaS_=1.0;
    eps_=1.0;
    prec_=1.0e-6;
}

CalcCGemLong::CalcCGemLong(double rc):Calculator(Calculator::Name::CGEM_CUT,rc){
    lambdaC_=1.0;
    lambdaS_=1.0;
    eps_=1.0;
    prec_=1.0e-6;
}

CalcCGemLong::CalcCGemLong(double rc, double lambdaC, double lambdaS):Calculator(Calculator::Name::CGEM_CUT,rc){
	lambdaC_=lambdaC;
    lambdaS_=lambdaS;
	if(lambdaC_<=0) throw std::invalid_argument("CalcCGemLong::CalcCGemLong(double,double,double): Invalid lambdaC\n");
    if(lambdaS_<=0) throw std::invalid_argument("CalcCGemLong::CalcCGemLong(double,double,double): Invalid lambdaS\n");
    eps_=1.0;
    prec_=1.0e-6;
}

//==== operator ====

std::ostream& operator<<(std::ostream& out, const CalcCGemLong& calc){
	return out<<static_cast<const Calculator&>(calc)
        <<" lambdaC "<<calc.lambdaC_<<" lambdaS "<<calc.lambdaS_<<" eps "<<calc.eps_<<" prec "<<calc.prec_;
}

//==== member functions ====

void CalcCGemLong::resize(int ntypes){
    if(CALC_CGEM_LONG_PRINT_FUNC>0) std::cout<<"CalcCGemLong::resize(int):\n";
	if(ntypes<0) throw std::invalid_argument("CalcCGemLong::resize(int): Invalid number of types.");
	ntypes_=ntypes;
	if(ntypes_>0){
        radius_=Eigen::VectorXd::Zero(ntypes_);
		aOver_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
		aRep_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
		muS_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
		muC_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
        rmuC_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
	}
}

void CalcCGemLong::init(){
    if(CALC_CGEM_LONG_PRINT_FUNC>0) std::cout<<"CalcCGemLong::init():\n";
    Eigen::VectorXd alphaS=Eigen::VectorXd::Zero(ntypes_);
    Eigen::VectorXd alphaC=Eigen::VectorXd::Zero(ntypes_);
    for(int i=0; i<ntypes_; ++i){
        alphaS[i]=lambdaS_/(2.0*radius_[i]*radius_[i]);
        alphaC[i]=lambdaC_/(2.0*radius_[i]*radius_[i]);
    }
    for(int i=0; i<ntypes_; ++i){
		for(int j=0; j<ntypes_; ++j){
            //aOver_(i,j) = 0.5*(aOver_(i,i)+aOver_(j,j));
			//aOver_(i,j) = sqrt(aOver_(i,i)*aOver_(j,j));
			if(aOver_(i,i)<ZERO && aOver_(j,j)<ZERO) aOver_(i,j)=0.0;
			else aOver_(i,j) = 2.0*(aOver_(i,i)*aOver_(j,j))/(aOver_(i,i)+aOver_(j,j));
			if(aRep_(i,i)<ZERO && aRep_(j,j)<ZERO) aRep_(i,j)=0.0;
			else aRep_(i,j) = 2.0*(aRep_(i,i)*aRep_(j,j))/(aRep_(i,i)+aRep_(j,j));
			muS_(i,j) = alphaS[i]*alphaS[j]/(alphaS[i]+alphaS[j]);
            muC_(i,j) = alphaC[i]*alphaC[j]/(alphaC[i]+alphaC[j]);
            rmuC_(i,j)=sqrt(muC_(i,j));
		}
	}
    coul_.prec()=prec_;
	coul_.rc()=rc_;
	coul_.econst()=false;
}

void CalcCGemLong::read(Token& token){
    static_cast<Calculator&>(*this).read(token);
    lambdaC_=std::atof(token.next().c_str());
    lambdaS_=std::atof(token.next().c_str());
    eps_=std::atof(token.next().c_str());
    prec_=std::atof(token.next().c_str());
    if(lambdaC_<=0) throw std::invalid_argument("CalcCGemLong::read(Token&): Invalid lambdaC\n");
    if(lambdaS_<=0) throw std::invalid_argument("CalcCGemLong::read(Token&): Invalid lambdaS\n");
    if(prec_<=0) throw std::invalid_argument("CalcCGemLong::read(Token&): Invalid prec\n");
    if(eps_<=0) throw std::invalid_argument("CalcCGemLong::read(Token&): Invalid eps\n");
}

void CalcCGemLong::coeff(Token& token){
    if(CALC_CGEM_LONG_PRINT_FUNC>0) std::cout<<"CalcCGemLong::coeff(Token&):\n";
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

double CalcCGemLong::energy(Structure& struc, const NeighborList& nlist){
    if(CALC_CGEM_LONG_PRINT_FUNC>0) std::cout<<"CalcCGemLong::energy(const Structure&,const NeighborList&):\n";
    // k-space
	coul_.init(struc);
	const double energyK=coul_.energy(struc);
	const double a=coul_.alpha();
    // r-space
    const double eps0=units::Consts::eps0();
	const double pe=1.0/(2.0*0.05*0.05);
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
			const Eigen::Vector3d disp=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
            const double dr2=disp.squaredNorm();
			if(dr2<rc2_){
				const double dr=sqrt(dr2);
                // Coulomb
				double eCoul=0.0;
                if(dr>ZERO){
					const double pf=1.0/(4.0*PI*eps0)*qi*qj/dr;
                	const double ferfg=std::erf(muC_(ti,tj)*dr);
                    const double ferfp=std::erf(a*dr);
                    eCoul=pf*(ferfg-ferfp);
				} else eCoul=1.0/(4.0*PI*eps0)*qi*qj*2.0/RadPI*(rmuC_(ti,tj)-a);
                // overlap
                const double eOver=aOver_(ti,tj)*zi*zj*fmexp(-muS_(ti,tj)*dr2);
				// repulsion
				const double eRep=aRep_(ti,tj)*fmexp(-pe*dr);
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

double CalcCGemLong::compute(Structure& struc, const NeighborList& nlist){
    if(CALC_CGEM_LONG_PRINT_FUNC>0) std::cout<<"CalcCGemLong::compute(const Structure&,const NeighborList&):\n";
    // k-space
    coul_.init(struc);
    const double energyK=coul_.compute(struc);
    const double a=coul_.alpha();
    // r-space
    const double eps0=units::Consts::eps0();
	const double pe=1.0/(2.0*0.05*0.05);
    double energyR=0;
	double energyCoul=0;
	double energyOver=0;
	double energyRep=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const int ti=struc.type(i);
        const double qi=struc.charge(i);
        const double zi=std::fabs(qi);
		for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
			const int tj=struc.type(jj);
            const double qj=struc.charge(jj);
            const double zj=std::fabs(qj);
			const Eigen::Vector3d disp=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
            const double dr2=disp.squaredNorm();
			if(dr2<rc2_){
				const double dr=sqrt(dr2);
                // Coulomb
				double eCoul=0.0;
                double fCoul=0.0;
				if(dr>ZERO){
                    const double pf=1.0/(4.0*PI*eps0)*qi*qj/dr;
                	const double ferfg=std::erf(muC_(ti,tj)*dr);
                    const double ferfp=std::erf(a*dr);
                    eCoul=pf*(ferfg-ferfp);
                    fCoul=pf/dr2*(
                        (ferfg-ferfp)
                        +2.0/RadPI*dr*(
                            -rmuC_(ti,tj)*std::exp(-muC_(ti,tj)*dr2)
                            +a*std::exp(-a*a*dr2)
                        )
                    );
				} else eCoul=1.0/(4.0*PI*eps0)*qi*qj*2.0/RadPI*(rmuC_(ti,tj)-a);
				// overlap
                const double eOver=aOver_(ti,tj)*zi*zj*fmexp(-muS_(ti,tj)*dr2);
				const double fOver=2.0*muS_(ti,tj)*eOver;
				// repulsion
				const double eRep=aRep_(ti,tj)*fmexp(-pe*dr);
				const double fRep=pe*eRep;
				//compute energy
				energyR+=eCoul+eOver+eRep;
				energyCoul+=eCoul;
				energyOver+=eOver;
				energyRep+=eRep;
                //compute force
				struc.force(i).noalias()+=(fCoul+fOver+fRep)*disp;
			}
		}
	}
    energyR*=0.5;
	energyCoul*=0.5;
	energyOver*=0.5;
	energyRep*=0.5;
    struc.pe()+=energyR;
	//return total energy
	//std::cout<<"eKspace "<<energyK<<" eCoul "<<energyCoul<<" eOver "<<energyOver<<" eRep "<<energyRep<<"\n";
	return energyR+energyK;
}

//**********************************************
// serialization
//**********************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const CalcCGemLong& obj){
		if(CALC_CGEM_LONG_PRINT_FUNC>0) std::cout<<"nbytes(const CalcCGemLong&):\n";
		int size=0;
		const int nt=obj.ntypes();
		size+=nbytes(static_cast<const Calculator&>(obj));
		size+=sizeof(int);//ntypes_
        size+=sizeof(double);//lambdaC_
        size+=sizeof(double);//lambdaS_
        size+=sizeof(double);//eps_
        size+=sizeof(double);//prec_
        size+=nt*sizeof(double);//radius_
		size+=nt*nt*sizeof(double);//aOver_
        size+=nt*nt*sizeof(double);//aRep_
        size+=nt*nt*sizeof(double);//muC_
        size+=nt*nt*sizeof(double);//muS_
		return size;
	}
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcCGemLong& obj, char* arr){
		if(CALC_CGEM_LONG_PRINT_FUNC>0) std::cout<<"pack(const CalcCGemLong&,char*):\n";
		int pos=0;
		pos+=pack(static_cast<const Calculator&>(obj),arr+pos);
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
            std::memcpy(arr+pos,obj.muC().data(),nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//muC_
            std::memcpy(arr+pos,obj.muS().data(),nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//muS_
		}
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcCGemLong& obj, const char* arr){
		if(CALC_CGEM_LONG_PRINT_FUNC>0) std::cout<<"unpack(CalcCGemLong&,const char*):\n";
		int pos=0,nt=0;
		pos+=unpack(static_cast<Calculator&>(obj),arr+pos);
		if(obj.name()!=Calculator::Name::CGEM_CUT) throw std::invalid_argument("serialize::unpack(CalcCGemLong&,const char*): Invalid name.");
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
            std::memcpy(obj.muC().data(),arr+pos,nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//muC_
            std::memcpy(obj.muS().data(),arr+pos,nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//muS_
		}
		obj.init();
		return pos;
	}
	
}