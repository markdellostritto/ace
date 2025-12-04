// math
#include "chem/units.hpp"
// math
#include "math/const.hpp"
#include "math/special.hpp"
// string
#include "str/string.hpp"
// sim
#include "sim/calc_cgemm_cut.hpp"

using math::constants::PI;
using math::constants::RadPI;
using math::constants::ZERO;
using math::special::fmexp;

//==== contructors/destructors ====

CalcCGemmCut::CalcCGemmCut(double rc, double lambdaC, double lambdaS):Calculator(Calculator::Name::CGEMM_CUT,rc){
    lambdaC_=lambdaC;
    lambdaS_=lambdaS;
	if(lambdaC_<=0) throw std::invalid_argument("CalcCGemmCut::CalcCGemmCut(double,double,double): Invalid lambdaC\n");
    if(lambdaS_<=0) throw std::invalid_argument("CalcCGemmCut::CalcCGemmCut(double,double,double): Invalid lambdaS\n");
}

//==== operator ====

std::ostream& operator<<(std::ostream& out, const CalcCGemmCut& calc){
	return out<<static_cast<const Calculator&>(calc)
        <<" mean "<<calc.mix_<<" lambdaC "<<calc.lambdaC_<<" lambdaS "<<calc.lambdaS_<<" rRep "<<calc.rRep_;
}

//==== member functions ====

void CalcCGemmCut::resize(int ntypes){
    if(CALC_CGEMM_CUT_PRINT_FUNC>0) std::cout<<"CalcCGemmCut::resize(int):\n";
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

void CalcCGemmCut::init(){
    if(CALC_CGEMM_CUT_PRINT_FUNC>0) std::cout<<"CalcCGemmCut::init():\n";
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
                    aOver_(i,j) = 0.5*(aOver_(i,i)+aOver_(j,j));
                break;
                case Calculator::Mix::GEOMETRIC:
                    aOver_(i,j) = sqrt(aOver_(i,i)*aOver_(j,j));
                    aRep_(i,j) = sqrt(aRep_(i,i)*aRep_(j,j));
                break;
                case Calculator::Mix::HARMONIC:
                    if(aOver_(i,i)<ZERO && aOver_(j,j)<ZERO) aOver_(i,j)=0.0;
                    else aOver_(i,j) = 2.0*(aOver_(i,i)*aOver_(j,j))/(aOver_(i,i)+aOver_(j,j));
                    if(aRep_(i,i)<ZERO && aRep_(j,j)<ZERO) aRep_(i,j)=0.0;
                    else aRep_(i,j) = 2.0*(aRep_(i,i)*aRep_(j,j))/(aRep_(i,i)+aRep_(j,j));
                break;
                default:
                    throw std::invalid_argument("CalcCGemmCut::init(): Invalid mean.");
                break;
            }
			gammaS_(i,j) = 2.0*alphaS[i]*alphaS[j]/(alphaS[i]+alphaS[j]);
            gammaC_(i,j) = 2.0*alphaC[i]*alphaC[j]/(alphaC[i]+alphaC[j]);
            rgammaC_(i,j)=sqrt(0.5*gammaC_(i,j));
		}
	}
}

void CalcCGemmCut::read(Token& token){
    static_cast<Calculator&>(*this).read(token);
    lambdaC_=std::atof(token.next().c_str());
    lambdaS_=std::atof(token.next().c_str());
	if(!token.end()){
        rRep_=std::atof(token.next().c_str());
    }
    if(!token.end()){
        mix_=Calculator::Mix::read(string::to_upper(token.next()).c_str());
    }
    if(lambdaC_<=0) throw std::invalid_argument("CalcCGemmCut::read(Token&): Invalid lambdaC.");
    if(lambdaS_<=0) throw std::invalid_argument("CalcCGemmCut::read(Token&): Invalid lambdaS.");
	if(rRep_<=0) throw std::invalid_argument("CalcCGemmCut::read(Token&): Invalid rRep.");
    if(mix_==Calculator::Mix::NONE) throw std::invalid_argument("CalcCGemmCut::read(Token&): Invalid mix type.");
}

void CalcCGemmCut::coeff(Token& token){
    if(CALC_CGEMM_CUT_PRINT_FUNC>0) std::cout<<"CalcCGemmCut::coeff(Token&):\n";
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

double CalcCGemmCut::energy(Structure& struc, const NeighborList& nlist)const{
    if(CALC_CGEMM_CUT_PRINT_FUNC>0) std::cout<<"CalcCGemmCut::energy(const Structure&,const NeighborList&):\n";
    const double ke=units::Consts::ke();
	const double cRep_=1.0/rRep_;
	double energy=0;
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
					eCoul=ke*qi*qj/dr*std::erf(rgammaC_(ti,tj)*dr);
				} else eCoul=ke*qi*qj*2.0/RadPI*rgammaC_(ti,tj);
                // overlap
                const double eOver=aOver_(ti,tj)*zi*zj*fmexp(-0.5*gammaS_(ti,tj)*dr2);
				// repulsion
				const double eRep=aRep_(ti,tj)*fmexp(-cRep_*dr);
				//const double eRep=aRep_(ti,tj)*cRep_*fmexp(-cRep_*dr);
				//compute energy
				energy+=eCoul+eOver+eRep;
			}
		}
	}
	energy*=0.5;
    struc.pe()+=energy;
	return energy;
}

double CalcCGemmCut::energy(Structure& struc)const{
    if(CALC_CGEMM_CUT_PRINT_FUNC>0) std::cout<<"CalcCGemmCut::energy(const Structure&):\n";
    const double ke=units::Consts::ke();
	const double cRep_=1.0/rRep_;
	double energy=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const int ti=struc.type(i);
        const double qi=struc.charge(i);
        const double zi=std::fabs(qi);
        for(int j=i+1; j<struc.nAtoms(); ++j){
			const int tj=struc.type(j);
            const double qj=struc.charge(j);
            const double zj=std::fabs(qj);
			const double dr2=(struc.posn(i)-struc.posn(j)).squaredNorm();
			if(dr2<rc2_){
				const double dr=sqrt(dr2);
                // Coulomb
				double eCoul=0.0;
                if(dr>ZERO){
					eCoul=ke*qi*qj/dr*std::erf(rgammaC_(ti,tj)*dr);
				} else eCoul=ke*qi*qj*2.0/RadPI*rgammaC_(ti,tj);
                // overlap
                const double eOver=aOver_(ti,tj)*zi*zj*fmexp(-0.5*gammaS_(ti,tj)*dr2);
				// repulsion
				const double eRep=aRep_(ti,tj)*fmexp(-cRep_*dr);
				//const double eRep=aRep_(ti,tj)*cRep_*fmexp(-cRep_*dr);
				//compute energy
				energy+=eCoul+eOver+eRep;
			}
		}
	}
	struc.pe()+=energy;
	return energy;
}

double CalcCGemmCut::compute(Structure& struc, const NeighborList& nlist)const{
    if(CALC_CGEMM_CUT_PRINT_FUNC>0) std::cout<<"CalcCGemmCut::compute(const Structure&,const NeighborList&):\n";
    const double ke=units::Consts::ke();
	const double cRep_=1.0/rRep_;
	double energy=0;
	double energyCoul=0;
	double energyOver=0;
	double energyRep=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		//store atom data
		const int ti=struc.type(i);
        const double qi=struc.charge(i);
        const double zi=std::fabs(qi);
        for(int j=0; j<nlist.size(i); ++j){
			//store atom data
			const int jj=nlist.index(i,j);
			const int tj=struc.type(jj);
            const double qj=struc.charge(jj);
            const double zj=std::fabs(qj);
			//compute distance
			const Eigen::Vector3d drv=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
            const double dr2=drv.squaredNorm();
			//compute interaction
			if(dr2<rc2_){
				const double dr=sqrt(dr2);
                // Coulomb
				double eCoul=0.0;
                double fCoul=0.0;
				if(dr>ZERO){
					const double pf=ke*qi*qj/dr;
                	const double ferf=std::erf(rgammaC_(ti,tj)*dr);
					eCoul=pf*ferf;
					fCoul=pf/dr2*(
						ferf-dr*2.0/RadPI*rgammaC_(ti,tj)*fmexp(-0.5*gammaC_(ti,tj)*dr2)
					);
				} else eCoul=ke*qi*qj*2.0/RadPI*rgammaC_(ti,tj);
                // overlap
                const double eOver=aOver_(ti,tj)*zi*zj*fmexp(-0.5*gammaS_(ti,tj)*dr2);
				const double fOver=gammaS_(ti,tj)*eOver;
				// repulsion
				const double eRep=aRep_(ti,tj)*fmexp(-cRep_*dr);
				//const double eRep=aRep_(ti,tj)*cRep_*fmexp(-cRep_*dr);
				const double fRep=cRep_*eRep;
				//compute energy
				energy+=eCoul+eOver+eRep;
                energyCoul+=eCoul;
				energyOver+=eOver;
				energyRep+=eRep;
                //compute force
				struc.force(i).noalias()+=(fCoul+fOver+fRep)*drv;
			}
		}
	}
	energy*=0.5;
    struc.pe()+=energy;
	//std::cout<<" eCoul "<<energyCoul<<" eOver "<<energyOver<<" eRep "<<energyRep<<"\n";
	return energy;
}

double CalcCGemmCut::compute(Structure& struc)const{
    if(CALC_CGEMM_CUT_PRINT_FUNC>0) std::cout<<"CalcCGemmCut::compute(const Structure&):\n";
    const double ke=units::Consts::ke();
	const double cRep_=1.0/rRep_;
	double energy=0;
	double energyCoul=0;
	double energyOver=0;
	double energyRep=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		//store atom data
		const int ti=struc.type(i);
        const double qi=struc.charge(i);
        const double zi=std::fabs(qi);
        for(int j=i+1; j<struc.nAtoms(); ++j){
			//store atom data
			const int tj=struc.type(j);
            const double qj=struc.charge(j);
            const double zj=std::fabs(qj);
			//compute distance
			const Eigen::Vector3d drv=struc.posn(i)-struc.posn(j);
            const double dr2=drv.squaredNorm();
			//compute interaction
			if(dr2<rc2_){
				const double dr=sqrt(dr2);
                // Coulomb
				double eCoul=0.0;
                double fCoul=0.0;
				if(dr>ZERO){
					const double pf=ke*qi*qj/dr;
                	const double ferf=std::erf(rgammaC_(ti,tj)*dr);
					eCoul=pf*ferf;
					fCoul=pf/dr2*(
						ferf-dr*2.0/RadPI*rgammaC_(ti,tj)*fmexp(-0.5*gammaC_(ti,tj)*dr2)
					);
				} else eCoul=ke*qi*qj*2.0/RadPI*rgammaC_(ti,tj);
                // overlap
                const double eOver=aOver_(ti,tj)*zi*zj*fmexp(-0.5*gammaS_(ti,tj)*dr2);
				const double fOver=gammaS_(ti,tj)*eOver;
				// repulsion
				const double eRep=aRep_(ti,tj)*fmexp(-cRep_*dr);
				//const double eRep=aRep_(ti,tj)*cRep_*fmexp(-cRep_*dr);
				const double fRep=cRep_*eRep;
				//compute energy
				energy+=eCoul+eOver+eRep;
                energyCoul+=eCoul;
				energyOver+=eOver;
				energyRep+=eRep;
                //compute force
				const Eigen::Vector3d fij=(fCoul+fOver+fRep)*drv;
				struc.force(i).noalias()+=fij;
				struc.force(j).noalias()-=fij;
			}
		}
	}
	struc.pe()+=energy;
	//std::cout<<" eCoul "<<energyCoul<<" eOver "<<energyOver<<" eRep "<<energyRep<<"\n";
	return energy;
}

//**********************************************
// serialization
//**********************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const CalcCGemmCut& obj){
		if(CALC_CGEMM_CUT_PRINT_FUNC>0) std::cout<<"nbytes(const CalcCGemmCut&):\n";
		int size=0;
		const int nt=obj.ntypes();
		size+=nbytes(static_cast<const Calculator&>(obj));
		size+=sizeof(obj.mix());//mix_
		size+=sizeof(int);//ntypes_
        size+=sizeof(double);//lambdaC_
        size+=sizeof(double);//lambdaS_
		size+=sizeof(double);//rRep_
        size+=nt*sizeof(double);//radius_
		size+=nt*nt*sizeof(double);//aOver_
        size+=nt*nt*sizeof(double);//aRep_
		return size;
	}
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcCGemmCut& obj, char* arr){
		if(CALC_CGEMM_CUT_PRINT_FUNC>0) std::cout<<"pack(const CalcCGemmCut&,char*):\n";
		int pos=0;
		pos+=pack(static_cast<const Calculator&>(obj),arr+pos);
		std::memcpy(arr+pos,&obj.mix(),sizeof(obj.mix())); pos+=sizeof(obj.mix());//mix_
		std::memcpy(arr+pos,&obj.ntypes(),sizeof(int)); pos+=sizeof(int);//ntypes_
        std::memcpy(arr+pos,&obj.lambdaC(),sizeof(double)); pos+=sizeof(double);//lambdaC_
        std::memcpy(arr+pos,&obj.lambdaS(),sizeof(double)); pos+=sizeof(double);//lambdaS_
		std::memcpy(arr+pos,&obj.rRep(),sizeof(double)); pos+=sizeof(double);//rRep_
		const int nt=obj.ntypes();
		if(nt>0){
            std::memcpy(arr+pos,obj.radius().data(),nt*sizeof(double)); pos+=nt*sizeof(double);//radius_
			std::memcpy(arr+pos,obj.aOver().data(),nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//aOver_
            std::memcpy(arr+pos,obj.aRep().data(),nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//aOver_
		}
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcCGemmCut& obj, const char* arr){
		if(CALC_CGEMM_CUT_PRINT_FUNC>0) std::cout<<"unpack(CalcCGemmCut&,const char*):\n";
		int pos=0,nt=0;
		pos+=unpack(static_cast<Calculator&>(obj),arr+pos);
		if(obj.name()!=Calculator::Name::CGEMM_CUT) throw std::invalid_argument("serialize::unpack(CalcCGemmCut&,const char*): Invalid name.");
		std::memcpy(&obj.mix(),arr+pos,sizeof(obj.mix())); pos+=sizeof(obj.mix());//mix_
		std::memcpy(&nt,arr+pos,sizeof(int)); pos+=sizeof(int);//ntypes_
		std::memcpy(&obj.lambdaC(),arr+pos,sizeof(double)); pos+=sizeof(double);//lambdaC_
        std::memcpy(&obj.lambdaS(),arr+pos,sizeof(double)); pos+=sizeof(double);//lambdaS_
		std::memcpy(&obj.rRep(),arr+pos,sizeof(double)); pos+=sizeof(double);//rRep_
		obj.resize(nt);
		if(nt>0){
            std::memcpy(obj.radius().data(),arr+pos,nt*sizeof(double)); pos+=nt*sizeof(double);//radius_
			std::memcpy(obj.aOver().data(),arr+pos,nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//aOver_
            std::memcpy(obj.aRep().data(),arr+pos,nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//aOver_
		}
		obj.init();
		return pos;
	}
	
}