// math
#include "chem/units.hpp"
// math
#include "math/const.hpp"
#include "math/special.hpp"
// sim
#include "sim/calc_bgem_cut.hpp"

using math::constants::PI;
using math::constants::RadPI;
using math::constants::ZERO;
using math::special::fmexp;

//==== operator ====

std::ostream& operator<<(std::ostream& out, const CalcBGemCut& calc){
	return out<<static_cast<const Calculator&>(calc);
}

//==== member functions ====

void CalcBGemCut::resize(int ntypes){
    if(CALC_BGEM_CUT_PRINT_FUNC>0) std::cout<<"CalcBGemCut::resize(int):\n";
	if(ntypes<0) throw std::invalid_argument("CalcBGemCut::resize(int): Invalid number of types.");
	ntypes_=ntypes;
	if(ntypes_>0){
        alpha_=Eigen::VectorXd::Zero(ntypes_);
        amp_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
		rep_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
        mu_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
		rmu_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
	}
}

void CalcBGemCut::init(){
    if(CALC_BGEM_CUT_PRINT_FUNC>0) std::cout<<"CalcBGemCut::init():\n";
    for(int i=0; i<ntypes_; ++i){
		for(int j=0; j<ntypes_; ++j){
            //amp_(i,j) = 0.5*(amp_(i,i)+amp_(j,j));
			//amp_(i,j) = sqrt(amp_(i,i)*amp_(j,j));
			if(amp_(i,i)<ZERO && amp_(j,j)<ZERO) amp_(i,j)=0.0;
			else amp_(i,j) = 2.0*(amp_(i,i)*amp_(j,j))/(amp_(i,i)+amp_(j,j));
			if(rep_(i,i)<ZERO && rep_(j,j)<ZERO) rep_(i,j)=0.0;
			else rep_(i,j) = 2.0*(rep_(i,i)*rep_(j,j))/(rep_(i,i)+rep_(j,j));
			mu_(i,j) = alpha_[i]*alpha_[j]/(alpha_[i]+alpha_[j]);
            rmu_(i,j)=sqrt(mu_(i,j));
		}
	}
}

void CalcBGemCut::read(Token& token){
    static_cast<Calculator&>(*this).read(token);
}

void CalcBGemCut::coeff(Token& token){
    if(CALC_BGEM_CUT_PRINT_FUNC>0) std::cout<<"CalcBGemCut::coeff(Token&):\n";
	//coeff lj_cut type radius amplitude
    //coeff cgem_cut type1 type2 alpha amp
	const int type=std::atoi(token.next().c_str())-1;
	const double alpha=std::atof(token.next().c_str());
    const double amp=std::atof(token.next().c_str());
	const double rep=std::atof(token.next().c_str());
	
    //check type
	if(type<0 || type>=ntypes_) throw std::invalid_argument("Invalid type.");
	
	alpha_[type]=alpha;
    amp_(type,type)=amp;
	rep_(type,type)=rep;
}

double CalcBGemCut::energy(Structure& struc, const NeighborList& nlist){
    if(CALC_BGEM_CUT_PRINT_FUNC>0) std::cout<<"CalcBGemCut::energy(const Structure&,const NeighborList&):\n";
    const double eps0=units::Consts::eps0();
	const double pe=1.0/(2.0*0.05*0.05);
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
            const Eigen::Vector3d disp=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
            const double dr2=disp.squaredNorm();
			if(dr2<rc2_){
				const double dr=sqrt(dr2);
                //Coulomb
				double eCoul=0.0;
                if(dr>ZERO){
					const double pf=1.0/(4.0*PI*eps0)*qi*qj/dr;
					eCoul=pf*std::erf(rmu_(ti,tj)*dr);
				} else eCoul=1.0/(4.0*PI*eps0)*qi*qj*2.0/RadPI*rmu_(ti,tj);
				// overlap
                const double eOver=amp_(ti,tj)*zi*zj*fmexp(-mu_(ti,tj)*dr2);
				// repulsion
				const double eRep=rep_(ti,tj)*fmexp(-pe*dr);
				//compute energy
				energy+=eCoul+eOver+eRep;
			}
		}
	}
	energy*=0.5;
    struc.pe()+=energy;
	return energy;
}

double CalcBGemCut::compute(Structure& struc, const NeighborList& nlist){
    if(CALC_BGEM_CUT_PRINT_FUNC>0) std::cout<<"CalcBGemCut::compute(const Structure&,const NeighborList&):\n";
    const double eps0=units::Consts::eps0();
	const double pe=1.0/(2.0*0.05*0.05);
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
            const Eigen::Vector3d disp=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
            const double dr2=disp.squaredNorm();
			if(dr2<rc2_){
				const double dr=sqrt(dr2);
                //Coulomb
				double eCoul=0.0;
                double fCoul=0.0;
				if(dr>ZERO){
					const double pf=1.0/(4.0*PI*eps0)*qi*qj/dr;
					const double ferf=std::erf(rmu_(ti,tj)*dr);
                	eCoul=pf*ferf;
					fCoul=pf/dr2*(
						ferf-2.0/RadPI*dr*rmu_(ti,tj)*fmexp(-mu_(ti,tj)*dr2)
					);
				} else eCoul=1.0/(4.0*PI*eps0)*qi*qj*2.0/RadPI*rmu_(ti,tj);
				// overlap
                const double eOver=amp_(ti,tj)*zi*zj*fmexp(-mu_(ti,tj)*dr2);
				const double fOver=2.0*mu_(ti,tj)*eOver;
				// repulsion
				const double eRep=rep_(ti,tj)*fmexp(-pe*dr);
				const double fRep=pe*eRep;
				//compute energy
				energy+=eCoul+eOver+eRep;
                //compute force
				struc.force(i).noalias()+=(fCoul+fOver+fRep)*disp;
			}
		}
	}
	energy*=0.5;
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
	
	template <> int nbytes(const CalcBGemCut& obj){
		if(CALC_BGEM_CUT_PRINT_FUNC>0) std::cout<<"nbytes(const CalcBGemCut&):\n";
		int size=0;
		const int nt=obj.ntypes();
		size+=nbytes(static_cast<const Calculator&>(obj));
		size+=sizeof(int);//ntypes_
        size+=nt*nt*sizeof(double);//amp_
		return size;
	}
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcBGemCut& obj, char* arr){
		if(CALC_BGEM_CUT_PRINT_FUNC>0) std::cout<<"pack(const CalcBGemCut&,char*):\n";
		int pos=0;
		pos+=pack(static_cast<const Calculator&>(obj),arr+pos);
		std::memcpy(arr+pos,&obj.ntypes(),sizeof(int)); pos+=sizeof(int);//ntypes_
        const int nt=obj.ntypes();
		if(nt>0){
            std::memcpy(arr+pos,obj.amp().data(),nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//amp_
		}
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcBGemCut& obj, const char* arr){
		if(CALC_BGEM_CUT_PRINT_FUNC>0) std::cout<<"unpack(CalcBGemCut&,const char*):\n";
		int pos=0,nt=0;
		pos+=unpack(static_cast<Calculator&>(obj),arr+pos);
		if(obj.name()!=Calculator::Name::CGEM_CUT) throw std::invalid_argument("serialize::unpack(CalcBGemCut&,const char*): Invalid name.");
		std::memcpy(&nt,arr+pos,sizeof(int)); pos+=sizeof(int);//ntypes_
		obj.resize(nt);
		if(nt>0){
            std::memcpy(obj.amp().data(),arr+pos,nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//amp_
		}
		obj.init();
		return pos;
	}
	
}