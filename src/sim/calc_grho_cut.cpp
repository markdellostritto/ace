// math
#include "chem/units.hpp"
// math
#include "math/const.hpp"
#include "math/special.hpp"
// string
#include "str/string.hpp"
// sim
#include "sim/calc_grho_cut.hpp"

using math::constants::RadPI;
using math::constants::ZERO;
using math::special::fmexp;

//==== contructors/destructors ====

CalcGRhoCut::CalcGRhoCut(double rc, Calculator::Mix mix):Calculator(Calculator::Name::GRHO_CUT,rc){
    mix_=mix;
	if(mix_==Calculator::Mix::NONE) throw std::invalid_argument("CalcGRhoCut::CalcGRhoCut(double,Calculator::Mix): Invalid mixing.\n");
}

//==== operator ====

std::ostream& operator<<(std::ostream& out, const CalcGRhoCut& calc){
	return out<<static_cast<const Calculator&>(calc)<<" mix "<<calc.mix_<<" eps "<<calc.eps_;
}

//==== member functions ====

void CalcGRhoCut::resize(int ntypes){
    if(CALC_GRHO_CUT_PRINT_FUNC>0) std::cout<<"CalcGRhoCut::resize(int):\n";
	Calculator::resize(ntypes);
	if(ntypes_>0){
        radius_=Eigen::VectorXd::Zero(ntypes_);
		gamma_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
		rgamma_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
	}
}

void CalcGRhoCut::init(){
    if(CALC_GRHO_CUT_PRINT_FUNC>0) std::cout<<"CalcGRhoCut::init():\n";
    Eigen::VectorXd alpha=Eigen::VectorXd::Zero(ntypes_);
    for(int i=0; i<ntypes_; ++i){
        alpha[i]=1.0/(2.0*radius_[i]*radius_[i]);
    }
    for(int i=0; i<ntypes_; ++i){
		for(int j=0; j<ntypes_; ++j){
            gamma_(i,j)=2.0*alpha[i]*alpha[j]/(alpha[i]+alpha[j]);
            rgamma_(i,j)=sqrt(0.5*gamma_(i,j));
		}
	}
}

void CalcGRhoCut::read(Token& token){
    static_cast<Calculator&>(*this).read(token);
    mix_=Calculator::Mix::read(string::to_upper(token.next()).c_str());
    if(mix_==Calculator::Mix::NONE) throw std::invalid_argument("CalcGRhoCut::read(Token&): Invalid mix type.");
}

void CalcGRhoCut::coeff(Token& token){
    if(CALC_GRHO_CUT_PRINT_FUNC>0) std::cout<<"CalcGRhoCut::coeff(Token&):\n";
	//coeff type radius 
	const int type=std::atoi(token.next().c_str())-1;
	const double radius=std::atof(token.next().c_str());
    //check type
	if(type>=ntypes_) throw std::invalid_argument("Invalid type.");
    if(radius<=0) throw std::invalid_argument("Invalid radius.");
    //assign
    radius_[type]=radius;
}

double CalcGRhoCut::energy(Structure& struc, const NeighborList& nlist)const{
    if(CALC_GRHO_CUT_PRINT_FUNC>0) std::cout<<"CalcGRhoCut::energy(const Structure&,const NeighborList&):\n";
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
				const double dr=sqrt(dr2);
                double eCoul=0.0;
                if(dr>ZERO){
					eCoul=ke*qi*qj/dr*std::erf(rgamma_(ti,tj)*dr);
				} else eCoul=ke*qi*qj*2.0/RadPI*rgamma_(ti,tj);
                //compute energy
				energy+=eCoul;
			}
		}
	}
	energy*=0.5;
    struc.pe()+=energy;
	return energy;
}

double CalcGRhoCut::energy(Structure& struc)const{
    if(CALC_GRHO_CUT_PRINT_FUNC>0) std::cout<<"CalcGRhoCut::energy(const Structure&):\n";
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
				const double dr=sqrt(dr2);
                double eCoul=0.0;
                if(dr>ZERO){
					eCoul=ke*qi*qj/dr*std::erf(rgamma_(ti,tj)*dr);
				} else eCoul=ke*qi*qj*2.0/RadPI*rgamma_(ti,tj);
                //compute energy
				energy+=eCoul;
			}
		}
	}
	struc.pe()+=energy;
	return energy;
}

double CalcGRhoCut::compute(Structure& struc, const NeighborList& nlist)const{
    if(CALC_GRHO_CUT_PRINT_FUNC>0) std::cout<<"CalcGRhoCut::compute(const Structure&,const NeighborList&):\n";
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
                double eCoul=0.0;
                double fCoul=0.0;
				if(dr>ZERO){
					const double pf=ke*qi*qj/dr;
                	const double ferf=std::erf(rgamma_(ti,tj)*dr);
					eCoul=pf*ferf;
					fCoul=pf/dr2*(
						ferf-dr*2.0/RadPI*rgamma_(ti,tj)*fmexp(-0.5*gamma_(ti,tj)*dr2)
					);
				} else eCoul=ke*qi*qj*2.0/RadPI*rgamma_(ti,tj);
                //compute energy
				energy+=eCoul;
                //compute force
				struc.force(i).noalias()+=fCoul*drv;
			}
		}
	}
	energy*=0.5;
    struc.pe()+=energy;
	return energy;
}

double CalcGRhoCut::compute(Structure& struc)const{
    if(CALC_GRHO_CUT_PRINT_FUNC>0) std::cout<<"CalcGRhoCut::compute(const Structure&):\n";
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
                double eCoul=0.0;
                double fCoul=0.0;
				if(dr>ZERO){
					const double pf=ke*qi*qj/dr;
                	const double ferf=std::erf(rgamma_(ti,tj)*dr);
					eCoul=pf*ferf;
					fCoul=pf/dr2*(
						ferf-dr*2.0/RadPI*rgamma_(ti,tj)*fmexp(-0.5*gamma_(ti,tj)*dr2)
					);
				} else eCoul=ke*qi*qj*2.0/RadPI*rgamma_(ti,tj);
                //compute energy
				energy+=eCoul;
                //compute force
				const Eigen::Vector3d fij=fCoul*drv;
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
	
	template <> int nbytes(const CalcGRhoCut& obj){
		if(CALC_GRHO_CUT_PRINT_FUNC>0) std::cout<<"nbytes(const CalcGRhoCut&):\n";
		int size=0;
		const int nt=obj.ntypes();
		size+=nbytes(static_cast<const Calculator&>(obj));
		size+=sizeof(obj.mix());//mix_
		size+=sizeof(int);//ntypes_
        size+=nt*sizeof(double);//radius_
		size+=sizeof(double);//eps_
		return size;
	}
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcGRhoCut& obj, char* arr){
		if(CALC_GRHO_CUT_PRINT_FUNC>0) std::cout<<"pack(const CalcGRhoCut&,char*):\n";
		int pos=0;
		pos+=pack(static_cast<const Calculator&>(obj),arr+pos);
		std::memcpy(arr+pos,&obj.mix(),sizeof(obj.mix())); pos+=sizeof(obj.mix());//mix_
		std::memcpy(arr+pos,&obj.ntypes(),sizeof(int)); pos+=sizeof(int);//ntypes_
        std::memcpy(arr+pos,&obj.eps(),sizeof(double)); pos+=sizeof(double);//eps_
		const int nt=obj.ntypes();
		if(nt>0){
            std::memcpy(arr+pos,obj.radius().data(),nt*sizeof(double)); pos+=nt*sizeof(double);//radius_
		}
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcGRhoCut& obj, const char* arr){
		if(CALC_GRHO_CUT_PRINT_FUNC>0) std::cout<<"unpack(CalcGRhoCut&,const char*):\n";
		int pos=0,nt=0;
		pos+=unpack(static_cast<Calculator&>(obj),arr+pos);
		if(obj.name()!=Calculator::Name::GRHO_CUT) throw std::invalid_argument("serialize::unpack(CalcGRhoCut&,const char*): Invalid name.");
		std::memcpy(&obj.mix(),arr+pos,sizeof(obj.mix())); pos+=sizeof(obj.mix());//mix_
		std::memcpy(&nt,arr+pos,sizeof(int)); pos+=sizeof(int);//ntypes_
		std::memcpy(&obj.eps(),arr+pos,sizeof(double)); pos+=sizeof(double);//eps_
        obj.resize(nt);
		if(nt>0){
            std::memcpy(obj.radius().data(),arr+pos,nt*sizeof(double)); pos+=nt*sizeof(double);//radius_
		}
		obj.init();
		return pos;
	}
	
}