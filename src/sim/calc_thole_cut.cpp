// chem
#include "chem/units.hpp"
// string
#include "str/string.hpp"
// math
#include "math/const.hpp"
// sim
#include "sim/calc_thole_cut.hpp"

using math::constants::RadPI;

//*****************************************
// IDD
//*****************************************

CalcTholeCut::IDD CalcTholeCut::IDD::read(const char* str){
	if(std::strcmp(str,"IDEAL")==0) return CalcTholeCut::IDD::IDEAL;
	else if(std::strcmp(str,"EXP")==0) return CalcTholeCut::IDD::EXP;
	else if(std::strcmp(str,"ERF")==0) return CalcTholeCut::IDD::ERF;
	else return CalcTholeCut::IDD::NONE;
}

const char* CalcTholeCut::IDD::name(const CalcTholeCut::IDD& name){
	switch(name){
		case CalcTholeCut::IDD::IDEAL: return "IDEAL";
		case CalcTholeCut::IDD::EXP: return "EXP";
		case CalcTholeCut::IDD::ERF: return "ERF";
		default: return "NONE";
	}
}

std::ostream& operator<<(std::ostream& out, const CalcTholeCut::IDD& name){
	switch(name){
		case CalcTholeCut::IDD::IDEAL: out<<"IDEAL"; break;
		case CalcTholeCut::IDD::EXP: out<<"EXP"; break;
		case CalcTholeCut::IDD::ERF: out<<"ERF"; break;
		default: out<<"NONE"; break;
	}
	return out;
}

//*****************************************
// CalcTholeCut
//*****************************************

//==== operator ====

std::ostream& operator<<(std::ostream& out, const CalcTholeCut& calc){
	return out<<static_cast<const Calculator&>(calc)<<" a "<<calc.a_<<" idd "<<calc.idd_;
}

//==== member functions ====

void CalcTholeCut::resize(int ntypes){
    if(CALC_THOLE_CUT_PRINT_FUNC>0) std::cout<<"CalcTholeCut::resize(int):\n";
	if(ntypes<0) throw std::invalid_argument("CalcTholeCut::resize(int): Invalid number of types.");
	ntypes_=ntypes;
	alpha_=Eigen::VectorXd::Zero(ntypes_);
}

void CalcTholeCut::read(Token& token){
    static_cast<Calculator&>(*this).read(token);
    a_=std::atof(token.next().c_str());
    idd_=CalcTholeCut::IDD::read(string::to_upper(token.next()).c_str());
}

void CalcTholeCut::coeff(Token& token){
    if(CALC_THOLE_CUT_PRINT_FUNC>0) std::cout<<"CalcTholeCut::coeff(Token&):\n";
    //coeff lj_cut type alpha
	const int type=std::atoi(token.next().c_str())-1;
	const double alpha=std::atof(token.next().c_str());
    //check type
	if(type>=ntypes_) throw std::invalid_argument("Invalid type.");
    if(alpha<=0) throw std::invalid_argument("Invalid alpha.");
    //assign
    alpha_[type]=alpha;
}

double CalcTholeCut::energy(Structure& struc, const NeighborList& nlist)const{
    if(CALC_THOLE_CUT_PRINT_FUNC>0) std::cout<<"CalcTholeCut::energy(const Structure&,const NeighborList&):\n";
    const double ke=units::Consts::ke();
	double energy=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const int ti=struc.type(i);
		const Eigen::Vector3d& mui=struc.dipole(i);
		for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
			const int tj=struc.type(jj);
            const Eigen::Vector3d& muj=struc.dipole(jj);
			const Eigen::Vector3d drv=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
			const double dr2=drv.squaredNorm();
			if(dr2<rc2_){
				const double dr=std::sqrt(dr2);
                switch(idd_){
                    case IDD::IDEAL:{
                        energy+=ke*(
							mui.dot(muj)
							-3.0*mui.dot(drv)*muj.dot(drv)/dr2
						)/(dr2*dr);
                    }; break;
                    case IDD::EXP:{
                        const double b=a_*dr;
						const double b2=b*b;
						energy+=ke*std::pow(alpha_[ti]*alpha_[tj],-1.0/6.0)*(
							(1.0-(0.5*b2+b+1)*std::exp(-b))*mui.dot(muj)-
							(1.0-((1.0/6.0)*b2*b+0.5*b2+b+1.0)*std::exp(-b))*3.0*mui.dot(drv)*muj.dot(drv)/dr2
						)/(dr2*dr);
                    }; break;
                    case IDD::ERF:{
						const double a2=a_*a_;
                        const double fexp=2.0/RadPI*a_*std::exp(-dr2*a2);
                        energy+=ke/(std::sqrt(std::pow(alpha_[ti]*2.0/3.0/RadPI,2.0/3.0)+std::pow(alpha_[tj]*2.0/3.0/RadPI,2.0/3.0)))*(
							2.0*a2*fexp*mui.dot(muj)/dr2-
							(3.0*mui.dot(drv)*muj.dot(drv)-mui.dot(muj)*dr2)*(std::erf(dr*a_)-dr*fexp)/(dr2*dr2*dr)
						);
                    }; break;
                    default: throw std::invalid_argument(
						"CalcTholeCut::energy(const Structure&,const NeighborList&): Invalid dipole-dipole interaction."
					);
                }
			}
		}
	}
	energy*=0.5;
    struc.pe()+=energy;
	return energy;
}

double CalcTholeCut::compute(Structure& struc, const NeighborList& nlist)const{
    if(CALC_THOLE_CUT_PRINT_FUNC>0) std::cout<<"CalcTholeCut::compute(const Structure&,const NeighborList&):\n";
    const double ke=units::Consts::ke();
	double energy=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const int ti=struc.type(i);
		const Eigen::Vector3d& mui=struc.dipole(i);
		for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
			const int tj=struc.type(jj);
            const Eigen::Vector3d& muj=struc.dipole(jj);
			const Eigen::Vector3d drv=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
			const double dr2=drv.squaredNorm();
			if(dr2<rc2_){
				const double dr=std::sqrt(dr2);
				const double dr4=dr2*dr2;
				switch(idd_){
                    case IDD::IDEAL:{
                        energy+=ke*(
							mui.dot(muj)
							-3.0*mui.dot(drv)*muj.dot(drv)/dr2
						)/(dr2*dr);
                    }; break;
                    case IDD::EXP:{
                        const double b=a_*dr;
						const double b2=b*b;
						energy+=ke*std::pow(alpha_[ti]*alpha_[tj],-1.0/6.0)*(
							(1.0-(0.5*b2+b+1)*std::exp(-b))*mui.dot(muj)-
							(1.0-((1.0/6.0)*b2*b+0.5*b2+b+1.0)*std::exp(-b))*3.0*mui.dot(drv)*muj.dot(drv)/dr2
						)/(dr2*dr);
                    }; break;
                    case IDD::ERF:{
						const double a2=a_*a_;
                        const double fexp=2.0/RadPI*a_*std::exp(-dr2*a2);
						energy+=ke/(std::sqrt(std::pow(alpha_[ti]*2.0/3.0/RadPI,2.0/3.0)+std::pow(alpha_[tj]*2.0/3.0/RadPI,2.0/3.0)))*(
							2.0*a2*fexp*mui.dot(muj)/dr2-
							(3.0*mui.dot(drv)*muj.dot(drv)-mui.dot(muj)*dr2)*(std::erf(dr*a_)-dr*fexp)/(dr2*dr2*dr)
						);
                    }; break;
                    default: throw std::invalid_argument(
						"CalcTholeCut::energy(const Structure&,const NeighborList&): Invalid dipole-dipole interaction."
					);
                }
				/*struc.force(i).noalias()+=ke*(
					(mui.dot(muj)*drv+mui*(muj.dot(drv))+muj*(mui.dot(drv)))*3.0
					-mui.dot(drv)*muj.dot(drv)*15.0/dr2*drv
				)/(dr4*dr);*/
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
	
	template <> int nbytes(const CalcTholeCut& obj){
		if(CALC_THOLE_CUT_PRINT_FUNC>0) std::cout<<"nbytes(const CalcTholeCut&):\n";
		int size=0;
		const int nt=obj.ntypes();
		size+=nbytes(static_cast<const Calculator&>(obj));
		size+=sizeof(int);//ntypes_
		return size;
	}
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcTholeCut& obj, char* arr){
		if(CALC_THOLE_CUT_PRINT_FUNC>0) std::cout<<"pack(const CalcTholeCut&,char*):\n";
		int pos=0;
		pos+=pack(static_cast<const Calculator&>(obj),arr+pos);
		std::memcpy(arr+pos,&obj.ntypes(),sizeof(int)); pos+=sizeof(int);//ntypes_
		const int nt=obj.ntypes();
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcTholeCut& obj, const char* arr){
		if(CALC_THOLE_CUT_PRINT_FUNC>0) std::cout<<"unpack(CalcTholeCut&,const char*):\n";
		int pos=0,nt=0;
		pos+=unpack(static_cast<Calculator&>(obj),arr+pos);
		if(obj.name()!=Calculator::Name::THOLE_CUT) throw std::invalid_argument("serialize::unpack(CalcTholeCut&,const char*): Invalid name.");
		std::memcpy(&nt,arr+pos,sizeof(int)); pos+=sizeof(int);//ntypes_
		obj.resize(nt);
		obj.init();
		return pos;
	}
	
}
