// chem
#include "chem/units.hpp"
// math
#include "math/const.hpp"
// sim
#include "sim/calc_thole_long.hpp"

//==== using statements ====

using math::constants::RadPI;
using math::constants::ZERO;

//*****************************************
// IDD
//*****************************************

CalcTholeLong::IDD CalcTholeLong::IDD::read(const char* str){
	if(std::strcmp(str,"IDEAL")==0) return CalcTholeLong::IDD::IDEAL;
	else if(std::strcmp(str,"EXP")==0) return CalcTholeLong::IDD::EXP;
	else if(std::strcmp(str,"ERF")==0) return CalcTholeLong::IDD::ERF;
	else return CalcTholeLong::IDD::NONE;
}

const char* CalcTholeLong::IDD::name(const CalcTholeLong::IDD& name){
	switch(name){
		case CalcTholeLong::IDD::IDEAL: return "IDEAL";
		case CalcTholeLong::IDD::EXP: return "EXP";
		case CalcTholeLong::IDD::ERF: return "ERF";
		default: return "NONE";
	}
}

std::ostream& operator<<(std::ostream& out, const CalcTholeLong::IDD& name){
	switch(name){
		case CalcTholeLong::IDD::IDEAL: out<<"IDEAL"; break;
		case CalcTholeLong::IDD::EXP: out<<"EXP"; break;
		case CalcTholeLong::IDD::ERF: out<<"ERF"; break;
		default: out<<"NONE"; break;
	}
	return out;
}

//*****************************************
// CalcTholeLong
//*****************************************

//==== operator ====

std::ostream& operator<<(std::ostream& out, const CalcTholeLong& calc){
	return out<<static_cast<const Calculator&>(calc)<<" prec "<<calc.prec_<<" a "<<calc.a_<<" idd "<<calc.idd_;
}
	
//==== member functions ====

void CalcTholeLong::resize(int ntypes){
    if(CALC_THOLE_LONG_PRINT_FUNC>0) std::cout<<"CalcTholeLong::resize(int):\n";
	if(ntypes<0) throw std::invalid_argument("CalcTholeLong::resize(int): Invalid number of types.");
	ntypes_=ntypes;
	alpha_=Eigen::VectorXd::Zero(ntypes_);
}

void CalcTholeLong::read(Token& token){
	static_cast<Calculator&>(*this).read(token);
	//calc dipole_long 6.0 1e-12
	prec_=std::atof(token.next().c_str());
	if(prec_<=0) throw std::invalid_argument("CalcTholeLong::read(Token&): invalid precision.");
	if(!token.end()){
		eps_=std::atof(token.next().c_str());
		if(eps_<=0.0) throw std::invalid_argument("CalcTholeLong::read(Token&): Invalid epsilon.");
		kDipole_.eps()=eps_;
	}
}

void CalcTholeLong::coeff(Token& token){
    if(CALC_THOLE_LONG_PRINT_FUNC>0) std::cout<<"CalcTholeLong::coeff(Token&):\n";
    //coeff lj_cut type alpha
	const int type=std::atoi(token.next().c_str())-1;
	const double alpha=std::atof(token.next().c_str());
    //check type
	if(type>=ntypes_) throw std::invalid_argument("Invalid type.");
    if(alpha<=0) throw std::invalid_argument("Invalid alpha.");
    //assign
    alpha_[type]=alpha;
}

void CalcTholeLong::init(){
	kDipole_.prec()=prec_;
	kDipole_.rc()=rc_;
}

void CalcTholeLong::init(const Structure& struc){
	kDipole_.init(struc);
}

double CalcTholeLong::energy(Structure& struc, const NeighborList& nlist)const{
	const double ke=units::Consts::ke()*eps_;
	// k-space
	const double energyK=kDipole_.energy(struc);
	const double a=kDipole_.alpha();
	const double a2=a*a;
	// r-space
	double energyR=0;
	for(int i=0; i<struc.nAtoms(); ++i){
        const int ti=struc.type(i);
		const Eigen::Vector3d& mui=struc.dipole(i);
		for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
            const int tj=struc.type(jj);
			const Eigen::Vector3d& muj=struc.dipole(jj);
			const Eigen::Vector3d drv=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
			const double& dr2=drv.squaredNorm();
			if(dr2<rc2_){
				const double dr=std::sqrt(dr2);
				//subtract the ideal dipole energy
                energyR-=ke*(
                    mui.dot(muj)
                    -3.0*mui.dot(drv)*muj.dot(drv)/dr2
                )/(dr2*dr);
                //add the thole energy
                switch(idd_){
                    case IDD::IDEAL:{
                        energyR+=ke*(
							mui.dot(muj)
							-3.0*mui.dot(drv)*muj.dot(drv)/dr2
						)/(dr2*dr);
                    }; break;
                    case IDD::EXP:{
                        const double s=a_*dr;
						const double s2=s*s;
						energyR+=ke*std::pow(alpha_[ti]*alpha_[tj],-1.0/6.0)*(
							(1.0-(0.5*s2+s+1)*std::exp(-s))*mui.dot(muj)-
							(1.0-((1.0/6.0)*s2*s+0.5*s2+s+1.0)*std::exp(-s))*3.0*mui.dot(drv)*muj.dot(drv)/dr2
						)/(dr2*dr);
                    }; break;
                    case IDD::ERF:{
						const double a2=a_*a_;
                        const double fexp=2.0/RadPI*a_*std::exp(-dr2*a2);
                        energyR+=ke/(std::sqrt(std::pow(alpha_[ti]*2.0/3.0/RadPI,2.0/3.0)+std::pow(alpha_[tj]*2.0/3.0/RadPI,2.0/3.0)))*(
							2.0*a2*fexp*mui.dot(muj)/dr2-
							(3.0*mui.dot(drv)*muj.dot(drv)-mui.dot(muj)*dr2)*(std::erf(dr*a_)-dr*fexp)/(dr2*dr2*dr)
						);
                    }; break;
                    default: throw std::invalid_argument("CalcTholeLong::energy(const Structure&,const NeighborList&): Invalid dipole-dipole interaction.");
                }
                //add the ewald contribution
                const double ferfc=erfc(a*dr);
				const double fexp=std::exp(-a*a*dr2);
				const double b=(ferfc+2.0*a*dr/RadPI*fexp)/(dr2*dr);
				const double c=(3.0*ferfc+2.0*a*dr/RadPI*(3.0+2.0*a2*dr2)*fexp)/(dr2*dr2*dr);
				energyR+=ke*(mui.dot(muj)*b-mui.dot(drv)*muj.dot(drv)*c);
			}
		}
	}
	energyR*=0.5;
	if(CALC_THOLE_LONG_PRINT_FUNC>1){
		std::cout<<"energyR = "<<energyR<<"\n";
		std::cout<<"energyK = "<<energyK<<"\n";
	}
	//total
	struc.pe()+=energyR;
	return energyR+energyK;
}

double CalcTholeLong::compute(Structure& struc, const NeighborList& nlist)const{
	const double ke=units::Consts::ke()*eps_;
	// k-space
	//kDipole_.init(struc);
	const double energyK=kDipole_.compute(struc);
	const double a=kDipole_.alpha();
	const double a2=a*a;
	const double a4=a2*a2;
	// r-space
	double energyR=0;
	for(int i=0; i<struc.nAtoms(); ++i){
        const int ti=struc.type(i);
		const Eigen::Vector3d& mui=struc.dipole(i);
		for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
            const int tj=struc.type(jj);
			const Eigen::Vector3d& muj=struc.dipole(jj);
			const Eigen::Vector3d drv=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
			const double& dr2=drv.squaredNorm();
			if(dr2<rc2_){
				const double dr=std::sqrt(dr2);
				//subtract the ideal dipole energy
                energyR-=ke*(
                    mui.dot(muj)
                    -3.0*mui.dot(drv)*muj.dot(drv)/dr2
                )/(dr2*dr);
                //add the thole energy
                switch(idd_){
                    case IDD::IDEAL:{
                        energyR+=ke*(
							mui.dot(muj)
							-3.0*mui.dot(drv)*muj.dot(drv)/dr2
						)/(dr2*dr);
                    }; break;
                    case IDD::EXP:{
                        const double s=a_*dr;
						const double s2=s*s;
                        energyR+=ke*std::pow(alpha_[ti]*alpha_[tj],-1.0/6.0)*(
							(1.0-(0.5*s2+s+1)*std::exp(-s))*mui.dot(muj)-
							(1.0-((1.0/6.0)*s2*s+0.5*s2+s+1.0)*std::exp(-s))*3.0*mui.dot(drv)*muj.dot(drv)/dr2
						)/(dr2*dr);
                    }; break;
                    case IDD::ERF:{
						const double a2=a_*a_;
                        const double fexp=2.0/RadPI*a_*std::exp(-dr2*a2);
                        energyR+=ke/(std::sqrt(std::pow(alpha_[ti]*2.0/3.0/RadPI,2.0/3.0)+std::pow(alpha_[tj]*2.0/3.0/RadPI,2.0/3.0)))*(
							2.0*a2*fexp*mui.dot(muj)/dr2-
							(3.0*mui.dot(drv)*muj.dot(drv)-mui.dot(muj)*dr2)*(std::erf(dr*a_)-dr*fexp)/(dr2*dr2*dr)
						);
                    }; break;
                    default: throw std::invalid_argument("CalcTholeLong::energy(const Structure&,const NeighborList&): Invalid dipole-dipole interaction.");
                }
                //add the ewald contribution
                const double ferfc=erfc(a*dr);
				const double fexp=std::exp(-a*a*dr2);
				const double b=(ferfc+2.0*a*dr/RadPI*fexp)/(dr2*dr);
				const double c=(3.0*ferfc+2.0*a*dr/RadPI*(3.0+2.0*a2*dr2)*fexp)/(dr2*dr2*dr);
				energyR+=ke*(mui.dot(muj)*b-mui.dot(drv)*muj.dot(drv)*c);
				/*struc.force(i).noalias()+=ke*(
					(mui.dot(muj)*drv+mui*(muj.dot(drv))+muj*(mui.dot(drv)))*c
					-mui.dot(drv)*muj.dot(drv)*d*drv
				);*/
			}
		}
	}
	energyR*=0.5;
	if(CALC_THOLE_LONG_PRINT_FUNC>1){
		std::cout<<"energyR = "<<energyR<<"\n";
		std::cout<<"energyK = "<<energyK<<"\n";
	}
	//total
	struc.pe()+=energyR;
	return energyR+energyK;
}

Eigen::Matrix3d& CalcTholeLong::interMat(const Eigen::Vector3d& drv, Eigen::Matrix3d& mat, double a){
	const double ke=units::Consts::ke()*eps_;
	//compute distance
	const double dr=drv.norm();
	const double dr2=dr*dr;
	const double dr4=dr2*dr2;
	//zero the matrix
	mat.setZero();
	//subtract the ideal interaction
	mat.noalias()-=(3.0*drv*drv.transpose()-Eigen::Matrix3d::Identity()*dr2)/(dr4*dr);
	//add the thole interaction
	switch(idd_){
		case IDD::IDEAL:{
			mat.noalias()+=(3.0*drv*drv.transpose()-Eigen::Matrix3d::Identity()*dr2)/(dr4*dr);
		}; break;
		case IDD::EXP:{
			const double s=a*dr;
			const double s2=s*s;
			mat.noalias()+=3.0*drv*drv.transpose()*(1.0-((1.0/6.0)*s2*s+0.5*s2+s+1.0)*std::exp(-s))/(dr4*dr)
				-Eigen::Matrix3d::Identity()*(1.0-(0.5*s2+s+1)*std::exp(-s))/(dr2*dr);
		}; break;
		case IDD::ERF:{
			const double a2=a*a;
			const double fexp=2.0/RadPI*a*std::exp(-dr2*a2);
			mat.noalias()+=(3.0*drv*drv.transpose()-Eigen::Matrix3d::Identity()*dr2)*(std::erf(dr*a)-dr*fexp)/(dr4*dr);
				-2.0*a2/dr2*fexp*drv*drv.transpose();
		}; break;
		default: throw std::invalid_argument("Invalid dipole-dipole interaction.");
	}
	//add the ewald interaction - rspace
	const double A=kDipole_.alpha();
	const double A2=A*A;
	const double ferfc=erfc(A*dr);
	const double fexp=std::exp(-A2*dr2);
	const double b=(ferfc+2.0*A*dr/RadPI*fexp)/(dr2*dr);
	const double c=(3.0*ferfc+2.0*A*dr/RadPI*(3.0+2.0*A2*dr2)*fexp)/(dr2*dr2*dr);
	mat.noalias()+=(3.0*drv*drv.transpose()*c-Eigen::Matrix3d::Identity()*dr2*b)/(dr2*dr2*dr);
	//scale by Coulomb constant
	mat*=ke;
	//add the ewald interaction - kspace
	Eigen::Matrix3d kMat;
	kDipole_.interMat(drv,kMat);
	mat.noalias()+=kMat;
	//return matrix
	return mat;
}

//**********************************************
// serialization
//**********************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const CalcTholeLong& obj){
		if(CALC_THOLE_LONG_PRINT_FUNC>0) std::cout<<"nbytes(const CalcTholeLong&):\n";
		int size=0;
		size+=nbytes(static_cast<const Calculator&>(obj));
		size+=sizeof(double);//eps_
		size+=sizeof(double);//prec_
		return size;
	}
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcTholeLong& obj, char* arr){
		if(CALC_THOLE_LONG_PRINT_FUNC>0) std::cout<<"pack(const CalcTholeLong&,char*):\n";
		int pos=0;
		pos+=pack(static_cast<const Calculator&>(obj),arr+pos);
		std::memcpy(arr+pos,&obj.eps(),sizeof(double)); pos+=sizeof(double);//eps_
		std::memcpy(arr+pos,&obj.prec(),sizeof(double)); pos+=sizeof(double);//prec_
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcTholeLong& obj, const char* arr){
		if(CALC_THOLE_LONG_PRINT_FUNC>0) std::cout<<"unpack(CalcTholeLong&,const char*):\n";
		int pos=0;
		pos+=unpack(static_cast<Calculator&>(obj),arr+pos);
		if(obj.name()!=Calculator::Name::THOLE_LONG) throw std::invalid_argument("serialize::unpack(CalcTholeLong&,const char*): Invalid name.");
		std::memcpy(&obj.eps(),arr+pos,sizeof(double)); pos+=sizeof(double);//eps_
		std::memcpy(&obj.prec(),arr+pos,sizeof(double)); pos+=sizeof(double);//prec_
		obj.init();
		return pos;
	}
	
}