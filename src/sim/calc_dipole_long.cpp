// chem
#include "chem/units.hpp"
// math
#include "math/const.hpp"
// sim
#include "sim/calc_dipole_long.hpp"

//==== using statements ====

using math::constants::RadPI;
using math::constants::ZERO;

//==== operator ====

std::ostream& operator<<(std::ostream& out, const CalcDipoleLong& calc){
	return out<<static_cast<const Calculator&>(calc)<<" "<<calc.prec_;
}
	
//==== member functions ====

void CalcDipoleLong::read(Token& token){
	static_cast<Calculator&>(*this).read(token);
	//calc dipole_long 6.0 1e-12
	prec_=std::atof(token.next().c_str());
	if(prec_<=0) throw std::invalid_argument("CalcDipoleLong::read(Token&): invalid precision.");
	if(!token.end()){
		eps_=std::atof(token.next().c_str());
		if(eps_<=0.0) throw std::invalid_argument("CalcDipoleLong::read(Token&): Invalid epsilon.");
		kDipole_.eps()=eps_;
	}
}

void CalcDipoleLong::init(){
	kDipole_.prec()=prec_;
	kDipole_.rc()=rc_;
}

void CalcDipoleLong::init(const Structure& struc){
	kDipole_.init(struc);
}

double CalcDipoleLong::energy(Structure& struc, const NeighborList& nlist)const{
	const double ke=units::Consts::ke()*eps_;
	// k-space
	const double energyK=kDipole_.energy(struc);
	const double a=kDipole_.alpha();
	const double a2=a*a;
	// r-space
	double energyR=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const Eigen::Vector3d& mui=struc.dipole(i);
		for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
			const Eigen::Vector3d& muj=struc.dipole(jj);
			const Eigen::Vector3d drv=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
			const double& dr2=drv.squaredNorm();
			if(dr2<rc2_){
				const double dr=std::sqrt(dr2);
				const double ferfc=erfc(a*dr);
				const double fexp=std::exp(-a*a*dr2);
				const double b=(ferfc+2.0*a*dr/RadPI*fexp)/(dr2*dr);
				const double c=(3.0*ferfc+2.0*a*dr/RadPI*(3.0+2.0*a2*dr2)*fexp)/(dr2*dr2*dr);
				energyR+=ke*(mui.dot(muj)*b-mui.dot(drv)*muj.dot(drv)*c);
			}
		}
	}
	energyR*=0.5;
	if(CALC_DIPOLE_LONG_PRINT_FUNC>1){
		std::cout<<"energyR = "<<energyR<<"\n";
		std::cout<<"energyK = "<<energyK<<"\n";
	}
	//total
	struc.pe()+=energyR;
	return energyR+energyK;
}

double CalcDipoleLong::compute(Structure& struc, const NeighborList& nlist)const{
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
		const Eigen::Vector3d& mui=struc.dipole(i);
		for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
			const Eigen::Vector3d& muj=struc.dipole(jj);
			const Eigen::Vector3d drv=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
			const double& dr2=drv.squaredNorm();
			if(dr2<rc2_){
				const double dr=std::sqrt(dr2);
				const double dr4=dr2*dr2;
				const double dr6=dr4*dr2;
				const double ferfc=erfc(a*dr);
				const double fexp=std::exp(-a*a*dr2);
				const double b=(ferfc+2.0*a*dr/RadPI*fexp)/(dr2*dr);
				const double c=(3.0*ferfc+2.0*a*dr/RadPI*(3.0+2.0*a*a*dr2)*fexp)/(dr4*dr);
				const double d=(15.0*ferfc+2.0*a*dr/RadPI*(15.0+10.0*a2*dr2+4.0*a4*dr4)*exp(-a2*dr2))/(dr6*dr);
				energyR+=ke*(mui.dot(muj)*b-mui.dot(drv)*muj.dot(drv)*c);
				struc.force(i).noalias()+=ke*(
					(mui.dot(muj)*drv+mui*(muj.dot(drv))+muj*(mui.dot(drv)))*c
					-mui.dot(drv)*muj.dot(drv)*d*drv
				);
			}
		}
	}
	energyR*=0.5;
	if(CALC_DIPOLE_LONG_PRINT_FUNC>1){
		std::cout<<"energyR = "<<energyR<<"\n";
		std::cout<<"energyK = "<<energyK<<"\n";
	}
	//total
	struc.pe()+=energyR;
	return energyR+energyK;
}

//**********************************************
// serialization
//**********************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const CalcDipoleLong& obj){
		if(CALC_DIPOLE_LONG_PRINT_FUNC>0) std::cout<<"nbytes(const CalcDipoleLong&):\n";
		int size=0;
		size+=nbytes(static_cast<const Calculator&>(obj));
		size+=sizeof(double);//eps_
		size+=sizeof(double);//prec_
		return size;
	}
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcDipoleLong& obj, char* arr){
		if(CALC_DIPOLE_LONG_PRINT_FUNC>0) std::cout<<"pack(const CalcDipoleLong&,char*):\n";
		int pos=0;
		pos+=pack(static_cast<const Calculator&>(obj),arr+pos);
		std::memcpy(arr+pos,&obj.eps(),sizeof(double)); pos+=sizeof(double);//eps_
		std::memcpy(arr+pos,&obj.prec(),sizeof(double)); pos+=sizeof(double);//prec_
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcDipoleLong& obj, const char* arr){
		if(CALC_DIPOLE_LONG_PRINT_FUNC>0) std::cout<<"unpack(CalcDipoleLong&,const char*):\n";
		int pos=0;
		pos+=unpack(static_cast<Calculator&>(obj),arr+pos);
		if(obj.name()!=Calculator::Name::DIPOLE_LONG) throw std::invalid_argument("serialize::unpack(CalcDipoleLong&,const char*): Invalid name.");
		std::memcpy(&obj.eps(),arr+pos,sizeof(double)); pos+=sizeof(double);//eps_
		std::memcpy(&obj.prec(),arr+pos,sizeof(double)); pos+=sizeof(double);//prec_
		obj.init();
		return pos;
	}
	
}