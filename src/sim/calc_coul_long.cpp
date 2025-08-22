// chem
#include "chem/units.hpp"
// math
#include "math/const.hpp"
// sim
#include "sim/calc_coul_long.hpp"

//==== using statements ====

using math::constants::RadPI;
using math::constants::ZERO;

//==== constructors/destructors ====

CalcCoulLong::CalcCoulLong():Calculator(Calculator::Name::COUL_LONG){
    eps_=1.0;
    prec_=1.0e-6;
}

CalcCoulLong::CalcCoulLong(double rc):Calculator(Calculator::Name::COUL_LONG,rc){
    eps_=1.0;
    prec_=1.0e-6;
}

//==== operator ====

std::ostream& operator<<(std::ostream& out, const CalcCoulLong& calc){
	return out<<static_cast<const Calculator&>(calc)<<" "<<calc.prec_;
}
	
//==== member functions ====

void CalcCoulLong::read(Token& token){
	static_cast<Calculator&>(*this).read(token);
	//calc coul_long 6.0 1e-12
	prec_=std::atof(token.next().c_str());
	if(prec_<=0) throw std::invalid_argument("CalcCoulLong::read(Token&): invalid precision.");
	if(!token.end()){
		eps_=std::atof(token.next().c_str());
		if(eps_<=0.0) throw std::invalid_argument("CalcCoulLong::read(Token&): Invalid epsilon.");
		coul_.eps()=eps_;
	}
}

void CalcCoulLong::init(){
	coul_.prec()=prec_;
	coul_.rc()=rc_;
}

void CalcCoulLong::init(const Structure& struc){
	coul_.init(struc);
}

double CalcCoulLong::energy(Structure& struc, const NeighborList& nlist)const{
	const double ke=units::Consts::ke()*eps_;
	// k-space
	//coul_.init(struc);
	const double energyK=coul_.energy(struc);
	const double a=coul_.alpha();
	// r-space
	double energyR=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const double& qi=struc.charge(i);
		for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
			const double& qj=struc.charge(jj);
			const Eigen::Vector3d disp=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
			const double& dr2=disp.squaredNorm();
			if(dr2<rc2_){
				const double dr=std::sqrt(dr2);
				const double pf=ke*qi*qj/dr;
				const double ferfc=erfc(a*dr);
				energyR+=pf*ferfc;
			}
		}
	}
	energyR*=0.5;
	if(CALC_COUL_LONG_PRINT_FUNC>1){
		std::cout<<"energyR = "<<energyR<<"\n";
		std::cout<<"energyK = "<<energyK<<"\n";
	}
	//total
	struc.pe()+=energyR;
	return energyR+energyK;
}

double CalcCoulLong::compute(Structure& struc, const NeighborList& nlist)const{
	const double ke=units::Consts::ke()*eps_;
	// k-space
	//coul_.init(struc);
	const double energyK=coul_.compute(struc);
	const double a=coul_.alpha();
	// r-space
	double energyR=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const double& qi=struc.charge(i);
		for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
			const double& qj=struc.charge(jj);
			const Eigen::Vector3d disp=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
			const double& dr2=disp.squaredNorm();
			if(dr2<rc2_){
				const double dr=std::sqrt(dr2);
				const double pf=ke*qi*qj/dr;
				const double ferfc=erfc(a*dr);
				const double expf=exp(-a*a*dr2);
				energyR+=pf*ferfc;
				struc.force(i).noalias()+=pf*(ferfc+2.0*a*dr/RadPI*expf)/dr2*disp;
			}
		}
	}
	energyR*=0.5;
	if(CALC_COUL_LONG_PRINT_FUNC>1){
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
	
	template <> int nbytes(const CalcCoulLong& obj){
		if(CALC_COUL_LONG_PRINT_FUNC>0) std::cout<<"nbytes(const CalcCoulLong&):\n";
		int size=0;
		size+=nbytes(static_cast<const Calculator&>(obj));
		size+=sizeof(double);//eps_
		size+=sizeof(double);//prec_
		return size;
	}
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcCoulLong& obj, char* arr){
		if(CALC_COUL_LONG_PRINT_FUNC>0) std::cout<<"pack(const CalcCoulLong&,char*):\n";
		int pos=0;
		pos+=pack(static_cast<const Calculator&>(obj),arr+pos);
		std::memcpy(arr+pos,&obj.eps(),sizeof(double)); pos+=sizeof(double);//eps_
		std::memcpy(arr+pos,&obj.prec(),sizeof(double)); pos+=sizeof(double);//prec_
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcCoulLong& obj, const char* arr){
		if(CALC_COUL_LONG_PRINT_FUNC>0) std::cout<<"unpack(CalcCoulLong&,const char*):\n";
		int pos=0;
		pos+=unpack(static_cast<Calculator&>(obj),arr+pos);
		if(obj.name()!=Calculator::Name::COUL_LONG) throw std::invalid_argument("serialize::unpack(CalcCoulLong&,const char*): Invalid name.");
		std::memcpy(&obj.eps(),arr+pos,sizeof(double)); pos+=sizeof(double);//eps_
		std::memcpy(&obj.prec(),arr+pos,sizeof(double)); pos+=sizeof(double);//prec_
		obj.init();
		return pos;
	}
	
}