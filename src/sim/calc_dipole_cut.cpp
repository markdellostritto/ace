// chem
#include "chem/units.hpp"
// sim
#include "sim/calc_dipole_cut.hpp"

//==== operator ====

std::ostream& operator<<(std::ostream& out, const CalcDipoleCut& calc){
	return out<<static_cast<const Calculator&>(calc);
}

//==== member functions ====

double CalcDipoleCut::energy(Structure& struc, const NeighborList& nlist)const{
    if(CALC_DIPOLE_CUT_PRINT_FUNC>0) std::cout<<"CalcDipoleCut::energy(const Structure&,const NeighborList&):\n";
    const double ke=units::Consts::ke();
	double energy=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const Eigen::Vector3d& mui=struc.dipole(i);
		for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
            const Eigen::Vector3d& muj=struc.dipole(jj);
			const Eigen::Vector3d drv=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
			const double dr2=drv.squaredNorm();
			if(dr2<rc2_){
				const double dr=std::sqrt(dr2);
				energy+=ke*(mui.dot(muj)-mui.dot(drv)*muj.dot(drv)*3.0/dr2)/(dr2*dr);
			}
		}
	}
	energy*=0.5;
    struc.pe()+=energy;
	return energy;
}

double CalcDipoleCut::compute(Structure& struc, const NeighborList& nlist)const{
    if(CALC_DIPOLE_CUT_PRINT_FUNC>0) std::cout<<"CalcDipoleCut::compute(const Structure&,const NeighborList&):\n";
    const double ke=units::Consts::ke();
	double energy=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const Eigen::Vector3d& mui=struc.dipole(i);
		for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
            const Eigen::Vector3d& muj=struc.dipole(jj);
			const Eigen::Vector3d drv=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
			const double dr2=drv.squaredNorm();
			if(dr2<rc2_){
				const double dr=std::sqrt(dr2);
				const double dr4=dr2*dr2;
				energy+=ke*(mui.dot(muj)-mui.dot(drv)*muj.dot(drv)*3.0/dr2)/(dr2*dr);
				struc.force(i).noalias()+=ke*(
					(mui.dot(muj)*drv+mui*(muj.dot(drv))+muj*(mui.dot(drv)))*3.0
					-mui.dot(drv)*muj.dot(drv)*15.0/dr2*drv
				)/(dr4*dr);
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
	
	template <> int nbytes(const CalcDipoleCut& obj){
		if(CALC_DIPOLE_CUT_PRINT_FUNC>0) std::cout<<"nbytes(const CalcDipoleCut&):\n";
		int size=0;
		const int nt=obj.ntypes();
		size+=nbytes(static_cast<const Calculator&>(obj));
		size+=sizeof(int);//ntypes_
		return size;
	}
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcDipoleCut& obj, char* arr){
		if(CALC_DIPOLE_CUT_PRINT_FUNC>0) std::cout<<"pack(const CalcDipoleCut&,char*):\n";
		int pos=0;
		pos+=pack(static_cast<const Calculator&>(obj),arr+pos);
		std::memcpy(arr+pos,&obj.ntypes(),sizeof(int)); pos+=sizeof(int);//ntypes_
		const int nt=obj.ntypes();
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcDipoleCut& obj, const char* arr){
		if(CALC_DIPOLE_CUT_PRINT_FUNC>0) std::cout<<"unpack(CalcDipoleCut&,const char*):\n";
		int pos=0,nt=0;
		pos+=unpack(static_cast<Calculator&>(obj),arr+pos);
		if(obj.name()!=Calculator::Name::DIPOLE_CUT) throw std::invalid_argument("serialize::unpack(CalcDipoleCut&,const char*): Invalid name.");
		std::memcpy(&nt,arr+pos,sizeof(int)); pos+=sizeof(int);//ntypes_
		obj.resize(nt);
		obj.init();
		return pos;
	}
	
}