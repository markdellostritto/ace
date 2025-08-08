// sim
#include "sim/calc_coul_cut.hpp"

//==== operator ====

std::ostream& operator<<(std::ostream& out, const CalcCoulCut& calc){
	return out<<static_cast<const Calculator&>(calc);
}

//==== member functions ====

void CalcCoulCut::resize(int ntypes){
    if(CALC_COUL_CUT_PRINT_FUNC>0) std::cout<<"CalcCoulCut::resize(int):\n";
	if(ntypes<0) throw std::invalid_argument("CalcCoulCut::resize(int): Invalid number of types.");
	ntypes_=ntypes;
}

void CalcCoulCut::init(){
    if(CALC_COUL_CUT_PRINT_FUNC>0) std::cout<<"CalcCoulCut::init():\n";
}

void CalcCoulCut::read(Token& token){
    static_cast<Calculator&>(*this).read(token);
}

void CalcCoulCut::coeff(Token& token){
    if(CALC_COUL_CUT_PRINT_FUNC>0) std::cout<<"CalcCoulCut::coeff(Token&):\n";
	//coeff coul_cut type1 type2 
}

double CalcCoulCut::energy(Structure& struc, const NeighborList& nlist){
    if(CALC_COUL_CUT_PRINT_FUNC>0) std::cout<<"CalcCoulCut::energy(const Structure&,const NeighborList&):\n";
	double energy=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const double qi=struc.charge(i);
		for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
            const double qj=struc.charge(jj);
			const Eigen::Vector3d disp=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
			const double dr2=disp.squaredNorm();
			if(dr2<rc2_){
				const double dri=1.0/sqrt(dr2);
				energy+=qi*qj*dri;
			}
		}
	}
	energy*=0.5;
    struc.pe()+=energy;
	return energy;
}

double CalcCoulCut::compute(Structure& struc, const NeighborList& nlist){
    if(CALC_COUL_CUT_PRINT_FUNC>0) std::cout<<"CalcCoulCut::compute(const Structure&,const NeighborList&):\n";
	double energy=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const double qi=struc.charge(i);
		for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
            const double qj=struc.charge(jj);
			const Eigen::Vector3d disp=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
			const double dr2=disp.squaredNorm();
			if(dr2<rc2_){
				const double dri=1.0/sqrt(dr2);
				energy+=qi*qj*dri;
				struc.force(i).noalias()+=qi*qj*dri/dr2*disp;
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
	
	template <> int nbytes(const CalcCoulCut& obj){
		if(CALC_COUL_CUT_PRINT_FUNC>0) std::cout<<"nbytes(const CalcCoulCut&):\n";
		int size=0;
		const int nt=obj.ntypes();
		size+=nbytes(static_cast<const Calculator&>(obj));
		size+=sizeof(int);//ntypes_
		return size;
	}
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcCoulCut& obj, char* arr){
		if(CALC_COUL_CUT_PRINT_FUNC>0) std::cout<<"pack(const CalcCoulCut&,char*):\n";
		int pos=0;
		pos+=pack(static_cast<const Calculator&>(obj),arr+pos);
		std::memcpy(arr+pos,&obj.ntypes(),sizeof(int)); pos+=sizeof(int);//ntypes_
		const int nt=obj.ntypes();
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcCoulCut& obj, const char* arr){
		if(CALC_COUL_CUT_PRINT_FUNC>0) std::cout<<"unpack(CalcCoulCut&,const char*):\n";
		int pos=0,nt=0;
		pos+=unpack(static_cast<Calculator&>(obj),arr+pos);
		if(obj.name()!=Calculator::Name::LJ_CUT) throw std::invalid_argument("serialize::unpack(CalcCoulCut&,const char*): Invalid name.");
		std::memcpy(&nt,arr+pos,sizeof(int)); pos+=sizeof(int);//ntypes_
		obj.resize(nt);
		obj.init();
		return pos;
	}
	
}