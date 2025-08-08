// sim
#include "sim/calc_lj_cut.hpp"

//==== operator ====

std::ostream& operator<<(std::ostream& out, const CalcLJCut& calc){
	return out<<static_cast<const Calculator&>(calc);
}

//==== member functions ====

void CalcLJCut::resize(int ntypes){
    if(CALC_LJ_CUT_PRINT_FUNC>0) std::cout<<"CalcLJCut::resize(int):\n";
	if(ntypes<0) throw std::invalid_argument("CalcLJCut::resize(int): Invalid number of types.");
	ntypes_=ntypes;
	if(ntypes_>0){
		f_=Eigen::MatrixXi::Zero(ntypes_,ntypes_);
		s_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
		e_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
	}
}

void CalcLJCut::init(){
    if(CALC_LJ_CUT_PRINT_FUNC>0) std::cout<<"CalcLJCut::init():\n";
	for(int i=0; i<ntypes_; ++i){
		for(int j=i; j<ntypes_; ++j){
			if(f_(i,j)==0){
				//s_(i,j)=0.5*(s_(i,i)+s_(j,j));
				s_(i,j)=sqrt(s_(i,i)*s_(j,j));
				s_(j,i)=s_(i,j);
				e_(i,j)=sqrt(e_(i,i)*e_(j,j));
				e_(j,i)=e_(i,j);
				f_(i,j)=1; f_(j,i)=1;
			}			
		}
	}
}

void CalcLJCut::read(Token& token){
    static_cast<Calculator&>(*this).read(token);
}

void CalcLJCut::coeff(Token& token){
    if(CALC_LJ_CUT_PRINT_FUNC>0) std::cout<<"CalcLJCut::coeff(Token&):\n";
	//coeff lj_cut type1 type2 eps sigma
	const int t1=std::atoi(token.next().c_str())-1;
	const int t2=std::atoi(token.next().c_str())-1;
	const double eps=std::atof(token.next().c_str());
	const double sigma=std::atof(token.next().c_str());
	
	if(t1>=ntypes_) throw std::invalid_argument("Invalid type.");
	if(t2>=ntypes_) throw std::invalid_argument("Invalid type.");
	
	int t1min=t1,t1max=t1;
	int t2min=t2,t2max=t2;
	if(t1<0){t1min=0;t1max=ntypes_-1;}
	if(t2<0){t2min=0;t2max=ntypes_-1;}
	for(int i=t1min; i<=t1max; ++i){
		for(int j=t2min; j<=t2max; ++j){
			s_(i,j)=sigma; s_(j,i)=sigma;
			e_(i,j)=eps; e_(j,i)=eps;
			f_(i,j)=1; f_(j,i)=1;
		}
	}
}

double CalcLJCut::energy(Structure& struc, const NeighborList& nlist){
    if(CALC_LJ_CUT_PRINT_FUNC>0) std::cout<<"CalcLJCut::energy(const Structure&,const NeighborList&):\n";
	double energy=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const int ti=struc.type(i);
		for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
			const int tj=struc.type(jj);
			const Eigen::Vector3d disp=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
            const double dr2=disp.squaredNorm();
			if(dr2<rc2_){
				const double du2=s_(ti,tj)*s_(ti,tj)/dr2;
				const double du6=du2*du2*du2;
				energy+=4.0*e_(ti,tj)*du6*(du6-1.0);
			}
		}
	}
	energy*=0.5;
    struc.pe()+=energy;
	return energy;
}

double CalcLJCut::compute(Structure& struc, const NeighborList& nlist){
    if(CALC_LJ_CUT_PRINT_FUNC>0) std::cout<<"CalcLJCut::compute(const Structure&,const NeighborList&):\n";
	double energy=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const int ti=struc.type(i);
		for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
			const int tj=struc.type(jj);
			const Eigen::Vector3d disp=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
            const double dr2=disp.squaredNorm();
			if(dr2<rc2_){
				const double dri=1.0/sqrt(dr2);
				const double du=s_(ti,tj)*dri;
				const double du3=du*du*du;
				const double du6=du3*du3;
				energy+=4.0*e_(ti,tj)*du6*(du6-1.0);
				const Eigen::Vector3d fij=24.0*e_(ti,tj)*du6*(2.0*du6-1.0)*dri*dri*disp;
				struc.force(i).noalias()+=fij;
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
	
	template <> int nbytes(const CalcLJCut& obj){
		if(CALC_LJ_CUT_PRINT_FUNC>0) std::cout<<"nbytes(const CalcLJCut&):\n";
		int size=0;
		const int nt=obj.ntypes();
		size+=nbytes(static_cast<const Calculator&>(obj));
		size+=sizeof(int);//ntypes_
		size+=nt*nt*sizeof(int);//f_
		size+=nt*nt*sizeof(double);//s_
		size+=nt*nt*sizeof(double);//e_
		return size;
	}
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcLJCut& obj, char* arr){
		if(CALC_LJ_CUT_PRINT_FUNC>0) std::cout<<"pack(const CalcLJCut&,char*):\n";
		int pos=0;
		pos+=pack(static_cast<const Calculator&>(obj),arr+pos);
		std::memcpy(arr+pos,&obj.ntypes(),sizeof(int)); pos+=sizeof(int);//ntypes_
		const int nt=obj.ntypes();
		if(nt>0){
			std::memcpy(arr+pos,obj.f().data(),nt*nt*sizeof(int)); pos+=nt*nt*sizeof(int);//f_
			std::memcpy(arr+pos,obj.s().data(),nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//s_
			std::memcpy(arr+pos,obj.e().data(),nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//e_
		}
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcLJCut& obj, const char* arr){
		if(CALC_LJ_CUT_PRINT_FUNC>0) std::cout<<"unpack(CalcLJCut&,const char*):\n";
		int pos=0,nt=0;
		pos+=unpack(static_cast<Calculator&>(obj),arr+pos);
		if(obj.name()!=Calculator::Name::LJ_CUT) throw std::invalid_argument("serialize::unpack(CalcLJCut&,const char*): Invalid name.");
		std::memcpy(&nt,arr+pos,sizeof(int)); pos+=sizeof(int);//ntypes_
		obj.resize(nt);
		if(nt>0){
			std::memcpy(obj.f().data(),arr+pos,nt*nt*sizeof(int)); pos+=nt*nt*sizeof(int);//f_
			std::memcpy(obj.s().data(),arr+pos,nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//s_
			std::memcpy(obj.e().data(),arr+pos,nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//e_
		}
		obj.init();
		return pos;
	}
	
}