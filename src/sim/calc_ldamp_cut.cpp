// sim
#include "sim/calc_ldamp_cut.hpp"

//==== operator ====

std::ostream& operator<<(std::ostream& out, const CalcLDampCut& calc){
	return out<<static_cast<const Calculator&>(calc);
}

//==== member functions ====

void CalcLDampCut::resize(int ntypes){
    if(CALC_LDAMP_CUT_PRINT_FUNC>0) std::cout<<"CalcLDampCut::resize(int):\n";
	Calculator::resize(ntypes);
	if(ntypes_>0){
		f_=Eigen::MatrixXi::Zero(ntypes_,ntypes_);
		rvdw_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
        rvdw6_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
		c6_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
	}
}

void CalcLDampCut::init(){
    if(CALC_LDAMP_CUT_PRINT_FUNC>0) std::cout<<"CalcLDampCut::init():\n";
	for(int i=0; i<ntypes_; ++i){
		for(int j=i+1; j<ntypes_; ++j){
			if(f_(i,j)==0){
				c6_(i,j)=sqrt(c6_(i,i)*c6_(j,j));
				c6_(j,i)=c6_(i,j);
				rvdw_(i,j)=0.5*(rvdw_(i,i)+rvdw_(j,j));
				rvdw_(j,i)=rvdw_(i,j);
				f_(i,j)=1; f_(j,i)=1;
			}
		}
	}
	for(int i=0; i<ntypes_; ++i){
		for(int j=i; j<ntypes_; ++j){
			const double rvdw=rvdw_(i,j);
			const double rvdw3=rvdw*rvdw*rvdw;
			const double rvdw6=rvdw3*rvdw3;
			rvdw6_(i,j)=rvdw6;
			rvdw6_(j,i)=rvdw6;
		}
	}
}

void CalcLDampCut::coeff(Token& token){
    if(CALC_LDAMP_CUT_PRINT_FUNC>0) std::cout<<"CalcLDampCut::coeff(Token&):\n";
	//coeff type1 type2 c6 rvdw
	const int t1=std::atof(token.next().c_str())-1;
	const int t2=std::atof(token.next().c_str())-1;
	const double c6=std::atof(token.next().c_str());
	const double rvdw=std::atof(token.next().c_str());
	
	if(t1>=ntypes_) throw std::invalid_argument("Invalid type.");
	if(t2>=ntypes_) throw std::invalid_argument("Invalid type.");
	
	int t1min=t1,t1max=t1;
	int t2min=t2,t2max=t2;
	if(t1<0){t1min=0;t1max=ntypes_-1;}
	if(t2<0){t2min=0;t2max=ntypes_-1;}
	for(int i=t1min; i<=t1max; ++i){
		for(int j=t2min; j<=t2max; ++j){
			const double rvdw3=rvdw*rvdw*rvdw;
			const double rvdw6=rvdw3*rvdw3;
			c6_(i,j)=c6; c6_(j,i)=c6;
			rvdw_(i,j)=rvdw; rvdw_(j,i)=rvdw;
			rvdw6_(i,j)=rvdw6; rvdw6_(j,i)=rvdw6;
			f_(i,j)=1; f_(j,i)=1;
		}
	}
}

double CalcLDampCut::energy(Structure& struc, const NeighborList& nlist)const{
    if(CALC_LDAMP_CUT_PRINT_FUNC>0) std::cout<<"CalcLDampCut::energy(const Structure&,const NeighborList&):\n";
	double energy=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const int ti=struc.type(i);
		for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
			const int tj=struc.type(jj);
			const double dr2=(struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j))).squaredNorm();
			if(dr2<rc2_){
				const double dr6=dr2*dr2*dr2;
				energy-=c6_(ti,tj)/(dr6+rvdw6_(ti,tj));
			}
		}
	}
	energy*=0.5;
    struc.pe()+=energy;
	return energy;
}

double CalcLDampCut::energy(Structure& struc)const{
    if(CALC_LDAMP_CUT_PRINT_FUNC>0) std::cout<<"CalcLDampCut::energy(const Structure&):\n";
	double energy=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const int ti=struc.type(i);
		for(int j=i+1; j<struc.nAtoms(); ++j){
			const int tj=struc.type(j);
			const double dr2=(struc.posn(i)-struc.posn(j)).squaredNorm();
			if(dr2<rc2_){
				const double dr6=dr2*dr2*dr2;
				energy-=c6_(ti,tj)/(dr6+rvdw6_(ti,tj));
			}
		}
	}
	struc.pe()+=energy;
	return energy;
}

double CalcLDampCut::compute(Structure& struc, const NeighborList& nlist)const{
    if(CALC_LDAMP_CUT_PRINT_FUNC>0) std::cout<<"CalcLDampCut::compute(const Structure&,const NeighborList&):\n";
	double energy=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const int ti=struc.type(i);
		for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
			const int tj=struc.type(jj);
			const Eigen::Vector3d drv=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
            const double dr2=drv.squaredNorm();
			if(dr2<rc2_){
				const double dr6=dr2*dr2*dr2;
				const double den=1.0/(dr6+rvdw6_(ti,tj));
				energy-=c6_(ti,tj)*den;
				struc.force(i).noalias()-=6.0*c6_(ti,tj)*dr2*dr2*den*den*drv;
			}
		}
	}
	energy*=0.5;
    struc.pe()+=energy;
	return energy;
}

double CalcLDampCut::compute(Structure& struc)const{
    if(CALC_LDAMP_CUT_PRINT_FUNC>0) std::cout<<"CalcLDampCut::compute(const Structure&):\n";
	double energy=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const int ti=struc.type(i);
		for(int j=i+1; j<struc.nAtoms(); ++j){
			const int tj=struc.type(j);
			const Eigen::Vector3d drv=struc.posn(i)-struc.posn(j);
            const double dr2=drv.squaredNorm();
			if(dr2<rc2_){
				const double dr6=dr2*dr2*dr2;
				const double den=1.0/(dr6+rvdw6_(ti,tj));
				energy-=c6_(ti,tj)*den;
				const Eigen::Vector3d fij=-6.0*c6_(ti,tj)*dr2*dr2*den*den*drv;
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
	
	template <> int nbytes(const CalcLDampCut& obj){
		if(CALC_LDAMP_CUT_PRINT_FUNC>0) std::cout<<"nbytes(const CalcLDampCut&):\n";
		int size=0;
		const int nt=obj.ntypes();
		size+=nbytes(static_cast<const Calculator&>(obj));
		size+=sizeof(int);//ntypes_
		size+=nt*nt*sizeof(int);//f_
		size+=nt*nt*sizeof(double);//rvdw_
		size+=nt*nt*sizeof(double);//c6_
		return size;
	}
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcLDampCut& obj, char* arr){
		if(CALC_LDAMP_CUT_PRINT_FUNC>0) std::cout<<"pack(const CalcLDampCut&,char*):\n";
		int pos=0;
		pos+=pack(static_cast<const Calculator&>(obj),arr+pos);
		std::memcpy(arr+pos,&obj.ntypes(),sizeof(int)); pos+=sizeof(int);//ntypes_
		const int nt=obj.ntypes();
		if(nt>0){
			std::memcpy(arr+pos,obj.f().data(),nt*nt*sizeof(int)); pos+=nt*nt*sizeof(int);//f_
			std::memcpy(arr+pos,obj.rvdw().data(),nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//s_
			std::memcpy(arr+pos,obj.c6().data(),nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//e_
		}
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcLDampCut& obj, const char* arr){
		if(CALC_LDAMP_CUT_PRINT_FUNC>0) std::cout<<"unpack(CalcLDampCut&,const char*):\n";
		int pos=0,nt=0;
		pos+=unpack(static_cast<Calculator&>(obj),arr+pos);
		if(obj.name()!=Calculator::Name::LDAMP_CUT) throw std::invalid_argument("serialize::unpack(CalcLDampCut&,const char*): Invalid name.");
		std::memcpy(&nt,arr+pos,sizeof(int)); pos+=sizeof(int);//ntypes_
		obj.resize(nt);
		if(nt>0){
			std::memcpy(obj.f().data(),arr+pos,nt*nt*sizeof(int)); pos+=nt*nt*sizeof(int);//f_
			std::memcpy(obj.rvdw().data(),arr+pos,nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//s_
			std::memcpy(obj.c6().data(),arr+pos,nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//e_
		}
		obj.init();
		return pos;
	}
	
}