// sim
#include "sim/calc_ldamp_long.hpp"

//==== operator ====

std::ostream& operator<<(std::ostream& out, const CalcLDampLong& calc){
	return out<<static_cast<const Calculator&>(calc);
}

//==== member functions ====

void CalcLDampLong::read(Token& token){
	static_cast<Calculator&>(*this).read(token);
	//calc lj_long 6.0 1e-12
	prec_=std::atof(token.next().c_str());
	if(prec_<=0) throw std::invalid_argument("CalcLDampLong::read(Token&): invalid precision.");
}

void CalcLDampLong::resize(int ntypes){
    if(CALC_LDAMP_LONG_PRINT_FUNC>0) std::cout<<"CalcLDampLong::resize(int):\n";
	Calculator::resize(ntypes);
	if(ntypes_>0){
		f_=Eigen::MatrixXi::Zero(ntypes_,ntypes_);
		rvdw_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
        rvdw6_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
		c6_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
	}
}

void CalcLDampLong::init(){
    if(CALC_LDAMP_LONG_PRINT_FUNC>0) std::cout<<"CalcLDampLong::init():\n";
	//compute interaction coeffecients
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
    //compute rvdw^6
	for(int i=0; i<ntypes_; ++i){
		for(int j=i; j<ntypes_; ++j){
			const double rvdw=rvdw_(i,j);
			const double rvdw3=rvdw*rvdw*rvdw;
			const double rvdw6=rvdw3*rvdw3;
			rvdw6_(i,j)=rvdw6;
			rvdw6_(j,i)=rvdw6;
		}
	}
	//set k-space coefficients
	kLondon_.prec()=prec_;
	kLondon_.rc()=rc_;
	kLondon_.b()=c6_;
}

void CalcLDampLong::init(const Structure& struc){
	kLondon_.init(struc);
}

void CalcLDampLong::coeff(Token& token){
    if(CALC_LDAMP_LONG_PRINT_FUNC>0) std::cout<<"CalcLDampLong::coeff(Token&):\n";
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

double CalcLDampLong::energy(Structure& struc, const NeighborList& nlist)const{
    if(CALC_LDAMP_LONG_PRINT_FUNC>0) std::cout<<"CalcLDampLong::energy(const Structure&,const NeighborList&):\n";
	//compute k-space energy
	//kLondon_.init(struc);
	const double energyK=kLondon_.energy(struc);
	const double a2=kLondon_.alpha()*kLondon_.alpha();
	//compute r-space energy
	double energyR=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const int ti=struc.type(i);
		for(int j=0; j<nlist.size(i); ++j){
            const int jj=nlist.index(i,j);
			const int tj=struc.type(jj);
            const double dr2=(struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j))).squaredNorm();
			if(dr2<rc2_){
                const double dr6=dr2*dr2*dr2;
				const double b2=a2*dr2;
				energyR+=c6_(ti,tj)*(
					(1.0-exp(-b2)*(1.0+b2*(1.0+0.5*b2)))/dr6
					-1.0/(dr6+rvdw6_(ti,tj))
				);
			}
		}
	}
	energyR*=0.5;
	//total
    struc.pe()+=energyR;
	return energyR+energyK;
}

double CalcLDampLong::compute(Structure& struc, const NeighborList& nlist)const{
    if(CALC_LDAMP_LONG_PRINT_FUNC>0) std::cout<<"CalcLDampLong::compute(const Structure&,const NeighborList&):\n";
	// k-space
	//kLondon_.init(struc);
	const double energyK=kLondon_.compute(struc);
	const double a2=kLondon_.alpha()*kLondon_.alpha();
	//compute r-space energy
	double energyR=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const int ti=struc.type(i);
		for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
			const int tj=struc.type(jj);
			const Eigen::Vector3d drv=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
            const double dr2=drv.squaredNorm();
			if(dr2<rc2_){
				const double dr6=dr2*dr2*dr2;
				const double b2=a2*dr2;
				const double expf=exp(-b2);
				const double den=1.0/(dr6+rvdw6_(ti,tj));
                energyR+=c6_(ti,tj)*(
					(1.0-expf*(1.0+b2*(1.0+0.5*b2)))/dr6
					-den
				);
                struc.force(i).noalias()+=drv*c6_(ti,tj)*(
                    (6.0-expf*(6.0+b2*(6.0+b2*(3.0+b2))))/(dr6*dr2)
					-6.0*dr2*dr2*den*den
                );
			}
		}
	}
    energyR*=0.5;
	// total
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
	
	template <> int nbytes(const CalcLDampLong& obj){
		if(CALC_LDAMP_LONG_PRINT_FUNC>0) std::cout<<"nbytes(const CalcLDampLong&):\n";
		int size=0;
		const int nt=obj.ntypes();
		size+=nbytes(static_cast<const Calculator&>(obj));
		size+=sizeof(int);//ntypes_
		size+=sizeof(double);//prec_
		size+=nt*nt*sizeof(int);//f_
		size+=nt*nt*sizeof(double);//rvdw_
		size+=nt*nt*sizeof(double);//c6_
		return size;
	}
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcLDampLong& obj, char* arr){
		if(CALC_LDAMP_LONG_PRINT_FUNC>0) std::cout<<"pack(const CalcLDampLong&,char*):\n";
		int pos=0;
		const int nt=obj.ntypes();
		pos+=pack(static_cast<const Calculator&>(obj),arr+pos);
		std::memcpy(arr+pos,&obj.ntypes(),sizeof(int)); pos+=sizeof(int);//ntypes_
		std::memcpy(arr+pos,&obj.prec(),sizeof(double)); pos+=sizeof(double);//prec_
		if(nt>0){
			std::memcpy(arr+pos,obj.f().data(),nt*nt*sizeof(int)); pos+=nt*nt*sizeof(int);//f_
			std::memcpy(arr+pos,obj.rvdw().data(),nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//rvdw_
			std::memcpy(arr+pos,obj.c6().data(),nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//c6_
		}
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(CalcLDampLong& obj, const char* arr){
		if(CALC_LDAMP_LONG_PRINT_FUNC>0) std::cout<<"unpack(CalcLDampLong&,const char*):\n";
		int pos=0,nt=0;
		pos+=unpack(static_cast<Calculator&>(obj),arr+pos);
		if(obj.name()!=Calculator::Name::LDAMP_LONG) throw std::invalid_argument("serialize::unpack(CalcLDampLong&,const char*): Invalid name.");
		std::memcpy(&nt,arr+pos,sizeof(int)); pos+=sizeof(int);//ntypes_
		std::memcpy(&obj.prec(),arr+pos,sizeof(double)); pos+=sizeof(double);//prec_
		obj.resize(nt);
		if(nt>0){
			std::memcpy(obj.f().data(),arr+pos,nt*nt*sizeof(int)); pos+=nt*nt*sizeof(int);//f_
			std::memcpy(obj.rvdw().data(),arr+pos,nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//rvdw_
			std::memcpy(obj.c6().data(),arr+pos,nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//c6_
		}
		obj.init();
		return pos;
	}
	
}