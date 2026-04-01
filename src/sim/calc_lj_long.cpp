// sim
#include "sim/calc_lj_long.hpp"

//==== operator ====

std::ostream& operator<<(std::ostream& out, const CalcLJLong& calc){
	return out<<static_cast<const Calculator&>(calc);
}

//==== member functions ====

void CalcLJLong::read(Token& token){
	static_cast<Calculator&>(*this).read(token);
	//calc lj_long 6.0 1e-12
	prec_=std::atof(token.next().c_str());
	if(prec_<=0) throw std::invalid_argument("CalcLJLong::read(Token&): invalid precision.");
}

void CalcLJLong::resize(int ntypes){
    if(CALC_LJ_LONG_PRINT_FUNC>0) std::cout<<"CalcLJLong::resize(int):\n";
	Calculator::resize(ntypes);
	if(ntypes_>0){
		f_=Eigen::MatrixXi::Zero(ntypes_,ntypes_);
		s_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
		e_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
	}
}

void CalcLJLong::init(){
    if(CALC_LJ_LONG_PRINT_FUNC>0) std::cout<<"CalcLJLong::init():\n";
	//compute interaction coeffecients
	for(int i=0; i<ntypes_; ++i){
		for(int j=i; j<ntypes_; ++j){
			if(f_(i,j)==0){
				s_(i,j)=sqrt(s_(i,i)*s_(j,j));
				s_(j,i)=s_(i,j);
				e_(i,j)=sqrt(e_(i,i)*e_(j,j));
				e_(j,i)=e_(i,j);
				f_(i,j)=1; f_(j,i)=1;
			}			
		}
	}
	//compute s^6
	Eigen::MatrixXd s6=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
	for(int i=0; i<ntypes_; ++i){
		for(int j=i; j<ntypes_; ++j){
			const double s3=s_(i,j)*s_(i,j)*s_(i,j);
			s6(i,j)=s3*s3;
			s6(j,i)=s6(i,j);
		}
	}
	//set k-space coefficients
	kLondon_.prec()=prec_;
	kLondon_.rc()=rc_;
	kLondon_.b()=4.0*e_.cwiseProduct(s6);
}

void CalcLJLong::init(const Structure& struc){
	kLondon_.init(struc);
}

void CalcLJLong::coeff(Token& token){
    if(CALC_LJ_LONG_PRINT_FUNC>0) std::cout<<"CalcLJLong::coeff(Token&):\n";
	//coeff type1 type2 eps sigma
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

double CalcLJLong::energy(Structure& struc, const NeighborList& nlist)const{
    if(CALC_LJ_LONG_PRINT_FUNC>0) std::cout<<"CalcLJLong::energy(const Structure&,const NeighborList&):\n";
	// k-space
	//kLondon_.init(struc);
	const double energyK=kLondon_.energy(struc);
	const double a=kLondon_.alpha();
	// r-space
	double energyR=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const int ti=struc.type(i);
		for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
			const int tj=struc.type(jj);
			const double dr2=(struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j))).squaredNorm();
			if(dr2<rc2_){
				const double du2=s_(ti,tj)*s_(ti,tj)/dr2;
				const double du6=du2*du2*du2;
				const double b2=dr2*a*a;
				energyR+=4.0*e_(ti,tj)*du6*(du6-exp(-b2)*(1.0+b2*(1.0+0.5*b2)));
			}
			
		}
	}
	energyR*=0.5;
	//total
    struc.pe()+=energyR;
	return energyR+energyK;
}

double CalcLJLong::compute(Structure& struc, const NeighborList& nlist)const{
    if(CALC_LJ_LONG_PRINT_FUNC>0) std::cout<<"CalcLJLong::compute(const Structure&,const NeighborList&):\n";
	// k-space
	//kLondon_.init(struc);
	const double energyK=kLondon_.compute(struc);
	const double a=kLondon_.alpha();
	// r-space
	double energyR=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		const int ti=struc.type(i);
		for(int j=0; j<nlist.size(i); ++j){
			const int jj=nlist.index(i,j);
			const int tj=struc.type(jj);
			const Eigen::Vector3d drv=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
            const double dr2=drv.squaredNorm();
			if(dr2<rc2_){
				const double dri2=1.0/dr2;
				const double du2=s_(ti,tj)*s_(ti,tj)*dri2;
				const double du6=du2*du2*du2;
				const double b2=dr2*a*a;
				const double expf=exp(-b2);
				energyR+=4.0*e_(ti,tj)*du6*(du6-expf*(1.0+b2*(1.0+0.5*b2)));
				const Eigen::Vector3d fij=24.0*e_(ti,tj)*du6*(2.0*du6-1.0/6.0*expf*(6.0+b2*(6.0+b2*(3.0+b2))))*dri2*drv;
				struc.force(i).noalias()+=fij;
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
	
	template <> int nbytes(const CalcLJLong& obj){
		if(CALC_LJ_LONG_PRINT_FUNC>0) std::cout<<"nbytes(const CalcLJLong&):\n";
		int size=0;
		const int nt=obj.ntypes();
		size+=nbytes(static_cast<const Calculator&>(obj));
		size+=sizeof(int);//ntypes_
		size+=sizeof(double);//prec_
		size+=nt*nt*sizeof(int);//f_
		size+=nt*nt*sizeof(double);//s_
		size+=nt*nt*sizeof(double);//e_
		return size;
	}
	
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const CalcLJLong& obj, char* arr){
		if(CALC_LJ_LONG_PRINT_FUNC>0) std::cout<<"pack(const CalcLJLong&,char*):\n";
		int pos=0;
		const int nt=obj.ntypes();
		pos+=pack(static_cast<const Calculator&>(obj),arr+pos);
		std::memcpy(arr+pos,&obj.ntypes(),sizeof(int)); pos+=sizeof(int);//ntypes_
		std::memcpy(arr+pos,&obj.prec(),sizeof(double)); pos+=sizeof(double);//prec_
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
	
	template <> int unpack(CalcLJLong& obj, const char* arr){
		if(CALC_LJ_LONG_PRINT_FUNC>0) std::cout<<"unpack(CalcLJLong&,const char*):\n";
		int pos=0,nt=0;
		pos+=unpack(static_cast<Calculator&>(obj),arr+pos);
		if(obj.name()!=Calculator::Name::LJ_LONG) throw std::invalid_argument("serialize::unpack(CalcLJLong&,const char*): Invalid name.");
		std::memcpy(&nt,arr+pos,sizeof(int)); pos+=sizeof(int);//ntypes_
		std::memcpy(&obj.prec(),arr+pos,sizeof(double)); pos+=sizeof(double);//prec_
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