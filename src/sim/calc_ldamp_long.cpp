// string
#include "str/string.hpp"
// sim
#include "sim/calc_ldamp_long.hpp"

//==== operator ====

std::ostream& operator<<(std::ostream& out, const CalcLDampLong& calc){
	return out<<static_cast<const Calculator&>(calc)<<" mix "<<calc.mix_<<" prec "<<calc.prec_;
}

//==== member functions ====

void CalcLDampLong::read(Token& token){
	static_cast<Calculator&>(*this).read(token);
	//calc ldamp_long 6.0 mix 1e-12
	mix_=Calculator::Mix::read(string::to_upper(token.next()).c_str());
	prec_=std::atof(token.next().c_str());
	if(mix_==Calculator::Mix::NONE) throw std::invalid_argument("CalcLDampLong::read(Token&): invalid mix mode.");
	if(prec_<=0) throw std::invalid_argument("CalcLDampLong::read(Token&): invalid precision.");
}

void CalcLDampLong::resize(int ntypes){
    if(CALC_LDAMP_LONG_PRINT_FUNC>0) std::cout<<"CalcLDampLong::resize(int):\n";
	Calculator::resize(ntypes);
	if(ntypes_>0){
		ie_=Eigen::VectorXd::Zero(ntypes_);
		alpha_=Eigen::VectorXd::Zero(ntypes_);
		rvdw_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
        rvdw6_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
		c6_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
	}
}

void CalcLDampLong::init(){
    if(CALC_LDAMP_LONG_PRINT_FUNC>0) std::cout<<"CalcLDampLong::init():\n";
	//compute interaction coeffecients
	for(int i=0; i<ntypes_; ++i){
		for(int j=0; j<ntypes_; ++j){
			c6_(i,j)=alpha_[i]*alpha_[j]*sqrt(ie_[i]*ie_[j]);
		}
	}
	//compute vdW radii
    for(int i=0; i<ntypes_; ++i){
		for(int j=0; j<ntypes_; ++j){
			switch(mix_){
				case Calculator::Mix::ARITHMETIC:
					rvdw_(i,j)=0.5*(rvdw_(i,i)+rvdw_(j,j));
				break;
				case Calculator::Mix::GEOMETRIC:
					rvdw_(i,j)=sqrt(rvdw_(i,i)*rvdw_(j,j));
				break;
				case Calculator::Mix::HARMONIC:
					rvdw_(i,j)=2.0*rvdw_(i,i)*rvdw_(j,j)/(rvdw_(i,i)+rvdw_(j,j));
				break;
				default:
					throw std::invalid_argument("CalcLDampLong::init(): Invalid mixing mode.");
				break;
			}
		}
	}
    //compute rvdW^6
	for(int i=0; i<ntypes_; ++i){
		for(int j=0; j<ntypes_; ++j){
			const double rvdw=rvdw_(i,j);
			const double rvdw3=rvdw*rvdw*rvdw;
			const double rvdw6=rvdw3*rvdw3;
			rvdw6_(i,j)=rvdw6;
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
	//coeff type rvdw alpha ie
	const int type=std::atoi(token.next().c_str())-1;
	const double rvdw=std::atof(token.next().c_str());
    const double alpha=std::atof(token.next().c_str());
	const double ie=std::atof(token.next().c_str());

	if(type>=ntypes_) throw std::invalid_argument("Invalid type.");

	rvdw_(type,type)=rvdw;
	alpha_[type]=alpha;
	ie_[type]=ie;
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
		size+=sizeof(obj.mix());//mix_
		size+=sizeof(int);//ntypes_
		size+=sizeof(double);//prec_
		size+=nt*sizeof(double);//ie_
		size+=nt*sizeof(double);//alpha_
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
		std::memcpy(arr+pos,&obj.mix(),sizeof(obj.mix())); pos+=sizeof(obj.mix());//mix_
		std::memcpy(arr+pos,&obj.ntypes(),sizeof(int)); pos+=sizeof(int);//ntypes_
		std::memcpy(arr+pos,&obj.prec(),sizeof(double)); pos+=sizeof(double);//prec_
		if(nt>0){
			std::memcpy(arr+pos,obj.ie().data(),nt*sizeof(double)); pos+=nt*sizeof(double);//ie_
			std::memcpy(arr+pos,obj.alpha().data(),nt*sizeof(double)); pos+=nt*sizeof(double);//alpha_
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
		std::memcpy(&obj.mix(),arr+pos,sizeof(obj.mix())); pos+=sizeof(obj.mix());//mix_
		std::memcpy(&nt,arr+pos,sizeof(int)); pos+=sizeof(int);//ntypes_
		std::memcpy(&obj.prec(),arr+pos,sizeof(double)); pos+=sizeof(double);//prec_
		obj.resize(nt);
		if(nt>0){
			std::memcpy(obj.ie().data(),arr+pos,nt*sizeof(double)); pos+=nt*sizeof(double);//ie_
			std::memcpy(obj.alpha().data(),arr+pos,nt*sizeof(double)); pos+=nt*sizeof(double);//alpha_
			std::memcpy(obj.rvdw().data(),arr+pos,nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//rvdw_
			std::memcpy(obj.c6().data(),arr+pos,nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//c6_
		}
		obj.init();
		return pos;
	}
	
}