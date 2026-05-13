// string
#include "str/string.hpp"
// sim
#include "sim/calc_ldamp_cut.hpp"

//==== operator ====

std::ostream& operator<<(std::ostream& out, const CalcLDampCut& calc){
	return out<<static_cast<const Calculator&>(calc)<<" mix "<<calc.mix_;
}

//==== member functions ====

void CalcLDampCut::read(Token& token){
	static_cast<Calculator&>(*this).read(token);
	//calc ldamp_long 6.0 mix 
	mix_=Calculator::Mix::read(string::to_upper(token.next()).c_str());
	if(mix_==Calculator::Mix::NONE) throw std::invalid_argument("CalcLDampCut::read(Token&): invalid mix mode.");
}

void CalcLDampCut::resize(int ntypes){
    if(CALC_LDAMP_CUT_PRINT_FUNC>0) std::cout<<"CalcLDampCut::resize(int):\n";
	Calculator::resize(ntypes);
	if(ntypes_>0){
		ie_=Eigen::VectorXd::Zero(ntypes_);
		alpha_=Eigen::VectorXd::Zero(ntypes_);
		rvdw_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
        rvdw6_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
		c6_=Eigen::MatrixXd::Zero(ntypes_,ntypes_);
	}
}

void CalcLDampCut::init(){
    if(CALC_LDAMP_CUT_PRINT_FUNC>0) std::cout<<"CalcLDampCut::init():\n";
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
}

void CalcLDampCut::coeff(Token& token){
    if(CALC_LDAMP_CUT_PRINT_FUNC>0) std::cout<<"CalcLDampCut::coeff(Token&):\n";
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
		size+=sizeof(obj.mix());//mix_
		size+=sizeof(int);//ntypes_
		size+=nt*sizeof(double);//ie_
		size+=nt*sizeof(double);//alpha_
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
		std::memcpy(arr+pos,&obj.mix(),sizeof(obj.mix())); pos+=sizeof(obj.mix());//mix_
		std::memcpy(arr+pos,&obj.ntypes(),sizeof(int)); pos+=sizeof(int);//ntypes_
		const int nt=obj.ntypes();
		if(nt>0){
			std::memcpy(arr+pos,obj.ie().data(),nt*sizeof(double)); pos+=nt*sizeof(double);//ie_
			std::memcpy(arr+pos,obj.alpha().data(),nt*sizeof(double)); pos+=nt*sizeof(double);//alpha_
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
		std::memcpy(&obj.mix(),arr+pos,sizeof(obj.mix())); pos+=sizeof(obj.mix());//mix_
		std::memcpy(&nt,arr+pos,sizeof(int)); pos+=sizeof(int);//ntypes_
		obj.resize(nt);
		if(nt>0){
			std::memcpy(obj.ie().data(),arr+pos,nt*sizeof(double)); pos+=nt*sizeof(double);//ie_
			std::memcpy(obj.alpha().data(),arr+pos,nt*sizeof(double)); pos+=nt*sizeof(double);//alpha_
			std::memcpy(obj.rvdw().data(),arr+pos,nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//s_
			std::memcpy(obj.c6().data(),arr+pos,nt*nt*sizeof(double)); pos+=nt*nt*sizeof(double);//e_
		}
		obj.init();
		return pos;
	}
	
}