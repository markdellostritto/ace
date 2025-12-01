// str
#include "str/print.hpp"
// chem
#include "chem/units.hpp"
#include "chem/ptable.hpp"
// sim
#include "sim/thole.hpp"

//==== operators ====
std::ostream& operator<<(std::ostream& out, Thole& thole){
	char* str=new char[print::len_buf];
	out<<print::buf(str)<<"\n";
	out<<print::title("THOLE",str)<<"\n";
	out<<"inter  = "<<thole.inter_<<"\n";
	out<<"alphar = "<<thole.alphar_<<"\n";
	out<<"calc   = "<<thole.calc_<<"\n";
	for(int i=0; i<thole.calc_.alpha().size(); ++i){
		out<<"alpha["<<i<<"] = "<<thole.calc_.alpha()[i]<<"\n";
	}
	for(int i=0; i<thole.alpha_.size(); ++i){
		out<<"alpha["<<i<<"] = "<<thole.alpha_[i]<<"\n";
	}
	out<<print::buf(str);
	delete[] str;
	return out;
}

//==== member functions ====

Eigen::Matrix3d& Thole::compute(const Structure& struc, Eigen::Matrix3d& alphaT){
	const double ke=units::Consts::ke();
    //set the utility matrices
    const int nAtoms=struc.nAtoms();
    const int size=3*nAtoms;
	A_=Eigen::MatrixXd::Zero(size,size);
	alphaC_.resize(size,3);
	identityC_.resize(size,3);
    r_.resize(nAtoms);
	r0_.resize(nAtoms);
	drAtom_.resize(nAtoms);

	//units
	double rscale=0.0;
	if(units::Consts::system()==units::System::AU) rscale=units::Bohr2Ang;
	else if(units::Consts::system()==units::System::METAL) rscale=1.0;
	else throw std::runtime_error("Invalid units.");

    //set the gas-phase atomic radius
	if(DEBUG_THOLE>1) std::cout<<"Setting gas-phase atomic radii.\n";
	for(int i=0; i<nAtoms; ++i) r0_[i]=ptable::radius_covalent(struc.an(i))*rscale;
	for(int i=0; i<nAtoms; ++i) r_[i]=ptable::radius_covalent(struc.an(i))*rscale;
	
    //calculate the atom-in-molecule radius
	if(DEBUG_THOLE>1) std::cout<<"Calculating the atom-in-molecule radius...\n";
	for(int i=0; i<nAtoms; ++i){
		double rMin=struc.R().squaredNorm();
		for(int j=0; j<nAtoms; ++j){
			if(i==j) continue;
            Eigen::Vector3d rVec;
			struc.diff(struc.posn(j),struc.posn(i),rVec);
			double dr=rVec.norm();
			double rTemp=r0_[i]-0.5*(r0_[i]+r0_[j]-dr);
			if(rTemp<rMin){
                Eigen::Vector3d r00,rNew;
				r_[i]=rTemp;
				rMin=rTemp;
				rVec/=rVec.norm();
				r00.noalias()=r0_[i]*rVec;
				rNew.noalias()=r_[i]*rVec;
				drAtom_[i][0]=std::fabs(rNew[0]/r00[0]);
				drAtom_[i][1]=std::fabs(rNew[1]/r00[1]);
				drAtom_[i][2]=std::fabs(rNew[2]/r00[2]);
				if(rVec[0]==0 && r00[0]==0) drAtom_[i][0]=1;
				if(rVec[1]==0 && r00[1]==0) drAtom_[i][1]=1;
				if(rVec[2]==0 && r00[2]==0) drAtom_[i][2]=1;
			}
		}
	}
	
    //populate the solution vector
	if(DEBUG_THOLE>1) std::cout<<"Populating solution vector...\n";
	for(int i=0; i<nAtoms; ++i) identityC_.block<3,3>(i*3,0)=Eigen::Matrix3d::Identity();
    
    //populate the diagonal of the A-matrix
	if(DEBUG_THOLE>1) std::cout<<"Populating diagonal of A-matrix...\n";
	if(alphar_){
		for(int i=0; i<nAtoms; ++i){
			Eigen::Matrix3d alpha=Eigen::Matrix3d::Identity()*alpha_[struc.type(i)];
			alpha(0,0)*=drAtom_[i][0];
			alpha(1,1)*=drAtom_[i][1];
			alpha(2,2)*=drAtom_[i][2];
			A_.block<3,3>(i*3,i*3).noalias()=alpha.inverse()*ke;
		}
	} else {
		for(int i=0; i<nAtoms; ++i){
			Eigen::Matrix3d alpha=Eigen::Matrix3d::Identity()*alpha_[struc.type(i)];
			A_.block<3,3>(i*3,i*3).noalias()=alpha.inverse()*ke;
		}
	}

    //calculate the interaction matrices
	if(DEBUG_THOLE>1) std::cout<<"Calculating the off-diagonal elements...\n";
	if(inter_){
		Eigen::Vector3d drv=Eigen::Vector3d::Zero();
        Eigen::Matrix3d interMat=Eigen::Matrix3d::Identity()*calc_.kDipole().vc()*ke/3.0;
		for(int i=0; i<nAtoms; ++i){
			A_.block<3,3>(i*3,i*3).noalias()-=interMat;
		}
		for(int i=0; i<nAtoms; ++i){
			for(int j=i+1; j<nAtoms; ++j){
                //calculate distance
				struc.diff(struc.posn(j),struc.posn(i),drv);
                //compute interaction matrix
				A_.block<3,3>(j*3,i*3).noalias()-=calc_.interMat(drv,interMat,calc_.kDipole().alpha());
			}
		}
	}

    //set the total matrix
	A_=A_.selfadjointView<Eigen::Lower>();

    //calculate the effective polarizabilities
	if(DEBUG_THOLE>1) std::cout<<"Calculating effective polarizabilities...\n";
    alphaC_.noalias()=A_.llt().solve(identityC_);

    //record the polarizability
    alphaT.setZero();
	for(int i=0; i<nAtoms; ++i){
		alphaT.noalias()+=alphaC_.block<3,3>(i*3,0);
	}

    //return the polarizability
    return alphaT;
}