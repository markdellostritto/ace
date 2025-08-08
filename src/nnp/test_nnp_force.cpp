// c
#include <cstdlib>
#include <ctime>
// c++
#include <iostream>
// struc
#include "struc/structure.hpp"
// nnp
#include "nnp/nnp.hpp"

Structure& make_random(Structure& struc, const std::vector<NNH::Type>& types){
	//resize struc
	const int natoms=20;
	Atom atom; 
	atom.name=true; atom.an=true; atom.type=true;
	atom.charge=false; atom.posn=true; atom.force=true; 
    atom.symm=true;
	struc.resize(natoms,atom);
	
	//resize the lattice
	const double a0=20.0;
	Eigen::Matrix3d lv=Eigen::Matrix3d::Identity()*a0;
	struc.init(lv);
	
	//set the atomic properties
	const int ntypes=types.size();
	for(int i=0; i<natoms; ++i){
		const int type=std::rand()%ntypes;
		//const int type=1;
		struc.type(i)=type;
		struc.name(i)=types[type].name();
		struc.posn(i)=Eigen::Vector3d::Random()*a0;
		Cell::returnToCell(struc.posn(i),struc.posn(i),struc.R(),struc.RInv());
	}
	
	//return 
	return struc;
}

void compute_force_num(Structure& struc, NNP& nnp){
	const double eps=1.0e-6;
	Structure struc_p=struc;
	Structure struc_m=struc;
	
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.force(i).setZero();
	}
	
	for(int i=0; i<struc.nAtoms(); ++i){
		//reset positions
		for(int n=0; n<struc.nAtoms(); ++n){
			struc_p.posn(n)=struc.posn(n);
			struc_m.posn(n)=struc.posn(n);
		}
        //compute force on atom i
		for(int j=0; j<3; ++j){
			struc_p.posn(i)=struc.posn(i); struc_p.posn(i)[j]+=eps;
			struc_m.posn(i)=struc.posn(i); struc_m.posn(i)[j]-=eps;
            nnp.compute_symm(struc_p);
			nnp.compute_symm(struc_m);
            const double ep=nnp.compute_energy(struc_p);
            const double em=nnp.compute_energy(struc_m);
			struc.force(i)[j]=-0.5*(ep-em)/eps;
		}
	}
}

int main(int argc, char* argv[]){
	
	//set the types
	const int ntypes=3;
	std::vector<NNH::Type> types(ntypes);
	types[0].name()="Ar";
	types[1].name()="Ne";
	types[2].name()="Xe";
	std::cout<<"Types = \n";
	for(int n=0; n<ntypes; ++n){
		std::cout<<types[n]<<"\n";
	}

	//==== make the structure ====
	std::cout<<"making structures\n";
	Structure struc_a;
	Structure struc_n;
	make_random(struc_a,types);
	struc_n=struc_a;
	std::cout<<struc_a<<"\n";
	for(int i=0; i<struc_a.nAtoms(); ++i){
		std::cout<<"atom "<<i<<" "<<struc_a.name(i)<<" "<<struc_a.type(i)<<" "<<struc_a.posn(i).transpose()<<"\n";
	}
    
    //nnp
	const double rcut=6.0;
	NNP nnp(types);
	nnp.rc()=rcut;
	NN::ANNP init;
	std::vector<int> nh(3);
	nh[0]=5; nh[1]=3; nh[2]=2;
	init.init()=NN::Init::HE;
	init.seed()=-1;
	init.neuron()=NN::Neuron::SWISH;

	//initialize basis
	std::cout<<"initializing basis\n";
	const int nR=6;
	BasisR basisR(rcut,Cutoff::Name::COS,nR,BasisR::Name::LOGCOSH);
    const int na=2;
	const int nA=2*na;
	BasisA basisA(rcut,Cutoff::Name::COS,nA,BasisA::Name::GAUSS);
	for(int i=0; i<nR; ++i){
		basisR.rs(i)=1.54;
		basisR.eta(i)=i+1;
	}
	for(int i=0; i<na; ++i){
		basisA.eta(i)=2.5279;
		basisA.zeta(i)=1.1390;
		basisA.lambda(i)=1;
	}
    for(int i=na; i<nA; ++i){
		basisA.eta(i)=2.5279;
		basisA.zeta(i)=1.1390;
		basisA.lambda(i)=-1;
	}

    //initialize neural network
	std::cout<<"initializing neural network\n";
	for(int i=0; i<ntypes; ++i){
		for(int j=0; j<ntypes; ++j){
			nnp.nnh(i).basisR(j)=basisR;
		}
		for(int j=0; j<ntypes; ++j){
			for(int k=j; k<ntypes; ++k){
				nnp.nnh(i).basisA(j,k)=basisA;
			}
		}
		nnp.nnh(i).init_input();
		nnp.nnh(i).nn().resize(init,nnp.nnh(i).nInput(),nh,1);
        nnp.nnh(i).dOdZ().resize(nnp.nnh(i).nn());
    }
    std::cout<<nnp<<"\n";
	std::cout<<basisR<<"\n";
	std::cout<<basisA<<"\n";

    //==== compute the analytical forces ====
	std::cout<<"computing forces - analytical\n";
	nnp.init(struc_a);
	nnp.compute_symm(struc_a);
	std::cout<<"pe = "<<nnp.compute_energy(struc_a)<<"\n";
	nnp.compute_force(struc_a);
    for(int i=0; i<struc_a.nAtoms(); ++i){
		std::cout<<"fa["<<i<<"] = "<<struc_a.force(i).transpose()<<"\n";
	}

    //==== compute the numerical forces ====
	std::cout<<"computing forces - numerical\n";
	nnp.init(struc_n);
	compute_force_num(struc_n,nnp);
	for(int i=0; i<struc_n.nAtoms(); ++i){
		std::cout<<"fn["<<i<<"] = "<<struc_n.force(i).transpose()<<"\n";
	}

	double error=0;
	for(int i=0; i<struc_a.nAtoms(); ++i){
		error+=(struc_a.force(i)-struc_n.force(i)).squaredNorm();
	}
	error=std::sqrt(error/struc_a.nAtoms());
	std::cout<<"error = "<<error<<"\n";
	
    return 0;
}