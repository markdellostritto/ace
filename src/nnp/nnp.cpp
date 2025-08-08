// c libraries
#if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
#include <cmath>
#elif (defined __ICC || defined __INTEL_COMPILER)
#include <mathimf.h> //intel math library
#else
#include <cmath>
#endif
// c++ libraries
#include <iostream>
// structure
#include "struc/structure.hpp"
// math
#include "math/const.hpp"
// str
#include "str/print.hpp"
#include "str/token.hpp"
// nnp
#include "nnp/nnh.hpp"
#include "nnp/nnp.hpp"

//************************************************************
// NNP - Neural Network Potential
//************************************************************

//==== operators ====

/**
* print the nnp to screen
* @param out - output stream
* @param nnp - the neural network potential
*/
std::ostream& operator<<(std::ostream& out, const NNP& nnp){
	char* str=new char[print::len_buf];
	out<<print::buf(str)<<"\n";
	out<<print::title("NN - POT",str)<<"\n";
	out<<"r_cut  = "<<nnp.rc_<<"\n";
	out<<"ntypes = "<<nnp.ntypes_<<"\n";
	for(int i=0; i<nnp.ntypes_; ++i) std::cout<<nnp.nnh_[i]<<"\n";
	out<<print::buf(str);
	delete[] str;
	return out;
}

//==== member functions ====

int NNP::size()const{
	int size=0;
	for(int i=0; i<ntypes_; ++i){
		size+=nnh_[i].nn().size();
	}
	return size;
}

/**
* set defaults for the neural network potential
*/
void NNP::defaults(){
	if(NNP_PRINT_FUNC>0) std::cout<<"NNP::defaults():\n";
	//types
		ntypes_=0;
		map_.clear();
		nnh_.clear();
	//cutoff
		rc_=0;
}

//==== resizing ====

/**
* resize the number of types and each NNH
* @param types - the types of the neural network potential
*/
void NNP::resize(const std::vector<NNH::Type>& types){
	if(NNP_PRINT_FUNC>0) std::cout<<"NNP::resize(const std::vector<NNH::Type>&)\n";
	if(types.size()<0) throw std::invalid_argument("NNP::resize(const std::vector<NNH::Type>&): invalid number of types");
	ntypes_=types.size();
	nnh_.resize(ntypes_);
	for(int i=0; i<ntypes_; ++i){
		nnh_[i].resize(ntypes_);
		nnh_[i].type()=types[i];
		map_.add(string::hash(types[i].name()),i);
	}
}

//==== static functions ====

//read/write basis

/**
* Read the basis for a given types from file.
* @param file - the name of the file from which the object will be read
* @param nnp - the neural network potential to be written
* @param name - the types for which we will read the basis
*/
void NNP::read_basis(const char* file, NNP& nnp, const char* name){
	if(NNP_PRINT_FUNC>0) std::cout<<"NNP::read_basis(const char*,NNP&,const char*):\n";
	FILE* reader=NULL;
	reader=fopen(file,"r");
	if(reader!=NULL){
		NNP::read_basis(reader,nnp,name);
		fclose(reader);
		reader=NULL;
	} else throw std::runtime_error(std::string("NNP::read_basis(const char*,NNP&,const char*): Could not open nnp file: \"")+std::string(file)+std::string("\""));
}

/**
* Read the basis for a given types from file.
* @param reader - the file pointer from which the object will be read
* @param nnp - the neural network potential to be written
* @param name - the types for which we will read the basis
*/
void NNP::read_basis(FILE* reader, NNP& nnp, const char* name){
	if(NNP_PRINT_FUNC>0) std::cout<<"NNP::read_basis(FILE*,NNP&,const char*):\n";
	//==== local function variables ====
	Token token;
	char* input=new char[string::M];
	//==== get atom types ====
	const int atomIndex=nnp.index(name);
	//==== global cutoff ====
	token.read(fgets(input,string::M,reader),string::WS); token.next();
	const double rc=std::atof(token.next().c_str());
	if(rc!=nnp.rc()) throw std::invalid_argument("NNP::read_basis(FILE*,NNP&,const char*): invalid cutoff.");
	//==== number of types ====
	token.read(fgets(input,string::M,reader),string::WS); token.next();
	const int ntypes=std::atoi(token.next().c_str());
	if(ntypes!=nnp.ntypes()) throw std::invalid_argument("NNP::read_basis(FILE*,NNP&,const char*): invalid number of types.");
	//==== central types ====
	token.read(fgets(input,string::M,reader),string::WS); token.next();
	const int II=nnp.index(token.next());
	//==== check indices ====
	if(atomIndex!=II) throw std::invalid_argument("NNP::read_basis(FILE*,NNP&,const char*): invalid central types.\n");
	//==== basis - radial ====
	for(int j=0; j<ntypes; ++j){
		//read types
		token.read(fgets(input,string::M,reader),string::WS); token.next();
		const int JJ=nnp.index(token.next());
		//read basis
		BasisR::read(reader,nnp.nnh(II).basisR(JJ));
	}
	//==== basis - angular ====
	for(int j=0; j<ntypes; ++j){
		for(int k=j; k<ntypes; ++k){
			//read types
			token.read(fgets(input,string::M,reader),string::WS); token.next();
			const int JJ=nnp.index(token.next());
			const int KK=nnp.index(token.next());
			//read basis
			BasisA::read(reader,nnp.nnh(II).basisA(JJ,KK));
		}
	}
	//==== initialize the inputs ====
	nnp.nnh(II).init_input();
	//==== clear local variables ====
	delete[] input;
}

//read/write nnp

/**
* Write the neural network to file
* @param file - the name of the file to which the object will be written
* @param nnp - the neural network potential to be written
*/
void NNP::write(const char* file, const NNP& nnp){
	if(NNP_PRINT_FUNC>0) std::cout<<"NNP::write(const char*,const NNP&):\n";
	FILE* writer=NULL;
	writer=fopen(file,"w");
	if(writer!=NULL){
		NNP::write(writer,nnp);
		fclose(writer);
		writer=NULL;
	} else throw std::runtime_error(std::string("NNP::write(const char*,const NNP&): Could not write to nnh file: \"")+std::string(file)+std::string("\""));
}

/**
* Read the neural network from file
* @param file - the name of the file fro
* @param nnp - stores the neural network potential to be read
*/
void NNP::read(const char* file, NNP& nnp){
	if(NNP_PRINT_FUNC>0) std::cout<<"NNP::read(const char*,NNP&):\n";
	FILE* reader=NULL;
	reader=fopen(file,"r");
	if(reader!=NULL){
		NNP::read(reader,nnp);
		fclose(reader);
		reader=NULL;
	} else throw std::runtime_error(std::string("NNP::read(const char*,NNP&): Could not open nnp file: \"")+std::string(file)+std::string("\""));
}

/**
* Write the neural network to file
* @param writer - the file pointer used to write the object to file
* @param nnp - the neural network potential to be written
*/
void NNP::write(FILE* writer, const NNP& nnp){
	if(NNP_PRINT_FUNC>0) std::cout<<"NNP::write(FILE*,const NNP&):\n";
	//==== header ====
	fprintf(writer,"ann\n");
	//==== types ====
	fprintf(writer, "ntypes %i\n",nnp.ntypes());
	for(int n=0; n<nnp.ntypes(); ++n){
		nnp.nnh(n).type().write(writer);
	}
	//==== cutoff ====
	fprintf(writer,"rc %f\n",nnp.rc());
	//==== basis ====
	for(int i=0; i<nnp.ntypes(); ++i){
		//write central types
		fprintf(writer,"basis %s\n",nnp.nnh(i).type().name().c_str());
		//write basis - radial
		for(int j=0; j<nnp.ntypes(); ++j){
			//write types
			fprintf(writer,"basis_radial %s\n",nnp.nnh(j).type().name().c_str());
			//write basis
			BasisR::write(writer,nnp.nnh(i).basisR(j));
		}
		//write basis - angular
		for(int j=0; j<nnp.ntypes(); ++j){
			for(int k=j; k<nnp.ntypes(); ++k){
				//write types
				fprintf(writer,"basis_angular %s %s\n",nnp.nnh(j).type().name().c_str(),nnp.nnh(k).type().name().c_str());
				//write basis
				BasisA::write(writer,nnp.nnh(i).basisA(j,k));
			}
		}
	}
	//==== neural network ====
	for(int n=0; n<nnp.ntypes(); ++n){
		//write central types
		fprintf(writer,"nn %s\n",nnp.nnh(n).type().name().c_str());
		//write the network
		NN::ANN::write(writer,nnp.nnh(n).nn());
	}
}

/**
* Read the neural network from file
* @param reader - the file pointer used to read the object from file
* @param nnp - stores the neural network potential to be read
*/
void NNP::read(FILE* reader, NNP& nnp){
	if(NNP_PRINT_FUNC>0) std::cout<<"NNP::read(FILE*,NNP&):\n";
	//==== local function variables ====
	Token token;
	char* input=new char[string::M];
	//==== header ====
	fgets(input,string::M,reader);
	//==== number of types ====
	token.read(fgets(input,string::M,reader),string::WS); token.next();
	const int ntypes=std::atoi(token.next().c_str());
	if(ntypes<=0) throw std::invalid_argument("NNP::read(FILE*,NNP&): invalid number of types.");
	//==== types ====
	std::vector<NNH::Type> types(ntypes);
	for(int n=0; n<ntypes; ++n){
		types[n].read(fgets(input,string::M,reader));
	}
	//==== resize ====
	nnp.resize(types);
	//==== global cutoff ====
	token.read(fgets(input,string::M,reader),string::WS); token.next();
	const double rc=std::atof(token.next().c_str());
	if(rc<=0) throw std::invalid_argument("NNP::read(FILE*,NNP&): invalid cutoff.");
	else nnp.rc()=rc;
	//==== basis ====
	for(int i=0; i<ntypes; ++i){
		//read central types
		token.read(fgets(input,string::M,reader),string::WS); token.next();
		const int II=nnp.index(token.next());
		//read basis - radial
		for(int j=0; j<ntypes; ++j){
			//read types
			token.read(fgets(input,string::M,reader),string::WS); token.next();
			const int JJ=nnp.index(token.next());
			//read basis
			BasisR::read(reader,nnp.nnh(II).basisR(JJ));
		}
		//read basis - angular
		for(int j=0; j<ntypes; ++j){
			for(int k=j; k<ntypes; ++k){
				//read types
				token.read(fgets(input,string::M,reader),string::WS); token.next();
				const int JJ=nnp.index(token.next());
				const int KK=nnp.index(token.next());
				//read basis
				BasisA::read(reader,nnp.nnh(II).basisA(JJ,KK));
			}
		}
	}
	//==== initialize inputs ====
	for(int i=0; i<ntypes; ++i){
		nnp.nnh(i).init_input();
	}
	//==== neural network ====
	for(int n=0; n<ntypes; ++n){
		//read types
		token.read(fgets(input,string::M,reader),string::WS); token.next();
		const int II=nnp.index(token.next());
		//read network
		NN::ANN::read(reader,nnp.nnh(II).nn());
		//resize gradient object
		nnp.nnh(II).dOdZ().resize(nnp.nnh(II).nn());
	}
	//==== clear local variables ====
	delete[] input;
}

//calculation

/**
* resize the symmetry function vectors to store the inputs
* @param nnp - the neural network potential
* @param struc - the structure which we will compute
*/
void NNP::init(Structure& struc)const{
	if(NNP_PRINT_FUNC>0) std::cout<<"NNP::init(Structure&):\n";
	for(int n=0; n<struc.nAtoms(); ++n){
		struc.symm(n).resize(nnh(index(struc.name(n))).nInput());
	}
}

/**
* execute all atomic networks and return energy
* @param struc - the structure for which we will compute the energy
* @return total energy
* it is assumed that the symmetry functions have been computed
*/
double NNP::compute_energy(const Structure& struc){
	if(NNP_PRINT_FUNC>0) std::cout<<"NNP::energy(const NNP&,const Structure&):\n";
	double pe=0;
	//loop over atoms
	for(int i=0; i<struc.nAtoms(); ++i){
		const int II=index(struc.name(i));
		pe+=nnh_[II].energy(struc.symm(i));
	}
	return pe;
}

/**
* Compute the symmetry functions for a given structure.
* @param nnp - the neural network potential
* @param struc - the structure which we will compute
* @param nlist - neighbor list for each atom (includes periodic images)
* The symmetry functions are computed by looping over all nearest-neighbor pairs
* and all unique nearest-neighbor triples.  Thus, the neighbor list must be set 
* for the structure, and the type must correspond to the index of the atomic types
* in the NNP.
*/
void NNP::compute_symm(Structure& struc, const NeighborList& nlist){
	if(NNP_PRINT_FUNC>0) std::cout<<"NNP::compute_symm(Structure&,const NeighborList&):\n";
	for(int i=0; i<struc.nAtoms(); ++i){
		//reset the inputs
		struc.symm(i).setZero();
		//get the index of types i
		const int II=struc.type(i);
		//loop over all neighbors
		for(int j=0; j<nlist.size(i); ++j){
			//get the index of types j
			const int jj=nlist.index(i,j);
			const int JJ=struc.type(jj);
			//get the distance from J to I
			const Eigen::Vector3d rIJ=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
			const double drIJ=rIJ.norm();
			//compute the IJ contribution to all radial basis functions
			nnh_[II].basisR(JJ).symm(drIJ);
			const int oR=nnh_[II].offsetR(JJ);
			const int nR=nnh_[II].basisR(JJ).size();
			for(int nr=0; nr<nR; ++nr){
				struc.symm(i)[oR+nr]+=nnh_[II].basisR(JJ).symm()[nr];
			}
			//loop over all unique triplets
			for(int k=j+1; k<nlist.size(i); ++k){
				//find the index of the types of atom k
				const int kk=nlist.index(i,k);
				const int KK=struc.type(kk);
				//get the distance from K to I
				const Eigen::Vector3d rIK=struc.posn(i)-(struc.posn(kk)+struc.R()*nlist.img(i,k));
				const double drIK=rIK.norm();
				//compute the cosIJK angle and store the distances
				const double cosIJK=rIJ.dot(rIK)/(drIJ*drIK);
				//compute the IJ,IK,JK contribution to all angular basis functions
				nnh_[II].basisA(JJ,KK).symm(drIJ,drIK,cosIJK);
				const int oA=nnh_[II].nInputR()+nnh_[II].offsetA(JJ,KK);
				const int nA=nnh_[II].basisA(JJ,KK).size();
				for(int na=0; na<nA; ++na){
					struc.symm(i)[oA+na]+=nnh_[II].basisA(JJ,KK).symm()[na];
				}
			}
		}
	}
}

void NNP::compute_symm(Structure& struc){
	NeighborList nlist(rc_);
	nlist.build(struc);
	compute_symm(struc,nlist);
}

/**
* Compute the forces on the atoms for a given structure
* @param nnp - the neural network potential
* @param struc - the structure which we will compute
* @param nlist - neighbor list for each atom (includes periodic images)
* The forces are computed by looping over all nearest-neighbor pairs
* and all unique nearest-neighbor triples.  Thus, the neighbor list must be set 
* for the structure, and the type must correspond to the index of the atomic types
* in the NNP.  In addition, the index for each neighbor must be set to -1 if it is
* a periodic image or within [0,natoms] in which case the Newton's third law force
* pair must be added to the neighbor atom.
*/
double NNP::compute_force(Structure& struc, const NeighborList& nlist){
	if(NNP_PRINT_FUNC>0) std::cout<<"NNP::compute_force(Structure&,const NeighborList&):\n";
	const int nAtoms=struc.nAtoms();
	//reset the energy
	struc.pe()=0;
	//reset the force
	for(int i=0; i<nAtoms; ++i) struc.force(i).setZero();
	//compute the forces
	for(int i=0; i<nAtoms; ++i){
		//get the index of types i
		const int II=struc.type(i);
		//execute the appropriate network
		nnh(II).nn().fpbp(struc.symm(i));
		//add to the energy 
		struc.pe()+=nnh(II).nn().out()[0];
		//calculate the network gradient
		nnh(II).dOdZ().grad(nnh(II).nn());
		//set the gradient (n.b. dodi - do/di - deriv. of out w.r.t. in)
		const Eigen::VectorXd& dEdG=nnh(II).dOdZ().dodi().row(0);
		//loop over all neighbors
		for(int j=0; j<nlist.size(i); ++j){
			//get the indices of the jth neighbor
			const int jj=nlist.index(i,j);
			const int JJ=struc.type(jj);
			//get the distance from J to I
			const Eigen::Vector3d rIJ=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
			const double drIJ=rIJ.norm();
			const double dIJi=1.0/drIJ;
			//compute the IJ contribution to the radial force
			const int offsetR_=nnh(II).offsetR(JJ);
			const double amp=nnh(II).basisR(JJ).force(drIJ,dEdG.data()+offsetR_)*dIJi;
			struc.force(i).noalias()+=amp*rIJ;
			//if(nlist.img(i,j).squaredNorm()==0) struc.force(jj).noalias()-=amp*rIJ;
			struc.force(jj).noalias()-=amp*rIJ;
			//loop over all unique triplets
			for(int k=j+1; k<nlist.size(i); ++k){
				//find the index of the types of atom k
				const int kk=nlist.index(i,k);
				const int KK=struc.type(kk);
				//get the distance from K to I
				const Eigen::Vector3d rIK=struc.posn(i)-(struc.posn(kk)+struc.R()*nlist.img(i,k));
				const double drIK=rIK.norm();
				const double dIKi=1.0/drIK;
				//compute the IJ,IK,JK contribution to the angular force
				const int offsetA_=nnh(II).nInputR()+nnh(II).offsetA(JJ,KK);
				double phi=0; double eta[2]={0,0};
				const double cosIJK=rIJ.dot(rIK)*dIJi*dIKi;
				nnh(II).basisA(JJ,KK).force(drIJ,drIK,cosIJK,phi,eta,dEdG.data()+offsetA_);
				struc.force(i).noalias()+=(phi*(dIKi-cosIJK*dIJi)+eta[0])*rIJ*dIJi;
				struc.force(i).noalias()+=(phi*(dIJi-cosIJK*dIKi)+eta[1])*rIK*dIKi;
				/*if(nlist.img(i,j).sum()==0){
					struc.force(jj).noalias()-=(-phi*cosIJK*dIJi+eta[0])*rIJ*dIJi+phi*dIJi*rIK*dIKi;
					//struc.force(jj).noalias()-=eta[2]*rJK*dJKi;
				}*/
				struc.force(jj).noalias()-=(-phi*cosIJK*dIJi+eta[0])*rIJ*dIJi+phi*dIJi*rIK*dIKi;
				/*if(nlist.img(i,k).sum()==0){
					struc.force(kk).noalias()-=(-phi*cosIJK*dIKi+eta[1])*rIK*dIKi+phi*dIKi*rIJ*dIJi;
					//struc.force(kk).noalias()+=eta[2]*rJK*dJKi;
				}*/
				struc.force(kk).noalias()-=(-phi*cosIJK*dIKi+eta[1])*rIK*dIKi+phi*dIKi*rIJ*dIJi;
			}
		}
	}
	return struc.pe();
}

double NNP::compute_force(Structure& struc){
	NeighborList nlist(rc_);
	nlist.build(struc);
	return compute_force(struc,nlist);
}

double NNP::compute(Structure& struc, const NeighborList& nlist){
	compute_symm(struc,nlist);
	return compute_force(struc,nlist);
}

double NNP::compute(Structure& struc){
	NeighborList nlist(rc_);
	nlist.build(struc);
	compute_symm(struc,nlist);
	return compute_force(struc,nlist);
}

//************************************************************
// serialization
//************************************************************

namespace serialize{
	
//**********************************************
// byte measures
//**********************************************

template <> int nbytes(const NNP& obj){
	if(NNP_PRINT_FUNC>0) std::cout<<"nbytes(const NNP&):\n";
	int size=0;
	//types
	size+=nbytes(obj.ntypes());
	size+=nbytes(obj.map());
	for(int i=0; i<obj.ntypes(); ++i){
		size+=nbytes(obj.nnh(i));
	}
	//cutoff
	size+=nbytes(obj.rc());
	//return the size
	return size;
}

//**********************************************
// packing
//**********************************************

template <> int pack(const NNP& obj, char* arr){
	if(NNP_PRINT_FUNC>0) std::cout<<"pack(const NNP&,char*):\n";
	int pos=0;
	//types
	pos+=pack(obj.ntypes(),arr+pos);
	pos+=pack(obj.map(),arr+pos);
	for(int i=0; i<obj.ntypes(); ++i){
		pos+=pack(obj.nnh(i),arr+pos);
	}
	//cutoff
	pos+=pack(obj.rc(),arr+pos);
	//return bytes written
	return pos;
}

//**********************************************
// unpacking
//**********************************************

template <> int unpack(NNP& obj, const char* arr){
	if(NNP_PRINT_FUNC>0) std::cout<<"unpack(NNP&,const char*):\n";
	int pos=0;
	//types
	int ntypes=0;
	Map<int,int> map;
	pos+=unpack(ntypes,arr+pos);
	pos+=unpack(map,arr+pos);
	std::vector<NNH> nnh(ntypes);
	std::vector<NNH::Type> types(ntypes);
	for(int i=0; i<ntypes; ++i){
		pos+=unpack(nnh[i],arr+pos);
		types[i]=nnh[i].type();
	}
	obj.resize(types);
	for(int i=0; i<obj.ntypes(); ++i){
		obj.nnh(i)=nnh[i];
	}
	obj.map()=map;
	//cutoff
	pos+=unpack(obj.rc(),arr+pos);
	//return bytes read
	return pos;
}

}
