//c++ libraries
#include <iostream>
//c libraries
#include <stdexcept>
// ann - strings
#include "str/string.hpp"
// ann - chemistry
#include "chem/ptable.hpp"
// ann - eigen
#include "math/const.hpp"
#include "math/eigen.hpp"
// ann - print
#include "str/print.hpp"
// ann - structure
#include "struc/structure.hpp"

//**********************************************************************************************
//AtomData
//**********************************************************************************************

//==== operators ====

std::ostream& operator<<(std::ostream& out, const AtomData& obj){
	out<<"natoms = "<<obj.nAtoms_<<"\n";
	out<<"type   = "<<obj.atom_;
	return out;
}

//==== member functions ====

void AtomData::clear(){
	if(STRUC_PRINT_FUNC>0) std::cout<<"AtomData::clear():\n";
	//basic properties
	name_.clear();
	an_.clear();
	type_.clear();
	//serial properties
	mass_.clear();
	charge_.clear();
	radius_.clear();
	//vector properties
	image_.clear();
	posn_.clear();
	vel_.clear();
	force_.clear();
	dipole_.clear();
	//symm
	symm_.clear();
}

//==== resizing ====

void AtomData::resize(int nAtoms, const Atom& atom){
	if(STRUC_PRINT_FUNC>0) std::cout<<"AtomData::resize(int,const Atom&):\n";
	//check arguments
	if(nAtoms<0) throw std::runtime_error("AtomData::resize(int,const Atom&): invalid number of atoms");
	//set atom info
	atom_=atom;
	nAtoms_=nAtoms;
	if(nAtoms_>0){
		//basic properties
		if(atom.name())   name_.resize(nAtoms_);
		if(atom.an())     an_.resize(nAtoms_,0);
		if(atom.type())   type_.resize(nAtoms_,-1);
		//serial properties
		if(atom.mass())   mass_.resize(nAtoms_,0.0);
		if(atom.charge()) charge_.resize(nAtoms_,0.0);
		if(atom.radius()) radius_.resize(nAtoms_,0.0);
		if(atom.eta())    eta_.resize(nAtoms_,0.0);
		//vector properties
		if(atom.image())  image_.resize(nAtoms_,Eigen::Vector3i::Zero());
		if(atom.posn())   posn_.resize(nAtoms_,Eigen::Vector3d::Zero());
		if(atom.vel())    vel_.resize(nAtoms_,Eigen::Vector3d::Zero());
		if(atom.force())  force_.resize(nAtoms_,Eigen::Vector3d::Zero());
		if(atom.dipole()) dipole_.resize(nAtoms_,Eigen::Vector3d::Zero());
		//nnp
		if(atom.symm())	symm_.resize(nAtoms_,Eigen::Vector3d::Zero());
	}
}

//**********************************************************************************************
//Structure
//**********************************************************************************************

//==== operators ====

std::ostream& operator<<(std::ostream& out, const Structure& struc){
	char* str=new char[print::len_buf];
	out<<print::buf(str)<<"\n";
	out<<print::title("STRUCTURE",str)<<"\n";
	out<<static_cast<const AtomData&>(struc)<<"\n";
	out<<static_cast<const Cell&>(struc)<<"\n";
	out<<static_cast<const State&>(struc)<<"\n";
	out<<print::buf(str);
	delete[] str;
	return out;
}

//==== member functions ====

void Structure::clear(){
	if(STRUC_PRINT_FUNC>0) std::cout<<"Structure::clear():\n";
	AtomData::clear();
	Cell::clear();
	State::clear();
}

Structure& Structure::super(const Structure& struc, Structure& superc, const Eigen::Vector3i nlat){
	if(nlat[0]<=0 || nlat[1]<=0 || nlat[2]<=0) throw std::invalid_argument("Invalid lattice.");
	const int np=nlat.prod();
	const int nAtomsT=struc.nAtoms()*np;
	superc.resize(nAtomsT,struc.atom());
	//set the atomic properties
	int c=0;
	const Atom& atom=struc.atom();
	for(int i=0; i<nlat[0]; ++i){
		for(int j=0; j<nlat[1]; ++j){
			for(int k=0; k<nlat[2]; ++k){
				const Eigen::Vector3d R=i*struc.R().col(0)+j*struc.R().col(1)+k*struc.R().col(2);
				for(int n=0; n<struc.nAtoms(); ++n){
					//basic properties
					if(atom.name())  superc.name(c)=struc.name(n);
					if(atom.an())    superc.an(c)=struc.an(n);
					if(atom.type())  superc.type(c)=struc.type(n);
					//serial properties
					if(atom.mass())  superc.mass(c)=struc.mass(n);
					if(atom.charge())superc.charge(c)=struc.charge(n);
					if(atom.radius())superc.radius(c)=struc.radius(n);
					if(atom.eta())   superc.eta(c)=struc.eta(n);
					//vector properties
					if(atom.image()) superc.image(c)=struc.image(n);
					if(atom.posn())  superc.posn(c)=struc.posn(n)+R;
					if(atom.vel())   superc.vel(c)=struc.vel(n);
					if(atom.force()) superc.force(c)=struc.force(n);
					if(atom.dipole())superc.dipole(c)=struc.dipole(n);
					//nnp
					if(atom.symm())  superc.symm(c)=struc.symm(n);
					//increment
					c++;
				}
			}
		}
	}
	Eigen::MatrixXd Rnew=struc.R();
	Rnew.col(0)*=nlat[0];
	Rnew.col(1)*=nlat[1];
	Rnew.col(2)*=nlat[2];
	static_cast<Cell&>(superc).init(Rnew);
	return superc;
}

//==== static functions ====

void Structure::write_binary(const Structure& struc, const char* file){
	if(STRUC_PRINT_FUNC>0) std::cout<<"Structure::write_binary(const char*):\n";
	//local variables
	char* arr=NULL;
	FILE* writer=NULL;
	bool error=false;
	int nWrite=-1;
	try{
		//open file
		writer=fopen(file,"wb");
		if(writer==NULL) throw std::runtime_error(std::string("write_binary(Structure&,const char*): Could not open file: ")+std::string(file));
		//allocate buffer
		const int nBytes=serialize::nbytes(struc);
		arr=new char[nBytes];
		if(arr==NULL) throw std::runtime_error("write_binary(Structure&,const char*): Could not allocate memory.");
		//write to buffer
		serialize::pack(struc,arr);
		//write to file
		nWrite=fwrite(&nBytes,sizeof(int),1,writer);
		if(nWrite!=1) throw std::runtime_error("write_binary(Structure&,const char*): Write error.");
		nWrite=fwrite(arr,sizeof(char),nBytes,writer);
		if(nWrite!=nBytes) throw std::runtime_error("write_binary(Structure&,const char*): Write error.");
		//close the file, free memory
		delete[] arr; arr=NULL;
		fclose(writer); writer=NULL;
	}catch(std::exception& e){
		std::cout<<"ERROR in write_binary(Structure& struc,const char*):\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	//free local variables
	if(arr!=NULL) delete[] arr;
	if(writer!=NULL) fclose(writer);
	if(error) throw std::runtime_error("Failed to write");
}

void Structure::read_binary(Structure& struc, const char* file){
	if(STRUC_PRINT_FUNC>0) std::cout<<"Structure::read_binary(const char*):\n";
	//local variables
	char* arr=NULL;
	FILE* reader=NULL;
	bool error=false;
	int nRead=-1;
	try{
		//open file
		reader=fopen(file,"rb");
		if(reader==NULL) throw std::runtime_error(std::string("read_binary(Structure&,const char*): Could not open file: ")+std::string(file));
		//find size
		int nBytes=0;
		nRead=fread(&nBytes,sizeof(int),1,reader);
		if(nRead!=1) throw std::runtime_error("read_binary(Structure&,const char*): Read error.");
		//allocate buffer
		arr=new char[nBytes];
		if(arr==NULL) throw std::runtime_error("read_binary(Structure&,const char*): Could not allocate memory.");
		//read from file
		nRead=fread(arr,sizeof(char),nBytes,reader);
		if(nRead!=nBytes) throw std::runtime_error("read_binary(Structure&,const char*): Read error.");
		//read from buffer
		serialize::unpack(struc,arr);
		//close the file, free memory
		delete[] arr; arr=NULL;
		fclose(reader); reader=NULL;
	}catch(std::exception& e){
		std::cout<<"ERROR in read_binary(Structure& struc,const char*):\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	//free local variables
	if(arr!=NULL) delete[] arr;
	if(reader!=NULL) fclose(reader);
	if(error) throw std::runtime_error("Failed to read");
}

//**********************************************************************************************
// serialization
//**********************************************************************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const AtomData& obj){
		if(STRUC_PRINT_FUNC>0) std::cout<<"nbytes(const AtomData&)\n";
		int size=0;
		//atom type
		size+=nbytes(obj.atom());
		//number of atoms
		size+=sizeof(obj.nAtoms());
		//basic properties
		if(obj.atom().name())   size+=nbytes(obj.name());
		if(obj.atom().an())     size+=nbytes(obj.an());
		if(obj.atom().type())   size+=nbytes(obj.type());
		//serial properties
		if(obj.atom().mass())   size+=nbytes(obj.mass());
		if(obj.atom().charge()) size+=nbytes(obj.charge());
		if(obj.atom().radius()) size+=nbytes(obj.radius());
		if(obj.atom().eta())    size+=nbytes(obj.eta());
		//vector properties
		if(obj.atom().image())  size+=nbytes(obj.image());
		if(obj.atom().posn())   size+=nbytes(obj.posn());
		if(obj.atom().vel())    size+=nbytes(obj.vel());
		if(obj.atom().force())  size+=nbytes(obj.force());
		if(obj.atom().dipole()) size+=nbytes(obj.dipole());
		//nnp
		if(obj.atom().symm())   size+=nbytes(obj.symm());
		//return
		return size;
	}
	template <> int nbytes(const Structure& obj){
		if(STRUC_PRINT_FUNC>0) std::cout<<"nbytes(const Structure&)\n";
		int size=0;
		size+=nbytes(static_cast<const Cell&>(obj));
		size+=nbytes(static_cast<const State&>(obj));
		size+=nbytes(static_cast<const AtomData&>(obj));
		return size;
	}
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const AtomData& obj, char* arr){
		if(STRUC_PRINT_FUNC>0) std::cout<<"pack(const AtomData&,char*):\n";
		int pos=0;
		//atom type
		pos+=pack(obj.atom(),arr+pos);
		//natoms
		std::memcpy(arr+pos,&obj.nAtoms(),sizeof(int)); pos+=sizeof(int);
		//basic properties
		if(obj.atom().name())   pos+=pack(obj.name(),arr+pos);
		if(obj.atom().an())     pos+=pack(obj.an(),arr+pos);
		if(obj.atom().type())   pos+=pack(obj.type(),arr+pos);
		//serial properties
		if(obj.atom().mass())   pos+=pack(obj.mass(),arr+pos);
		if(obj.atom().charge()) pos+=pack(obj.charge(),arr+pos);
		if(obj.atom().radius()) pos+=pack(obj.radius(),arr+pos);
		if(obj.atom().eta())    pos+=pack(obj.eta(),arr+pos);
		//vector properties
		if(obj.atom().image())  pos+=pack(obj.image(),arr+pos);
		if(obj.atom().posn())   pos+=pack(obj.posn(),arr+pos);
		if(obj.atom().vel())    pos+=pack(obj.vel(),arr+pos);
		if(obj.atom().force())  pos+=pack(obj.force(),arr+pos);
		if(obj.atom().dipole()) pos+=pack(obj.dipole(),arr+pos);
		//nnp
		if(obj.atom().symm())   pos+=pack(obj.symm(),arr+pos);
		//return
		return pos;
	}
	template <> int pack(const Structure& obj, char* arr){
		if(STRUC_PRINT_FUNC>0) std::cout<<"pack(const Structure&,char*):\n";
		int pos=0;
		pos+=pack(static_cast<const Cell&>(obj),arr+pos);
		pos+=pack(static_cast<const State&>(obj),arr+pos);
		pos+=pack(static_cast<const AtomData&>(obj),arr+pos);
		return pos;
	}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(AtomData& obj, const char* arr){
		if(STRUC_PRINT_FUNC>0) std::cout<<"unpack(AtomData&,const char*):\n";
		int pos=0;
		//atom type
		Atom atom;
		pos+=unpack(atom,arr+pos);
		//natoms
		int nAtoms=0;
		std::memcpy(&nAtoms,arr+pos,sizeof(int)); pos+=sizeof(int);
		//resize
		obj.resize(nAtoms,atom);
		//basic properties
		if(obj.atom().name())   pos+=unpack(obj.name(),arr+pos);
		if(obj.atom().an())     pos+=unpack(obj.an(),arr+pos);
		if(obj.atom().type())   pos+=unpack(obj.type(),arr+pos);
		//serial properties
		if(obj.atom().mass())   pos+=unpack(obj.mass(),arr+pos);
		if(obj.atom().charge()) pos+=unpack(obj.charge(),arr+pos);
		if(obj.atom().radius()) pos+=unpack(obj.radius(),arr+pos);
		if(obj.atom().eta())    pos+=unpack(obj.eta(),arr+pos);
		//vector properties
		if(obj.atom().image())  pos+=unpack(obj.image(),arr+pos);
		if(obj.atom().posn())   pos+=unpack(obj.posn(),arr+pos);
		if(obj.atom().vel())    pos+=unpack(obj.vel(),arr+pos);
		if(obj.atom().force())  pos+=unpack(obj.force(),arr+pos);
		if(obj.atom().dipole()) pos+=unpack(obj.dipole(),arr+pos);
		//nnp
		if(obj.atom().symm())   pos+=unpack(obj.symm(),arr+pos);
		//return
		return pos;
	}
	template <> int unpack(Structure& obj, const char* arr){
		if(STRUC_PRINT_FUNC>0) std::cout<<"unpack(Structure&,const char*):\n";
		int pos=0;
		pos+=unpack(static_cast<Cell&>(obj),arr+pos);
		pos+=unpack(static_cast<State&>(obj),arr+pos);
		pos+=unpack(static_cast<AtomData&>(obj),arr+pos);
		return pos;
	}
	
}
