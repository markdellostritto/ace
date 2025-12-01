// c libraries
#include <cstdio>
#include <ctime>
// c++ libraries
#include <string>
#include <stdexcept>
#include <iostream>
// ann - structure
#include "struc/structure.hpp"
// ann - strings
#include "str/string.hpp"
// ann - chem
#include "chem/units.hpp"
#include "chem/ptable.hpp"
// ann - format
#include "format/xyz_struc.hpp"

namespace XYZ{

//*****************************************************
//reading
//*****************************************************

void read(const char* xyzfile, const Atom& atom, Structure& struc){
	//open file
	if(XYZ_PRINT_STATUS>0) std::cout<<"opening file\n";
	FILE* reader=fopen(xyzfile,"r");
	if(reader==NULL) throw std::runtime_error(std::string("ERROR in XYZ::read(const char*,const Atom&,Structure&): Could not open file: ")+std::string(xyzfile));
	//read
	read(reader,atom,struc);
	//close the file
	if(XYZ_PRINT_STATUS>0) std::cout<<"closing file\n";
	fclose(reader);
	reader=NULL;
}

void read(FILE* reader, const Atom& atom, Structure& struc){
	const std::string funcame=std::string("XYZ::read(FILE*,const Atom&,Structure&)");
	if(XYZ_PRINT_FUNC>0) std::cout<<funcame<<":\n";
	//==== local function variables ====
	//file i/o
		char* input=new char[string::M];
		char* name=new char[string::M];
		Token token;
	//atom info
		int nAtoms=0;
		double pe=0;
		Eigen::Matrix3d lv=Eigen::Matrix3d::Zero();
	//units
		double s_len=0.0,s_energy=0.0,s_mass=0.0;
		if(units::Consts::system()==units::System::LJ){
			s_len=1.0;
			s_energy=1.0;
			s_mass=1.0;
		} else if(units::Consts::system()==units::System::AU){
			s_len=units::Ang2Bohr;
			s_energy=units::Ev2Eh;
			s_mass=units::MPoME;
		} else if(units::Consts::system()==units::System::METAL){
			s_len=1.0;
			s_energy=1.0;
			s_mass=1.0;
		} 
		else throw std::runtime_error("Invalid units.");
		
	//read natoms
	if(XYZ_PRINT_STATUS>0) std::cout<<"reading natoms\n";
	fgets(input,string::M,reader);
	nAtoms=std::atoi(input);
	if(nAtoms<=0) throw std::runtime_error(funcame+": found zero atoms.");
	
	//read header
	int ndata=0;
	int ni=-1; int qi=-1; int mi=-1;
	Eigen::Vector3i ri=Eigen::Vector3i::Constant(-1);
	Eigen::Vector3i fi=Eigen::Vector3i::Constant(-1);
	Eigen::Vector3i di=Eigen::Vector3i::Constant(-1);
	fgets(input,string::M,reader);
	string::to_upper(input);
	if(std::strstr(input,"PROPERTIES")!=NULL){
		token.read(std::strstr(input,"PROPERTIES")," \r\t\n="); token.next();
		const std::string propstr=token.next();
		Token proptok=Token(propstr,":");
		while(!proptok.end()){
			const std::string tag=proptok.next();
			if(tag=="SPECIES"){
				if(proptok.next()!="S") throw std::runtime_error(funcame+": invalid name data type.");
				else if(std::atoi(proptok.next().c_str())!=1) throw std::runtime_error(funcame+": invalid name length.");
				else ni=ndata++;
			} else if(tag=="CHARGE"){
				if(proptok.next()!="R") throw std::runtime_error(funcame+": invalid charge data type.");
				else if(std::atoi(proptok.next().c_str())!=1) throw std::runtime_error(funcame+": invalid charge length.");
				else qi=ndata++;
			} else if(tag=="MASS"){
				if(proptok.next()!="R") throw std::runtime_error(funcame+": invalid charge data type.");
				else if(std::atoi(proptok.next().c_str())!=1) throw std::runtime_error(funcame+": invalid charge length.");
				else mi=ndata++;
			} else if(tag=="POS" || tag=="POSITION"){
				if(proptok.next()!="R") throw std::runtime_error(funcame+": invalid position data type.");
				else if(std::atoi(proptok.next().c_str())!=3) throw std::runtime_error(funcame+": invalid position length.");
				else for(int i=0; i<3; ++i) ri[i]=ndata++;
			} else if(tag=="FORCES" || tag=="FORCE"){
				if(proptok.next()!="R") throw std::runtime_error(funcame+": invalid force data type.");
				else if(std::atoi(proptok.next().c_str())!=3) throw std::runtime_error(funcame+": invalid force length.");
				else for(int i=0; i<3; ++i) fi[i]=ndata++;
			} else if(tag=="DIPOLE" || tag=="DIPOLES"){
				if(proptok.next()!="R") throw std::runtime_error(funcame+": invalid dipole data type.");
				else if(std::atoi(proptok.next().c_str())!=3) throw std::runtime_error(funcame+": invalid dipole length.");
				else for(int i=0; i<3; ++i) di[i]=ndata++;
			} 
		}
	}
	if(std::strstr(input,"POTENTIAL_ENERGY")!=NULL){
		token.read(std::strstr(input,"POTENTIAL_ENERGY")," \r\t\n="); token.next();
		pe=std::atof(token.next().c_str());
	}
	if(std::strstr(input,"LATTICE")!=NULL){
		token.read(std::strstr(input,"LATTICE")," \r\t\n=\"");
		token.next();
		lv(0,0)=std::atof(token.next().c_str());
		lv(1,0)=std::atof(token.next().c_str());
		lv(2,0)=std::atof(token.next().c_str());
		lv(0,1)=std::atof(token.next().c_str());
		lv(1,1)=std::atof(token.next().c_str());
		lv(2,1)=std::atof(token.next().c_str());
		lv(0,2)=std::atof(token.next().c_str());
		lv(1,2)=std::atof(token.next().c_str());
		lv(2,2)=std::atof(token.next().c_str());
		lv*=s_len;
	}
	
	//resize the structure
	if(XYZ_PRINT_STATUS>0) std::cout<<"resizing structure\n";
	struc.resize(nAtoms,atom);
	
	//read iatoms
	if(XYZ_PRINT_STATUS>0) std::cout<<"reading names and posns\n";
	std::vector<std::string> sarr(ndata);
	for(int i=0; i<nAtoms; ++i){
		int c=0;
		token.read(fgets(input,string::M,reader),string::WS);
		while(!token.end()) sarr[c++]=token.next();
		if(struc.atom().name() && ni>=0){
			struc.name(i)=sarr[ni];
		} 
		if(struc.atom().charge() && qi>=0){
			struc.charge(i)=std::atof(sarr[qi].c_str());
		}
		if(struc.atom().mass() && mi>=0){
			struc.mass(i)=std::atof(sarr[mi].c_str());
		}
		if(struc.atom().posn() && ri.minCoeff()>=0){
			struc.posn(i)[0]=std::atof(sarr[ri[0]].c_str())*s_len;
			struc.posn(i)[1]=std::atof(sarr[ri[1]].c_str())*s_len;
			struc.posn(i)[2]=std::atof(sarr[ri[2]].c_str())*s_len;
		}
		if(struc.atom().force() && fi.minCoeff()>=0){
			struc.force(i)[0]=std::atof(sarr[fi[0]].c_str())*s_energy/s_len;
			struc.force(i)[1]=std::atof(sarr[fi[1]].c_str())*s_energy/s_len;
			struc.force(i)[2]=std::atof(sarr[fi[2]].c_str())*s_energy/s_len;
		}
		if(struc.atom().dipole() && di.minCoeff()>=0){
			struc.dipole(i)[0]=std::atof(sarr[di[0]].c_str())*s_len;
			struc.dipole(i)[1]=std::atof(sarr[di[1]].c_str())*s_len;
			struc.dipole(i)[2]=std::atof(sarr[di[2]].c_str())*s_len;
		}
	}

	//set the cell
	if(lv.norm()>0){
		if(XYZ_PRINT_STATUS>0) std::cout<<"setting cell\n";
		static_cast<Cell&>(struc).init(lv);
		for(int i=0; i<nAtoms; ++i){
			Cell::returnToCell(struc.posn(i),struc.posn(i),struc.R(),struc.RInv());
		}
	}
	
	//set the energy
	if(XYZ_PRINT_STATUS>0) std::cout<<"setting energy\n";
	struc.pe()=s_energy*pe;
	
	//set an
	if(atom.an() && atom.name()){
		for(int i=0; i<nAtoms; ++i){
			struc.an(i)=ptable::an(struc.name(i).c_str());
		}
	}
	
	//set mass
	if(mi<0){
		if(atom.an() && atom.mass()){
			for(int i=0; i<nAtoms; ++i){
				struc.mass(i)=ptable::mass(struc.an(i))*s_mass;
			}
		} else if(atom.name() && atom.mass()){
			for(int i=0; i<nAtoms; ++i){
				const int an=ptable::an(struc.name(i).c_str());
				struc.mass(i)=ptable::mass(an)*s_mass;
			}
		}
	}
	
	//set radius
	if(atom.an() && atom.radius()){
		for(int i=0; i<nAtoms; ++i){
			struc.radius(i)=ptable::radius_covalent(struc.an(i));
		}
	} else if(atom.name() && atom.radius()){
		for(int i=0; i<nAtoms; ++i){
			const int an=ptable::an(struc.name(i).c_str());
			struc.radius(i)=ptable::radius_covalent(an);
		}
	}

	//set the type
	if(atom.name() && atom.type()){
		std::vector<std::string> names;
		int ntypes=0;
		for(int i=0; i<nAtoms; ++i){
			bool match=false;
			for(int j=0; j<names.size(); ++j){
				if(struc.name(i)==names[j]){
					struc.type(i)=j;
					match=true; break;
				}
			}
			if(!match){
				struc.type(i)=ntypes++;
				names.push_back(struc.name(i));
			}
		}
	}
	
	//free memory
	delete[] input;
	delete[] name;
}

//*****************************************************
//writing
//*****************************************************

void write(const char* file, const Atom& atom, const Structure& struc){
	//open file
	if(XYZ_PRINT_STATUS>0) std::cout<<"opening file\n";
	FILE* writer=fopen(file,"w");
	if(writer==NULL) throw std::runtime_error("Runtime Error: Could not open file: \""+std::string(file)+"\"");
	//write xyz
	if(XYZ_PRINT_STATUS>0) std::cout<<"write xyz\n";
	write(writer,atom,struc);
	//close file
	if(XYZ_PRINT_STATUS>0) std::cout<<"closing file\n";
	fclose(writer);
	writer=NULL;
}
void write(FILE* writer, const Atom& atom, const Structure& struc){
	if(XYZ_PRINT_FUNC>0) std::cout<<"write(const char*,const Atom&,const Structure&):\n";
	fprintf(writer,"%i\n",struc.nAtoms());
	const Eigen::Matrix3d& R=struc.R();
	if(atom.force()){
		fprintf(writer,"Properties=species:S:1:pos:R:3:forces:R:3 potential_energy=%f pbc=\"T T T\" Lattice=\"%f %f %f %f %f %f %f %f %f\"\n",
			struc.pe(),R(0,0),R(1,0),R(2,0),R(0,1),R(1,1),R(1,2),R(0,2),R(1,2),R(2,2)
		);
		for(int i=0; i<struc.nAtoms(); ++i){
			fprintf(writer,"%-2s %19.10f %19.10f %19.10f %19.10f %19.10f %19.10f\n",struc.name(i).c_str(),
				struc.posn(i)[0],struc.posn(i)[1],struc.posn(i)[2],
				struc.force(i)[0],struc.force(i)[1],struc.force(i)[2]
			);
		}
	} else {
		fprintf(writer,"Properties=species:S:1:pos:R:3 potential_energy=%f pbc=\"T T T\" Lattice=\"%f %f %f %f %f %f %f %f %f\"\n",
			struc.pe(),R(0,0),R(1,0),R(2,0),R(0,1),R(1,1),R(1,2),R(0,2),R(1,2),R(2,2)
		);
		for(int i=0; i<struc.nAtoms(); ++i){
			fprintf(writer,"%-2s %19.10f %19.10f %19.10f\n",struc.name(i).c_str(),
				struc.posn(i)[0],struc.posn(i)[1],struc.posn(i)[2]
			);
		}
	}
}

	
}