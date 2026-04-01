// c libraries
#include <cstdio>
#include <ctime>
// c++ libraries
#include <string>
#include <stdexcept>
#include <iostream>
// str
#include "str/string.hpp"
// chem
#include "chem/units.hpp"
#include "chem/ptable.hpp"
// format
#include "format/xyz_traj.hpp"

namespace XYZ{

void read(const char* file, const Interval& interval, const Atom& atom, Trajectory& traj){
	if(XYZ_PRINT_FUNC>0) std::cout<<"XYZ::read(const char*,Interval&,const Atom&,Trajectory&):\n";
	//==== local function variables ====
	//file i/o
		FILE* reader=NULL;
		char* input=new char[string::M];
		char* name=new char[string::M];
		Token token;
	//atom info	
		int nAtoms=0;
		Eigen::Matrix3d lv=Eigen::Matrix3d::Zero();
	//positions
		Eigen::Vector3d r;
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
	
	//open file
	if(XYZ_PRINT_STATUS>0) std::cout<<"opening file\n";
	reader=fopen(file,"r");
	if(reader==NULL) throw std::runtime_error("Runtime Error: Could not open file.");
	
	//read natoms
	if(XYZ_PRINT_STATUS>0) std::cout<<"reading natoms\n";
	fgets(input,string::M,reader);
	nAtoms=std::atoi(input);
	if(nAtoms<=0) throw std::runtime_error("Runtime Error: found zero atoms.");
	if(XYZ_PRINT_DATA>0) std::cout<<"natoms = "<<nAtoms<<"\n";
	
	//find the total number of timesteps
	if(XYZ_PRINT_STATUS>0) std::cout<<"reading timesteps\n";
	std::rewind(reader);
	int nlines=0;
	while(fgets(input,string::M,reader)) ++nlines;
	int ts=nlines/(nAtoms+2);//natoms + natoms-line + comment-line
	if(XYZ_PRINT_DATA>0) std::cout<<"ts = "<<ts<<"\n";
	
	//set the interval
	if(XYZ_PRINT_STATUS>0) std::cout<<"setting interval\n";
	const int ibeg=Interval::index(interval.beg(),ts);
	const int iend=Interval::index(interval.end(),ts);
	const int len=iend-ibeg+1;
	const int nsteps=len/interval.stride();
	
	//resize the simulation
	if(XYZ_PRINT_STATUS>0) std::cout<<"resizing simulation\n";
	traj.resize(nsteps,nAtoms,atom);
	
	//read the simulation
	if(XYZ_PRINT_STATUS>0) std::cout<<"reading simulation\n";
	std::rewind(reader);
	for(int t=0; t<ibeg; ++t){
		fgets(input,string::M,reader);//natoms
		const int nAtoms_=std::atoi(input);
		fgets(input,string::M,reader);//comment line
		for(int n=0; n<nAtoms_; ++n){
			fgets(input,string::M,reader);
		}
	}
	for(int t=0; t<traj.timesteps(); ++t){
		//read natoms
		fgets(input,string::M,reader);//natoms
		const int nAtoms_=std::atoi(input);
		if(traj.frame(t).nAtoms()!=nAtoms_) throw std::invalid_argument("Error in XYZ::read(const char*,Interval&,const Atom&,Trajectory&): Invalid number of atoms.");
		//read header
		double pe=0;
		int ndata=0;
		int ni=-1; int qi=-1; int mi=-1;
		Eigen::Vector3i ri=Eigen::Vector3i::Constant(-1);
		Eigen::Vector3i fi=Eigen::Vector3i::Constant(-1);
		fgets(input,string::M,reader);
		string::to_upper(input);
		if(std::strstr(input,"PROPERTIES")!=NULL){
			token.read(std::strstr(input,"PROPERTIES")," \r\t\n=");
			const std::string propstr=token.next();
			Token proptok=Token(propstr,":");
			while(!proptok.end()){
				const std::string tag=proptok.next();
				if(tag=="SPECIES"){
					if(token.next()!="S") throw std::runtime_error("ERROR in XYZ::read(const char*,const Atom&,Structure&): invalid name data type.");
					else if(std::atoi(token.next().c_str())!=1) throw std::runtime_error("ERROR in XYZ::read(const char*,const Atom&,Structure&): invalid name length.");
					else ni=ndata++;
				} else if(tag=="POS" || tag=="POSN" || tag=="POSITION"){
					if(token.next()!="R") throw std::runtime_error("ERROR in XYZ::read(const char*,const Atom&,Structure&): invalid position data type.");
					else if(std::atoi(token.next().c_str())!=3) throw std::runtime_error("ERROR in XYZ::read(const char*,const Atom&,Structure&): invalid position length.");
					else for(int i=0; i<3; ++i) ri[i]=ndata++;
				} else if(tag=="FORCES" || tag=="FORCE"){
					if(token.next()!="R") throw std::runtime_error("ERROR in XYZ::read(const char*,const Atom&,Structure&): invalid force data type.");
					else if(std::atoi(token.next().c_str())!=3) throw std::runtime_error("ERROR in XYZ::read(const char*,const Atom&,Structure&): invalid force length.");
					else for(int i=0; i<3; ++i) fi[i]=ndata++;
				} else if(tag=="CHARGE" || tag=="CHG" || tag=="Q"){
					if(proptok.next()!="R") throw std::runtime_error(funcame+": invalid charge data type.");
					else if(std::atoi(proptok.next().c_str())!=1) throw std::runtime_error(funcame+": invalid charge length.");
					else qi=ndata++;
				} else if(tag=="MASS"){
					if(proptok.next()!="R") throw std::runtime_error(funcame+": invalid charge data type.");
					else if(std::atoi(proptok.next().c_str())!=1) throw std::runtime_error(funcame+": invalid charge length.");
					else mi=ndata++;
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
		traj.frame(t).pe()=pe*s_energy;
		//read atom data
		std::vector<std::string> sarr(ndata);
		for(int n=0; n<nAtoms_; ++n){
			/*fgets(input,string::M,reader);
			std::sscanf(input,"%s %lf %lf %lf",name,&r[0],&r[1],&r[2]);
			traj.frame(t).posn(n).noalias()=s_len*r;
			if(atom.name()) traj.frame(t).name(n)=name;*/
			token.read(fgets(input,string::M,reader),string::WS);
			while(!token.end()) sarr[c++]=token.next();
			if(atom.name() && ni>=0){
				traj.frame(t).name(i)=sarr[ni];
			} 
			if(atom.charge() && qi>=0){
				traj.frame(t).charge(i)=std::atof(sarr[qi].c_str());
			}
			if(atom.mass() && mi>=0){
				traj.frame(t).mass(i)=std::atof(sarr[mi].c_str());
			}
			if(atom.posn() && ri.minCoeff()>=0){
				traj.frame(t).posn(i)[0]=std::atof(sarr[ri[0]].c_str())*s_len;
				traj.frame(t).posn(i)[1]=std::atof(sarr[ri[1]].c_str())*s_len;
				traj.frame(t).posn(i)[2]=std::atof(sarr[ri[2]].c_str())*s_len;
			}
			if(atom.force() && fi.minCoeff()>=0){
				traj.frame(t).force(i)[0]=std::atof(sarr[fi[0]].c_str())*s_energy/s_len;
				traj.frame(t).force(i)[1]=std::atof(sarr[fi[1]].c_str())*s_energy/s_len;
				traj.frame(t).force(i)[2]=std::atof(sarr[fi[2]].c_str())*s_energy/s_len;
			}
		}
		//set the cell
		if(XYZ_PRINT_STATUS>0) std::cout<<"setting cell\n";
		if(lv.norm()>0) static_cast<Cell&>(traj.frame(t)).init(lv);
		//skip "stride-1" steps
		for(int tt=0; tt<interval.stride()-1; ++tt){
			fgets(input,string::M,reader);//natoms
			const int NN=std::atoi(input);
			fgets(input,string::M,reader);//comment line
			for(int n=0; n<NN; ++n) fgets(input,string::M,reader);
		}
	}
	
	//set an
	if(atom.name() && atom.an()){
		for(int t=0; t<traj.timesteps(); ++t){
			Structure& struc=traj.frame(t);
			for(int i=0; i<nAtoms; ++i){
				struc.an(i)=ptable::an(struc.name(i).c_str());
			}
		}
	}
	
	//set mass
	if(atom.mass()){
		if(atom.an()){
			for(int t=0; t<traj.timesteps(); ++t){
				Structure& struc=traj.frame(t);
				for(int i=0; i<nAtoms; ++i){
					struc.mass(i)=ptable::mass(struc.an(i))*s_mass;
				}
			}
		} else if(atom.name()){
			for(int t=0; t<traj.timesteps(); ++t){
				Structure& struc=traj.frame(t);
				for(int i=0; i<nAtoms; ++i){
					const int an=ptable::an(struc.name(i).c_str());
					struc.mass(i)=ptable::mass(an)*s_mass;
				}
			}
		}
	}
	
	//close file
	if(XYZ_PRINT_STATUS>0) std::cout<<"closing file\n";
	fclose(reader);
	reader=NULL;
	
	//free memory
	delete[] input;
	delete[] name;
}

void write(const char* file, const Interval& interval, const Atom& atom, const Trajectory& traj){
	if(XYZ_PRINT_FUNC>0) std::cout<<"write(const char*,Interval&,const Atom&,Trajectory&):\n";
	FILE* writer=NULL;
	
	//open file
	if(XYZ_PRINT_STATUS>0) std::cout<<"opening file\n";
	writer=fopen(file,"w");
	if(writer==NULL) throw std::runtime_error("Runtime Error: Could not open file.");
	
	//check timing info
	if(XYZ_PRINT_STATUS>0) std::cout<<"setting interval\n";
	const int ibeg=Interval::index(interval.beg(),traj.timesteps());
	const int iend=Interval::index(interval.end(),traj.timesteps());
	
	//write simulation
	if(XYZ_PRINT_STATUS>0) std::cout<<"writing simulation\n";
	if(atom.force()){
		for(int t=ibeg; t<=iend; ++t){
			Structure struc=traj.frame(t);
			fprintf(writer,"%i\n",struc.nAtoms());
			if(struc.R().norm()>0){
				const Eigen::Matrix3d& R=struc.R();
				fprintf(writer,"Properties=species:S:1:pos:R:3:forces:R:3 potential_energy=%f pbc=\"T T T\" Lattice=\"%f %f %f %f %f %f %f %f %f\"\n",
					struc.pe(),R(0,0),R(1,0),R(2,0),R(0,1),R(1,1),R(1,2),R(0,2),R(1,2),R(2,2)
				);
			} else {
				fprintf(writer,"%s\n",traj.name().c_str());
			}
			for(int i=0; i<struc.nAtoms(); ++i){
				fprintf(writer,"  %-2s %19.10f %19.10f %19.10f %19.10f %19.10f %19.10f\n",struc.name(i).c_str(),
					struc.posn(i)[0],struc.posn(i)[1],struc.posn(i)[2],
					struc.force(i)[0],struc.force(i)[1],struc.force(i)[2]
				);
			}
		}
	} else {
		for(int t=ibeg; t<=iend; ++t){
			Structure struc=traj.frame(t);
			fprintf(writer,"%i\n",struc.nAtoms());
			if(struc.R().norm()>0){
				const Eigen::Matrix3d& R=struc.R();
				fprintf(writer,"Properties=species:S:1:pos:R:3 potential_energy=%f pbc=\"T T T\" Lattice=\"%f %f %f %f %f %f %f %f %f\"\n",
					struc.pe(),R(0,0),R(1,0),R(2,0),R(0,1),R(1,1),R(1,2),R(0,2),R(1,2),R(2,2)
				);
			} else {
				fprintf(writer,"%s\n",traj.name().c_str());
			}
			for(int i=0; i<struc.nAtoms(); ++i){
				fprintf(writer,"  %-2s %19.10f %19.10f %19.10f\n",struc.name(i).c_str(),
					struc.posn(i)[0],struc.posn(i)[1],struc.posn(i)[2]
				);
			}
		}
	}
	
	//close file
	if(XYZ_PRINT_STATUS>0) std::cout<<"closing file\n";
	fclose(writer);
	writer=NULL;
}
	
}