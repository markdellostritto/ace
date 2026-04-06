// c libraries
#include <cstdio>
#include <ctime>
// c++ libraries
#include <string>
#include <stdexcept>
#include <iostream>
// string
#include "str/string.hpp"
// chem
#include "chem/units.hpp"
// math
#include "math/const.hpp"
//struc
#include "struc/structure.hpp"
//cp2k
#include "format/cp2k_struc.hpp"

namespace CP2K{

//*****************************************************
//reading
//*****************************************************

void read(const char* file, const Atom& atom, Structure& struc){
	//file i/o
		FILE* reader=NULL;
		char* input=new char[string::M];
		//std::vector<std::string> strlist;
		Token token;
	//units
		double s_len=0.0,s_posn=0.0,s_energy=0.0,s_mass=0.0;
		if(units::Consts::system()==units::System::LJ){
			s_len=1.0;
			s_posn=1.0;
			s_energy=1.0;
			s_mass=1.0;
		} else if(units::Consts::system()==units::System::AU){
			s_len=1.0;
			s_posn=units::Ang2Bohr;
			s_energy=1.0;
			s_mass=units::MPoME;
		} else if(units::Consts::system()==units::System::METAL){
			s_len=units::Bohr2Ang;
			s_posn=1.0;
			s_energy=units::Eh2Ev;
			s_mass=1.0;
		}
		else throw std::runtime_error("Invalid units.");
	//structure
		int natomst=0;
		std::vector<int> natoms;
		int nspecies=0;
		std::vector<std::string> species;
	//flags
		const char* flag_atom="ATOMIC COORDINATES";
		const char* flag_forces="ATOMIC FORCES in";
		const char* flag_cell="CELL|";
		const char* flag_energy="ENERGY|";
		const char* flag_natoms="Atoms:";
	//misc
		bool error=false;
		
	try{
		//open the file
		if(CP2K_PRINT_STATUS>0) std::cout<<"opening the file: "<<file<<"\n";
		reader=fopen(file,"r");
		if(reader==NULL) throw std::runtime_error(std::string("ERROR: Could not open file: \"")+std::string(file)+std::string("\"\n"));
		
		//read data
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strstr(input,flag_natoms)!=NULL){
				natomst=std::atoi(std::strstr(input,flag_natoms)+6);
				if(natomst==0) throw std::runtime_error("Found zero atoms");
				struc.resize(natomst,atom);
			} else if(std::strstr(input,flag_atom)!=NULL){
				//skip three lines
				fgets(input,string::M,reader);
				fgets(input,string::M,reader);
				fgets(input,string::M,reader);
				for(int i=0; i<natomst; ++i){
					std::vector<std::string> tokens;
					token.read(fgets(input,string::M,reader),string::WS);
					while(!token.end()) tokens.push_back(token.next());
					if(atom.name()){
						struc.name(i)=tokens[2];
					}
					if(atom.posn()){
						struc.posn(i)[0]=std::atof(tokens[4].c_str());
						struc.posn(i)[1]=std::atof(tokens[5].c_str());
						struc.posn(i)[2]=std::atof(tokens[6].c_str());
					}
					if(atom.mass()){
						struc.mass(i)=std::atof(tokens[8].c_str())*s_mass;
					}
				}
			} else if(std::strstr(input,flag_forces)!=NULL){
				//skip two lines
				fgets(input,string::M,reader);
				fgets(input,string::M,reader);
				for(int i=0; i<natomst; ++i){
					token.read(fgets(input,string::M,reader),string::WS).next(3);
					struc.force(i)[0]=std::atof(token.next().c_str())*s_energy/s_len;
					struc.force(i)[1]=std::atof(token.next().c_str())*s_energy/s_len;
					struc.force(i)[2]=std::atof(token.next().c_str())*s_energy/s_len;
				}
			} else if(std::strstr(input,flag_energy)!=NULL){
				token.read(input,string::WS);
				while(!token.end()) token.next();
				struc.pe()=std::atof(token.token().c_str())*s_energy;
			} else if(std::strstr(input,flag_cell)!=NULL){
				Eigen::Matrix3d R;
				//vector - a
				token.read(fgets(input,string::M,reader),string::WS).next(4);
				R(0,0)=std::atof(token.next().c_str());
				R(1,0)=std::atof(token.next().c_str());
				R(2,0)=std::atof(token.next().c_str());
				//vector - b
				token.read(fgets(input,string::M,reader),string::WS).next(4);
				R(0,1)=std::atof(token.next().c_str());
				R(1,1)=std::atof(token.next().c_str());
				R(2,1)=std::atof(token.next().c_str());
				//vector - c
				token.read(fgets(input,string::M,reader),string::WS).next(4);
				R(0,2)=std::atof(token.next().c_str());
				R(1,2)=std::atof(token.next().c_str());
				R(2,2)=std::atof(token.next().c_str());
				static_cast<Cell&>(struc).init(R);
				//skip further "CELL|" lines
				fgets(input,string::M,reader);
				fgets(input,string::M,reader);
				fgets(input,string::M,reader);
				fgets(input,string::M,reader);
			} 
		}
	}catch(std::exception& e){
		std::cout<<"Error in CP2K::read()\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
}

}