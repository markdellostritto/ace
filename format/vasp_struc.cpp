// c libraries
#include <cstdio>
#include <ctime>
// c++ libraries
#include <iostream>
#include <string>
#include <stdexcept>
// eigen libraries
#include <Eigen/Dense>
// ann - structure
#include "struc/structure.hpp"
// ann - strings
#include "str/string.hpp"
// ann - chem
#include "chem/units.hpp"
#include "chem/ptable.hpp"
// ann - vasp
#include "format/vasp_struc.hpp"

namespace VASP{

//*****************************************************
//POSCAR
//*****************************************************

namespace POSCAR{

void read(const char* file, const Atom& atom, Structure& struc){
	const char* funcName="read(const char*,const Atom&,Structure&)";
	if(VASP_PRINT_FUNC>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<funcName<<":\n";
	//open file 
	FILE* reader=fopen(file,"r");
	if(reader==NULL) throw std::runtime_error(std::string("I/O Error: Could not open file: ")+file);
	//read vasp
	read(reader,atom,struc);
	//close file
	fclose(reader);
	reader=NULL;
}
void read(FILE* reader, const Atom& atom, Structure& struc){
	const char* funcName="read(FILE*,const Atom&,Structure&)";
	if(VASP_PRINT_FUNC>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<funcName<<":\n";
	
	//====  local function variables ==== 
	//file i/o
		char* input=new char[string::M];
		char* str_name=new char[string::M];
		char* str_number=new char[string::M];
		Token token;
	//simulation flags
		bool direct;//whether the coordinates are in direct or Cartesian coordinates
	//cell
		double scale=0.0;
		Cell cell;
		Eigen::Matrix3d lv;
	//atom info
		int nAtomsT=0;
		std::vector<int> nAtoms;//the number of atoms in each species
		std::vector<std::string> names;//the names of each species
	//units
		double s_len=0.0,s_mass=0.0;
		if(units::Consts::system()==units::System::LJ){
			s_len=1.0;
			s_mass=1.0;
		} else if(units::Consts::system()==units::System::AU){
			s_len=units::Ang2Bohr;
			s_mass=units::MPoME;
		} else if(units::Consts::system()==units::System::METAL){
			s_len=1.0;
			s_mass=1.0;
		} else throw std::runtime_error("Invalid units.");
	//misc
		bool error=false;
		
	try{
		//==== read header ====
		if(VASP_PRINT_STATUS>0) std::cout<<"read lattice vectors\n";
		fgets(input,string::M,reader);//name
		std::sscanf(fgets(input,string::M,reader),"%lf",&scale);
		std::sscanf(fgets(input,string::M,reader),"%lf %lf %lf",&lv(0,0),&lv(1,0),&lv(2,0));
		std::sscanf(fgets(input,string::M,reader),"%lf %lf %lf",&lv(0,1),&lv(1,1),&lv(2,1));
		std::sscanf(fgets(input,string::M,reader),"%lf %lf %lf",&lv(0,2),&lv(1,2),&lv(2,2));
		lv*=s_len*scale;
		static_cast<Cell&>(struc).init(lv);
		
		//==== read species ====
		if(VASP_PRINT_STATUS>0) std::cout<<"read species\n";
		//read number of species
		fgets(str_name, string::M, reader);
		fgets(str_number, string::M, reader);
		//read the species names and numbers
		token.read(str_name,string::WS);
		while(!token.end()) names.push_back(token.next());
		token.read(str_number,string::WS);
		while(!token.end()) nAtoms.push_back(std::atof(token.next().c_str()));
		if(names.size()==0 || nAtoms.size()==0 || names.size()!=nAtoms.size()){
			throw std::runtime_error("Invalid number of species");
		}
		const int nSpecies=names.size();
		//compute the total number
		nAtomsT=0; for(int i=0; i<nSpecies; ++i) nAtomsT+=nAtoms[i];
		
		//==== read coord ====
		if(VASP_PRINT_STATUS>0) std::cout<<"read coord\n";
		fgets(input, string::M, reader);
		if(input[0]=='D') direct=true;
		else direct=false;
		
		//====  print data to screen ==== 
		if(VASP_PRINT_DATA>0){
			std::cout<<"CELL    = \n"<<static_cast<Cell&>(struc)<<"\n";
			std::cout<<"DIRECT  = "<<(direct?"T":"F")<<"\n";
			std::cout<<"SPECIES = "; for(int i=0; i<names.size(); ++i) std::cout<<names[i]<<" "; std::cout<<"\n";
			std::cout<<"NUMBERS = "; for(int i=0; i<nAtoms.size(); ++i) std::cout<<nAtoms[i]<<" "; std::cout<<"\n";
		}
		
		//====  resize the simulation ==== 
		if(VASP_PRINT_STATUS>0) std::cout<<"allocating memory\n";
		Atom atomTl=atom;
		atomTl.frac=direct;
		struc.resize(nAtomsT,atomTl);
		
		//====  read positions ==== 
		if(VASP_PRINT_STATUS>0) std::cout<<"reading positions\n";
		for(int n=0; n<struc.nAtoms(); ++n){
			std::sscanf(
				fgets(input,string::M,reader),"%lf %lf %lf",
				&struc.posn(n)[0],&struc.posn(n)[1],&struc.posn(n)[2]
			);
		}
		
		//====  convert to cartesian coordinates (if necessary) ==== 
		if(VASP_PRINT_STATUS>0) std::cout<<"Converting to Cartesian coordinates\n";
		if(direct){
			for(int n=0; n<struc.nAtoms(); ++n){
				struc.posn(n)=struc.R()*struc.posn(n);
			}
		} else {
			for(int n=0; n<struc.nAtoms(); ++n){
				struc.posn(n)*=s_len;
			}
		}
		
		//==== set species ====
		if(atom.name){
			int count=0;
			for(int n=0; n<nSpecies; ++n){
				for(int m=0; m<nAtoms[n]; ++m){
					struc.name(count++)=names[n];
				}
			}
		}
		
		//==== set type ====
		if(atom.type){
			int count=0;
			for(int n=0; n<nSpecies; ++n){
				for(int m=0; m<nAtoms[n]; ++m){
					struc.type(count++)=n;
				}
			}
		}
		
		//==== set an ====
		if(atom.an && atom.name){
			int count=0;
			for(int n=0; n<nSpecies; ++n){
				for(int m=0; m<nAtoms[n]; ++m){
					struc.an(count)=ptable::an(struc.name(count).c_str());
					++count;
				}
			}
		}
		
		//==== set mass ====
		if(atom.mass && atom.an){
			int count=0;
			for(int n=0; n<nSpecies; ++n){
				for(int m=0; m<nAtoms[n]; ++m){
					struc.mass(count)=ptable::mass(struc.an(count))*s_mass;
					++count;
				}
			}
		} else if(atom.name && atom.mass){
			int count=0;
			for(int n=0; n<nSpecies; ++n){
				for(int m=0; m<nAtoms[n]; ++m){
					const int an=ptable::an(struc.name(count).c_str());
					struc.mass(count)=ptable::mass(an)*s_mass;
					++count;
				}
			}
		}
	}catch(std::exception& e){
		std::cout<<"ERROR in "<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	delete[] input;
	delete[] str_name;
	delete[] str_number;
	
	if(error) throw std::runtime_error("I/O Exception: Could not read data.");
}

void write(const char* file, const Atom& atom, const Structure& struc){
	static const char* funcName="read<AtomT>(const char*,Structure&)";
	if(VASP_PRINT_FUNC>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
	//open the file
	if(VASP_PRINT_STATUS>0) std::cout<<"opening file\n";
	FILE* writer=fopen(file,"w");
	if(writer==NULL) throw std::runtime_error(std::string("I/O Error: Could not open file: ")+std::string(file));
	//writer poscar
	write(writer,atom,struc);
	//close the file
	if(VASP_PRINT_STATUS>0) std::cout<<"closing file\n";
	fclose(writer);
	writer=NULL;
}
void write(FILE* writer, const Atom& atom, const Structure& struc){
	static const char* funcName="read<AtomT>(FILE*,const Atom&,Structure&)";
	if(VASP_PRINT_FUNC>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
	bool error=false;
	
	try{
		//find species names and numbers
		std::vector<std::string> names;
		std::vector<int> nAtoms;
		for(int i=0; i<struc.nAtoms(); ++i){
			int index=-1;
			for(int n=0; n<names.size(); ++n){
				if(struc.name(i)==names[n]){index=n;break;}
			}
			if(index<0){
				names.push_back(struc.name(i));
				nAtoms.push_back(1);
			} else ++nAtoms[index];
		}
		
		//writer header
		if(VASP_PRINT_STATUS>0) std::cout<<"writing header\n";
		fprintf(writer,"system\n");
		fprintf(writer,"%10.4f\n",1.0);
		for(int i=0; i<3; ++i){
			for(int j=0; j<3; ++j){
				fprintf(writer,"%11.6f ",struc.R()(j,i));
			}
			fprintf(writer,"\n");
		}
		for(int i=0; i<names.size(); ++i) fprintf(writer,"%s ",names[i].c_str()); fprintf(writer,"\n");
		for(int i=0; i<nAtoms.size(); ++i) fprintf(writer,"%i ",nAtoms[i]); fprintf(writer,"\n");
		
		//write positions
		if(VASP_PRINT_STATUS>0) std::cout<<"writing posns\n";
		if(!atom.frac){
			fprintf(writer,"Cart\n");
			for(int n=0; n<struc.nAtoms(); ++n){
				fprintf(writer, "%.8f %.8f %.8f\n",
					struc.posn(n)[0],
					struc.posn(n)[1],
					struc.posn(n)[2]
				);
			}
		} else {
			fprintf(writer,"Direct\n");
			for(int n=0; n<struc.nAtoms(); ++n){
				Eigen::Vector3d posn=struc.RInv()*struc.posn(n);
				fprintf(writer, "%.8f %.8f %.8f\n",posn[0],posn[1],posn[2]);
			}
		}
		
		if(VASP_PRINT_STATUS>0) std::cout<<"closing file\n";
		fclose(writer);
		writer=NULL;
	}catch(std::exception& e){
		std::cout<<"ERROR in "<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	if(error) throw std::runtime_error("I/O Exception: Could not write data.");
}

}

}