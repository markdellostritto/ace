// c libraries
#include <cstring>
#include <cmath>
// c++ libraries
#include <iostream>
// Eigen
#include <Eigen/Dense>
// chem
#include "chem/ptable.hpp"
#include "chem/units.hpp"
// str
#include "str/string.hpp"
#include "str/token.hpp"
// format
#include "format/cube_struc.hpp"

namespace CUBE{
		
void read(const char* file, const Atom& atom, Structure& struc, Grid& grid){
	if(CUBE_PRINT_FUNC>0) std::cout<<"read(const char*,const Atom&,Structure&,Grid&):\n";
	//==== local variables ====
	//file i/o
		FILE* reader=NULL;
		char* input=new char[string::M];
	//structure
		Eigen::Vector3d origin;
		int natoms=0;
		int np=0;
		Eigen::Vector3i gridsize=Eigen::Vector3i::Zero();
		Eigen::Matrix3d lv,lvv;//lattice vectors
	//units
		double s_len=0.0;
		if(units::Consts::system()==units::System::LJ){
			s_len=1.0;
		} else if(units::Consts::system()==units::System::AU){
			s_len=1.0;
		} else if(units::Consts::system()==units::System::METAL){
			s_len=units::Bohr2Ang;
		} else throw std::runtime_error("Invalid units.");
		
	//==== open file ====
	if(CUBE_PRINT_STATUS>0) std::cout<<"opening file\n";
	reader=fopen(file,"r");
	if(reader==NULL) throw std::runtime_error("Could not open cube file");
	
	//==== read grid ====
	if(CUBE_PRINT_STATUS>0) std::cout<<"reading grid\n";
	fgets(input,string::M,reader);//skip line
	fgets(input,string::M,reader);//skip line
	//read natoms, origin
    std::sscanf(fgets(input,string::M,reader),"%i %lf %lf %lf",&natoms,&origin[0],&origin[1],&origin[2]);
	origin*=s_len;
	if(CUBE_PRINT_DATA>0) std::cout<<"natoms = "<<natoms<<"\n";
	if(CUBE_PRINT_DATA>0) std::cout<<"origin = "<<origin.transpose()<<"\n";
	//read grid
	std::sscanf(fgets(input,string::M,reader),"%i %lf %lf %lf",&gridsize[0],&lvv(0,0),&lvv(1,0),&lvv(2,0));
	std::sscanf(fgets(input,string::M,reader),"%i %lf %lf %lf",&gridsize[1],&lvv(0,1),&lvv(1,1),&lvv(2,1));
	std::sscanf(fgets(input,string::M,reader),"%i %lf %lf %lf",&gridsize[2],&lvv(0,2),&lvv(1,2),&lvv(2,2));
	if(CUBE_PRINT_DATA>0) std::cout<<"n = "<<gridsize.transpose()<<"\n";
	if(CUBE_PRINT_DATA>0) std::cout<<"lvv = \n"<<lvv<<"\n";
	//set units
	if(gridsize[0]<0){
		//lengths are in Angstroms
		if(units::Consts::system()==units::System::AU) s_len=units::Ang2Bohr;
		else if(units::Consts::system()==units::System::METAL) s_len=1.0;
		else throw std::runtime_error("Invalid units.");
	} else {
		//lengths are in Bohr
		if(units::Consts::system()==units::System::AU) s_len=1.0;
		else if(units::Consts::system()==units::System::METAL) s_len=units::Bohr2Ang;
		else throw std::runtime_error("Invalid units.");
	}
	lvv*=s_len;
	gridsize[0]=std::fabs(gridsize[0]);
	gridsize[1]=std::fabs(gridsize[1]);
	gridsize[2]=std::fabs(gridsize[2]);
	np=gridsize[0]*gridsize[1]*gridsize[2];
	//set lattice vector
	lv.col(0)=lvv.col(0)*gridsize[0];
	lv.col(1)=lvv.col(1)*gridsize[1];
	lv.col(2)=lvv.col(2)*gridsize[2];
	if(CUBE_PRINT_DATA>0) std::cout<<"s_len = "<<s_len<<"\n";
	if(CUBE_PRINT_DATA>1) std::cout<<"lv = \n"<<lv<<"\n";
	
	//=== read atoms ===
	if(CUBE_PRINT_STATUS>0) std::cout<<"reading atoms\n";
	//resize structure
	struc.resize(natoms,atom);
	//read atoms
	for(int i=0; i<natoms; ++i){
		int an; double q;
		Eigen::Vector3d p;
		std::sscanf(fgets(input,string::M,reader),"%i %lf %lf %lf %lf",&an,&q,&p[0],&p[1],&p[2]);
		const std::string name=ptable::name(an);
		if(atom.posn) struc.posn(i)=p*s_len;
		if(atom.an) struc.an(i)=an;
		if(atom.an) struc.charge(i)=q;
		if(atom.name) struc.name(i)=name;
	}
	
	//==== read grid ====
	if(CUBE_PRINT_STATUS>0) std::cout<<"reading grid\n";
	grid.resize(gridsize);
	grid.voxel()=lvv;
	grid.origin()=origin;
	int c=0;
	Token token;
	while(fgets(input,string::M,reader)){
		token.read(input,string::WS);
		while(!token.end()) grid.data()[c++]=std::atof(token.next().c_str());
	}
	
	delete[] input;
}

void write(const char* file, const Atom& atom, const Structure& struc, const Grid& grid){
	if(CUBE_PRINT_FUNC>0) std::cout<<"write(const char*,const Atom&,const Structure&,const Grid&):\n";
	//==== local variables ====
	//file i/o
		FILE* writer=NULL;
	//units
		double s_len=0.0,s_mass=1.0;
		if(units::Consts::system()==units::System::LJ){
			s_len=1.0;
		} else if(units::Consts::system()==units::System::AU){
			s_len=1.0;
		} else if(units::Consts::system()==units::System::METAL){
			s_len=units::Bohr2Ang;
		} else throw std::runtime_error("Invalid units.");
		
	//==== open file ====
	if(CUBE_PRINT_STATUS>0) std::cout<<"opening file\n";
	writer=fopen(file,"w");
	if(writer==NULL) throw std::runtime_error("Could not open cube file");
	
	//==== write grid ====
	//commment
	fprintf(writer,"test\n");
	fprintf(writer,"Density\n");
	//natoms origin
	const Eigen::Vector3d origin=grid.origin()*1.0/s_len;
	fprintf(writer,"%5i %12.6f %12.6f %12.6f\n",struc.nAtoms(),origin[0],origin[1],origin[2]);
	//print grid
	const Eigen::Matrix3d voxel=grid.voxel()*1.0/s_len;
	fprintf(writer,"%5i %12.6f %12.6f %12.6f\n",grid.n(0),voxel(0,0),voxel(0,1),voxel(0,2));
	fprintf(writer,"%5i %12.6f %12.6f %12.6f\n",grid.n(1),voxel(1,0),voxel(1,1),voxel(1,2));
	fprintf(writer,"%5i %12.6f %12.6f %12.6f\n",grid.n(2),voxel(2,0),voxel(2,1),voxel(2,2));
	
	//=== write atoms ===
	for(int i=0; i<struc.nAtoms(); ++i){
		const Eigen::Vector3d posn=struc.posn(i)*1.0/s_len;
		fprintf(writer,"%5i %12.6f %12.6f %12.6f %12.6f\n",
			struc.an(i),1.0*struc.an(i),
			posn[0],posn[1],posn[2]
		);
	}
	
	//==== write grid ====
	int c=0;
	for(int i=0; i<grid.n(0); ++i){
		for(int j=0; j<grid.n(1); ++j){
			for(int k=0; k<grid.n(2); ++k){
				fprintf(writer,"%13.5E ",grid.data(i,j,k));
				if(k%6==5) fprintf(writer,"\n");
			}
			fprintf(writer,"\n");
		}
	}
	
	//==== close file ====
	fclose(writer);
	writer=NULL;
}

}