#include <iostream>
#include <cstring>
#include "format/format.hpp"

//**********************************************************************************************
//FILE_FORMAT struct
//**********************************************************************************************

FILE_FORMAT::type FILE_FORMAT::read(const char* str){
	if(std::strcmp(str,"POSCAR")==0) return FILE_FORMAT::POSCAR;
	else if(std::strcmp(str,"XDATCAR")==0) return FILE_FORMAT::XDATCAR;
	else if(std::strcmp(str,"LAMMPS")==0) return FILE_FORMAT::LAMMPS;
	else if(std::strcmp(str,"XYZ")==0) return FILE_FORMAT::XYZ;
	else if(std::strcmp(str,"CUBE")==0) return FILE_FORMAT::CUBE;
	else if(std::strcmp(str,"QE")==0) return FILE_FORMAT::QE;
	else if(std::strcmp(str,"CP2K")==0) return FILE_FORMAT::CP2K;
	else return FILE_FORMAT::NONE;
}

static const char* name(const FILE_FORMAT::type& format){
	switch(format){
		case FILE_FORMAT::POSCAR: return "POSCAR";
		case FILE_FORMAT::XDATCAR: return "XDATCAR";
		case FILE_FORMAT::LAMMPS: return "LAMMPS";
		case FILE_FORMAT::XYZ: return "XYZ";
		case FILE_FORMAT::CUBE: return "CUBE";
		case FILE_FORMAT::QE: return "QE";
		case FILE_FORMAT::CP2K: return "CP2K";
		default: return "NONE";
	}
}

std::ostream& operator<<(std::ostream& out, FILE_FORMAT::type& format){
	switch(format){
		case FILE_FORMAT::POSCAR: out<<"POSCAR"; break;
		case FILE_FORMAT::XDATCAR: out<<"XDATCAR"; break;
		case FILE_FORMAT::LAMMPS: out<<"LAMMPS"; break;
		case FILE_FORMAT::XYZ: out<<"XYZ"; break;
		case FILE_FORMAT::CUBE: out<<"CUBE"; break;
		case FILE_FORMAT::QE: out<<"QE"; break;
		case FILE_FORMAT::CP2K: out<<"CP2K"; break;
		default: out<<"NONE"; break;
	}
	return out;
}
