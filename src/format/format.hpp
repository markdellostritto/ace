#pragma once
#ifndef FORMAT_HPP
#define FORMAT_HPP

#include <iosfwd>

//**********************************************************************************************
//FILE_FORMAT struct
//**********************************************************************************************

struct FILE_FORMAT{
	enum type{
		POSCAR,//VASP poscar file
		XDATCAR,//VASP xdatcar file
		LAMMPS,//LAMMPS dump file
		XYZ,//XYZ file
		CUBE,//cube file
		QE,//qe output
		CP2K,//cp2k output
		NONE//Unknown format
	};
	static FILE_FORMAT::type read(const char* str);
	static const char* name(const FILE_FORMAT::type& format);
};

std::ostream& operator<<(std::ostream& out, FILE_FORMAT::type& format);

#endif