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
		XYZ,//XYZ file
		CUBE,//cube file
		QE,//qe output
		NONE//Unknown format
	};
	static FILE_FORMAT::type read(const char* str);
	static const char* name(const FILE_FORMAT::type& format);
};

std::ostream& operator<<(std::ostream& out, FILE_FORMAT::type& format);

#endif