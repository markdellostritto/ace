#pragma once
#ifndef XYZ_STRUC_HPP
#define XYZ_STRUC_HPP

// eigen libraries
#include <Eigen/Dense>
// format
#include "format/xyz.hpp"
// structure
#include "struc/structure.hpp"

namespace XYZ{
	
//*****************************************************
//reading
//*****************************************************

void read(const char* file, const Atom& atom, Structure& struc);
void read(FILE* writer, const Atom& atom, Structure& struc);

//*****************************************************
//writing
//*****************************************************

void write(const char* file, const Atom& atom, const Structure& struc);
void write(FILE* writer, const Atom& atom, const Structure& struc);

}

#endif
