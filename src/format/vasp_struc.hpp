#pragma once
#ifndef VASP_STRUC_HPP
#define VASP_STRUC_HPP

// c++ libraries
#include <vector>
#include <string>
// format
#include "format/vasp.hpp"
// structure
#include "format/format.hpp"
#include "struc/structure.hpp"

namespace VASP{

//*****************************************************
//POSCAR
//*****************************************************

namespace POSCAR{

static const char* NAMESPACE_LOCAL="POSCAR";
void read(const char* file, const Atom& atom, Structure& struc);
void read(FILE* reader, const Atom& atom, Structure& struc);
void write(const char* file, const Atom& atom, const Structure& struc);
void write(FILE* writer, const Atom& atom, const Structure& struc);

}

}


#endif
