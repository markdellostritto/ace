#pragma once
#ifndef CP2K_STRUC_HPP
#define CP2K_STRUC_HPP

// eigen libraries
#include <Eigen/Dense>
// format
#include "format/cp2k.hpp"
// structure
#include "struc/structure.hpp"

#ifndef __cplusplus
	#error A C++ compiler is required
#endif

namespace CP2K{

//*****************************************************
//reading
//*****************************************************

void read(const char* file, const Atom& atom, Structure& struc);

}

#endif