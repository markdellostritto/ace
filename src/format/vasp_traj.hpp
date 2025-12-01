#pragma once
#ifndef VASP_SIM_HPP
#define VASP_SIM_HPP

// c++ libraries
#include <vector>
#include <string>
// format
#include "format/vasp.hpp"
// structure
#include "format/format.hpp"
#include "struc/trajectory.hpp"

namespace VASP{

namespace XDATCAR{

static const char* NAMESPACE_LOCAL="XDATCAR";
void read(const char* file, const Interval& interval, const Atom& atom, Trajectory& traj);
void write(const char* file, const Interval& interval, const Atom& atom, const Trajectory& traj);

}

}

#endif