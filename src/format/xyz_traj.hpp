#pragma once
#ifndef XYZ_TRAJ_HPP
#define XYZ_TRAJ_HPP

// eigen libraries
#include <Eigen/Dense>
//format
#include "format/xyz.hpp"
// structure
#include "struc/trajectory.hpp"
#include "struc/interval.hpp"

namespace XYZ{
	
//unwrapping

void unwrap(Structure& struc);

void read(const char* file, const Interval& interval, const Atom& atom, Trajectory& traj);
void write(const char* file, const Interval& interval, const Atom& atom, const Trajectory& traj);

}

#endif