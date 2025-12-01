#pragma once
#ifndef FILE_SIM_HPP
#define FILE_SIM_HPP

// structure
#include "struc/trajectory.hpp"
#include "struc/interval.hpp"
// file i/o
#include "format/format.hpp"

Trajectory& read_sim(const char* file, FILE_FORMAT::type format, const Interval& interval, const Atom& atom, Trajectory& traj);
const Trajectory& write_sim(const char* file, FILE_FORMAT::type format, const Interval& interval, const Atom& atom, const Trajectory& traj);

#endif