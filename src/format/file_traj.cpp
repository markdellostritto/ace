// c++ libraries
#include <iostream>
#include <exception>
// file i/o
#include "format/file_traj.hpp"
#include "format/vasp_traj.hpp"
#include "format/lammps_traj.hpp"
#include "format/xyz_traj.hpp"

Trajectory& read_sim(const char* file, FILE_FORMAT::type format, const Interval& interval, const Atom& atom, Trajectory& traj){
	switch(format){
		case FILE_FORMAT::XDATCAR:
			VASP::XDATCAR::read(file,interval,atom,traj);
		break;
		case FILE_FORMAT::LAMMPS:
			LAMMPS::DUMP::read(file,interval,atom,traj);
		break;
		case FILE_FORMAT::XYZ:
			XYZ::read(file,interval,atom,traj);
		break;
		default:
			throw std::invalid_argument("ERROR in read_sim(const char*,FILE_Format::type,const Atom&,Trajectory&): invalid file format.");
		break;
	}
	return traj;
}

const Trajectory& write_sim(const char* file, FILE_FORMAT::type format, const Interval& interval, const Atom& atom, const Trajectory& traj){
	switch(format){
		case FILE_FORMAT::XDATCAR:
			VASP::XDATCAR::write(file,interval,atom,traj);
		break;
		case FILE_FORMAT::LAMMPS:
			LAMMPS::DUMP::write(file,interval,atom,traj);
		break;
		case FILE_FORMAT::XYZ:
			XYZ::write(file,interval,atom,traj);
		break;
		default:
			throw std::invalid_argument("ERROR in write_sim(const char*,FILE_Format::type,const Atom&,const Trajectory&): invalid file format.");
		break;
	}
	return traj;
}