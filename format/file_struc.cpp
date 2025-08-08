// c++ libraries
#include <iostream>
#include <exception>
// ann - structure
#include "struc/structure.hpp"
// ann - file i/o
#include "format/file_struc.hpp"
#include "format/vasp_struc.hpp"
#include "format/xyz_struc.hpp"
#include "format/qe_struc.hpp"

Structure& read_struc(const char* file, FILE_FORMAT::type format, const Atom& atom, Structure& struc){
	switch(format){
		case FILE_FORMAT::POSCAR:
			VASP::POSCAR::read(file,atom,struc);
		break;
		case FILE_FORMAT::XYZ:
			XYZ::read(file,atom,struc);
		break;
		case FILE_FORMAT::QE:
			QE::OUT::read(file,atom,struc);
		break;
		default:
			throw std::invalid_argument("ERROR in read(const char*,FILE_Format::type,const Atom&,Structure&): invalid file format.");
		break;
	}
	return struc;
}

const Structure& write_struc(const char* file, FILE_FORMAT::type format, const Atom& atom, const Structure& struc){
	switch(format){
		case FILE_FORMAT::POSCAR:
			VASP::POSCAR::write(file,atom,struc);
		break;
		case FILE_FORMAT::XYZ:
			XYZ::write(file,atom,struc);
		break;
		default:
			throw std::invalid_argument("ERROR in write(const char*,FILE_Format::type,const Atom&,const Structure&): invalid file format.");
		break;
	}
	return struc;
}