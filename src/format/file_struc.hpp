#ifndef FILE_STRUC_HPP
#define FILE_STRUC_HPP

// ann - structure
#include "struc/structure.hpp"
// ann - file i/o
#include "format/format.hpp"

Structure& read_struc(const char* file, FILE_FORMAT::type format, const Atom& atom, Structure& struc);
const Structure& write_struc(const char* file, FILE_FORMAT::type format, const Atom& atom, const Structure& struc);

#endif