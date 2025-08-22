#pragma once
#ifndef CALC_FACTORY_HPP
#define CALC_FACTORY_HPP

// c++
#include <memory>
// sim
#include "sim/calc.hpp"

//==== constructors ====

std::shared_ptr<Calculator> make_calc(const Calculator::Name& name, double rc);
std::shared_ptr<Calculator> copy(const std::shared_ptr<Calculator>& ptr);

//==== reading ====

std::shared_ptr<Calculator>& read_calc(std::shared_ptr<Calculator>& pot, Token& token);
void coeff(std::vector<std::shared_ptr<Calculator> >& calcs, Token& token);

//==== output ====

std::ostream& operator<<(std::ostream& out, const std::shared_ptr<Calculator>& calc);

//==== memory ====

namespace serialize{

    template <> int nbytes(const std::shared_ptr<Calculator>& ptr);

    template <> int pack(const std::shared_ptr<Calculator>& ptr, char* arr);

    template <> int unpack(std::shared_ptr<Calculator>& ptr, const char* arr);
}

#endif