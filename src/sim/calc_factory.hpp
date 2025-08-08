#pragma once
#ifndef CALC_FACTORY_HPP
#define CALC_FACTORY_HPP

// c++
#include <memory>
// sim
#include "sim/calc.hpp"

std::shared_ptr<Calculator> make_calc(const Calculator::Name& name, double rc);
std::shared_ptr<Calculator>& read_calc(std::shared_ptr<Calculator>& pot, Token& token);
void coeff(std::vector<std::shared_ptr<Calculator> >& calcs, Token& token);
std::ostream& operator<<(std::ostream& out, const std::shared_ptr<Calculator>& pot);

#endif