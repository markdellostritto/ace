#pragma once
#ifndef CONSTRAINT_FACTORY_HPP
#define CONSTRAINT_FACTORY_HPP

// c++
#include <memory>
// sim
#include "sim/constraint.hpp"

std::shared_ptr<Constraint> make_constraint(const Constraint::Name& name);
std::shared_ptr<Constraint>& read_constraint(std::shared_ptr<Constraint>& pot, Token& token);
std::ostream& operator<<(std::ostream& out, const std::shared_ptr<Constraint>& constraint);

#endif