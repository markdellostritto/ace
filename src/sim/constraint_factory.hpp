#pragma once
#ifndef CONSTRAINT_FACTORY_HPP
#define CONSTRAINT_FACTORY_HPP

// c++
#include <memory>
// sim
#include "sim/constraint.hpp"

//==== constructors ====

std::shared_ptr<Constraint> make_constraint(const Constraint::Name& name);
std::shared_ptr<Constraint> copy(const std::shared_ptr<Constraint>& ptr);

//==== reading ====

std::shared_ptr<Constraint>& read_constraint(std::shared_ptr<Constraint>& pot, Token& token);

//==== output ====

std::ostream& operator<<(std::ostream& out, const std::shared_ptr<Constraint>& constraint);

//==== memory ====

namespace serialize{

    template <> int nbytes(const std::shared_ptr<Constraint>& ptr);

    template <> int pack(const std::shared_ptr<Constraint>& ptr, char* arr);

    template <> int unpack(std::shared_ptr<Constraint>& ptr, const char* arr);
}

#endif