// string
#include "str/string.hpp"
// sim
#include "sim/constraint_factory.hpp"
// constraints
#include "sim/constraint_freeze.hpp"

std::shared_ptr<Constraint> make_constraint(const Constraint::Name& name){
    std::shared_ptr<Constraint> ptr;
    switch(name){
        case Constraint::Name::FREEZE:{
            ptr.reset(new ConstraintFreeze());
        } break;
        default:{
            throw std::invalid_argument("make_constraint(const Constraint::Name&): Invalid constraint.");
        }
    }
    return ptr;
}

std::shared_ptr<Constraint>& read_constraint(std::shared_ptr<Constraint>& calc, Token& token){
    const Constraint::Name name=Constraint::Name::read(string::to_upper(token.next()).c_str());
    switch(name){
        case Constraint::Name::FREEZE:{
            calc.reset(new ConstraintFreeze());
            static_cast<ConstraintFreeze&>(*calc).read(token);
        }break;
        default:{
            throw std::invalid_argument("read(std::shared_ptr<Constraint>&,Token&): Invalid potential.");
        }break;
    }
    return calc;
}

std::ostream& operator<<(std::ostream& out, const std::shared_ptr<Constraint>& calc){
    switch(calc->name()){
        case Constraint::Name::FREEZE: out<<static_cast<const Constraint&>(*calc); break;
        default: break;
    }
    return out;
}