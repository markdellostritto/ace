// string
#include "str/string.hpp"
// sim
#include "sim/constraint_factory.hpp"
// constraints
#include "sim/constraint_freeze.hpp"

//==== constructors ====

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

std::shared_ptr<Constraint> copy(const std::shared_ptr<Constraint>& ptr){
    switch(ptr->name()){
        case Constraint::Name::FREEZE:{
            return std::make_shared<ConstraintFreeze>(static_cast<const ConstraintFreeze&>(*ptr));
        } break;
        default:{
            throw std::invalid_argument("make_constraint(const Constraint::Name&): Invalid constraint.");
        }
    }
}

//==== reading ====

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

//==== output ====

std::ostream& operator<<(std::ostream& out, const std::shared_ptr<Constraint>& calc){
    switch(calc->name()){
        case Constraint::Name::FREEZE: out<<static_cast<const Constraint&>(*calc); break;
        default: break;
    }
    return out;
}

//==== memory ====

namespace serialize{

    template <> int nbytes(const std::shared_ptr<Constraint>& ptr){
        int size=0;
        size+=nbytes(ptr->name());
        switch(ptr->name()){
            case Constraint::Name::FREEZE:{
                size+=nbytes(static_cast<const ConstraintFreeze&>(*ptr));
            }break;
            default: break;
        }
        return size;
    }

    template <> int pack(const std::shared_ptr<Constraint>& ptr, char* arr){
        int pos=0;
        Constraint::Name name=ptr->name();
        pos+=pack(name,arr+pos);
        switch(name){
            case Constraint::Name::FREEZE:{
                pos+=pack(static_cast<const ConstraintFreeze&>(*ptr),arr+pos);
            }break;
            default: break;
        }
        return pos;
    }

    template <> int unpack(std::shared_ptr<Constraint>& ptr, const char* arr){
        int pos=0;
        Constraint::Name name=Constraint::Name::NONE;
        pos+=unpack(name,arr+pos);
        switch(name){
            case Constraint::Name::FREEZE:{
                ptr=std::make_shared<ConstraintFreeze>();
                pos+=unpack(static_cast<ConstraintFreeze&>(*ptr),arr+pos);
            }break;
            default: break;
        }
        return pos;
    }

}