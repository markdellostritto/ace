// string
#include "str/string.hpp"
// sim
#include "sim/calc_factory.hpp"
// calculators
#include "sim/calc_lj_cut.hpp"
#include "sim/calc_coul_cut.hpp"
#include "sim/calc_coul_long.hpp"
#include "sim/calc_cgem_cut.hpp"
#include "sim/calc_cgem_long.hpp"

std::shared_ptr<Calculator> make_calc(const Calculator::Name& name, double rc){
    std::shared_ptr<Calculator> ptr;
    switch(name){
        case Calculator::Name::LJ_CUT:{
            ptr.reset(new CalcLJCut(rc));
        } break;
        case Calculator::Name::COUL_CUT:{
            ptr.reset(new CalcCoulCut(rc));
        } break;
        case Calculator::Name::COUL_LONG:{
            ptr.reset(new CalcCoulLong(rc));
        } break;
        case Calculator::Name::CGEM_CUT:{
            ptr.reset(new CalcCGemCut(rc));
        } break;
        case Calculator::Name::CGEM_LONG:{
            ptr.reset(new CalcCGemLong(rc));
        } break;
        default:{
            throw std::invalid_argument("make_calc(const Calculator::Name&): Invalid potential.");
        }
    }
    return ptr;
}

std::shared_ptr<Calculator>& read_calc(std::shared_ptr<Calculator>& calc, Token& token){
    const Calculator::Name name=Calculator::Name::read(string::to_upper(token.next()).c_str());
    switch(name){
        case Calculator::Name::LJ_CUT:{
            calc.reset(new CalcLJCut());
            static_cast<CalcLJCut&>(*calc).read(token);
        }break;
        case Calculator::Name::COUL_CUT:{
            calc.reset(new CalcCoulCut());
            static_cast<CalcCoulCut&>(*calc).read(token);
        }break;
        case Calculator::Name::COUL_LONG:{
            calc.reset(new CalcCoulLong());
            static_cast<CalcCoulLong&>(*calc).read(token);
        }break;
        case Calculator::Name::CGEM_CUT:{
            calc.reset(new CalcCGemCut());
            static_cast<CalcCGemCut&>(*calc).read(token);
        }break;
        case Calculator::Name::CGEM_LONG:{
            calc.reset(new CalcCGemLong());
            static_cast<CalcCGemLong&>(*calc).read(token);
        }break;
        default:{
            throw std::invalid_argument("read(std::shared_ptr<Calculator>&,Token&): Invalid potential.");
        }break;
    }
    return calc;
}

void coeff(std::vector<std::shared_ptr<Calculator> >& calcs, Token& token){
    Calculator::Name name=Calculator::Name::read(string::to_upper(token.next()).c_str());
    for(int i=0; i<calcs.size(); ++i){
        if(calcs[i]->name()==name){
            calcs[i]->coeff(token);
            break;
        }
    }
}

std::ostream& operator<<(std::ostream& out, const std::shared_ptr<Calculator>& calc){
    switch(calc->name()){
        case Calculator::Name::LJ_CUT: out<<static_cast<const CalcLJCut&>(*calc); break;
        case Calculator::Name::COUL_CUT: out<<static_cast<const CalcCoulCut&>(*calc); break;
        case Calculator::Name::COUL_LONG: out<<static_cast<const CalcCoulLong&>(*calc); break;
        case Calculator::Name::CGEM_CUT: out<<static_cast<const CalcCGemCut&>(*calc); break;
        case Calculator::Name::CGEM_LONG: out<<static_cast<const CalcCGemLong&>(*calc); break;
        default: break;
    }
    return out;
}