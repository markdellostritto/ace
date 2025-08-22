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
#include "sim/calc_cgemm_cut.hpp"
#include "sim/calc_cgemm_long.hpp"

//==== constructors ====

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
        case Calculator::Name::CGEMM_CUT:{
            ptr.reset(new CalcCGemmCut(rc));
        } break;
        case Calculator::Name::CGEMM_LONG:{
            ptr.reset(new CalcCGemmLong(rc));
        } break;
        default:{
            throw std::invalid_argument("make_calc(const Calculator::Name&): Invalid potential.");
        }
    }
    return ptr;
}

std::shared_ptr<Calculator> copy(const std::shared_ptr<Calculator>& ptr){
    std::shared_ptr<Calculator> ptrnew;
    switch(ptr->name()){
        case Calculator::Name::LJ_CUT:{
            ptrnew=std::make_shared<CalcLJCut>(static_cast<const CalcLJCut&>(*ptr));
        } break;
        case Calculator::Name::COUL_CUT:{
            ptrnew=std::make_shared<CalcCoulCut>(static_cast<const CalcCoulCut&>(*ptr));
        } break;
        case Calculator::Name::COUL_LONG:{
            ptrnew=std::make_shared<CalcCoulLong>(static_cast<const CalcCoulLong&>(*ptr));
        } break;
        case Calculator::Name::CGEM_CUT:{
            ptrnew=std::make_shared<CalcCGemCut>(static_cast<const CalcCGemCut&>(*ptr));
        } break;
        case Calculator::Name::CGEM_LONG:{
            ptrnew=std::make_shared<CalcCGemLong>(static_cast<const CalcCGemLong&>(*ptr));
        } break;
        case Calculator::Name::CGEMM_CUT:{
            ptrnew=std::make_shared<CalcCGemmCut>(static_cast<const CalcCGemmCut&>(*ptr));
        } break;
        case Calculator::Name::CGEMM_LONG:{
            ptrnew=std::make_shared<CalcCGemmLong>(static_cast<const CalcCGemmLong&>(*ptr));
        } break;
        default:{
            throw std::invalid_argument("copy(const std::shared_ptr<Calculator>&): Invalid potential.");
        }
    }
    return ptrnew;
}

//==== reading ====

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
        case Calculator::Name::CGEMM_CUT:{
            calc.reset(new CalcCGemmCut());
            static_cast<CalcCGemmCut&>(*calc).read(token);
        }break;
        case Calculator::Name::CGEMM_LONG:{
            calc.reset(new CalcCGemmLong());
            static_cast<CalcCGemmLong&>(*calc).read(token);
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

//==== output ====

std::ostream& operator<<(std::ostream& out, const std::shared_ptr<Calculator>& calc){
    switch(calc->name()){
        case Calculator::Name::LJ_CUT: out<<static_cast<const CalcLJCut&>(*calc); break;
        case Calculator::Name::COUL_CUT: out<<static_cast<const CalcCoulCut&>(*calc); break;
        case Calculator::Name::COUL_LONG: out<<static_cast<const CalcCoulLong&>(*calc); break;
        case Calculator::Name::CGEM_CUT: out<<static_cast<const CalcCGemCut&>(*calc); break;
        case Calculator::Name::CGEM_LONG: out<<static_cast<const CalcCGemLong&>(*calc); break;
        case Calculator::Name::CGEMM_CUT: out<<static_cast<const CalcCGemmCut&>(*calc); break;
        case Calculator::Name::CGEMM_LONG: out<<static_cast<const CalcCGemmLong&>(*calc); break;
        default: break;
    }
    return out;
}

//==== memory ====

namespace serialize{

    template <> int nbytes(const std::shared_ptr<Calculator>& ptr){
        int size=0;
        size+=nbytes(ptr->name());
        switch(ptr->name()){
            case Calculator::Name::LJ_CUT:{
                size+=nbytes(static_cast<const CalcLJCut&>(*ptr));
            }break;
            case Calculator::Name::COUL_CUT:{
                size+=nbytes(static_cast<const CalcCoulCut&>(*ptr));
            }break;
            case Calculator::Name::COUL_LONG:{
                size+=nbytes(static_cast<const CalcCoulLong&>(*ptr));
            }break;
            case Calculator::Name::CGEM_CUT:{
                size+=nbytes(static_cast<const CalcCGemCut&>(*ptr));
            }break;
            case Calculator::Name::CGEM_LONG:{
                size+=nbytes(static_cast<const CalcCGemLong&>(*ptr));
            }break;
            case Calculator::Name::CGEMM_CUT:{
                size+=nbytes(static_cast<const CalcCGemmCut&>(*ptr));
            }break;
            case Calculator::Name::CGEMM_LONG:{
                size+=nbytes(static_cast<const CalcCGemmLong&>(*ptr));
            }break;
            default: break;
        }
        return size;
    }
    
    template <> int pack(const std::shared_ptr<Calculator>& ptr, char* arr){
        int pos=0;
        Calculator::Name name=ptr->name();
        pos+=pack(name,arr+pos);
        switch(name){
            case Calculator::Name::LJ_CUT:{
                pos+=pack(static_cast<const CalcLJCut&>(*ptr),arr+pos);
            }break;
            case Calculator::Name::COUL_CUT:{
                pos+=pack(static_cast<const CalcCoulCut&>(*ptr),arr+pos);
            }break;
            case Calculator::Name::COUL_LONG:{
                pos+=pack(static_cast<const CalcCoulLong&>(*ptr),arr+pos);
            }break;
            case Calculator::Name::CGEM_CUT:{
                pos+=pack(static_cast<const CalcCGemCut&>(*ptr),arr+pos);
            }break;
            case Calculator::Name::CGEM_LONG:{
                pos+=pack(static_cast<const CalcCGemLong&>(*ptr),arr+pos);
            }break;
            case Calculator::Name::CGEMM_CUT:{
                pos+=pack(static_cast<const CalcCGemmCut&>(*ptr),arr+pos);
            }break;
            case Calculator::Name::CGEMM_LONG:{
                pos+=pack(static_cast<const CalcCGemmLong&>(*ptr),arr+pos);
            }break;
            default: break;
        }
        return pos;
    }
    
    template <> int unpack(std::shared_ptr<Calculator>& ptr, const char* arr){
        int pos=0;
        Calculator::Name name=Calculator::Name::NONE;
        pos+=unpack(name,arr+pos);
        switch(name){
            case Calculator::Name::LJ_CUT:{
                ptr=std::make_shared<CalcLJCut>();
                pos+=unpack(static_cast<const CalcLJCut&>(*ptr),arr+pos);
            }break;
            case Calculator::Name::COUL_CUT:{
                ptr=std::make_shared<CalcCoulCut>();
                pos+=unpack(static_cast<const CalcCoulCut&>(*ptr),arr+pos);
            }break;
            case Calculator::Name::COUL_LONG:{
                ptr=std::make_shared<CalcCoulLong>();
                pos+=unpack(static_cast<const CalcCoulLong&>(*ptr),arr+pos);
            }break;
            case Calculator::Name::CGEM_CUT:{
                ptr=std::make_shared<CalcCGemCut>();
                pos+=unpack(static_cast<const CalcCGemCut&>(*ptr),arr+pos);
            }break;
            case Calculator::Name::CGEM_LONG:{
                ptr=std::make_shared<CalcCGemLong>();
                pos+=unpack(static_cast<const CalcCGemLong&>(*ptr),arr+pos);
            }break;
            case Calculator::Name::CGEMM_CUT:{
                ptr=std::make_shared<CalcCGemmCut>();
                pos+=unpack(static_cast<const CalcCGemmCut&>(*ptr),arr+pos);
            }break;
            case Calculator::Name::CGEMM_LONG:{
                ptr=std::make_shared<CalcCGemmLong>();
                pos+=unpack(static_cast<const CalcCGemmLong&>(*ptr),arr+pos);
            }break;
            default: break;
        }
        return pos;
    }
    
}