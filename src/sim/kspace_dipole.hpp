#pragma once
#ifndef KSPACE_DIPOLE_HPP
#define KSPACE_DIPOLE_HPP

//mem
#include "mem/serialize.hpp"
// structure
#include "struc/structure.hpp"
// sim
#include "sim/kspace.hpp"

#ifndef KSPACED_PRINT_FUNC
#define KSPACED_PRINT_FUNC 0
#endif

#ifndef KSPACED_PRINT_STATUS
#define KSPACED_PRINT_STATUS 0
#endif

#ifndef KSPACED_PRINT_DATA
#define KSPACED_PRINT_DATA 0
#endif

namespace KSpace{

class Dipole: public Base{
private:
    double eps_{1.0};
    double m2_{0.0};//sum of squares of dipoles
    double vc_{0.0};//constant term
	mutable std::vector<double> c_;
	mutable std::vector<double> s_;
public:
    //==== constructors/destructors ====
    Dipole(){}
    virtual ~Dipole(){}
    
    //==== operators ====
    friend std::ostream& operator<<(std::ostream& out, const Dipole& d);
    
    //==== access ====
    double& eps(){return eps_;}
    const double& eps()const{return eps_;}
    const double& m2()const{return m2_;}
    const double& vc()const{return vc_;}
    
    //==== member functions ====
    void init(const Structure& struc);
    double energy(Structure& struc)const;
    double compute(Structure& struc)const;
    Eigen::Matrix3d& interMat(const Eigen::Vector3d& drv, Eigen::Matrix3d& mat);
};
//Estimate of the cutoff errors in the Ewald summation for dipolar systems
//doi: 10.1063/1.1398588

}

//**********************************************
// serialization
//**********************************************

namespace serialize{
    
    //**********************************************
    // byte measures
    //**********************************************
    
    template <> int nbytes(const KSpace::Dipole& obj);
    
    
    //**********************************************
    // packing
    //**********************************************
    
    template <> int pack(const KSpace::Dipole& obj, char* arr);
    
    //**********************************************
    // unpacking
    //**********************************************
    
    template <> int unpack(KSpace::Dipole& obj, const char* arr);
    
}

#endif