#pragma once
#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

// c++
#include <iosfwd>
#include <memory>
// chem
#include "chem/units.hpp"
// struc
#include "struc/structure.hpp"
// str
#include "str/token.hpp"
// sim
#include "sim/engine.hpp"

#ifndef INTEGRATOR_PRINT_FUNC
#define INTEGRATOR_PRINT_FUNC 0
#endif

//********************************************************************
// Integrator
//********************************************************************

class Integrator{
public:
	class Name{
	public:
		enum Type{
			QUICKMIN,
			FIRE,
			CG,
			VERLET,
			VSCALE,
			BERENDSEN,
			LANGEVIN,
			NONE
		};
		//constructor
		Name():t_(Type::NONE){}
		Name(Type t):t_(t){}
		//operators
		friend std::ostream& operator<<(std::ostream& out, const Name& name);
		operator Type()const{return t_;}
		//member functions
		static Name read(const char* str);
		static const char* name(const Name& name);
	private:
		Type t_;
	};
protected:
	static const double eps_;//divergence protection
	Name name_;//thermostat name
	double dt_;
	double dmax_;
public:
	//==== constructors/destructors ===
	Integrator():name_(Name::NONE),dmax_(0.1){}
	Integrator(const Name& name):name_(name),dmax_(0.1){}
	virtual ~Integrator(){}
	
	//==== operators ====
	friend std::ostream& operator<<(std::ostream& out, const Integrator& thermo);
	
	//==== access ====
	const Name& name()const{return name_;}
	double& dt(){return dt_;}
	const double& dt()const{return dt_;}
	double& dmax(){return dmax_;}
	const double& dmax()const{return dmax_;}
	
	//==== member functions ====
	void clear(){};
	
	//==== virtual functions ====
	virtual void read(Token& token);
	virtual void compute(Structure& struc, Engine& engine)=0;
	virtual void compute(Structure& struc, Engine& engine, const NeighborList& nlist)=0;
	virtual std::ostream& print(std::ostream& out)const{}
	
	//==== static functions ====
	static double ke(const Structure& struc);
	static double temp(const Structure& struc);
	static void compute(Structure& struc);
	static std::shared_ptr<Integrator>& make(const Name& name, std::shared_ptr<Integrator>& therm);
	static std::shared_ptr<Integrator>& read(std::shared_ptr<Integrator>& therm, Token& token);
};

//***********************************************************
// Quickmin
//***********************************************************

class Quickmin: public Integrator{
public:
	//==== constructors/destructors ====
	Quickmin():Integrator(Name::QUICKMIN){}
	~Quickmin(){}
	
	//==== operators ====
	friend std::ostream& operator<<(std::ostream& out, const Quickmin& thermo);
	
	//==== member functions ====
	void read(Token& token);
	void compute(Structure& struc, Engine& engine);
	void compute(Structure& struc, Engine& engine, const NeighborList& nlist);
	std::ostream& print(std::ostream& out)const;	
};

//***********************************************************
// Fire
//***********************************************************

class Fire: public Integrator{
private:
	int ndelay_;
	int npos_,nneg_,nmax_;
	double alpha0_;
	double alpha_;
	double falpha_;
	double fdtd_,fdti_;
	double dtmin_,dtmax_;
public:
	//==== constructors/destructors ====
	Fire():Integrator(Name::FIRE),npos_(0){}
	~Fire(){}
	
	//==== operators ====
	friend std::ostream& operator<<(std::ostream& out, const Fire& thermo);
	
	//==== access ====
	int& ndelay(){return ndelay_;}
	const int& ndelay()const{return ndelay_;}
	int& nmax(){return nmax_;}
	const int& nmax()const{return nmax_;}
	int& npos(){return npos_;}
	const int& npos()const{return npos_;}
	int& nneg(){return nneg_;}
	const int& nneg()const{return nneg_;}
	double& alpha0(){return alpha0_;}
	const double& alpha0()const{return alpha0_;}
	double& alpha(){return alpha_;}
	const double& alpha()const{return alpha_;}
	double& falpha(){return falpha_;}
	const double& falpha()const{return falpha_;}
	double& fdtd(){return fdtd_;}
	const double& fdtd()const{return fdtd_;}
	double& fdti(){return fdti_;}
	const double& fdti()const{return fdti_;}
	double& dtmin(){return dtmin_;}
	const double& dtmin()const{return dtmin_;}
	double& dtmax(){return dtmax_;}
	const double& dtmax()const{return dtmax_;}

	//==== member functions ====
	void read(Token& token);
	void compute(Structure& struc, Engine& engine);
	void compute(Structure& struc, Engine& engine, const NeighborList& nlist);
	std::ostream& print(std::ostream& out)const;	
};

//***********************************************************
// CG
//***********************************************************

class CG: public Integrator{
private:
	std::vector<Eigen::Vector3d> d_;
	std::vector<Eigen::Vector3d> f_;
public:
	//==== constructors/destructors ====
	CG():Integrator(Name::CG){}
	~CG(){}
	
	//==== operators ====
	friend std::ostream& operator<<(std::ostream& out, const CG& cg);
	
	//==== member functions ====
	void read(Token& token);
	void compute(Structure& struc, Engine& engine);
	void compute(Structure& struc, Engine& engine, const NeighborList& nlist);
	std::ostream& print(std::ostream& out)const;	
};

//***********************************************************
// Verlet
//***********************************************************

class Verlet: public Integrator{
public:
	//==== constructors/destructors ====
	Verlet():Integrator(Name::VERLET){}
	~Verlet(){}
	
	//==== operators ====
	friend std::ostream& operator<<(std::ostream& out, const Verlet& thermo);
	
	//==== member functions ====
	void read(Token& token);
	void compute(Structure& struc, Engine& engine);
	void compute(Structure& struc, Engine& engine, const NeighborList& nlist);
	std::ostream& print(std::ostream& out)const;	
};

//***********************************************************
// VScale
//***********************************************************

class VScale: public Integrator{
private:
	double T_;
	int tau_;
public:
	//==== constructors/destructors ====
	VScale():Integrator(Name::VSCALE){}
	~VScale(){}
	
	//==== operators ====
	friend std::ostream& operator<<(std::ostream& out, const VScale& thermo);
	
	//==== access ====
	double& T(){return T_;}
	const double& T()const{return T_;}
	int& tau(){return tau_;}
	const int& tau()const{return tau_;}
	
	//==== member functions ====
	void read(Token& token);
	void compute(Structure& struc, Engine& engine);
	void compute(Structure& struc, Engine& engine, const NeighborList& nlist);
	std::ostream& print(std::ostream& out)const;	
};

//***********************************************************
// Berendsen
//***********************************************************

class Berendsen: public Integrator{
private:
	double T_;
	double tau_;
public:
	//==== constructors/destructors ====
	Berendsen():Integrator(Name::BERENDSEN){}
	~Berendsen(){}
	
	//==== operators ====
	friend std::ostream& operator<<(std::ostream& out, const Berendsen& thermo);
	
	//==== access ====
	double& T(){return T_;}
	const double& T()const{return T_;}
	double& tau(){return tau_;}
	const double& tau()const{return tau_;}
	
	//==== member functions ====
	void read(Token& token);
	void compute(Structure& struc, Engine& engine);
	void compute(Structure& struc, Engine& engine, const NeighborList& nlist);
	std::ostream& print(std::ostream& out)const;	
};

//***********************************************************
// Langevin
//***********************************************************

class Langevin: public Integrator{
private:
	double T_;
	double gamma_;//damping constant
public:
	//==== constructors/destructors ====
	Langevin():Integrator(Name::LANGEVIN){}
	~Langevin(){}
	
	//==== operators ====
	friend std::ostream& operator<<(std::ostream& out, const Langevin& thermo);
	
	//==== access ====
	double& T(){return T_;}
	const double& T()const{return T_;}
	double& gamma(){return gamma_;}
	const double& gamma()const{return gamma_;}
	
	//==== member functions ====
	void read(Token& token);
	void compute(Structure& struc, Engine& engine);
	void compute(Structure& struc, Engine& engine, const NeighborList& nlist);
	std::ostream& print(std::ostream& out)const;	
};

#endif