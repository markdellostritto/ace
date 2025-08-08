// c
#include <cstring>
// c++
#include <iostream>
#include <string>
// math
#include "math/const.hpp"
// str
#include "str/string.hpp"
// sim
#include "sim/integrator.hpp"

using math::constants::ZERO;

//***********************************************************
// Name (Integrator)
//***********************************************************

Integrator::Name Integrator::Name::read(const char* str){
	if(std::strcmp(str,"QUICKMIN")==0) return Integrator::Name::QUICKMIN;
	else if(std::strcmp(str,"CG")==0) return Integrator::Name::CG;
	else if(std::strcmp(str,"FIRE")==0) return Integrator::Name::FIRE;
	else if(std::strcmp(str,"VERLET")==0) return Integrator::Name::VERLET;
	else if(std::strcmp(str,"VSCALE")==0) return Integrator::Name::VSCALE;
	else if(std::strcmp(str,"BERENDSEN")==0) return Integrator::Name::BERENDSEN;
	else if(std::strcmp(str,"LANGEVIN")==0) return Integrator::Name::LANGEVIN;
	else return Integrator::Name::NONE;
}

const char* Integrator::Name::name(const Integrator::Name& t){
	switch(t){
		case Integrator::Name::QUICKMIN: return "QUICKMIN";
		case Integrator::Name::CG: return "CG";
		case Integrator::Name::FIRE: return "FIRE";
		case Integrator::Name::VERLET: return "VERLET";
		case Integrator::Name::VSCALE: return "VSCALE";
		case Integrator::Name::BERENDSEN: return "BERENDSEN";
		case Integrator::Name::LANGEVIN: return "LANGEVIN";
		default: return "NONE";
	}
}

std::ostream& operator<<(std::ostream& out, const Integrator::Name& t){
	switch(t){
		case Integrator::Name::QUICKMIN: out<<"QUICKMIN"; break;
		case Integrator::Name::CG: out<<"CG"; break;
		case Integrator::Name::FIRE: out<<"FIRE"; break;
		case Integrator::Name::VERLET: out<<"VERLET"; break;
		case Integrator::Name::VSCALE: out<<"VSCALE"; break;
		case Integrator::Name::BERENDSEN: out<<"BERENDSEN"; break;
		case Integrator::Name::LANGEVIN: out<<"LANGEVIN"; break;
		default: out<<"NONE"; break;
	}
	return out;
}

//***********************************************************
// Integrator
//***********************************************************

//==== constants ====

const double Integrator::eps_=1e-14;

//==== operators ====

std::ostream& operator<<(std::ostream& out, const Integrator& intg){
	return out;
}

//==== static functions ====

double Integrator::ke(const Structure& struc){
	double energy=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		energy+=struc.mass(i)*struc.vel(i).squaredNorm();
	}
	return 0.5*energy;
}

double Integrator::temp(const Structure& struc){
	return struc.ke()*(2.0/3.0)/(struc.nAtoms()*units::Consts::kb());
}

void Integrator::read(Token& token){
	dt_=std::atof(token.next().c_str());
}

static std::shared_ptr<Integrator>& make(const Integrator::Name& name, std::shared_ptr<Integrator>& intg){
	switch(name){
		case Integrator::Name::QUICKMIN:{
			intg.reset(new Quickmin());
		} break;
		case Integrator::Name::FIRE:{
			intg.reset(new Fire());
		} break;
		case Integrator::Name::VERLET:{
			intg.reset(new Verlet());
		} break;
		case Integrator::Name::VSCALE:{
			intg.reset(new VScale());
		} break;
		case Integrator::Name::BERENDSEN:{
			intg.reset(new Berendsen());
		} break;
		default:{
			intg.reset();
		} break;
	}
	return intg;
}

std::shared_ptr<Integrator>& Integrator::read(std::shared_ptr<Integrator>& intg, Token& token){
	//read name
	Integrator::Name name=Integrator::Name::read(string::to_upper(token.next()).c_str());
	switch(name){
		case Integrator::Name::QUICKMIN: {
			intg.reset(new Quickmin());
			static_cast<Quickmin&>(*intg).read(token);
		} break;
		case Integrator::Name::CG: {
			intg.reset(new CG());
			static_cast<CG&>(*intg).read(token);
		} break;
		case Integrator::Name::FIRE: {
			intg.reset(new Fire());
			static_cast<Fire&>(*intg).read(token);
		} break;
		case Integrator::Name::VERLET: {
			intg.reset(new Verlet());
			static_cast<Verlet&>(*intg).read(token);
		} break;
		case Integrator::Name::VSCALE: {
			intg.reset(new VScale());
			static_cast<VScale&>(*intg).read(token);
		} break;
		case Integrator::Name::BERENDSEN: {
			intg.reset(new Berendsen());
			static_cast<Berendsen&>(*intg).read(token);
		} break;
		case Integrator::Name::LANGEVIN: {
			intg.reset(new Langevin());
			static_cast<Langevin&>(*intg).read(token);
		} break;
		default:{
			throw std::invalid_argument("Integrator::read(Token&): invalid integrator name.");
		} break;
	}
	return intg;
}

//***********************************************************
// Quickmin
//***********************************************************

void Quickmin::compute(Structure& struc, Engine& engine){
	//compute force
	for(int i=0; i<struc.nAtoms(); ++i) struc.force(i).setZero();
	struc.pe()=engine.compute(struc);
	//project velocity along force
	for(int i=0; i<struc.nAtoms(); ++i){
		//project velocity along force
		Eigen::Vector3d fhat=Eigen::Vector3d::Zero();
		const double fnorm=struc.force(i).norm();
		if(fnorm>math::constants::ZERO) fhat=struc.force(i)/struc.force(i).norm();
		const double vdotf=struc.vel(i).dot(fhat);
		struc.vel(i)=fhat*vdotf;
		if(struc.vel(i).dot(fhat)<0.0) struc.vel(i).setZero();
		//euler step
		Eigen::Vector3d drv=struc.vel(i)*dt_;
		const double drn=drv.norm();
		if(drn>dmax_) drv*=dmax_/drn;
		struc.posn(i).noalias()+=drv;
		struc.vel(i).noalias()+=struc.force(i)*dt_;
	}
	//compute KE, T
	struc.ke()=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.ke()+=struc.mass(i)*struc.vel(i).squaredNorm();
	}
	struc.ke()*=0.5;
	struc.temp()=struc.ke()*(2.0/3.0)/(struc.nAtoms()*units::Consts::kb());
	//increment
	++struc.t();
}

void Quickmin::read(Token& token){
	Integrator::read(token);
}

std::ostream& Quickmin::print(std::ostream& out)const{
	return out;
}

std::ostream& operator<<(std::ostream& out, const Quickmin& intg){
	return out<<intg.name()<<" "<<intg.dt();
}

//***********************************************************
// Fire
//***********************************************************

void Fire::compute(Structure& struc, Engine& engine){
	if(struc.t()==0){
		for(int i=0; i<struc.nAtoms(); ++i) struc.vel(i).setZero();
		for(int i=0; i<struc.nAtoms(); ++i) struc.force(i).setZero();
		struc.pe()=engine.compute(struc);
		alpha_=alpha0_;
		npos_=0;
		nneg_=0;
	}

	//**** FIRE STEP ****

	//compute P
	double P=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		P+=struc.force(i).dot(struc.vel(i));
	}
	//std::cout<<"P = "<<P<<"\n";
	if(P>0.0){
		npos_++;
		nneg_=0;
		if(npos_>ndelay_){
			dt_=std::min(dt_*fdti_,dtmax_);
			alpha_*=falpha_;
		}
	} else {
		nneg_++;
		npos_=0;
		//if(nneg_>nmax_) break;
		if(struc.t()>ndelay_){
			dt_=std::max(dt_*fdtd_,dtmin_);
			alpha_=alpha0_;
		}
		for(int i=0; i<struc.nAtoms(); ++i){
			struc.posn(i).noalias()-=0.5*dt_*struc.vel(i);
			struc.vel(i).setZero();
		}
	}
	//std::cout<<"P "<<P<<" dt "<<dt_<<" alpha "<<alpha_<<"\n";
	
	//**** MD STEP ****

	//Velocity Verlet
	/*
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.vel(i).noalias()+=0.5*dt_*struc.force(i)/struc.mass(i);
		const double vnorm=struc.vel(i).norm();
		struc.vel(i)*=(1.0-alpha_);
		struc.vel(i).noalias()+=alpha_*vnorm*struc.force(i)/(1.0e-16+struc.force(i).norm());
		struc.posn(i).noalias()+=struc.vel(i)*dt_;
		struc.force(i).setZero();
	}
	struc.pe()=engine.compute(struc);
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.vel(i).noalias()+=0.5*dt_*struc.force(i)/struc.mass(i);
	}
	*/
	
	//Euler (explicit)
	/*
	for(int i=0; i<struc.nAtoms(); ++i){
		const double vnorm=struc.vel(i).norm();
		struc.vel(i)*=(1.0-alpha_);
		struc.vel(i).noalias()+=alpha_*vnorm*struc.force(i)/(1.0e-16+struc.force(i).norm());
		struc.posn(i).noalias()+=struc.vel(i)*dt_;
		struc.vel(i).noalias()+=dt_*struc.force(i)/struc.mass(i);
		struc.force(i).setZero();
	}
	struc.pe()=engine.compute(struc);
	*/

	//Euler (semi-implicit)
	for(int i=0; i<struc.nAtoms(); ++i){
		const double vnorm=struc.vel(i).norm();
		struc.vel(i).noalias()+=dt_*struc.force(i)/struc.mass(i);
		struc.vel(i)*=(1.0-alpha_);
		struc.vel(i).noalias()+=alpha_*vnorm*struc.force(i)/(1.0e-16+struc.force(i).norm());
		Eigen::Vector3d drv=struc.vel(i)*dt_;
		const double drn=drv.norm();
		if(drn>dmax_) drv*=dmax_/drn;
		struc.posn(i).noalias()+=drv;
		struc.force(i).setZero();
	}
	struc.pe()=engine.compute(struc);
	
	//compute KE, T
	struc.ke()=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.ke()+=struc.mass(i)*struc.vel(i).squaredNorm();
	}
	struc.ke()*=0.5;
	struc.temp()=struc.ke()*(2.0/3.0)/(struc.nAtoms()*units::Consts::kb());
	//increment
	++struc.t();
}

void Fire::read(Token& token){
	//int fire dt dtmax fdt alpha falpha ndelay 
	Integrator::read(token);
	dtmin_=std::atof(token.next().c_str());
	dtmax_=std::atof(token.next().c_str());
	fdtd_=std::atof(token.next().c_str());
	fdti_=std::atof(token.next().c_str());
	alpha_=std::atof(token.next().c_str());
	falpha_=std::atof(token.next().c_str());
	ndelay_=std::atoi(token.next().c_str());
	nmax_=std::atoi(token.next().c_str());
	alpha0_=alpha_;
	npos_=0;
	nneg_=0;
}

std::ostream& Fire::print(std::ostream& out)const{
	return out;
}

std::ostream& operator<<(std::ostream& out, const Fire& intg){
	return out<<intg.name()<<" "<<intg.dt()<<" "<<intg.dtmin()<<" "<<intg.dtmax()
		<<" "<<intg.fdtd()<<" "<<intg.fdti()
		<<" "<<intg.alpha()<<" "<<intg.falpha()
		<<" "<<intg.ndelay()<<" "<<intg.nmax();
}

//***********************************************************
// CG
//***********************************************************

void CG::compute(Structure& struc, Engine& engine){
	if(struc.t()==0){
		//compute the force
		for(int i=0; i<struc.nAtoms(); ++i) struc.force(i).setZero();
		struc.pe()=engine.compute(struc);
		//store the force
		f_.resize(struc.nAtoms());
		for(int i=0; i<struc.nAtoms(); ++i){
			f_[i]=struc.force(i);
		}
		//set the search direction
		d_.resize(struc.nAtoms());
		for(int i=0; i<struc.nAtoms(); ++i){
			d_[i]=struc.force(i);
			const double norm=d_[i].norm();
			if(norm>ZERO) d_[i]/=norm;
		}
	}

	//compute new position
	std::vector<double> alpha(struc.nAtoms(),0.0);
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.posn(i).noalias()+=dt_*d_[i];
	}
	
	//store force
	for(int i=0; i<struc.nAtoms(); ++i){
		f_[i]=struc.force(i);
	}
	//compute the force
	for(int i=0; i<struc.nAtoms(); ++i) struc.force(i).setZero();
	struc.pe()=engine.compute(struc);

	//compute new direction
	for(int i=0; i<struc.nAtoms(); ++i){
		const double fnorm=f_[i].squaredNorm();
		if(fnorm>0) d_[i]*=struc.force(i).dot(struc.force(i)-f_[i])/fnorm;
		else d_[i]*=0;
		d_[i].noalias()+=struc.force(i);
	}
	
	//compute KE, T
	struc.ke()=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.ke()+=struc.mass(i)*struc.vel(i).squaredNorm();
	}
	struc.ke()*=0.5;
	struc.temp()=struc.ke()*(2.0/3.0)/(struc.nAtoms()*units::Consts::kb());
	//increment
	++struc.t();
}

void CG::read(Token& token){
	Integrator::read(token);
}

std::ostream& CG::print(std::ostream& out)const{
	return out;
}

std::ostream& operator<<(std::ostream& out, const CG& intg){
	return out<<intg.name()<<" "<<intg.dt();
}

//***********************************************************
// Verlet
//***********************************************************

void Verlet::compute(Structure& struc, Engine& engine){
	if(INTEGRATOR_PRINT_FUNC>0) std::cout<<"Verlet::compute(Structure&,Engine&):\n";
	//first half-step
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.vel(i).noalias()+=0.5*dt_*struc.force(i)/struc.mass(i);
		struc.posn(i).noalias()+=struc.vel(i)*dt_;
		Cell::returnToCell(struc.posn(i),struc.posn(i),struc.R(),struc.RInv());
	}
	//compute force
	struc.pe()=engine.compute(struc);
	//second half-step
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.vel(i).noalias()+=0.5*dt_*struc.force(i)/struc.mass(i);
	}
	//compute KE, T
	struc.ke()=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.ke()+=struc.mass(i)*struc.vel(i).squaredNorm();
	}
	struc.ke()*=0.5;
	struc.temp()=struc.ke()*(2.0/3.0)/(struc.nAtoms()*units::Consts::kb());
	//increment
	++struc.t();
}

void Verlet::read(Token& token){
	Integrator::read(token);
}

std::ostream& Verlet::print(std::ostream& out)const{
	return out;
}

std::ostream& operator<<(std::ostream& out, const Verlet& intg){
	return out<<intg.name()<<" "<<intg.dt();
}

//***********************************************************
// VScale
//***********************************************************

void VScale::compute(Structure& struc, Engine& engine){
	//first half-step
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.vel(i).noalias()+=0.5*dt_*struc.force(i)/struc.mass(i);
		struc.posn(i).noalias()+=struc.vel(i)*dt_;
		Cell::returnToCell(struc.posn(i),struc.posn(i),struc.R(),struc.RInv());
	}
	//compute force
	struc.pe()=engine.compute(struc);
	//second half-step
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.vel(i).noalias()+=0.5*dt_*struc.force(i)/struc.mass(i);
	}
	//compute KE, T
	struc.ke()=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.ke()+=struc.mass(i)*struc.vel(i).squaredNorm();
	}
	struc.ke()*=0.5;
	struc.temp()=struc.ke()*(2.0/3.0)/(struc.nAtoms()*units::Consts::kb());
	//alter velocities
	if(struc.t()%tau_==0){
		const double fac=sqrt(T_/(struc.temp()+1e-6));
		for(int i=0; i<struc.nAtoms(); ++i){
			struc.vel(i)*=fac;
		}
	}
	//compute KE, T
	struc.ke()=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.ke()+=struc.mass(i)*struc.vel(i).squaredNorm();
	}
	struc.ke()*=0.5;
	struc.temp()=struc.ke()*(2.0/3.0)/(struc.nAtoms()*units::Consts::kb());
	//increment
	++struc.t();
}

void VScale::read(Token& token){
	Integrator::read(token);
	T_=std::atof(token.next().c_str());
	tau_=std::atoi(token.next().c_str());
}

std::ostream& VScale::print(std::ostream& out)const{
	return out;
}

std::ostream& operator<<(std::ostream& out, const VScale& intg){
	return out<<intg.name()<<" "<<intg.dt()<<" "<<intg.T()<<" "<<intg.tau();
}

//***********************************************************
// Berendsen
//***********************************************************

void Berendsen::compute(Structure& struc, Engine& engine){
	//first half-step
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.vel(i).noalias()+=0.5*dt_*struc.force(i)/struc.mass(i);
		struc.posn(i).noalias()+=struc.vel(i)*dt_;
		Cell::returnToCell(struc.posn(i),struc.posn(i),struc.R(),struc.RInv());
	}
	//compute force
	struc.pe()=engine.compute(struc);
	//second half-step
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.vel(i).noalias()+=0.5*dt_*struc.force(i)/struc.mass(i);
	}
	//compute KE, T
	struc.ke()=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.ke()+=struc.mass(i)*struc.vel(i).squaredNorm();
	}
	struc.ke()*=0.5;
	struc.temp()=struc.ke()*(2.0/3.0)/(struc.nAtoms()*units::Consts::kb());
	//alter velocities
	const double T=struc.temp();
	const double t=struc.t();
	const double dt=struc.dt();
	const double fac=sqrt(1.0+dt/tau_*(T_/(T*(t-0.5*dt))-1.0));
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.vel(i)*=fac;
	}
	//compute KE, T
	struc.ke()=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.ke()+=struc.mass(i)*struc.vel(i).squaredNorm();
	}
	struc.ke()*=0.5;
	struc.temp()=struc.ke()*(2.0/3.0)/(struc.nAtoms()*units::Consts::kb());
	//increment
	++struc.t();
}

void Berendsen::read(Token& token){
	Integrator::read(token);
	T_=std::atof(token.next().c_str());
	tau_=std::atof(token.next().c_str());
}

std::ostream& Berendsen::print(std::ostream& out)const{
	return out;
}

std::ostream& operator<<(std::ostream& out, const Berendsen& intg){
	return out<<intg.name()<<" "<<intg.dt()<<" "<<intg.T()<<" "<<intg.tau();
}

//***********************************************************
// Berendsen
//***********************************************************

void Langevin::compute(Structure& struc, Engine& engine){
	//first half-step
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.vel(i).noalias()+=0.5*dt_*struc.force(i)/struc.mass(i);
		struc.posn(i).noalias()+=struc.vel(i)*dt_;
		Cell::returnToCell(struc.posn(i),struc.posn(i),struc.R(),struc.RInv());
	}
	//compute force
	struc.pe()=engine.compute(struc);
	//add damping and random forces
	const double c=sqrt((T_*units::Consts::kb())/(dt_*gamma_));
	for(int i=0; i<struc.nAtoms(); ++i){
		//damping force
		struc.force(i).noalias()-=struc.mass(i)/gamma_*struc.vel(i);
		//random force
		const Eigen::Vector3d rv=Eigen::Vector3d::Random();
		struc.force(i).noalias()+=c*sqrt(struc.mass(i))*rv/rv.norm();
	}
	//second half-step
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.vel(i).noalias()+=0.5*dt_*struc.force(i)/struc.mass(i);
	}
	//compute KE, T
	struc.ke()=0;
	for(int i=0; i<struc.nAtoms(); ++i){
		struc.ke()+=struc.mass(i)*struc.vel(i).squaredNorm();
	}
	struc.ke()*=0.5;
	struc.temp()=struc.ke()*(2.0/3.0)/(struc.nAtoms()*units::Consts::kb());
	//increment
	++struc.t();
}

void Langevin::read(Token& token){
	Integrator::read(token);
	T_=std::atof(token.next().c_str());
	gamma_=std::atof(token.next().c_str());
}

std::ostream& Langevin::print(std::ostream& out)const{
	return out;
}

std::ostream& operator<<(std::ostream& out, const Langevin& intg){
	return out<<intg.name()<<" "<<intg.dt()<<" "<<intg.T()<<" "<<intg.gamma();
}