// c++
#include <iostream>
#include <stdexcept>
#include <chrono>
#include <random>
#include <memory>
#include <map>
// math
#include "math/const.hpp"
// util
#include "util/time.hpp"
// chem
#include "chem/units.hpp"
#include "chem/ptable.hpp"
// struc
#include "struc/structure.hpp"
// format
#include "format/format.hpp"
#include "format/file_struc.hpp"
#include "format/xyz_struc.hpp"
// str
#include "str/string.hpp"
#include "str/token.hpp"
#include "str/print.hpp"
// sim
#include "sim/job.hpp"
#include "sim/engine.hpp"
#include "sim/integrator.hpp"
#include "sim/calc_factory.hpp"
#include "sim/constraint_factory.hpp"

//****************************************************************************
// Main
//****************************************************************************

int main(int argc, char* argv[]){
	//units
		units::System unitsys=units::System::NONE;
	//files
		std::string fparam;
		std::string fstruc;
		char* input=new char[string::M];
		char* strbuf=new char[print::len_buf];
		Token token;
	//struc
		Atom atom;
		Structure struc;
		FILE_FORMAT::type format;
		int ntypes=0;
		std::map<std::string,int> types;
		Eigen::Vector3d shift=Eigen::Vector3d::Zero();
	//md
		Job job;
		Engine engine;
		std::shared_ptr<Integrator> intg;
		int nstep=0;
        int nprint=-1;
        int nwrite=-1;
		double ftol=math::constants::ZERO;
	//rand
		std::srand(std::time(NULL));
		int seed=std::chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator(seed);
	//misc
		bool error=false;
		Clock clock;
	
	try{
		//==== check the arguments ====
		if(argc!=2) throw std::invalid_argument("Torch::main(int,char**): Invalid number of arguments.");
		
		//==== open the parameter file ==== 
		fparam=argv[1];
		FILE* reader=fopen(fparam.c_str(),"r");
		if(reader==NULL) throw std::runtime_error("Torch::main(int,char**): Could not open parameter file.");
		
		//==== read the parameter file ==== 
		std::cout<<"reading general parameters\n";
		while(fgets(input,string::M,reader)!=NULL){
			string::trim_right(input,string::COMMENT);//trim comments
			Token token(input,string::WS); //split line into tokens
			if(token.end()) continue; //skip empty lines
			std::string tag=string::to_upper(token.next());
			if(tag=="JOB"){
				job=Job::read(string::to_upper(token.next()).c_str());
			} else if(tag=="UNITS"){
				unitsys=units::System::read(string::to_upper(token.next()).c_str());
				units::Consts::init(unitsys);
			} else if(tag=="ATOM"){
				atom.read(token);
				if(atom.posn()==false) throw std::invalid_argument("md::main(int,char**): Atom type missing position.");
				if(atom.mass()==false) throw std::invalid_argument("md::main(int,char**): Atom type missing mass.");
				if(atom.type()==false) throw std::invalid_argument("md::main(int,char**): Atom type missing type.");
			} else if(tag=="FORMAT"){//simulation format
				format=FILE_FORMAT::read(string::to_upper(token.next()).c_str());
			} else if(tag=="FSTRUC"){
				fstruc=token.next();
			} else if(tag=="ENGINE"){
				engine.read(token);
			} else if(tag=="INTEGRATOR" || tag=="INT"){
				Integrator::read(intg,token);
			} else if(tag=="CALC"){
				std::shared_ptr<Calculator> calc;
				read_calc(calc,token);
				engine.calcs().push_back(calc);
			} else if(tag=="CONSTRAINT"){
				std::shared_ptr<Constraint> constraint;
				read_constraint(constraint,token);
				engine.constraints().push_back(constraint);
			} else if(tag=="NSTEP"){
				nstep=std::atoi(token.next().c_str());
			} else if(tag=="NPRINT"){
				nprint=std::atoi(token.next().c_str());
			} else if(tag=="NWRITE"){
				nwrite=std::atoi(token.next().c_str());
			} else if(tag=="TYPE"){
				types[token.next()]=ntypes++;
			} else if(tag=="SHIFT"){
				shift[0]=std::atof(token.next().c_str());
				shift[1]=std::atof(token.next().c_str());
				shift[2]=std::atof(token.next().c_str());
			} else if(tag=="FTOL"){
				ftol=std::atof(token.next().c_str());
				if(ftol<0.0) throw std::invalid_argument("Invalid force tolerance.");
			}
		}
		
		//==== print ====
		std::cout<<print::buf(strbuf)<<"\n";
		std::cout<<print::title("MATHEMATICAL CONSTANTS",strbuf)<<"\n";
		std::printf("PI     = %.15f\n",math::constants::PI);
		std::printf("RadPI  = %.15f\n",math::constants::RadPI);
		std::printf("Rad2   = %.15f\n",math::constants::Rad2);
		std::printf("Log2   = %.15f\n",math::constants::LOG2);
		std::printf("Eps<D> = %.15e\n",std::numeric_limits<double>::epsilon());
		std::printf("Min<D> = %.15e\n",std::numeric_limits<double>::min());
		std::printf("Max<D> = %.15e\n",std::numeric_limits<double>::max());
		std::printf("Min<I> = %i\n",std::numeric_limits<int>::min());
		std::printf("Max<I> = %i\n",std::numeric_limits<int>::max());
		std::cout<<print::buf(strbuf)<<"\n";
		std::cout<<print::title("PHYSICAL CONSTANTS",strbuf)<<"\n";
		std::cout<<"ke = "<<units::Consts::ke()<<"\n";
		std::cout<<"kb = "<<units::Consts::kb()<<"\n";
		std::cout<<print::buf(strbuf)<<"\n";
		std::cout<<print::title("MD",strbuf)<<"\n";
		std::cout<<"job    = "<<job<<"\n";
		std::cout<<"units  = "<<unitsys<<"\n";
		std::cout<<"atom   = "<<atom<<"\n";
		std::cout<<"nstep  = "<<nstep<<"\n";
        std::cout<<"nprint = "<<nprint<<"\n";
        std::cout<<"nwrite = "<<nwrite<<"\n";
		std::cout<<"fstruc = "<<fstruc<<"\n";
		std::cout<<"format = "<<format<<"\n";
		std::cout<<"ftol   = "<<ftol<<"\n";
		std::cout<<"shift  = "<<shift.transpose()<<"\n";
		std::cout<<print::buf(strbuf)<<"\n";
		if(intg!=NULL){
			switch(intg->name()){
				case Integrator::Name::QUICKMIN: std::cout<<static_cast<const Quickmin&>(*intg)<<"\n"; break;
				case Integrator::Name::FIRE: std::cout<<static_cast<const Fire&>(*intg)<<"\n"; break;
				case Integrator::Name::VERLET: std::cout<<static_cast<const Verlet&>(*intg)<<"\n"; break;
				case Integrator::Name::VSCALE: std::cout<<static_cast<const VScale&>(*intg)<<"\n"; break;
				case Integrator::Name::BERENDSEN: std::cout<<static_cast<const Berendsen&>(*intg)<<"\n"; break;
				case Integrator::Name::LANGEVIN: std::cout<<static_cast<const Langevin&>(*intg)<<"\n"; break;
				default: std::cout<<"INTEGRATOR = NONE\n";
			}
		}
		std::cout<<print::buf(strbuf)<<"\n";
		std::cout<<print::title("TYPES",strbuf)<<"\n";
		for (const auto& [key, value] : types){
        	std::cout<<key<<" = "<<value<<"\n";
		}
		std::cout<<print::buf(strbuf)<<"\n";
		
		//==== read the structure ====
		std::cout<<"reading the structure\n";
		read_struc(fstruc.c_str(),format,atom,struc);
		
		//==== print structure ====
		std::cout<<struc<<"\n";
		for(int i=0; i<struc.nAtoms(); ++i){
			struc.vel(i).setZero();
		}

		//==== shift atoms ====
		if(shift.norm()>math::constants::ZERO){
			std::cout<<"shifting atoms\n";
			for(int i=0; i<struc.nAtoms(); ++i){
				struc.posn(i).noalias()+=shift;
				Cell::returnToCell(struc.posn(i),struc.posn(i),struc.R(),struc.RInv());
			}
		}
		
		//==== set the types ====
		std::cout<<"setting the types\n";
		for(int i=0; i<struc.nAtoms(); ++i) struc.type(i)=-1;
		for(int i=0; i<struc.nAtoms(); ++i){
			struc.type(i)=types.at(struc.name(i));
		}
		for(int i=0; i<struc.nAtoms(); ++i){
			if(struc.type(i)<0){
				throw std::invalid_argument("md::main(int,char**): No type for name "+struc.name(i)+"\n");
			}
		}

		//==== resize the engine ====
		std::cout<<"resizing the engine\n";
		/*int nt=-1;
		for(int i=0; i<struc.nAtoms(); ++i){
			std::cout<<struc.name(i)<<" "<<struc.type(i)<<"\n";
			if(struc.type(i)>nt) nt=struc.type(i);
		}
		nt++;
		std::cout<<"ntypes = "<<nt<<"\n";
		engine.resize(nt);*/
		engine.resize(ntypes);
		
		//==== read the coefficients ====
		std::cout<<"reading coefficients\n";
		std::rewind(reader);
		while(fgets(input,string::M,reader)!=NULL){
			string::trim_right(input,string::COMMENT);//trim comments
			Token token(input,string::WS); //split line into tokens
			if(token.end()) continue; //skip empty lines
			std::string tag=string::to_upper(token.next());
			if(tag=="COEFF"){
				coeff(engine.calcs(),token);
			} 
		}
		
		//==== initialize the engine ====
		std::cout<<"initializing the engine\n";
		if(intg!=NULL) struc.dt()=intg->dt();
		engine.init();
		NeighborList nlist(engine.rcmax());
		std::cout<<engine<<"\n";
		
		//==== close parameter file ==== 
		std::fclose(reader);
		reader=NULL;
		
		//==== compute ==== 
		clock.start();
		switch(job){
			case Job::SP:{
				std::cout<<"JOB - SP\n";
				std::cout<<"building nlist\n";
				nlist.build(struc);
				std::cout<<"computing energy\n";
				engine.init(struc);
				const double energy=engine.compute(struc,nlist);
				printf("energy = %.10f\n",energy);
				XYZ::write("out.xyz",atom,struc);
			}break;
			case Job::MD:{
				std::cout<<"JOB - MD\n";
				FILE* writer=fopen("out.xyz","w");
				if(writer==NULL) throw std::runtime_error("Could not open output file.");
				printf("N T KE PE TE\n");
				engine.init(struc);
				for(int t=0; t<nstep; ++t){
					struc.t()=t;
					if(t%engine.stride()==0) nlist.build(struc);
					//step
					intg->compute(struc,engine,nlist);
					//print
					if(t%nprint==0) printf("%i %4.5f %4.5f %4.5f %4.5f\n",t,struc.temp(),struc.ke(),struc.pe(),struc.ke()+struc.pe());
					//write
					if(t%nwrite==0) XYZ::write(writer,atom,struc);
				}
				fclose(writer); writer=NULL;
			}break;
			case Job::MIN:{
				std::cout<<"JOB - MIN\n";
				FILE* writer=fopen("out.xyz","w");
				if(writer==NULL) throw std::runtime_error("Could not open output file.");
				printf("N T KE PE TE\n");
				engine.init(struc);
				for(int t=0; t<nstep; ++t){
					struc.t()=t;
					if(t%engine.stride()==0) nlist.build(struc);
					//step
					intg->compute(struc,engine,nlist);
					//print
					if(t%nprint==0) printf("%i %4.5f %4.5f %4.5f %4.5f\n",t,struc.temp(),struc.ke(),struc.pe(),struc.ke()+struc.pe());
					//write
					if(t%nwrite==0) XYZ::write(writer,atom,struc);
					//compute total force
					double ftot=0;
					for(int i=0; i<struc.nAtoms(); ++i){
						ftot+=struc.force(i).squaredNorm();
					}
					ftot=std::sqrt(ftot/struc.nAtoms());
					if(ftot<ftol){
						std::cout<<"Convergence achieved.\n";
						break;
					}
				}
				fclose(writer); writer=NULL;
			}break;
			case Job::NUMDIFF:{
				std::cout<<"JOB - NUMDIFF\n";
				const double eps=1.0e-6;
				Structure strucP=struc;
				Structure strucM=struc;
				Structure strucA=struc;
				Engine engineP=engine;
				Engine engineM=engine;
				Engine engineA=engine;
				engineP.init(strucP);
				engineM.init(strucM);
				engineA.init(strucA);
				NeighborList nlistP(engine.rcmax());
				NeighborList nlistM(engine.rcmax());
				NeighborList nlistA(engine.rcmax());
				//compute force - analytical
				std::cout<<"computing force - analytical\n";
				nlistA.build(strucA);
				engineA.compute(strucA,nlistA);
				//compute force - numerical
				std::cout<<"computing force - numerical\n";
				for(int n=0; n<struc.nAtoms(); ++n){
					std::cout<<"atom "<<n<<"\n";
					//reset positions
					for(int i=0; i<struc.nAtoms(); ++i){
						strucP.posn(i)=struc.posn(i);
						strucM.posn(i)=struc.posn(i);
					}
					//store position
					for(int i=0; i<3; ++i){
						//perturb position
						strucP.posn(n)=struc.posn(n); strucP.posn(n)[i]+=eps;
						strucM.posn(n)=struc.posn(n); strucM.posn(n)[i]-=eps;
						//compute energy - plus
						nlistP.build(strucP);
						const double energyP=engineP.energy(strucP,nlistP);
						//compute energy - minus
						nlistM.build(strucP);
						const double energyM=engineM.energy(strucM,nlistM);
						//compute force
						struc.force(n)[i]=-0.5*(energyP-energyM)/eps;
					}
				}
				for(int n=0; n<struc.nAtoms(); ++n){
					std::cout<<struc.name(n)<<" "<<struc.force(n).transpose()<<"\n";
				}
				//compute error
				double error=0;
				for(int n=0; n<struc.nAtoms(); ++n){
					error+=(strucA.force(n)-struc.force(n)).squaredNorm();
				}
				error=sqrt(error/struc.nAtoms());
				//print forces
				for(int n=0; n<strucA.nAtoms(); ++n){
					std::cout<<"fa "<<strucA.name(n)<<" "<<strucA.force(n).transpose()<<"\n";
					std::cout<<"fn "<<struc.name(n)<<" "<<struc.force(n).transpose()<<"\n";
				}
				//print error
				std::cout<<"error - numdiff = "<<error<<"\n";
			}break;
			default:{
				std::cout<<"WARNING: Invalid job.";
			}break;
		}
		clock.stop();
		std::cout<<"time = "<<clock.duration()<<"\n";
		
	}catch(std::exception& e){
		std::cout<<"ERROR in md::main(int,char**):\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	//free memory
	std::cout<<"freeing memory\n";
	delete[] input;
	delete[] strbuf;
	
	if(error) return 1;
	else return 0;
}