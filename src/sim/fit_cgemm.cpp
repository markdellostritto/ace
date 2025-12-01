// c++
#include <iostream>
#include <string>
#include <random>
#include <chrono>
#include <map>
#include <random>
// math
#include "math/const.hpp"
// chem
#include "chem/units.hpp"
// str
#include "str/print.hpp"
#include "str/string.hpp"
#include "str/token.hpp"
// calcs
#include "sim/calc_cgemm_cut.hpp"
// format
#include "format/cube_struc.hpp"
//sim
#include "sim/fit_cgemm.hpp"
#include "sim/constraint_freeze.hpp"
#include "sim/calc_factory.hpp"
#include "sim/constraint_factory.hpp"
//nlopt
#include <nlopt.hpp>

//************************************************************
// Constants
//************************************************************

using math::constants::PI;

//************************************************************
// Global Variables
//************************************************************

int step=0;
int count=0;

bool opt_lambda=false;
bool opt_radius=false;
bool opt_aOver=false;
bool opt_rRep=false;
bool opt_weight=false;
double sigma=1.0;
double min_lambda=1.0e-6;
double min_radius=1.0e-6;
double min_aOver=1.0e-6;
double max_lambda=HUGE_VAL;
double max_radius=HUGE_VAL;
double max_aOver=HUGE_VAL;

//*****************************************************
// CGemmType
//*****************************************************

//==== operators ====

std::ostream& operator<<(std::ostream& out, const CGemmType& type){
    return out<<type.name()<<" "<<type.index()<<" "<<type.radius()<<" "<<type.rcut()<<" "<<type.aOver()<<" "<<type.aRep();
}

//==== member functions ====

void CGemmType::read(Token& token){
    name_=token.next();
    mass_=std::atof(token.next().c_str());
    radius_=std::atof(token.next().c_str());
    aOver_=std::atof(token.next().c_str());
    aRep_=std::atof(token.next().c_str());
    if(radius_<0) throw std::invalid_argument("CGemmType::read(Token&): invalid mass.");
    if(radius_<0) throw std::invalid_argument("CGemmType::read(Token&): invalid radius.");
    if(aOver_<0) throw std::invalid_argument("CGemmType::read(Token&): invalid overlap amplitude.");
    if(aRep_<0) throw std::invalid_argument("CGemmType::read(Token&): invalid repulsive amplitude.");
}

//*****************************************************
// Objective Function
//*****************************************************

double objf(const std::vector<double> &x, std::vector<double> &grad, void* fData){
    FunctionData* fData_=static_cast<FunctionData*>(fData);
    const int ntypes=fData_->types.size();
    //std::cout<<"x = "; for(int i=0; i<x.size(); ++i) std::cout<<x[i]<<" "; std::cout<<"\n";
    int c;
    
    //==== unpack the parameters - CGemm ====
    CalcCGemmCut& calcCGemmCut=static_cast<CalcCGemmCut&>(*fData_->engine.calcs().front());
    calcCGemmCut.lambdaC()=1.0;
    calcCGemmCut.lambdaS()=1.0;
    c=0;
    if(opt_lambda){
        calcCGemmCut.lambdaC()=x[c++];
        calcCGemmCut.lambdaS()=x[c++];
    }
    if(opt_radius) for(int i=0; i<ntypes; ++i) calcCGemmCut.radius()[i]=x[c++];
    else for(int i=0; i<ntypes; ++i) calcCGemmCut.radius()[i]=fData_->types[i].radius();
    if(opt_aOver) for(int i=0; i<ntypes; ++i) calcCGemmCut.aOver()(i,i)=x[c++];
    else for(int i=0; i<ntypes; ++i) calcCGemmCut.aOver()(i,i)=fData_->types[i].aOver();
    for(int i=0; i<ntypes; ++i) calcCGemmCut.aRep()(i,i)=fData_->types[i].aRep();
    if(opt_rRep) calcCGemmCut.rRep()=x[c++];
    calcCGemmCut.init();

    //==== unpack the parameters - types ====
    c=0;
    if(opt_lambda) c+=2;
    if(opt_radius) for(int i=0; i<ntypes; ++i) fData_->types[i].radius()=x[c++];
    if(opt_aOver) for(int i=0; i<ntypes; ++i) fData_->types[i].aOver()=x[c++];
    if(opt_rRep) c+=1;
    
    //==== loop over each structure ====
    Engine& engine=fData_->engine;
    double t_avg=0;
    double t_max=0;
    double f_avg=0;
    double f_max=0;
    double error=0;
    double np=0;
    for(int n=0; n<fData_->data.size(); ++n){
        Structure& strucA=fData_->data[n].strucA();
        Structure& strucN=fData_->data[n].strucN();
        Grid& espA=fData_->data[n].espA();
        Grid& espN=fData_->data[n].espN();

        //set the constraints
        std::vector<int> indices;
        for(int m=0; m<strucN.nAtoms(); ++m){
            if(strucN.name(m)!="X"){
                indices.push_back(m);
            }
        }
        engine.constraints().front()->indices()=indices;
        
        //reset the electron positions
        const int nCores=strucA.nAtoms();
        const int nTotal=strucN.nAtoms();
        const int nShells=nTotal-nCores;
        //find average posn
        Eigen::Vector3d avg=Eigen::Vector3d::Zero();
        for(int j=0; j<nCores; ++j){
            avg.noalias()+=strucA.posn(j);
        }
        avg*=1.0/nCores;
        double dev=0;
        for(int j=0; j<nCores; ++j){
            dev+=(avg-strucA.posn(j)).squaredNorm();
        }
        dev=std::sqrt(dev/nCores);
        //std::cout<<"struc "<<n<<"\n";
        //std::cout<<"avg = "<<avg.transpose()<<"\n";
        //std::cout<<"dev = "<<dev<<"\n";
        //set the electron positions
        if(nShells<nCores){
            for(int j=0; j<nShells; ++j){
                strucN.posn(nCores+j)=strucA.posn(j);
            }
        } else {
            //add electrons to atoms
            for(int j=0; j<nCores; ++j){
                strucN.posn(nCores+j)=strucA.posn(j);
            }
            //add extra electrons to random posns near the center
            for(int j=nCores; j<nShells; ++j){
                strucN.posn(nCores+j)=avg;
                strucN.posn(nCores+j).noalias()+=0.1*dev*Eigen::Vector3d::Random();
            }
            /*const int nExtra=nShells-nCores;
            const int nAttempts=10;
            //std::cout<<"nExtra = "<<nExtra<<"\n";
            std::vector<int> cores(nExtra,-1);
            int icore=0;
            for(int j=nCores; j<nShells; ++j){
                //choose core
                bool match=false;
                for(int ii=0; ii<nAttempts; ++ii){
                    cores[icore]=std::rand()%nCores;
                    for(int i=0; i<icore; ++i){
                        if(cores[icore]==cores[i]){
                            match=true; break;
                        }
                    }
                    if(!match) break;
                }
                strucN.posn(nCores+j)=strucA.posn(icore);
                strucN.posn(nCores+j).noalias()+=0.1*Eigen::Vector3d::Random();
                ++icore;
            }*/
        }

        //relax the system
        if(PRINT_STRUC>0){
            std::cout<<strucN.nAtoms()<<"\nmolecule\n";
            for(int m=0; m<strucN.nAtoms(); ++m){
                std::cout<<strucN.name(m)
                    <<" "<<strucN.type(m)
                    <<" "<<strucN.charge(m)
                    <<" "<<fData_->types[strucN.type(m)].radius()
                    <<" "<<fData_->types[strucN.type(m)].rcut()
                    <<" "<<fData_->types[strucN.type(m)].aOver()
                    <<" "<<fData_->types[strucN.type(m)].aRep()
                    <<" "<<strucN.posn(m).transpose()
                    <<"\n";
            }
        }
        int t=0;
        double ftot=0;
        strucN.t()=0;
        for(t=0; t<fData_->nsteps; ++t){
            //compute step
            fData_->intg->compute(strucN,engine);
            //compute total force
            ftot=0;
            for(int m=0; m<strucN.nAtoms(); ++m){
                ftot+=strucN.force(m).squaredNorm();
            }
            ftot=std::sqrt(ftot/strucN.nAtoms());
            if(ftot<fData_->ftol) break;
        }
        t_avg+=t;
        f_avg+=ftot;
        if(t>t_max) t_max=t;
        if(ftot>f_max) f_max=ftot;
        if(PRINT_STRUC>0){
            std::cout<<strucN.nAtoms()<<"\nmolecule\n";
            for(int m=0; m<strucN.nAtoms(); ++m){
                std::cout<<strucN.name(m)
                    <<" "<<strucN.type(m)
                    <<" "<<strucN.charge(m)
                    <<" "<<fData_->types[strucN.type(m)].radius()
                    <<" "<<fData_->types[strucN.type(m)].rcut()
                    <<" "<<fData_->types[strucN.type(m)].aOver()
                    <<" "<<fData_->types[strucN.type(m)].aRep()
                    <<" "<<strucN.posn(m).transpose()
                    <<"\n";
            }
        }

        //compute the potential
        const Eigen::Vector3i& npts=espA.n();
        const double ke=units::Consts::ke();
        for(int i=0; i<npts[0]; ++i){
            for(int j=0; j<npts[1]; ++j){
                for(int k=0; k<npts[2]; ++k){
                    const Eigen::Vector3d indexD=(Eigen::Vector3d()<<i,j,k).finished();
                    const Eigen::Vector3d r=espA.origin()+espA.voxel()*indexD;
                    double& data=espN.data(i,j,k);
                    data=0.0;
                    for(int m=0; m<strucN.nAtoms(); ++m){
                        const double dr=(r-strucN.posn(m)).norm();
                        const double R=fData_->types[strucN.type(m)].radius();
                        const double alpha=calcCGemmCut.lambdaC()/(2.0*R*R);
                        data+=ke*strucN.charge(m)/dr*std::erf(std::sqrt(alpha)*dr);
                    }
                }
            }
        }

        //compute the error
        for(int i=0; i<npts[0]; ++i){
            for(int j=0; j<npts[1]; ++j){
                for(int k=0; k<npts[2]; ++k){
                    //compute the grid position
                    const Eigen::Vector3d indexD=(Eigen::Vector3d()<<i,j,k).finished();
                    const Eigen::Vector3d r=espA.origin()+espA.voxel()*indexD;
                    //compute the difference
                    const double diff=espA.data(i,j,k)-espN.data(i,j,k);
                    //compute the weight
                    double w=1;
                    if(opt_weight){
                        for(int m=0; m<strucN.nAtoms(); ++m){
                            const double dr=(r-strucN.posn(m)).norm();
                            const double Rvdw=fData_->types[strucN.type(m)].rcut();
                            w*=0.5*(tanh((dr-Rvdw)/sigma)+1.0);
                        }
                    }
                    //compute the error
                    error+=std::fabs(w*diff);
                }
            }
        }
        np+=espA.np();
    }
    error/=np;
    t_avg/=fData_->data.size();
    f_avg/=fData_->data.size();
    
    //==== print error ====
    //std::cout<<"step "<<step<<" error "<<error<<" f_avg "<<f_avg<<" f_max "<<f_max<<" t_avg "<<t_avg<<" t_max "<<t_max<<"\n";
    std::printf(
        "%4i %10.6e %10.6e %10.6e %7i %7i",
        step,error,f_avg,f_max,
        static_cast<int>(std::round(t_avg)),
        static_cast<int>(std::round(t_max))
    );
    for(int i=0; i<x.size(); ++i) std::printf("%12.8f ",x[i]);
    std::printf("\n");
    step++;
    
    //==== return error ====
    return error;
}

//*****************************************************
// Main
//*****************************************************

int main(int argc, char* argv[]){
    //units
        units::System unitsys=units::System::NONE;
    //files
        std::string fparam;
        std::string fstruc;
        FILE* reader=NULL;
        char* input=new char[string::M];
        char* strbuf=new char[print::len_buf];
        Token token;
    //atom
        Atom atom;
        atom.name()=true; 
        atom.type()=true; 
        atom.an()=true;
        atom.charge()=true; 
        atom.mass()=true; 
        atom.posn()=true; 
        atom.force()=true; 
        atom.vel()=true;
    //cube files
        std::string cubelist;
        std::vector<std::string> cubefiles;
        std::vector<int> qtot;
    //fit
        double rc=0;//cutoff
        double lambdaC=1.0;//initial guess for lambdaC
        double lambdaS=1.0;//initial guess for lambdaS
        double rRep=0.05;
        double tol=0.0;
        int miter=0;
        nlopt::algorithm algo;
        double srcut=1.0;
    //function data
        FunctionData functionData;
        std::vector<CubeData>& data=functionData.data;
        Engine& engine=functionData.engine;
        std::vector<CGemmType>& types=functionData.types;
        std::shared_ptr<Integrator>& intg=functionData.intg;
    //rand
		std::srand(std::time(NULL));
		int seed=std::chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator(seed);
	//misc
		bool error=false;
	
    try{
        //==== check the arguments ====
        if(argc!=2) throw std::invalid_argument("Torch::main(int,char**): Invalid number of arguments.");

        //==== open the parameter file ==== 
		fparam=argv[1];
        reader=fopen(fparam.c_str(),"r");
		if(reader==NULL) throw std::runtime_error("Torch::main(int,char**): Could not open parameter file.");
		
        //==== read the parameter file ==== 
		std::cout<<"reading general parameters\n";
		while(fgets(input,string::M,reader)!=NULL){
            string::trim_right(input,string::COMMENT);//trim comments
			Token token(input,string::WS); //split line into tokens
			if(token.end()) continue; //skip empty lines
			std::string tag=string::to_upper(token.next());
			if(tag=="TYPE"){
                types.push_back(CGemmType());
                types.back().read(token);
            } else if(tag=="UNITS"){
				unitsys=units::System::read(string::to_upper(token.next()).c_str());
				units::Consts::init(unitsys);
			} else if(tag=="CUBEFILES"){
                cubelist=token.next();
            } else if(tag=="ENGINE"){
				engine.read(token);
			} else if(tag=="RCUT" || tag=="RC"){
				rc=std::atof(token.next().c_str());
			} else if(tag=="LAMBDAC"){
				lambdaC=std::atof(token.next().c_str());
			} else if(tag=="LAMBDAS"){
				lambdaS=std::atof(token.next().c_str());
			} else if(tag=="TOL" || tag=="TOLERANCE"){
				tol=std::atof(token.next().c_str());
			} else if(tag=="MITER" || tag=="MAX_ITER"){
				miter=std::atoi(token.next().c_str());
			} else if(tag=="NSTEPS"){
				functionData.nsteps=std::atoi(token.next().c_str());
			} else if(tag=="FTOL"){
				functionData.ftol=std::atof(token.next().c_str());
			} else if(tag=="INTEGRATOR" || tag=="INTG"){
				Integrator::read(intg,token);
			} else if(tag=="ALGO"){
                std::string salgo=string::to_upper(token.next());
                if(salgo=="CRS") algo=nlopt::GN_CRS2_LM;
                else if(salgo=="ISRES") algo=nlopt::GN_ISRES;
                else if(salgo=="ESCH") algo=nlopt::GN_ESCH;
                else if(salgo=="SBPLX") algo=nlopt::LN_SBPLX;
                else if(salgo=="BOBYQA") algo=nlopt::LN_BOBYQA;
                else if(salgo=="COBYLA") algo=nlopt::LN_COBYLA;
                else if(salgo=="PRAXIS") algo=nlopt::LN_PRAXIS;
                else if(salgo=="NELDERMEAD") algo=nlopt::LN_NELDERMEAD;
                else throw std::invalid_argument("Invalid optimization algorithm.");
            } else if(tag=="OPT_LAMBDA"){
				opt_lambda=string::boolean(token.next().c_str());
			} else if(tag=="OPT_RADIUS"){
				opt_radius=string::boolean(token.next().c_str());
			} else if(tag=="OPT_AOVER"){
				opt_aOver=string::boolean(token.next().c_str());
            } else if(tag=="OPT_RREP"){
				opt_rRep=string::boolean(token.next().c_str());
			} else if(tag=="OPT_WEIGHT" || tag=="OPT_WEIGHTED"){
                opt_weight=string::boolean(token.next().c_str());
            } else if(tag=="MIN_LAMBDA"){
				min_lambda=std::atof(token.next().c_str());
			} else if(tag=="MIN_RADIUS"){
                min_radius=std::atof(token.next().c_str());
			} else if(tag=="MIN_AOVER"){
                min_aOver=std::atof(token.next().c_str());
			} else if(tag=="MAX_LAMBDA"){
				max_lambda=std::atof(token.next().c_str());
			} else if(tag=="MAX_RADIUS"){
                max_radius=std::atof(token.next().c_str());
			} else if(tag=="MAX_AOVER"){
                max_aOver=std::atof(token.next().c_str());
			} else if(tag=="RREP"){
                rRep=std::atof(token.next().c_str());
			} else if(tag=="SIGMA"){
                sigma=std::atof(token.next().c_str());
			} else if(tag=="SRCUT"){
                srcut=std::atof(token.next().c_str());
			} 
        }

        //=== close the parameter file ====
        fclose(reader);
        reader=NULL;

        //=== check the parameters ====
        if(min_lambda<0) throw std::invalid_argument("Error in fit_cgemm(int,char**): Invalid min lambda.");
        if(min_radius<0) throw std::invalid_argument("Error in fit_cgemm(int,char**): Invalid min radius.");
        if(min_aOver<0) throw std::invalid_argument("Error in fit_cgemm(int,char**): Invalid min aOver.");
        if(max_lambda<0) throw std::invalid_argument("Error in fit_cgemm_omp(int,char**): Invalid max lambda.");
        if(max_radius<0) throw std::invalid_argument("Error in fit_cgemm_omp(int,char**): Invalid max radius.");
        if(max_aOver<0) throw std::invalid_argument("Error in fit_cgemm_omp(int,char**): Invalid max aOver.");
        if(rRep<0) throw std::invalid_argument("Error in fit_cgemm_omp(int,char**): Invalid rRep.");
        if(sigma<=0.0) throw std::invalid_argument("Error in fit_cgemm_omp(int,char**): Invalid sigma.");
        if(srcut<=0.0) throw std::invalid_argument("Error in fit_cgemm_omp(int,char**): Invalid srcut.");
        
        //==== add CGemm cut to engine ====
        std::cout<<"initializing the engine\n";
        const int ntypes=types.size();
        engine.calcs().push_back(
            std::make_shared<CalcCGemmCut>(rc,lambdaC,lambdaS)
        );
        engine.constraints().push_back(
            std::make_shared<ConstraintFreeze>()
        );
        engine.resize(ntypes);
        engine.init();
        static_cast<CalcCGemmCut&>(*engine.calcs().back()).rRep()=rRep;

        //==== make the type map ====
        std::cout<<"making the type map\n";
        std::map<std::string,int> type_map;
        for(int i=0; i<ntypes; ++i){
            type_map[types[i].name()]=i;
            types[i].index()=i;
        }
        
        //==== read cube files ====
        reader=fopen(cubelist.c_str(),"r");
        if(reader==NULL) throw std::runtime_error("Could not open cube list file.");
        while(fgets(input,string::M,reader)!=NULL){
            string::trim_right(input,string::COMMENT);
            Token token(input,string::WS); //split line into tokens
			if(token.end()) continue; //skip empty lines
            cubefiles.push_back(token.next());
            if(!token.end()) qtot.push_back(std::atoi(token.next().c_str()));
            else qtot.push_back(0);
        }
        fclose(reader);
        reader=NULL;
        
        //==== print ====
        int dim=0;
        if(opt_lambda) dim+=2;
        if(opt_radius) dim+=ntypes;
        if(opt_aOver) dim+=ntypes;
        if(opt_rRep) dim+=1;
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
        std::cout<<atom<<"\n";
		std::cout<<print::buf(strbuf)<<"\n";
		std::cout<<print::title("TYPES",strbuf)<<"\n";
        for(int i=0; i<types.size(); ++i){
            std::cout<<types[i]<<"\n";
        }
		std::cout<<print::buf(strbuf)<<"\n";
		std::cout<<print::title("CUBE FILES",strbuf)<<"\n";
        for(int i=0; i<cubefiles.size(); ++i){
            std::cout<<i<<" "<<cubefiles[i]<<" "<<qtot[i]<<"\n";
        }
        std::cout<<print::buf(strbuf)<<"\n";
        std::cout<<engine<<"\n";
        std::cout<<print::title("MD",strbuf)<<"\n";
        std::cout<<"rcut     = "<<rc<<"\n";
        std::cout<<"lambdaC  = "<<lambdaC<<"\n";
        std::cout<<"lambdaS  = "<<lambdaS<<"\n";
        std::cout<<"rrep     = "<<rRep<<"\n";
        std::cout<<"nsteps   = "<<functionData.nsteps<<"\n";
        std::cout<<"ftol     = "<<functionData.ftol<<"\n";
        std::cout<<print::buf(strbuf)<<"\n";
		std::cout<<print::title("OPT",strbuf)<<"\n";
        std::cout<<"dim        = "<<dim<<"\n";
        std::cout<<"tol        = "<<tol<<"\n";
        std::cout<<"miter      = "<<miter<<"\n";
        std::cout<<"opt_lambda = "<<opt_lambda<<"\n";
        std::cout<<"opt_radius = "<<opt_radius<<"\n";
        std::cout<<"opt_aOver  = "<<opt_aOver<<"\n";
        std::cout<<"opt_rRep   = "<<opt_rRep<<"\n";
        std::cout<<"opt_weight = "<<opt_weight<<"\n";
        std::cout<<"min_lambda = "<<min_lambda<<"\n";
        std::cout<<"min_radius = "<<min_radius<<"\n";
        std::cout<<"min_aOver  = "<<min_aOver<<"\n";
        std::cout<<"max_lambda = "<<max_lambda<<"\n";
        std::cout<<"max_radius = "<<max_radius<<"\n";
        std::cout<<"max_aOver  = "<<max_aOver<<"\n";
        std::cout<<"sigma      = "<<sigma<<"\n";
        std::cout<<"srcut      = "<<srcut<<"\n";
        if(algo==nlopt::LN_SBPLX) std::cout<<"algo  = SBPLX\n";
        else if(algo==nlopt::LN_BOBYQA) std::cout<<"algo  = BOBYQA\n";
        else if(algo==nlopt::LN_COBYLA) std::cout<<"algo  = COBYLA\n";
        else if(algo==nlopt::LN_PRAXIS) std::cout<<"algo  = PRAXIS\n";
        else if(algo==nlopt::LN_NELDERMEAD) std::cout<<"algo  = NELDERMEAD\n";
        std::cout<<print::buf(strbuf)<<"\n";
        if(intg!=NULL){
			switch(intg->name()){
				case Integrator::Name::QUICKMIN: std::cout<<static_cast<const Quickmin&>(*intg)<<"\n"; break;
                case Integrator::Name::CG: std::cout<<static_cast<const CG&>(*intg)<<"\n"; break;
				case Integrator::Name::FIRE: std::cout<<static_cast<const Fire&>(*intg)<<"\n"; break;
				case Integrator::Name::VERLET: std::cout<<static_cast<const Verlet&>(*intg)<<"\n"; break;
				case Integrator::Name::VSCALE: std::cout<<static_cast<const VScale&>(*intg)<<"\n"; break;
				case Integrator::Name::BERENDSEN: std::cout<<static_cast<const Berendsen&>(*intg)<<"\n"; break;
				case Integrator::Name::LANGEVIN: std::cout<<static_cast<const Langevin&>(*intg)<<"\n"; break;
				default: std::cout<<"INTEGRATOR = NONE\n";
			}
		} else throw std::invalid_argument("No integrator.");
        std::cout<<print::buf(strbuf)<<"\n";
		
        //==== set up random number generators ====
        std::cout<<"setting random\n";
        std::mt19937 rngen=std::mt19937(std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_real_distribution<double> dist_pos(0.1,1.0);
        
        //==== check the types ====
        std::cout<<"checking the types\n";
        int eIndex=-1;
        for(int i=0; i<types.size(); ++i){
            if(types[i].name()=="X"){
                eIndex=i; break;
            }
            if(types[i].radius()<=0.0) throw std::invalid_argument("Invalid radius.");
            if(types[i].rcut()<=0.0) throw std::invalid_argument("Invalid rcut.");
            if(types[i].aOver()<0.0) throw std::invalid_argument("Invalid overlap amplitude.");
            if(types[i].aRep()<0.0) throw std::invalid_argument("Invalid repulsive amplitude.");
        }
        for(int i=0; i<types.size(); ++i){
            types[i].rcut()*=srcut;
        }
        if(eIndex<0) throw std::invalid_argument("No electron X found in types.");

        //==== read the cube files ====
        std::cout<<"reading the cube files\n";
        data.resize(cubefiles.size());
        for(int i=0; i<cubefiles.size(); ++i){
            CUBE::read(
                cubefiles[i].c_str(),atom,
                data[i].strucA(),
                data[i].espA()
            );
            data[i].espN()=data[i].espA();
        }
        double fpot=0.0;
        if(units::Consts::system()==units::System::AU){
            fpot=1.0;
        } else if(units::Consts::system()==units::System::METAL){
            fpot=units::Eh2Ev;
        }
        for(int n=0; n<data.size(); ++n){
            const Eigen::Vector3i& npts=data[n].espA().n();
            for(int i=0; i<npts[0]; ++i){
                for(int j=0; j<npts[1]; ++j){
                    for(int k=0; k<npts[2]; ++k){
                        const Eigen::Vector3i index=(Eigen::Vector3i()<<i,j,k).finished();
                        data[n].espA().data(index)*=fpot;
                    }
                }
            }
        }

        //==== add electrons ====
        std::cout<<"adding electrons\n";
        for(int i=0; i<cubefiles.size(); ++i){
            Structure& strucA=data[i].strucA();
            Structure& strucN=data[i].strucN();
            const int nCores=strucA.nAtoms();
            const int nShells=nCores-qtot[i];
            const int nTotal=nCores+nShells;
            strucN=Structure(nTotal,strucA.atom());
            //set each nucleus to be a particle with charge +1
            for(int j=0; j<nCores; ++j){
                strucN.name(j)=strucA.name(j);
                strucN.an(j)=strucA.an(j);
                strucN.type(j)=-1;
                strucN.charge(j)=1.0;
                strucN.mass(j)=strucA.mass(j);
                strucN.posn(j)=strucA.posn(j);
                strucN.vel(j)=Eigen::Vector3d::Zero();
                strucN.force(j)=Eigen::Vector3d::Zero();
            }
            //set each electron to be a particle with charge -1
            for(int j=nCores; j<nTotal; ++j){
                strucN.name(j)="X";
                strucN.an(j)=0;
                strucN.type(j)=-1;
                strucN.charge(j)=-1.0;
                strucN.mass(j)=0.1;
                strucN.vel(j)=Eigen::Vector3d::Zero();
                strucN.force(j)=Eigen::Vector3d::Zero();
                strucN.posn(j)=Eigen::Vector3d::Zero();
            }
            //find average posn
            Eigen::Vector3d avg=Eigen::Vector3d::Zero();
            for(int j=0; j<nCores; ++j){
                avg.noalias()+=strucA.posn(j);
            }
            avg*=1.0/nCores;
            double dev=0;
            for(int j=0; j<nCores; ++j){
                dev+=(avg-strucA.posn(j)).squaredNorm();
            }
            dev=std::sqrt(dev/nCores);
            //set electron positions
            if(nShells<nCores){
                //add electrons to atoms
                for(int j=0; j<nShells; ++j){
                    strucN.posn(nCores+j)=strucA.posn(j);
                }
            } else {
                //add electrons to atoms
                for(int j=0; j<nCores; ++j){
                    strucN.posn(nCores+j)=strucA.posn(j);
                }
                //add extra electrons to random posns near the center
                for(int j=nCores; j<nShells; ++j){
                    strucN.posn(nCores+j)=avg;
                    strucN.posn(nCores+j).noalias()+=dev*Eigen::Vector3d::Random();
                }
            }
        }

        //==== set the types ====
        std::cout<<"setting the types\n";
        for(int i=0; i<data.size(); ++i){
            for(int j=0; j<data[i].strucA().nAtoms(); ++j){
                data[i].strucA().type(j)=type_map[data[i].strucA().name(j)];
            }
        }
        for(int i=0; i<data.size(); ++i){
            for(int j=0; j<data[i].strucN().nAtoms(); ++j){
                data[i].strucN().type(j)=type_map[data[i].strucN().name(j)];
            }
        }

        //==== set the mass ====
        std::cout<<"setting the masses\n";
        for(int i=0; i<data.size(); ++i){
            for(int j=0; j<data[i].strucA().nAtoms(); ++j){
                data[i].strucA().mass(j)=types[type_map[data[i].strucA().name(j)]].mass();
            }
        }
        for(int i=0; i<data.size(); ++i){
            for(int j=0; j<data[i].strucN().nAtoms(); ++j){
                data[i].strucN().mass(j)=types[type_map[data[i].strucN().name(j)]].mass();
            }
        }

        /*for(int i=0; i<data.size(); ++i){
            std::cout<<"struc "<<i<<"\n";
            std::cout<<data[i].strucN().R()<<"\n";
            std::cout<<"qtot = "<<qtot[i]<<"\n";
            std::cout<<"natoms = "<<data[i].strucN().nAtoms()<<"\n";
            std::cout<<"name an type charge mass x y z\n";
            for(int j=0; j<data[i].strucN().nAtoms(); ++j){
                std::cout
                    <<data[i].strucN().name(j)<<" "
                    <<data[i].strucN().an(j)<<" "
                    <<data[i].strucN().type(j)<<" "
                    <<data[i].strucN().charge(j)<<" "
                    <<data[i].strucN().mass(j)<<" "
                    <<data[i].strucN().posn(j).transpose()<<" "
                    <<"\n";
            }
        }*/
        
        //==== set up nlopt calculation ====
        std::cout<<"setting up nlopt\n";
        const double MAX_VAL=1.0e2;
        int c;
        std::vector<double> lb(dim),ub(dim),x(dim);
        c=0;
        if(opt_radius) for(int i=0; i<2; ++i) lb[c++]=min_lambda;
        if(opt_radius) for(int i=0; i<ntypes; ++i) lb[c++]=min_radius;
        if(opt_aOver) for(int i=0; i<ntypes; ++i) lb[c++]=min_aOver;
        if(opt_rRep) lb[c++]=min_radius;
        c=0;
        if(opt_lambda) for(int i=0; i<2; ++i) ub[c++]=max_lambda;
        if(opt_radius) for(int i=0; i<ntypes; ++i) ub[c++]=max_radius;
        if(opt_aOver) for(int i=0; i<ntypes; ++i) ub[c++]=max_aOver;
        if(opt_rRep) ub[c++]=max_radius;
        c=0;
        if(opt_lambda) for(int i=0; i<2; ++i) x[c++]=dist_pos(rngen);
        if(opt_radius) for(int i=0; i<ntypes; ++i) x[c++]=dist_pos(rngen);
        if(opt_aOver) for(int i=0; i<ntypes; ++i) x[c++]=dist_pos(rngen);
        if(opt_rRep) x[c++]=dist_pos(rngen);
        
        nlopt::opt opt(algo,(unsigned int)dim);
        opt.set_min_objective(objf,&functionData);
        opt.set_ftol_abs(tol);
        opt.set_maxeval(miter);
		opt.set_lower_bounds(lb);
		opt.set_upper_bounds(ub);

        //==== optimize ====
        std::cout<<"optimizing\n";
		double minf;
        //print head
        std::printf("step error f_avg f_max t_avg t_max ");
        if(opt_lambda) std::printf("lambdaC lambdaS ");
        if(opt_radius) for(int i=0; i<types.size(); ++i) std::printf("R[%s] ",types[i].name().c_str());
        if(opt_aOver) for(int i=0; i<types.size(); ++i) std::printf("A[%s] ",types[i].name().c_str());
        if(opt_rRep) std::printf("rRep ");
        std::printf("\n");
        //optimize
		nlopt::result result=opt.optimize(x,minf);
        //print tail
        std::printf("step error f_avg f_max t_avg t_max ");
        if(opt_lambda) std::printf("lambdaC lambdaS ");
        if(opt_radius) for(int i=0; i<types.size(); ++i) std::printf("R[%s] ",types[i].name().c_str());
        if(opt_aOver) for(int i=0; i<types.size(); ++i) std::printf("A[%s] ",types[i].name().c_str());
        if(opt_rRep) std::printf("rRep ");
        std::printf("\n");
        //print status
		if(result>=0) std::cout<<"optimization successful\n";
		
        //compute the potential
        std::cout<<"computing the potential\n";
        for(int n=0; n<data.size(); ++n){
            Grid& espA=data[n].espA();
            Grid& espN=data[n].espN();
            Structure& strucN=data[n].strucN();
            CalcCGemmCut& calcCGemmCut=static_cast<CalcCGemmCut&>(*engine.calcs().back());
            const Eigen::Vector3i& npts=espA.n();
            const double pf=1.0/(4.0*PI*units::Consts::eps0());
            for(int i=0; i<npts[0]; ++i){
                for(int j=0; j<npts[1]; ++j){
                    for(int k=0; k<npts[2]; ++k){
                        const Eigen::Vector3i indexI=(Eigen::Vector3i()<<i,j,k).finished();
                        const Eigen::Vector3d indexD=(Eigen::Vector3d()<<i,j,k).finished();
                        const Eigen::Vector3d r=espN.origin()+espN.voxel()*indexD;
                        espN.data(indexI)=0.0;
                        for(int m=0; m<strucN.nAtoms(); ++m){
                            const double dr=(r-strucN.posn(m)).norm();
                            const double R=types[strucN.type(m)].radius();
                            const double alpha=calcCGemmCut.lambdaC()/(2.0*R*R);
                            espN.data(indexI)+=pf*strucN.charge(m)/dr*std::erf(std::sqrt(alpha)*dr);
                        }
                    }
                }
            }
            double error=0;
            for(int i=0; i<npts[0]; ++i){
                for(int j=0; j<npts[1]; ++j){
                    for(int k=0; k<npts[2]; ++k){
                        error+=std::fabs(espN(i,j,k)-espA(i,j,k));
                    }
                }
            }
            error/=espA.np();
            std::cout<<"error["<<n<<"] = "<<error<<"\n";
            for(int i=0; i<npts[0]; ++i){
                for(int j=0; j<npts[1]; ++j){
                    for(int k=0; k<npts[2]; ++k){
                        espN(i,j,k)*=1.0/fpot;
                    }
                }
            }
            std::string filename="out_n"+std::to_string(n)+".cube";
            CUBE::write(filename.c_str(),atom,strucN,espN);
        }
    }catch(std::exception& e){
		std::cout<<"ERROR in fit_CGemm::main(int,char**):\n";
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