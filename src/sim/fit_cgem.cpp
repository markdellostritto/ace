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
#include "sim/calc_cgem_cut.hpp"
// format
#include "format/cube_struc.hpp"
//sim
#include "sim/fit_cgem.hpp"
#include "sim/constraint_freeze.hpp"
//nlopt
#include <nlopt.hpp>

using math::constants::PI;

int step=0;
int count=0;

bool opt_lambda=false;
bool opt_radius=false;
bool opt_aOver=false;

//*****************************************************
// CGemType
//*****************************************************

//==== constructors/destructors ====

CGemType::CGemType(){
    name_="NONE";
    radius_=0;
    aOver_=0;
    aRep_=0;
    index_=-1;    
}

//==== operators ====

std::ostream& operator<<(std::ostream& out, const CGemType& type){
    return out<<type.name()<<" "<<type.index()<<" "<<type.radius()<<" "<<type.aOver()<<" "<<type.aRep();
}

//==== member functions ====

void CGemType::read(Token& token){
    name_=token.next();
    mass_=std::atof(token.next().c_str());
    radius_=std::atof(token.next().c_str());
    aOver_=std::atof(token.next().c_str());
    aRep_=std::atof(token.next().c_str());
    if(radius_<0) throw std::invalid_argument("CGemType::read(Token&): invalid mass.");
    if(radius_<0) throw std::invalid_argument("CGemType::read(Token&): invalid radius.");
    if(aOver_<0) throw std::invalid_argument("CGemType::read(Token&): invalid overlap amplitude.");
    if(aRep_<0) throw std::invalid_argument("CGemType::read(Token&): invalid repulsive amplitude.");
}

//*****************************************************
// Objective Function
//*****************************************************

double objf(const std::vector<double> &x, std::vector<double> &grad, void* fData){
    FunctionData& fData_=*((FunctionData*)fData);
    const int ntypes=fData_.types.size();
    std::cout<<"x = "; for(int i=0; i<x.size(); ++i) std::cout<<x[i]<<" "; std::cout<<"\n";
    double error=0;
    int c;
    
    //==== unpack the parameters ====
    CalcCGemCut& calcCGemCut=static_cast<CalcCGemCut&>(*fData_.engine.calcs().front());
    calcCGemCut.lambdaC()=1.0;
    calcCGemCut.lambdaS()=1.0;
    c=0;
    if(opt_lambda){
        calcCGemCut.lambdaC()=x[c++];
        calcCGemCut.lambdaS()=x[c++];
    }
    if(opt_radius) for(int i=0; i<ntypes; ++i) calcCGemCut.radius()[i]=x[c++];
    else for(int i=0; i<ntypes; ++i) calcCGemCut.radius()[i]=fData_.types[i].radius();
    if(opt_aOver) for(int i=0; i<ntypes; ++i) calcCGemCut.aOver()(i,i)=x[c++];
    else for(int i=0; i<ntypes; ++i) calcCGemCut.aOver()(i,i)=fData_.types[i].aOver();
    for(int i=0; i<ntypes; ++i) calcCGemCut.aRep()(i,i)=fData_.types[i].aRep();
    calcCGemCut.init();
    c=0;
    if(opt_lambda){
        c++;
        c++;
    }
    if(opt_radius) for(int i=0; i<ntypes; ++i) fData_.types[i].radius()=x[c++];
    if(opt_aOver) for(int i=0; i<ntypes; ++i) fData_.types[i].aOver()=x[c++];
    //std::cout<<calcCGemCut<<"\n";
    //std::cout<<"radius = "<<calcCGemCut.radius()<<"\n";
    //std::cout<<"aOver  = "<<calcCGemCut.aOver()<<"\n";
    //std::cout<<"aRep   = "<<calcCGemCut.aRep()<<"\n";

    //==== loop over each structure ====
    Engine& engine=fData_.engine;
    double np=0;
    double t_avg=0;
    double f_avg=0;
    int t_max=0;
    for(int n=0; n<fData_.data.size(); ++n){
        Structure& strucA=fData_.data[n].strucA();
        Structure& strucN=fData_.data[n].strucN();
        Grid& espA=fData_.data[n].espA();
        Grid& espN=fData_.data[n].espN();

        //set the constraints
        std::vector<int> indices;
        for(int m=0; m<strucN.nAtoms(); ++m){
            if(strucN.name(m)!="X"){
                indices.push_back(m);
            }
        }
        engine.constraints().front()->indices()=indices;
        //for(int i=0; i<indices.size(); ++i) std::cout<<"indices["<<i<<"] = "<<indices[i]<<"\n";

        //reset the electron positions
        const int nAtoms=strucA.nAtoms();
        for(int i=0; i<nAtoms; ++i){
            strucN.posn(i+nAtoms)=strucA.posn(i);
        }

        //relax the system
        if(PRINT_STRUC>0){
            std::cout<<strucN.nAtoms()<<"\nmolecule\n";
            for(int m=0; m<strucN.nAtoms(); ++m){
                std::cout<<strucN.name(m)
                    <<" "<<strucN.type(m)
                    <<" "<<strucN.charge(m)
                    <<" "<<fData_.types[strucN.type(m)].radius()
                    <<" "<<fData_.types[strucN.type(m)].aOver()
                    <<" "<<fData_.types[strucN.type(m)].aRep()
                    <<" "<<strucN.posn(m).transpose()
                    <<"\n";
            }
        }
        int t=0;
        double ftot=0;
        for(t=0; t<fData_.nsteps; ++t){
            //compute step
            fData_.intg->compute(strucN,engine);
            //compute total force
            ftot=0;
            for(int m=0; m<strucN.nAtoms(); ++m){
                ftot+=strucN.force(m).squaredNorm();
            }
            ftot=std::sqrt(ftot/strucN.nAtoms());
            if(ftot<fData_.ftol) break;
        }
        t_avg+=t;
        f_avg+=ftot;
        if(t>t_max) t_max=t;
        if(PRINT_STRUC>0){
            std::cout<<strucN.nAtoms()<<"\nmolecule\n";
            for(int m=0; m<strucN.nAtoms(); ++m){
                std::cout<<strucN.name(m)
                    <<" "<<strucN.type(m)
                    <<" "<<strucN.charge(m)
                    <<" "<<fData_.types[strucN.type(m)].radius()
                    <<" "<<fData_.types[strucN.type(m)].aOver()
                    <<" "<<fData_.types[strucN.type(m)].aRep()
                    <<" "<<strucN.posn(m).transpose()
                    <<"\n";
            }
        }

        //compute the potential
        const Eigen::Vector3i& npts=espA.n();
        const double pf=1.0/(4.0*PI*units::Consts::eps0());
        for(int i=0; i<npts[0]; ++i){
            for(int j=0; j<npts[1]; ++j){
                for(int k=0; k<npts[2]; ++k){
                    const Eigen::Vector3d indexD=(Eigen::Vector3d()<<i,j,k).finished();
                    const Eigen::Vector3d r=espA.origin()+espA.voxel()*indexD;
                    //std::cout<<"indexI = "<<indexI.transpose()<<"\n";
                    //std::cout<<"indexD = "<<indexD.transpose()<<"\n";
                    //std::cout<<"r = "<<r.transpose()<<"\n";
                    double& data=espN.data(i,j,k);
                    data=0.0;
                    for(int m=0; m<strucN.nAtoms(); ++m){
                        const double dr=(r-strucN.posn(m)).norm();
                        const double R=fData_.types[strucN.type(m)].radius();
                        const double alpha=calcCGemCut.lambdaC()/(2.0*R*R);
                        data+=pf*strucN.charge(m)/dr*std::erf(std::sqrt(alpha)*dr);
                        //std::cout<<strucN.name(m)<<" "<<strucN.type(m)<<" "<<strucN.charge(m)<<" "<<strucN.posn(m).transpose()<<" "<<pf*strucN.charge(m)/dr*std::erf(std::sqrt(alpha)*dr)<<"\n";
                    }
                    //std::cout<<"espN = "<<espN.data(indexI)<<"\n";
                    //std::cout<<"espA = "<<espA.data(indexI)<<"\n";
                }
            }
        }

        //compute the error
        for(int i=0; i<npts[0]; ++i){
            for(int j=0; j<npts[1]; ++j){
                for(int k=0; k<npts[2]; ++k){
                    //std::cout<<"espA = "<<espA.data(index)<<"\n";
                    //std::cout<<"espN = "<<espN.data(index)<<"\n";
                    error+=std::fabs(espA.data(i,j,k)-espN.data(i,j,k));
                }
            }
        }
        np+=espA.np();
    }
    error/=np;
    t_avg/=fData_.data.size();
    f_avg/=fData_.data.size();
    
    //==== print error ====
    std::cout<<"step "<<step<<" error "<<error<<" f_avg "<<f_avg<<" t_avg "<<t_avg<<" t_max "<<t_max<<"\n";
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
        atom.name=true; atom.type=true; atom.an=true;
        atom.charge=true; atom.mass=true; 
        atom.posn=true; atom.force=true; atom.vel=true;
    //cube files
        std::string cubelist;
        std::vector<std::string> cubefiles;
        std::vector<int> qtot;
    //fit
        double rc=0;//cutoff
        double lambdaC=1.0;//initial guess for lambdaC
        double lambdaS=1.0;//initial guess for lambdaS
        double tol=0.0;
        int miter=0;
        nlopt::algorithm algo;
    //function data
        FunctionData functionData;
        std::vector<CubeData>& data=functionData.data;
        Engine& engine=functionData.engine;
        std::vector<CGemType>& types=functionData.types;
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
                types.push_back(CGemType());
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
			} 
        }

        //=== close the parameter file ====
        fclose(reader);
        reader=NULL;
        
        //==== add cgem cut to engine ====
        std::cout<<"initializing the engine\n";
        const int ntypes=types.size();
        engine.calcs().push_back(
            std::make_shared<CalcCGemCut>(rc,lambdaC,lambdaS)
        );
        engine.constraints().push_back(
            std::make_shared<ConstraintFreeze>()
        );
        engine.resize(ntypes);
        engine.init();

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
            std::cout<<cubefiles[i]<<" "<<qtot[i]<<"\n";
        }
        std::cout<<print::buf(strbuf)<<"\n";
        std::cout<<engine<<"\n";
        std::cout<<print::title("MD",strbuf)<<"\n";
        std::cout<<"rcut     = "<<rc<<"\n";
        std::cout<<"lambdaC  = "<<lambdaC<<"\n";
        std::cout<<"lambdaS  = "<<lambdaS<<"\n";
        std::cout<<"nsteps   = "<<functionData.nsteps<<"\n";
        std::cout<<"ftol     = "<<functionData.ftol<<"\n";
        std::cout<<print::buf(strbuf)<<"\n";
		std::cout<<print::title("OPT",strbuf)<<"\n";
        std::cout<<"dim   = "<<dim<<"\n";
        std::cout<<"tol   = "<<tol<<"\n";
        std::cout<<"miter = "<<miter<<"\n";
        std::cout<<"opt_lambda = "<<opt_lambda<<"\n";
        std::cout<<"opt_radius = "<<opt_radius<<"\n";
        std::cout<<"opt_aOver    = "<<opt_aOver<<"\n";
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
        std::uniform_real_distribution<double> dist_uni(-1.0,1.0);

        //==== check the types ====
        std::cout<<"checking the types\n";
        int eIndex=-1;
        for(int i=0; i<types.size(); ++i){
            if(types[i].name()=="X"){
                eIndex=i; break;
            }
            if(types[i].radius()<=0.0) throw std::invalid_argument("Invalid radius.");
            if(types[i].aOver()<0.0) throw std::invalid_argument("Invalid overlap amplitude.");
            if(types[i].aRep()<0.0) throw std::invalid_argument("Invalid repulsive amplitude.");
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
            const int natoms=strucA.nAtoms();
            const int nelectrons=natoms-qtot[i];
            const int ntot=natoms+nelectrons;
            strucN=Structure(ntot,strucA.atom());
            //set each nucleus to be a particle with charge +1
            for(int j=0; j<natoms; ++j){
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
            for(int j=natoms; j<ntot; ++j){
                strucN.name(j)="X";
                strucN.an(j)=0;
                strucN.type(j)=-1;
                strucN.charge(j)=-1.0;
                strucN.mass(j)=0.1;
                strucN.vel(j)=Eigen::Vector3d::Zero();
                strucN.force(j)=Eigen::Vector3d::Zero();
                strucN.posn(j)=Eigen::Vector3d::Zero();
            }
            //set electron positions
            Eigen::Vector3d avg=Eigen::Vector3d::Zero();
            for(int j=0; j<natoms; ++j){
                avg.noalias()+=strucA.posn(j);
            }
            avg*=1.0/natoms;
            double dev=0;
            for(int j=0; j<natoms; ++j){
                dev+=(avg-strucA.posn(j)).squaredNorm();
            }
            dev=std::sqrt(dev/natoms);
            //set electron positions
            if(nelectrons<natoms){
                //add electrons to atoms
                for(int j=0; j<nelectrons; ++j){
                    strucN.posn(natoms+j)=strucA.posn(j);
                }
            } else {
                //add electrons to atoms
                for(int j=0; j<natoms; ++j){
                    strucN.posn(natoms+j)=strucA.posn(j);
                }
                //add extra electrons to random posns near the center
                for(int j=natoms; j<nelectrons; ++j){
                    strucN.posn(natoms+j)=avg;
                    strucN.posn(natoms+j).noalias()+=dev*Eigen::Vector3d::Random();
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
        if(opt_lambda){
            lb[c++]=1.0e-3;
            lb[c++]=1.0e-3;
        }
        if(opt_radius) for(int i=0; i<ntypes; ++i) lb[c++]=1.0e-6;
        if(opt_aOver) for(int i=0; i<ntypes; ++i) lb[c++]=1.0e-6;
        //if(opt_aOver) for(int i=0; i<ntypes; ++i) lb[c++]=-HUGE_VAL;
        c=0;
        if(opt_lambda){
            ub[c++]=HUGE_VAL;
            ub[c++]=HUGE_VAL;
        }
        if(opt_radius) for(int i=0; i<ntypes; ++i) ub[c++]=HUGE_VAL;
        if(opt_aOver) for(int i=0; i<ntypes; ++i) ub[c++]=HUGE_VAL;
        c=0;
        if(opt_lambda){
            x[c++]=dist_pos(rngen);
            x[c++]=dist_pos(rngen);
        }
        if(opt_radius) for(int i=0; i<ntypes; ++i) x[c++]=dist_pos(rngen);
        if(opt_aOver) for(int i=0; i<ntypes; ++i) x[c++]=dist_pos(rngen);
        //if(opt_aOver) for(int i=0; i<ntypes; ++i) x[c++]=dist_uni(rngen);
        
        nlopt::opt opt(algo,(unsigned int)dim);
        opt.set_min_objective(objf,&functionData);
        opt.set_ftol_abs(tol);
        opt.set_maxeval(miter);
		opt.set_lower_bounds(lb);
		opt.set_upper_bounds(ub);

        //==== optimize ====
        std::cout<<"optimizing\n";
		double minf;
		nlopt::result result=opt.optimize(x,minf);
		if(result>=0) std::cout<<"optimization successful\n";
		
        //compute the potential
        std::cout<<"computing the potential\n";
        for(int n=0; n<data.size(); ++n){
            Grid& espA=data[n].espA();
            Grid& espN=data[n].espN();
            Structure& strucN=data[n].strucN();
            CalcCGemCut& calcCGemCut=static_cast<CalcCGemCut&>(*engine.calcs().back());
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
                            const double alpha=calcCGemCut.lambdaC()/(2.0*R*R);
                            espN.data(indexI)+=pf*strucN.charge(m)/dr*std::erf(std::sqrt(alpha)*dr);
                        }
                    }
                }
            }
            for(int i=0; i<npts[0]; ++i){
                for(int j=0; j<npts[1]; ++j){
                    for(int k=0; k<npts[2]; ++k){
                        espN(i,j,k)*=1.0/fpot;
                        //espA(i,j,k)*=1.0/fpot;
                    }
                }
            }
            std::string filename="out_n"+std::to_string(n)+".cube";
            CUBE::write(filename.c_str(),atom,strucN,espN);
            //CUBE::write(filename.c_str(),atom,strucA,espA);
        }
    }catch(std::exception& e){
		std::cout<<"ERROR in fit_cgem::main(int,char**):\n";
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