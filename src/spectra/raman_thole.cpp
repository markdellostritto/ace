// c++
#include <map>
// structure
#include "struc/trajectory.hpp"
#include "struc/interval.hpp"
#include "struc/group.hpp"
// chem
#include "chem/units.hpp"
// format
#include "format/format.hpp"
#include "format/file_traj.hpp"
// str
#include "str/print.hpp"
// sim
#include "sim/calc.hpp"
#include "sim/calc_thole_long.hpp"
#include "sim/thole.hpp"
// math
#include "math/const.hpp"
// signal
#include "signal/fft.hpp"
#include "signal/window.hpp"

int main(int argc, char* argv[]){
	//==== typdefs ====
		typedef fourier::FFT<1,fourier::DataT::COMPLEX,fourier::DataT::COMPLEX> FFT1D;
    //==== file i/o ====
		FILE* reader=NULL;
		FILE_FORMAT::type fileFormat;
		char* paramFile=new char[string::M];
		char* input    =new char[string::M];
		char* trajStr  =new char[string::M];
		char* strbuf   =new char[print::len_buf];
	//==== simulation variables ====
		Trajectory traj;
		Interval interval;
		double ts=0;
		double T=0;
		Group group;
		int ntypes=0;
		std::map<std::string,int> types;
	//==== signal variables ====
		double fvis=600.0;
		double fwidth=0.0;
		double fpmin=0.0,fpmax=0.0;
	//==== atom type ====
		traj.atom().an=true;
		traj.atom().name=true;
		traj.atom().mass=true;
		traj.atom().type=true;
		traj.atom().posn=true;
		traj.atom().dipole=true;
		traj.atom().image=true;
	//==== thole ====
		Thole thole;
		double a=0;
		double rc=0;
		double prec=0;
		CalcTholeLong::IDD idd=CalcTholeLong::IDD::NONE;
	//==== raman ====
		std::vector<Eigen::Matrix3d> alpha;
	//==== units ====
		units::System unitsys;
	//==== misc ====
		int nprint=0;
		bool error=false;
	
	try{
		if(argc!=2) throw std::invalid_argument("Invalid number of command-line arguments.");
		
		//======== copy the parameter file ========
		std::cout<<"reading parameter file\n";
		std::strcpy(paramFile,argv[1]);
		
		//======== read in the general parameters ========
		reader=fopen(paramFile,"r");
		if(reader==NULL) throw std::runtime_error("I/O Error: could not open parameter file.");
		std::cout<<"reading general parameters\n";
		while(fgets(input,string::M,reader)!=NULL){
			string::trim_right(input,string::COMMENT);
			Token token(input,string::WS);
			if(token.end()) continue;
			const std::string tag=string::to_upper(token.next());
			if(tag=="INTERVAL"){
				interval=Interval::read(token.next().c_str(),interval);
			} else if(tag=="TRAJ"){
				std::strcpy(trajStr,token.next().c_str());
			} else if(tag=="FORMAT"){
				fileFormat=FILE_FORMAT::read(string::to_upper(token.next()).c_str());
			} else if(tag=="UNITS"){
				unitsys=units::System::read(string::to_upper(token.next()).c_str());
			} else if(tag=="NPRINT"){
				nprint=std::atoi(token.next().c_str());
			} else if(tag=="TS" || tag=="TIMESTEP"){
				ts=std::atof(token.next().c_str());
			} else if(tag=="TEMP"){
				T=std::atof(token.next().c_str());
			} else if(tag=="FWIDTH"){
				fwidth=std::atof(token.next().c_str());
			} else if(tag=="GROUP"){
				group.read(token);
			} else if(tag=="FINT"){
				fpmin=std::atof(token.next().c_str());
				fpmax=std::atof(token.next().c_str());
			} else if(tag=="INTER"){
				thole.inter()=string::boolean(token.next().c_str());
			} else if(tag=="ALPHAR"){
				thole.alphar()=string::boolean(token.next().c_str());
			} else if(tag=="CALC"){
				const Calculator::Name name=Calculator::Name::read(string::to_upper(token.next()).c_str());
				if(name!=Calculator::Name::THOLE_LONG) throw std::invalid_argument("Invalid calculator.");
				thole.calc().read(token);
			} else if(tag=="TYPE"){
				std::string name=token.next();
				const double alpha=std::atof(token.next().c_str());
				types[name]=ntypes++;
				if(ntypes>thole.alpha().size()) thole.alpha().resize(ntypes);
				thole.alpha()[ntypes-1]=alpha;
			} else if(tag=="A"){
				a=std::atof(token.next().c_str());
			} else if(tag=="RC" || tag=="RCUT"){
				rc=std::atof(token.next().c_str());
				if(rc<=0) throw std::invalid_argument("Invalid rc.");
			} else if(tag=="PREC"){
				prec=std::atof(token.next().c_str());
				if(prec<=0) throw std::invalid_argument("Invalid prec.");
			} else if(tag=="IDD"){
				idd=CalcTholeLong::IDD::read(string::to_upper(token.next()).c_str());
			} else if(tag=="FVIS"){
				fvis=std::atof(token.next().c_str());
			}
		}
		//close the file
		fclose(reader);
		reader=NULL;
		
		//======== initialize the unit system ========
		std::cout<<"initializing the unit system\n";
		units::Consts::init(unitsys);
		const double hbar=units::Consts::hbar();
		const double kb=units::Consts::kb();
		const double mvv_to_e=units::Consts::mvv_to_e();
		
		//======== set calculator ========
		std::cout<<"setting calculator\n";
		thole.calc()=CalcTholeLong(rc);
		thole.calc().prec()=prec;
		thole.calc().idd()=idd;
		thole.calc().a()=a;
		thole.calc().init();
		thole.calc().resize(ntypes);
		for (std::map<std::string,int>::iterator it = types.begin(); it != types.end(); it++){
			thole.calc().alpha()[it->second]=thole.alpha()[it->second];
		}

		//======== print the parameters ========
		std::cout<<print::buf(strbuf,'*')<<"\n";
		std::cout<<print::title("PARAMETERS",strbuf)<<"\n";
		std::cout<<"units    = "<<unitsys<<"\n";
		std::cout<<"hbar     = "<<hbar<<"\n";
		std::cout<<"kb       = "<<kb<<"\n";
		std::cout<<"nprint   = "<<nprint<<"\n";
		std::cout<<"atomt    = "<<traj.atom()<<"\n";
		std::cout<<"traj     = \""<<trajStr<<"\"\n";
		std::cout<<"format   = "<<fileFormat<<"\n";
		std::cout<<"interval = "<<interval<<"\n";
		std::cout<<"timestep = "<<ts<<"\n";
		std::cout<<"temp     = "<<T<<"\n";
		std::cout<<"fwidth   = "<<fwidth<<"\n";
		std::cout<<"fint     = "<<fpmin<<" "<<fpmax<<"\n";
		std::cout<<"fvis     = "<<fvis<<"\n";
		std::cout<<"group    = "<<group<<"\n";
		std::cout<<print::buf(strbuf,'*')<<"\n";
		std::cout<<print::title("TYPES",strbuf)<<"\n";
		for (const auto& [key, value] : types){
        	std::cout<<key<<" = "<<value<<"\n";
		}
		std::cout<<print::buf(strbuf)<<"\n";
		std::cout<<thole<<"\n";
		std::cout<<print::buf(strbuf)<<"\n";
		
		//======== check the parameters ========
		if(unitsys==units::System::NONE) throw std::invalid_argument("Invalid unit system.");
		if(fileFormat==FILE_FORMAT::NONE) throw std::invalid_argument("Invalid file format.");
		if(ts<=0) throw std::invalid_argument("Invalid timestep.");
		
		//======== read the simulation ========
		std::cout<<"reading simulation\n";
		read_sim(trajStr,fileFormat,interval,traj.atom(),traj);
		
		//======== set the types ========
		std::cout<<"setting the types\n";
		for(int t=0; t<traj.timesteps(); ++t){
			for(int n=0; n<traj.frame(t).nAtoms(); ++n){
				traj.frame(t).type(n)=types.at(traj.frame(t).name(n));
			}
		}
		for(int t=0; t<traj.timesteps(); ++t){
			for(int n=0; n<traj.frame(t).nAtoms(); ++n){
				traj.frame(t).dipole(n)=Eigen::Vector3d::Random()/math::constants::Rad3/1000.0;
			}
		}

		//======== unwrap the coordinates ========
		std::cout<<"unwrapping the coordinates\n";
		Trajectory::set_image(traj);
		Trajectory::unwrap(traj);
		
		//======== print the data to screen ========
		std::cout<<"SIMULATION = \n"<<traj<<"\n";
		std::cout<<traj.frame(0)<<"\n";
		
		//======== build group ========
		std::cout<<"building group\n";
		group.build(traj.frame(0));
		std::cout<<"group size: "<<group.size()<<"\n";
		
		//======== set frequency data ========
		std::cout<<"setting frequency data\n";
		//printing
		if(nprint==0) nprint=traj.timesteps()/10;
		if(nprint==0) nprint=1;
		//number of steps
		const int Nt=traj.timesteps();
		//frequency - min/max/step
		const double fmin=1000.0/(ts*2.0*Nt);
		const double fmax=1000.0/ts;
		const double df=fmin;
		//frequency - min/max - index
		const int imin=(int)std::floor(fpmin/df);
		const int imax=(int)std::ceil(fpmax/df);
		std::cout<<"fint = ["<<fmin<<","<<fmax<<","<<df<<"]\n";
		std::cout<<"pint = ["<<imin<<","<<imax<<","<<1<<"]\n";
		
		//======== compute polarizability ========
		std::cout<<"computing polarizability\n";
		thole.calc().init(traj.frame(0));
		alpha.resize(Nt);
		for(int t=0; t<Nt; ++t){
			if(t%nprint==0) std::cout<<"Timestep "<<t<<"\n";
			thole.compute(traj.frame(t),alpha[t]);
			alpha[t]/=traj.frame(t).nAtoms();
			//std::cout<<"alpha["<<t<<"] = \n"<<alpha[t]<<"\n";
		}
		
		//======== compute isotropic polarizability ========
		std::cout<<"computing polarizability - isotropic\n";
		std::vector<double> alphaI(Nt);
		for(int t=0; t<Nt; ++t){
			alphaI[t]=1.0/3.0*alpha[t].trace();
		}
		std::vector<double> alphaID(Nt);
		alphaID[0]=(alphaI[1]-alphaI[0])/ts;
		for(int t=1; t<Nt-1; ++t){
			alphaID[t]=(alphaI[t+1]-alphaI[t-1])/(2.0*ts);
		}
		alphaID[Nt-1]=(alphaI[Nt-1]-alphaI[Nt-2])/ts;

		//======== compute raman spectrum ========
		std::cout<<"computing raman spectrum - isotropic\n";
		FFT1D fftRI(2*Nt,FFTW_FORWARD); fftRI.init();
		FFT1D fftb(2*Nt,FFTW_BACKWARD); fftb.init();
		FFT1D fftf(2*Nt,FFTW_FORWARD); fftf.init();
		std::function<double(double)> window=window::BlackmanHarris(Nt);
		std::vector<double> ramanI(2*Nt,0);
		//compute the forward transform
		for(int t=0; t<Nt; ++t){
			fftRI.in(t)[0]=alphaID[t];
			fftRI.in(t)[1]=0.0;
		}
		for(int t=Nt; t<2*Nt; ++t){
			fftRI.in(t)[0]=0.0;
			fftRI.in(t)[1]=0.0;
		}
		fftRI.transform();
		//compute power spectrum - freq space
		for(int t=0; t<2*Nt; ++t){
			fftb.in(t)[0]=fftRI.out(t)[0]*fftRI.out(t)[0]+fftRI.out(t)[1]*fftRI.out(t)[1];
			fftb.in(t)[1]=0.0;
		}
		//perform the reverse transform
		fftb.transform();
		//normalize, shift, and window the correlation function
		for(int t=0; t<Nt; ++t){
			const double fac=window(t)/(t+1);
			fftf.in(t)[0]=fftb.out(t+Nt)[0]*fac;
			fftf.in(t)[1]=fftb.out(t+Nt)[1]*fac;
		}
		for(int t=Nt; t<2*Nt; ++t){
			const double fac=window(t-Nt)/(2*Nt-t);
			fftf.in(t)[0]=fftb.out(t-Nt)[0]*fac;
			fftf.in(t)[1]=fftb.out(t-Nt)[1]*fac;
		}
		//perform a forward transform 
		fftf.transform();
		//copy the data into the vdos array
		const double norm=1.0/(2.0*Nt);
		for(int t=0; t<2*Nt; ++t){
			const double w=fmin*(t+0.5);
			ramanI[t]=std::fabs(fftf.out(t)[0])*norm/pow(w/fvis-1.0,4.0);
		}

		//======== smooth raman spectrum ========
		if(fwidth>0.0) signala::smooth(ramanI,fwidth);
		
		//======== write raman spectrum ========
		FILE* writer=fopen("raman.dat","w");
		if(writer!=NULL){
			fprintf(writer,"#freq(THz) ramanI ramanA\n");
			for(int t=imin; t<imax; ++t){
				//fprintf(writer,"%12.8f %12.8f %12.8f\n",t*df,vdosr[t],vdosi[t]);
				fprintf(writer,"%12.8f %12.8f\n",t*df,ramanI[t]);
			}
			fclose(writer);
			writer=NULL;
		}
		
	}catch(std::exception& e){
		std::cout<<e.what()<<"\n";
		std::cout<<"ANALYSIS FAILED.\n";
		error=true;
	}
	
	std::cout<<"freeing local variables\n";
	delete[] paramFile;
	delete[] input;
	delete[] trajStr;
	delete[] strbuf;
	
	std::cout<<"exiting program\n";
	if(error) return 1;
	else return 0;
}