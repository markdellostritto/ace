// c libraries
#include <cstdlib>
#include <ctime>
//c++
#include <iostream>
#include <vector>
// structure
#include "struc/trajectory.hpp"
#include "struc/group.hpp"
// signal
#include "signal/fft.hpp"
#include "signal/window.hpp"
// format
#include "format/file_traj.hpp"
// str
#include "str/string.hpp"
#include "str/print.hpp"
// chem
#include "chem/units.hpp"

//***********************************************************************
// Main
//***********************************************************************

int main(int argc, char* argv[]){
	//======== local function variables ========
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
	//==== signal variables ====
		double fwidth=0.0;
		double fpmin=0.0,fpmax=0.0;
	//==== atom type ====
		traj.atom().name=true;
		traj.atom().mass=true;
		traj.atom().type=true;
		traj.atom().posn=true;
		traj.atom().vel=true;
		traj.atom().image=true;
	//==== miscellaneous ====
		bool error=false;
	//==== units ====
		units::System unitsys;
	//==== misc ====
		int nprint=0;
		
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
			} else if(tag=="TS"){
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
		std::cout<<"group    = "<<group<<"\n";
		std::cout<<print::buf(strbuf,'*')<<"\n";
		
		//======== check the parameters ========
		if(unitsys==units::System::NONE) throw std::invalid_argument("Invalid unit system.");
		if(fileFormat==FILE_FORMAT::NONE) throw std::invalid_argument("Invalid file format.");
		if(ts<=0) throw std::invalid_argument("Invalid timestep.");
		
		//======== read the simulation ========
		std::cout<<"reading simulation\n";
		read_sim(trajStr,fileFormat,interval,traj.atom(),traj);
		
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
		
		//======== compute the requested properties ========
		std::cout<<"computing properties\n";
		clock_t start=std::clock();
		if(nprint==0) nprint=traj.timesteps()/10;
		if(nprint==0) nprint=1;
		
		//======== set frequency data ========
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
		
		FFT1D fftvx(2*Nt,FFTW_FORWARD); fftvx.init();
		FFT1D fftvy(2*Nt,FFTW_FORWARD); fftvy.init();
		FFT1D fftvz(2*Nt,FFTW_FORWARD); fftvz.init();
		FFT1D fftb(2*Nt,FFTW_BACKWARD); fftb.init();
		FFT1D fftf(2*Nt,FFTW_FORWARD); fftf.init();
		std::function<double(double)> window=window::BlackmanHarris(Nt);
		
		//compute velocity
		std::cout<<"computing velocity\n";
		const double c2[1]={0.5};
		const double c4[2]={-1.0/12.0,2.0/3.0};
		const double c6[3]={1.0/60,-3.0/20.0,3.0/4.0};
		const double c8[4]={-1.0/280.0,4.0/105.0,-1.0/5.0,4.0/5.0};
		const int tb=0;
		const int te=traj.timesteps()-1;
		{//t=tb
			//compute velocity
			for(int n=0; n<group.size(); ++n){
				const int nn=group.atom(n);
				Eigen::Vector3d v=Eigen::Vector3d::Zero();
				v.noalias()+=(traj.frame(tb+1).posn(nn)-traj.frame(tb).posn(nn));
				traj.frame(tb).vel(nn).noalias()=1.0/ts*v;
			}
		}
		{//t=tb+1
			//compute velocity
			for(int n=0; n<group.size(); ++n){
				const int nn=group.atom(n);
				Eigen::Vector3d v=Eigen::Vector3d::Zero();
				for(int i=0; i<1; i++) v.noalias()+=c2[i]*(traj.frame((tb+1)+1-i).posn(nn)-traj.frame((tb+1)-1+i).posn(nn));
				traj.frame(tb+1).vel(nn).noalias()=1.0/ts*v;
			}
		}
		{//t=tb+2
			//compute velocity
			for(int n=0; n<group.size(); ++n){
				const int nn=group.atom(n);
				Eigen::Vector3d v=Eigen::Vector3d::Zero();
				for(int i=0; i<2; i++) v.noalias()+=c4[i]*(traj.frame((tb+2)+2-i).posn(nn)-traj.frame((tb+2)-2+i).posn(nn));
				traj.frame(tb+2).vel(nn).noalias()=1.0/ts*v;
			}
		}
		{//t=tb+3
			//compute velocity
			for(int n=0; n<group.size(); ++n){
				const int nn=group.atom(n);
				Eigen::Vector3d v=Eigen::Vector3d::Zero();
				for(int i=0; i<3; i++) v.noalias()+=c6[i]*(traj.frame((tb+3)+3-i).posn(nn)-traj.frame((tb+3)-3+i).posn(nn));
				traj.frame(tb+3).vel(nn).noalias()=1.0/ts*v;
			}
		}
		//tb+4<=t<=te-4
		for(int t=4; t<traj.timesteps()-4; ++t){
			//compute velocity
			for(int n=0; n<group.size(); ++n){
				const int nn=group.atom(n);
				Eigen::Vector3d v=Eigen::Vector3d::Zero();
				for(int i=0; i<4; i++) v.noalias()+=c8[i]*(traj.frame(t+4-i).posn(nn)-traj.frame(t-4+i).posn(nn));
				traj.frame(t).vel(nn).noalias()=1.0/ts*v;
			}
		}
		{//t=te-3
			//compute velocity
			for(int n=0; n<group.size(); ++n){
				const int nn=group.atom(n);
				Eigen::Vector3d v=Eigen::Vector3d::Zero();
				for(int i=0; i<3; i++) v.noalias()+=c6[i]*(traj.frame((te-3)+3-i).posn(nn)-traj.frame((te-3)-3+i).posn(nn));
				traj.frame(te-3).vel(nn).noalias()=1.0/ts*v;
			}
		}
		{//t=te-2
			//compute velocity
			for(int n=0; n<group.size(); ++n){
				const int nn=group.atom(n);
				Eigen::Vector3d v=Eigen::Vector3d::Zero();
				for(int i=0; i<2; i++) v.noalias()+=c4[i]*(traj.frame((te-2)+2-i).posn(nn)-traj.frame((te-2)-2+i).posn(nn));
				traj.frame(te-2).vel(nn).noalias()=1.0/ts*v;
			}
		}
		{//t=te-1
			//compute velocity
			for(int n=0; n<group.size(); ++n){
				const int nn=group.atom(n);
				Eigen::Vector3d v=Eigen::Vector3d::Zero();
				for(int i=0; i<1; i++) v.noalias()+=c2[i]*(traj.frame((te-1)+1-i).posn(nn)-traj.frame((te-1)-1+i).posn(nn));
				traj.frame(te-1).vel(nn).noalias()=1.0/ts*v;
			}
		}
		{//t=te
			//compute velocity
			for(int n=0; n<group.size(); ++n){
				const int nn=group.atom(n);
				Eigen::Vector3d v=Eigen::Vector3d::Zero();
				v.noalias()+=(traj.frame(te).posn(nn)-traj.frame(te-1).posn(nn));
				traj.frame(te).vel(nn).noalias()=1.0/ts*v;
			}
		}
		
		//compute v2avg
		std::cout<<"computing v2avg\n";
		double v2avg=0.0;
		for(int t=0; t<Nt; ++t){
			double v2avgl=0;
			for(int n=0; n<group.size(); ++n){
				const int nn=group.atom(n);
				v2avgl+=traj.frame(t).vel(nn).squaredNorm();
			}
			v2avg+=v2avgl;
		}
		v2avg/=Nt;
		std::cout<<"v2avg = "<<v2avg<<"\n";
		
		//compute vdos
		std::cout<<"computing vdos\n";
		std::vector<double> vdosr(2*Nt,0);
		std::vector<double> vdosi(2*Nt,0);
		for(int n=0; n<group.size(); ++n){
			const int nn=group.atom(n);
			const double mass=traj.frame(0).mass(nn);
			//velocity transforms
			for(int t=0; t<Nt; ++t){
				fftvx.in(t)[0]=traj.frame(t).vel(nn)[0]; fftvx.in(t)[1]=0.0;
				fftvy.in(t)[0]=traj.frame(t).vel(nn)[1]; fftvy.in(t)[1]=0.0;
				fftvz.in(t)[0]=traj.frame(t).vel(nn)[2]; fftvz.in(t)[1]=0.0;
			}
			for(int t=Nt; t<2*Nt; ++t){
				fftvx.in(t)[0]=0.0; fftvx.in(t)[1]=0.0;
				fftvy.in(t)[0]=0.0; fftvy.in(t)[1]=0.0;
				fftvz.in(t)[0]=0.0; fftvz.in(t)[1]=0.0;
			}
			fftvx.transform();
			fftvy.transform();
			fftvz.transform();
			//power spectrum - freq space
			for(int t=0; t<2*Nt; ++t){
				fftb.in(t)[0]=1.0*(
					fftvx.out(t)[0]*fftvx.out(t)[0]+fftvx.out(t)[1]*fftvx.out(t)[1]+
					fftvy.out(t)[0]*fftvy.out(t)[0]+fftvy.out(t)[1]*fftvy.out(t)[1]+
					fftvz.out(t)[0]*fftvz.out(t)[0]+fftvz.out(t)[1]*fftvz.out(t)[1]
				);
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
			const double norm=1.0/(2.0*Nt)*1.0/(2.0*Nt)*mass*mvv_to_e*2.0/(kb*T)*2.0;
			for(int t=0; t<2*Nt; ++t){
				vdosr[t]+=std::fabs(fftf.out(t)[0])*norm;
				vdosi[t]+=std::fabs(fftf.out(t)[1])*norm;
			}
		}
		
		double real,imag=0;
		for(int t=0; t<2*Nt; ++t){
			real+=vdosr[t];
			imag+=vdosi[t];
		}
		std::cout<<"error - fourier = "<<imag/real<<"\n";
		
		//integrate the vdos - before smoothing (for testing)
		double intg_bs=0;
		intg_bs+=0.5*vdosr[0]*df;
		for(int t=1; t<2*Nt-1; ++t){
			intg_bs+=vdosr[t]*df;
		}
		intg_bs+=0.5*vdosr[2*Nt-1]*df;
		std::cout<<"integral - before smoothing = "<<intg_bs<<"\n";
		
		//smooth the vdos
		signala::smooth(vdosr,fwidth);
		
		//integrate the vdos - after smoothing
		double intg_as=0;
		intg_as+=0.5*vdosr[0]*df;
		for(int t=1; t<2*Nt-1; ++t){
			intg_as+=vdosr[t]*df;
		}
		intg_as+=0.5*vdosr[2*Nt-1]*df;
		std::cout<<"integral - after smoothing = "<<intg_as<<"\n";
		std::cout<<"error (%) - smoothing = "<<(intg_as-intg_bs)/intg_bs*100.0<<"\n";
		
		//write the vdos
		FILE* writer=fopen("vdos.dat","w");
		if(writer!=NULL){
			//fprintf(writer,"#freq(THz) vdosr vdosi\n");
			fprintf(writer,"#freq(THz) vdos\n");
			for(int t=imin; t<imax; ++t){
				//fprintf(writer,"%12.8f %12.8f %12.8f\n",t*df,vdosr[t],vdosi[t]);
				fprintf(writer,"%12.8f %12.8f\n",t*df,vdosr[t]);
			}
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
