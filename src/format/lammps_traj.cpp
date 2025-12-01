// c libraries
#include <ctime>
// c++ libraries
#include <iostream>
// str
#include "str/string.hpp"
#include "str/token.hpp"
// chem
#include "chem/ptable.hpp"
#include "chem/units.hpp"
// structure
#include "format/lammps_traj.hpp"

namespace LAMMPS{

//*****************************************************
//STYLE_ATOM struct
//*****************************************************

STYLE_ATOM::type STYLE_ATOM::read(const char* str){
	if(std::strcmp(str,"FULL")==0) return STYLE_ATOM::FULL;
	else if(std::strcmp(str,"BOND")==0) return STYLE_ATOM::BOND;
	else if(std::strcmp(str,"ATOMIC")==0) return STYLE_ATOM::ATOMIC;
	else if(std::strcmp(str,"CHARGE")==0) return STYLE_ATOM::CHARGE;
	else return STYLE_ATOM::UNKNOWN;
}

std::ostream& operator<<(std::ostream& out, STYLE_ATOM::type& t){
	if(t==STYLE_ATOM::FULL) out<<"FULL";
	else if(t==STYLE_ATOM::BOND) out<<"BOND";
	else if(t==STYLE_ATOM::ATOMIC) out<<"ATOMIC";
	else if(t==STYLE_ATOM::CHARGE) out<<"CHARGE";
	else out<<"UNKNOWN";
	return out;
}

//*****************************************************
//FORMAT_ATOM struct
//*****************************************************

std::ostream& operator<<(std::ostream& out, FORMAT_ATOM& f){
	out<<"index = "<<f.index<<"\n";
	out<<"mol   = "<<f.mol<<"\n";
	out<<"type  = "<<f.type<<"\n";
	out<<"x     = "<<f.x<<"\n";
	out<<"y     = "<<f.y<<"\n";
	out<<"z     = "<<f.z<<"\n";
	out<<"xu    = "<<f.xu<<"\n";
	out<<"yu    = "<<f.yu<<"\n";
	out<<"zu    = "<<f.zu<<"\n";
	out<<"q     = "<<f.q<<"\n";
	out<<"fx    = "<<f.fx<<"\n";
	out<<"fy    = "<<f.fy<<"\n";
	out<<"fz    = "<<f.fz;
	return out;
}

namespace DUMP{

void read(const char* file, const Interval& interval, const Atom& atom, Trajectory& traj){
	static const char* funcName="read(const char*,const Interval&,Atom&,Trajectory&,Format&)";
	if(LAMMPS_PRINT_FUNC>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<funcName<<":\n";
	//======== local function variables ========
	//==== file i/o ====
		FILE* reader=NULL;
		char* input=new char[string::M];
		char* temp=new char[string::M];
		Token token;
	//==== time info ====
		int ts=0;//number of timesteps
	//==== cell info ====
		Eigen::Matrix3d lv=Eigen::Matrix3d::Zero();
	//==== atom info ====
		int natoms=0;
		DATA_ATOM dataAtom;
		FORMAT_ATOM formatAtom;
		int minindex=-1;
	//==== timing ====
		clock_t start,stop;
		double time;
	//==== misc ====
		bool error=false;
	//==== units ====
		double s_len=1.0,s_energy=1.0,s_mass=1.0;
		
	try{
		//==== start the timer ====
		start=std::clock();
		
		//==== open the file ====
		if(LAMMPS_PRINT_STATUS>0) std::cout<<"opening file\n";
		reader=fopen(file,"r");
		if(reader==NULL) throw std::runtime_error("I/O Error: Could not open file.");
		
		//==== find the number of timesteps ====
		if(LAMMPS_PRINT_STATUS>0) std::cout<<"reading number of timesteps\n";
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strstr(input,"TIMESTEP")!=NULL) ++ts;
		}
		std::rewind(reader);
		
		//==== read in the atom format ====
		if(LAMMPS_PRINT_STATUS>0) std::cout<<"reading atom format\n";
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strstr(input,"ITEM: ATOMS")!=NULL){
				token.read(input,string::WS).next(2);
				int i=0;
				while(!token.end()){
					const std::string tmp=token.next();
					if(tmp=="id") formatAtom.index=i;
					else if(tmp=="mol") formatAtom.mol=i;
					else if(tmp=="type") formatAtom.type=i;
					else if(tmp=="q") formatAtom.q=i;
					else if(tmp=="mass") formatAtom.m=i;
					else if(tmp=="x") formatAtom.x=i;
					else if(tmp=="y") formatAtom.y=i;
					else if(tmp=="z") formatAtom.z=i;
					else if(tmp=="xu") formatAtom.xu=i;
					else if(tmp=="yu") formatAtom.yu=i;
					else if(tmp=="zu") formatAtom.zu=i;
					else if(tmp=="vx") formatAtom.vx=i;
					else if(tmp=="vy") formatAtom.vy=i;
					else if(tmp=="vz") formatAtom.vz=i;
					else if(tmp=="fx") formatAtom.fx=i;
					else if(tmp=="fy") formatAtom.fy=i;
					else if(tmp=="fz") formatAtom.fz=i;
					i++;
				}
				break;
			}
		}
		if(LAMMPS_PRINT_DATA>0) std::cout<<"formatAtom = \n"<<formatAtom<<"\n";
		std::rewind(reader);
		
		//==== read in the number of atoms ====
		if(LAMMPS_PRINT_STATUS>0) std::cout<<"reading number of atoms\n";
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strstr(input,"NUMBER OF ATOMS")!=NULL){
				natoms=std::atoi(fgets(input,string::M,reader));
				break;
			}
		}
		
		//==== find the min index ===
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strstr(input,"ITEM: ATOMS")!=NULL){
				minindex=natoms;
				for(int i=0; i<1; ++i){
					token.read(fgets(input,string::M,reader),string::WS).next(formatAtom.index);
					int index=-1;
					if(formatAtom.index>=0) index=std::atoi(token.next().c_str())-1;
					minindex=index;
				}
				for(int i=1; i<natoms; ++i){
					token.read(fgets(input,string::M,reader),string::WS).next(formatAtom.index);
					//find index
					int index=-1;
					if(formatAtom.index>=0) index=std::atoi(token.next().c_str())-1;
					if(index<minindex) minindex=index;
				}
				break;
			}
		}
		
		std::rewind(reader);
		
		//==== set the timesteps ====
		const int ibeg=Interval::index(interval.beg(),ts);
		const int iend=Interval::index(interval.end(),ts);
		const int tsint=iend-ibeg+1;
		const int nsteps=tsint/interval.stride();
		if(LAMMPS_PRINT_DATA>1){
			std::cout<<"interval  = "<<interval<<"\n";
			std::cout<<"ts        = "<<ts<<"\n";
			std::cout<<"(beg,end) = ("<<ibeg<<","<<iend<<")\n";
			std::cout<<"tsint     = "<<tsint<<"\n";
			std::cout<<"nsteps    = "<<nsteps<<"\n";
		}
		
		//==== resize the simulation ====
		if(LAMMPS_PRINT_STATUS>0) std::cout<<"allocating memory\n";
		traj.resize(nsteps,natoms,atom);
		
		//==== read atoms ====
		if(LAMMPS_PRINT_STATUS>0) std::cout<<"reading atoms\n";
		int timestep=0;
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strstr(input,"BOX")!=NULL){
				//check timestep
				if(timestep<ibeg){continue;}
				if(timestep%interval.stride()!=0){continue;}
				if(LAMMPS_PRINT_DATA>1) std::cout<<"Cell: "<<timestep<<"\n";
				//local variables
				Token token;
				lv.setZero();
				//x
				token.read(fgets(input,string::M,reader),string::WS);
				const double xlob=std::atof(token.next().c_str());//xlo
				const double xhib=std::atof(token.next().c_str());//xhi
				double xy=0;
				if(!token.end()) xy=std::atof(token.next().c_str());//xy
				//y
				token.read(fgets(input,string::M,reader),string::WS);
				const double ylob=std::atof(token.next().c_str());//ylo
				const double yhib=std::atof(token.next().c_str());//yhi
				double xz=0;
				if(!token.end()) xz=std::atof(token.next().c_str());//xz
				//z
				token.read(fgets(input,string::M,reader),string::WS);
				const double zlob=std::atof(token.next().c_str());//zlo
				const double zhib=std::atof(token.next().c_str());//zhi
				double yz=0;
				if(!token.end()) yz=std::atof(token.next().c_str());//yz
				//set xlo,xhi,ylo,yhi,zlo,zhi
				const double xlo=xlob-std::min(0.0,std::min(xy,std::min(xz,xy+xz)));
				const double xhi=xhib-std::max(0.0,std::max(xy,std::max(xz,xy+xz)));
				const double ylo=ylob-std::min(0.0,yz);
				const double yhi=yhib-std::max(0.0,yz);
				const double zlo=zlob;
				const double zhi=zhib;
				//set lv
				lv(0,0)=xhi-xlo;
				lv(1,1)=yhi-ylo;
				lv(2,2)=zhi-zlo;
				lv(0,1)=xy;
				lv(0,2)=xz;
				lv(1,2)=yz;
				//set cell
				static_cast<Cell&>(traj.frame(timestep/interval.stride()-ibeg)).init(lv*s_len);
			}
			if(std::strstr(input,"ITEM: ATOMS")!=NULL){
				if(LAMMPS_PRINT_DATA>1) std::cout<<"Atoms: "<<timestep<<"\n";
				if(timestep<ibeg){++timestep; continue;}
				if(timestep%interval.stride()!=0){++timestep; continue;}
				const int lts=timestep/interval.stride()-ibeg;//local time step
				for(int i=0; i<natoms; ++i){
					std::vector<std::string> tokens;
					token.read(fgets(input,string::M,reader),string::WS);
					while(!token.end()) tokens.push_back(token.next());
					//read in the data
					if(formatAtom.q>=0) dataAtom.q=std::atof(tokens[formatAtom.q].c_str());
					if(formatAtom.m>=0) dataAtom.m=std::atof(tokens[formatAtom.m].c_str());
					if(formatAtom.x>=0) dataAtom.posn[0]=std::atof(tokens[formatAtom.x].c_str());
					if(formatAtom.y>=0) dataAtom.posn[1]=std::atof(tokens[formatAtom.y].c_str());
					if(formatAtom.z>=0) dataAtom.posn[2]=std::atof(tokens[formatAtom.z].c_str());
					if(formatAtom.xu>=0) dataAtom.posn[0]=std::atof(tokens[formatAtom.xu].c_str());
					if(formatAtom.yu>=0) dataAtom.posn[1]=std::atof(tokens[formatAtom.yu].c_str());
					if(formatAtom.zu>=0) dataAtom.posn[2]=std::atof(tokens[formatAtom.zu].c_str());
					if(formatAtom.vx>=0) dataAtom.vel[0]=std::atof(tokens[formatAtom.vx].c_str());
					if(formatAtom.vy>=0) dataAtom.vel[1]=std::atof(tokens[formatAtom.vy].c_str());
					if(formatAtom.vz>=0) dataAtom.vel[2]=std::atof(tokens[formatAtom.vz].c_str());
					if(formatAtom.fx>=0) dataAtom.force[0]=std::atof(tokens[formatAtom.fx].c_str());
					if(formatAtom.fy>=0) dataAtom.force[1]=std::atof(tokens[formatAtom.fy].c_str());
					if(formatAtom.fz>=0) dataAtom.force[2]=std::atof(tokens[formatAtom.fz].c_str());
					if(formatAtom.index>=0) dataAtom.index=std::atoi(tokens[formatAtom.index].c_str())-1;
					if(formatAtom.type>=0) dataAtom.type=std::atoi(tokens[formatAtom.type].c_str())-1;
					//set the simulation data
					const int index=dataAtom.index-minindex;
					if(atom.type()) traj.frame(lts).type(index)=dataAtom.type;
					if(atom.posn()) traj.frame(lts).posn(index)=dataAtom.posn*s_len;
					if(atom.charge()) traj.frame(lts).charge(index)=dataAtom.q;
					if(atom.mass()) traj.frame(lts).mass(index)=dataAtom.m;
					if(atom.vel()) traj.frame(lts).vel(index)=dataAtom.vel*s_len;
					if(atom.force()) traj.frame(lts).force(index)=dataAtom.force*s_energy/s_len;
					if(atom.name()) traj.frame(lts).name(index)=std::string("X")+std::to_string(dataAtom.type);
				}
				++timestep;
			}
			if(timestep>iend) break;
		}
		
		//==== move the atoms within the cell ====
		for(int t=0; t<traj.timesteps(); ++t){
			for(int n=0; n<traj.frame(t).nAtoms(); ++n){
				traj.frame(t).modv(traj.frame(t).posn(n),traj.frame(t).posn(n));
			}
		}
		
		//==== set the mass ====
		if(formatAtom.m<0 && atom.mass()){
			if(atom.an()){
				for(int t=0; t<traj.timesteps(); ++t){
					for(int i=0; i<traj.frame(t).nAtoms(); ++i){
						traj.frame(t).mass(i)=ptable::mass(traj.frame(t).an(i))*s_mass;
					}
				}
			} else if(atom.name()){
				for(int t=0; t<traj.timesteps(); ++t){
					for(int i=0; i<traj.frame(t).nAtoms(); ++i){
						const int an=ptable::an(traj.frame(t).name(i).c_str());
						traj.frame(t).mass(i)=ptable::mass(an)*s_mass;
					}
				}
			}
		}

		//==== stop the timer ====
		stop=std::clock();
		
		//==== print the time ====
		time=((double)(stop-start))/CLOCKS_PER_SEC;
		std::cout<<"positions read in "<<time<<" seconds\n";
	}catch(std::exception& e){
		std::cout<<"ERROR in "<<NAMESPACE_GLOBAL<<"::"<<funcName<<":\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	//free local variables
	if(reader!=NULL) fclose(reader);
	delete[] input;
	delete[] temp;
	
	if(error) throw std::runtime_error("I/O Exception Occurred.");
}

void write(const char* file, const Interval& interval, const Atom& atom, const Trajectory& traj){
	static const char* funcName="write<AtomT>(const char*,const Interval&,const Atom&,Trajectory&)";
	if(LAMMPS_PRINT_FUNC>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<funcName<<":\n";
	//local variables
	FILE* writer=NULL;
	bool error=false;
	
	try{
		//open the file
		writer=fopen(file,"w");
		if(writer==NULL) throw std::runtime_error("Unable to open file.");
		
		//set the beginning and ending timesteps
		const int ibeg=Interval::index(interval.beg(),traj.timesteps());
		const int iend=Interval::index(interval.end(),traj.timesteps());
		
		for(int t=ibeg; t<=iend; ++t){
			fprintf(writer,"ITEM: TIMESTEP\n");
			fprintf(writer,"%i\n",t);
			fprintf(writer,"ITEM: NUMBER OF ATOMS\n");
			fprintf(writer,"%i\n",traj.frame(t).nAtoms());
			fprintf(writer,"ITEM: BOX BOUNDS pp pp pp\n");
			const Eigen::Vector3d ra=traj.frame(t).R().col(0);
			const Eigen::Vector3d rb=traj.frame(t).R().col(1);
			const Eigen::Vector3d rc=traj.frame(t).R().col(2);
			const double A=ra.norm();
			const double B=rb.norm();
			const double C=rc.norm();
			const double cosA=rb.dot(rc)/(B*C);
			const double cosB=ra.dot(rc)/(A*C);
			const double cosC=ra.dot(rb)/(A*B);
			const double lx=A;
			const double xy=B*cosC;
			const double xz=C*cosB;
			const double ly=sqrt(B*B-xy*xy);
			const double yz=(B*C*cosA-xy*xz)/ly;
			const double lz=sqrt(C*C-xz*xz-yz*yz);
			fprintf(writer,"%f %f %f\n",0.0,lx,xy);
			fprintf(writer,"%f %f %f\n",0.0,ly,xz);
			fprintf(writer,"%f %f %f\n",0.0,lz,yz);
			fprintf(writer,"ITEM: ATOMS id type x y z\n");
			for(int n=0; n<traj.frame(t).nAtoms(); ++n){
				fprintf(writer,"%i %i %f %f %f\n",n+1,traj.frame(t).type(n),
					traj.frame(t).posn(n)[0],traj.frame(t).posn(n)[1],traj.frame(t).posn(n)[2]
				);
			}
		}
		
		fclose(writer);
		writer=NULL;
	}catch(std::exception& e){
		std::cout<<"ERROR in "<<NAMESPACE_GLOBAL<<"::"<<funcName<<":\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	if(error) throw std::runtime_error("I/O Exception Occurred.");
}
	
}

}
