// c libraries
#include <cstdlib>
#include <cstdio>
// c++ libraries
#include <iostream>
// io
#include "str/parse.hpp"
#include "str/string.hpp"
// structure
#include "struc/structure.hpp"
// format
#include "format/file_struc.hpp"
#include "format/format.hpp"
// math
#include "math/const.hpp"

int main(int argc, char* argv[]){
	//simulation
		Structure struc;
		FILE_FORMAT::type format_in;
		FILE_FORMAT::type format_out;
	//file i/o
		std::string file_in;
		std::string file_out;
	//arguments
		std::vector<input::Arg> args;
	//atom type
		Atom atom;
		atom.name()=true;
		atom.an()=true;
		atom.type()=true;
		atom.posn()=true;
		atom.force()=true;
		atom.frac()=true;
	//misc
		bool error=false;
		
	try{
		//check the number of arguments
		if(argc==1) throw std::invalid_argument("No arguments provided.");
		
		//parse the arguments
		input::parse(argc,argv,args);
		for(int i=0; i<args.size(); ++i){
			if(args[i].key()=="fin") format_in=FILE_FORMAT::read(args[i].val(0).c_str());
			else if(args[i].key()=="fout") format_out=FILE_FORMAT::read(args[i].val(0).c_str());
			else if(args[i].key()=="in") file_in=args[i].val(0).c_str();
			else if(args[i].key()=="out") file_out=args[i].val(0).c_str();
			else if(args[i].key()=="cart") atom.frac()=false;
			else if(args[i].key()=="frac") atom.frac()=true;
		}
		
		//print parameters
		std::cout<<"format-in  = "<<format_in<<"\n";
		std::cout<<"format-out = "<<format_out<<"\n";
		std::cout<<"file-in    = "<<file_in<<"\n";
		std::cout<<"file-out   = "<<file_out<<"\n";
		std::cout<<"frac       = "<<atom.frac()<<"\n";
		
		//check parameters
		if(format_in==FILE_FORMAT::NONE) throw std::invalid_argument("Invalid input format.");
		if(format_out==FILE_FORMAT::NONE) throw std::invalid_argument("Invalid output format.");
		
		//read
		std::cout<<"reading structure\n";
		read_struc(file_in.c_str(),format_in,atom,struc);
		
		//write
		std::cout<<"writing structure\n";
		write_struc(file_out.c_str(),format_out,atom,struc);
		
	}catch(std::exception& e){
		std::cout<<"ERROR in convert_struc::main(int,char**):\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	if(error) return 1;
	else return 0;
}