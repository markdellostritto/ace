// c++ libraries
#include <iostream>
// string
#include "str/print.hpp"
// math
#include "math/special.hpp"
// structure
#include "struc/structure.hpp"
#include "struc/clist.hpp"

//==== operators ====

std::ostream& operator<<(std::ostream& out, const CellList& c){
	char* str=new char[print::len_buf];
	out<<print::buf(str)<<"\n";
	out<<print::title("CELL LIST",str)<<"\n";
	out<<"DIM  = "<<c.dim_[0]<<" "<<c.dim_[1]<<" "<<c.dim_[2]<<"\n";
	out<<"FLEN = "<<c.flen_[0]<<" "<<c.flen_[1]<<" "<<c.flen_[2]<<"\n";
	out<<print::title("CELL LIST",str)<<"\n";
	out<<print::buf(str);
	delete[] str;
	return out;
}

//==== member functions ====

const std::vector<int>& CellList::atoms(int i, int j, int k)const{
	return atoms_[
		index(
			math::special::mod(i,dim_[0]),
			math::special::mod(j,dim_[1]),
			math::special::mod(k,dim_[2])
		)
	];
}

const std::vector<int>& CellList::atoms(const Eigen::Vector3i& i)const{
	return atoms_[
		index(
			math::special::mod(i[0],dim_[0]),
			math::special::mod(i[1],dim_[1]),
			math::special::mod(i[2],dim_[2])
		)
	];
}
	
void CellList::clear(){
	dim_[0]=0; dim_[1]=0; dim_[2]=0;
	flen_[0]=0.0; flen_[1]=0.0; flen_[2]=0.0;
	cell_.clear();
	atoms_.clear();
}

int CellList::index(int i, int j, int k)const{
	return i*dim_[1]*dim_[2]+j*dim_[2]+k;
}

void CellList::compute(double rc, const Structure& struc){
	//compute the dimension
	dim_[0]=std::ceil(struc.R().row(0).lpNorm<Eigen::Infinity>()/rc);
	dim_[1]=std::ceil(struc.R().row(1).lpNorm<Eigen::Infinity>()/rc);
	dim_[2]=std::ceil(struc.R().row(2).lpNorm<Eigen::Infinity>()/rc);
	flen_[0]=1.0/(1.0*dim_[0]);
	flen_[1]=1.0/(1.0*dim_[1]);
	flen_[2]=1.0/(1.0*dim_[2]);
	//resize
	atoms_.clear();
	atoms_.resize(dim_[0]*dim_[1]*dim_[2]);
	cell_.resize(struc.nAtoms());
	//iterate over all atoms
	for(int n=0; n<struc.nAtoms(); ++n){
		//convert to frac coordinates
		const Eigen::Vector3d p=struc.RInv()*struc.posn(n);
		//bin the position
		cell_[n][0]=math::special::mod((int)(p[0]/flen_[0]),dim_[0]);
		cell_[n][1]=math::special::mod((int)(p[1]/flen_[1]),dim_[1]);
		cell_[n][2]=math::special::mod((int)(p[2]/flen_[2]),dim_[2]);
		//add to the cell list
		atoms_[index(cell_[n][0],cell_[n][1],cell_[n][2])].push_back(n);
	}
}