// c++
#include <stdexcept>
// math
#include "math/const.hpp"
// struc
#include "struc/neighbor.hpp"

using math::constants::ZERO;

//***********************************************************
// NeighborList
//***********************************************************

//==== constructors/destructors ====

NeighborList::NeighborList(double rc){
    if(rc<0) throw std::invalid_argument("NeighborList::NeighborList(double): invalid cutoff radius.");
    rc_=rc;
    rc2_=rc_*rc_;
    self_=true;
}

NeighborList::NeighborList(double rc, const Structure& struc){
    if(rc<0) throw std::invalid_argument("NeighborList::NeighborList(double,const Structure&): invalid cutoff radius.");
    rc_=rc;
    rc2_=rc_*rc_;
    self_=true;
    build(struc);
}

//==== operators ====

std::ostream& operator<<(std::ostream& out, const NeighborList& nlist){
    return out<<"neighbor rc "<<nlist.rc()<<" self "<<nlist.self();
}

//==== member functions ====

void NeighborList::clear(){
    rc_=0;
    indices_.clear();
    img_.clear();
}

void NeighborList::build(const Structure& struc){
    //compute image vectors
    Eigen::Vector3i cellf=Eigen::Vector3i::Zero();
    if(struc.R().norm()>ZERO){
        for(int n=0; n<3; ++n){
            cellf[n]=std::ceil(3.0*rc_/struc.R().col(n).norm());
        }
    }
    const int Rsize=(2*cellf[0]+1)*(2*cellf[1]+1)*(2*cellf[2]+1);
    std::vector<Eigen::Vector3d> ivecs(Rsize);
    int c=0;
    for(int i=-cellf[0]; i<=cellf[0]; ++i){
        for(int j=-cellf[1]; j<=cellf[1]; ++j){
            for(int k=-cellf[2]; k<=cellf[2]; ++k){
                ivecs[c++]<<1.0*i,1.0*j,1.0*k;
            }
        }
    }
    //resize neighbor lists
    indices_.resize(struc.nAtoms());
    img_.resize(struc.nAtoms());
    //compute
    for(int n=0; n<struc.nAtoms(); ++n){
        indices_[n].clear();
        img_[n].clear();
        for(int m=0; m<struc.nAtoms(); ++m){
            const Eigen::Vector3d dr=struc.posn(n)-struc.posn(m);
            //loop over all nearby cells
            for(int i=0; i<ivecs.size(); ++i){
                //compute distance
                const Eigen::Vector3d ds=dr-struc.R()*ivecs[i];//add to m : subtract from dr
                const double ds2=ds.squaredNorm();
                //check for self-image
                if(m!=n || (self_ && ds2>ZERO)){
                    if(ds2<rc2_){
                        indices_[n].push_back(m);
                        img_[n].push_back(ivecs[i]);
                    }
                }
            }
        }
    }
}

//***********************************************************
// PairList
//***********************************************************

void PairList::compute(const Structure& struc, const NeighborList& nlist){
    const int natoms=nlist.size();
    indices_.resize(natoms);
    dr2_.resize(natoms);
    disp_.resize(natoms);
    for(int i=0; i<natoms; ++i){
        indices_[i].resize(nlist.size(i));
        dr2_[i].resize(nlist.size(i));
        disp_[i].resize(nlist.size(i));
        for(int j=0; j<nlist.size(i); ++j){
            const int jj=nlist.index(i,j);
            indices_[i][j]=jj;
            disp_[i][j]=struc.posn(i)-(struc.posn(jj)+struc.R()*nlist.img(i,j));
            dr2_[i][j]=disp_[i][j].squaredNorm();
        }
    }
}