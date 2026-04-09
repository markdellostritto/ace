// c++
#include <stdexcept>
// math
#include "math/const.hpp"
// struc
#include "struc/neighbor.hpp"
#include "struc/clist.hpp"

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
                const double ds2=(dr-struc.R()*ivecs[i]).squaredNorm();//add to m : subtract from dr
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

void NeighborList::build_clist(const Structure& struc){
    //compute image vectors
    Eigen::Vector3i cellf=Eigen::Vector3i::Zero();
    if(struc.R().norm()>ZERO){
        for(int n=0; n<3; ++n){
            cellf[n]=std::ceil(3.0*rc_/struc.R().col(n).norm());
        }
    }
    const int Rsize=(2*cellf[0]+1)*(2*cellf[1]+1)*(2*cellf[2]+1);
    //resize neighbor lists
    indices_.resize(struc.nAtoms());
    img_.resize(struc.nAtoms());
    if(cellf.maxCoeff()>1){
        //make neighbor cell vectors
        std::vector<Eigen::Vector3d> ivecs(Rsize);
        int c=0;
        for(int i=-cellf[0]; i<=cellf[0]; ++i){
            for(int j=-cellf[1]; j<=cellf[1]; ++j){
                for(int k=-cellf[2]; k<=cellf[2]; ++k){
                    ivecs[c++]<<1.0*i,1.0*j,1.0*k;
                }
            }
        }
        //compute
        for(int n=0; n<struc.nAtoms(); ++n){
            indices_[n].clear();
            img_[n].clear();
            for(int m=0; m<struc.nAtoms(); ++m){
                const Eigen::Vector3d dr=struc.posn(n)-struc.posn(m);
                //loop over all nearby cells
                for(int i=0; i<ivecs.size(); ++i){
                    //compute distance
                    const double ds2=(dr-struc.R()*ivecs[i]).squaredNorm();//add to m : subtract from dr
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
    } else {
        Eigen::Vector3d drv;
        //compute nnc
        std::vector<Eigen::Vector3i> nnc_(27);
        int count=0;
        for(int i=-1; i<=1; ++i){
            for(int j=-1; j<=1; ++j){
                for(int k=-1; k<=1; ++k){
                    nnc_[count++]<<i,j,k;
                }
            }
        }
        //make the cell list
        CellList cellList(rc_,struc);
        //compute
        for(int n=0; n<struc.nAtoms(); ++n){
            //clear neighbor list
            indices_[n].clear();
            img_[n].clear();
            //get cell for atom n
            const Eigen::Vector3i& cell=cellList.cell(n);
            //loop over all nearest-neighbor-cells
            for(int i=0; i<nnc_.size(); ++i){
                const Eigen::Vector3i ncell=cell+nnc_[i];
                for(int j=0; j<cellList.atoms(ncell).size(); ++j){
                    const int m=cellList.atoms(ncell)[j];
                    const double dr2=struc.dist2(struc.posn(n),struc.posn(m),drv);
                    if(math::constants::ZERO<dr2 && dr2<rc2_){
                        indices_[n].push_back(m);
                        img_[n].push_back(Eigen::Vector3d::Zero());
                    }
                }
            }
        }
    }
}
