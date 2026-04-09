#pragma once
#ifndef NEIGHBOR_HPP
#define NEIGHBOR_HPP

// c++
#include <iostream>
// struc
#include "struc/structure.hpp"

//***********************************************************
// NeighborList
//***********************************************************

class NeighborList{
private:
    bool self_;//self-interaction
    double rc_,rc2_;//cutoff radius
    std::vector<std::vector<int> > indices_;
    std::vector<std::vector<Eigen::Vector3d> > img_;
public:
    //==== constructors/destructors ====
    NeighborList():self_(true),rc_(0),rc2_(0){}
    NeighborList(double rc);
    NeighborList(double rc, const Structure& struc);
    ~NeighborList(){}

    //==== operator ====
    friend std::ostream& operator<<(std::ostream& out, const NeighborList& nlist);

    //==== access ====
    bool& self(){return self_;}
    const bool& self()const{return self_;}
    const double& rc()const{return rc_;}
    int size()const{return indices_.size();}
    int size(int i)const{return indices_[i].size();}
    const std::vector<int>& indices(int i)const{return indices_[i];}
    const std::vector<Eigen::Vector3d>& img(int i)const{return img_[i];}
    const int& index(int i, int j)const{return indices_[i][j];}
    const Eigen::Vector3d& img(int i, int j)const{return img_[i][j];}

    //==== member functions ====
    void clear();
    void build(const Structure& struc);
    void build_clist(const Structure& struc);
};

#endif