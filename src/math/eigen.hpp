#pragma once
#ifndef EIGEN_HPP
#define EIGEN_HPP

#define EIGEN_NO_DEBUG

//c++ libraries
#include <iosfwd>
//eigen
#include <Eigen/Dense>
//ann serialize
#include "mem/serialize.hpp"

namespace serialize{

//**********************************************
// byte measures
//**********************************************

template <> int nbytes(const Eigen::Vector3i& obj);
template <> int nbytes(const Eigen::Vector3d& obj);
template <> int nbytes(const Eigen::VectorXd& obj);
template <> int nbytes(const Eigen::Matrix3d& obj);
template <> int nbytes(const Eigen::MatrixXd& obj);
template <> int nbytes(const std::vector<Eigen::Vector3d>& obj);
template <> int nbytes(const std::vector<Eigen::VectorXd>& obj);
template <> int nbytes(const std::vector<Eigen::MatrixXd>& obj);

//**********************************************
// packing
//**********************************************

template <> int pack(const Eigen::Vector3i& obj, char* arr);
template <> int pack(const Eigen::Vector3d& obj, char* arr);
template <> int pack(const Eigen::VectorXd& obj, char* arr);
template <> int pack(const Eigen::Matrix3d& obj, char* arr);
template <> int pack(const Eigen::MatrixXd& obj, char* arr);
template <> int pack(const std::vector<Eigen::Vector3d>& obj, char* arr);
template <> int pack(const std::vector<Eigen::VectorXd>& obj, char* arr);
template <> int pack(const std::vector<Eigen::MatrixXd>& obj, char* arr);

//**********************************************
// unpacking
//**********************************************

template <> int unpack(Eigen::Vector3i& obj, const char* arr);
template <> int unpack(Eigen::Vector3d& obj, const char* arr);
template <> int unpack(Eigen::VectorXd& obj, const char* arr);
template <> int unpack(Eigen::Matrix3d& obj, const char* arr);
template <> int unpack(Eigen::MatrixXd& obj, const char* arr);
template <> int unpack(std::vector<Eigen::Vector3d>& obj, const char* arr);
template <> int unpack(std::vector<Eigen::VectorXd>& obj, const char* arr);
template <> int unpack(std::vector<Eigen::MatrixXd>& obj, const char* arr);

}

#endif