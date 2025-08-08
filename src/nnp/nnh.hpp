#pragma once
#ifndef NNH_HPP
#define NNH_HPP

// c++ libraries
#include <iosfwd>
// structure
#include "struc/structure.hpp"
#include "struc/neighbor.hpp"
#include "struc/type.hpp"
// ml
#include "ml/nn.hpp"
// string
#include "str/string.hpp"
// nnp
#include "nnp/tensortri.hpp"
#include "nnp/basis_radial.hpp"
#include "nnp/basis_angular.hpp"

//***********************************************************************
// COMPILER DIRECTIVES
//***********************************************************************

#ifndef NNH_PRINT_FUNC
#define NNH_PRINT_FUNC 0
#endif

#ifndef NNH_PRINT_STATUS
#define NNH_PRINT_STATUS 0
#endif

#ifndef NNH_PRINT_DATA
#define NNH_PRINT_DATA 0
#endif

//************************************************************
// NEURAL NETWORK HAMILTONIAN (NNH)
//************************************************************

/**
* Class defining a Neural Network Hamiltonian
* This class contains all the data and methods necessary to define a Hamiltonian yielding the
* potential energy surface for a specific atom type given a set of inputs to the neural network, 
* the output of which yields the atomic energy for each set of inputs.
*/
class NNH{
public:
	class Type{
	private:
		std::string name_;
		double energy_;
		double mass_;
		double radius_;
		double amp_;
	public:
		//===== constructors/destructors ====
		Type(){clear();}
		Type(std::string& name):name_(name){};
		~Type(){}

		//==== operators ====
		friend std::ostream& operator<<(std::ostream& out, const Type& type);

		//==== access ====
		std::string& name(){return name_;}
		const std::string& name()const{return name_;}
		double& energy(){return energy_;}
		const double& energy()const{return energy_;}
		double& mass(){return mass_;}
		const double& mass()const{return mass_;}
		double& radius(){return radius_;}
		const double& radius()const{return radius_;}
		double& amp(){return amp_;}
		const double& amp()const{return amp_;}

		//==== member functions ====
		void clear();
		Type& read(const char* str);
		Type& read(Token& token);
		void write(FILE* out)const;
	};
private:
	//network configuration
	int nInput_;//number of radial + angular symmetry functions
	int nInputR_;//number of radial symmetry functions
	int nInputA_;//number of angular symmetry functions
	
	//hamiltonian
	int ntypes_;//the total number of types
	Type type_;//type associated with NNH
	NN::ANN nn_;//neural network hamiltonian
	NN::DODZ dOdZ_;//gradient of the output w.r.t. node values
	
	//basis for pair/triple interactions
	TensorTri<1,BasisR> basisR_;//radial basis functions (ntypes_)
	TensorTri<1,int> offsetR_;//offset for the given radial basis (ntypes_)
	TensorTri<2,BasisA> basisA_;//angular basis functions (ntypes x (ntypes+1)/2)
	TensorTri<2,int> offsetA_;//offset for the given radial basis (ntypes x (ntypes+1)/2)
public:
	//==== constructors/destructors ====
	NNH(){defaults();}
	~NNH(){}
	
	//operators
	friend std::ostream& operator<<(std::ostream& out, const NNH& nmh);
	
	//==== access ====
	//hamiltonian
		const int& ntypes()const{return ntypes_;}
		Type& type(){return type_;}
		const Type& type()const{return type_;}
		NN::ANN& nn(){return nn_;}
		const NN::ANN& nn()const{return nn_;}
		NN::DODZ& dOdZ(){return dOdZ_;}
		const NN::DODZ& dOdZ()const{return dOdZ_;}
	//basis for pair/triple interactions
		BasisR& basisR(int i){return basisR_[i];}
		const BasisR& basisR(int i)const{return basisR_[i];}
		BasisA& basisA(int i, int j){return basisA_(i,j);}
		const BasisA& basisA(int i, int j)const{return basisA_(i,j);}
	//network configuration
		const int& nInput()const{return nInput_;}
		const int& nInputR()const{return nInputR_;}
		const int& nInputA()const{return nInputA_;}
		const int& offsetR(int i)const{return offsetR_[i];}
		const int& offsetA(int i, int j)const{return offsetA_(i,j);}
	
	//==== member functions ====
	//misc
		void defaults();//set defaults
		void clear(){defaults();}//clear the potential
	//resizing
		void resize(int ntypes);//resize
		void init_input();//initialize the inputs
	//output
		double energy(const Eigen::VectorXd& symm);//compute energy of atom
};

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> int nbytes(const NNH::Type& obj);
	template <> int nbytes(const NNH& obj);
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> int pack(const NNH::Type& obj, char* arr);
	template <> int pack(const NNH& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> int unpack(NNH::Type& obj, const char* arr);
	template <> int unpack(NNH& obj, const char* arr);
	
}

#endif