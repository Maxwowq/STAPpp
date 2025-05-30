/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Material.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//	Read material data from stream Input
bool CBarMaterial::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> Area;	// Young's modulus and section area

	return true;
}

//	Write material data to Stream
void CBarMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Area << endl;
}

//	Read material data from stream Input
bool CT3Material::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> nu >> t >> PlaneStress;	// Young's modulus, Poisson's ratio, thickness and plane stress indicator

	this->ComputeElasticityMatrix(); // Compute and store Elasticity matrix after read info
	return true;
}

//	Write material data to Stream
void CT3Material::Write(COutputter& output)
{
	output << setw(15) << E << setw(15) << nu << setw(15) << t << setw(15);
	if(PlaneStress){
		output << "PlaneStress" << endl;
	}
	else
	{
		output << "PlaneStrain" << endl;
	}
}

// Compute elasticity matrix
void CT3Material::ComputeElasticityMatrix(){
	// For simplicity, store all data of the matrix
	double E_, nu_;
	if(!PlaneStress){
		// If the material is plane strain problem
		// convert to elasticity property into plane stress problem
		E_ = E/(1-nu*nu);
		nu_ = nu/(1-nu);
	}
	else{
		E_ = E;
		nu_ = nu;
	}

	// compute coefficient
	double coefficient = E_/(1-nu_*nu_);

	// compute matrix
	D[0][0] = coefficient * 1;       // D11
	D[0][1] = coefficient * nu_;      // D12
	D[0][2] = 0;                     // D13
	D[1][0] = coefficient * nu_;      // D21
	D[1][1] = coefficient * 1;       // D22
	D[1][2] = 0;                     // D23
	D[2][0] = 0;                     // D31
	D[2][1] = 0;                     // D32
	D[2][2] = coefficient * (1 - nu_) / 2; // D33 (剪切项)
}