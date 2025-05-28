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
	output << setw(16) << E << setw(16) << t << endl;
}

// Compute elasticity matrix
void CT3Material::ComputeElasticityMatrix(){
	// For simplicity, store all data of the matrix
	if(!PlaneStress){
		// If the material is plane strain problem
		// convert to elasticity property into plane stress problem
		E = E/(1-nu*nu);
		nu = nu/(1-nu);
	}

	// compute coefficient
	double coefficient = E/(1-nu*nu);

	// compute matrix
	D[0][0] = coefficient * 1;       // D11
	D[0][1] = coefficient * nu;      // D12
	D[0][2] = 0;                     // D13
	D[1][0] = coefficient * nu;      // D21
	D[1][1] = coefficient * 1;       // D22
	D[1][2] = 0;                     // D23
	D[2][0] = 0;                     // D31
	D[2][1] = 0;                     // D32
	D[2][2] = coefficient * (1 - nu) / 2; // D33 (剪切项)
}