/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Node.h"
#include "ElementGroup.h"
#include "Outputter.h"
#include "Solver.h"
#include "LoadCaseData.h"
#include "SkylineMatrix.h"
#include "EssentialBoudary.h"

using namespace std;

//!	Clear an array
template <class type> void clear( type* a, unsigned int N );

//!	Domain class : Define the problem domain
/*!	Only a single instance of Domain class can be created */
class CDomain
{
private:

//!	The instance of the Domain class
	static CDomain* _instance;

//!	Input file stream for reading data from input data file
	ifstream Input;

//!	Heading information for use in labeling the outpu
	char Title[256]; 

//!	Solution MODEX
/*!		0 : Data check only;
		1 : Execution */
	unsigned int MODEX;

//!	Total number of nodal points
	unsigned int NUMNP;

//!	List of all nodes in the domain
	CNode* NodeList;

//!	Total number of element groups.
/*! An element group consists of a convenient collection of elements with same type */
	unsigned int NUMEG;

//! Element group list
    CElementGroup* EleGrpList;
    
//!	Number of load cases
	unsigned int NLCASE;

//!	List of all load cases
	CLoadCaseData* LoadCases;

//! Essential boundary data
	CEBData* EBdata;

//!	Number of concentrated loads applied in each load case
	unsigned int* NLOAD;

//!	Total number of equations in the system
	unsigned int NEQ;

//! Total number of boudary conditions (bcode = 2)
	unsigned int NBC;

//!	Banded stiffness matrix
/*! A one-dimensional array storing only the elements below the	skyline of the 
    global stiffness matrix. */
    CSkylineMatrix<double>* StiffnessMatrix;

//! Banded global stiffness matrix including dof on essential boudary conditions
/*! A one-dimensional array storing only the elements below the	skyline of the 
    global stiffness matrix. */
	CSkylineMatrix<double>* GlobalStiffnessMatrix;

//!	Global nodal force/displacement vector
	double* Force;

//! Global nodal force (including boundary condition)
	double* GlobalNodalForce;

//! Displacement at essential boundary
	double* EssentialDisplacement;

//! Global displacement
	double* GlobalDisplacement;

private:

//!	Constructor
	CDomain();

//!	Desconstructor
	~CDomain();

public:

//!	Return pointer to the instance of the Domain class
	static CDomain* GetInstance();

//!	Read domain data from the input data file
	bool ReadData(string FileName, string OutFile);

//!	Read nodal point data
	bool ReadNodalPoints();

//!	Read load case data
	bool ReadLoadCases();

//! Read essential boundary data
	bool ReadEBData();

//!	Read element data
	bool ReadElements();

//!	Calculate global equation numbers corresponding to every degree of freedom of each node
	void CalculateEquationNumber();

//!	Calculate column heights
	void CalculateColumnHeights();

//! Calculate global column heights
	void CalculateGlobalColumnHeights();

//! Allocate storage for matrices
/*!	Allocate Force, ColumnHeights, DiagonalAddress and StiffnessMatrix and 
    calculate the column heights and address of diagonal elements */
	void AllocateMatrices();

//! Allocate storage for global matrices
/*!	Allocate GlobalForce, ColumnHeights, DiagonalAddress and GlobalStiffnessMatrix and 
    calculate the column heights and address of diagonal elements */
	void AllocateGlobalMatrices();

//!	Assemble the banded gloabl stiffness matrix
	void AssembleStiffnessMatrix();

//! Assemble global stiffness matrix
	void AssembleGlobalStiffnessMatrix();

//!	Assemble the global nodal force vector for load case LoadCase
	bool AssembleForce(unsigned int LoadCase); 

//! Assemble the essential boundary displacement vector
	bool AssmbleEssentialDisplacement();

//! Substitute contribution of essential boundary of force vector: \hat{f_x} = f_x - K_{ux}d_u
	bool SubstitueEssentialBoundary();

//! Compute nodal force (including nodes on essential boudndary)
	bool ComputeNodalForce();

//! Concate to get global displacement
	void ConcateGlobalDisplacement();

//!	Return solution mode
	inline unsigned int GetMODEX() { return MODEX; }

//!	Return the title of problem
	inline string GetTitle() { return Title; }

//!	Return the total number of equations
	inline unsigned int GetNEQ() { return NEQ; }

//! Return the number of essential boundary conditon
	inline unsigned int GetNBC() { return NBC; }

//!	Return the total number of nodal points
	inline unsigned int GetNUMNP() { return NUMNP; }

//!	Return the node list
	inline CNode* GetNodeList() { return NodeList; }

//!	Return total number of element groups
	inline unsigned int GetNUMEG() { return NUMEG; }

//! Return element group list
    inline CElementGroup* GetEleGrpList() { return EleGrpList; }

//!	Return pointer to the global nodal force vector
	inline double* GetForce() { return Force; }

//!	Return pointer to the global nodal displacement vector
	inline double* GetDisplacement() { return Force; }

//! Return pointer to the essential displacement vector
	inline double* GetEssentialDisplacement() { return EssentialDisplacement; }

//! Return pointer to the global displacement vector
	inline double* GetGlobalDisplacement() { return GlobalDisplacement; }

//!	Return the total number of load cases
	inline unsigned int GetNLCASE() { return NLCASE; }

//!	Return the number of concentrated loads applied in each load case
	inline unsigned int* GetNLOAD() { return NLOAD; }

//!	Return the list of load cases
	inline CLoadCaseData* GetLoadCases() { return LoadCases; }

//! Return the EBdata
	inline CEBData* GetEBdata() { return EBdata; }

//! Return pointer to global nodal force
	inline double* GetNodalForce() { return GlobalNodalForce; }

//!	Return pointer to the banded stiffness matrix
	inline CSkylineMatrix<double>* GetStiffnessMatrix() { return StiffnessMatrix; }

};
