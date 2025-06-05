/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Domain.h"
#include "Material.h"

using namespace std;

//	Clear an array
template <class type> void clear( type* a, unsigned int N )
{
	for (unsigned int i = 0; i < N; i++)
		a[i] = 0;
}

CDomain* CDomain::_instance = nullptr;

//	Constructor
CDomain::CDomain()
{
	Title[0] = '0';
	MODEX = 0;

	NUMNP = 0;
	NodeList = nullptr;
	
	NUMEG = 0;
	EleGrpList = nullptr;
	
	NLCASE = 0;
	NLOAD = nullptr;
	LoadCases = nullptr;
	
	NEQ = 0;

	Force = nullptr;
	StiffnessMatrix = nullptr;
}

//	Desconstructor
CDomain::~CDomain()
{
	delete [] NodeList;

	delete [] EleGrpList;

	delete [] NLOAD;
	delete [] LoadCases;

	delete [] Force;
	delete StiffnessMatrix;
}

//	Return pointer to the instance of the Domain class
CDomain* CDomain::GetInstance()
{
	if (!_instance) 
		_instance = new CDomain();
	
	return _instance;
}

//	Read domain data from the input data file
bool CDomain::ReadData(string FileName, string OutFile)
{
	Input.open(FileName);

	if (!Input) 
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}

	COutputter* Output = COutputter::GetInstance(OutFile);

//	Read the heading line
	Input.getline(Title, 256);
	Output->OutputHeading();

//	Read the control line
	Input >> NUMNP >> NUMEG >> NLCASE >> MODEX;

//	Read nodal point data
	if (ReadNodalPoints())
        Output->OutputNodeInfo();
    else
        return false;

//	Update equation number
	CalculateEquationNumber();
	Output->OutputEquationNumber();

// Read essential boundray data
	if (ReadEBData())
		Output->OutputEBinfo();
	else
		return false;

//	Read load data
	if (ReadLoadCases())
        Output->OutputLoadInfo();
    else
        return false;

//	Read element data
	if (ReadElements())
        Output->OutputElementInfo();
    else
        return false;

	return true;
}

//	Read nodal point data
bool CDomain::ReadNodalPoints()
{

//	Read nodal point data lines
	NodeList = new CNode[NUMNP];

//	Loop over for all nodal points
	for (unsigned int np = 0; np < NUMNP; np++)
    {
		if (!NodeList[np].Read(Input))
			return false;
    
        if (NodeList[np].NodeNumber != np + 1)
        {
            cerr << "*** Error *** Nodes must be inputted in order !" << endl
            << "   Expected node number : " << np + 1 << endl
            << "   Provided node number : " << NodeList[np].NodeNumber << endl;
        
            return false;
        }
    }

	return true;
}

//	Calculate global equation numbers corresponding to every degree of freedom of each node
void CDomain::CalculateEquationNumber()
{
	NEQ = 0;
	NBC = 0;
	for (unsigned int np = 0; np < NUMNP; np++)	// Loop over for all node
	{
		for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
		{
			if (NodeList[np].bcode[dof] == 1)
			{ 
				NodeList[np].bcode[dof] = 0;
				NodeList[np].gbcode[dof] = 0;
			}
			else if(NodeList[np].bcode[dof] == 0)
			{
				NEQ++;
				NodeList[np].bcode[dof] = NEQ;
				NodeList[np].gbcode[dof] = NEQ;
			}
		}
	}
	for (unsigned int np = 0; np < NUMNP; np++)	// Loop over for all node
	{
		for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
		{
			// for essential boundary
			if (NodeList[np].bcode[dof] == 2)
			{ 
				// set bcode to 0, so the bcode is consistent with former definition
				NodeList[np].bcode[dof] = 0;
				// set gbcode to NEQ + NBC
				NBC++;
				NodeList[np].gbcode[dof] = NEQ + NBC;
			}
		}
	}
}

//	Read load case data
bool CDomain::ReadLoadCases()
{
//	Read load data lines
	LoadCases = new CLoadCaseData[NLCASE];	// List all load cases

//	Loop over for all load cases
	for (unsigned int lcase = 0; lcase < NLCASE; lcase++)
    {
        unsigned int LL;
        Input >> LL;
        
        if (LL != lcase + 1)
        {
            cerr << "*** Error *** Load case must be inputted in order !" << endl
            << "   Expected load case : " << lcase + 1 << endl
            << "   Provided load case : " << LL << endl;
            
            return false;
        }

        LoadCases[lcase].Read(Input);
    }

	return true;
}

// Read essential boundary data
bool CDomain::ReadEBData()
{
	EBdata = new CEBData;
	EBdata->Read(Input);
}

// Read element data
bool CDomain::ReadElements()
{
    EleGrpList = new CElementGroup[NUMEG];

//	Loop over for all element group
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
        if (!EleGrpList[EleGrp].Read(Input))
            return false;
    
    return true;
}

//	Calculate column heights
void CDomain::CalculateColumnHeights()
{
#ifdef _DEBUG_
    COutputter* Output = COutputter::GetInstance();
    *Output << setw(9) << "Ele = " << setw(22) << "Location Matrix" << endl;
#endif

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)		//	Loop over for all element groups
    {
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();
        
		for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp
        {
            CElement& Element = ElementGrp[Ele];

            // Generate location matrix
            Element.GenerateLocationMatrix();
            
#ifdef _DEBUG_
            unsigned int* LocationMatrix = Element.GetLocationMatrix();
            
            *Output << setw(9) << Ele+1;
            for (int i=0; i<Element.GetND(); i++)
                *Output << setw(5) << LocationMatrix[i];
            *Output << endl;
#endif

            StiffnessMatrix->CalculateColumnHeight(Element.GetLocationMatrix(), Element.GetND());
        }
    }
    
    StiffnessMatrix->CalculateMaximumHalfBandwidth();
    
#ifdef _DEBUG_
    *Output << endl;
	Output->PrintColumnHeights();
#endif

}

//	Calculate global column heights
void CDomain::CalculateGlobalColumnHeights()
{
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)		//	Loop over for all element groups
    {
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();
        
		for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp
        {
            CElement& Element = ElementGrp[Ele];

            // Generate location matrix
            Element.GenerateGlobalLocationMatrix();
            
            GlobalStiffnessMatrix->CalculateColumnHeight(Element.GetGlobalLocationMatrix(), Element.GetND());
        }
    }
    
    GlobalStiffnessMatrix->CalculateMaximumHalfBandwidth();
}

//    Allocate storage for matrices Force, ColumnHeights, DiagonalAddress and StiffnessMatrix
//    and calculate the column heights and address of diagonal elements
void CDomain::AllocateMatrices()
{
    //    Allocate for global force/displacement vector
    Force = new double[NEQ];
    
    //  Create the banded stiffness matrix
    StiffnessMatrix = new CSkylineMatrix<double>(NEQ);
    
    //    Calculate column heights
    CalculateColumnHeights();
    
    //    Calculate address of diagonal elements in banded matrix
    StiffnessMatrix->CalculateDiagnoalAddress();
    
    //    Allocate for banded global stiffness matrix
    StiffnessMatrix->Allocate();
    
    COutputter* Output = COutputter::GetInstance();
    Output->OutputTotalSystemData();
}

//	Allocate GlobalForce, ColumnHeights, DiagonalAddress and GlobalStiffnessMatrix and 
//  calculate the column heights and address of diagonal elements
void CDomain::AllocateGlobalMatrices()
{
	// Allocate gloabl force
	Force = new double[NEQ + NBC];

	// Create the banded global stiffness matrix
	GlobalStiffnessMatrix = new CSkylineMatrix<double>(NEQ + NBC);

	// Calculate column heights
	CalculateGlobalColumnHeights();

    //    Calculate address of diagonal elements in banded matrix
    GlobalStiffnessMatrix->CalculateDiagnoalAddress();
    
    //    Allocate for banded global stiffness matrix
    GlobalStiffnessMatrix->Allocate();
}

//	Assemble the banded gloabl stiffness matrix
void CDomain::AssembleStiffnessMatrix()
{
//	Loop over for all element groups
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();

		unsigned int size = ElementGrp[0].SizeOfStiffnessMatrix();
		double* Matrix = new double[size];

//		Loop over for all elements in group EleGrp
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
        {
            CElement& Element = ElementGrp[Ele];
            Element.ElementStiffness(Matrix);
            StiffnessMatrix->Assembly(Matrix, Element.GetLocationMatrix(), Element.GetND());
        }

		delete[] Matrix;
		Matrix = nullptr;
	}

#ifdef _DEBUG_
	COutputter* Output = COutputter::GetInstance();
	Output->PrintStiffnessMatrix();
#endif

}

//	Assemble global stiffness matrix
void CDomain::AssembleGlobalStiffnessMatrix()
{
//	Loop over for all element groups
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();

		unsigned int size = ElementGrp[0].SizeOfStiffnessMatrix();
		double* Matrix = new double[size];

//		Loop over for all elements in group EleGrp
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
        {
            CElement& Element = ElementGrp[Ele];
            Element.ElementStiffness(Matrix);
            GlobalStiffnessMatrix->Assembly(Matrix, Element.GetGlobalLocationMatrix(), Element.GetND());
        }

		delete[] Matrix;
		Matrix = nullptr;
	}
}

//	Assemble the global nodal force vector for load case LoadCase
bool CDomain::AssembleForce(unsigned int LoadCase)
{
	if (LoadCase > NLCASE) 
		return false;

	CLoadCaseData* LoadData = &LoadCases[LoadCase - 1];

    clear(Force, NEQ);

//	Loop over for all concentrated loads in load case LoadCase
	for (unsigned int lnum = 0; lnum < LoadData->nloads; lnum++)
	{
		unsigned int dof = NodeList[LoadData->node[lnum] - 1].bcode[LoadData->dof[lnum] - 1];
        
        if(dof) // The DOF is activated
            Force[dof - 1] += LoadData->load[lnum];
	}

	return true;
}

// Assemble the essential boundary condition vector
bool CDomain::AssmbleEssentialDisplacement()
{
	// Allocate and initialize
	EssentialDisplacement = new double[NBC];
	for(int i=0; i<NBC; i++) EssentialDisplacement[i] = 0;

	// Loop over non-zero essential boundary and assign data
	for(int i=0; i<EBdata->nebcs; i++){
		unsigned int rank = NodeList[EBdata->node[i]-1].gbcode[EBdata->dof[i]-1] - NEQ - 1;
		EssentialDisplacement[rank] = EBdata->strain[i];
	}

	return true;
}

// Substitute contribution of essential boundary of force vector: \hat{f_x} = f_x - K_{ux}d_u
bool CDomain::SubstitueEssentialBoundary()
{
	unsigned int *ColumnHeights = GlobalStiffnessMatrix->GetColumnHeights();

	// Loop in column to get best cache hit
	for(unsigned int j = NEQ; j< NEQ + NBC; j++){
		// Start from first non-zero row
		for(unsigned int i = j - ColumnHeights[j]; i< NEQ; i++){
			Force[i] -= GlobalStiffnessMatrix->operator()(i,j) * EssentialDisplacement[j-NEQ];
		}
	}

	return true;
}

// Compute nodal force (including nodes on essential boudndary) f = Kd
bool CDomain::ComputeNodalForce()
{
	// Allocate and initialize
	GlobalNodalForce = new double[NEQ + NBC];
	for(int i=0; i<NEQ+NBC; i++) GlobalNodalForce[i]=0;
	unsigned int *ColumnHeights = GlobalStiffnessMatrix->GetColumnHeights();

	// Loop in column to get best cache hit
	for(unsigned int j = 0; j < NEQ+NBC; j++){
		// Start from first non-zero row
		for(unsigned int i = j - ColumnHeights[j]; i< NEQ; i++){
			double delta;
			if(j <NEQ){
				delta = (*GlobalStiffnessMatrix)(i, j) * Force[j];
			}
			else{
				delta = (*GlobalStiffnessMatrix)(i, j) * EssentialDisplacement[j - NEQ];
			}
			GlobalNodalForce[i] += delta;
		}
	}
}