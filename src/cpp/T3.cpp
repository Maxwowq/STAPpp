#include "T3.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//  Constructor
CT3::CT3()
{
    NEN_ = 3; // Each element has 3 nodes
    nodes_ = new CNode*[NEN_];

    ND_ = 9;
    LocationMatrix_ = new unsigned int[ND_];

    ElementMaterial_ = nullptr;
}

//  Destructor
CT3::~CT3()
{
}

// Read element data from stram Input
bool CT3::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int MSet;  // Material property set number
    unsigned int N1, N2, N3;    // T3 element has three nodes that are counterclockwise numbered

    Input >> N1 >> N2 >> N3 >> MSet;
    ElementMaterial_ = dynamic_cast<CT3Material*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
    nodes_[2] = &NodeList[N3 - 1];

    return true;
}

//	Write element data to stream
void CT3::Write(COutputter& output)
{
    output  << setw(11) << nodes_[0]->NodeNumber 
            << setw(9)  << nodes_[1]->NodeNumber 
            << setw(9)  << nodes_[2]->NodeNumber 
            << setw(12) << ElementMaterial_->nset 
            << endl;
}

// Calculate elasticity matrix


//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CT3::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());

// Calculate elasticity matrix


// Calcuclate a, b, c and element area A
}