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

    ND_ = 6; // 2D element in 3D, z axis is constrained
    LocationMatrix_ = new unsigned int[ND_];

    ElementMaterial_ = nullptr;
}

//  Destructor
CT3::~CT3()
{
}

// Generate location matrix of T3 element(ignoring z axis)
void CT3::GenerateLocationMatrix()
{
    int i=0;
    for(int N=0; N < NEN_; N++){
        for(int D=0; D < 2; D++){
            LocationMatrix_[i++] = nodes_[N]->bcode[D];
        }
    }
}

// Generate global location matrix of T3 element (ignoring z axis)
void CT3::GenerateGlobalLocationMatrix()
{
    int i=0;
    for(int N=0; N < NEN_; N++){
        for(int D=0; D < 2; D++){
            GlobalLocationMatrix_[i++] = nodes_[N]->gbcode[D];
        }
    }
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

    GenerateLocationMatrix();

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

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CT3::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());

// Calculate elasticity matrix(finished in material reading period)
    CT3Material* material_ = dynamic_cast<CT3Material*>(ElementMaterial_);

// Calcuclate a, b, c and element area A
    double a[3], b[3], c[3];

    for(int i=0; i<3; i++){
        int j = (i+1)%3;
        int k = (i+2)%3;
        a[i] = nodes_[j]->XYZ[0]*nodes_[k]->XYZ[1] - nodes_[k]->XYZ[0]*nodes_[j]->XYZ[1];
        b[i] = nodes_[j]->XYZ[1] - nodes_[k]->XYZ[1];
        c[i] = nodes_[k]->XYZ[0] - nodes_[j]->XYZ[0];
    }

    double A = (a[0] + a[1] + a[2])/2;

// Calculate B Matrix
    for (int i = 0; i < 3; i++) {
        int node_idx = i;
        int col_u = 2 * node_idx;
        int col_v = 2 * node_idx + 1;

        B[0][col_u] = b[i] / (2 * A); // ε_xx
        B[0][col_v] = 0;

        B[1][col_u] = 0;
        B[1][col_v] = c[i] / (2 * A); // ε_yy

        B[2][col_u] = c[i] / (2 * A); // γ_xy
        B[2][col_v] = b[i] / (2 * A);
    }

// Calculate DB[3][6] matrix
    double DB[3][6];
    for(int i=0; i<3; i++){
        for(int j=0; j<6; j++){
            DB[i][j] = 0;
            for(int k=0; k<3; k++){
                DB[i][j] += material_->D[i][k] * B[k][j];
            }
        }
    }

// Calculate final K^e (Matrix)
    for(int i=0; i<6; i++){
        // only compute upper half
        for(int j=i; j<6; j++){
            // convert i,j to matrix rank
            int rank = (j+1)*(j+2)/2 - (i+1);
            Matrix[rank] = 0;
            for(int k=0; k<3; k++){
                Matrix[rank] += B[k][i] * DB[k][j];
            }
            Matrix[rank] = Matrix[rank] * A * material_->t;
        }
    }
}

//	Calculate element stress
void CT3::ElementStress(double* stress, double* Displacement)
{
    CT3Material* material_ = dynamic_cast<CT3Material*>(ElementMaterial_);

    double epsilon[3];

    // epsilon = B * d
    for(int i=0; i<3; i++){
        epsilon[i] = 0;
        for(int k=0; k<6; k++){
            if(LocationMatrix_[k]){
                epsilon[i] += B[i][k] * Displacement[LocationMatrix_[k]-1];
            }
        }
    }

    // s = D * epsilon
    for(int i=0; i<3; i++){
        stress[i] = 0;
        for(int k=0; k<3; k++){
            stress[i] += material_->D[i][k] * epsilon[k];
        }
    }

    if(material_->PlaneStress){
        stress[3] = 0;
    }
    else{
        // s22 = lambda * (epsilon_x + epsilon_y)
        double E = material_->E;
        double nu = material_->nu;
        double lambda = E*nu/(1+nu)/(1-2*nu);
        stress[3] = lambda * (epsilon[0] + epsilon[1]);
    }
}