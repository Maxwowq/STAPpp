#pragma once

#include "Outputter.h"

using namespace std;

//! Class EBData is used to store non-zero essential boundary
class CEBData
{
public:
    unsigned int nebcs; //!< Number of non zero essential boudary conditions
    unsigned int* node; //!< Node number to which this ebc is applied
    unsigned int* dof;  //!< Degree of freedom number for this ebc
    double* strain;     //!< Magnitude of displacement

public:
    CEBData() : nebcs(0), node(NULL), dof(NULL), strain(NULL) {};
    ~CEBData();

//!	Set nebcs, and new array node, dof and strain
	void Allocate(unsigned int num);

//!	Read ebc data from stream Input
	bool Read(ifstream& Input);

//!	Write ebc data to stream
	void Write(COutputter& output);
};