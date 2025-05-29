#pragma once

#include "Element.h"

using namespace std;

//! T3 element class
class CT3 : public CElement
{
private:
	double B[3][6];

public:

//! COnstructor
    CT3();

//! Destructor
    ~CT3();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream
	virtual void Write(COutputter& output);

//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);
};