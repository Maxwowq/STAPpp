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