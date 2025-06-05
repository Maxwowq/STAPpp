#include "EssentialBoudary.h"

#include <iomanip>
#include <iostream>

using namespace std;

CEBData :: ~CEBData()
{
    delete [] node;
    delete [] dof;
    delete [] strain;
}

void CEBData::Allocate(unsigned int num)
{
    nebcs = num;
    node = new unsigned int[nebcs];
    dof = new unsigned int[nebcs];
    strain = new double[nebcs];
}

// Read essential boudary data from stream Input
bool CEBData::Read(ifstream& Input)
{
    unsigned int NEBC;

    Input >> NEBC;

    Allocate(NEBC);

    for (unsigned int i=0; i<NEBC; i++)
        Input >> node[i] >> dof[i] >> strain[i];

    return true;
}

//	Write essential boudary data data to stream
void CEBData::Write(COutputter& output)
{
	for (unsigned int i = 0; i < nebcs; i++)
		output << setw(7) << node[i] << setw(13) << dof[i]  << setw(19) << strain[i] << endl;
}
