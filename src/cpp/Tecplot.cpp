#include "Domain.h"
#include "Tecplot.h"

using namespace std;

//	Constructor
CTecplot::CTecplot(string FileName)
{
    string datFileName = FileName + "_tecplot.dat";
    DatOutputFile.open(datFileName);

    if(!DatOutputFile)
    {
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
    }
}

CTecplot* CTecplot::_instance = nullptr;

//	Return the single instance of the class
CTecplot* CTecplot::GetInstance(string FileName)
{
    if(!_instance){
        _instance = new CTecplot(FileName);
    }

    return _instance;
}

// Output filename_tecplot.dat file
void CTecplot::OutputDatFile()
{
    CDomain* FEMData = CDomain::GetInstance();
    
    // output the title line
    DatOutputFile << "Title=\"" << FEMData->GetTitle() << "\"" <<endl;

    // Output the variables
    DatOutputFile << "Variables=\"X\"   \"Y\"   \"U\"   \"V\"   \"FX\"   \"FY\" " << endl;

    unsigned int NUMEG = FEMData->GetNUMEG();

    for(unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++){
        ElementTypes ElementType = FEMData->GetEleGrpList()[EleGrp].GetElementType();
        switch (ElementType)
        {
            case ElementTypes::T3:
                OutputT3Zone(EleGrp);
                break;
            default:
                cout << ElementType << "has no Tecplot output implementation." << endl;
                break;
        }
    }
}

// Output T3 zone in .dat file
void CTecplot::OutputT3Zone(unsigned int EleGrp)
{
    CDomain* FEMData = CDomain::GetInstance();

    unsigned int NUMNP = FEMData->GetNUMNP();
    unsigned int NUME = FEMData->GetEleGrpList()[EleGrp].GetNUME();

    // Output Zone info
    DatOutputFile << "Zone N=" << NUMNP << " E=" << NUME 
                  << " F=FEPOINT ET=TRIANGLE" << endl;

    // Output nodal info
    CNode* Nodelist = FEMData->GetNodeList();
    double* Displacement = FEMData->GetGlobalDisplacement();
    double* NodalForce = FEMData->GetNodalForce();
    for(unsigned int np = 0; np < NUMNP; np++){
        // Output X, Y
        DatOutputFile << Nodelist[np].XYZ[0] << setw(20) <<Nodelist[np].XYZ[1] << setw(20);

        // Output U, V
        double u=0, v=0;
        if(Nodelist[np].gbcode[0]){
            u = Displacement[Nodelist[np].gbcode[0] - 1];
        }
        if(Nodelist[np].gbcode[1]){
            v = Displacement[Nodelist[np].gbcode[1] -1];
        }
        DatOutputFile << u << setw(20) << v << setw(20);

        // Output FX, FY
        double fx=0, fy=0;
        if(Nodelist[np].gbcode[0]){
            fx = NodalForce[Nodelist[np].gbcode[0] - 1];
        }
        if(Nodelist[np].gbcode[1]){
            fy = NodalForce[Nodelist[np].gbcode[1] -1];
        }
        DatOutputFile << fx << setw(20) << fy << endl;
    }

    // Output element info
    CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
    for(unsigned int Ele=0; Ele < NUME; Ele++){
        unsigned int NodeNumber[3];
        for(int i=0; i<3; i++){
            NodeNumber[i] = ElementGroup[Ele].GetNodes()[i]->NodeNumber;
        }
        DatOutputFile << NodeNumber[0] << setw(20) << NodeNumber[1] << setw(20) << NodeNumber[2] << endl;
    }
}