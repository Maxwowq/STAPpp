#pragma once

#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

//! Outputer class is used to output results
class CTecplot
{
private:

//!	File stream for output
	ofstream DatOutputFile;

//!	Designed as a single instance class
	static CTecplot* _instance;

//! Constructor
    CTecplot(string FileName);

public:

//!	Return the single instance of the class
	static CTecplot* GetInstance(string FileName = " ");

//! Output of .dat file with x, y, u, v
    void OutputDatFile();

//! Output T3 zone in .dat file
	void OutputT3Zone(unsigned int EleGrp);
};