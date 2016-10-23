#ifndef TOOLIO_H_
#define TOOLIO_H_

#include "ToolCalc.h"
//#include <lapackpp.h>
#include <Eigen/Core>

using namespace std;

class ToolIO
{
public:
	ToolIO();
	virtual ~ToolIO();
	void parseCalctype();
	void parseSETUP(ToolCalc &A);

	void printMol(ToolCalc &A, bool shel=false);
	void printFrac(ToolCalc &a, bool shel=false, bool peridoic=false);
	void fractoFile(ToolCalc &a, bool shel=false, bool peridoic=false);
	void printGrad(ToolCalc &A);
	void moltoFile(ToolCalc &A);
	void moltoFile2(int &atnumber, double &energy, string atoms[], ToolCalc::VectorXd &xyz,  string core_shell[]=NULL, int nr_shells=0, bool shel=false, string type="LJ");
	void fractoFile2(int &atnumber, double &energy, string atoms[], ToolCalc::VectorXd &xyz,  string core_shell[]=NULL, int nr_shells=0, bool shel=false, double Lx=10, double Ly=10, double Lz=10);
	void printRDF(int &atnumber, double &energy, string atoms[], ToolCalc::VectorXd &xyz,double cutoff=25, double Lx=10, double Ly=10, double Lz=10);
	void datatoFile(double energy, double dipnorm, bool erase=false);
	void genom2File(double lenergy, int lgenom_compressed[], int atnumber, bool erase=false);
	void gridout(int lgenom[],int lgridp);
	void genomout(double lenergy, int lgenom_compressed[], int atnumber);
	void printParameters(ToolCalc &A);
	void printTiming(timeval &start, timeval &end);
};

#endif /*TOOLIO_H_*/
