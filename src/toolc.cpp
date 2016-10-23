//============================================================================
// Name        : toolc.cpp
// Author      : Christoph Loschen
// Version     : v1.0
// Copyright   : GNU Public licence
// Description : Pair Potential Global Optimisation & Surface Constructor
//============================================================================
//TO DO:
// find residue error in shell modell, in long range part...
//Increase performance of code
//Implement Ewald for VdW
//introduce strip O shells and spawn O shells function!!!!
//Optimierungspotential: Individual cutoffs, less steps optimizer, pol=false setup frueher
//Shell model for periodic approach
//Introduce Grid & Grid interface
//Predict ForceConstants with NN
//Predict surface reconstruction with basin hopping
//bond order potentials?
#include <sys/time.h>
//#include <laprefs.h>
#include <iostream>
//#include <omp.h>
#include "toolc.h"
#include "ToolCalc.h"
#include "ToolCalcPP.h"
#include "ToolCalcLJ.h"
#include "ToolIO.h"

using namespace std;

int main() {
	//Ran myRan(188);
	cout << "###TOOL++ V1.1, C. Loschen 2008###\n";
	ToolIO *newIO = new ToolIO;

	//2x parsen...
	(*newIO).parseCalctype();
	ToolCalc *myECalc;
	if (toolc::ctype == "PP") {
		myECalc = new ToolCalcPP();
		(*myECalc).maxiter = 300;
		(*myECalc).ewald_thresh = 0.01;
		(*myECalc).threshhold = 0.01;
		(*myECalc).max_displacement = 0.5;
		(*myECalc).step = 0.1;
		(*myECalc).cutoff = 35;
		(*myECalc).dipole_cor = true;

	} else {
		myECalc = new ToolCalcLJ;
		//(*myECalc).max_displacement = 0.1;
		(*myECalc).step = 0.002;
	}
	(*newIO).parseSETUP(*myECalc);

	if (toolc::ctype == "PP") {
		(*myECalc).setDipole(true);
	}
	Ran myRan((*myECalc).rseed);
	//new better timing
	timeval t1, t2;
	gettimeofday(&t1, NULL);
//	timespec res, start, stop;
//	clock_getres(CLOCK_REALTIME, &res);
//	clock_gettime(CLOCK_REALTIME, &start);
	//#############################

	//######Global OPT##########

	(*myECalc).moveRandom(myRan, .5);

	//(*myECalc).periodic = true;
	(*myECalc).setUP_periodic((*myECalc).periodic, false, true);
	//(*myECalc).E_periodic((*myECalc).shell_modell,true,true);
	//(*newIO).fractoFile((*myECalc), false, true);
	//(*newIO).printMol((*myECalc),(*myECalc).shell_modell);
	//WARNING IF FROZEN ATOMS
	//(*myECalc).xyz2grid(16,16,16,true);
	//(*newIO).gridout((*myECalc).genom,(*myECalc).gridp);

	//(*myECalc).grid2genom(true);
	//
	//
	//(*myECalc).orderOx();

	//(*myECalc).monte_carlo_sampling(myRan, 2);
	//(*newIO).genomout((*myECalc).energy,(*myECalc).genom_compressed,(*myECalc).atnumber);
	//(*myECalc).grid2xyz(true);
	(*myECalc).basin_hopping(myRan, 0.5, false);

	//(*myECalc).shell_modell = false;
	//(*myECalc).E((*myECalc).shell_modell, true, true);
	(*myECalc).opt((*myECalc).shell_modell, true);
	//(*myECalc).shell_modell=true;

	//(*myECalc).basin_hopping(myRan, 0.5, false);
	//(*myECalc).append_shells(false);
	//(*myECalc).shell_modell=true;
	//(*myECalc).setcoreshell(false,true);
	//(*myECalc).opt((*myECalc).shell_modell, true);
	//(*myECalc).strip_shells(false);
	//(*myECalc).opt((*myECalc).shell_modell, true);
	//(*myECalc).E_periodic((*myECalc).shell_modell,true,true);
	//(*myECalc).E_periodic_Pol((*myECalc).shell_modell,true,true);
	//#############################
	//(*newIO).printMol((*myECalc),(*myECalc).shell_modell);
	//(*newIO).printGrad((*myECalc));
	//(*newIO).fractoFile((*myECalc), false, true);
	gettimeofday(&t2, NULL);
	cout << "\n\nTIMING: ";
	(*newIO).printTiming(t1, t2);
	(*newIO).printParameters(*myECalc);
	//(*newIO).gridtoFile((*myECalc).energy, (*myECalc).genom, (*myECalc).gridp, false);
	//(*newIO).printMol((*myECalc), true);
	//(*newIO).moltoFile((*myECalc));
	//(*newIO).printGrad((*myECalc));
	//cout <<ToolCalc::erfc_new(2.0)<< endl;
	delete newIO;
	delete myECalc;
}
//Konstruktor
toolc::toolc() {
}

string toolc::ctype = "PP";

//Destruktor, has to be defined
toolc::~toolc() {
	//delete *xyz;
}

