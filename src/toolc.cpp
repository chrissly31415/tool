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
#include <iostream>
#include "toolc.h"
#include "ToolCalc.h"
#include "ToolCalcPP.h"
#include "ToolCalcLJ.h"
#include "ToolCalcZE.h"
#include "ToolIO.h"

using namespace std;

int main(int argc, char *argv[]) {

	timeval t1, t2;
	gettimeofday(&t1, NULL);
	cout << "###TOOL++ V2.0, C. Loschen ###\n";
	ToolIO *newIO = new ToolIO;

	const char* setupfile = "SETUP";
	if (argc > 1) {
		setupfile = argv[1];
	}

	newIO->parseCalctype(setupfile);
	ToolCalc *myECalc;

	//some defaults for pair potentials
	if (toolc::ctype == "PP") {
		myECalc = new ToolCalcPP();
		myECalc->maxiter = 300;
		myECalc->ewald_thresh = 0.01;
		myECalc->threshhold = 0.01;
		myECalc->max_displacement = 0.5;
		myECalc->step = 0.1;
		myECalc->cutoff = 35;
		myECalc->dipole_cor = true;
	} else if (toolc::ctype == "ZERNICKE") {
		myECalc = new ToolCalcZE();

	} else {
		myECalc = new ToolCalcLJ;
		//(*myECalc).max_displacement = 0.1;
		myECalc->step = 0.002;
	}
	//read geometry
	newIO->parseSETUP(*myECalc, setupfile);
	Ran myRan((*myECalc).rseed);

	if (toolc::ctype == "PP") {
		myECalc->setDipole(true);
		myECalc->setUP_periodic((*myECalc).periodic, false, true);
		myECalc->moveRandom(myRan, .5);
		myECalc->basin_hopping(myRan, 0.5, false);
		myECalc->opt((*myECalc).shell_modell, true);
		//(*myECalc).periodic = true;
		//(*myECalc).E_periodic((*myECalc).shell_modell,true,true);
	} else if (toolc::ctype == "LJ") {
		myECalc->moveRandom(myRan, .5);
		myECalc->basin_hopping(myRan, 0.5, false);
		//(*myECalc).monte_carlo_sampling(myRan, 2);
	} else if (toolc::ctype == "ZERNICKE") {
		ToolCalcZE* zeCalc = dynamic_cast<ToolCalcZE*>(myECalc);
		zeCalc->printSegments();
		zeCalc->seg2grid(16,16,16,true);

		//(*newIO).gridout((*myECalc).genom,(*myECalc).gridp);

		//(*myECalc).grid2genom(true);
		//
		//
		//(*myECalc).orderOx();

		//(*newIO).genomout((*myECalc).energy,(*myECalc).genom_compressed,(*myECalc).atnumber);
		//(*myECalc).grid2xyz(true);
	}

	//(*newIO).fractoFile((*myECalc), false, true);
	//(*newIO).printMol((*myECalc),(*myECalc).shell_modell);
	//WARNING IF FROZEN ATOMS

	//(*myECalc).shell_modell = false;
	//(*myECalc).E((*myECalc).shell_modell, true, true);

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
	newIO->printParameters(*myECalc);
	gettimeofday(&t2, NULL);
	cout << "\n\nTIMING: ";
	newIO->printTiming(t1, t2);

	//(*newIO).gridtoFile((*myECalc).energy, (*myECalc).genom, (*myECalc).gridp, false);
	//(*newIO).printMol((*myECalc), true);
	//(*newIO).moltoFile((*myECalc));
	//(*newIO).printGrad((*myECalc));

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

