#include "ToolCalc.h"
#include "ToolCalcPP.h"
#include "ToolIO.h"

//#include <omp.h>

using namespace std;

//Constructor with variable initialization
ToolCalc::ToolCalc() {
	atnumber = 0;
	calculated = false;
	energy = 0.;
	gradnorm = 0.;
	gradinfnorm = 0.;
	converged = false;
	//periodic variables
	periodic = false;
	Lx = 100., Ly = 100., Lz = 100.;
	alpha = 0.2;
	//general
	maxatoms = 10000;
	verbose = false;
	nproc = 1;
	shell_modell = false;
	scutoff = 0.6;
	rseed = 123;
	rotation_angle = 0.0;
	rotation_axis = -1;
	voxelstep = 0.1;
	voxelmode=0;

}
//PES
string ToolCalc::elements[87] = { "X", "H", "He", "Li", "Be", "B", "C", "N",
		"O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
		"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
		"As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
		"Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba",
		"La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
		"Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
		"Tl", "Pb", "Bi", "Po", "At", "Rn" };

//atomic masses
double ToolCalc::emass[87] = { 0.0, 1.0079, 4.00260, 6.941, 9.01218, 10.81,
		12.011, 14.0067, 16.0, 18.998403, 20.1797, 22.98977, 24.305, 26.98154,
		28.0855, 30.97376, 32.066, 35.4527, 39.948, 39.0983, 40.078, 44.9559,
		47.867, 50.9415, 51.9961, 54.9380, 55.845, 58.9332, 58.6934, 63.546,
		65.39, 69.723, 72.61, 74.9216, 78.96, 79.904, 83.80, 85.4678, 87.62,
		88.9059, 91.224, 92.9064, 95.94, 98.0, 101.07, 102.9055, 106.42,
		107.868, 112.41, 114.818, 118.69, 121.760, 127.60, 126.9045, 131.29,
		132.9054, 137.327, 138.9055, 140.12, 140.9077, 144.24, 145.0, 150.36,
		151.964, 157.25, 158.9254, 162.50, 164.9303, 167.26, 168.9342, 173.04,
		174.967, 178.49, 180.9479, 183.85, 186.207, 190.23, 192.217, 195.078,
		196.9665, 200.59, 204.3883, 207.2, 208.9804, 209.0, 210.0, 222.0 };

//makes random move
void ToolCalc::moveRandom(Ran &myRan, double stepsize, double zscale) {
	double j;
	double z;
	//double stepsize =0.3;
	for (int i = 0; i < atnumber * 3; i++) {
		j = myRan.doub();
		j = (2 * j - 1.0) * stepsize * (*moveMat)(i);
		//cout << j <<endl;
		(*xyz)(i) = ((*xyz)(i) + j);
		if (i % 3 == 2) {
			//this is z-coordinate
			z = myRan.doub();
			z = (2 * z - 1.0) * zscale * (*moveMat)(i);
			(*xyz)(i) = ((*xyz)(i) + z);
		}
	}
}

void ToolCalc::moveRandomPol(Ran &myRan, double stepsize) {
	//cores should be positioned at the end of the input
	double j;
	int atom = 0;
	//double stepsize =0.3;
	for (int i = 0; i < atnumber * 3; i++) {
		atom = i / 3;
		j = myRan.doub();
		j = (2 * j - 1.0) * stepsize * (*moveMat)(i);
		(*xyz)(i) = ((*xyz)(i) + j);
		if (core_shell[atom] == "SHEL") {
			(*xyz)(i + nr_cores * 3) = (*xyz)(i);
		}
		//cout << j <<endl;

	}
}

//creates random structures
void ToolCalc::createRandom(Ran &myRan, double radius) {
	double j;
	for (int i = 0; i < atnumber * 3; i++) {
		j = myRan.doub();
		j = (2 * j - 1.0) * radius;
		//cout << j <<endl;
		//complicated expression: moveMat=0: Partikel frozen, moveMat=1: can be moved / scaled
		(*xyz)(i) = (*xyz)(i) * (1 - (*moveMat)(i)) + j * (*moveMat)(i);
	}
}

//rotate structure
void ToolCalc::rotateCoordinates(int axis, double phi) {
	cout << "Rotation of coordinates around axis:" << axis << " using phi="
			<< phi << "\n";

	phi = phi * 3.14159265 / 180.0;
	Eigen::MatrixXd R(3, 3);
	if (axis == 0) {
		R << 1.0, 0.0, 0.0, 0.0, cos(phi), -sin(phi), 0.0, sin(phi), cos(phi);
	} else if (axis == 1) {
		R << cos(phi), 0.0, sin(phi), 0.0, 1.0, 0.0, -sin(phi), 0.0, cos(phi);
	} else if (axis == 2) {
		R << cos(phi), -sin(phi), 0.0, sin(phi), cos(phi), 0.0, 0.0, 0.0, 1.0;
	}
	//rotate atoms
	double x, y, z;
	for (int i = 0; i < atnumber; ++i) {
		x = (*xyz)(3 * i);
		y = (*xyz)(3 * i + 1);
		z = (*xyz)(3 * i + 2);
		(*xyz)(3 * i) = R(0, 0) * x + R(0, 1) * y + R(0, 2) * z;
		(*xyz)(3 * i + 1) = R(1, 0) * x + R(1, 1) * y + R(1, 2) * z;
		(*xyz)(3 * i + 2) = R(2, 0) * x + R(2, 1) * y + R(2, 2) * z;
	}
}

//sets xyz-coordinates
void ToolCalc::setXYZ(VectorXd xyz_new) {
	for (int i = 0; i < atnumber * 3; i++) {
		(*xyz)(i) = xyz_new(i);
	}
}

//includes all 1x tasks
void ToolCalc::setUP_periodic(bool periodic, bool pol, bool verbose) {

	if (verbose == true) {
		cout << endl;
		cout << "Setup of periodic calculation..." << endl;
	}

	if (periodic != true) {
		cout << "Should set the keyword >periodic<..." << endl;
		return;
	}
	double LLx = 0, LLy = 0, LLz = 0;
	//V=LLx*LLy*LLz;
	LLx = Lx * ANG2BOHR;
	LLy = Ly * ANG2BOHR;
	LLz = Lz * ANG2BOHR;
	V = LLx * LLy * LLz;
	double tr2tf = 3.6;
	if (verbose == true) {
		cout << "Cell volume: " << V << " bohr³" << endl << "Optimum alpha: "
				<< pow((tr2tf * pow(PI, 3) * atnumber / pow(V, 2)), 0.16667)
				<< endl;
	}
	double tmp;
	//set up of real and reciprocal lattice vectors
	//real lattice vectors
	ar[0] = LLx;
	ar[1] = 0;
	ar[2] = 0;
	br[0] = 0;
	br[1] = LLy;
	br[2] = 0;
	cr[0] = 0;
	cr[1] = 0;
	cr[2] = LLz;
	//reciprocal lattice vectors
	//a*
	af[0] = br[1] * cr[2] - br[2] * cr[1];
	af[1] = br[2] * cr[0] - br[0] * cr[2];
	af[2] = br[0] * cr[1] - br[2] * cr[0];
	tmp = ar[0] * af[0] + ar[1] * af[1] + ar[2] * af[2];
	af[0] = af[0] / tmp;
	af[1] = af[1] / tmp;
	af[2] = af[2] / tmp;
	//b*
	bf[0] = ar[1] * cr[2] - ar[2] * cr[1];
	bf[1] = ar[2] * cr[0] - ar[0] * cr[2];
	bf[2] = ar[0] * cr[1] - ar[2] * cr[0];
	tmp = br[0] * bf[0] + br[1] * bf[1] + br[2] * bf[2];
	bf[0] = bf[0] / tmp;
	bf[1] = bf[1] / tmp;
	bf[2] = bf[2] / tmp;
	//c*
	cf[0] = ar[1] * br[2] - ar[2] * br[1];
	cf[1] = ar[2] * br[0] - ar[0] * br[2];
	cf[2] = ar[0] * br[1] - ar[2] * br[0];
	tmp = cr[0] * cf[0] + cr[1] * cf[1] + cr[2] * cf[2];
	cf[0] = cf[0] / tmp;
	cf[1] = cf[1] / tmp;
	cf[2] = cf[2] / tmp;

	//project atoms into cell
	if (verbose == true) {
		cout << "Shifting " << atnumber << " atoms into the cell..." << endl;
	}
	for (int i = 0; i < atnumber * 3; i++) {
		////			//			cout << "!";
		if (i % 3 == 0) {
			if ((*xyz)(i) < 0.0) {
				(*xyz)(i) = (*xyz)(i) + Lx;
			}
			if ((*xyz)(i) > Lx) {
				(*xyz)(i) = (*xyz)(i) - Lx;
			}
		}
		if (i % 3 == 1) {
			if ((*xyz)(i) < 0.0) {
				(*xyz)(i) = (*xyz)(i) + Ly;
			}
			if ((*xyz)(i) > Ly) {
				(*xyz)(i) = (*xyz)(i) - Ly;
			}
		}
		if (i % 3 == 2) {
			if ((*xyz)(i) < 0.0) {
				(*xyz)(i) = (*xyz)(i) + Lz;
			}
			if ((*xyz)(i) > Lz) {
				(*xyz)(i) = (*xyz)(i) - Lz;
			}
		}
	}
}

//Lennard Jones Energy optimisation
void ToolCalc::monte_carlo_sampling(Ran &myRan, double distortion) {
	ToolIO newIO;
	VectorXd xyz_min(3 * atnumber);
	double min_energy;
	VectorXd xyz_global(3 * atnumber);
	double global_energy = 0;
	VectorXd xyz_old(3 * atnumber);
	double energy_old = 0;
	double k = 8.617343 / 10000.0; //eV/K
	double random = 0;
	double boltzf = 0;
	cout << "###Monte Carlo Sampling at T=" << temperature << " K##" << endl;
	//outer loop, creating new random structures
	for (int j = 0; j < nr_reset; j++) {
		cout << "###Creating Random Distortion, iteration " << (j + 1)
				<< " out of " << nr_reset << " ###" << endl;
		moveRandom(myRan, distortion);
		opt();
		cout << endl;
		//our minimum is our first structure
		min_energy = energy;
		//xyz_min.copy(*xyz);
		xyz_min = xyz->replicate(1, 1);
		if (j == 0) {
			global_energy = energy;
			xyz_global = xyz->replicate(1, 1);
		}
		//inner loop
		for (int i = 0; i < global_maxiter; i++) {
			cout << "MC Move " << i << ": ";
			xyz_old = xyz->replicate(1, 1);
			energy_old = energy;
			moveRandom(myRan, distortion);
			opt(false, false);

			//project structure to grid
			setUP_periodic(periodic, false, false);
			///xyz2grid(16, 16, 16, false);
			///grid2genom(false);
			//cout<<endl;
			//newIO.genomout(energy, genom_compressed, atnumber);
			newIO.printRDF(atnumber, energy, atoms, *xyz, cutoff, Lx, Ly, Lz);
			newIO.genom2File(energy, genom_compressed, atnumber);
			//higher energy
			if (energy > (energy_old)) {
				//test for Boltzmann Faktor
				random = myRan.doub();
				boltzf = exp(
						-1 * ((energy - energy_old) * HARTREE2EV)
								/ (k * temperature));
				//				cout << energy << " " << min_energy << endl;
				//
				//				cout << "BF: " << boltzf << endl;
				//				cout << "Random: " << random << endl;
				//check acceptance criterion
				if (boltzf > random) {
					//accept step
					cout << " -";
					min_energy = energy;
				} else {
					//Do not accept MC step
					cout << " #";
					*xyz = xyz_old.replicate(1, 1);
				}
				//lower energy anyway
			} else {
				if (energy < (min_energy - 0.01)) {
					min_energy = energy;
					xyz_min = xyz->replicate(1, 1);
					//cout << energy << " " << min_energy << endl;
					cout << " +";
				} else {
					//no significant change
					cout << " 0";
				}
				if (energy < (global_energy - 0.01)) {
					cout << " GM";
					global_energy = energy;
					xyz_global = xyz->replicate(1, 1);
					min_energy = energy;
					xyz_min = xyz->replicate(1, 1);
					E();
					newIO.fractoFile2(atnumber, energy, atoms, *xyz, core_shell,
							nr_shells, true, Lx, Ly, Lz);

				}
			}
			cout << endl;
		}
	}
	*xyz = xyz_global.replicate(1, 1);
}
//new try for basin hoppin
void ToolCalc::basin_hopping(Ran &myRan, double distortion, double verbose) {
	ToolIO newIO;
	VectorXd xyz_min(3 * atnumber);
	VectorXd xyz_global(3 * atnumber);
	VectorXd xyz_actual(3 * atnumber);
	VectorXd xyz_opt(3 * atnumber);
	double min_energy;
	double global_energy = 0;
	double random;
	double boltzf;
	//Qpart=1.0;

	cout << endl << "###Basin Hopping Global Optimization with T="
			<< temperature << " K##" << endl;
	//opt();
	global_energy = energy;
	xyz_global = xyz->replicate(1, 1);
	//outer loop, creating new random structures
	for (int j = 0; j < nr_reset; j++) {
		//
		if (j == 0) {
			cout << endl << "###Creating Random Structure, scale: "
					<< distortion << ", zscale: " << distortion << " iteration "
					<< (j + 1) << " out of " << nr_reset << " ###" << endl;
			//createRandom(myRan, distortion);
			moveRandom(myRan, distortion);
		} else {
			cout << endl << "###Creating Random Distortion, scale: "
					<< distortion << ", zscale: " << distortion * 3
					<< " iteration " << (j + 1) << " out of " << nr_reset
					<< " ###" << endl;
			moveRandom(myRan, distortion, distortion * 3);
		}
		//new, create VdW cluster first
		//createRandom(myRan, 0.1);
		//periodic=false;
		//		cout <<"### Refining Random Structure by Cluster Optimisation"
		//				<< " ###" <<endl;
		opt(false, false);
		//The actual minimum as a base
		min_energy = energy;
		xyz_min = xyz->replicate(1, 1);
		//our minimum is our first structure
		//inner loop
		cout << endl;
		for (int i = 0; i < global_maxiter; i++) {
			//Important, otherwise code slows down!
			if (periodic == true) {
				setUP_periodic(periodic, false, false);
			}
			//timing
			timeval t1, t2;
			gettimeofday(&t1, NULL);
			cout << "MC Move " << i << ": ";
			moveRandom(myRan, distortion);
			E();
			//save random starting structure
			xyz_actual = xyz->replicate(1, 1);

			opt(false, false);
			//save optimized structure
			xyz_opt = xyz->replicate(1, 1);
			//cout << energy << " "<<(k+temperature);
			//4 cases:
			//GM: New global minimum
			//- : Accepted upwards MC step
			//# : rejected MC step
			//+ : new minimum
			//0 : No significant change

			if (energy > (min_energy + 0.01)) {
				//Testing for Boltzmann-criterion
				random = myRan.doub();
				boltzf = exp(
						-1 * ((energy - min_energy) * HARTREE2EV)
								/ (k * temperature));
				if (boltzf > random) {
					//Accept upwards MC step
					min_energy = energy;
					xyz_min = xyz_opt.replicate(1, 1);
					cout << " -";
					//calculate Q
					//Qpart = Qpart + exp(-1* energy / (k * temperature));
					//prob = exp(-1* global_energy / (k * temperature)) / Qpart;

				} else {
					//Do not accept MC step
					setXYZ(xyz_min);
					xyz_opt = xyz_min.replicate(1, 1);
					E();
					cout << " #";

				}
			} else {
				if ((energy < (min_energy - 0.01))) {
					//New local minimum
					min_energy = energy;
					xyz_min = xyz_opt.replicate(1, 1);
					cout << " +";
					//calculate Q
					//Qpart = Qpart + exp(-1* energy / (k * temperature));
					//prob = exp(-1* global_energy / (k * temperature)) / Qpart;
				} else {
					//no significant change
					setXYZ(xyz_actual);
					//(*xyz).copy(xyz_actual);
					cout << " 0";
				}

			}
			if (energy < (global_energy - 0.01)) {
				//If new global minimum
				global_energy = energy;
				xyz_global = xyz_opt.replicate(1, 1);
				cout << " GM";
				//cout << '\a';
				newIO.fractoFile2(atnumber, energy, atoms, xyz_global,
						core_shell, nr_shells, true, Lx, Ly, Lz);
				newIO.moltoFile2(atnumber, energy, atoms, xyz_global,
						core_shell, nr_shells, shell_modell);

			}
			//cout << "\tQpart:" << Qpart;
			//Helmholtz free energy
			//F = -1* k * temperature * log (Qpart);
			//cout << setiosflags(ios::fixed) << setprecision(2)<< "\tF:" << F;
			//cout<<"\tQ:"<<Qpart;
			cout << endl;
			setDipole(false);
			E();
			newIO.datatoFile(min_energy, dipnorm, false);
			//			newIO.fractoFile2(atnumber, energy, atoms, xyz_opt, core_shell,
			//					nr_shells, true, Lx, Ly, Lz);
			//setXYZ(xyz_actual);
			//(*xyz).copy(xyz_actual);
			gettimeofday(&t2, NULL);
			if (verbose == true) {
				cout << "TIMING: ";
				newIO.printTiming(t1, t2);
			}
		}
		(*xyz) = xyz_global.replicate(1, 1);

		//opt();
	}
	if (verbose == true) {
		cout << setiosflags(ios::fixed) << setprecision(2) << "\nE(global_min):"
				<< global_energy << endl;
	}
}

//get center of mass
void ToolCalc::setCOM() {
	VectorXd mxyz(3);
	double mass = 0;
	double total_mass = 0;
	for (int i = 0; i < atnumber; i++) {
		mass = ToolCalc::emass[atom_nr[i]];
		total_mass = +mass;
		mxyz(0) = (*xyz)(3 * i) * mass;
		mxyz(1) = (*xyz)(3 * i + 1) * mass;
		mxyz(2) = (*xyz)(3 * i + 2) * mass;
	}
	mxyz(0) = mxyz(0) / total_mass;
	mxyz(1) = mxyz(1) / total_mass;
	mxyz(2) = mxyz(2) / total_mass;
	cout << setiosflags(ios::fixed) << setprecision(8);
	if (verbose == true) {
		cout << "###Center of mass: ###";
		cout << "x: " << mxyz(0) << " y: " << mxyz(1) << " z: " << mxyz(2)
				<< endl;
	}
	COM[0] = mxyz(0);
	COM[1] = mxyz(1);
	COM[2] = mxyz(2);
}
//sets gradient elemets zero which should not be moved
void ToolCalc::correctGrad() {
	for (int i = 0; i < atnumber * 3; i++) {
		(*grad)(i) = (*grad)(i) * (*moveMat)(i);
		//cout << grad [i] << " " << moveMat[i] << endl;
		//cout << (*grad)(i) << endl;
	}
}
//set dipole moment
void ToolCalc::setDipole(bool verbose) {
	//double dist;
	dipole[0] = 0;
	dipole[1] = 0;
	dipole[2] = 0;
	dipnorm = 0;
	for (int i = 0; i < atnumber * 3; i += 3) {
		dipole[0] = dipole[0] + (*xyz)(i) * q[i / 3];
		dipole[1] = dipole[1] + (*xyz)(i + 1) * q[i / 3];
		dipole[2] = dipole[2] + (*xyz)(i + 2) * q[i / 3];
	}
	if (verbose == true) {
		cout << endl << "Dipole moment: " << dipole[0] << " " << dipole[1]
				<< " " << dipole[2] << " e.angs" << endl;
	}
	dipnorm = sqrt(pow(dipole[0], 2) + pow(dipole[1], 2) + pow(dipole[2], 2));
	if (verbose == true) {
		cout << "||Dipole moment||: " << dipnorm << endl << endl;
	}

}

//used to get interaction during energy calculation
int ToolCalc::getInteraction(int type1, int type2) {
	int interact = 0;
	for (int i = 0; i < interactions; i++) {
		if ((interMat[i][0] == type1) && (interMat[i][1] == type2)) {
			interact = i;
		}
		if ((interMat[i][0] == type2) && (interMat[i][1] == type1)) {
			interact = i;
		}
	}
	return interact;
}
//Validity check: have all shells their corresponding cores?
void ToolCalc::setcoreshell(bool verbose, bool reset) {
	cout << "\nValidiy check for Catlow shell modell, spring cutoff: "
			<< scutoff << endl;
	int numpos = 3 * (atnumber - nr_cores);
	int numposall = atnumber * 3;
	double b[3];
	double dist = 0;
	int atom1 = 0, atom2 = 0;
	bool shell_exists;
	if (reset == true) {
		cout << "Placing oxygen cores at corresponding shells." << endl;
		for (int i = 0; i < numpos; i += 3) {
			atom1 = i / 3;
			if (atom_nr[atom1] == 8 && core_shell[atom1] == "SHEL") {
				(*xyz)(i + 3 * nr_cores) = (*xyz)(i);
				(*xyz)(i + 1 + 3 * nr_cores) = (*xyz)(i + 1);
				(*xyz)(i + 2 + 3 * nr_cores) = (*xyz)(i + 2);
			}
		}
	}
	for (int i = 0; i < numpos; i += 3) {
		atom1 = i / 3;
		if (atom_nr[atom1] == 8 && core_shell[atom1] == "SHEL") {
			shell_exists = false;
			if (verbose == true) {
				cout << "Found oxygen shell at position " << atom1 + 1 << " ";
			}
			for (int j = numpos; j < numposall; j += 3) {
				atom2 = j / 3;
				//cout << "core_shell:" <<core_shell[atom2] << endl;
				if (atom_nr[atom2] == 8 && core_shell[atom2] == "CORE") {
					b[0] = (*xyz)(i) - (*xyz)(j);
					b[1] = (*xyz)(i + 1) - (*xyz)(j + 1);
					b[2] = (*xyz)(i + 2) - (*xyz)(j + 2);
					dist = sqrt(pow(b[0], 2) + pow(b[1], 2) + pow(b[2], 2));
					if (dist < scutoff) {
						shell_exists = true;
						if (verbose == true) {
							cout << ", distance to nearest core at position "
									<< atom2 + 1 << " is: " << dist << endl;
						}
						break;
					}		//end if
				}		//end if
			}		//end for j
					//could should not get here
			if (shell_exists == false) {
				cout << "\nWARNING: No core found for shell atom #" << atom1 + 1
						<< " within cutoff radius, EXIT..." << endl;
				exit(1);
			}

		}			// end if

	}			// end for i
	cout << "Shell specifications are ok.\n";
}

//append shells to structure
void ToolCalc::append_shells(bool verbose) {
	//counting oxygen atoms
	//cout << "nr_cores:" << nr_cores << endl;
	//cout << "nr_shell:" << nr_cores << endl;
	cout << "\nAppending shells to oxygen atoms." << endl;
	if (shell_modell == false && nr_cores == 0) {
		cout << "shell_modell=FALSE" << endl;
	} else {
		cout << "Can not append, shells have been defined already!" << endl;
		exit(1);
	}
	int numposall = 0;
	int numpos = atnumber * 3;
	nr_cores = 0;
	nr_shells = atnumber;
	int atom1 = 0;
	//counting oxygen
	for (int i = 0; i < numpos; i += 3) {
		atom1 = i / 3;
		if (atom_nr[atom1] == 8) {
			//cout << "core_shell:" <<core_shell[atom1]<< endl;
			nr_cores++;
		}
	}
	if (verbose == true) {
		cout << "We found " << nr_cores << " oxygen atoms, " << atnumber
				<< " atoms in total.\n";
		cout
				<< "Shells are appended to end of xyz-matrix, oxygen atoms should have been put to the end of zmatrix.";
	}
	//new number of atoms
	atnumber = atnumber + nr_cores;
	numposall = (atnumber) * 3;

	//atomtype&position
	VectorXd lxyz(3 * atnumber);
	VectorXd lgrad(3 * atnumber);
	VectorXd lmovemat(3 * atnumber);
	//Go through shells
	for (int i = 0; i < (atnumber - nr_cores) * 3; i++) {
		atom1 = i / 3;
		//cout << "xyz:" << (*xyz)(i) << endl;
		if (i % 3 == 0) {
			//			latoms[k] = atoms[k];
			//cout << endl;
			//cout << atoms[atom1];

			core_shell[atom1] = "SHEL";
			//setting charge
			if (atom_nr[atom1] == 8) {
				q[i / 3] = -2.869;
			}
		}
		(lxyz)(i) = (*xyz)(i);
		(lgrad)(i) = 0.0;
		(lmovemat)(i) = (*moveMat)(i);
		//cout << " " << (lxyz)(i);
	}

	//go through cores
	for (int i = (atnumber - nr_cores) * 3; i < (atnumber) * 3; i++) {
		atom1 = i / 3;
		if (i % 3 == 0) {
			//cout << endl;
			atoms[atom1] = "O";
			core_shell[atom1] = "CORE";
			atom_nr[atom1] = 8;
			q[atom1] = 0.869;
			//cout << atom1;
			//cout << atoms[atom1];
		}
		(lxyz)(i) = (*xyz)(i - nr_cores * 3);
		(lgrad)(i) = 0.0;
		(lmovemat)(i) = 1.0;
		;
		//cout << " " << (lxyz)(i);
	}

	(*xyz).resize(atnumber * 3);
	(*grad).resize(atnumber * 3);
	(*moveMat).resize(atnumber * 3);
	(*xyz) = lxyz.replicate(1, 1);
	(*grad) = lgrad.replicate(1, 1);
	(*moveMat) = lmovemat.replicate(1, 1);
	//charges
	//	calculation.q = new double[countatom];
	//	calculation.core_shell = new string[countatom];
	if (verbose == true) {
		for (int i = 0; i < atnumber * 3; i++) {
			atom1 = i / 3;
			if (i % 3 == 0) {
				cout << endl;
				cout << atoms[atom1] << "\t";
			}
			cout.precision(8);
			cout << (*xyz)(i) << "\t";
		}

	}
	TotalQ = 0;
	//cout << endl;
	for (int i = 0; i < atnumber * 3; i++) {
		atom1 = i / 3;
		if (i % 3 == 0) {
			TotalQ = TotalQ + q[atom1];
			if (verbose == true) {
				cout << endl;
				cout.precision(4);
				cout << fixed << atoms[atom1];
				if (core_shell[atom1] == "CORE") {
					cout << "\tc" << "\t";
				} else {
					cout << "\ts" << "\t";
				}

				cout << setw(8) << q[atom1] << "\t";
			}
		}

	}
	if (verbose == true || (TotalQ != 0)) {
		cout << "Total Charge: " << TotalQ << endl;
	}
	//interactions?
}

//strip shells from structure
void ToolCalc::strip_shells(bool verbose) {
	//nr_cores = 0;
	//nr_shells = atnumber;
	int atom1 = 0;
	cout << "\nStripping shells from oxygen atoms." << endl;
	if (shell_modell == true && nr_cores > 0) {
		cout << "Found shell_modell=TRUE" << endl;
	} else {
		cout << "Can not strip, shells have not been defined." << endl;
		exit(1);
	}
	atnumber = atnumber - nr_cores;
	VectorXd lxyz(3 * atnumber);
	//counting oxygen shells && setting charges
	for (int i = 0; i < atnumber * 3; i++) {
		atom1 = i / 3;
		if (i % 3 == 0) {

			core_shell[atom1] == "NORM";
			if (atom_nr[atom1] == 8 && core_shell[atom1] == "SHEL") {
				//cout << "core_shell:" <<core_shell[atom1]<< endl;
				q[i / 3] = -2.0;
				nr_cores--;
				nr_shells--;
			}
		}
		(lxyz)(i) = (*xyz)(i);
		//cout << "lxyz:" << (lxyz)(i) << endl;
	}
	(*xyz).resize(atnumber * 3);
	(*grad).resize(atnumber * 3);
	(*xyz) = lxyz.replicate(1, 1);
	//setting charges
	TotalQ = 0;
	//cout << endl;
	for (int i = 0; i < atnumber * 3; i++) {
		atom1 = i / 3;
		if (i % 3 == 0) {
			TotalQ = TotalQ + q[atom1];
			if (verbose == true) {

				cout.precision(4);
				cout << fixed << atoms[atom1];
				cout << setw(8) << q[atom1] << "\n";
			}
		}

	}
	if (verbose == true || (TotalQ != 0)) {
		cout << "Total Charge: " << TotalQ << endl;
	}
	cout << "Setting shell_modell=FALSE" << endl;
	shell_modell = false;
}

void ToolCalc::opt(bool pol, bool verbose) {
	VectorXd h(3 * atnumber);
	VectorXd old_xyz(3 * atnumber);
	VectorXd old_grad(3 * atnumber);
	VectorXd dxyz;
	if (verbose == true && pol == true) {
		cout << "Shell model switched ON...\n";
	}
	//function pointer to suitable energy function
	p_E = &ToolCalc::E;
	if (periodic == true) {
		p_E = &ToolCalc::E_periodic;
	}
	//Do we use shell-modell?
	if (nr_shells > 0) {
		p_E = &ToolCalc::E_Pol;
	}
	if (nr_shells > 0 && periodic == true) {
		p_E = &ToolCalc::E_periodic_Pol;
	}
	//open file for output
	ToolIO newIO;
	double lstep = step;
	double lambda = 0;
	double opt_lambda = 0;
	//using heavily nr3.h
	cout << setiosflags(ios::fixed) << setprecision(8);
	//cout << endl;
	if (verbose == true) {
		cout
				<< "###Optimization, using Polak-Ribiere Conjugate-Gradient Algorithm...###"
				<< endl;
	}
	for (int i = 0; i < maxiter; i++) {
		//writing to file
		//cout << "TEST" << endl;
		if (verbose == true) {
			newIO.moltoFile2(atnumber, energy, atoms, (*xyz), core_shell,
					nr_shells, shell_modell);
			newIO.fractoFile2(atnumber, energy, atoms, (*xyz), core_shell,
					nr_shells, shell_modell, Lx, Ly, Lz);
			setDipole(false);
			newIO.datatoFile(energy, dipnorm, false);
		}
		//only for first iteration
		if (i == 0) {
			(this->*p_E)(pol, false, true);
		}
		if (gradnorm < threshhold) {
			//			if (calctype == "PP" || calctype == "POL") {
			cout << "Optimisation CONVERGED in ";
			cout << "Cycle " << i << ", ||grad||< " << threshhold
					<< " # Energy: " << energy << " a.u. "
					<< energy * HARTREE2EV << " eV";
			//			} else {
			//				cout << "CONVERGED in ";
			//				cout <<"Cycle "<< i << ". Energy:" << energy << ", ||grad||: "
			//						<< gradnorm;
			//			}
			converged = true;
			break;
		} else {
			if (i == maxiter - 1) {
				cout << "NOT converged after " << maxiter << " iterations..."
						<< endl;
				converged = false;
				break;
			}
		}
		if (i == 0) {
			h = grad->replicate(1, 1);
			h = h * -1.0;
		}
		double old_gradn;
		old_gradn = grad->norm();
		old_grad = grad->replicate(1, 1);
		old_xyz = xyz->replicate(1, 1);
		//we are starting with line search
		//step=0.001;
		lambda = 0;
		VecDoub ls_energy(ls_steps);
		VecDoub ls_x(ls_steps);
		//VecDoub opt_step(4);
		//LINE SEARCH
		ls_x[0] = 0;
		ls_energy[0] = energy;
		//cout << "("<<0<<"), lambda:  "<<ls_x[0]<<" Energy:" << setw(12) << ls_energy[0] <<endl;
		for (int k = 1; k < ls_steps; k++) {
			//makes use of optimum line search element after some time
			if (gradnorm < 100 * threshhold && k == 1) {
				lambda = opt_lambda - (ls_steps - 1) * lstep;
				//cout << "Switching on opt_lambda..."<< lambda << endl;
			}
			lambda = (lambda + lstep * k);
			ls_x[k] = lambda;
			//cout << "ls_x:" <<ls_x[k] << " ";
			//Blas_Add_Mult(*xyz, lambda, h);
			*xyz = (*xyz) + lambda * h;
			//only rough energy, no gradients!
			//ewald_thresh = ewald_thresh * 10;
			//cutoff = cutoff - 10;
			(this->*p_E)(pol, false, false);
			//ewald_thresh = ewald_thresh * 0.1;
			//cutoff = cutoff + 10;
			ls_energy[k] = energy;
			//cout << "("<<k<<"), lambda:  "<<ls_x[k]<<" Energy:" << setw(12) << ls_energy[k] <<endl;

			//resetting coordinates
			//setXYZ(old_xyz);
			(*xyz) = old_xyz.replicate(1, 1);
			//(*xyz).copy(old_xyz);
		}
		opt_lambda = fitenergy(ls_x, ls_energy, fpoly_np);
		//cout << "opt_lambda: " << opt_lambda << "  "<< endl;
		//Final CG scale
		//Blas_Add_Mult(*xyz, opt_lambda, h);
		*xyz = (*xyz) + opt_lambda * h;

		//Filtering max elements
		dxyz = xyz->replicate(1, 1);
		//Blas_Add_Mult(dxyz, -1, old_xyz);
		dxyz = dxyz - 1.0 * old_xyz;
		//check for max displacement:
		for (int j = 0; j < (dxyz).size(); j++) {
			if (fabs(dxyz(j)) > max_displacement) {
				//cout << "Max. Displacement gt. treshhold " << (dxyz)(j) << "at iteration:"<<i <<endl;
				if (dxyz(j) > 0) {
					(dxyz)(j) = max_displacement;
				} else {
					(dxyz)(j) = -1 * max_displacement;
				}
			}
		}
		//		for (int j=0; j<(*xyz).size(); j++) {
		//						cout << (*xyz)(j) << endl;
		//				}
		//cout << "Max. Displacement:" << (dxyz)(Blas_Index_Max(dxyz)) << " at iteration: " <<i << endl;
		//Blas_Add_Mult(dxyz, 1, old_xyz);
		dxyz = dxyz + old_xyz;
		setXYZ(dxyz);
		//(*xyz).copy(dxyz);
		//Blas_Add_Mult(old_xyz, -1, *xyz);
		old_xyz = old_xyz - 1 * (*xyz);

		(this->*p_E)(pol, false, true);
		if (verbose == true) {
			cout << "Cycle " << i + 1 << " Energy:" << energy << ", ||grad||: "
					<< gradnorm << ", ||displacement||: " << old_xyz.norm()
					<< endl;
		}

		//Polak-Ribiere, p.517 numerical recipes
		double temp;
		//Blas_Add_Mult(old_grad, -1, *grad);
		old_grad = old_grad - 1.0 * (*grad);

		//Blas_Scale(-1., old_grad);
		old_grad = -1.0 * old_grad;
		temp = (old_grad.dot(*grad)) / (old_gradn * old_gradn);
		if (temp > 1) {
			temp = 1;
		}
		//cout << "temp:" << temp << endl;
		//h.scale(-1* temp );
		h = h * -1.0 * temp;

		//Blas_Add_Mult(h, 1., *grad);
		h = h + 1.0 * (*grad);

		//h.scale(-1.);
		h = h * -1.0;
		//convergence checks
		//
	}
	//newIO.printGrad(calculation);
}

//implementation in inherited classes
void ToolCalc::E(bool pol, bool verbose, bool gcalc) {
	cout << "This function should never be called..." << endl;
}
void ToolCalc::E_Pol(bool pol, bool verbose, bool gcalc) {
	cout << "This function should never be called..." << endl;
}
void ToolCalc::E_periodic(bool pol, bool verbose, bool gcalc) {
	cout << "This function should never be called..." << endl;
}
void ToolCalc::E_periodic_Pol(bool pol, bool verbose, bool gcalc) {
	cout << "This function should never be called..." << endl;
}

ToolCalc::~ToolCalc() {
}
