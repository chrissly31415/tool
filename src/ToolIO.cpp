#include <boost/regex.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>
#include <string.h>
#include <fstream>
#include <iomanip>
#include "ToolIO.h"
#include "ToolCalc.h"
#include "toolc.h"

using namespace std;

ToolIO::ToolIO() {
}

void ToolIO::parseCalctype(const char* setupfile) {
	ifstream myfile;
	string line;

	cout << "Open setup file: " << setupfile << endl;
	myfile.open(setupfile);

	boost::regex type("lj|pp|pol|zernicke", boost::regex::icase);

	if (myfile.is_open()) {
		while (!myfile.eof()) {
			getline(myfile, line);
			boost::smatch matches;
			if (boost::regex_search(line, matches, type)) {
				if (matches.str() == "PP" || matches.str() == "pp") {
					toolc::ctype = "PP";
					cout << "Type of calculation: Empirical Pair potential"
							<< endl;
				} else if (matches.str() == "LJ" || matches.str() == "lj") {
					toolc::ctype = "LJ";
					cout << "Type of calculation: Lennard-Jones potential"
							<< endl;
				} else if (matches.str() == "zernicke"
						|| matches.str() == "ZERNICKE") {
					toolc::ctype = "ZERNICKE";
					cout << "Type of calculation: ZERNICKE SHAPE EXPANSION"
							<< endl;
				}
			}
		}
		myfile.close();
	} else {
		cout << "Unable to open file " << setupfile << endl;
		;
	}
}

void ToolIO::parseSETUP(ToolCalc &calculation, const char* setupfile) {
	ifstream myfile;
	string line;
	char buffer[256];
	std::string temp;
	int maxatom = calculation.maxatoms;
	double lxyz[3 * maxatom];
	double lmoveMat[3 * maxatom];
	double lq[maxatom];
	string lcore_shell[maxatom];
	//given that nr of interaction is never larger than maxatoms...
	string specpot[maxatom][2];
	double lpot[maxatom][4];
	//initalizing all atoms are movable
	for (int i = 0; i < 3 * maxatom; i++) {
		lmoveMat[i] = 1.;
	}
	//initializing species mat
	for (int i = 0; i < maxatom; i++) {
		lq[i] = 0.;
		lcore_shell[i] = "NORM";
	}
	string latoms[maxatom];
	int latom_nr[maxatom];

	int countatom = 0;
	int countpos = 0;

	myfile.open(setupfile);

	boost::regex celldat(
			"(CELL|cell|Cell)[[:blank:]]+([0-9]{1,}\\.[0-9]*)[[:blank:]]+([0-9]{1,}\\.?[0-9]*)[[:blank:]]+([0-9]{1,}\\.?[0-9]*)");
	boost::regex coord("COORD|coord|Coord");
	boost::regex cosmofile(".{1,}\\.cosmo", boost::regex::icase);
	boost::regex coordfile("\\.cosmo|\\.xyz|\\.sdf", boost::regex::icase);
	boost::regex end("END", boost::regex::icase);
	boost::regex type("LJ|PP|POL", boost::regex::icase);
	boost::regex at("[A-z|a-z]{1,2}");
	boost::regex parallel("(par)[[:blank:]]+([0-9]{1,2})", boost::regex::icase);
	boost::regex temperat("(temp)[[:blank:]]+([0-9]{1,}\\.?[0-9]*)",
			boost::regex::icase);
	boost::regex alpha("(alpha)[[:blank:]]+([0-9]{1,}\\.?[0-9]*)",
			boost::regex::icase);
	boost::regex ran("(ran)[[:blank:]]+([0-9]{1,}\\.?[0-9]*)",
			boost::regex::icase);
	boost::regex mcsteps("(mcsteps)[[:blank:]]+([0-9]{1,9})",
			boost::regex::icase);
	boost::regex resets("(resets)[[:blank:]]+([0-9]{1,9})",
			boost::regex::icase);

	//boost::regex cs("[cC|sS]{1}");
	boost::regex frozen("([F|T]{1})[[:blank:]]+([F|T])[[:blank:]]+([F|T])",
			boost::regex::icase);
	boost::regex atpos(
			"(-?[0-9]{1,}\\.?[0-9]*)[[:blank:]]+(-?[0-9]{1,}\\.?[0-9]*)[[:blank:]]+(-?[0-9]{1,}\\.?[0-9]*)");
	boost::regex species("SPECIE", boost::regex::icase);
	boost::regex specdat(
			"([A-z|a-z]{1,2})[[:blank:]]+([cC|sS|nN])?[[:blank:]]+(-?[0-9]{0,}\\.?[0-9]*)");
	boost::regex potentials("POTEN", boost::regex::icase);
	boost::regex potdat(
			"([A-z|a-z]{1,2})[[:blank:]]+([A-z|a-z]{1,2})[[:blank:]]+(-?[0-9]{1,}\\.?[0-9]*)[[:blank:]]+(-?[0-9]{1,}\\.?[0-9]*)[[:blank:]]+(-?[0-9]{1,}\\.?[0-9]*)[[:blank:]]+(-?[0-9]{1,}\\.?[0-9]*)");
	cout << "\ndir: " << getcwd(buffer, 256) << "\n";
	bool xyzstart = false;
	bool specstart = false;
	bool potstart = false;
	if (myfile.is_open()) {
		while (!myfile.eof()) {
			getline(myfile, line);
			if (boost::starts_with(line, "#"))
				continue;
			boost::smatch matches;
			//parsing periodic cell info
			if (boost::regex_search(line, matches, celldat)) {
				//getline(myfile, line);
				if (boost::regex_search(line, matches, celldat)) {
					cout << "Reading periodic cell information...\n";
					calculation.Lx = string2double(matches[2].str());
					calculation.Ly = string2double(matches[3].str());
					calculation.Lz = string2double(matches[4].str());
				}
				getline(myfile, line);
				if (boost::regex_search(line, matches, alpha)) {
					calculation.alpha = string2double(matches[2].str());
				}
			}
			//parsing parallel information
			if (boost::regex_search(line, matches, parallel)) {
				calculation.nproc = string2integer(matches[2].str());
				cout << "Using " << calculation.nproc << " processor(s).\n";
			}
			if (boost::regex_search(line, matches, ran)) {
				calculation.rseed = string2integer(matches[2].str());
				cout << "Using " << calculation.rseed << " as random seed.\n";
			}
			if (boost::regex_search(line, matches, temperat)) {
				calculation.temperature = string2double(matches[2].str());
				cout << "Temperature " << calculation.temperature << " K\n";
			}

			if (boost::regex_search(line, matches, mcsteps)) {
				calculation.global_maxiter = string2integer(matches[2].str());
				cout << "MC Iterations: " << calculation.global_maxiter << "\n";
			}
			if (boost::regex_search(line, matches, resets)) {
				calculation.nr_reset = string2integer(matches[2].str());
				cout << "Nr. of MC resets: " << calculation.nr_reset << "\n";
			}
			//parsing species information
			if (boost::regex_search(line, matches, species)) {

				specstart = true;
				cout
						<< "Reading species information, please specify core atoms at the end...\n";
				getline(myfile, line);
				int l = 0;
				while (boost::regex_search(line, matches, specdat)
						&& (specstart == true)) {
					if (matches[2].str() == "c" || matches[2].str() == "C"
							|| matches[2].str() == "s"
							|| matches[2].str() == "S"
							|| matches[2].str() == "n"
							|| matches[2].str() == "N") {
						if (matches[2].str() == "c"
								|| matches[2].str() == "C") {
							lcore_shell[l] = "CORE";
							calculation.nr_cores++;
						} else {
							if (matches[2].str() == "s"
									|| matches[2].str() == "S") {
								lcore_shell[l] = "SHEL";
								calculation.nr_shells++;
							} else {
								lcore_shell[l] = "NORM";
							}
						}
						lq[l] = string2double(matches[3].str());
						l++;
					} else {
						for (int i = 0; i < countatom; i++) {
							if (latoms[i] == matches[1].str()) {
								lq[i] = string2double(matches[3].str());
							}

						}

					}
					calculation.nr_species++;
					getline(myfile, line);
				}
				//continue;
				if (boost::regex_search(line, matches, end)) {
					//cout << "Found " << matches <<endl;
					specstart = false;
				}
			}
			if (boost::regex_search(line, matches, potentials)) {
				potstart = true;
				cout << "Reading potential information...\n";
				getline(myfile, line);
				while (boost::regex_search(line, matches, potdat)
						&& (potstart == true)) {
					specpot[calculation.interactions][0] = matches[1].str();
					specpot[calculation.interactions][1] = matches[2].str();
					lpot[calculation.interactions][0] = string2double(
							matches[3].str());
					lpot[calculation.interactions][1] = string2double(
							matches[4].str());
					lpot[calculation.interactions][2] = string2double(
							matches[5].str());
					lpot[calculation.interactions][3] = string2double(
							matches[6].str());
					calculation.interactions++;
					getline(myfile, line);
					if (boost::regex_search(line, matches, end)) {
						//cout << "Found " << matches <<endl;
						potstart = false;
					}
				}

			}

			if (boost::regex_search(line, matches, coord)) {
				xyzstart = true;
				cout << "Found xyz-section, reading molecular structure...\n";
				continue;
			}
			if (boost::regex_search(line, matches, end)) {
				xyzstart = false;
			}
			if (xyzstart) {
				if (boost::regex_search(line, matches, at)) {
					//cout << "Found atom " << matches << endl;
					latoms[countatom] = matches.str();
					for (int i = 0; i < 87; i++) {
						if (matches.str() == ToolCalc::elements[i]) {
							latom_nr[countatom] = i;
						}
					}
					countatom++;
				}
				if (boost::regex_search(line, matches, atpos)) {
					//converting from std::string to char*
					//matches[0] contains all results
					lxyz[countpos] = string2double(matches[1].str());
					countpos++;
					lxyz[countpos] = string2double(matches[2].str());
					countpos++;
					lxyz[countpos] = string2double(matches[3].str());
					countpos++;
				}
				if (boost::regex_search(line, matches, frozen)) {
					//resetting counter
					countpos = countpos - 3;
					temp = matches[1].str();
					if (temp == "F" || temp == "f") {
						lmoveMat[countpos] = 0.0;
					} else {
						lmoveMat[countpos] = 1.;
					}
					countpos++;
					temp = matches[2].str();
					if (temp == "F" || temp == "f") {
						lmoveMat[countpos] = 0.;
					} else {
						lmoveMat[countpos] = 1.;
					}
					countpos++;
					temp = matches[3].str();
					if (temp == "F" || temp == "f") {
						lmoveMat[countpos] = 0.;
					} else {
						lmoveMat[countpos] = 1.;
					}
					countpos++;
				}
				if (boost::regex_search(line, matches, cosmofile)) {
					temp = matches[0].str();
					cout << "Reading coordinates from file: " << matches[0]
							<< endl;
					parseCOSMO(calculation, temp);
				}
				//
			}
		}
		myfile.close();
	} else
		cout << "Unable to open file ";

	//in case coordinates have not been defined in file
	if (calculation.atnumber == 0) {
		//definition of calculation object
		calculation.atnumber = countatom;
		calculation.xyz = new ToolCalc::VectorXd(3 * calculation.atnumber);
		for (int i = 0; i < 3 * calculation.atnumber; ++i) {
			(*calculation.xyz)(i) = lxyz[i];
		}

		//2* countatom:to allow for appending shells
		calculation.atoms = new string[2 * calculation.atnumber];
		calculation.atom_nr = new int[2 * calculation.atnumber];
		for (int i = 0; i < calculation.atnumber; i++) {
			calculation.atoms[i] = latoms[i];
			calculation.atom_nr[i] = latom_nr[i];
		}
	}

	calculation.grad = new ToolCalc::VectorXd(3 * calculation.atnumber);
	calculation.moveMat = new ToolCalc::VectorXd(3 * calculation.atnumber);
	for (int i = 0; i < calculation.atnumber * 3; i++) {
		(*calculation.moveMat)(i) = lmoveMat[i];
		(*calculation.grad)(i) = 0.0;
	}

	//setting potentials
	if (toolc::ctype == "PP") {
		cout << "#NR of species: " << calculation.nr_species;
		cout << " -- Cores: " << calculation.nr_cores << " -- Shells: "
				<< calculation.nr_shells << endl;
		if (calculation.nr_shells > 0
				&& (calculation.nr_shells + calculation.nr_cores
						!= calculation.atnumber)) {
			cout << "###WARNING!Shells or cores definition is incomplete...###"
					<< endl;
		}

		if (calculation.atnumber != calculation.nr_species) {
			cout
					<< "###Warning! Number of species does not coincide with number of atoms specified!###"
					<< "atoms  :" << calculation.atnumber << endl << "species:"
					<< calculation.nr_species << endl << endl;
		}

		calculation.q = new double[2 * countatom];
		calculation.core_shell = new string[2 * countatom];
		for (int i = 0; i < countatom; i++) {
			calculation.q[i] = lq[i];
			calculation.TotalQ = calculation.TotalQ + calculation.q[i];
			calculation.core_shell[i] = lcore_shell[i];
			cout.precision(4);
			cout << fixed << calculation.atoms[i] << "\t"
					<< calculation.core_shell[i] << " species\tcharge ";
			cout << setw(8) << calculation.q[i] << endl;
		}
		cout << "Total charge: " << calculation.TotalQ << " e" << endl;
	}

	if (toolc::ctype == "PP" || toolc::ctype == "LJ") {
		cout << "#NR of interactions: " << calculation.interactions << endl;
		calculation.paraMat = new double*[calculation.interactions];
		calculation.interMat = new int*[calculation.interactions];
		for (int i = 0; i < calculation.interactions; i++) {
			calculation.interMat[i] = new int[2];
			calculation.paraMat[i] = new double[4];
			calculation.paraMat[i][0] = lpot[i][0];
			calculation.paraMat[i][1] = lpot[i][1];
			calculation.paraMat[i][2] = lpot[i][2];
			calculation.paraMat[i][3] = lpot[i][3];
			for (int j = 0; j < 87; j++) {
				if (specpot[i][0] == ToolCalc::elements[j]) {
					calculation.interMat[i][0] = j;
				}
			}
			for (int j = 0; j < 87; j++) {
				if (specpot[i][1] == ToolCalc::elements[j]) {
					calculation.interMat[i][1] = j;
				}
			}
			cout.precision(4);
			cout << fixed << calculation.interMat[i][0] << "\t";
			cout << fixed << calculation.interMat[i][1] << "\t";
			cout << setw(10) << calculation.paraMat[i][0] << " ";
			cout << setw(10) << calculation.paraMat[i][1] << " ";
			cout << setw(10) << calculation.paraMat[i][2] << " ";
			cout << setw(10) << calculation.paraMat[i][3] << endl;
		}
	}
}

void ToolIO::parseCOSMO(ToolCalc &calculation, string filename) {
	int maxatom = calculation.maxatoms;
	double lxyz[3 * maxatom];
	double segment_pos[3 * maxatom * 10];
	int seg_atom[maxatom];
	double sq[maxatom * 10]; //charge
	double sa[maxatom * 10]; //area
	//double pot[maxatom * 10]; //potential
	string latoms[maxatom];
	int latom_nr[maxatom];

	string atreg = "(-?[0-9]{1,}\\.[0-9]*)[[:blank:]]*";
	string atname = "([A-z|a-z]{1,2})";
	int pos = 0;
	int nseg = 0;
	boost::regex atpos_header("\\$coord_rad", boost::regex::icase);
	boost::regex atpos(atreg + atreg + atreg + atname);
	boost::regex segment_header("\\$segment_information", boost::regex::icase);
	boost::regex segments_re(
			"([0-9]{1,})[[:blank:]]*" + atreg + atreg + atreg + atreg + atreg);
	bool atoms_section = false;
	bool segment_section = false;
	cout << "open:" << filename << endl;
	ifstream myfile;
	string line;
	myfile.open(filename.c_str());
	if (myfile.is_open()) {
		while (!myfile.eof()) {
			getline(myfile, line);
			if (boost::starts_with(line, "#"))
				continue;
			boost::smatch matches;
			if (boost::regex_search(line, matches, atpos_header)) {
				//cout << "line:" << line << endl;
				atoms_section = true;
				continue;
			} else if (boost::regex_search(line, matches, segment_header)) {
				segment_section = true;
				continue;
			}

			if (atoms_section && boost::regex_search(line, matches, atpos)) {
				for (int i = 1; i < 4; ++i) {
					lxyz[pos] = string2double(matches[i].str()) / ANG2BOHR;
					pos++;
				}
				line = matches[4].str();

				latoms[(pos - 1) / 3] = line;
				for (int i = 0; i < 87; i++) {
					if (boost::iequals(line, ToolCalc::elements[i])) {
						latom_nr[(pos - 1) / 3] = i;
					}
				}
				continue;

			} else if (atoms_section) {
				atoms_section = false;
			}

			if (segment_section
					&& boost::regex_search(line, matches, segments_re)) {
				seg_atom[nseg / 3] = string2integer(matches[1].str());
				for (int i = 2; i < 5; ++i) {
					segment_pos[nseg] = string2double(
							matches[i].str())/ANG2BOHR;
					nseg++;
				}
				sq[nseg / 3] = string2double(matches[5].str());
				sa[nseg / 3] = string2double(matches[6].str());

			} else if (segment_section) {
				segment_section = false;
			}

		}
	} else
		cout << "Unable to open file: " << filename << endl;
	myfile.close();

	calculation.atnumber = pos / 3;

	calculation.xyz = new ToolCalc::VectorXd(pos);
	for (int i = 0; i < pos; ++i) {
		(*calculation.xyz)(i) = lxyz[i];
	}

	calculation.atoms = new string[calculation.atnumber];
	calculation.atom_nr = new int[calculation.atnumber];
	for (int i = 0; i < calculation.atnumber; i++) {
		calculation.atoms[i] = latoms[i];
		calculation.atom_nr[i] = latom_nr[i];
	}

	//segments
	nseg = nseg / 3;
	cout << "Number of segments:" << (nseg) << endl;
	calculation.nseg = nseg;
	calculation.segments = ToolCalc::VectorXd::Zero(nseg * 3);
	calculation.scharge = ToolCalc::VectorXd::Zero(nseg);
	calculation.sarea = ToolCalc::VectorXd::Zero(nseg);
	calculation.satom = ToolCalc::VectorXi::Zero(nseg);
	for (int i = 0; i < nseg * 3; i += 3) {
		calculation.segments(i) = segment_pos[i];
		calculation.segments(i + 1) = segment_pos[i + 1];
		calculation.segments(i + 2) = segment_pos[i + 2];

		calculation.satom(i / 3) = seg_atom[i / 3] - 1;
		calculation.scharge(i / 3) = sq[i / 3];
		calculation.sarea(i / 3) = sa[i / 3];

	}

	cout << "Parsing .cosmo file - no atoms:" << calculation.atnumber << endl;

}

int ToolIO::string2integer(string str) {
	char* temp = new char[str.length() + 1];
	strcpy(temp, str.c_str());
	return atoi(temp);
}

double ToolIO::string2double(string str) {
	char* temp = new char[str.length() + 1];
	strcpy(temp, str.c_str());
	return atof(temp);
}

//prints fractional coordinates in Gulp readable format
void ToolIO::fractoFile(ToolCalc &calculation, bool shel, bool periodic) {
	ofstream f;
	f.open("tool.gin", ios::app);
	f << endl << "single conp" << endl;
	f << endl << "name tool" << endl;
	f << endl << "cell" << endl;
	f << calculation.Lx << " " << calculation.Ly << " " << calculation.Lz << " "
			<< "90.000000 90.000000 90.000000" << endl;
	f << "fractional" << endl;
	int k = 0;
	for (int i = 0; i < calculation.atnumber * 3; i++) {
		if (calculation.core_shell[i / 3] == "SHEL" && shel == false) {
			continue;
		}
		if (i % 3 == 0) {
			f << fixed << calculation.atoms[k] << "\t";
			k++;
		}
		//f << "Positions: " <<endl << *calculation.xyz;

		if (i % 3 == 0) {
			f << setw(14) << (*calculation.xyz)(i) / calculation.Lx << "\t";
		}
		if (i % 3 == 1) {
			f << setw(14) << (*calculation.xyz)(i) / calculation.Ly << "\t";
		}
		if (i % 3 == 2) {
			f << setw(12) << (*calculation.xyz)(i) / calculation.Lz << "\t";
		}
		if ((i + 1) % 3 == 0) {
			f << endl;
		}
	}

	f << endl << "space" << endl << "P1" << endl;
	f << endl << "print 1" << endl;
	f.close();
}

//prints fractional coordinates in Gulp readable format
void ToolIO::printFrac(ToolCalc &calculation, bool shel, bool periodic) {
	//cout << endl <<  "single conp" << endl;
	//cout << endl << "name tool" <<  endl;
	cout << endl << "cell" << endl;
	cout << calculation.Lx << " " << calculation.Ly << " " << calculation.Lz
			<< " " << "90.000000 90.000000 90.000000" << endl;
	cout << "fractional" << endl;
	int k = 0;
	for (int i = 0; i < calculation.atnumber * 3; i++) {
		if (calculation.core_shell[i / 3] == "SHEL" && shel == false) {
			continue;
		}
		if (i % 3 == 0) {
			cout << fixed << calculation.atoms[k] << "\t";
			k++;
		}
		//cout << "Positions: " <<endl << *calculation.xyz;

		if (i % 3 == 0) {
			cout << setw(14) << (*calculation.xyz)(i) / calculation.Lx << "\t";
		}
		if (i % 3 == 1) {
			cout << setw(14) << (*calculation.xyz)(i) / calculation.Ly << "\t";
		}
		if (i % 3 == 2) {
			cout << setw(12) << (*calculation.xyz)(i) / calculation.Lz << "\t";
		}
		if ((i + 1) % 3 == 0) {
			cout << endl;
		}
	}

	//cout << endl << "space"<<endl << "P1"<<endl;
	//cout << endl << "print 1" << endl;

}

//Function prints xyz data
void ToolIO::printMol(ToolCalc &calculation, bool shel) {
	//Ausgabe
	cout.precision(8);
	cout << endl << " " << calculation.atnumber << endl;
	if (calculation.calculated == true) {
		if (toolc::ctype == "PP") {
			cout << setw(12) << calculation.energy * 27.211396132 << endl;
		} else {
			cout << setw(12) << calculation.energy << endl;
		}

	}
	//Ausgabe
	cout.precision(8);
	int k = 0;
	for (int i = 0; i < calculation.atnumber * 3; i++) {
		//cout <<i<< endl;
		//cout << calculation.core_shell[i];
		if (calculation.core_shell[i / 3] == "CORE" && shel == false) {

			continue;
		}
		if (i % 3 == 0) {
			//somehow i is changed when given here????

			cout << fixed << calculation.atoms[k] << "\t";
			k++;
		}

		//cout << "Positions: " <<endl << *calculation.xyz;
		cout << setw(12) << (*calculation.xyz)(i) << "\t";
		if ((i + 1) % 3 == 0) {
			cout << endl;
		}
	}
}

//Function prints gradient
void ToolIO::printGrad(ToolCalc &calculation) {
	calculation.periodic = true;
	if (calculation.calculated == true) {
		cout << "Gradient:\tx\t\ty\t\tz" << endl;
		cout.precision(8);
		int k = 0;
		for (int i = 0; i < calculation.atnumber * 3; i++) {
			if (i % 3 == 0) {
				//somehow i is changed when given here????
				cout << fixed << calculation.atoms[k] << "\t";
				k++;
			}
			//cout << "Positions: " <<endl << *calculation.xyz;
			if (calculation.periodic == true) {
				if (i % 3 == 0) {
					cout << setw(14)
							<< (*calculation.grad)(i) * HARTREE2EV * ANG2BOHR
									* calculation.Lx << "\t";
				}
				if (i % 3 == 1) {
					cout << setw(14)
							<< (*calculation.grad)(i) * HARTREE2EV * ANG2BOHR
									* calculation.Ly << "\t";
				}
				if (i % 3 == 2) {
					cout << setw(12)
							<< (*calculation.grad)(i) * HARTREE2EV * ANG2BOHR
									* calculation.Lz << "\t";
				}
			} else {
				cout << setw(14)
						<< (*calculation.grad)(i) * HARTREE2EV * ANG2BOHR
						<< "\t";
			}
			if ((i + 1) % 3 == 0) {
				cout << endl;
			}
		}

	}
	cout << "||grad||:" << calculation.gradnorm << endl;
	cout << "Maximum abs. gradient component at position "
			<< calculation.gradmaxelem << ": " << calculation.gradinfnorm
			<< endl;
}

//static needs only to be declared in header
void ToolIO::eigencoord2file(Eigen::MatrixXd coords, const char* filename) {
	//Eigen::IOFormat NumpyFmt(Eigen::FullPrecision, 0, ", ", ",\n", "[", "]", "[", "]");
	const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
			Eigen::DontAlignCols, ", ", "\n");
	std::ofstream file(filename);
	if (file.is_open()) {
		file << coords.format(CSVFormat) << endl;
	}
	file.close();
}

//print to file
void ToolIO::moltoFile(ToolCalc* calculation) {
	ofstream f;
	//cout << "being called..."<< endl;
	f.open("movie.xyz", ios::app);
	f << " " << calculation->atnumber << endl;
	if (calculation->calculated == true) {
		if (toolc::ctype == "PP") {
			f.precision(12);
			f << setw(12) << calculation->energy * 27.211 << endl;
		} else {
			f.precision(12);
			f << setw(12) << calculation->energy << endl;
		}

	}
	//Ausgabe
	f.precision(8);
	int k = 0;
	for (int i = 0; i < calculation->atnumber * 3; i++) {
		if (i % 3 == 0) {
			//somehow i is changed when given here????
			f << fixed << calculation->atoms[k] << "\t";
			k++;
		}
		//cout << "Positions: " <<endl << *calculation.xyz;
		f << setw(12) << (*calculation->xyz)(i) << "\t";
		if ((i + 1) % 3 == 0) {
			f << endl;
		}
	}
	f.close();
}
//second moltofile for call during optimisation without call to object
void ToolIO::moltoFile2(int &atnumber, double &energy, string atoms[],
		ToolCalc::VectorXd &xyz, string core_shell[], int nr_shells, bool shel,
		string type) {
	//cout << "atnumber:" <<atnumber<<endl;
	ofstream f;
	f.open("movie.xyz", ios::app);
	f << " " << atnumber << endl;

	if (type == "LJ") {
		f.precision(12);
		f << setw(12) << energy << endl;
	} else {
		f.precision(12);
		f << setw(12) << energy * 27.211 << endl;
	}
	//Ausgabe
	f.precision(8);
	int k = 0;
	for (int i = 0; i < atnumber * 3; i++) {
		if (core_shell[i / 3] == "CORE" && shel == false) {
			continue;
		}
		if (i % 3 == 0) {
			//somehow i is changed when given here????
			f << fixed << atoms[k] << "\t";
			k++;
		}
		//cout << "Positions: " <<endl << *calculation.xyz;
		f << setw(12) << (xyz)(i) << "\t";
		if ((i + 1) % 3 == 0) {
			f << endl;
		}
	}
	f.close();
}
//second version for call from opt function
void ToolIO::fractoFile2(int &atnumber, double &energy, string atoms[],
		ToolCalc::VectorXd &xyz, string core_shell[], int nr_shells, bool shel,
		double Lx, double Ly, double Lz) {
	//cout << "I am being called!!!!!!!!!"<<endl;
	ofstream f;
	f.open("tool.gin", ios::app);
	f << endl << "single conp" << endl;
	f << endl << "name tool" << endl;
	f << endl << "cell" << endl;
	f << Lx << " " << Ly << " " << Lz << " " << "90.000000 90.000000 90.000000"
			<< endl;
	f << "fractional" << endl;
	int atom1 = 0;
	//cout << "atnumber:"<<atnumber<<endl;
	for (int i = 0; i < atnumber * 3; i++) {
		atom1 = i / 3;
		if (core_shell[i / 3] == "CORE" && shel == false) {
			continue;
		}
		if (i % 3 == 0) {
			f << fixed << atoms[i / 3] << "\t";
		}
		//f << "Positions: " <<endl << *calculation.xyz;

		if (i % 3 == 0) {
			f << setw(14) << xyz(i) / Lx << "\t";
		}
		if (i % 3 == 1) {
			f << setw(14) << xyz(i) / Ly << "\t";
		}
		if (i % 3 == 2) {
			f << setw(12) << xyz(i) / Lz << "\t";
		}
		if ((i + 1) % 3 == 0) {
			f << endl;
		}
	}

	f << endl << "space" << endl << "P1" << endl;
	f << endl << "print 1" << endl;
	f.close();
}

void ToolIO::printRDF(int &atnumber, double &energy, string atoms[],
		ToolCalc::VectorXd &xyz, double cutoff, double Lx, double Ly,
		double Lz) {
	int nmax = 0;
	double L[3], a[3], b[3];
	double dist = 0, rsq = 0;
	int atom1 = 0, atom2 = 0;
	for (int nact = 0; nact <= nmax; nact++) {
		for (int nx = -nact; nx <= nact; nx++) {
			for (int ny = -nact; ny <= nact; ny++) {
				for (int nz = -nact; nz <= nact; nz++) {
					if (abs(nx) == nact || abs(ny) == nact || abs(nz) == nact) {
						L[0] = nx * Lx;
						L[1] = ny * Ly;
						L[2] = nz * Lz;
						for (int i = 0; i < atnumber * 3; i += 3) {
							atom1 = i / 3;
							a[0] = (xyz)(i) + L[0];
							a[1] = (xyz)(i + 1) + L[1];
							a[2] = (xyz)(i + 2) + L[2];
							for (int j = 0; j < atnumber * 3; j += 3) {
								//skipo inner loop
								if (nx == 0 && ny == 0 && nz == 0) {
									//cout << "We are in unit cell..." << endl;
									if (i == j) {
										//skipping selfinteraction...
										continue;
									}
								}
								atom2 = j / 3;
								b[0] = a[0] - (xyz)(j);
								b[1] = a[1] - (xyz)(j + 1);
								b[2] = a[2] - (xyz)(j + 2);
								rsq = pow(b[0], 2) + pow(b[1], 2)
										+ pow(b[2], 2);
								dist = sqrt(rsq);
								cout << "Distance between atom ["
										<< atoms[atom1] << atom1 << "] and ["
										<< atoms[atom2] << atom2 << "]: "
										<< dist << endl;
								//distance cut off, avoiding also Coulomb explosion??
								if (dist > cutoff) {
									//if (dist > cutoff) {
									//cout << "Cut off too big...";
									continue;
								}
							}
						}

					}
				}
			}
		}
	}
}
void ToolIO::datatoFile(double energy, double dipnorm, bool erase) {
	ofstream f;
	f.open("tool.data", ios::app);
	if (erase == true) {
		f.close();
		f.open("tool.data", ios::out);
	}
	f.precision(12);
	f << setw(14) << energy << "\t" << setw(4) << dipnorm << endl;
	f.close();
}
//prints grid to file, csv format
void ToolIO::genom2File(double lenergy, int lgenom_compressed[], int atnumber,
		bool erase) {
	ofstream f;
	f.open("grid.csv", ios::app);
	if (erase == true) {
		f.close();
		f.open("tool.grid", ios::out);
	}
	f.precision(12);
	f << setw(14) << lenergy << ",";

	for (int i = 0; i < (atnumber); i++) {
		f << setw(4) << lgenom_compressed[i] << ",";
	}
	f << endl;
	f.close();
}

// gives out grid
void ToolIO::gridout(int lgenom[], int lgridp) {
	for (int i = 0; i < lgridp; i++) {
		if (i % 50 == 0) {
			cout << endl;
		}
		cout << setw(2) << lgenom[i];
		//if (lgenom[i] != 0) {
		//cout <<" "<<i;
		//}
	}
	cout << endl;
}

//prints out genom with energy
void ToolIO::genomout(double lenergy, int lgenom_compressed[], int atnumber) {
	cout << lenergy << " ";
	for (int i = 0; i < (atnumber); i++) {
		cout << setw(2) << lgenom_compressed[i] << " ";
	}
	cout << endl;
}

//give out parameters
void ToolIO::printParameters(ToolCalc &calc) {
	cout << endl;
	cout.precision(3);

	if (toolc::ctype == "PP") {
		if (calc.periodic == true) {
			cout << "Ewald alpha: " << calc.alpha << endl;
			cout << "Threshhold Ewald sums: " << calc.ewald_thresh << endl;
		}
		if (calc.dipole_cor == true) {
			cout << "Dipole correction is ON " << endl;
		} else {
			cout << "Dipole correction is OFF " << endl;
		}
	}

	if (toolc::ctype == "PP" || toolc::ctype == "LJ") {
		cout << "Interaction cutoff: " << calc.cutoff << endl;
		cout << "Line search parameter (step): " << calc.step << endl;
		cout << "Maximum displacement in geometry optimization: "
				<< calc.max_displacement << endl;
		cout << "||Grad|| threshhold geometry optimization: " << calc.threshhold
				<< endl;
		cout << "Maximum iteration geometry optimization: " << calc.maxiter
				<< endl;
		cout << "Maximum iteration MC sampling: " << calc.global_maxiter
				<< endl;
		cout << "MC temperature: " << calc.temperature << endl;
	}

	cout << "Number of processors: " << calc.nproc << endl;
}

//Function prints timing
void ToolIO::printTiming(timeval &start, timeval &end) {
	int difsec;
	int difusec;
	difsec = end.tv_sec - start.tv_sec;
	difusec = end.tv_usec - start.tv_usec;
	//cout << start.tv_sec << ':' << start.tv_usec << endl;
	//cout << end.tv_sec << ':' << end.tv_usec << endl;
	//cout <<"Simulation ended after: " << setprecision(6) << dif << " sec.";
	if (difusec < 0) {
		difsec--;
		difusec = 1000000 + difusec;
		//cout <<"modifying";
	}
	//cout.precision(8);
	cout << difsec << "." << difusec << " sec" << endl;
}

//Destruktor
ToolIO::~ToolIO() {
}
