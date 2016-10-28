#include "ToolCalcZE.h"
#include "ToolIO.h"
#include <iostream>
#include <Eigen/Core>

using namespace std;

ToolCalcZE::ToolCalcZE() {
	Lx = 0.0;
	Ly = 0.0;
	Lz = 0.0;

}

void ToolCalcZE::printSegments() {
	cout << "atom: " << "seg:" << "x:       " << "y:       " << "z:       "
			<< endl;
	cout << fixed << setprecision(4);
	for (int i = 0; i < nseg * 3; i += 3) {
		cout << setw(5) << satom(i / 3) << setw(4) << (i / 3);
		cout << setw(10) << segments(i) << setw(10) << segments(i + 1)
				<< setw(10) << segments(i + 2) << endl;
	}

}

void ToolCalcZE::defineCubeSize() {
	//reshape matrix - just a view
	Eigen::MatrixXd M1(*xyz);

	Eigen::Map<
			Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
					Eigen::RowMajor>> M2(M1.data(), 3, M1.cols() / 3);
	Eigen::Vector3d maxv = M2.colwise().maxCoeff();
	Eigen::Vector3d minv = M2.colwise().minCoeff();
	double delta = 1.0;
	if (Lx < 1E-5) {
		Lx = maxv(0) - minv(0) + delta;
	}
	if (Ly < 1E-5) {
		Ly = maxv(1) - minv(1) + delta;
	}
	if (Lz < 1E-5) {
		Lz = maxv(2) - minv(2) + delta;
	}
	cout << "Lx,Ly,Lz: ";
	cout << setprecision(4);
	cout << setw(10);
	cout << Lx << " " << Ly << " " << Lz << endl;
}

void ToolCalcZE::seg2grid(int na, int nb, int nc, bool verbose) {
	defineCubeSize();
	gx = na;
	gy = nb;
	gz = nc;
	gridp = na * nb * nc;

	double cubex = 0, cubey = 0, cubez = 0;
	double spacex = 0, spacey = 0, spacez = 0;
	spacex = Lx / (double) na;
	spacey = Ly / (double) nb;
	spacez = Lz / (double) nc;
	int hit = 0;
	int segment = 0;
	int segment_old = 0;
	int atom1;
	//cout << atnumber << endl;
	genom = new int[gridp];

	//segment=0;
	for (int a = 0; a < na; a++) {

		cubex = a * spacex;
		for (int b = 0; b < nb; b++) {

			cubey = b * spacey;
			for (int c = 0; c < nc; c++) {
				cubez = c * spacez;
				genom[segment] = 0;
				for (int i = 0; i < atnumber * 3; i += 3) {
					atom1 = i / 3;
					if (((*xyz)(i) > cubex) && ((*xyz)(i) < (cubex + spacex))
							&& (((*xyz)(i + 1) > cubey)
									&& ((*xyz)(i + 1) < (cubey + spacey)))
							&& (((*xyz)(i + 2) > cubez)
									&& ((*xyz)(i + 2) < (cubez + spacez)))) {
						//cout << "Hit:"<< (*xyz)(i) << " " <<(*xyz)(i+1) << " "
						//		<<(*xyz)(i+2)<< endl;
						//cout << "Actual seg: "<< segment<< " coords x:" <<cubex
						//	<<" y:"<<cubey <<" z:" <<cubez<<endl;
						if (segment == segment_old) {
							cout << "\nWARNING: Two atoms at one grid_point:"
									<< segment << endl;
							segment++;
						}
						genom[segment] = (int) q[atom1];
						segment_old = segment;
						//cout << "Genom["<<segment<<"]"<<genom[segment] << endl;
						//cout << " "<< segment;
						hit++;
						//cout << hit<< endl;
					}

					//(*xyz)(i)=(*xyz)(i)*(1-(*moveMat)(i))+j*(*moveMat)(i);
				}
				segment++;

			}
		}

	}
	//cout << "Actual seg: "<< segment<< " coords x:" <<cubex+spacex<<" y:"
	//		<<cubey+spacey <<" z:" <<cubez+spacez<<endl;
	if (verbose == true) {
		cout << endl << "xyz2grid:" << endl << "Grid size x:" << na << " y:"
				<< nb << " z:" << nc << endl;
		cout << "New grid created with " << gridp << " positions.\n";
		cout << "Projected " << hit << " out of " << atnumber
				<< " total atoms into Grid..." << endl;
	}
}

//projects coordinates to grid
void ToolCalcZE::xyz2grid(int na, int nb, int nc, bool verbose) {
	gx = na;
	gy = nb;
	gz = nc;
	gridp = na * nb * nc;

	double cubex = 0, cubey = 0, cubez = 0;
	double spacex = 0, spacey = 0, spacez = 0;
	spacex = Lx / (double) na;
	spacey = Ly / (double) nb;
	spacez = Lz / (double) nc;
	int hit = 0;
	int segment = 0;
	int segment_old = 0;
	int atom1;
	//cout << atnumber << endl;
	genom = new int[gridp];

	//segment=0;
	for (int a = 0; a < na; a++) {

		cubex = a * spacex;
		for (int b = 0; b < nb; b++) {

			cubey = b * spacey;
			for (int c = 0; c < nc; c++) {
				cubez = c * spacez;
				genom[segment] = 0;
				for (int i = 0; i < atnumber * 3; i += 3) {
					atom1 = i / 3;
					if (((*xyz)(i) > cubex) && ((*xyz)(i) < (cubex + spacex))
							&& (((*xyz)(i + 1) > cubey)
									&& ((*xyz)(i + 1) < (cubey + spacey)))
							&& (((*xyz)(i + 2) > cubez)
									&& ((*xyz)(i + 2) < (cubez + spacez)))) {
						//cout << "Hit:"<< (*xyz)(i) << " " <<(*xyz)(i+1) << " "
						//		<<(*xyz)(i+2)<< endl;
						//cout << "Actual seg: "<< segment<< " coords x:" <<cubex
						//	<<" y:"<<cubey <<" z:" <<cubez<<endl;
						if (segment == segment_old) {
							cout << "\nWARNING: Two atoms at one grid_point:"
									<< segment << endl;
							segment++;
						}
						genom[segment] = (int) q[atom1];
						segment_old = segment;
						//cout << "Genom["<<segment<<"]"<<genom[segment] << endl;
						//cout << " "<< segment;
						hit++;
						//cout << hit<< endl;
					}

					//(*xyz)(i)=(*xyz)(i)*(1-(*moveMat)(i))+j*(*moveMat)(i);
				}
				segment++;

			}
		}

	}
	//cout << "Actual seg: "<< segment<< " coords x:" <<cubex+spacex<<" y:"
	//		<<cubey+spacey <<" z:" <<cubez+spacez<<endl;
	if (verbose == true) {
		cout << endl << "xyz2grid:" << endl << "Grid size x:" << na << " y:"
				<< nb << " z:" << nc << endl;
		cout << "New grid created with " << gridp << " positions.\n";
		cout << "Projected " << hit << " out of " << atnumber
				<< " total atoms into Grid..." << endl;
	}
}

// compresses grid to "genom"
void ToolCalcZE::grid2genom(bool verbose) {
	ToolIO newIO;
	genom_compressed = new int[atnumber];
	int gridp_compressed = 0;
	//vaccum shifted to account for zero vaccum = 1
	int vacuum_counter = 1;
	int total_vacuum = 0;
	int counter = 0;
	for (int i = 0; i < gridp; i++) {
		if (genom[i] != 0) {
			//cout <<setw(4)<< " V:" << vacuum_counter << " C:" << genom[i];
			//cout << " "<<i;
			//			if (counter % 10 == 0) {
			//				cout << endl;
			//			}
			if (genom[i] > 0) {
				genom_compressed[counter] = vacuum_counter;
			} else {
				genom_compressed[counter] = vacuum_counter * -1;
			}
			//cout << " "<<atoms[gridp_compressed];
			counter++;
			total_vacuum = total_vacuum + vacuum_counter;
			vacuum_counter = 1;
			gridp_compressed++;
		} else {
			vacuum_counter++;
		}
	}
	genom_compressed[counter] = vacuum_counter;
	//cout << "counter:" << counter << endl;
	//cout << " V:" << vacuum_counter << endl;
	total_vacuum = total_vacuum + vacuum_counter;
	if (gridp_compressed != atnumber) {
		cout << "\ngrid2genom: something went terribly wrong: number of atoms:"
				<< atnumber << " does not coincide with genes:"
				<< gridp_compressed << " Exiting..." << endl;
		newIO.fractoFile2(atnumber, energy, atoms, *xyz, core_shell, nr_shells,
				true, Lx, Ly, Lz);
		exit(1);
	}
	if (verbose == true) {
		cout << "grid2genom: atoms:" << gridp_compressed << " total vacuum:"
				<< total_vacuum << endl;
	}
}

//decompresses genom
void ToolCalcZE::genom2grid() {
	int vacuum_counter = 0;
	int counter = 0;

	for (int i = 0; i < (atnumber); i++) {
		//cout << atoms[i]<<endl;
		vacuum_counter = genom_compressed[i];
		if (vacuum_counter > 0) {
			while (vacuum_counter > +1) {
				cout << "0 ";
				genom[counter] = 0;
				counter++;
				vacuum_counter -= 1;
			}
			genom[counter] = 2;
			counter++;
			//cout << "2 ";
		}
		if (vacuum_counter < 0) {
			while (vacuum_counter < -1) {
				//cout << "0 ";
				genom[counter] = 0;
				counter++;
				vacuum_counter += 1;
			}
			genom[counter] = -2;
			counter++;
			//cout << "-2 ";
		}
	}
	cout << endl;

}

// creates xyzmatrix and orders it, oxygens last
void ToolCalcZE::grid2xyz(bool verbose) {
	double cubex = 0, cubey = 0, cubez = 0;
	double spacex = 0, spacey = 0, spacez = 0;
	int segment = 0;
	int pos = 0;
	int atom1 = 0;
	spacex = Lx / (double) gx;
	spacey = Ly / (double) gy;
	spacez = Lz / (double) gz;
	//check if frozen
	//	for (int i = 0; i < atnumber * 3; i++) {
	//		if ((*moveMat)(i) == 1.0) {
	//		cout<< "Frozen atoms specified, will lead to problems with grid projections!"
	//				<< endl;
	//		cout << (*moveMat)(i) << endl;
	//		}
	//	}

	for (int a = 0; a < gx; a++) {
		cubex = a * spacex;
		for (int b = 0; b < gy; b++) {
			cubey = b * spacey;
			for (int c = 0; c < gz; c++) {
				cubez = c * spacez;
				//cout << "Actual seg: "<< segment<< " coords x:" <<cubex
				//							<<" y:"<<cubey <<" z:" <<cubez;
				if (genom[segment] != 0) {
					atom1 = pos / 3;
					//cout << " - Found atom: "<< genom[segment];
					//cout << pos <<endl;
					if (genom[segment] == 2) {
						atom_nr[atom1] = 30;
						q[atom1] = 3.0;
						atoms[atom1] = "Zn";
					}
					if (genom[segment] == 3) {
						atom_nr[atom1] = 13;
						q[atom1] = 3.0;
						atoms[atom1] = "Al";
					}
					if (genom[segment] == -2) {
						atom_nr[atom1] = 8;
						q[atom1] = -2.0;
						atoms[atom1] = "O";
					}
					(*xyz)(pos) = cubex;
					(*xyz)(pos + 1) = cubey;
					(*xyz)(pos + 2) = cubez;
					pos += 3;
				}
				//cout << endl;
				segment++;
			}
		}

	}
	if (verbose == true) {
		cout << "grid2xyz: XYZ-matrix created with " << pos / 3 << " atoms.\n";
	}
	if (pos != atnumber * 3) {
		cout << "Something went terribly wrong, exiting..." << endl;
		exit(1);
	}
}

//sets xyz-coordinates
void ToolCalcZE::orderOx() {
	int atom1;
	int n;
	int temp_atom_nr;
	double temp_q;
	string temp_atoms;
	double temp_x = 0, temp_y = 0, temp_z = 0;
	n = (atnumber - 1) * 3;
	while (n > 1) {
		for (int i = 0; i < n; i++) {
			atom1 = i / 3;
			// if oxygen swap with next atom
			if (i % 3 == 0) {
				if (atom_nr[atom1] == 8) {
					//temp_save
					//					cout << "Swapping...";
					//					cout << atoms[atom1] << " x" << (*xyz)(i) << "y "
					//							<< (*xyz)(i + 1) << "z " << (*xyz)(i + 2)
					//							<< " with: ";
					//					cout << atoms[atom1 + 1] << " x" << (*xyz)(i + 3) << "y "
					//							<< (*xyz)(i + 4) << "z " << (*xyz)(i + 5)
					//							<< endl;
					temp_atom_nr = atom_nr[atom1];
					temp_q = q[atom1];
					temp_atoms = atoms[atom1];
					temp_x = (*xyz)(i);
					temp_y = (*xyz)(i + 1);
					temp_z = (*xyz)(i + 2);
					//get next atom
					atom_nr[atom1] = atom_nr[atom1 + 1];
					q[atom1] = q[atom1 + 1];
					atoms[atom1] = atoms[atom1 + 1];
					(*xyz)(i) = (*xyz)(i + 3);
					(*xyz)(i + 1) = (*xyz)(i + 4);
					(*xyz)(i + 2) = (*xyz)(i + 5);
					//swap oxygen
					atom_nr[atom1 + 1] = temp_atom_nr;
					q[atom1 + 1] = temp_q;
					atoms[atom1 + 1] = temp_atoms;
					(*xyz)(i + 3) = temp_x;
					(*xyz)(i + 4) = temp_y;
					(*xyz)(i + 5) = temp_z;
					//(*xyz)(i)
				}
			}
		}
		n = n - 1;
	}

}

ToolCalcZE::~ToolCalcZE() {

}
