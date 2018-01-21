#include "ToolCalcZE.h"
#include "ToolIO.h"
#include <iostream>
#include <Eigen/Core>
#include <cmath>

using namespace std;

ToolCalcZE::ToolCalcZE() {
	Lx = 0.0;
	Ly = 0.0;
	Lz = 0.0;

}

void ToolCalcZE::printSegments(bool unitsBohr) {
	cout << setw(5) << "atom" << setw(8) << "segment" << setw(10) << "x"
			<< setw(10) << "y" << setw(10) << "z" << endl;
	cout << fixed << setprecision(4);
	double conv = 1.0;
	if (unitsBohr)
		conv = ANG2BOHR;
	for (int i = 0; i < nseg * 3; i += 3) {
		cout << setw(5) << satom(i / 3) << setw(8) << (i / 3);
		cout << setw(10) << segments(i) * conv << setw(10)
				<< segments(i + 1) * conv << setw(10) << segments(i + 2) * conv
				<< endl;
	}

}

void ToolCalcZE::defineCubeSize(bool roundit) {
	//reshape matrix - just a view
	Eigen::MatrixXd M1(segments);

	EigenCoords M2(M1.data(), M1.cols() / 3, 3);
	Eigen::Vector3d maxv = M2.colwise().maxCoeff();
	Eigen::Vector3d minv = M2.colwise().minCoeff();
	double delta = 1.0;
	if (Lx < 1E-15) {
		Lx = maxv(0) - minv(0) + delta;
	}
	if (Ly < 1E-15) {
		Ly = maxv(1) - minv(1) + delta;
	}
	if (Lz < 1E-15) {
		Lz = maxv(2) - minv(2) + delta;
	}
	if (roundit) {
		Lx = ceil(Lx);
		Ly = ceil(Ly);
		Lz = ceil(Lz);
	}
	//make it square
	Eigen::Vector3f L;
	L << Lx, Ly, Lz;
	Lx = L.maxCoeff();
	Ly = L.maxCoeff();
	Lz = L.maxCoeff();

}

Eigen::Vector3d ToolCalcZE::getCubeOrigin() {
	//??
	Eigen::MatrixXd M1(segments);
	EigenCoords M2(M1.data(), M1.cols() / 3, 3);
	Eigen::Vector3d minv = M2.colwise().minCoeff();
	double delta = 1E-1;
	int n = minv.size();
	return minv - delta * Eigen::Vector3d::Ones(n);
}

void ToolCalcZE::showgrid() {
	EigenCoords M2(gridseg.data(), gridseg.rows() / 3, 3);

	Eigen::IOFormat NumpyFmt(Eigen::FullPrecision, 0, ", ", ",\n", "[", "]",
			"[", "]");
	cout << M2.format(NumpyFmt) << endl;
}

void ToolCalcZE::seg2qhull() {
	cout << "3 cosmo segments" << endl;
	cout << nseg << endl;
	Eigen::MatrixXd M1(segments);
	EigenCoords M2(M1.data(), M1.cols() / 3, 3);
	//Eigen::IOFormat QHullFmt(Eigen::FullPrecision, 0, ", ", ",\n", "", "", "", "");
	//cout<<M2.format(QHullFmt)<<endl;
	cout << M2 << endl;
}

void ToolCalcZE::moveSystem2octantI() {
	//move segments and atom coordinates to first octant
	double x, y, z;
	double delta = voxelstep;
	Eigen::MatrixXd M1(segments);
	EigenCoords M2(M1.data(), M1.cols() / 3, 3);
	//find most negative x,y,z
	Eigen::Vector3d minv = M2.colwise().minCoeff();
	for (int i = 0; i < nseg * 3; i += 3) {
		segments(i) = segments(i) - minv(0) + delta;
		segments(i + 1) = segments(i + 1) - minv(1) + delta;
		segments(i + 2) = segments(i + 2) - minv(2) + delta;

	}

	for (int i = 0; i < atnumber; ++i) {
		x = (*xyz)(3 * i);
		y = (*xyz)(3 * i + 1);
		z = (*xyz)(3 * i + 2);
		(*xyz)(3 * i) = x - minv(0) + delta;
		(*xyz)(3 * i + 1) = y - minv(1) + delta;
		(*xyz)(3 * i + 2) = z - minv(2) + delta;
	}
}

void ToolCalcZE::rotateSegments(int axis, double phi) {
	cout << "Rotation of segments around axis:" << axis << " using phi=" << phi
			<< "\n";
	phi = phi * 3.14159265 / 180.0;
	Eigen::MatrixXd R(3, 3);
	if (axis == 0) {
		R << 1.0, 0.0, 0.0, 0.0, cos(phi), -sin(phi), 0.0, sin(phi), cos(phi);
	} else if (axis == 1) {
		R << cos(phi), 0.0, sin(phi), 0.0, 1.0, 0.0, -sin(phi), 0.0, cos(phi);
	} else if (axis == 2) {
		R << cos(phi), -sin(phi), 0.0, sin(phi), cos(phi), 0.0, 0.0, 0.0, 1.0;
	}
	//cout << R << endl;
	//rotate segments
	double x, y, z;
	for (int i = 0; i < nseg * 3; i += 3) {
		x = R(0, 0) * segments(i) + R(0, 1) * segments(i + 1)
				+ R(0, 2) * segments(i + 2);
		y = R(1, 0) * segments(i) + R(1, 1) * segments(i + 1)
				+ R(1, 2) * segments(i + 2);
		z = R(2, 0) * segments(i) + R(2, 1) * segments(i + 1)
				+ R(2, 2) * segments(i + 2);

		segments(i) = x;
		segments(i + 1) = y;
		segments(i + 2) = z;
	}
}

void ToolCalcZE::seg2voxel(bool verbose) {
	int mode = 0;
	cout << "\nCreating voxel grid..." << endl;
	//cout << " atnumber:" << atnumber << endl;
	int at;
	double x, y, z, sx, sy, sz, dist, q, sigma, area;
	double meandist[atnumber];
	vector<double> maxsig(atnumber);

	int segcount[atnumber];
	memset(meandist, 0.0, atnumber * sizeof(double));
	memset(segcount, 0, sizeof segcount);
	//preparation loop, only necessary for fast algo without checking for segments
	for (int i = 0; i < nseg * 3; i += 3) {
		//get segment positions
		sx = segments(i);
		sy = segments(i + 1);
		sz = segments(i + 2);
		//get atom number
		at = satom(i / 3);
		//get atom position
		x = (*xyz)(3 * at);
		y = (*xyz)(3 * at + 1);
		z = (*xyz)(3 * at + 2);
		//mean distance between segments and atom
		dist = sqrt(pow(x - sx, 2) + pow(y - sy, 2) + pow(z - sz, 2));
		meandist[at] = meandist[at] + dist;
		segcount[at] += 1;

		sigma = scharge(i / 3);
		if (abs(sigma) > abs(maxsig[at]))
			maxsig[at] = sigma;
	}
	for (int i = 0; i < atnumber; ++i) {
		meandist[i] = meandist[i] / segcount[i];
		//cout << "atom:" << i << " MeanD:" << meandist[i] << endl;
		//cout << "max.sigma:" << i << " :" << maxsig[i] << endl;

	}
	//ToolIO::moltoFile(this);
	moveSystem2octantI();
	//ToolIO::moltoFile(this);
	Eigen::MatrixXd M1(segments);
	EigenCoords M2(M1.data(), M1.cols() / 3, 3);

	//Eigen::IOFormat NumpyFmt(Eigen::FullPrecision, 0, ", ", ",\n", "[", "]",
	//		"[", "]");
	//cout << M2.format(NumpyFmt) << endl;

	//ToolIO::eigencoord2file(M2,"coords.csv");
	defineCubeSize();
	double side_length = voxelstep;
	int maxi = (int) Lx / side_length;
	int maxj = (int) Ly / side_length;
	int maxk = (int) Lz / side_length;

	//voxel 1D vector
	Eigen::VectorXd voxels = Eigen::VectorXd::Zero(maxi * maxj * maxk);
	//voxel occupations for sigma averaging
	Eigen::VectorXi occ = Eigen::VectorXi::Zero(maxi * maxj * maxk);

	//cout << "no. of voxels in cube:" << voxels.rows() << endl;
	double xpos = 0.0;
	double ypos = 0.0;
	double zpos = 0.0;

	int pos_counter = 0;
	//for every grid point check if voxels are on/off
	//fast mode 0, we shoud iterate only locally e.g. in 10x10x10 cubes starting with icheck!
	if (mode == 0) {
		cout<<"atoms:"<<atnumber<<endl;
		for (int n = 0; n < atnumber; ++n) {
			//atom by atom
			for (int l = 0; l < nseg * 3; l += 3) {
				//get segment positions
				sx = segments(l);
				sy = segments(l + 1);
				sz = segments(l + 2);
				at = satom(l / 3);
				if (at!=n) continue;
//			x = (*xyz)(3 * at);
//			y = (*xyz)(3 * at + 1);
//			z = (*xyz)(3 * at + 2);
//			//get norm vector
//			double nx = sx - x;
//			double ny = sy - y;
//			double nz = sz - z;
				int icheck = sx / side_length - 0.5;
				int jcheck = sy / side_length - 0.5;
				int kcheck = sz / side_length - 0.5;
				for (int m = 0; m < nseg * 3; m += 3) {
					sx = segments(m);
					sy = segments(m + 1);
					sz = segments(m + 2);
					at = satom(m / 3);
					if (at!=n) continue;
					//make little cubes around central segment
					int csize = 5;
					for (int i = -csize;
							(icheck + i) < maxi && (icheck + i) >= 0
									&& i <= csize; ++i) {

						xpos = (icheck + i) * side_length + 0.5 * side_length;
						//cout<<"icheck: "<<icheck<<"i:"<<i<<"xpos:"<<xpos<<"maxi:"<<maxi<<endl;
						for (int j = -csize;
								(jcheck + j) < maxj && (jcheck + j) >= 0
										&& j <= csize; ++j) {
							ypos = (jcheck + j) * side_length
									+ 0.5 * side_length;
							for (int k = -csize;
									(kcheck + k) < maxk && (kcheck + k) >= 0
											&& k <= csize; ++k) {
								zpos = (kcheck + k) * side_length
										+ 0.5 * side_length;
								pos_counter = ((icheck + i) * maxi
										+ (jcheck + j)) * maxj + (kcheck + k);
								dist = sqrt(
										pow(xpos - sx, 2) + pow(ypos - sy, 2)
												+ pow(zpos - sz, 2));
								//now should average again!
								if (dist <= 3 * side_length) {
									occ(pos_counter) += 1;
									q = scharge(l / 3);
									area = sarea(l / 3);
									sigma = 1.0 + 1.0 * q / area;
									//moving average
									voxels(pos_counter) = voxels(pos_counter)
											- voxels(pos_counter)
													/ (double) occ(pos_counter);
									voxels(pos_counter) = voxels(pos_counter)
											+ sigma / (double) occ(pos_counter);
									voxels(pos_counter) = pos_counter;
									voxels(pos_counter) = n+1;
								}
							}
						}
					}
				}
			}
		}
	//simple segment mode
	} else if (mode == 1) {
		for (int i = 0; i < maxi; ++i) {
			xpos = i * side_length + 0.5 * side_length;
			for (int j = 0; j < maxj; ++j) {
				ypos = j * side_length + 0.5 * side_length;
				for (int k = 0; k < maxk; ++k) {
					zpos = k * side_length + 0.5 * side_length;
					//check which segment is there
					sigma = 0.0;
					for (int l = 0; l < nseg * 3; l += 3) {
						//get segment positions
						sx = segments(l);
						sy = segments(l + 1);
						sz = segments(l + 2);
						at = satom(l / 3);
						dist = sqrt(
								pow(xpos - sx, 2) + pow(ypos - sy, 2)
										+ pow(zpos - sz, 2));
						if (dist <= 3 * side_length) {
							q = scharge(l / 3);
							area = sarea(l / 3);
							sigma = 1.0 + 10.0 * q / area;
							sigma = pos_counter;
							voxels(pos_counter) = sigma;
							break;
						}
					}
					voxels(pos_counter) = sigma;
					pos_counter += 1;
				}
			}
		}
		//mean dist mode
	} else {
		for (int i = 0; i < maxi; ++i) {
			xpos = i * side_length + 0.5 * side_length;
			for (int j = 0; j < maxj; ++j) {
				ypos = j * side_length + 0.5 * side_length;
				for (int k = 0; k < maxk; ++k) {
					zpos = k * side_length + 0.5 * side_length;
					for (int l = 0; l < atnumber; ++l) {
						x = (*xyz)(3 * l);
						y = (*xyz)(3 * l + 1);
						z = (*xyz)(3 * l + 2);
						dist = sqrt(
								pow(x - xpos, 2) + pow(y - ypos, 2)
										+ pow(z - zpos, 2));
						//check if we are at a grid position outside or inside of molecule
						// surface
						if (dist < meandist[l]) {
							if (abs(dist - meandist[l]) < side_length
									&& xpos < 2.5) {
								//
								//cout << "JUHU surface:" << 1.0 + maxsig[l]
								//		<< " dist:" << dist1 << " mean dist:"
								//		<< meandist[l] << endl;
								voxels(pos_counter) = 1.0 + maxsig[l];
							} else {
								//switch on if inside
								voxels(pos_counter) = 0.0;
							}
							break;
						} else {
							//switch off if outside
							voxels(pos_counter) = 0.0;
						}
					}
					pos_counter += 1;
				}
			}
		}
	}

	cout << "Lx,Ly,Lz: ";
	cout << setprecision(4);
	cout << Lx << " " << Ly << " " << Lz << endl;
	cout << "Stepsize: " << side_length << endl;
	cout << "Counter : " << pos_counter << endl;
	cout << "Voxels  : " << voxels.sum() << endl;

	//cout << "Density : " << voxels.sum() / ((double) voxels.size()) << endl;
	cout << "Name    :" << jobname << endl;
	ToolIO::eigencoord2file(voxels, jobname + ".csv");
	//save array dimension
	//std::ofstream file(jobname+"_dim.csv");
	//if (file.is_open()) {
	//	file << "# voxel dimension" << endl;
	//	file << maxi << "," << maxj << "," << maxk << endl;
	//}
	//file.close();

}

void ToolCalcZE::seg2grid(int na, int nb, int nc, bool verbose) {
	//need list of grid points?
	//C   2 4 0
	//C   1 4 9
	gridseg = Eigen::MatrixXd::Zero(nseg * 3, 1);
	defineCubeSize();
	gx = na;
	gy = nb;
	gz = nc;
	gridp = na * nb * nc;
	double cubex = 0.0, cubey = 0.0, cubez = 0.0;
	double spacex = 0.0, spacey = 0.0, spacez = 0.0;
	spacex = Lx / (double) na;
	spacey = Ly / (double) nb;
	spacez = Lz / (double) nc;
	int hit = 0;

	Eigen::Vector3d cubeo = getCubeOrigin();
	for (int i = 0; i < nseg * 3; i += 3) {
		cubex = cubeo(0) - spacex;
		for (int a = 0; a < na; a++) {
			cubex = cubex + spacex;
			cubey = cubeo(1) - spacey;
			for (int b = 0; b < nb; b++) {
				cubey = cubey + spacey;
				cubez = cubeo(2) - spacez;
				for (int c = 0; c < nc; c++) {
					cubez = cubez + spacez;

					if (((segments)(i) > cubex)
							&& ((segments)(i) < (cubex + spacex))
							&& (((segments)(i + 1) > cubey)
									&& ((segments)(i + 1) < (cubey + spacey)))
							&& (((segments)(i + 2) > cubez)
									&& ((segments)(i + 2) < (cubez + spacez)))) {
						hit++;
						gridseg(i) = a;
						gridseg(i + 1) = b;
						gridseg(i + 2) = c;
					}
				}
			}
		}

	}
	//cout << "Actual seg: "<< segment<< " coords x:" <<cubex+spacex<<" y:"
	//		<<cubey+spacey <<" z:" <<cubez+spacez<<endl;
	if (verbose == true) {
		cout << endl << "seg2grid:" << endl << "Grid size x:" << na << " y:"
				<< nb << " z:" << nc << endl;
		cout << "New grid created with " << gridp << " positions.\n";
		cout << "Projected " << hit << " out of " << nseg
				<< " total segments into grid coordinates..." << endl;
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

//sets xyz-coordinates, puts oxygen to the end??
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

