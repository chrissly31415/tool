//#include <omp.h>
#include "ToolCalcPP.h"
#include "ToolIO.h"
#include <time.h>

using namespace std;

ToolCalcPP::ToolCalcPP() {
	step = 0.1;//PP calculations
	cutoff = 50.;
	max_displacement = 0.2;
	threshhold = 0.01;
	maxiter = 500;
	global_maxiter = 20;
	nr_reset = 2;
	temperature = 273;
	k = 8.617343 / 10000.0; //for PP
	nr_species = 0;
	interactions = 0;
	TotalQ = 0;
	nr_shells = 0;
	nr_cores = 0;
	ewald_thresh = 0.01;
	scutoff = 0.6;
	springc = 74.92;
}

//energies and gradients according to Catlow-Lewis potential
void ToolCalcPP::E(bool pol, bool verbose, bool gcalc) {
	double lxyz[atnumber * 3];
	double lgrad[atnumber * 3];
	double a[3];
	double b[3];
	double dist = 0;
	double en = 0;
	double enq = 0;
	double esum = 0;
	double esumq = 0;
	double ngrad = 0;
	double ncgrad = 0;
	double lA = 0;
	double lrho = 1;
	double q1 = 0;
	double q2 = 0;
	int atom1 = 0;
	int atom2 = 0;
	//energyFunc EF1;
	int numpos = atnumber * 3;
	//cout <<"nupos:" << numpos;
	for (int i = 0; i < numpos; i++) {
		lxyz[i] = (*xyz)(i);
		//initialisation of grad!!
		lgrad[i] = 0;
	}
	for (int i = 0; i < numpos; i += 3) {
		a[0] = lxyz[i];
		a[1] = lxyz[i + 1];
		a[2] = lxyz[i + 2];
		for (int j = i + 3; j < numpos; j += 3) {
			b[0] = lxyz[j];
			b[1] = lxyz[j + 1];
			b[2] = lxyz[j + 2];
			b[0] = a[0] - b[0];
			b[1] = a[1] - b[1];
			b[2] = a[2] - b[2];
			dist = pow(b[0], 2) + pow(b[1], 2) + pow(b[2], 2);
			dist = sqrt(dist);
			if (dist > cutoff) {
				//cout << "Cut off too big...";
				continue;
			}
			atom1 = i / 3;
			atom2 = j / 3;
			//ENERGY, these arrays are definitily too big
			lA = paraMat[getInteraction(atom_nr[atom1], atom_nr[atom2])][0];
			lrho = paraMat[getInteraction(atom_nr[atom1], atom_nr[atom2])][1];
			//			lA = A[atom1][atom2];
			//			lrho = rho[atom1][atom2];
			en = lA * exp(-1 * dist / lrho);
			en = en / HARTREE2EV;
			//cout << "E: " <<en << endl;
			esum = esum + en;
			q1 = q[atom1];
			q2 = q[atom2];
			//1.889684482 should be changed to 1.88972598
			enq = q1 * q2 / (dist * ANG2BOHR);
			//cout << "Enq: " <<enq << endl;
			esumq = esumq + enq;
			//GRADIENT
			ngrad = lA * -1 / lrho * exp(-dist / lrho);
			// grad in a.u.
			ngrad = ngrad / (HARTREE2EV * ANG2BOHR);
			ncgrad = -1 * q1 * q2 / pow(dist * ANG2BOHR, 2);
			ngrad = ngrad + ncgrad;
			//cout << "ngrad: " <<ngrad<<endl;
			//cout << "dist:" <<dist<<endl;
			b[0] = b[0] / dist * ngrad;
			b[1] = b[1] / dist * ngrad;
			b[2] = b[2] / dist * ngrad;
			//			cout << b[0]<<endl;
			//			cout << b[1]<<endl;
			//			cout << b[2]<<endl;
			lgrad[i] = lgrad[i] + b[0];
			lgrad[i + 1] = lgrad[i + 1] + b[1];
			lgrad[i + 2] = lgrad[i + 2] + b[2];
			lgrad[j] = lgrad[j] - b[0];
			lgrad[j + 1] = lgrad[j + 1] - b[1];
			lgrad[j + 2] = lgrad[j + 2] - b[2];
			//			cout <<grad[i]<<endl;
			//			cout <<grad[i+1]<<endl;
			//			cout <<grad[i+2]<<endl;
		}
	}
	energy = esum + esumq;
	//cout << "Energy: " << esum*27.211 <<endl;
	//cout << "Coulomb: " << esumq*27.211 <<endl;
	//grad = new LaVectorDouble(lgrad,numpos);
	for (int i = 0; i < numpos; i++) {
		(*grad)(i) = lgrad[i];
		//which gradients are frozen??
		(*grad)(i) = lgrad[i] * (*moveMat)(i);
		//cout << "grad:"<< (*grad)(i)<<endl;
		//cout << lgrad[i]<<endl;
	}
	//gradnorm = Blas_Norm2((*grad));
	gradnorm = grad->norm();
	grad->maxCoeff(&gradmaxelem);
	//gradmaxelem = Blas_Index_Max((*grad));
	//gradinfnorm = Blas_Norm_Inf((*grad));
	//cout << "Energy: "<<esum <<endl;
	calculated = true;
}

void ToolCalcPP::E_Pol(bool pol, bool verbose, bool gcalc) {
	double lxyz[atnumber * 3];
	double lgrad[atnumber * 3];
	double a[3];
	double b[3];
	double dist = 0, en = 0, enq = 0, esum = 0, esumq = 0, es = 0, sgrad = 0,
			ngrad = 0, ncgrad = 0, lA = 0, lrho = 1;
	double q1 = 0;
	double q2 = 0;
	int atom1 = 0;
	int atom2 = 0;
	double Ospring = 74.92;
	double scutoff = 0.6;
	//energyFunc EF1;
	int numpos = atnumber * 3;
	for (int i = 0; i < numpos; i++) {
		lxyz[i] = (*xyz)(i);
		//initialisation of grad!!
		lgrad[i] = 0;
	}
	if (verbose == true && pol == true) {
		cout << "Using shell-modell for energy calculation...\n";
	}
	bool core1 = false;
	bool core2 = false;
	bool norm1 = true;
	bool norm2 = true;
	bool shel1 = false;
	bool shel2 = false;
	for (int i = 0; i < numpos; i += 3) {
		atom1 = i / 3;
		q1 = q[atom1];
		if (core_shell[atom1] == "CORE" && pol == false) {
			//skip core atoms if no shell model is wanted
			continue;
		} else {
			if (core_shell[atom1] == "CORE") {
				//cout << "O_Ci"<<endl;
				core1 = true;
				norm1 = false;
			}
			if (core_shell[atom1] == "SHEL") {
				//cout << "O_Ci"<<endl;
				shel1 = true;
				norm1 = false;
			}
		}
		if (atom_nr[atom1] == 8 && pol == false) {
			//cout << "O_Ci"<<endl;
			q1 = -2.0;
		}
		a[0] = lxyz[i];
		a[1] = lxyz[i + 1];
		a[2] = lxyz[i + 2];
		for (int j = i + 3; j < numpos; j += 3) {
			//resetting of energy terms
			en = 0;
			es = 0;
			ngrad = 0;
			sgrad = 0;
			atom2 = j / 3;
			q2 = q[atom2];
			if (core_shell[atom2] == "CORE" && pol == false) {
				//cout << "O_Ci"<<endl;
				continue;
			} else {
				if (core_shell[atom2] == "CORE") {
					//cout << "O_Ci"<<endl;
					core2 = true;
					norm2 = false;
				}
				if (core_shell[atom2] == "SHEL") {
					//cout << "O_Ci"<<endl;
					shel2 = true;
					norm2 = false;
				}
			}
			if (atom_nr[atom2] == 8 && pol == false) {
				//cout << "O_Ci"<<endl;
				//cout << atom_nr[atom1]<< endl;
				q2 = -2.0;
			}
			b[0] = lxyz[j];
			b[1] = lxyz[j + 1];
			b[2] = lxyz[j + 2];
			b[0] = a[0] - b[0];
			b[1] = a[1] - b[1];
			b[2] = a[2] - b[2];
			dist = pow(b[0], 2) + pow(b[1], 2) + pow(b[2], 2);
			dist = sqrt(dist);
			if (dist > cutoff) {
				//cout << "Cut off too big...";
				continue;
			}
			//SPRING POTENTIAL
			if (((core1 == true && shel2 == true) || (shel1 == true && core2
					== true)) && (dist < scutoff)) {
				//energy
				es = 0.5 * Ospring * pow(dist, 2);
				es = es / HARTREE2EV;
				en = 0;
				enq = 0;
				//gradient
				sgrad = Ospring * dist / (HARTREE2EV * ANG2BOHR);
				ngrad = 0;
				ncgrad = 0;
				//cout << "Spring energy" << es*HARTREE2EV << endl;
			}
			//PAIR POTENTIAL
			if (shel1 == true && shel2 == true) {
				//Both are SHELLs!
				lA = paraMat[getInteraction(atom_nr[atom1], atom_nr[atom2])][0];
				lrho
						= paraMat[getInteraction(atom_nr[atom1], atom_nr[atom2])][1];
				en = lA * exp(-1 * dist / lrho);
				en = en / HARTREE2EV;
				es = 0;
				//gradient
				ngrad = lA * -1 / lrho * exp(-dist / lrho);
				ngrad = ngrad / (HARTREE2EV * ANG2BOHR);
				sgrad = 0;
			}
			// a few conditiona if there is no potential
			if (core1 == true && core2 == true) {
				en = 0;
				es = 0;
				ngrad = 0;
				sgrad = 0;
			}
			//remote core-shell interaction
			if (((core1 == true && shel2 == true) || (shel1 == true && core2
					== true)) && (dist > scutoff)) {
				en = 0;
				es = 0;
				ngrad = 0;
				sgrad = 0;
			}
			en = en + es;
			esum = esum + en;
			//COULOMB INTERACTION
			//NO Coulomb if core interacts with closest shell
			if (((core1 == true && shel2 == true) || (core2 == true && shel1
					== true)) && (dist < scutoff)) {
				enq = 0;
				ncgrad = 0;
			} else {

				enq = q1 * q2 / (dist * ANG2BOHR);
				ncgrad = -1 * q1 * q2 / pow(dist * ANG2BOHR, 2);
				//cout << atom_nr[atom1] << " " << atom_nr[atom2] << " " << enq << endl;
			}
			esumq = esumq + enq;
			//cout << "esumq : " << esumq << endl;
			ngrad = ngrad + ncgrad + sgrad;
			if (dist > 0.0001) {
				b[0] = b[0] / dist * ngrad;
				b[1] = b[1] / dist * ngrad;
				b[2] = b[2] / dist * ngrad;
			} else {
				b[0] = 0;
				b[1] = 0;
				b[2] = 0;
			}
			lgrad[i] = lgrad[i] + b[0];
			lgrad[i + 1] = lgrad[i + 1] + b[1];
			lgrad[i + 2] = lgrad[i + 2] + b[2];
			lgrad[j] = lgrad[j] - b[0];
			lgrad[j + 1] = lgrad[j + 1] - b[1];
			lgrad[j + 2] = lgrad[j + 2] - b[2];
			core2 = false;
			shel2 = false;
			norm2 = true;
		}
		core1 = false;
		shel1 = false;
		norm1 = true;
	}
	if (verbose == true) {
		cout.precision(8);
		cout << "Interatomic potentials:\t" << esum * 27.211 << " eV" << endl;
		cout << "Coulomb interaction:\t" << esumq * 27.211 << " eV" << endl;
	}
	energy = esum + esumq;
	//cout << "Energy: " << esum*27.211 <<endl;
	//cout << "Coulomb: " << esumq*27.211 <<endl;
	//grad = new LaVectorDouble(lgrad,numpos);
	for (int i = 0; i < numpos; i++) {
		(*grad)(i) = lgrad[i];
		//which gradients are frozen??
		(*grad)(i) = lgrad[i] * (*moveMat)(i);
		//cout << "grad:"<< (*grad)(i)<<endl;
		//cout << lgrad[i]<<endl;
	}
	//gradnorm = Blas_Norm2((*grad));
	gradnorm = grad->norm();
	grad->maxCoeff(&gradmaxelem);
	//gradmaxelem = Blas_Index_Max((*grad));
	//gradinfnorm = Blas_Norm_Inf((*grad));
	//cout << "Energy: "<<esum <<endl;
	calculated = true;
}
//PARALLEL NOT WORKING 100% correctly!
//emp. pair potential under periodic boundary conditions
//TO DO FOR SPEED UP: Symmetry, linked list, FFT, Nearest image trick, combine loops  ...
void ToolCalcPP::E_periodic(bool pol, bool verbose, bool gcalc) {
	bool ewald_converged = false;
	//cout<<"+++++++OPENMP switched ON++++++++++++++\n";
	#ifdef OMP_H
	omp_set_num_threads(1);
	//#pragma omp parallel
	//cout<<"OPENMP switched ON, NUMTHREADS:"<<omp_get_num_threads()<<"\n";
	#endif
	//ToolIO newIO;
	//shell counters
	int nmax_reci = 50;
	int nmax_real = 1;
	//Only orthorombic cells at the moment
	//Optimisations: prepare cell parameters in separate function
	double lxyz[atnumber * 3];
	double lgrad[atnumber * 3];
	double a[3];
	double b[3];
	double dist = 0, rsq = 0;
	double en = 0, enq_sr = 0, esum = 0, esumq_sr = 0, enq_f = 0, esumq_f = 0,
			enq_sc = 0, esumq_sc = 0, esum_dipole = 0;
	double ngrad = 0, ncgrad = 0;
	double lA = 0., lrho = 1.;
	double q1 = 0., q2 = 0.;
	int atom1 = 0, atom2 = 0, tmpi = 0;
	//	//orthogonal box
	double L[3];
	//width of gaussian
	double epsilon = 1.0;
	//reciprocal space variables
	double tmp = 0;
	double ksq = 0;
	double rhok_re = 0;
	double rhok_im = 0;
	double rhoki_re = 0;
	double rhoki_im = 0;
	double phr = 0;
	double phi = 0;
	double grad1 = 0, grad2 = 0, grad3 = 0;
	//lattice vectors in real space
	//double ar[3], br[3], cr[3];
	//lattice vector in Fourier space
	//double af[3], bf[3], cf[3];
	double G[3];
	if (verbose == true) {
		cout.precision(2);
		cout << endl << endl << "Periodic calculation with Lx=" << Lx
				<< ", Ly=" << Ly << ", Lz=" << Lz << " angstrom" << endl;
		cout.precision(8);
		cout << "Coulomb interaction calculated via Ewald approach, alpha="
				<< alpha << endl << endl;
	}
	if (nr_shells > 0 || pol == true) {
		cout
				<< "WARNING! Shells have been defined or polarization switched on, wrong energy function used!\n";
		exit(1);
	}
	int numpos = atnumber * 3;
	//setting up local arrays
	#pragma omp parallel for
	for (int i = 0; i < numpos; i++) {
		lxyz[i] = (*xyz)(i) * ANG2BOHR;
		lgrad[i] = 0;
	}
	//Gradient dipole correction + self interaction
	double tmpQX = 0;
	double tmpQY = 0;
	double tmpQZ = 0;
	if (dipole_cor == true) {
		//Energy dipole correction
		setDipole(false);
		esum_dipole = 2*PI * dipnorm * dipnorm * ANG2BOHR * ANG2BOHR / ((2
				+ epsilon) * V);
		//SUM Q*R
		#pragma omp parallel for private (atom2,q2,a) reduction(+:tmpQX,tmpQY,tmpQZ)
		for (int j = 0; j < numpos; j += 3) {
			atom2 = j / 3;
			q2 = q[atom2];
			//atoms in cell, right dipole?
			a[0] = lxyz[j];
			a[1] = lxyz[j + 1];
			a[2] = lxyz[j + 2];
			tmpQX = tmpQX + q2 * a[0];
			tmpQY = tmpQY + q2 * a[1];
			tmpQZ = tmpQZ + q2 * a[2];
		}
	}
	#pragma omp parallel for private (atom1,q1,enq_sc) reduction(+:esumq_sc)
	for (int i = 0; i < numpos; i += 3) {
		//correction for self interaction
		atom1 = i / 3;
		q1 = q[atom1];
		enq_sc = pow(q1, 2);
		esumq_sc = esumq_sc + enq_sc;
		//gradient dipole correction
		if (dipole_cor == true) {
			lgrad[i] = 4*PI * q1 * tmpQX / ((2 + epsilon) * V);
			lgrad[i + 1] = 4*PI * q1 * tmpQY / ((2 + epsilon) * V);
			lgrad[i + 2] = 4*PI * q1 * tmpQZ / ((2 + epsilon) * V);
		}
	}

	//correction for self interaction
	esumq_sc = -1 * alpha / sqrt(PI) * esumq_sc;

	//timing
	timespec res, t1, t2, t3;
	clock_getres(CLOCK_REALTIME, &res);
	clock_gettime(CLOCK_REALTIME, &t1);
	//RECIPROCAL SPACE - Symmetry!!!!?
	double shell_contr = 0.0;
	//iteration of cell counter
	for (int nact = 1; nact <= nmax_reci; nact++) {
		if (ewald_converged == false) {
			//#pragma omp parallel for private (tmp,G,ksq,a,b,atom1,atom2,q1,q2,phr,phi,rhoki_re,rhoki_im) reduction(+:rhok_re,rhok_im,esumq_f,shell_contr,enq_f) schedule(static,CHUNKSIZE)
			#pragma omp parallel for private (tmp,G,ksq,a,b,atom1,atom2,q1,q2,phr,phi,rhoki_re,rhoki_im) reduction(+:rhok_re,rhok_im,esumq_f,shell_contr,enq_f)
			for (int nx = -nact; nx <= nact; nx++) {
				for (int ny = -nact; ny <= nact; ny++) {
					for (int nz = -nact; nz <= nact; nz++) {
						//skip central cell
						if (nx == 0 && ny == 0 && nz == 0) {
							enq_f = 0;
							continue;
						}
						//skip cells previously summed over
						if (abs(nx) == nact || abs(ny) == nact || abs(nz)
								== nact) {
							//Reciprocal space
							G[0] = 2* PI * nx * af[0] + 2* PI * ny * bf[0] + 2
									* PI * nz * cf[0];
							G[1] = 2* PI * nx * af[1] + 2* PI * ny * bf[1] + 2
									* PI * nz * cf[1];
							G[2] = 2* PI * nx * af[2] + 2* PI * ny * bf[2] + 2
									* PI * nz * cf[2];
							//k²
							ksq = G[0] * G[0] + G[1] * G[1] + G[2] * G[2];
							//real and imaginary part
							rhok_re = 0;
							rhok_im = 0;
							//ENERGY
							#pragma omp critical
							for (int i = 0; i < numpos; i += 3) {
								a[0] = lxyz[i];
								a[1] = lxyz[i + 1];
								a[2] = lxyz[i + 2];
								atom1 = i / 3;
								q1 = q[atom1];
								//calculate k * ri
								tmp = G[0] * a[0] + G[1] * a[1] + G[2] * a[2];
								rhoki_re = q1 * cos(tmp);
								//rhoki_im = rhoki_re - PI*0.5;
								rhoki_im = q1 * sin(tmp);
								//calculate rhok
								rhok_re = rhok_re + rhoki_re;
								rhok_im = rhok_im + rhoki_im;
							}
							tmp = 2* PI * (rhok_re * rhok_re + rhok_im
									* rhok_im) / (V * ksq);
							enq_f = tmp * exp(-1* ksq / (4* pow (alpha, 2)));
							esumq_f = esumq_f + enq_f;
							//Gradient Screening
							//#pragma omp critical
							if (gcalc == true && abs((enq_f) * HARTREE2EV)
									< 10E-8) {
								//printf ("ENQ_F: %2.4E\n", enq_f);
							}
							//GRADIENT, second summation neeed
							else {
								if (gcalc == true) {
									//#pragma omp critical
									for (int i = 0; i < numpos; i += 3) {
										phr = 0;
										phi = 0;
										//only phase important
										a[0] = lxyz[i];
										a[1] = lxyz[i + 1];
										a[2] = lxyz[i + 2];
										atom1 = i / 3;
										q1 = q[atom1];
										tmp = G[0] * a[0] + G[1] * a[1] + G[2]
												* a[2];
										phr = cos(tmp);
										phi = -sin(tmp);
										tmp = (phi * rhok_re + phr * rhok_im)
												* q1 * 4*PI / (V * ksq) * exp(
												-1 * ksq / (4 * pow(alpha, 2)));
										//cout << "tmp:" << tmp << endl;
										#pragma omp critical
										{
											lgrad[i] = lgrad[i] + tmp * G[0];
											lgrad[i + 1] = lgrad[i + 1] + tmp
													* G[1];
											lgrad[i + 2] = lgrad[i + 2] + tmp
													* G[2];
										}
									}
								}
							}
						}//end if
						//cout.precision(8);
						//cout << "nx: " << nx << " ny: " << ny << " nz: " << nz << " enq_f:\t" << enq_f << endl;
					} // nz for loop
				} // ny for loop
			}// nx for loop

			shell_contr = esumq_f - shell_contr;

			if (verbose == true) {
				cout << "Shell #" << nact << " contribution, Q recip. part: "
						<< shell_contr * HARTREE2EV << endl;
			}
			////#pragma omp critical
			if (abs(shell_contr * HARTREE2EV < ewald_thresh)) {
				//cout << "We could leave now...";
				ewald_converged = true;
				//#pragma omp flush (ewald_converged)
			}
		}//end if ewald_converged
		shell_contr = esumq_f;
	}// nact for loop

	//timing
	clock_gettime(CLOCK_REALTIME, &t2);
	if (verbose == true) {
		cout << "Reciprocal Part: " << (t2.tv_sec - t1.tv_sec)
				+ (float) (t2.tv_nsec - t1.tv_nsec) / 1000000.0 << " ms\n";
	}
	//REAL SPACE
	shell_contr = 0.0;
	ewald_converged = false;
	for (int nact = 0; nact <= nmax_real; nact++) {
		if (!ewald_converged) {
			//#pragma omp parallel for private (tmp,tmpi,a,b,L,atom1,atom2,q1,q2,enq_sr,en,rsq,dist,lA,lrho) reduction(+:shell_contr,esum,esumq_sr,ngrad,grad1,grad2,grad3) schedule(static,CHUNKSIZE)
			#pragma omp parallel for private (tmp,tmpi,a,b,L,atom1,atom2,q1,q2,enq_sr,en,rsq,dist,lA,lrho) reduction(+:shell_contr,esum,esumq_sr,ngrad,grad1,grad2,grad3)
			for (int nx = -nact; nx <= nact; nx++) {
				////#pragma omp parallel for  reduction(+:esum,esumq_sr,esumq_e,ngrad) schedule(static,CHUNKSIZE)
				for (int ny = -nact; ny <= nact; ny++) {
					////#pragma omp parallel for  reduction(+:esum,esumq_sr,esumq_e,ngrad) schedule(static,CHUNKSIZE)
					for (int nz = -nact; nz <= nact; nz++) {
						//skip cells previously summed over
						if (abs(nx) == nact || abs(ny) == nact || abs(nz)
								== nact) {
							L[0] = nx * ar[0] + ny * br[0] + nz * cr[0];
							L[1] = nx * ar[1] + ny * br[1] + nz * cr[1];
							L[2] = nx * ar[2] + ny * br[2] + nz * cr[2];
							//cout << "L: <"<<L[0]<< " " << L[1] << " " <<L[2] <<">"
							//<< endl;
							//real space part
							for (int i = 0; i < numpos; i += 3) {
								a[0] = lxyz[i] + L[0];
								a[1] = lxyz[i + 1] + L[1];
								a[2] = lxyz[i + 2] + L[2];
								atom1 = i / 3;
								q1 = q[atom1];
								// inner loop
								//#pragma omp critical
								for (int j = 0; j < numpos; j += 3) {
									if (nx == 0 && ny == 0 && nz == 0) {
										//cout << "We are in unit cell..." << endl;
										if (i == j) {
											//cout << "skipping selfinteraction..." << endl;
											continue;
										}
									}

									//hier kommt der teuerste Teil des codes
									atom2 = j / 3;
									//small optimization, shift L to a!!!
									//get distances
									b[0] = a[0] - lxyz[j];
									b[1] = a[1] - lxyz[j + 1];
									b[2] = a[2] - lxyz[j + 2];
									rsq = pow(b[0], 2) + pow(b[1], 2) + pow(
											b[2], 2);
									dist = sqrt(rsq);


									//distance cut off
									if (dist > cutoff) {
										//cout << "Cut off too big...";
										continue;
									}
									//if this comes earliers codes slows down
									tmpi = getInteraction_fast(atom_nr[atom1],atom_nr[atom2]);
									lA = paraMat[tmpi][0];
									lrho = paraMat[tmpi][1];
									//############################################################
									//ENERGY
									//PP Part

									q2 = q[atom2];
									en = 0.5 * lA * exp(-1 * dist / (ANG2BOHR
											* lrho));
									en = en * 1 / HARTREE2EV;
									esum = esum + en;
									//Coulomb Ewald Part
									//Evaluate SHORT RANGE part of Ewald, seems correct
									tmp = erfc_new(alpha * dist);
									enq_sr = 0.5 * q1 * q2 * tmp / dist;
									esumq_sr = esumq_sr + enq_sr;
									//Add Ewald components
									//esumq_e = esumq_sr + esumq_f + esumq_sc;
									//Evaluate conventional charge
									//enq = q1 * q2 / (dist);
									//esumq = esumq + enq;
									//#############################################################
									//Gradient Screening
									//#pragma omp critical
									if (gcalc == true && abs((enq_sr + en)
											* HARTREE2EV) < 10E-8) {
										//printf ("ENQ_SR: %2.4E   EN: %2.4E\n", enq_sr, en);
//																				grad1 = 0;
//																				grad2 = 0;
//																				grad3 = 0;
										//continue;
									}
									//GRADIENTS
									else {
										//#pragma omp critical
										if (gcalc == true) {
											ncgrad = q1 * q2 * ((-2.0 * sqrt(1
													/ PI) * exp(-pow(alpha, 2)
													* rsq) * alpha * dist)
													- tmp) / rsq;
											ngrad = lA * -1 / lrho * exp(-dist
													/ (ANG2BOHR * lrho));
											ngrad = ngrad / (HARTREE2EV
													* ANG2BOHR) + ncgrad;

											grad1 = grad1 + b[0] / dist * ngrad;
											grad2 = grad2 + b[1] / dist * ngrad;
											grad3 = grad3 + b[2] / dist * ngrad;
										}
									}
								}//end loop j
								#pragma omp critical
								if (gcalc == true) {
									lgrad[i] = lgrad[i] + grad1;
									lgrad[i + 1] = lgrad[i + 1] + grad2;
									lgrad[i + 2] = lgrad[i + 2] + grad3;
									grad1 = 0;
									grad2 = 0;
									grad3 = 0;
								}
							}//end loop i
						} //end if abs(nx)==nact || abs(ny)==nact || abs(nz)==nact
						//cout.precision(8);
						//if (nx < 1 && ny < 1 && nz < 1) {
						//					cout << "nx: " << nx << " ny: " << ny << " nz: " << nz
						//							<< " enq_sr:\t" << enq_sr << "en:\t" << en
						//							*HARTREE2EV<< endl;
						//					}
					}// nz for loop
				}// ny for loop
			}// nx for loop
			shell_contr = esumq_sr - shell_contr;
			if (verbose == true) {
				cout.precision(8);
				cout << "Shell #" << nact << " contribution, Q real part: "
						<< setw(12) << shell_contr * HARTREE2EV;
				cout << "\tPP part: " << esum * HARTREE2EV << endl;
			}
			if (abs(shell_contr * HARTREE2EV) < ewald_thresh) {
				ewald_converged = true;
				#pragma omp flush (ewald_converged)
			}
			shell_contr = esumq_sr;
		}
	}// nact for loop

	energy = esum + esumq_sr + esumq_f + esumq_sc;
	//cout << "esum_cor : " << esum_cor*HARTREE2EV <<endl;
	//	cout << endl << " Energy - Conventional:\t\t" << esumq *HARTREE2EV << endl;


	//timing
	clock_gettime(CLOCK_REALTIME, &t3);
	if (verbose == true) {
		cout << "Real Part: " << (t3.tv_sec - t1.tv_sec) + (float) (t3.tv_nsec
				- t1.tv_nsec) / 1000000.0 << " ms" << endl;
	}
	if (verbose == true) {
		cout.precision(8);
		cout << endl << "Energy - pair potential\t\t:	" << setw(12) << esum
				* HARTREE2EV << " eV" << endl;
		cout << endl << "Energy - Coulomb\t\t:	" << setw(12) << (esumq_sr
				+ esumq_f + esumq_sc) * HARTREE2EV << " eV";
		cout << endl << "Ewald - real space\t\t:	" << setw(12) << esumq_sr
				* HARTREE2EV << " eV" << endl;
		cout << "Ewald - reciprocal space\t:	" << setw(12) << (esumq_f)
				* HARTREE2EV << " eV" << endl;
		cout << "Ewald - self correction\t\t:	" << setw(12) << esumq_sc
				* HARTREE2EV << " eV" << endl << endl;
		cout << "#Total energy\t\t\t:	" << setw(12) << energy * HARTREE2EV
				<< " eV" << endl << endl;
		//energy=esum+esumq_sr+esumq_f+esumq_sc+esum_dipole;
		cout << "Energy - dipole correction\t:	" << setw(12) << esum_dipole
				* HARTREE2EV << " eV" << endl;
		cout << "#Total energy (+dipole corr.)\t:	" << setw(12) << (energy
				+ esum_dipole) * HARTREE2EV << " eV" << endl;
	}
	if (gcalc == true) {
		//#pragma omp parallel for
		for (int i = 0; i < numpos; i++) {
			(*grad)(i) = lgrad[i];
			//which gradients are frozen??
			(*grad)(i) = lgrad[i] * (*moveMat)(i);
			//cout << "grad:"<< (*grad)(i)<<endl;
			//cout << lgrad[i]*27.211*1.88972598<<endl;
		}
		//gradnorm = Blas_Norm2((*grad));
		gradnorm = grad->norm();
		//(gradmaxelem = Blas_Index_Max((*grad));
		grad->maxCoeff(&gradmaxelem);
		//gradinfnorm = Blas_Norm_Inf((*grad));
	}
	calculated = true;
	//delete newIO;
}
//this is a large cell version, with cutoff < L/2
//Should give significant speed up for reciprocal part, not implemented et
void ToolCalcPP::E_periodic_LC(bool pol, bool verbose, bool gcalc) {
	bool ewald_converged = false;
	//omp_set_num_threads(nproc);
	ToolIO newIO;
	//image counters
	int nmax = 50;
	//Only orthorombic cells at the moment
	//Optimisations: prepare cell paramters in separate function
	double lxyz[atnumber * 3];
	double lgrad[atnumber * 3];
	double a[3];
	double b[3];
	double dist = 0, rsq = 0;
	double en = 0, enq_sr = 0, esum = 0, esumq_sr = 0, enq_f = 0, esumq_f = 0,
			enq_sc = 0, esumq_sc = 0, esum_dipole = 0;
	double ngrad = 0, ncgrad = 0;
	double lA = 0., lrho = 1.;
	double q1 = 0., q2 = 0.;
	int atom1 = 0, atom2 = 0, tmpi = 0;
	//	//orthogonal box
	double L[3];
	//width of gaussian
	double epsilon = 1;
	//reciprocal space variables
	double tmp = 0;
	double ksq = 0;
	double rhok_re = 0;
	double rhok_im = 0;
	double rhoki_re = 0;
	double rhoki_im = 0;
	double phr = 0;
	double phi = 0;
	double grad1 = 0, grad2 = 0, grad3 = 0;
	//lattice vectors in real space
	//double ar[3], br[3], cr[3];
	//lattice vector in Fourier space
	//double af[3], bf[3], cf[3];
	double G[3];
	if (verbose == true) {
		cout.precision(2);
		cout << endl << endl << "Periodic calculation with Lx=" << Lx
				<< ", Ly=" << Ly << ", Lz=" << Lz << " angstrom" << endl;
		cout.precision(8);
		cout << "Coulomb interaction calculated via Ewald approach, alpha="
				<< alpha << endl << endl;
	}
	if (nr_shells > 0 || pol == true) {
		cout
				<< "WARNING! Shells have been defined or polarization switched on, wrong energy function used!\n";
		exit(1);
	}
	int numpos = atnumber * 3;
	//setting up local arrays
	//#pragma omp parallel for schedule(static,CHUNKSIZE)
	for (int i = 0; i < numpos; i++) {
		lxyz[i] = (*xyz)(i) * ANG2BOHR;
		lgrad[i] = 0;
	}
	//Gradient dipole correction + self interaction
	double tmpQX = 0;
	double tmpQY = 0;
	double tmpQZ = 0;
	if (dipole_cor == true) {
		//Energy dipole correction
		setDipole(false);
		esum_dipole = 2*PI * dipnorm * dipnorm * ANG2BOHR * ANG2BOHR / ((2
				+ epsilon) * V);
		//SUM Q*R
		//#pragma omp parallel for private (atom2,q2,a) reduction(+:tmpQX,tmpQY,tmpQZ) schedule(static,CHUNKSIZE)
		for (int j = 0; j < numpos; j += 3) {
			atom2 = j / 3;
			q2 = q[atom2];
			//atoms in cell, right dipole?
			a[0] = lxyz[j];
			a[1] = lxyz[j + 1];
			a[2] = lxyz[j + 2];
			tmpQX = tmpQX + q2 * a[0];
			tmpQY = tmpQY + q2 * a[1];
			tmpQZ = tmpQZ + q2 * a[2];
		}
	}
	//#pragma omp parallel for private (atom1,q1,enq_sc) reduction(+:esumq_sc) schedule(static,CHUNKSIZE)
	for (int i = 0; i < numpos; i += 3) {
		//correction for self interaction
		atom1 = i / 3;
		q1 = q[atom1];
		enq_sc = pow(q1, 2);
		esumq_sc = esumq_sc + enq_sc;
		//gradient dipole correction
		//cout << tmpQ[0] << endl;
		//Not stable during optimization
		//Are all atoms inside the cell?!!!!!!!!project them back...
		if (dipole_cor == true) {
			lgrad[i] = 4*PI * q1 * tmpQX / ((2 + epsilon) * V);
			lgrad[i + 1] = 4*PI * q1 * tmpQY / ((2 + epsilon) * V);
			lgrad[i + 2] = 4*PI * q1 * tmpQZ / ((2 + epsilon) * V);
		}
	}

	//correction for self interaction
	esumq_sc = -1 * alpha / sqrt(PI) * esumq_sc;

	//timing
	timespec res, t1, t2, t3;
	clock_getres(CLOCK_REALTIME, &res);
	clock_gettime(CLOCK_REALTIME, &t1);

	//RECIPROCAL SPACE - Symmetry!!!!?
	double shell_contr = 0.0;
	//iteration of cell counter
	for (int nact = 1; nact <= nmax; nact++) {
		if (ewald_converged == false) {
			//#pragma omp parallel for private (tmp,G,ksq,a,b,atom1,atom2,q1,q2,phr,phi,rhoki_re,rhoki_im) reduction(+:rhok_re,rhok_im,esumq_f,shell_contr,enq_f) schedule(static,CHUNKSIZE)
			//cout << "hello\n";
			for (int nx = -nact; nx <= nact; nx++) {
				for (int ny = -nact; ny <= nact; ny++) {
					for (int nz = -nact; nz <= nact; nz++) {
						//skip central cell
						if (nx == 0 && ny == 0 && nz == 0) {
							enq_f = 0;
							continue;
						}
						//skip cells previously summed over
						if (abs(nx) == nact || abs(ny) == nact || abs(nz)
								== nact) {
							//Reciprocal space
							G[0] = 2* PI * nx * af[0] + 2* PI * ny * bf[0] + 2
									* PI * nz * cf[0];
							G[1] = 2* PI * nx * af[1] + 2* PI * ny * bf[1] + 2
									* PI * nz * cf[1];
							G[2] = 2* PI * nx * af[2] + 2* PI * ny * bf[2] + 2
									* PI * nz * cf[2];
							//k²
							ksq = G[0] * G[0] + G[1] * G[1] + G[2] * G[2];
							//real and imaginary part
							rhok_re = 0;
							rhok_im = 0;
							//ENERGY
							////#pragma omp critical
							for (int i = 0; i < numpos; i += 3) {
								a[0] = lxyz[i];
								a[1] = lxyz[i + 1];
								a[2] = lxyz[i + 2];
								atom1 = i / 3;
								q1 = q[atom1];
								//calculate k * ri
								tmp = G[0] * a[0] + G[1] * a[1] + G[2] * a[2];
								rhoki_re = q1 * cos(tmp);
								rhoki_im = q1 * sin(tmp);
								//calculate rhok

								rhok_re = rhok_re + rhoki_re;
								rhok_im = rhok_im + rhoki_im;
							}
							tmp = 2* PI * (rhok_re * rhok_re + rhok_im
									* rhok_im) / (V * ksq);
							enq_f = tmp * exp(-1* ksq / (4* pow (alpha, 2)));
							esumq_f = esumq_f + enq_f;
							//Gradient Screening
							if (gcalc == true && abs((enq_f) * HARTREE2EV)
									< 10E-6) {
								//printf ("ENQ_F: %2.4E\n", enq_f);
								//										grad1 = 0;
								//										grad2 = 0;
								//										grad3 = 0;
								continue;
							}
							//GRADIENT, second summation neeed
							////#pragma omp critical
							if (gcalc == true) {
								////#pragma omp serial
								for (int i = 0; i < numpos; i += 3) {
									phr = 0;
									phi = 0;
									//only phase important
									a[0] = lxyz[i];
									a[1] = lxyz[i + 1];
									a[2] = lxyz[i + 2];
									atom1 = i / 3;
									q1 = q[atom1];
									tmp = G[0] * a[0] + G[1] * a[1] + G[2]
											* a[2];
									phr = cos(tmp);
									phi = -sin(tmp);
									tmp = (phi * rhok_re + phr * rhok_im) * q1
											* 4*PI / (V * ksq) * exp(-1 * ksq
											/ (4 * pow(alpha, 2)));
									//cout << "tmp:" << tmp << endl;
									//#pragma omp critical
									{
										lgrad[i] = lgrad[i] + tmp * G[0];
										lgrad[i + 1] = lgrad[i + 1] + tmp
												* G[1];
										lgrad[i + 2] = lgrad[i + 2] + tmp
												* G[2];
									}
								}
							}

						}//end if
						//cout.precision(8);
						//cout << "nx: " << nx << " ny: " << ny << " nz: " << nz << " enq_f:\t" << enq_f << endl;
					} // nz for loop
				} // ny for loop
			}// nx for loop
			shell_contr = esumq_f - shell_contr;

			if (verbose == true) {
				cout << "Shell #" << nact << " contribution, Q recip. part: "
						<< shell_contr * HARTREE2EV << endl;
			}
			////#pragma omp critical
			if (abs(shell_contr * HARTREE2EV < ewald_thresh)) {
				//cout << "We could leave now...";
				ewald_converged = true;
				//#pragma omp flush (ewald_converged)
			}
		}//end if ewald_converged
		shell_contr = esumq_f;
	}// nact for loop

	//timing
	clock_gettime(CLOCK_REALTIME, &t2);
	if (verbose == true) {
		cout << "Reciprocal Part: " << (t2.tv_sec - t1.tv_sec)
				+ (float) (t2.tv_nsec - t1.tv_nsec) / 1000000.0 << " ms\n";
	}
	//REAL SPACE
	shell_contr = 0.0;
	ewald_converged = false;
	for (int nact = 0; nact <= nmax; nact++) {
		if (!ewald_converged) {
			//#pragma omp parallel for private (tmp,tmpi,a,b,L,atom1,atom2,q1,q2,enq_sr,en,rsq,dist,lA,lrho) reduction(+:shell_contr,esum,esumq_sr,ngrad,grad1,grad2,grad3) schedule(static,CHUNKSIZE)
			for (int nx = -nact; nx <= nact; nx++) {
				////#pragma omp parallel for  reduction(+:esum,esumq_sr,esumq_e,ngrad) schedule(static,CHUNKSIZE)
				for (int ny = -nact; ny <= nact; ny++) {
					////#pragma omp parallel for  reduction(+:esum,esumq_sr,esumq_e,ngrad) schedule(static,CHUNKSIZE)
					for (int nz = -nact; nz <= nact; nz++) {
						//skip cells previously summed over
						if (abs(nx) == nact || abs(ny) == nact || abs(nz)
								== nact) {
							L[0] = nx * ar[0] + ny * br[0] + nz * cr[0];
							L[1] = nx * ar[1] + ny * br[1] + nz * cr[1];
							L[2] = nx * ar[2] + ny * br[2] + nz * cr[2];
							//cout << "L: <"<<L[0]<< " " << L[1] << " " <<L[2] <<">"
							//<< endl;
							//real space part
							for (int i = 0; i < numpos; i += 3) {
								a[0] = lxyz[i] + L[0];
								a[1] = lxyz[i + 1] + L[1];
								a[2] = lxyz[i + 2] + L[2];
								atom1 = i / 3;
								q1 = q[atom1];
								// inner loop
								////#pragma omp critical
								for (int j = 0; j < numpos; j += 3) {
									if (nx == 0 && ny == 0 && nz == 0) {
										//cout << "We are in unit cell..." << endl;
										if (i == j) {
											//cout << "skipping selfinteraction..." << endl;
											continue;
										}
									}

									//hier kommt der teuerste Teil des codes
									atom2 = j / 3;
									//small optimization, shift L to a!!!
									//get distances
									b[0] = a[0] - lxyz[j];
									b[1] = a[1] - lxyz[j + 1];
									b[2] = a[2] - lxyz[j + 2];
									rsq = pow(b[0], 2) + pow(b[1], 2) + pow(
											b[2], 2);
									dist = sqrt(rsq);
									//distance cut off

									if (dist > cutoff) {
										//cout << "Cut off too big...";
										continue;
									}
									//if this comes earliers codes slows down
									tmpi = getInteraction(atom_nr[atom1],
											atom_nr[atom2]);
									//############################################################
									//ENERGY
									//PP Part
									lA = paraMat[tmpi][0];
									lrho = paraMat[tmpi][1];
									q2 = q[atom2];
									en = 0.5 * lA * exp(-1 * dist / (ANG2BOHR
											* lrho));
									en = en * 1 / HARTREE2EV;
									esum = esum + en;
									//Coulomb Ewald Part
									//Evaluate SHORT RANGE part of Ewald, seems correct
									tmp = erfc_new(alpha * dist);
									//cout << tmp;
									enq_sr = 0.5 * q1 * q2 * tmp / dist;
									esumq_sr = esumq_sr + enq_sr;
									//Add Ewald components
									//esumq_e = esumq_sr + esumq_f + esumq_sc;
									//Evaluate conventional charge
									//enq = q1 * q2 / (dist);
									//esumq = esumq + enq;
									//#############################################################
									//Gradient Screening
									if (gcalc == true && abs((enq_sr + en)
											* HARTREE2EV) < 10E-6) {
										//printf ("ENQ_SR: %2.4E   EN: %2.4E\n", enq_sr, en);
										//										grad1 = 0;
										//										grad2 = 0;
										//										grad3 = 0;
										continue;
									}
									//GRADIENTS
									//#pragma omp critical
									if (gcalc == true) {
										ncgrad = q1 * q2
												* ((-2.0 * sqrt(1 / PI) * exp(
														-pow(alpha, 2) * rsq)
														* alpha * dist) - tmp)
												/ rsq;
										ngrad = lA * -1 / lrho * exp(-dist
												/ (ANG2BOHR * lrho));
										//ngrad = lA * -1 / lrho * exp(-dist / (lrho));
										ngrad = ngrad / (HARTREE2EV * ANG2BOHR)
												+ ncgrad;

										grad1 = grad1 + b[0] / dist * ngrad;
										grad2 = grad2 + b[1] / dist * ngrad;
										grad3 = grad3 + b[2] / dist * ngrad;
									}
								}//end loop j
								//#pragma omp critical
								if (gcalc == true) {
									lgrad[i] = lgrad[i] + grad1;
									lgrad[i + 1] = lgrad[i + 1] + grad2;
									lgrad[i + 2] = lgrad[i + 2] + grad3;
									grad1 = 0;
									grad2 = 0;
									grad3 = 0;
								}
							}//end loop i
						} //end if abs(nx)==nact || abs(ny)==nact || abs(nz)==nact
						//cout.precision(8);
						//if (nx < 1 && ny < 1 && nz < 1) {
						//					cout << "nx: " << nx << " ny: " << ny << " nz: " << nz
						//							<< " enq_sr:\t" << enq_sr << "en:\t" << en
						//							*HARTREE2EV<< endl;
						//					}
					}// nz for loop
				}// ny for loop
			}// nx for loop
			shell_contr = esumq_sr - shell_contr;
			if (verbose == true) {
				cout.precision(8);
				cout << "Shell #" << nact << " contribution, Q real part: "
						<< setw(12) << shell_contr * HARTREE2EV;
				cout << "\tPP part: " << esum * HARTREE2EV << endl;
			}
			if (abs(shell_contr * HARTREE2EV) < ewald_thresh) {
				ewald_converged = true;
				//#pragma omp flush (ewald_converged)
			}
			shell_contr = esumq_sr;
		}
	}// nact for loop

	energy = esum + esumq_sr + esumq_f + esumq_sc;
	//cout << "esum_cor : " << esum_cor*HARTREE2EV <<endl;
	//	cout << endl << " Energy - Conventional:\t\t" << esumq *HARTREE2EV << endl;


	//timing
	clock_gettime(CLOCK_REALTIME, &t3);
	if (verbose == true) {
		cout << "Real Part: " << (t3.tv_sec - t1.tv_sec) + (float) (t3.tv_nsec
				- t1.tv_nsec) / 1000000.0 << " ms" << endl;
	}
	if (verbose == true) {
		cout.precision(8);
		cout << endl << "Energy - pair potential\t\t:	" << setw(12) << esum
				* HARTREE2EV << " eV" << endl;
		cout << endl << "Energy - Coulomb\t\t:	" << setw(12) << (esumq_sr
				+ esumq_f + esumq_sc) * HARTREE2EV << " eV";
		cout << endl << "Ewald - real space\t\t:	" << setw(12) << esumq_sr
				* HARTREE2EV << " eV" << endl;
		cout << "Ewald - reciprocal space\t:	" << setw(12) << (esumq_f)
				* HARTREE2EV << " eV" << endl;
		cout << "Ewald - self correction\t\t:	" << setw(12) << esumq_sc
				* HARTREE2EV << " eV" << endl << endl;
		cout << "#Total energy\t\t\t:	" << setw(12) << energy * HARTREE2EV
				<< " eV" << endl << endl;
		//energy=esum+esumq_sr+esumq_f+esumq_sc+esum_dipole;
		cout << "Energy - dipole correction\t:	" << setw(12) << esum_dipole
				* HARTREE2EV << " eV" << endl;
		cout << "#Total energy (+dipole corr.)\t:	" << setw(12) << (energy
				+ esum_dipole) * HARTREE2EV << " eV" << endl;
	}
	if (gcalc == true) {
		//#pragma omp parallel for schedule(static,CHUNKSIZE)
		for (int i = 0; i < numpos; i++) {
			(*grad)(i) = lgrad[i];
			//which gradients are frozen??
			(*grad)(i) = lgrad[i] * (*moveMat)(i);
			//cout << "grad:"<< (*grad)(i)<<endl;
			//cout << lgrad[i]*27.211*1.88972598<<endl;
		}
		//gradnorm = Blas_Norm2((*grad));
		gradnorm = grad->norm();
		//gradmaxelem = Blas_Index_Max((*grad));
		grad->maxCoeff(&gradmaxelem);
		//gradinfnorm = Blas_Norm_Inf((*grad));
	}
	//cout << "Energy: "<<esum <<endl;

	//timing
	//cout << "Rec. space:\t";
	//newIO.printTiming(rec_start, rec_end);
	//cout << "Real space:\t";
	//newIO.printTiming(real_start, real_end);
	calculated = true;
}
//To Do: check method if each O has corresponding shell
//To Do: Add back the interaction of each shell with its core, only for inner cell-> real space!!!
//periodic version of shell-modell
//Shell-Modell needs to be redefined
//Warning: non-polarizable species is defined as shell!
//					PP      Coulomb     Spring
//shell-shell   	X			X          -
//core-shell(close) -           -          X
//core-core     	-           X          -
//core-shell(remote)-			X		   -
void ToolCalcPP::E_periodic_Pol(bool pol, bool verbose, bool gcalc) {
	//Only orthorombic cells at the moment
	//dipole correction not checked...
	bool ewald_converged = false;
	//image counters
	int nmax = 50;
	double lxyz[atnumber * 3], lgrad[atnumber * 3];
	double a[3], b[3], L[3], G[3];
	double dist = 0, rsq = 0, en = 0, enq_sr = 0, ecorr = 0, esum = 0,
			esumq_sr = 0, enq_f = 0, esumq_f = 0, enq_sc = 0, esumq_sc = 0,
			esum_dipole = 0, es = 0, sgrad = 0;
	double ngrad = 0, ncgrad = 0, gulp_fac = 0;
	double lA = 0., lrho = 1.;
	double q1 = 0., q2 = 0.;
	int atom1 = 0, atom2 = 0, tmpi = 0;
	double tmp = 0;
	//width of gaussian
	double epsilon = 1.0;
	//reciprocal space variables
	double ksq = 0, rhok_re = 0, rhok_im = 0, rhoki_re = 0, rhoki_im = 0, phr =
			0, phi = 0, grad1 = 0, grad2 = 0, grad3 = 0;
	//Setup
	if (verbose == true) {
		cout.precision(2);
		cout << endl << endl << "Periodic calculation with Lx=" << Lx
				<< ", Ly=" << Ly << ", Lz=" << Lz << " angstrom" << endl;
		cout.precision(8);
		cout << "Coulomb interaction calculated via Ewald approach, alpha="
				<< alpha << endl;
	}
	int numpos = 3 * (atnumber - nr_cores);
	int numposall = atnumber * 3;
	if (verbose == true && pol == true) {
		cout << "Shell-modell switched ON...\n";
	}
	if (verbose == true && pol == false) {
		cout << "Shell-modell switched OFF...\n";
	}
	if (!(nr_shells > 0)) {
		cout << "WARNING! No shells defined, wrong energy function used!\n";
		exit(1);
	}
	bool core1 = false, core2 = false, shel1 = false, shel2 = false;
	//cout << "numpos: " << numpos << "\n";
	//cout << "numposall: " << numposall << "\n";
	//setting up local arrays
	for (int i = 0; i < numposall; i++) {
		lxyz[i] = (*xyz)(i) * ANG2BOHR;
		lgrad[i] = 0;
	}

	//Dipole correction + self interaction
	double tmpQX = 0;
	double tmpQY = 0;
	double tmpQZ = 0;
	//Energy dipole correction
	if (dipole_cor == true) {
		setDipole(false);
		esum_dipole = 2*PI * dipnorm * dipnorm * ANG2BOHR * ANG2BOHR / ((2
				+ epsilon) * V);
		//SUM Q*R
		for (int j = 0; j < numpos; j += 3) {
			atom2 = j / 3;
			q2 = q[atom2];
			//atoms in cell, right dipole?
			a[0] = lxyz[j];
			a[1] = lxyz[j + 1];
			a[2] = lxyz[j + 2];
			tmpQX = tmpQX + q2 * a[0];
			tmpQY = tmpQY + q2 * a[1];
			tmpQZ = tmpQZ + q2 * a[2];
		}
	}
	//correction for self interaction
	for (int i = 0; i < numposall; i += 3) {
		atom1 = i / 3;
		q1 = q[atom1];
		//Make flexible self interaction!!!
		if (atom_nr[atom1] == 8 && core_shell[atom1] == "CORE") {
			continue;
		}
		if (atom_nr[atom1] == 8 && core_shell[atom1] == "SHEL") {
			continue;
		}
		enq_sc = pow(q1, 2);
		esumq_sc = esumq_sc + enq_sc;
		//gradient dipole correction
		//Not stable during optimization...
		//Are all atoms inside the cell?!!!!!!!!project them back...
		if (dipole_cor == true) {
			lgrad[i] = 4*PI * q1 * tmpQX / ((2 + epsilon) * V);
			lgrad[i + 1] = 4*PI * q1 * tmpQY / ((2 + epsilon) * V);
			lgrad[i + 2] = 4*PI * q1 * tmpQZ / ((2 + epsilon) * V);
		}
	}

	//RECIPROCAL SPACE - exploit symmetry...?
	double shell_contr = 0.0;
	//iteration of cell counter
	for (int nact = 1; nact <= nmax; nact++) {
		if (ewald_converged == false) {
			for (int nx = -nact; nx <= nact; nx++) {
				for (int ny = -nact; ny <= nact; ny++) {
					for (int nz = -nact; nz <= nact; nz++) {
						//skip central cell
						if (nx == 0 && ny == 0 && nz == 0) {
							enq_f = 0;
							continue;
						}
						//skip cells previously summed over
						if (abs(nx) == nact || abs(ny) == nact || abs(nz)
								== nact) {
							//Reciprocal space
							G[0] = 2* PI * nx * af[0] + 2* PI * ny * bf[0] + 2
									* PI * nz * cf[0];
							G[1] = 2* PI * nx * af[1] + 2* PI * ny * bf[1] + 2
									* PI * nz * cf[1];
							G[2] = 2* PI * nx * af[2] + 2* PI * ny * bf[2] + 2
									* PI * nz * cf[2];
							//k²
							ksq = G[0] * G[0] + G[1] * G[1] + G[2] * G[2];
							//cout << "ksq: " << ksq << endl;
							//real and imaginary part
							rhok_re = 0;
							rhok_im = 0;
							//ENERGY
							for (int i = 0; i < numposall; i += 3) {
								atom1 = i / 3;
								if (core_shell[atom1] == "CORE" && pol == false) {
									//skip core atoms if no shell model is wanted
									continue;
								}
								a[0] = lxyz[i];
								a[1] = lxyz[i + 1];
								a[2] = lxyz[i + 2];
								q1 = q[atom1];
								if (atom_nr[atom1] == 8 && pol == false) {
									q1 = -2.0;
								}
								//calculate k * ri
								tmp = G[0] * a[0] + G[1] * a[1] + G[2] * a[2];
								rhoki_re = q1 * cos(tmp);
								rhoki_im = q1 * sin(tmp);
								//calculate rhok
								rhok_re = rhok_re + rhoki_re;
								rhok_im = rhok_im + rhoki_im;
							}
							//GRADIENT, second summation neeed
							if (gcalc == true) {
								for (int i = 0; i < numposall; i += 3) {
									atom1 = i / 3;
									if (core_shell[atom1] == "CORE" && pol
											== false) {
										//skip core atoms if no shell model is wanted
										continue;
									}
									phr = 0;
									phi = 0;
									//only phase important
									a[0] = lxyz[i];
									a[1] = lxyz[i + 1];
									a[2] = lxyz[i + 2];
									q1 = q[atom1];
									if (atom_nr[atom1] == 8 && pol == false) {
										q1 = -2.0;
									}
									tmp = G[0] * a[0] + G[1] * a[1] + G[2]
											* a[2];
									phr = cos(tmp);
									phi = -sin(tmp);
									tmp = (phi * rhok_re + phr * rhok_im) * q1
											* 4*PI / (V * ksq) * exp(-1 * ksq
											/ (4 * pow(alpha, 2)));
									//cout << "tmp:" << tmp << endl;
									lgrad[i] = lgrad[i] + tmp * G[0];
									lgrad[i + 1] = lgrad[i + 1] + tmp * G[1];
									lgrad[i + 2] = lgrad[i + 2] + tmp * G[2];

								}
							}
							tmp = 2* PI * (rhok_re * rhok_re + rhok_im
									* rhok_im) / (V * ksq);
							enq_f = tmp * exp(-1* ksq / (4* pow (alpha, 2)));
							esumq_f = esumq_f + enq_f;
						}//end if
						//cout.precision(8);
						//cout << "nx: " << nx << " ny: " << ny << " nz: " << nz << " enq_f:\t" << enq_f << endl;
					} // nz for loop
				} // ny for loop
			}// nx for loop
			shell_contr = esumq_f - shell_contr;
			if (verbose == true) {
				cout << "Shell #" << nact << " contribution, Q recip. part: "
						<< shell_contr * HARTREE2EV << endl;
			}
			if (abs(shell_contr * HARTREE2EV < ewald_thresh) && nact > 1) {
				//cout << "We could leave now...";
				ewald_converged = true;
			}
		}//end if ewald_converged
		shell_contr = esumq_f;
	}// nact for loop

	//REAL SPACE
	shell_contr = 0.0;
	ewald_converged = false;
	for (int nact = 0; nact <= nmax; nact++) {
		if (!ewald_converged) {
			for (int nx = -nact; nx <= nact; nx++) {
				for (int ny = -nact; ny <= nact; ny++) {
					for (int nz = -nact; nz <= nact; nz++) {
						//skip cells previously summed over
						if (abs(nx) == nact || abs(ny) == nact || abs(nz)
								== nact) {
							L[0] = nx * ar[0] + ny * br[0] + nz * cr[0];
							L[1] = nx * ar[1] + ny * br[1] + nz * cr[1];
							L[2] = nx * ar[2] + ny * br[2] + nz * cr[2];
							//cout << "L: <"<<L[0]<< " " << L[1] << " " <<L[2] <<">"
							//<< endl;
							//real space part
							for (int i = 0; i < numposall; i += 3) {
								core1 = false;
								shel1 = false;
								atom1 = i / 3;
								if (core_shell[atom1] == "CORE" && pol == false) {
									//skip core atoms if no shell model is wanted
									continue;
								} else {
									if (core_shell[atom1] == "CORE") {
										//cout << "O_Ci"<<endl;
										core1 = true;
									}
									if (core_shell[atom1] == "SHEL") {
										//cout << "O_Ci"<<endl;
										shel1 = true;
									}
								}
								a[0] = lxyz[i] + L[0];
								a[1] = lxyz[i + 1] + L[1];
								a[2] = lxyz[i + 2] + L[2];
								q1 = q[atom1];
								if (atom_nr[atom1] == 8 && pol == false) {
									q1 = -2.0;
								}
								// inner loop
								for (int j = 0; j < numposall; j += 3) {
									if (nx == 0 && ny == 0 && nz == 0) {
										//cout << "We are in unit cell..." << endl;
										if (i == j) {
											//skipping selfinteraction...
											continue;
										}
									}
									core2 = false;
									shel2 = false;
									atom2 = j / 3;
									if (core_shell[atom2] == "CORE" && pol
											== false) {
										//skip core atoms if no shell model is wanted
										continue;
									} else {
										if (core_shell[atom2] == "CORE") {
											//cout << "CORE"<<endl;
											core2 = true;
										}
										if (core_shell[atom2] == "SHEL") {
											//cout << "SHELL"<<endl;
											shel2 = true;
										}
									}
									//small optimization, shift L to a!!!
									//get distances
									b[0] = a[0] - lxyz[j];
									b[1] = a[1] - lxyz[j + 1];
									b[2] = a[2] - lxyz[j + 2];
									rsq = pow(b[0], 2) + pow(b[1], 2) + pow(
											b[2], 2);
									dist = sqrt(rsq);
									//distance cut off, avoiding also Coulomb explosion??
									if (dist > cutoff) {
										//if (dist > cutoff) {
										//cout << "Cut off too big...";
										continue;
									}
									//if this comes earlier codes slows down
									//perhaps paraMat can be decreased without shells...
									tmpi = getInteraction(atom_nr[atom1],
											atom_nr[atom2]);
									//############################################################
									//SPRING POTENTIAL
									if (((core1 == true && shel2 == true)
											|| (shel1 == true && core2 == true))
											&& (dist < scutoff * ANG2BOHR)) {
										//energy
										es = 0.25 * springc * pow(dist
												/ ANG2BOHR, 2);
										es = es / HARTREE2EV;
										//gradient * empirical corr
										gulp_fac = 1.0584;
										sgrad = 0.5 * springc * dist
												/ (ANG2BOHR * HARTREE2EV)
												* gulp_fac;
										//ngrad = 0;
										ncgrad = 0;
									}
									//PP Part
									if (shel1 == true && shel2 == true) {
										lA = paraMat[tmpi][0];
										lrho = paraMat[tmpi][1];

										en = 0.5 * lA * exp(-1 * dist
												/ (ANG2BOHR * lrho));
										en = en * 1 / HARTREE2EV;
										ngrad = lA * -1 / lrho * exp(-dist
												/ (ANG2BOHR * lrho));
										sgrad = 0;
									}
									esum = esum + en + es;
									//Coulomb Ewald Part
									q2 = q[atom2];
									if (atom_nr[atom2] == 8 && pol == false) {
										q2 = -2.0;
									}
									//Evaluate SHORT RANGE part of Ewald except for nearest C-S pairs
									if (((core1 == true && shel2 == true)
											|| (core2 == true && shel1 == true))
											&& (dist < scutoff)) {
										//Treat core-shell pair with distance=0 separately
										//Flexible self interaction
										enq_sc = 0.5 * pow(2.0, 2);
										ncgrad = 0;
										enq_sr = 0;
										ecorr = 0;

										if (dist > 0.0001) {
											ecorr = -0.5 * q1 * q2 / (dist);
											tmp = erfc_new(alpha * dist);
											enq_sr = 0.5 * q1 * q2 * tmp / dist;
											//											cout << "Coulomb-Interaction for C-S pair "
											//													<< atom1 + 1 << " and "
											//													<< atom2 + 1 << ": "
											//													<< ecorr * HARTREE2EV
											//													<< "& Q_sr: " << enq_sr
											//													* HARTREE2EV << endl;
											//Flexible self interaction
											enq_sc = 0.5 * (pow(q1, 2) + pow(
													q2, 2));
											//Gradient
											//											grad_corr =   q1 * q2 / pow(dist, 2);
											//											ncgrad = q1 * q2 * ((-2.0 * sqrt(1 / PI) * exp(	-pow(alpha, 2) * rsq) * alpha * dist) - tmp)
											//																						/ rsq;
											//											cout << "grad_corr: "<< grad_corr << " ncgrad: " << ncgrad << endl;
											//											ncgrad = ncgrad + grad_corr;
										}
										//Flexible self interaction
										esumq_sc = esumq_sc + enq_sc;
									} else {
										tmp = erfc_new(alpha * dist);
										//cout << tmp;
										enq_sr = 0.5 * q1 * q2 * tmp / dist;
										ncgrad = q1 * q2
												* ((-2.0 * sqrt(1 / PI) * exp(
														-pow(alpha, 2) * rsq)
														* alpha * dist) - tmp)
												/ rsq;
									}
									//cout << ncgrad << endl;
									esumq_sr = esumq_sr + enq_sr + ecorr;
									ecorr = 0;
									//#############################################################
									//GRADIENTS
									if (gcalc == true) {
										ngrad = ngrad / (HARTREE2EV * ANG2BOHR);
										//cout << "ngrad:" << ngrad;
										//cout << " sgrad: " << sgrad << endl;
										ngrad = ngrad + ncgrad + sgrad;
										if (dist > 0) {
											grad1 = grad1 + b[0] / dist * ngrad;
											grad2 = grad2 + b[1] / dist * ngrad;
											grad3 = grad3 + b[2] / dist * ngrad;

											//										cout << " sgradx: " <<b[0] / dist * sgrad *HARTREE2EV * ANG2BOHR*Lx ;
											//										cout << " sgrady: " <<b[1] / dist * sgrad *HARTREE2EV * ANG2BOHR*Ly ;
											//										cout << " sgradz " <<b[2] / dist * sgrad *HARTREE2EV * ANG2BOHR*Lz << endl;
										}
									}
									en = 0;
									es = 0;
									ngrad = 0;
									sgrad = 0;
								}//end loop j
								if (gcalc == true) {
									lgrad[i] = lgrad[i] + grad1;
									lgrad[i + 1] = lgrad[i + 1] + grad2;
									lgrad[i + 2] = lgrad[i + 2] + grad3;
									grad1 = 0;
									grad2 = 0;
									grad3 = 0;
								}
							}//end loop i
						} //end if abs(nx)==nact || abs(ny)==nact || abs(nz)==nact
						//cout.precision(8);
						//if (nx < 1 && ny < 1 && nz < 1) {
						//					cout << "nx: " << nx << " ny: " << ny << " nz: " << nz
						//							<< " enq_sr:\t" << enq_sr << "en:\t" << en
						//							*HARTREE2EV<< endl;
						//					}
					}// nz for loop
				}// ny for loop
			}// nx for loop
			shell_contr = esumq_sr - shell_contr;
			if (verbose == true) {
				cout.precision(8);
				cout << "Shell #" << nact << " contribution, Q real part: "
						<< setw(12) << shell_contr * HARTREE2EV;
				cout << "\tPP part: " << esum * HARTREE2EV << endl;
			}
			if ((abs(shell_contr * HARTREE2EV) < ewald_thresh) && nact > 0) {
				ewald_converged = true;
			}
			shell_contr = esumq_sr;
		}//end if ewald converged
	}// nact for loop

	//correction for self interaction
	esumq_sc = -1 * alpha / sqrt(PI) * esumq_sc;

	energy = esum + esumq_sr + esumq_f + esumq_sc;
	//cout << "esum_cor : " << esum_cor*HARTREE2EV <<endl;
	//	cout << endl << " Energy - Conventional:\t\t" << esumq *HARTREE2EV << endl;
	if (verbose == true) {
		cout.precision(8);
		cout << endl << "Energy - pair potential\t\t:	" << setw(12) << esum
				* HARTREE2EV << " eV" << endl;
		cout << endl << "Energy - Coulomb\t\t:	" << setw(12) << (esumq_sr
				+ esumq_f + esumq_sc) * HARTREE2EV << " eV";
		cout << endl << "Ewald - real space\t\t:	" << setw(12) << esumq_sr
				* HARTREE2EV << " eV" << endl;
		cout << "Ewald - reciprocal space\t:	" << setw(12) << (esumq_f)
				* HARTREE2EV << " eV" << endl;
		cout << "Ewald - self correction\t\t:	" << setw(12) << esumq_sc
				* HARTREE2EV << " eV" << endl << endl;
		cout << "#Total energy\t\t\t:	" << setw(12) << energy * HARTREE2EV
				<< " eV" << endl << endl;
		//energy=esum+esumq_sr+esumq_f+esumq_sc+esum_dipole;
		cout << "Energy - dipole correction\t:	" << setw(12) << esum_dipole
				* HARTREE2EV << " eV" << endl;
		cout << "#Total energy (+dipole corr.)\t:	" << setw(12) << (energy
				+ esum_dipole) * HARTREE2EV << " eV" << endl;
	}
	if (gcalc == true) {
		for (int i = 0; i < numposall; i++) {
			(*grad)(i) = lgrad[i];
			//which gradients are frozen??
			(*grad)(i) = lgrad[i] * (*moveMat)(i);
			//cout << "grad:"<< (*grad)(i)<<endl;
			//			cout << lgrad[i]* HARTREE2EV
			//			* ANG2BOHR<<endl;
		}
		//gradnorm = Blas_Norm2((*grad));
		gradnorm = grad->norm();
		//gradmaxelem = Blas_Index_Max((*grad));
		grad->maxCoeff(&gradmaxelem);
		//gradinfnorm = Blas_Norm_Inf((*grad));
	}
	//cout << "Energy: "<<esum <<endl;
	calculated = true;
}

//emp. pair potential under periodic boundary conditions
//TO DO: Symmetry, linked list, FFT, combine loops, parinline ...
//void ToolCalcPP::E_periodic(bool pol, bool verbose, bool gcalc) {
//
//	//image counters
//	int nmax=50;
//	//Only orthorombic cells at the moment
//	//Optimisations: prepare cell paramters in separate function
//	double lxyz[atnumber*3];
//	double lgrad[atnumber*3];
//	double a[3];
//	double b[3];
//	double dist=0, rsq=0;
//	double en=0, enq_sr=0, esum=0, esumq_sr=0, enq_f=0, esumq_f=0, enq_sc=0,
//			esumq_sc=0, esumq_e=0, esum_dipole=0;
//	double ngrad=0, ncgrad=0;
//	double lA =0., lrho=1.;
//	double q1=0., q2=0.;
//	int atom1=0, atom2=0, tmpi=0;
//	//	//orthogonal box
//	double L[3];
//	//width of gaussian
//	double epsilon = 1;
//	//reciprocal space variables
//	int nx=0;
//	int ny=0;
//	int nz=0;
//	double tmp=0;
//	double ksq=0;
//	double rhok_re=0;
//	double rhok_im=0;
//	double rhoki_re=0;
//	double rhoki_im=0;
//	double phr=0;
//	double phi=0;
//	double grad1=0, grad2=0, grad3=0;
//	//lattice vectors in real space
//	//double ar[3], br[3], cr[3];
//	//lattice vector in Fourier space
//	//double af[3], bf[3], cf[3];
//	double G[3];
//	if (verbose == true) {
//		cout.precision(2);
//		cout << endl <<endl << "Periodic calculation with Lx="<<Lx<< ", Ly="
//				<<Ly <<", Lz="<<Lz << " angstrom"<< endl;
//		cout.precision(8);
//		cout << "Coulomb interaction calculated via Ewald approach, alpha="
//				<< alpha <<endl << endl;
//	}
//	int numpos=atnumber*3;
//
//	//setting up local arrays
//	for (int i=0; i<numpos; i++) {
//		lxyz[i]=(*xyz)(i)*ANG2BOHR;
//		lgrad[i]=0;
//	}
//	//Gradient dipole correction + self interaction
//	double tmpQ[3]= { 0, 0, 0 };
//	if (dipole_cor == true) {
//		//Energy dipole correction
//		setDipole(false);
//		esum_dipole = 2*PI * dipnorm*dipnorm * ANG2BOHR *ANG2BOHR/((2+epsilon)
//				*V);
//		//SUM Q*R
//		for (int j=0; j<numpos; j += 3) {
//			atom2=j/3;
//			q2=q[atom2];
//			//atoms in cell, right dipole?
//			a[0]=lxyz[j];
//			a[1]=lxyz[j+1];
//			a[2]=lxyz[j+2];
//			tmpQ[0]=tmpQ[0]+q2*a[0];
//			tmpQ[1]=tmpQ[1]+q2*a[1];
//			tmpQ[2]=tmpQ[2]+q2*a[2];
//		}
//	}
//
//	for (int i=0; i<numpos; i +=3) {
//		//correction for self interaction
//		atom1=i/3;
//		q1=q[atom1];
//		enq_sc = pow(q1, 2);
//		esumq_sc = esumq_sc +enq_sc;
//		//gradient dipole correction
//		//cout << tmpQ[0] << endl;
//		//Not stable during optimization
//		//Are all atoms inside the cell?!!!!!!!!project them back...
//		if (dipole_cor == true) {
//			lgrad[i] =4*PI*q1*tmpQ[0]/((2+epsilon)*V);
//			lgrad[i + 1] =4*PI*q1*tmpQ[1]/((2+epsilon)*V);
//			lgrad[i + 2] =4*PI*q1*tmpQ[2]/((2+epsilon)*V);
//		}
//	}
//
//	//correction for self interaction
//	esumq_sc = -1 * alpha / sqrt(PI) * esumq_sc;
//
//	//RECIPROCAL SPACE - Symmetry!!!!?
//
//
//	double shell_contr = 0.0;
//	//iteration of cell counter
//	for (int nact=1; nact<=nmax; nact++) {
//		for (nx=-nact; nx<=nact; nx++) {
//			//id = omp_get_thread_num();
//			//printf("Hello World from thread %d\n", id);
//			for (ny=-nact; ny<=nact; ny++) {
//				for (nz=-nact; nz<=nact; nz++) {
//					//skip central cell
//					if (nx==0 && ny==0 && nz==0) {
//						enq_f=0;
//						continue;
//					}
//					//skip cells previously summed over
//					if (abs(nx)==nact || abs(ny)==nact || abs(nz)==nact) {
//						//Reciprocal space
//						G[0] = 2* PI * nx * af[0] + 2* PI * ny * bf[0] + 2* PI
//								* nz * cf[0];
//						G[1] = 2* PI * nx * af[1] + 2* PI * ny * bf[1] + 2* PI
//								* nz * cf[1];
//						G[2] = 2* PI * nx * af[2] + 2* PI * ny * bf[2] + 2* PI
//								* nz * cf[2];
//						//k²
//						ksq= G[0]*G[0]+G[1]*G[1]+G[2]*G[2];
//						//real and imaginary part
//						rhok_re=0;
//						rhok_im=0;
//						//ENERGY
//						for (int i=0; i<numpos; i += 3) {
//							a[0]=lxyz[i];
//							a[1]=lxyz[i+1];
//							a[2]=lxyz[i+2];
//							atom1=i/3;
//							q1=q[atom1];
//							//calculate k * ri
//							tmp = G[0]*a[0]+G[1]*a[1]+G[2]*a[2];
//							rhoki_re = q1 * cos(tmp);
//							rhoki_im = q1 * sin(tmp);
//							//calculate rhok
//							rhok_re = rhok_re + rhoki_re;
//							rhok_im = rhok_im + rhoki_im;
//						}
//						//GRADIENT, second summation neeed
//						if (gcalc == true) {
//							for (int i=0; i<numpos; i += 3) {
//								phr=0;
//								phi=0;
//								//only phase important
//								a[0]=lxyz[i];
//								a[1]=lxyz[i+1];
//								a[2]=lxyz[i+2];
//								atom1=i/3;
//								q1=q[atom1];
//								tmp = G[0]*a[0]+G[1]*a[1]+G[2]*a[2];
//								phr=cos(tmp);
//								phi=-sin(tmp);
//								tmp=(phi*rhok_re+phr*rhok_im)*q1*4*PI/(V*ksq)
//										*exp(-1 *ksq/(4*pow(alpha, 2)));
//								//cout << "tmp:" << tmp << endl;
//								lgrad[i] = lgrad[i] + tmp*G[0];
//								lgrad[i + 1] = lgrad[i + 1] +tmp*G[1];
//								lgrad[i + 2] = lgrad[i + 2] + tmp*G[2];
//							}
//						}
//						tmp = 2* PI * (rhok_re*rhok_re+rhok_im*rhok_im) / (V
//								* ksq);
//						enq_f = tmp * exp(-1*ksq/(4*pow(alpha, 2)));
//						//enq_f = 2* tmp * exp(-1*ksq/(4*alpha));
//						esumq_f=esumq_f + enq_f;
//					}//end if
//					//cout.precision(8);
//					//cout << "nx: " << nx << " ny: " << ny << " nz: " << nz << " enq_f:\t" << enq_f << endl;
//				} // nz for loop
//			} // ny for loop
//		}// nx for loop
//		shell_contr = esumq_f - shell_contr;
//		if (verbose == true) {
//			cout << "Shell #"<<nact<<" contribution, Q recip. part: "
//					<< shell_contr * HARTREE2EV << endl;
//		}
//		if (shell_contr * HARTREE2EV < ewald_thresh) {
//			//cout << "We could leave now...";
//			break;
//		}
//		shell_contr = esumq_f;
//	}// nact for loop
//	//
//	//Alternative reciproc. gradients, Deserno et al.
//	/*double tmpG[3]= { 0, 0, 0 };
//	 double tmpQ[3]= { 0, 0, 0 };
//	 for (int i=0; i<numpos; i += 3) {
//	 atom1=i/3;
//	 //Coulomb Part
//	 q1=q[atom1];
//	 a[0]=lxyz[i]*ANG2BOHR;
//	 a[1]=lxyz[i+1]*ANG2BOHR;
//	 a[2]=lxyz[i+2]*ANG2BOHR;
//	 tmpQ[0]=0;
//	 tmpQ[1]=0;
//	 tmpQ[2]=0;
//	 for (int j=0; j<numpos; j += 3) {
//	 if (j==i) {
//	 continue;
//	 }
//	 atom2=j/3;
//	 q2=q[atom2];
//	 b[0]=lxyz[j]*ANG2BOHR;
//	 b[1]=lxyz[j+1]*ANG2BOHR;
//	 b[2]=lxyz[j+2]*ANG2BOHR;
//	 //get distances
//	 b[0]=a[0]-b[0];
//	 b[1]=a[1]-b[1];
//	 b[2]=a[2]-b[2];
//	 //get G
//	 tmpG[0]=0;
//	 tmpG[1]=0;
//	 tmpG[2]=0;
//	 for (int nact=1; nact<=nmax; nact++) {
//	 for (nx=-nact; nx<=nact; nx++) {
//	 for (ny=-nact; ny<=nact; ny++) {
//	 for (nz=-nact; nz<=nact; nz++) {
//	 //skip central cell
//	 if (nx==0 && ny==0 && nz==0) {
//	 enq_f=0;
//	 continue;
//	 }
//	 if (abs(nx)==nact || abs(ny)==nact || abs(nz)==nact) {
//	 //Reciprocal space
//	 G[0] = 2* PI * nx * af[0] + 2* PI * ny * bf[0]
//	 + 2* PI * nz * cf[0];
//	 G[1] = 2* PI * nx * af[1] + 2* PI * ny * bf[1]
//	 + 2* PI * nz * cf[1];
//	 G[2] = 2* PI * nx * af[2] + 2* PI * ny * bf[2]
//	 + 2* PI * nz * cf[2];
//	 //k²
//	 ksq= G[0]*G[0]+G[1]*G[1]+G[2]*G[2];
//	 tmp = G[0]*b[0]+G[1]*b[1]+G[2]*b[2];
//	 //cout << "k*r: " << tmp << endl;
//	 tmp=sin(tmp)*4*PI/ksq *exp(-ksq/(4 *pow(alpha,
//	 2)));
//	 //cout << "grad:" << tmp << endl;
//	 //tmp=0;
//	 tmpG[0] = tmpG[0] + tmp*G[0];
//	 tmpG[1] = tmpG[1] + tmp*G[1];
//	 tmpG[2] = tmpG[2] + tmp*G[2];
//	 }
//	 } // nz for loop
//	 }// ny for loop
//	 }// nx for loop
//	 }// nact loop
//	 tmpQ[0]=tmpQ[0]+q2*tmpG[0];
//	 tmpQ[1]=tmpQ[1]+q2*tmpG[1];
//	 tmpQ[2]=tmpQ[2]+q2*tmpG[2];
//	 }//end qj loop
//	 lgrad[i] = -tmpQ[0]*q1/V;
//	 lgrad[i + 1] = -tmpQ[1]*q1/V;
//	 lgrad[i + 2] = -tmpQ[2]*q1/V;
//	 //
//	 //		cout << "lgrad[i]:" << lgrad[i]*Lx*HARTREE2EV*ANG2BOHR;
//	 //		cout << " lgrad[i+1]:" << lgrad[i+1]*Ly*HARTREE2EV*ANG2BOHR;
//	 //		cout << " lgrad[i+2]:" << lgrad[i+2]*Lz*HARTREE2EV*ANG2BOHR << endl;
//
//	 }//end qi loop
//	 */
//	//REAL SPACE
//	shell_contr = 0.0;
//	for (int nact=0; nact<=nmax; nact++) {
//		for (nx=-nact; nx<=nact; nx++) {
//			for (ny=-nact; ny<=nact; ny++) {
//				for (nz=-nact; nz<=nact; nz++) {
//					//skip cells previously summed over
//					if (abs(nx)==nact || abs(ny)==nact || abs(nz)==nact) {
//						L[0]= nx * ar[0] + ny * br[0] + nz * cr[0];
//						L[1]= nx * ar[1] + ny * br[1] + nz * cr[1];
//						L[2]= nx * ar[2] + ny * br[2] + nz * cr[2];
//						//cout << "L: <"<<L[0]<< " " << L[1] << " " <<L[2] <<">"
//						//<< endl;
//						//real space part
//						for (int i=0; i<numpos; i += 3) {
//							a[0]=lxyz[i]+L[0];
//							a[1]=lxyz[i+1]+L[1];
//							a[2]=lxyz[i+2]+L[2];
//							//							a[0]=lxyz[i];
//							//							a[1]=lxyz[i+1];
//							//							a[2]=lxyz[i+2];
//							atom1=i/3;
//							q1=q[atom1];
//							// inner loop
//							for (int j=0; j<numpos; j+=3) {
//								if (nx==0 && ny==0 && nz==0) {
//									//cout << "We are in unit cell..." << endl;
//									if (i==j) {
//										//cout << "skipping selfinteraction..." << endl;
//										continue;
//									}
//								}
//								//hier kommt der teuerste Teil des codes
//								atom2=j/3;
//								//small optimization, shift L to a!!!
//								//get coordinates
//								//	b[0]=lxyz[j]+L[0];
//								//	b[1]=lxyz[j+1]+L[1];
//								//	b[2]=lxyz[j+2]+L[2];
//								//get distances
//								b[0]=a[0]-lxyz[j];
//								b[1]=a[1]-lxyz[j+1];
//								b[2]=a[2]-lxyz[j+2];
//								rsq=pow(b[0], 2)+pow(b[1], 2)+pow(b[2], 2);
//								dist = sqrt(rsq);
//								//distance cut off
//								if (dist > cutoff) {
//									//cout << "Cut off too big...";
//									continue;
//								}
//								//if this comes earliers codes slows down
//								tmpi=getInteraction(atom_nr[atom1],
//										atom_nr[atom2]);
//								//############################################################
//								//ENERGY
//								//PP Part
//								lA = paraMat[tmpi][0];
//								lrho = paraMat[tmpi][1];
//								q2=q[atom2];
//								en = 0.5 * lA* exp(-1 * dist / (ANG2BOHR*lrho));
//								en = en *1 / HARTREE2EV;
//								esum=esum + en;
//								//Coulomb Ewald Part
//								//Evaluate SHORT RANGE part of Ewald, seems correct
//								tmp=erfc_new(alpha*dist);
//								//cout << tmp;
//								enq_sr = 0.5 * q1 * q2 * tmp / dist;
//								esumq_sr = esumq_sr + enq_sr;
//								//Add Ewald components
//								esumq_e = esumq_sr + esumq_f + esumq_sc;
//								//Evaluate conventional charge
//								//enq = q1 * q2 / (dist);
//								//esumq = esumq + enq;
//								//#############################################################
//								//GRADIENTS
//								if (gcalc == true) {
//									ncgrad = q1 *q2 * ( (-2.0 * sqrt(1/PI)
//											*exp(-pow(alpha, 2) * rsq)*alpha
//											*dist) -tmp)/rsq;
//									ngrad = lA * -1 / lrho * exp(-dist
//											/ (ANG2BOHR *lrho));
//									//ngrad = lA * -1 / lrho * exp(-dist / (lrho));
//									ngrad = ngrad / (HARTREE2EV * ANG2BOHR);
//									ngrad = ngrad + ncgrad;
//									grad1=grad1 + b[0] /dist * ngrad;
//									grad2=grad2 + b[1] /dist * ngrad;
//									grad3=grad3 + b[2] /dist * ngrad;
//								}
//							}//end loop j
//							if (gcalc==true) {
//								lgrad[i] = lgrad[i] + grad1;
//								lgrad[i + 1] = lgrad[i + 1] + grad2;
//								lgrad[i + 2] = lgrad[i + 2] + grad3;
//								grad1=0;
//								grad2=0;
//								grad3=0;
//							}
//						}//end loop i
//					} //end if abs(nx)==nact || abs(ny)==nact || abs(nz)==nact
//					//cout.precision(8);
//					//if (nx < 1 && ny < 1 && nz < 1) {
//					//					cout << "nx: " << nx << " ny: " << ny << " nz: " << nz
//					//							<< " enq_sr:\t" << enq_sr << "en:\t" << en
//					//							*HARTREE2EV<< endl;
//					//					}
//				}// nz for loop
//			}// ny for loop
//		}// nx for loop
//		shell_contr = esumq_sr - shell_contr;
//		if (verbose == true) {
//			cout.precision(8);
//			cout << "Shell #"<<nact<< " contribution, Q real part: "
//					<< setw(12) << shell_contr * HARTREE2EV;
//			cout << "\tPP part: " << esum* HARTREE2EV << endl;
//		}
//		if (abs(shell_contr * HARTREE2EV) < ewald_thresh) {
//			break;
//		}
//		shell_contr = esumq_sr;
//	}// nact for loop
//
//	energy=esum+esumq_sr+esumq_f+esumq_sc;
//	//cout << "esum_cor : " << esum_cor*HARTREE2EV <<endl;
//	//	cout << endl << " Energy - Conventional:\t\t" << esumq *HARTREE2EV << endl;
//	if (verbose == true) {
//		cout.precision(8);
//		cout << endl << "Energy - pair potential\t\t:	" << setw(12) << esum
//				*HARTREE2EV << " eV"<<endl;
//		cout << endl << "Energy - Coulomb\t\t:	" << setw(12) << esumq_e
//				*HARTREE2EV << " eV";
//		cout << endl << "Ewald - real space\t\t:	" << setw(12) << esumq_sr
//				*HARTREE2EV << " eV" <<endl;
//		cout << "Ewald - reciprocal space\t:	"<< setw(12) << (esumq_f)
//				*HARTREE2EV << " eV" << endl;
//		cout << "Ewald - self correction\t\t:	" << setw(12)<< esumq_sc
//				*HARTREE2EV << " eV"<<endl<<endl;
//		cout << "#Total energy\t\t\t:	"<< setw(12) << energy *HARTREE2EV
//				<< " eV" << endl<< endl;
//		//energy=esum+esumq_sr+esumq_f+esumq_sc+esum_dipole;
//		cout << "Energy - dipole correction\t:	" << setw(12)<< esum_dipole
//				*HARTREE2EV<< " eV"<<endl;
//		cout << "#Total energy (+dipole corr.)\t:	"<< setw(12) << (energy
//				+ esum_dipole) *HARTREE2EV<< " eV"<< endl;
//	}
//	if (gcalc == true) {
//		for (int i=0; i<numpos; i++) {
//			(*grad)(i)=lgrad[i];
//			//which gradients are frozen??
//			(*grad)(i)=lgrad[i]*(*moveMat)(i);
//			//cout << "grad:"<< (*grad)(i)<<endl;
//			//cout << lgrad[i]*27.211*1.88972598<<endl;
//		}
//		gradnorm=Blas_Norm2((*grad));
//		gradmaxelem=Blas_Index_Max((*grad));
//		gradinfnorm=Blas_Norm_Inf((*grad));
//	}
//	//cout << "Energy: "<<esum <<endl;
//	calculated=true;
//}


ToolCalcPP::~ToolCalcPP() {
}

