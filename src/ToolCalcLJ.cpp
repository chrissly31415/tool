#include "ToolCalcLJ.h"

ToolCalcLJ::ToolCalcLJ() {
	step =0.001; //LJ cluster
	cutoff = 40.;
	max_displacement = 0.1;
	threshhold = 0.1;
	maxiter = 500;
	global_maxiter = 20;
	nr_reset = 2;
	temperature=1;
	k=1;//for LJ
	dipole_cor=false;
}

//Energy calculation
void ToolCalcLJ::E(bool pol, bool verbose, bool gcalc) {
	double lxyz[atnumber*3];
	double lgrad[atnumber*3];
	double a[3];
	double b[3];
	double dist=0;
	double en=0;
	double esum=0;
	double ngrad=0;
	double sigma =1;
	double epsilon=1;
	//energyFunc EF1;
	int numpos=atnumber*3;
	//cout <<"nupos:" << numpos;
	for (int i=0; i<numpos; i++) {
		lxyz[i]=(*xyz)(i);
		//initialisation of grad!!
		lgrad[i]=0;
	}
	for (int i=0; i<numpos; i += 3) {
		a[0]=lxyz[i];
		a[1]=lxyz[i+1];
		a[2]=lxyz[i+2];
		for (int j=i+3; j<numpos; j+=3) {
			b[0]=lxyz[j];
			b[1]=lxyz[j+1];
			b[2]=lxyz[j+2];
			b[0]=a[0]-b[0];
			b[1]=a[1]-b[1];
			b[2]=a[2]-b[2];
			dist=pow(b[0], 2)+pow(b[1], 2)+pow(b[2], 2);
			dist = sqrt(dist);
			if (dist > cutoff) {
				//cout << "Cut off too big...";
				continue;
			}
			//en=energyFunc::enLJ(dist, epsilon, sigma);
			en = 4 * epsilon * (pow(sigma / dist, 12) - pow(sigma / dist, 6));
			esum=esum+en;
			ngrad=24* epsilon* ((-2 * pow(dist / sigma, -13) / sigma) + (pow(
					dist / sigma, -7) / sigma));
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
	energy=esum;
	if (pol == true) {
		cout << "energy:"<< energy << endl;
	}
	//grad = new LaVectorDouble(lgrad,numpos);
	for (int i=0; i<numpos; i++) {
		(*grad)(i)=lgrad[i];
		//which gradients are frozen??
		(*grad)(i)=lgrad[i]*(*moveMat)(i);
		//cout << "grad:"<< (*grad)(i)<<endl;
		//cout << lgrad[i]<<endl;
	}
	gradnorm=grad->norm();

	//gradmaxelem=Blas_Index_Max((*grad));
	grad->maxCoeff(&gradmaxelem);
	//gradinfnorm=Blas_Norm_Inf((*grad));
	//cout << "Energy: "<<esum <<endl;
	calculated=true;
}
//Energy calculation
void ToolCalcLJ::E_periodic(bool pol, bool verbose, bool gcalc) {
	//TODO include identical particle in other cells
	//change order of interactions ie image wise, perhaps two summation one for center cell and the rest for images 
	double lxyz[atnumber*3];
	double lgrad[atnumber*3];
	double a[3];
	double b[3];
	double L[3];
	double shell_contr = 0.0;
	double dist=0, rsq=0;;
	double en=0;
	double esum=0;
	double ngrad=0;
	double grad1=0, grad2=0, grad3=0;
	double sigma =1;
	double epsilon=1;
	int numpos=atnumber*3;
	//cubic box
	int nmax=1;
	int nx=0;
	int ny=0;
	int nz=0;
	//cout <<"nupos:" << numpos;
	for (int i=0; i<numpos; i++) {
		lxyz[i]=(*xyz)(i);
		//initialisation of grad!!
		lgrad[i]=0;
	}

	for (int nact=0; nact<=nmax; nact++) {
		for (nx=-nact; nx<=nact; nx++) {
			for (ny=-nact; ny<=nact; ny++) {
				for (nz=-nact; nz<=nact; nz++) {
					//skip cells previously summed over
					if (abs(nx)==nact || abs(ny)==nact || abs(nz)==nact) {
						L[0]= nx * ar[0] + ny * br[0] + nz * cr[0];
						L[1]= nx * ar[1] + ny * br[1] + nz * cr[1];
						L[2]= nx * ar[2] + ny * br[2] + nz * cr[2];
						//cout << "L: <"<<L[0]<< " " << L[1] << " " <<L[2] <<">"
						//<< endl;
						//real space part
						for (int i=0; i<numpos; i += 3) {
							//							a[0]=lxyz[i];
							//							a[1]=lxyz[i+1];
							//							a[2]=lxyz[i+2];
							a[0]=lxyz[i]+L[0];
							a[1]=lxyz[i+1]+L[1];
							a[2]=lxyz[i+2]+L[2];
							//							a[0]=lxyz[i];
							//							a[1]=lxyz[i+1];
							//							a[2]=lxyz[i+2];
							// inner loop
							for (int j=0; j<numpos; j+=3) {
								if (nx==0 && ny==0 && nz==0) {
									//cout << "We are in unit cell..." << endl;
									if (i==j) {
										//cout << "skipping selfinteraction..." << endl;
										continue;
									}
								}
								//small optimization, shift L to a!!!
								//get coordinates
								//								b[0]=lxyz[j]+L[0];
								//								b[1]=lxyz[j+1]+L[1];
								//								b[2]=lxyz[j+2]+L[2];
								//								b[0]=lxyz[j];
								//								b[1]=lxyz[j+1];
								//								b[2]=lxyz[j+2];
								//get distances
								b[0]=a[0]-lxyz[j];
								b[1]=a[1]-lxyz[j+1];
								b[2]=a[2]-lxyz[j+2];
								rsq=pow(b[0], 2)+pow(b[1], 2)+pow(b[2], 2);
								dist = sqrt(rsq);
								//distance cut off
								if (dist > cutoff) {
									//cout << "Cut off too big...";
									continue;
								}
								//############################################################
								//ENERGY
								en = 4 * epsilon * (pow(sigma / dist, 12)
										- pow(sigma / dist, 6));

								esum=esum+en;

								//#############################################################
								//GRADIENTS
								if (gcalc == true) {
									ngrad=24* epsilon* ((-2 * pow(dist / sigma,
											-13) / sigma) + (pow(dist / sigma,
											-7) / sigma));
									//cout << "ngrad: " <<ngrad<<endl;
									//cout << "dist:" <<dist<<endl;
									//			cout << b[0]<<endl;
									//			cout << b[1]<<endl;
									//			cout << b[2]<<endl;
									grad1=grad1 + b[0] /dist * ngrad;
									grad2=grad2 + b[1] /dist * ngrad;
									grad3=grad3 + b[2] /dist * ngrad;

								}
							}//end loop j
							if (gcalc==true) {
								lgrad[i] = lgrad[i] + grad1;
								lgrad[i + 1] = lgrad[i + 1] + grad2;
								lgrad[i + 2] = lgrad[i + 2] + grad3;
								grad1=0;
								grad2=0;
								grad3=0;
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
		shell_contr = esum - shell_contr;
		if (verbose == true) {
			cout.precision(8);
			cout << "Shell #"<<nact<< " contribution, Q real part: "
					<< setw(12) << shell_contr * HARTREE2EV;
			cout << "\tPP part: " << esum* HARTREE2EV << endl;
		}
//		if (abs(shell_contr * HARTREE2EV) < ewald_thresh) {
//			break;
//		}
		shell_contr = esum;
	}

	//cout <<"Energy: "<<esum <<" nx: " <<nx <<" ny: "<<ny << " nz: "<<nz <<endl;
	energy=esum;
	if (pol == true) {
		cout << "energy:"<< energy << endl;
	}
	//grad = new LaVectorDouble(lgrad,numpos);
	for (int i=0; i<numpos; i++) {
		(*grad)(i)=lgrad[i];
		//which gradients are frozen??
		(*grad)(i)=lgrad[i]*(*moveMat)(i);
		//cout << "grad:"<< (*grad)(i)<<endl;
		//cout << lgrad[i]<<endl;
	}
	//gradnorm=Blas_Norm2((*grad));
	gradnorm = grad->norm();

	//gradmaxelem=Blas_Index_Max((*grad));
	grad->maxCoeff(&gradmaxelem);

	//gradinfnorm=Blas_Norm_Inf((*grad));
	//cout << "Energy: "<<esum <<endl;
	calculated=true;
}

void ToolCalcLJ::E_Pol(bool pol, bool verbose, bool gcalc) {
	cout << "This function should never be called..." << endl;
}

ToolCalcLJ::~ToolCalcLJ() {
}
