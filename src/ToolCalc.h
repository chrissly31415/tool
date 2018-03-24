#ifndef TOOLCALC_H_
#define TOOLCALC_H_

#define PI 3.1415926535897
#define ANG2BOHR 1.8897162
#define HARTREE2EV 27.211396132
//Global variable defines degree of polynomial: robust settings 4/5 or 3/4
#define fpoly_np 3
#define ls_steps 3

#include <sys/time.h>
#include <Eigen/Core>
#include <math.h>

#include "fitlin.h"
#include "ran.h"

//Calc ---> 	Calc_PP   --> Periodic
//				Calc_LJ   --> Lennard Jones
//				...

using namespace std;

class ToolCalc {
public:
	ToolCalc();
	virtual ~ToolCalc();
	typedef Eigen::Matrix<double, 1, Eigen::Dynamic> VectorXd;
	typedef Eigen::Matrix<int, 1, Eigen::Dynamic> VectorXi;
	typedef Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 3,Eigen::RowMajor>> EigenCoords;

	string jobname;
	VectorXd* xyz;
	VectorXd* grad;
	//atoms which can move (=1)
	VectorXd* moveMat;

	//stores the cosmo segments
	VectorXd segments;
	VectorXi satom;
	VectorXd scharge;
	VectorXd sarea;
	VectorXd maxsig;
	VectorXd meansig;
	int nseg;

	//rotate system for debugging mainly
	double rotation_angle;
	double voxelstep;
	int voxelmode;
	int rotation_axis;

	//properties of calculation
	int nproc;
	int atnumber=0;
	bool calculated;
	double energy;
	//element names
	string* atoms;
	//element numbers
	int* atom_nr;
	//species/type numbers

	//periodic variables
	bool periodic;
	bool dipole_cor;
	double Lx, Ly, Lz, V;
	double alpha;
	double ewald_thresh;
	//lattice vectors in real space
	double ar[3];
	double br[3];
	double cr[3];
	//lattice vector in Fourier space
	double af[3];
	double bf[3];
	double cf[3];

	double gradnorm;
	double gradinfnorm;
	double cutoff;
	double threshhold;
	double max_displacement;

	//step for line search
	double step;
	double temperature;
	double k;		//for LJ
	double COM[3]; //center of mass

	int gradmaxelem;
	int maxiter;
	int global_maxiter;
	int nr_reset;
	//General
	int maxatoms;
	bool verbose;
	double rseed;

	//Parameters for pair potentials
	//species related
	double* q;
	double TotalQ;
	double dipole[3];
	double dipnorm;
	string* core_shell;
	bool shell_modell;
	int nr_shells;
	int nr_cores;
	double scutoff;
	double springc;
	// interaction terms

	double** paraMat;
	int** interMat;
	int nr_species;
	int interactions;

	//Optimization related
	bool converged;
	static string elements[87];
	static double emass[87];
	static double bondiradii[87];

	//Grid related
	int* genom;
	int* genom_compressed;
	int gridp;
	int gx, gy, gz;

	//These functions are implemented in inherited classes
	virtual void E(bool pol = false, bool verbose = false, bool gcalc = true);
	virtual void E_Pol(bool pol = false, bool verbose = false,
			bool gcalc = true);
	virtual void E_periodic(bool pol = false, bool verbose = false, bool gcalc =
			true);
	virtual void E_periodic_Pol(bool pol = false, bool verbose = false,
			bool gcalc = true);
	void setUP_periodic(bool periodic = true, bool pol = false, bool verbose =
			false);
	void append_shells(bool verbose = false);
	void strip_shells(bool verbose = false);
	//Variable energy function pointer
	void (ToolCalc::*p_E)(bool, bool, bool);
	void moveRandom(Ran &myRan, double scale = 0.5, double zscale = 0.0);
	void moveRandomPol(Ran &myRan, double scale = 0.5);
	void rotateCoordinates(int axis, double phi);
	void createRandom(Ran &myRan, double radius);

	void opt(bool pol = false, bool verbose = false);
	void monte_carlo_sampling(Ran &myRan, double distortion);
	void basin_hopping(Ran &myRan, double distortion, double verbose);
	void basin_hopping_par(Ran &myRan, double distortion);
	void setCOM();
	void correctGrad();
	void setDipole(bool);
	void setXYZ(VectorXd xyz_new);
	int getInteraction(int type1, int type2);

	int inline getInteraction_fast(int type1, int type2) {
		int interact = 0;
		//Customized version for binary oxides
		if (type1 != 8 && type2 != 8) {
			interact = 2;
		} else {
			interact = 0;
		}
		if (type1 == 8 && type2 == 8) {
			interact = 1;
		}
		//cout<< " "<<interact<< endl;
		return interact;
	}
	void setcoreshell(bool verbose, bool reset);

	//INLINED FUNCTIONS (does not bring too much gain)
	//this method generates polynomial of arbitrarily degree
	static inline VecDoub fpoly(const Doub x) {
		Int j;
		VecDoub p(fpoly_np);
		p[0] = 1.0;
		for (j = 1; j < fpoly_np; j++)
			p[j] = p[j - 1] * x;
		return p;
	}

	//fits energy to a cubic polynomial
	static inline double fitenergy(const VecDoub xx, const VecDoub EE,
			const int order) {
		const Int npts = 6;
		double sig[] = { 1, 1, 1, 1, 1, 1 };
		VecDoub ssig(npts, sig);
		Fitlin myfit(xx, EE, ssig, &fpoly);
		myfit.fit();
		//	for (int i=0; i<fpoly_np; i++) {
		//		cout << myfit.a[i] << endl;
		//	}
		double min = xx[0];
		double max = xx[ls_steps - 1];
		//cout << "Min: "<< min << " Max:" << max <<endl;
		double dist;
		double grid = 20;
		dist = (max - min);
		//cout << "Dist: "<< dist <<endl;
		double step;
		step = dist / grid;
		VecDoub xs(fpoly_np - 1);
		//xs =
		double minimum = 0;
		//cout << "Value: " << minimum << endl;
		for (int i = 0; i <= grid; i++) {
			min = min + step;
			xs = ToolCalc::fpoly(min);
			double fvalue = 0;
			for (int i = 0; i < xs.size(); i++) {
				//cout << xs[i] << " ";
				//cout << myfit.a[i] << endl;
				fvalue = fvalue + xs[i] * myfit.a[i];

			}
			if (fvalue < minimum) {
				minimum = fvalue;
				max = min;
			}
			//cout << min << " " << fvalue << endl;
		}
		return max;
	}

	static inline double erfc_new(double z) {
		return erf(z);
	}

};

#endif /*TOOLCALC_H_*/
