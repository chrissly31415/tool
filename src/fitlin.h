#ifndef FITLIN_H_
#define FITLIN_H_

#include "nr3.h"

using namespace std;

struct Fitlin {
	Int ndat, ma;
	VecDoub_I &x, &y, &sig;
	VecDoub (*funcs)(const Doub);
	VecBool ia;

	VecDoub a;
	MatDoub covar;
	Doub chisq;

//https://eigen.tuxfamily.org/dox-devel/group__LeastSquares.html
	Fitlin(VecDoub_I &xx, VecDoub_I &yy, VecDoub_I &ssig, VecDoub funks(const Doub)) :
		ndat(xx.size()), x(xx), y(yy), sig(ssig), funcs(funks) {
		ma = funcs(x[0]).size();
		a.resize(ma);
		covar.resize(ma, ma);
		ia.resize(ma);
		for (Int i=0; i<ma; i++)
			ia[i] = true;
	}

	void gaussj(MatDoub_IO &a, MatDoub_IO &b) {
		Int i, icol=0, irow=0, j, k, l, ll, n=a.nrows(), m=b.ncols();
		Doub big, dum, pivinv;
		VecInt indxc(n), indxr(n), ipiv(n);
		for (j=0; j<n; j++)
			ipiv[j]=0;
		for (i=0; i<n; i++) {
			big=0.0;
			for (j=0; j<n; j++)
				if (ipiv[j] != 1)
					for (k=0; k<n; k++) {
						if (ipiv[k] == 0) {
							if (abs(a[j][k]) >= big) {
								big=abs(a[j][k]);
								irow=j;
								icol=k;
							}
						}
					}
			++(ipiv[icol]);
			if (irow != icol) {
				for (l=0; l<n; l++)
					SWAP(a[irow][l], a[icol][l]);
				for (l=0; l<m; l++)
					SWAP(b[irow][l], b[icol][l]);
			}
			indxr[i]=irow;
			indxc[i]=icol;
			if (a[icol][icol] == 0.0)
				throw("gaussj: Singular Matrix");
			pivinv=1.0/a[icol][icol];
			a[icol][icol]=1.0;
			for (l=0; l<n; l++)
				a[icol][l] *= pivinv;
			for (l=0; l<m; l++)
				b[icol][l] *= pivinv;
			for (ll=0; ll<n; ll++)
				if (ll != icol) {
					dum=a[ll][icol];
					a[ll][icol]=0.0;
					for (l=0; l<n; l++)
						a[ll][l] -= a[icol][l]*dum;
					for (l=0; l<m; l++)
						b[ll][l] -= b[icol][l]*dum;
				}
		}
		for (l=n-1; l>=0; l--) {
			if (indxr[l] != indxc[l])
				for (k=0; k<n; k++)
					SWAP(a[k][indxr[l]], a[k][indxc[l]]);
		}
	}
	
	void hold(const Int i, const Doub val) {
		ia[i]=false;
		a[i]=val;
	}
	void free(const Int i) {
		ia[i]=true;
	}

	void fit() {
		Int i, j, k, l, m, mfit=0;
		Doub ym, wt, sum, sig2i;
		VecDoub afunc(ma);
		for (j=0; j<ma; j++)
			if (ia[j])
				mfit++;
		if (mfit == 0)throw("lfit: no parameters to be fitted");
		MatDoub temp(mfit, mfit, 0.), beta(mfit, 1, 0.);
		for (i=0; i<ndat; i++) {
			afunc = funcs(x[i]);
			ym=y[i];
			if (mfit < ma) {
				for (j=0; j<ma; j++)
					if (!ia[j])
						ym -= a[j]*afunc[j];
			}
			sig2i=1.0/SQR(sig[i]);
			for (j=0, l=0; l<ma; l++) {
				if (ia[l]) {
					wt=afunc[l]*sig2i;
					for (k=0, m=0; m<=l; m++)
						if (ia[m])
							temp[j][k++] += wt*afunc[m];
					beta[j++][0] += ym*wt;
				}
			}
		}
		for (j=1; j<mfit; j++)
			for (k=0; k<j; k++)
				temp[k][j]=temp[j][k];
		gaussj(temp, beta);
		for (j=0, l=0; l<ma; l++)
			if (ia[l])
				a[l]=beta[j++][0];
		chisq=0.0;
		for (i=0; i<ndat; i++) {
			afunc = funcs(x[i]);
			sum=0.0;
			for (j=0; j<ma; j++)
				sum += a[j]*afunc[j];
			chisq += SQR((y[i]-sum)/sig[i]);
		}
		for (j=0; j<mfit; j++)
			for (k=0; k<mfit; k++)
				covar[j][k]=temp[j][k];
		for (i=mfit; i<ma; i++)
			for (j=0; j<i+1; j++)
				covar[i][j]=covar[j][i]=0.0;
		k=mfit-1;
		for (j=ma-1; j>=0; j--) {
			if (ia[j]) {
				for (i=0; i<ma; i++)
					SWAP(covar[i][k], covar[i][j]);
				for (i=0; i<ma; i++)
					SWAP(covar[k][i], covar[j][i]);
				k--;
			}
		}
	}
};


#endif 
