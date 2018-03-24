#ifndef TOOLCALCZE_H_
#define TOOLCALCZE_H_

#include "ToolCalc.h"
#include "ToolParam.h"

class ToolCalcZE: public ToolCalc {
public:
	ToolCalcZE();
	virtual ~ToolCalcZE();

	Eigen::MatrixXd gridseg;

	Eigen::Vector3d getCubeOrigin();

	void moveSystem2octantI();
	void rotateSegments(int axis, double phi);
	void printSegments(bool unitsBohr=false);
	void defineCubeSize(bool roundit=true);
	void xyz2grid(int a, int b, int c, bool verbose = false);
	void seg2grid(int a, int b, int c, bool verbose = false);
	void seg2voxel(bool verbose = false);
	Eigen::VectorXd voxelize_cosmo(int maxi,int maxj,int maxk, double side_length);
	Eigen::VectorXd voxelize_cosmo_simple(int maxi,int maxj,int maxk,double side_length);
	Eigen::VectorXd voxelize_simple(int maxi,int maxj,int maxk,double side_length);
	void seg2qhull();
	void showgrid();
	void gridout();
	void grid2genom(bool verbose = false);
	void genomout();
	void genom2grid();
	void grid2xyz(bool verbose = false);
	void orderOx();

	static EigenCoords convert2eigen();

};

#endif /* TOOLCALCZE_H_ */
