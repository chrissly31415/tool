#ifndef TOOLCALCZE_H_
#define TOOLCALCZE_H_

#include "ToolCalc.h"

class ToolCalcZE: public ToolCalc {
public:
	ToolCalcZE();
	virtual ~ToolCalcZE();

	void printSegments(bool unitsBohr=false);
	void defineCubeSize();
	Eigen::Vector3d getCubeOrigin();
	void xyz2grid(int a, int b, int c, bool verbose = false);
	void seg2grid(int a, int b, int c, bool verbose = false);
	void seg2qhull();
	void showgrid();
	void gridout();
	void grid2genom(bool verbose = false);
	void genomout();
	void genom2grid();
	void grid2xyz(bool verbose = false);
	void orderOx();

	Eigen::MatrixXd gridseg;

};

#endif /* TOOLCALCZE_H_ */
