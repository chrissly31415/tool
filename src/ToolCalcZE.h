#ifndef TOOLCALCZE_H_
#define TOOLCALCZE_H_

#include "ToolCalc.h"

class ToolCalcZE: public ToolCalc {
public:
	ToolCalcZE();
	virtual ~ToolCalcZE();

	void printSegments();
	void defineCubeSize();
	void xyz2grid(int a, int b, int c, bool verbose = false);
	void seg2grid(int a, int b, int c, bool verbose = false);
	void gridout();
	void grid2genom(bool verbose = false);
	void genomout();
	void genom2grid();
	void grid2xyz(bool verbose = false);
	void orderOx();

};

#endif /* TOOLCALCZE_H_ */
