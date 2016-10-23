#ifndef ENERGYFUNC_H_
#define ENERGYFUNC_H_

using namespace std;

#include <math.h>

class energyFunc
{
public:
	energyFunc();
	virtual ~energyFunc();
public:
	static double enLJ(double r, double epsilon, double sigma);
};

#endif /*ENERGYFUNC_H_*/
