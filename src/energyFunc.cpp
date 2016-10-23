#include "energyFunc.h"

energyFunc::energyFunc()
{
}

energyFunc::~energyFunc()
{
}

double energyFunc::enLJ(double r, double epsilon, double sigma) {
	double energy = 0;
	energy = 4 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6));
	//energy=4*pow(sigma / r, 12);
	//cout << energy <<endl;
	return energy;
}