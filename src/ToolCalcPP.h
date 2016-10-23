#ifndef TOOLCALCPP_H_
#define TOOLCALCPP_H_

#define CHUNKSIZE 4


#include "ToolCalc.h"


using namespace std;

class ToolCalcPP : public ToolCalc
{
public:
	ToolCalcPP();
	virtual ~ToolCalcPP();
	void E(bool pol=false, bool verbose=false, bool gcalc = true);
	void E_Pol(bool pol=true, bool verbose=false, bool gcalc = true );
	void E_periodic(bool pol=false, bool verbose=false, bool gcalc = true);
	void E_periodic_Pol(bool pol=false, bool verbose=false, bool gcalc = true);
	void E_periodic_LC(bool pol=false, bool verbose=false, bool gcalc = true);
	//wird in Basisklasse deklariert
	//double alpha;
	//bool dipole_cor;

	//INLINE

};



#endif /*TOOLCALCPP_H_*/
