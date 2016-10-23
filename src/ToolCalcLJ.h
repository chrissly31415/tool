#ifndef TOOLCALCLJ_H_
#define TOOLCALCLJ_H_

#include "ToolCalc.h"

class ToolCalcLJ : public ToolCalc
{
public:
	ToolCalcLJ();
	virtual ~ToolCalcLJ();
	void E(bool pol=false, bool verbose=false, bool gcalc = true);
	void E_periodic(bool pol=false, bool verbose=false, bool gcalc = true);
	void E_Pol(bool pol=true, bool verbose=false, bool gcalc = true );
};

#endif /*TOOLCALCLJ_H_*/
