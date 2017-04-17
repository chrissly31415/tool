/*
 * ToolParam.h
 *
 *  Created on: Apr 8, 2017
 *      Author: loschen
 */

#ifndef TOOLPARAM_H_
#define TOOLPARAM_H_

#include <map>

class ToolParam
{
std::map<std::string, double> vdwradius;
public:
	ToolParam() {
//		vdwradius['C']=1.6;
//	    vdwradius['H']=1.1;
	}
	virtual ~ToolParam();

};



#endif /* TOOLPARAM_H_ */
