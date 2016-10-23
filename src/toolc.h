#ifndef TOOLC_H_
#define TOOLC_H_

//#include <omp.h>

using namespace std;

class toolc {
public:
	toolc();
	virtual ~toolc();
	//static variables
	static string ctype;
	int nproc;
};
#endif /*TOOLC_H_*/
