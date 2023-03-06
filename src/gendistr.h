#ifndef _GENDISTR_H__
#define _GENDISTR_H__

#include <math.h>
#include <mytypes.h>

class CGenDistr
{
public:
	CGenDistr(){};
	
	void GetGauss(  double n, double r0, double x0, double x1,
	                double y0, double y1, double& x, double& y );

	void GetGauss(  double brightness, double width, double x0, double y0,
	                double& x, double& y );

	double GetDistr( double (*func)(double,double) );
};


#endif
