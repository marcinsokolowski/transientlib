#include "gendistr.h"
#include "mathdefs.h"
#include <random.h>
#include <cexcp.h>

double Gauss_func( double r, double n, double s){
	Assert(s>0,"r0 must be greater then 0");
	double ret = double(n)*exp(-(r*r)/(2*s*s));
	return ret;
}

double Gauss_func( double x, double y, double r0, double max_value )
{
	Assert(r0>0,"r0 must be greater then 0");
   double ret = double(max_value)*
                   exp(-(sqrt(x*x+y*y))/r0);
   return ret;
}


void CGenDistr::GetGauss(  double brightness, double width, double x0, double y0,
                           double& x, double& y )
{
	CRandom rnd;
	double z=1;
	double val=0;
	double phi,r;
	double x_new,y_new;
	
	while(z>val){
		phi = 2*PI_VALUE*rnd.GetRandom();
		r = 10*width*rnd.GetRandom();		

		x_new = x0 + r*cos(phi);
		y_new = y0 + r*sin(phi);
		val = Gauss_func( r , brightness, width );
		z = brightness*rnd.GetRandom();
	}
	x = x_new;
	y = y_new;			
}

void   CGenDistr::GetGauss( double n, double r0, double x0, double x1,
                            double y0, double y1, double& x, double& y )
{
	CRandom rnd;
	double z=1;
	double val=0;
	double x_new,y_new;
	
	while(z>val){
		x_new = x0 + rnd.GetRandom()*(x1-x0);
		y_new = y0 + rnd.GetRandom()*(y1-y0);
		val = Gauss_func( x_new, y_new, r0, n );
		z = n*rnd.GetRandom();
	}
	x = x_new;
	y = y_new;	
}


double CGenDistr::GetDistr( double (*func)(double,double) )
{
	return 0;
}
