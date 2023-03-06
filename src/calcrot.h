#ifndef _MY_CALC_ROT_H
#define _MY_CALC_ROT_H

#include <mytypes.h>

class CMyCalcRot
{
public:
	static void CalcOrtoLine( double x1, double y1, double x2, double y2,
               	 	     	  double& a, double& b, double& c );

	static void CalcLine( double x1, double y1, double x2, double y2,
	                      double& a, double& b, double& c );

	static BOOL_T CalcCross( double a, double b, double c,
			                   double p, double q, double r,
         			          double& x, double& y );

	static double CalcRotAngle( double x1, double y1, double x2, double y2,
			                      double xs, double ys );

};

#endif
