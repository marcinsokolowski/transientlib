/***************************************************************************\
 **  transientlib : library for identification of short optical flashes 
 **						with the wide field cameras.
 **  This software was written by Marcin Sokolowski ( msok@fuw.edu.pl ) 
 **	it was a substantial part of analysis performed for PHD thesis 
 **  ( Investigation of astrophysical phenomena in short time scales with "Pi of the Sky" apparatus )
 **	it can be used under the terms of GNU General Public License.
 **	In case this software is used for scientific purposes and results are
 **	published, please refer to the PHD thesis submited to astro-ph :
 **
 **		http://arxiv.org/abs/0810.1179
 **
 ** Public distribution was started on 2008-10-31
 **
 ** 
 ** NOTE : some of the files (C files) were created by other developers and 
 **        they maybe distributed under different conditions.
 ** 

 ******************************************************************************
 ** This program is free software; you can redistribute it and/or modify it
 ** under the terms of the GNU General Public License as published by the
 ** Free Software Foundation; either version 2 of the License or any later
 ** version. 
 **
 ** This program is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY; without even the implied warranty of
 ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 ** General Public License for more details. 
 **
 *\**************************************************************************

*/           
#include "calcrot.h"
#include "mathfunc.h"
#include "mathdefs.h"

void CMyCalcRot::CalcOrtoLine( double x1, double y1, double x2, double y2,
               double& a, double& b, double& c )
{
   a = ( x2 - x1 );
   b = ( y2 - y1 );
   c = -0.5*( (x2+x1)*a + (y2+y1)*b );
}

void CMyCalcRot::CalcLine( double x1, double y1, double x2, double y2,
               double& a, double& b, double& c )
{
   a = ( y2 - y1 );
   b = ( x1 - x2 );
   c = (y1*x2-y2*x1);
}


BOOL_T CMyCalcRot::CalcCross( double a, double b, double c,
                double p, double q, double r,
                double& x, double& y )
{
	if(p==0 && q==0)
		return FALSE;

	int p_10 = (int)p*10;
	int q_10 = (int)q*10;
	int a_10 = (int)a*10;
	int b_10 = (int)b*10;

	if(p_10==a_10 && q_10==b_10)
		return FALSE;

	// original formula :
	// x = ( (r/q) - (c/b) )/( (a/b) - (p/q) );
   // y = ( (r/p) - (c/a) )/( (b/a) - (q/p) );


	if(p!=0 && q!=0){
		x = ( b*(r/q) - c )/( a - b*(p/q) );
		y = ( a*(r/p) - c )/( b - a*(q/p) );
	}else{
		if(p==0){
			if(a==0){
				// paralel :
				return FALSE;
			}else{
				x = ( b*(r/q) - c )/( a - b*(p/q) );

				// y = ( r - (c/a)*p )/( p*(b/a) - q );				
				// p=0 ->
				y = - (r)/(q);
			}
		}
		if(q==0){
			if(b==0){
				// paralel :
				return FALSE;
			}else{				
				// x = ( r - q*(c/b) )/( q*(a/b) - p );
				// q=0 ->
				x = -(r)/(p);
				y = ( a*(r/p) - c )/( b - a*(q/p) );
			}
		}
	}
	return TRUE;
}


double CMyCalcRot::CalcRotAngle( double x1, double y1, double x2, double y2,
                     double xs, double ys )
{

	double r = sqrt( CMyMathFunc::mysqr(x1-xs)+CMyMathFunc::mysqr(y1-ys) );

	double sqr_val = CMyMathFunc::mysqr(x1-x2)+CMyMathFunc::mysqr(y1-y2);
	double angle = 2*asin( 0.5*( sqrt( sqr_val )/r ) );
//   double angle = 2*asin( 0.5*( sqrt( (CMyMathFunc::mysqr(x1-x2)+CMyMathFunc::mysqr(y1-y2))/(CMyMathFunc::mysqr(x1-xs)+CMyMathFunc::mysqr(y1-ys)) ) ) );
	

	// convert to center of rotation frame :
	x1 = x1 - xs;
	y1 = y1 - ys;
	x2 = x2 - xs;
	y2 = y2 - ys;

	int sign = 1;

	double sin_alfa_1 = fabs( y1/r );
	double alfa_1_tmp = asin( sin_alfa_1 );
	double alfa_1 = alfa_1_tmp;

	if(x1>=0 && y1>=0)
		alfa_1 = alfa_1_tmp;
	if(x1>=0 && y1<0)
		alfa_1 = (2*PI_VALUE) - alfa_1;
	if(x1<0 && y1>=0)
		alfa_1 = PI_VALUE - alfa_1;
	if(x1<0 && y1<0)
		alfa_1 = PI_VALUE + alfa_1;
	
	double sin_alfa_2 = fabs( y2/r );
	double alfa_2_tmp = asin( sin_alfa_2 );
	double alfa_2 = alfa_2_tmp;

	if(x2>=0 && y2>=0)
		alfa_2 = alfa_2_tmp;
	if(x2>=0 && y2<0)
		alfa_2 = (2*PI_VALUE) - alfa_2;
	if(x2<0 && y2>=0)
		alfa_2 = PI_VALUE - alfa_2;
	if(x2<0 && y2<0)
		alfa_2 = PI_VALUE + alfa_2;
	


	// now we have alfa_1 and alfa_2 on the circle so :

	if(alfa_2 < alfa_1){
		// rotation was CLOCK-WISE - NEGATIVE
		sign = -1;
	}else{
		if( alfa_2==alfa_1 ){
			// means only - sinuses are equal , but not angles :
			if( y1>0 ){
				if( x2>x1 )
					sign = -1;
			}else{
				if( x2<x1 )
					sign = -1;
			}
		}
	}

	angle = angle*sign;

   return angle;
}

