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
		phi = 2*PI*rnd.GetRandom();
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
