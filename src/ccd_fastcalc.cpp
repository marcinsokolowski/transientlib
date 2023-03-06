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
#include "ccd_fastcalc.h"
#include <math.h>


CCDFastCalc::CCDFastCalc()
{}



int CCDFastCalc::GetPixelNoInRedial( double r )
{
	int count=0;
	register int rp1 = (int)(r+1);
	for(register int x=-rp1;x<=rp1;x++){
		for(register int y=-rp1;y<=rp1;y++){
			double dist = sqrt(x*x+y*y);
			if(dist<=r)
				count++;
		}
	}
	return count;
}

double CCDFastCalc::CalcSphericity( double r, double nPixels )
{
	//double nPixelsInSphere = GetPixelNoInRedial( r );
	double nPixelsInSphere = 3.1415*r*r;
	if(nPixelsInSphere==0)
		nPixelsInSphere = 1;
	double spher = ( nPixels/nPixelsInSphere );
	return spher;
}


double CCDFastCalc::FindMaxRedial( LONG_T* cluster, int cluster_cnt, double x0, double y0, int xSize )
{
	double max_redial = 0;
	for(register int i=0;i<cluster_cnt;i++){
		double xx = (cluster[i] % xSize);
		double yy = (cluster[i] / xSize);

		double x = xx-0.5;
		double y = yy-0.5;		
		double r = sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0) );
		if(r>max_redial)
			max_redial = r;

		x = xx+0.5;
		y = yy-0.5;
		r = sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0) );
		if(r>max_redial)
			max_redial = r;

		x = xx-0.5;
		y = yy+0.5;
		r = sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0) );
		if(r>max_redial)
			max_redial = r;

		x = xx+0.5;
		y = yy+0.5;
		r = sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0) );
		if(r>max_redial)
			max_redial = r;


	}
	return max_redial;
}

void CCDFastCalc::FindCenter( LONG_T* cluster, int cluster_cnt, int xSize,
			      double& x0, double& y0 )
{
	double x_sum=0.00,y_sum=0.00;
	for(int i=0;i<cluster_cnt;i++){
		int x = (cluster[i] % xSize);
		int y = (cluster[i] / xSize);
		x_sum += x;
		y_sum += y;
	}
	x0 = x_sum / cluster_cnt;
	y0 = y_sum / cluster_cnt;
}

double CCDFastCalc::CalcSphericity( LONG_T* cluster, int cluster_cnt, int xSize, double& r )
{
	double x0,y0;
	FindCenter( cluster, cluster_cnt, xSize, x0, y0 );
	r = CCDFastCalc::FindMaxRedial( cluster, cluster_cnt, x0, y0, xSize );
   double spher = CCDFastCalc::CalcSphericity( r, cluster_cnt );
	return spher;
}