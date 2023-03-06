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
#ifndef _CCD_ASTRO_FORMULAE_H__
#define _CCD_ASTRO_FORMULAE_H__

#include <math.h>

class AstroCCD;

class CCDAstroFormulas
{
public:
	static AstroCCD* m_pAstroFormulas;

	void InitAstroCCD( double focus, double pixel_size, 
							 int sizeX, int sizeY, 
							 int day, int month, int year, int time_zone );


	CCDAstroFormulas();
	~CCDAstroFormulas();
	
	static void xy2rd (double x, double y, double& Dec, double& RA, 
					double Dec_obs, double RA_obs, double orientation,
		  	      double x_center, double y_center, double focus, 
		  	      double pixel_size);

	static void rd2xy (double Dec, double RA, double& x, double& y, 
					double Dec_obs, double RA_obs, double orientation,
	            double x_center, double y_center, double focus, 
	            double pixel_size);
	            
	static void xyAfterTime(double x, double y, double timeInSec,
					double& x_new, double& y_new, 
					double Dec_obs, double RA_obs,
					double orientation, double x_center, double y_center, 
					double focus, double pixel_size);
	            

	static double ang2pix (double angle, double focus, double pixel_size);
	static double pix2ang (double length, double focus, double pixel_size);
	static double deg2rad (double degrees, double minutes, double seconds);
	static double deg2rad (double degrees);
	static double time2rad (double hours, double minutes, double seconds);
	static double time2rad (double hours);
	static double rad2deg (double rad, int& degrees, int& minutes, int& seconds);
	static double rad2deg (double rad);
	static double rad2time (double rad, int& hours, int& minutes, int& seconds);
	static double rad2time (double rad);

};



#endif

