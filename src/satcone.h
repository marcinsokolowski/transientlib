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
#ifndef _SAT_CONE_H__
#define _SAT_CONE_H__

#include <mydate.h>

// all in [km]
#define EARTH_RADIUS  6378 
#define SUN_RADIUS    696000
#define SUN_DIST      149600000

// in [deg] :
#define CONE_ANGLE 0.264121516

// ??? o co chodzi - jak dam R_s = 696000 to wychodzi to samo !!!???
// #define SUN_CONE_H    1383582
// #define SUN_CONE_C    0.0046097728

extern double gE;
extern double gC;

class CSatCone 
{
public:
	CSatCone(){}
	~CSatCone(){}
	
	static int calc_min_dist( double ra_in_deg, double dec_in_deg, time_t unix_time,
									  double& r_min );

	static int get_earth_centered_line( double ra_in_deg, double dec_in_deg, time_t unix_time, double d );

	static int get_local_line( double ra_in_deg, double dec_in_deg, time_t unix_time, double d );

	static int calc_cone_line_crossing();
}; 


#endif