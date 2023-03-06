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
#ifndef _ASTRO_BODY_H__
#define _ASTRO_BODY_H__

#include "AstroCCD.h"
#include <vector>
#include <mystring.h>
#include <sat_interface.h>

using namespace std;

class CBodyPosition
{
public:
	CBodyPosition( time_t t, double _ra, double _dec, double _az, double _h, 
						double _illum=0.00 );


	time_t unix_time;
	double ra;     // RADIANS
	double dec;    // RADIANS
	double h;      // RADIANS
	double az;     // RADIANS
	double illum;
	double mag;
};

class CBodyPositions : public vector<CBodyPosition>
{
public :
	CBodyPositions(){}

	void Add( time_t unix_time, double ra, double dec, double az, double h, 
				 double illum=0.00 );	
	CBodyPosition* FindRise( int& pos );
	CBodyPosition* FindSet( int& pos );
};

class CAstroBody : public AstroCCD
{
protected:
	static double calc_sidt( double JD_date );

public:
	// time step when calculating object positions
	static int m_CalcPosStepInSec;
	
	static double calcDistDeg( double ra1, double dec1, double ra2, double dec2 );
	

	CAstroBody();
	CAstroBody( const char* szName, double ra=0, double dec=0, double fov=0.1 );

	// specialy optimized for SWIFT :
	// getting coordinates in case satellite information is provided 
	// it checks pre-defined list of satellite positions ,
	// returns TRUE if satellite poinint was used and FALSE if just body_ra,body_dec
	// was used :
	BOOL_T GetRADEC( time_t unix_time, double& ra_in_rad, double& dec_in_rad,
						double& ra_in_deg, double& dec_in_deg, 
						double lat, double longit, double min_acceptable_alt,
						time_t min_track_time=600, double max_dist=70.00,
						time_t max_future_time=900 );

	// simple position retriving :						
	BOOL_T GetSimpleRADEC( time_t unix_time, double& ra_in_rad, double& dec_in_rad,
	                 double& ra_in_deg, double& dec_in_deg, time_t& start_time,
	                 time_t& end_time );

	void ShowSkyPath();

	mystring m_szName;
	CBodyPositions m_PosList;

	double body_ra;  // radians
	double body_dec; // radians

	double body_ra_in_deg;  // degrees 
	double body_dec_in_deg; // degrees
	
	
	// information about satellite from WWW :
	CSatInterface* m_pSatInfo;
	
	double illum;    // illum (mainly for MOON)

	double m_FOV; // Field Of View - in degrees 
	
	// times of rise / set :
	
	// RISE/SET above/below horizont ( h=0 ):
	time_t m_hor_rise_ut;
	time_t m_hor_set_ut;

	double m_hor_rise_jd;
	double m_hor_set_jd;
	
	// RISE/SET above/below minimum h for observation :
	time_t m_obs_rise_ut;
	time_t m_obs_set_ut;
	
	double m_obs_rise_jd;
	double m_obs_set_jd;
		
	// sun - returns ra,dec in DEGREES 
	static int calc_sun( double j_date, double &ra, double &dec );	
	static int calc_sun2( double j_date, double &ra, double &dec );

	// calculate height ( in DEGREES ) :
	static double calc_height( double JD_date , double lat, double lon,
	       			            double body_ra, double body_dec );

	static double when( double j_date, double height, double ra, double dec, double la, 
							  double lo, double& ah );
	static double when_new( double j_date, double height, double ra, double dec, double la, 
							  double lo, double& ah, double prec );
	static double when_sun( double j_date, double height, double ra, double dec, double la, double lo );

	CBodyPosition* FindRiseAbove( double h );
	CBodyPosition* FindSetBelow( double h );
	
	CBodyPosition* SetMidnightRADEC( time_t ut_time );
	CBodyPosition* GetPosition( time_t ut_time );
	
	BOOL_T CalcPositions( time_t dtm_start, time_t dtm_end,
								 double long_in_deg, double lat_in_deg );
								 
									 
};


#endif
