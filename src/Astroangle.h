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
 **        Astroangle.h was originally created in 2003 by 
 **        Bogumil Pilecki (pilecki@astrouw.edu.pl) and Dorota Szczygiel (dszczyg@astrouw.edu.pl)
 **        and latter developed by Marcin Sokolowski 
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


#ifndef ASTROANGLE_H
#define ASTROANGLE_H

#include <string>
#include <math.h>
#include "mathdefs.h"
#include "mystring.h"

using namespace std;

// #include "Header.h"


const double ARCSEC_IN_RAD = 206264.8062;	// 180*60*60/PI
const double SEC_IN_RAD = 13750.98708;		// 12*60*60/PI
#define RAD_TO_HOURS 3.81971863420548807255
#define HOURS_TO_RAD 0.26179938779914940783
// #define DEG_TO_RAD 0.017453293
// #define RAD_TO_DEG 57.29578
#define RAD_TO_ARCSEC 206264.8062

typedef signed char angle_type;

const angle_type ANGLE_NOT_RANGED = -1;	// this angle can NEVER be cut or shifted
const angle_type ANGLE_0_2PI      =  0;	// standard '0 - +2PI' range ('0 - +360' degress), periodic (can be shifted)
const angle_type ANGLE_TYPE_1     =  1;	// '-PI/2 - +PI/2' range ('-90 - +90' degrees), non-periodic (can be cut)
const angle_type ANGLE_TYPE_2     =  2;	// '-PI - +PI' range ('-180 - +180' degrees),  periodic (can be shifted)
const angle_type ANGLE_TYPE_3     =  3;	// '0 - +PI' range ('0 - +180' degrees), non-periodic (can be cut)
const angle_type ANGLE_TIME_TYPE  =  4;	// standard time '0 - 24h' range, periodic (can be shifted)
const angle_type ANGLE_HOURS      =  5;   // hours as 12.35454

const angle_type ANGLE_HA_TYPE    = 10;	// hour angle type ('0 - 24' hours, periodic)
const angle_type ANGLE_RA_TYPE    = 11;	// right ascension type ('0 - 24' hours, periodic)
const angle_type ANGLE_DEC_TYPE   = 12;	// declination type ('-90 - +90' degrees, non-periodic)
const angle_type ANGLE_LAT_TYPE   = 13;	// latitude type ('-90 - +90' degrees, non-periodic)
const angle_type ANGLE_LONG1_TYPE = 14;	// longitude type 1 ('0 - +360' degrees, periodic)
const angle_type ANGLE_LONG2_TYPE = 15;	// longitude type 2 ('-180 - +180' degrees, periodic)
const angle_type ANGLE_AZIM1_TYPE = 16;	// azimuth type 1 ('0 - +360' degrees, periodic)
const angle_type ANGLE_AZIM2_TYPE = 17;	// azimuth type 2 ('-180 - +180' degrees, periodic)
const angle_type ANGLE_ZD_TYPE    = 18;	// zenithal distance type ('0 - +180' degrees, non-periodic)
const angle_type ANGLE_ALT_TYPE   = 19;	// altitude type ('-90 - +90' degrees, non-periodic)
const angle_type ANGLE_DEG_TYPE   = 20;   // normal deg type ( 0 - 360 ) degrees 



/******************************************************************************************************/
/**********************                ASTRO-ANGLE CLASS                 ******************************/
/******************************************************************************************************/


class AstroAngle
{
  private:
	double     angle;	// standard radian representation
	angle_type type;	// type of angle
	//char       sign;	// sign of angle ('+' or '-', 'E' or 'W', 'N' or 'S', depending on angle type)

  public:

	AstroAngle (angle_type type_to_set = ANGLE_0_2PI, double angle_to_set = 0.0);
	AstroAngle (angle_type type_to_set, double main, double min, double sec);

	void setAngle (double main = 0.0, double min = 0.0, double sec = 0.0);	// remember to set correct angle type before setting angle !!!
	void setAngleType (angle_type type_to_set);
	//void setSign (char sign_to_set);

	void setInArcDegrees (int deg, int min, double sec);
	void setInArcDegrees (double deg);
	void setInArcMinutes (double min);
	void setInArcSeconds (double sec);

	void setInHours (int hour, int min, double sec);
	void setInHours (double hour);
	void setInMinutes (double min);
	void setInSeconds (double sec);
	
	void setInRadians (double rad);

	string toString();
	string toString(angle_type in_type);
	static string toString(double _angle, angle_type in_type);

	double inArcDeg();
	double inArcMin();
	double inArcSec();

	double inHours();
	double inMin();
	double inSec();

	double inRad();

	static int arcdeg( double _angle );
	static int arcmin( double _angle );
	static double arcsec( double _angle );

	inline int arcdeg() { return arcdeg(angle); }
	inline int arcmin() { return arcmin(angle); }
	inline double arcsec() { return arcsec(angle); }

	// conversions :
	static double degtime2deg( int deg, int min, double sec );
	static double ra2deg( const char* szRA );
	static double arcsec2rad( double arcsec );
	static double timeangle2rad( double time_angle );
	static double timeangle2rad( int hour, int min, double sec  );	
	static double timeangle2deg( int hour, int min, double sec  );
	static double deg2rad( double in_deg );
	static void deg2deg( double in_deg, int& deg, int& min, int& sec );
	static void deg2deg( double in_deg, int& deg, int& min, double& sec );
	static double degdeg2deg( int deg, int min, double sec );
	static double deg2hours( double in_deg );
	static double rad2deg( double in_rad );
	static double time2hours( int hour, int min, double sec );
	static void rad2timeangle( double in_rad, int& hour, int& min, int& sec );
	static void rad2timeangle_new( double in_rad, int& hour, int& min, double& sec );
	static void rad2timeangle( double in_rad, double& hour, double& min, 
										double& sec );
	static void deg2timeangle( double in_deg, int& hour, int& min, double& sec );										
	static double rad2hours( double in_rad );
	static double hours2rad( double in_hours );
	static double hours2deg( double in_hours );
	
	static double rad2arcsec( double in_rad );
	static mystring rad2degstring( double in_rad );

	inline int hour(){ return hour( angle ); }
	inline int min(){ return min( angle ); }
	inline double sec(){ return sec( angle ); }
	
	static int hour( double _angle );
	static int min( double _angle );
	static double sec( double _angle );

	char getSign ();	// sign of angle ('+' or '-', 'E' or 'W', 'N' or 'S', depending on angle type)-
	static char getSign( double _angle, angle_type _type );
	static string& getSign( double _angle, angle_type _type, string& szSign );
	

	static void cutToRange( double& _angle, angle_type _type );		// eg. for 0-2PI range cuts to 0 if less then 0 and to 2PI when greater
	inline void cutToRange(){ cutToRange( angle, type ); }


	// angular distances :
	static double getDistDeg( double angle1 , double angle2 );
	static double getDist( double angle1 , double angle2 );
	static double getDist( double ra1, double dec1, double ra2, double dec2 );
	
	// comparison
	static int CompareAngle( double angle1, double angle2 );

	static inline void shiftToRangeRAType( double& _angle )
	{
		if ( ( _angle < 0.0 ) || ( _angle >= 2*PI_VALUE ) ){
	   	_angle = _angle - 2*PI_VALUE * ( floor( _angle/(2*PI_VALUE) + 0.000000000001 ) );
	   }
	}	   	
	         	
	inline void shiftToRange (){ shiftToRange( angle, type ); }
	static void shiftToRange( double& _angle, angle_type _type );

	static void setDefaultForPole( double& _angle, angle_type _type, double& out_angle, angle_type _out_type );
	inline void setDefaultForPole ( AstroAngle& angle_to_check )
		{ setDefaultForPole( angle_to_check.angle, angle_to_check.type, angle, type ); }

}; // END OF ASTROANGLE CLASS


#endif

