#ifndef ASTROCCD_H
#define ASTROCCD_H

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


// #include "Header.h"
#include "Astroangle.h"
#include "ccd_common_defines.h"

	// Orientation constants:
	// (change it if image from CCD Camera is fliped horizontal or vertical)
extern	const double NORTHERN_HEMISPHERE_STANDARD; //   0 degrees
extern const double NORTHERN_HEMISPHERE_UP_WEST;	//  90 degrees
extern const double NORTHERN_HEMISPHERE_UP_EAST;	// -90 degrees
extern const double SOUTHERN_HEMISPHERE_STANDARD;	// 180 degrees
extern const double SOUTHERN_HEMISPHERE_UP_WEST;	// 270 degrees
extern const double SOUTHERN_HEMISPHERE_UP_EAST;	//  90 degrees

extern const double UP_NORTH; //   0 degrees
extern const double UP_SOUTH;	// 180 degrees
extern const double UP_WEST;	//  90 degrees
extern const double UP_EAST;	// -90 degrees

// Notes about convention:
// 1) When both angle coordinates are used in function calling,
//    first goes 2PI_VALUE range coordinate (eg. Longitude, Right Ascension, etc.),
//    then PI range coordinate (eg. Latitude, Altitude, Declination, etc.).
//

/******************************************************************************************************/
/**********************                 ASTRO-CCD CLASS                  ******************************/
/******************************************************************************************************/



class AstroCCD
{
  protected :
  	void initializeAstroAngles();

	eObservationMode_T m_ObsMode;

	// some parameters should only be set by special functions - sometimes 
	// additional calculations are required :
	double Lat;
	double Sin_Lat,Cos_Lat;
	
	double Long;
	double Sin_Long,Cos_Long;
	
	int x_resolution, y_resolution;
	double x_center;
	double y_center;	   
	double focus;
	double pixel_size;		
	double pixel_size_div_focus;
	
	
	// observation characteristics :
	// NOTE : all angles are kept in RADIANS !!!!	   
	int time_zone;
	time_t UT_time;
	double SID_time; // sideral time 

	double Azim_obs;
	double Sin_Azim_obs,Cos_Azim_obs;

	double Alt_obs;
	double Sin_Alt_obs,Cos_Alt_obs;
	
	double RA_obs;
	double Sin_RA_obs,Cos_RA_obs;
	
	double HA_obs;
	double Sin_HA_obs,Cos_HA_obs;
	
	double Dec_obs;
	double Sin_Dec_obs,Cos_Dec_obs;
	
	double orientation; // rotation of CCD with respect to azimutal meridian 
	double eq_orient;
	
	      
	// OBSOLATE :
	int ccd_orient; // ccd coordinates orientation +1 - normal , -1 - miror of Y axis
	   
	// better :
	int ccd_X_orient; // 1/-1 - decides of CCD X axis orientation 
	int ccd_Y_orient; // 1/-1 - decides of CCD Y axis orientation
	            
  public:
	// member variables :

	AstroCCD ();
	AstroCCD( double _focus, double _pixel_size, int _xSize, int _ySize,
				 double _Lat, double _Long, double _Azim_obs, double _Alt_obs,
				 double _orientation, int _X_orient=1, int _Y_orient=-1,
				 eObservationMode_T mode=eNoMovingMode, double _Dec_obs=0, double _RA_obs=0,
				 int _time_zone=0 );

	// changing oberved position :
	void ChangeObsCoordinates( double ra, double dec, double azim, double alt, eObservationMode_T obsMode );
	void GetObsCoordinates( double& ra, double& dec, double& azim, double& alt,
								   eObservationMode_T& obsMode );

	// update of coordinates - no movement change :
	void ChangeObsCoordinates( double ra, double dec, double azim, double alt );
	
	void SetObsHorizontalCoo( double _Azim_obs, double _Alt_obs );
	
	void SetObsEquatorialCoo( double dec, double ra );

	// after change of Lat,Long,Azim_obs, UT_time call :
	void ReCalcParams();				 	
	
	void CalculateAllCoord();


	// functions for updating members :
	void SetUT( time_t ut_time );
	inline time_t GetUT() { return UT_time; }

	inline eObservationMode_T GetObsMode(){ return m_ObsMode; }
	inline void SetObsMode( eObservationMode_T obsMode ){ m_ObsMode=obsMode; }
	
	inline void GetObsEq( double& ra, double& dec ){
		ra = RA_obs;
		dec = Dec_obs;
	}

	inline void GetObsEq( double& ra, double& dec, double& ha ){
		ra = RA_obs;
		dec = Dec_obs;
		ha = HA_obs;
	}
	
	inline void GetObsAzim( double& azim, double& alt ){
		azim = Azim_obs;
		alt  = Alt_obs;
	}
	
	inline void GetObsCoo( double& ra, double& dec, double& azim, double& alt,
								  eObservationMode_T& obsMode )
	{
		ra   = RA_obs;
		dec  = Dec_obs;
		azim = Azim_obs;
		alt  = Alt_obs;			
		obsMode = m_ObsMode;	
	}								  
	
	
	// calculating distances :

	// BOTH FUNCTIONS RETURN angular distance in DEGRESS !!!!!!!!!!!!!!!	

	static double CalcDistInDeg( double ra1, double dec1, double ra2, double dec2 );

	static double CalcDistInRad( double ra1, double dec1, double ra2, double dec2 );
	
	// ra - hours, dec - degrees 
	// RETURNS in DEGREES :	
	static double CalcDistRADEC( double rah1, double dec_deg1, double rah2, double dec_deg2 );
	
	// calculates RA_min and RA_max - range of ra in which frame is contained :
	static int calc_fov_range( double ra0, double dec0, double fov, 
										double& ra_min, double& ra_max,
									   double& dec_min, double& dec_max );
  
	// statics :
	static double getPointZenithalDistance( double alt );

	void setPrimeEquatorialCoordinates(double RA_obs_to_set, double Dec_obs_to_set);
	void setScondaryEquatorialCoordinates(double HA_obs_to_set, double Dec_obs_to_set);



	// HORIZONTAL COORDINATES :
	void calculateHorizontalCoordinatesBase( double dec, double ha, double& alt, double& azim );
	static void calculateHorizontalCoordinatesBaseSlow( double dec, double ha, double lat, double& alt, double& azim );
	
	
	void calculateHorizontalCoordinatesFromEq( double dec, double ra, time_t ut_time, 
													 double& alt, double& azim );
	void calculateHorizontalCoordinatesFromEq( double dec, double ra, 
													 double& alt, double& azim );
	static void calculateHorizontalCoordinatesFromEq( double dec, double ra, time_t ut_time, 
															double& geo_long, double& geo_lat,
													 		double& alt, double& azim );

	inline void calculateHorizontalCoordinatesObs()
		{ 
			calculateHourAngleObs();
			calculateHorizontalCoordinatesBase( Dec_obs, HA_obs, Alt_obs, Azim_obs ); 
			ReCalcParams();				
		}
	

	void calculatePointHorizonatal( double x, double y, 
											  double& h, double& az );

	void calculatePointFromHorizonatal( double h, double az, double& x, double& y );

	// TEST :
	void calculatePointFromEquatorial( double dec, double ra, double& x, double& y );


	// REFRACTION :
	void correctAltitudeForRefraction( double& Alt );
	// void correctForRefraction();

	// EQUATORIAL COORDINATES :
	// from x,y -> RA,DEC
	void calculatePointEquatorial( double x, double y,
	                               double& RA, double& Dec );

	void calculatePointEquatorialNEW( double x, double y,
	                               double& RA, double& Dec );
	                               
	void calculatePointEquatorial2( double x, double y,
											  double& RA, double& Dec, double& HA,
											  double& alt, double& azim );    

	// test - calculating RA, Dec,HA from x,y without first calculating az,h
	void calculatePointEquatorialTest( double x, double y, time_t ut_time,
											  double& RA, double& Dec, double& HA,
											  double& alt, double& azim );    

	void calculatePointEquatorialTest( double x, double y,
											  double& RA, double& Dec, double& HA,
											  double& alt, double& azim );    

	void calculatePointEquatorialTest( double x, double y, double& RA, double& Dec );

	                               
	// from azim,alt -> HA,RA,DEC	             
	// STUPID ORDER OF DEC,RA !!!!!!!!!!!!!!!!!		                  
	void calculateEquatorialCoordinates( double azim, double alt, time_t ut_time, double geo_long,
	                                     double& ha, double& dec, double& ra );

	// STUPID ORDER OF DEC,RA !!!!!!!!!!!!!!!!!
	static void calculateEquatorialCoordinates( double azim, double alt, time_t ut_time, 
													 double geo_long, double geo_lat,
	                                     double& ha, double& dec, double& ra );

	// from azim,alt -> HA,RA,DEC
	// STUPID ORDER OF DEC,RA !!!!!!!!!!!!!!!!!
	inline void calculateEquatorialCoordinates( double azim, double alt, double& ha, double& dec, double& ra )
		{ calculateEquatorialCoordinates( azim, alt, UT_time, Long, ha, dec, ra ); }	

	// STUPID ORDER OF DEC,RA !!!!!!!!!!!!!!!!!		
	inline void calculateEquatorialCoordinates( double azim, double alt, time_t ut_time, double& ha, double& dec, double& ra )
		{ calculateEquatorialCoordinates( azim, alt, ut_time, Long, ha, dec, ra ); }		

	// from Azim_obs,Alt_obs -> RA_obs, DEC_obs
	void calculateEquatorialCoordinatesObs();													 	

	void calculatePointXY( double RA, double Dec, double& x, double& y );

	void calculateEquatorialOrientation();	         
	
	static void calculateRightAscension( double ha, double dec, time_t ut_time, double geo_long, double& ra );
	inline void calculateRightAscension( double ha, double dec, double& ra )
		{ calculateRightAscension( ha, dec, UT_time, Long, ra ); }
	inline void calculateRightAscensionObs(){ 
		calculateRightAscension( HA_obs, Dec_obs, RA_obs ); 
		ReCalcParams();
	}

	static void calculateHourAngleBase( double ra, double dec, double& sid_time, double& ha );
	inline static void calculateHourAngle( double ra, double dec, time_t ut_time, double geo_long, double& ha )
	{
		 double sid_time = getSiderealTimeLocal( ut_time, geo_long );
		 calculateHourAngleBase( ra, dec, sid_time, ha );		    
	}
	inline void calculateHourAngle( double ra, double dec, time_t ut_time, double& ha )
	{
		double sid_time = getSiderealTimeLocal( ut_time, Long );
		calculateHourAngleBase( ra, dec, sid_time, ha );
	}
	inline void calculateHourAngle( double ra, double dec, double& ha )
		{ calculateHourAngle( ra, dec, UT_time, Long, ha ); }
	inline void calculateHourAngleObs(){ 
		calculateHourAngleBase( RA_obs, Dec_obs, SID_time, HA_obs ); 
		ReCalcParams();			
	}
	

	// getting center knowing coordinates for star :
	int CalcCenterFromStar( double ra, double dec, double azim, double alt,
									 time_t ut_time, double x, double y, 
									 double& ra_obs, double& dec_obs, double& azim_obs, double& alt_obs );

	// SHIFTS CALCULATIONS :
	void XYAfterTime (double sec);	// change X and Y to as it should be after some time
	void XYAfterTime (double x_start, double y_start, 
							double& end_x, double& end_y,
							double sec); // same as above, but you can set X and Y at start
							
	// void XYAfterTimeNEW(double x_start, double y_start, double& end_x, double& end_y, double sec);							  
	
	void XYAfterTimeTEST(double x_start, double y_start, double& end_x, double& end_y, double sec);	
	
	void XYAfterTimeTEST( double x, double y, int* prevX, int* prevY, int backStart, int backTo,
	                     const double* PrevFramesTime,int time_sign );
	
	void XYAfterTimeFromEq( double dec, double ra, double& end_x, double& end_y, double sec);
	

	// Standard coordinates :
	double calcStandardCoord( double ra, double dec, double& eta, double& teta );		
	
	//-------------------------------------------------------------------------------------
	// DATE TIME :
	static double getHJD( time_t ut_time, double ra, double dec );
	static double jd2hjd( double jd, double ra, double dec );
	static double getJulianDay( time_t ut_time );
	static double getJulianDay( int year, int month, int day, double hour );
	static double getJulianDay1( int year, int month, int day, double hour );
	
	// julian date to UTC :
	static time_t jd2ux( double jd );

	static double getTJD( time_t ut_time );
	static time_t unixtime_from_tjdsod( int tjd, double sod );
	
	static double JD_to_TJD( double jd );
	
	

	// sideral time :
	static double getSiderealTimeAtGreenwichFromJD( double JD_date );
	
	static double getSiderealTimeAtGreenwich( time_t ut_time );
	inline double getSiderealTimeAtGreenwich(){ return getSiderealTimeAtGreenwich(UT_time); }

	static double getSiderealTimeLocal( time_t ut_time, double geo_long );
	inline double getSiderealTimeLocal(){ return getSiderealTimeLocal( UT_time, Long ); }		
	//----------------------------------------------------------------------------------------	


	                                             
	// TYPICAL FOR CCD :
	// angular <-> pixel distance conversion
	inline double ang2pix(double angle){ return tan(angle)/pixel_size_div_focus; }
	inline double pix2ang(double length){ return atan( length*pixel_size_div_focus ); }
	inline double in_focal( double length ){ return length*pixel_size_div_focus; }


	void getCoordinates( time_t ut_time, double x, double y,
								double& azim, double& altit,
								double& dec, double& ra, double& ha );

	void getCoordinates( double x, double y,
								double& azim, double& altit,
								double& dec, double& ra, double& ha );


	// calculating horizontal from azimuthal :
	// IN :
	// ra - right acscension ( radians )
	// dec - declination ( radians )
	// OUT :
	// h - 
	// az - azimuth 
	static void getHorizontalCoo( double ra, double dec,
											time_t ut_time, double geo_long, double geo_lat,
											double& alt, double& azim );
											
											
	// some standard formulae :
	static void calcHorToEqHour( double h, double az, double phi, double& t, double& dec );
	static void calcEqHourToHoriz( double t, double dec, double phi, double& h, double& az );


	// calculation of galactic <-> equatorial coordinates :
	static int c__1;
	static double c_b11;
	static double c_b13;
	static int c__5;
	static double c_b35;
	
	static int gal2eq(double *gl, double *b, double *raepo,double *depo, double *epoch);
	static int eq2gal2k(double *ra2000, double *d2000, double *gl, double *b);	
	static int prenew(double *dje1, double *dje2, double* ra1, double *d1, double *ra2, double *d2);
	static int eq2gal0(double *raepo, double *depo, double* epoch, double *gl, double *b);
	static int eq2gal(double *raepo, double *depo, double* epoch, double *gl, double *b);
			
												
}; // END OF ASTROCCD CLASS


#endif

