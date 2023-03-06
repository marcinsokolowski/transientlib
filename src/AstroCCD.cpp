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


#include "AstroCCD.h"
#include <stdio.h>
#include <mathfunc.h>
#include <myutil.h>

extern "C" {
#include <astutil.h>
// extern double d_mod(double *, double *);
}

// #include <math.h>

// Orientation constants:
// (change it if image from CCD Camera is fliped horizontal or vertical)
const double NORTHERN_HEMISPHERE_STANDARD =     0.0;	//   0 degrees
const double NORTHERN_HEMISPHERE_UP_WEST  =  PI_2_VALUE;	//  90 degrees
const double NORTHERN_HEMISPHERE_UP_EAST  = -PI_2_VALUE;	// -90 degrees
const double SOUTHERN_HEMISPHERE_STANDARD =      PI_VALUE;	// 180 degrees
const double SOUTHERN_HEMISPHERE_UP_WEST  =  1.5*PI_VALUE;	// 270 degrees
const double SOUTHERN_HEMISPHERE_UP_EAST  =  PI_2_VALUE;	//  90 degrees

const double UP_NORTH =     0.0;	//   0 degrees
const double UP_SOUTH =      PI_VALUE;	// 180 degrees
const double UP_WEST  =  PI_2_VALUE;	//  90 degrees
const double UP_EAST  = -PI_2_VALUE;	// -90 degrees


AstroCCD::AstroCCD()
{
	initializeAstroAngles();
	CalculateAllCoord();
}

AstroCCD::AstroCCD( double _focus, double _pixel_size, int _xSize, int _ySize,
          		     double _Lat, double _Long, double _Azim_obs, double _Alt_obs,
		              double _orientation, int _X_orient, int _Y_orient,
						  eObservationMode_T mode, double _Dec_obs, double _RA_obs,
						  int _time_zone )
: focus( _focus ), pixel_size(_pixel_size), x_resolution(_xSize),
  y_resolution(_ySize),Lat(_Lat), Long(_Long),Azim_obs(_Azim_obs),
  Alt_obs( _Alt_obs ), orientation(_orientation), time_zone( _time_zone ), ccd_orient(-1),
  ccd_X_orient(_X_orient), ccd_Y_orient(_Y_orient), m_ObsMode( mode ), 
  Dec_obs( _Dec_obs ), RA_obs( _RA_obs ),UT_time(0),
  HA_obs(0),eq_orient(0.00),x_center(0),y_center(0),SID_time(0)
{
	CalculateAllCoord();
}


void AstroCCD::initializeAstroAngles ()
{
	// CCD angles and time
	Lat         = 0.00;
	Long        = 0.00;
	Azim_obs    = 0.00;
	Alt_obs     = 0.0;
	RA_obs      = 0.0;
	HA_obs      = 0.0;
	Dec_obs     = 0.0;
	orientation = 0.0;
	eq_orient   = 0.0;
	UT_time     = 0;
	ccd_orient  = -1;
	ccd_X_orient = 1;
	ccd_Y_orient = -1;
	x_center = 0;
	y_center = 0;
	m_ObsMode = eNoMovingMode;
	SID_time = 0;
}

void AstroCCD::ReCalcParams()
{
	x_center = ( x_resolution + 1.0 ) / 2.0;
	y_center = ( y_resolution + 1.0 ) / 2.0;

	Sin_Lat = sin(Lat);
	Cos_Lat = cos(Lat);
	
	Sin_Long = sin(Long);
	Cos_Long = cos(Long);

	Sin_Azim_obs = sin(Azim_obs);
	Cos_Azim_obs = cos(Azim_obs);

	Sin_Alt_obs = sin(Alt_obs);
	Cos_Alt_obs = cos(Alt_obs);

	Sin_RA_obs = sin(RA_obs);
	Cos_RA_obs = cos(RA_obs);

	Sin_HA_obs = sin(HA_obs);
	Cos_HA_obs = cos(HA_obs);

	Sin_Dec_obs = sin(Dec_obs);
	Cos_Dec_obs = cos(Dec_obs);	

	pixel_size_div_focus = (pixel_size/focus);
}

void AstroCCD::CalculateAllCoord()
{
	if( m_ObsMode==eNoMovingMode ){
		SetObsHorizontalCoo( Azim_obs, Alt_obs );
	}
	if( m_ObsMode==eEarthMovingMode ){
		SetObsEquatorialCoo( Dec_obs, RA_obs );
	}
}

void AstroCCD::GetObsCoordinates( double& ra, double& dec, double& azim, double& alt,
										   eObservationMode_T& obsMode )
{
	ra = RA_obs;
	dec = Dec_obs;
	azim = Azim_obs;
	alt = Alt_obs;
	obsMode = m_ObsMode;
}


void AstroCCD::ChangeObsCoordinates( double ra, double dec, double azim, double alt )
{
	RA_obs = ra;
	Dec_obs = dec;
	Alt_obs = alt;
	Azim_obs = azim;	

	if( m_ObsMode==eNoMovingMode ){
		SetObsHorizontalCoo( azim, alt );
	}
	if( m_ObsMode==eEarthMovingMode ){
		SetObsEquatorialCoo( dec, ra );
	}
	ReCalcParams();
}


void AstroCCD::ChangeObsCoordinates( double ra, double dec, double azim, double alt, eObservationMode_T obsMode )
{
	m_ObsMode = obsMode;
	if( m_ObsMode==eNoMovingMode ){
		SetObsHorizontalCoo( azim, alt );
	}
	if( m_ObsMode==eEarthMovingMode ){
		SetObsEquatorialCoo( dec, ra );
	}
	ReCalcParams();
}

void AstroCCD::SetObsHorizontalCoo( double _Azim_obs, double _Alt_obs )
{
	Azim_obs = _Azim_obs;
	Alt_obs  = _Alt_obs;	
	ReCalcParams();

	// calls also ReCalcParams after calculation of coordinates :
	if( UT_time>0 ){
		calculateEquatorialCoordinatesObs();
	}
}

void AstroCCD::SetObsEquatorialCoo( double dec, double ra )
{
	Dec_obs = dec;
	RA_obs  = ra;	
	ReCalcParams();
	if( UT_time>0 ){
		calculateHourAngleObs();

		// calls also ReCalcParams after calculation of coordinates :
		calculateHorizontalCoordinatesObs();
	}
}


void AstroCCD::SetUT( time_t ut_time )
{
	UT_time = ut_time;
	SID_time = getSiderealTimeLocal();

	// update also observation coordinates (RA_obs,Dec_obs - when Azim_obs,Alt_obs=constant
	// or oposite ) :
	if(m_ObsMode!=eEarthMovingMode){
		calculateEquatorialCoordinatesObs();
	}
}

void AstroCCD::getCoordinates( time_t ut_time, double x, double y,
 			     						 double& azim, double& altit,
										 double& dec, double& ra, double& ha )
{
	
}

void AstroCCD::getCoordinates( double x, double y,
										 double& azim, double& altit,
										 double& dec, double& ra, double& ha )
{
	getCoordinates( UT_time, x, y, azim, altit, dec, ra, ha );
}







/*----------------------------   HORIZONTAL COORDINATES   -------------------------------------*/
double AstroCCD::getPointZenithalDistance( double alt )
{
	double z = ( PI_2_VALUE - alt );
	AstroAngle::cutToRange( z, ANGLE_ZD_TYPE );
   return z;
}



/*----------------------------   EQUATORIAL COORDINATES   -------------------------------------*/
void AstroCCD::setPrimeEquatorialCoordinates( double RA_obs_to_set, 
														    double Dec_obs_to_set)
{
	RA_obs = RA_obs_to_set;
	Dec_obs = Dec_obs_to_set;

	calculateHourAngleObs();
	//calculateHorizontalCoordinates();
}


void AstroCCD::setScondaryEquatorialCoordinates( double HA_obs_to_set, 
																 double Dec_obs_to_set)
{
	HA_obs = HA_obs_to_set;
	Dec_obs = Dec_obs_to_set;

	calculateRightAscensionObs();
	//calculateHorizontalCoordinates();
}



/*----------------------------   TIME AND DATE   -------------------------------------*/





/*----------------------------   CALCULATIONS  -------------------------------------*/


void AstroCCD::calculateEquatorialCoordinates( double azim, double alt, time_t ut_time, double geo_long,
														     double& ha, double& dec, double& ra )
{
	AstroAngle::cutToRange( alt, ANGLE_ALT_TYPE );

	if( Lat > PI_2_VALUE ){
	  ha = azim;
	  dec = alt;
	}else{
		if( Lat < -PI_2_VALUE ){
	   	ha = azim;
	 		dec = - alt;
		}else{
			// Declination is...
			double sin_val = Sin_Lat*sin(alt) - Cos_Lat * cos(alt) * cos(azim);
			dec =  asin ( sin_val );

			// Hour Angle is...
			ha = atan2 ( cos(alt)*sin(azim)*cos(dec) * Cos_Lat,
		          ( sin(alt) - sin(dec) * Sin_Lat )* cos(dec) );

			//double tan_tmp = ( cos(alt)*sin(azim) )/( cos(Lat)*sin(alt)+sin(Lat)*cos(alt)*cos(azim) );

			//double ha = atan2 ( cos(alt)*sin(azim),
			//						(cos(Lat)*sin(alt)+sin(Lat)*cos(alt)*cos(azim)) );
			// printf("ha=%.5f, ha2=%.5f\n",ha,ha2);
	 	}
	}

	// but hour angle cannot be less then zero and greater then 24h, so...
	AstroAngle::shiftToRange( ha, ANGLE_HA_TYPE );

	// and if Declination is 90 or -90 degrees there is a default Hour Angle of 0h 0m 0s
	AstroAngle::setDefaultForPole( dec, ANGLE_DEC_TYPE, ha, ANGLE_HA_TYPE );

	// now, Right Ascension can be calculated...
	calculateRightAscension( ha, dec, ut_time, geo_long, ra );

	// equatorial orientation can be obtained now...
	// calculateEquatorialOrientation();
}

void AstroCCD::calculateEquatorialCoordinates( double azim, double alt, time_t ut_time, 
															  double geo_long, double geo_lat,
														     double& ha, double& dec, double& ra )
{
	AstroAngle::cutToRange( alt, ANGLE_ALT_TYPE );

	if( geo_lat > PI_2_VALUE ){
	  ha = azim;
	  dec = alt;
	}else{
		if( geo_lat < -PI_2_VALUE ){
	   	ha = azim;
	 		dec = - alt;
		}else{
			// Declination is...
			double sin_val = sin(geo_lat)*sin(alt) - cos(geo_lat) * cos(alt) * cos(azim);
			dec =  asin ( sin_val );

			// Hour Angle is...
			ha = atan2 ( cos(alt)*sin(azim)*cos(dec) * cos(geo_lat),
		          ( sin(alt) - sin(dec) * sin(geo_lat) )* cos(dec) );

	 	}
	}

	// but hour angle cannot be less then zero and greater then 24h, so...
	AstroAngle::shiftToRange( ha, ANGLE_HA_TYPE );

	// and if Declination is 90 or -90 degrees there is a default Hour Angle of 0h 0m 0s
	AstroAngle::setDefaultForPole( dec, ANGLE_DEC_TYPE, ha, ANGLE_HA_TYPE );

	// now, Right Ascension can be calculated...
	calculateRightAscension( ha, dec, ut_time, geo_long, ra );

	// equatorial orientation can be obtained now...
	// calculateEquatorialOrientation();
}



void AstroCCD::calculateEquatorialCoordinatesObs()
{
	calculateEquatorialCoordinates( Azim_obs, Alt_obs, UT_time, Long, HA_obs, Dec_obs, RA_obs );

	// first re-calculate sin,cos of Obs point :
	ReCalcParams();

	// then calc Equatorial orientation
	calculateEquatorialOrientation();
}

void AstroCCD::calculatePointEquatorial2( double x, double y,double& RA, 
														double& Dec, double& HA,
														double& alt, double& azim )
{
	calculatePointHorizonatal( x, y, alt, azim );
	correctAltitudeForRefraction( alt );
	calculateEquatorialCoordinates( azim, alt, UT_time, Long, 
											  HA, Dec, RA );	
}


void AstroCCD::calculatePointFromEquatorial( double dec, double ra, double& x, double& y )
{
	// to calculate only once :
	double sin_dec =  sin(dec);
	double cos_dec =  cos(dec);

	double cos_alfa = sin_dec*Sin_Dec_obs+cos_dec*Cos_Dec_obs*cos(ra-RA_obs);
	double alfa = acos(cos_alfa);
	double sin_alfa = sin(alfa);
	double sin_sum = (sin(ra-RA_obs)*cos_dec)/sin_alfa;
	double sum = asin(sin_sum);

	// now check quarter - as asin gives only results in range : -PI/2,PI/2 
	double cos_sum = -( ( sin_dec  - Sin_Dec_obs*cos_alfa )/(Cos_Dec_obs*sin_alfa) );
	if( cos_sum<0 && fabs(sum)<PI_2_VALUE){
		// in such case angle obtained from asin(sin_sum) is not correct
		// because cos is negative so this is another quarter :
		sum = PI_VALUE - sum;
	}
	double beta = sum - Lat;

	beta = beta - orientation;

	//if( dec-Dec_obs < 0 )	
	//	beta = TWO_PI_VALUE - beta;

	if(beta<0)
		beta = TWO_PI_VALUE + beta;
	
	double tg_beta = tan(beta);
	double rx = fabs( (focus/pixel_size*tan(alfa) )/(sqrt(1+tg_beta*tg_beta)) );
	double ry = fabs( rx*tg_beta );

	
	if(beta>PI_2_VALUE && beta<PI_VALUE){
		rx = - rx;
	}else{
		if(beta>PI_VALUE && beta<1.5*PI_VALUE){
			rx = -rx;
			ry = -ry;
		}else{
			if(beta>1.5*PI_VALUE && beta<TWO_PI_VALUE)
				ry = -ry;
		}	
	}

	ry = ccd_Y_orient*ry;
   rx = ccd_X_orient*rx;


	x = ( rx + x_center );
   y = ( ry + y_center );
	
}

void AstroCCD::calculatePointEquatorialTest( double x, double y, time_t ut_time,
					  double& RA, double& Dec, double& HA,
					  double& alt, double& azim )
{
	SetUT( ut_time );
	calculatePointEquatorialTest( x, y, RA, Dec, HA, alt, azim );
}


void AstroCCD::calculatePointEquatorialTest( double x, double y, double& RA, double& Dec )
{
double rx, ry, distance;
double frame_rotation_angle;
	
	//distance from center of image :
	rx = ( x - x_center );
	ry = ( y - y_center );

	// TODO :
	ry = ccd_Y_orient*ry;
	rx = ccd_X_orient*rx;

	distance = sqrt( rx*rx + ry*ry );

	//angular distance in radians:
	double ad = pix2ang( distance );
	double sin_ad = sin(ad);

	//there are two great cirlces that goe straight across the center
	//of the picture horizontaly and verticaly;
	//we rotate these circles clockwise until the "vertical" circle
	//eats our point of interest
	double beta = atan2(ry, rx) + orientation;


	// good for 20030202 - with orientation=0.6 
	// double sin_delta = Sin_Dec_obs*cos(ad) - sin_ad*Cos_Dec_obs*cos(Lat+beta );
	// TEST - should be for paralactic mount :
	double sin_delta = Sin_Dec_obs*cos(ad) + sin_ad*Cos_Dec_obs*cos(PI_2_VALUE-beta);

	Dec = asin(sin_delta);

	// good for 20030202 - with orientation=0.6
	// double sin_dRA = sin_ad*sin( Lat+beta )/cos(Dec);
	// TEST - should be for paralactic mount :
	double sin_dRA = sin_ad*sin( PI_2_VALUE-beta )/cos(Dec);

	double dRA = asin( sin_dRA );
	
	RA = RA_obs + dRA;
	
	//but Right Ascension cannot be less then zero and greater then 24h, so...
   AstroAngle::shiftToRangeRAType( RA );

   //and for Declination 90 or -90 degrees there is a default Right Ascension of 0h 0m 0s
	AstroAngle::setDefaultForPole(Dec, ANGLE_DEC_TYPE, RA, ANGLE_RA_TYPE );
}


double gDelta=0.00;
double gAlfa=0.00;
double gBeta=0.00;
double delta_func( double d_obs )
{
	double sin_delta = sin(gDelta);
	double minus = (sin(d_obs)*cos(gAlfa)+sin(gAlfa)*cos(d_obs)*sin(gBeta));
	double ret = (sin_delta-minus);
	return ret;
}

int AstroCCD::CalcCenterFromStar( double ra, double dec, double azim, double alt,
				   time_t ut_time, double x, double y, 
				   double& ra_obs, double& dec_obs, 
					double& azim_obs, double& alt_obs )
{
double rx, ry, distance;
double frame_rotation_angle;
	int ret=0;
	
	//distance from center of image :
	rx = ( x - x_center );
	ry = ( y - y_center );

	// TODO :
	ry = ccd_Y_orient*ry;
	rx = ccd_X_orient*rx;

	distance = sqrt( rx*rx + ry*ry );

	//angular distance in radians:
	double ad = pix2ang( distance );
	double sin_ad = sin(ad);

	//there are two great cirlces that goe straight across the center
	//of the picture horizontaly and verticaly;
	//we rotate these circles clockwise until the "vertical" circle
	//eats our point of interest
	double beta = atan2(ry, rx) + orientation;

	// calculating of RA_obs :
	double sin_delta = sin_ad*cos(beta)/cos(dec);
	double dRA_obs = asin( sin_delta );
	ra_obs = ra - dRA_obs;

	// calculationb of Dec_obs 
	gDelta = dec;
	gAlfa = ad;
	gBeta = beta;	
	double zero_place=0.00;
	if(CMyMathFunc::find_zero_place( delta_func, zero_place, 0.00, PI_2_VALUE )){
		ret=1;
		dec_obs = zero_place;		
	}

	calculateHorizontalCoordinatesFromEq( dec_obs, ra_obs, ut_time, alt_obs, azim_obs );

	return ret;
}


void AstroCCD::calculatePointEquatorialTest( double x, double y,
               			                    double& RA, double& Dec, double& HA,
                        			           double& alt, double& azim )
{
	calculatePointEquatorialTest( x, y, RA, Dec );

	calculateHourAngleBase( RA , Dec, SID_time, HA );
	calculateHorizontalCoordinatesBase( Dec, HA, alt, azim );
}


void AstroCCD::getHorizontalCoo( double ra, double dec,
										time_t ut_time, double geo_long, double geo_lat,
										double& alt, double& azim  )
{
	AstroCCD calc;
	// calc.setRightAscension( AstroAngle( ANGLE_0_2PI, ra ) );
	// calc.setDeclination( AstroAngle( ANGLE_0_2PI, dec ) );	
	
	calc.UT_time = ut_time;
	calc.Lat = geo_lat;
	calc.Long = geo_long;
	calc.calculateHorizontalCoordinatesFromEq( dec, ra, alt, azim );
}


void AstroCCD::calculateHorizontalCoordinatesFromEq( double dec, double ra,
															  		  double& alt, double& azim ) 
{
	double ha;
	calculateHourAngleBase( ra, dec, SID_time, ha );
	calculateHorizontalCoordinatesBase( dec, ha, alt, azim );
}


void AstroCCD::calculateHorizontalCoordinatesFromEq( double dec, double ra, time_t ut_time,
															  		  double& alt, double& azim ) 
{
	double ha;
	calculateHourAngle( ra, dec, ut_time, Long, ha );	
	calculateHorizontalCoordinatesBase( dec, ha, alt, azim );
}

void AstroCCD::calculateHorizontalCoordinatesFromEq( double dec, double ra, time_t ut_time, 
															double& geo_long, double& geo_lat,
													 		double& alt, double& azim )
{
	double ha;
	calculateHourAngle( ra, dec, ut_time, geo_long, ha );
	calculateHorizontalCoordinatesBaseSlow( dec, ha, geo_lat , alt, azim );
}


void AstroCCD::calculateHorizontalCoordinatesBase( double dec, double ha,
															  		double& alt, double& azim )
{
	AstroAngle::cutToRange( dec, ANGLE_DEC_TYPE );

	if( Lat > PI_2_VALUE )
	{
	  alt = dec;
	  azim = ha;
	}else{
		if( Lat < -PI_2_VALUE ){
		  alt = -dec;
		  azim = ha;
		}else{
			// Altitude is...
			alt = asin ( Sin_Lat * sin(dec)
                   + Cos_Lat * cos(dec) * cos(ha) );

			// Azimuth is...
			// azim = atan2 ( cos(dec) * sin(ha)*cos(alt) * cos(Lat),
         //       ( sin(alt) * sin(Lat) - sin(dec) )*cos(alt) );
			double up = cos(dec) * sin(ha);
			double bot1 = Sin_Lat*cos(dec)*cos(ha);
			double bot2 = Cos_Lat*sin(dec);
			
			azim = atan2 ( up , ( bot1-bot2 ) );
		 }
	}

	AstroAngle::cutToRange( azim, ANGLE_AZIM1_TYPE );

	//and for Altitude 90 or -90 degrees there is a default Azimuth of 0d 0m 0s
	AstroAngle::setDefaultForPole( alt, ANGLE_ALT_TYPE, azim, ANGLE_AZIM1_TYPE );
}

void AstroCCD::calculateHorizontalCoordinatesBaseSlow( double dec, double ha, double lat,
															  		double& alt, double& azim )
{
	AstroAngle::cutToRange( dec, ANGLE_DEC_TYPE );

	if( lat > PI_2_VALUE )
	{
	  alt = dec;
	  azim = ha;
	}else{
		if( lat < -PI_2_VALUE ){
		  alt = -dec;
		  azim = ha;
		}else{
			// Altitude is...
			alt = asin ( sin(lat) * sin(dec)
                   + cos(lat) * cos(dec) * cos(ha) );

			// Azimuth is...
			// azim = atan2 ( cos(dec) * sin(ha)*cos(alt) * cos(Lat),
         //       ( sin(alt) * sin(Lat) - sin(dec) )*cos(alt) );
			double up = cos(dec) * sin(ha);
			double bot1 = sin(lat)*cos(dec)*cos(ha);
			double bot2 = cos(lat)*sin(dec);
			
			azim = atan2 ( up , ( bot1-bot2 ) );
		 }
	}

	AstroAngle::cutToRange( azim, ANGLE_AZIM1_TYPE );

	//and for Altitude 90 or -90 degrees there is a default Azimuth of 0d 0m 0s
	AstroAngle::setDefaultForPole( alt, ANGLE_ALT_TYPE, azim, ANGLE_AZIM1_TYPE );
}


void AstroCCD::calculateRightAscension( double ha, double dec, time_t ut_time, double geo_long, double& ra )
{
	ra = getSiderealTimeLocal( ut_time, geo_long ) - ha;
	// ra = getSiderealTimeAtGreenwich( ut_time ) - ha;

	// but Right Ascension cannot be less then zero and greater then 24h, so...
	AstroAngle::shiftToRange( ra, ANGLE_RA_TYPE );

	//and if Declination is 90 or -90 degrees there is a default Right Ascension of 0h 0m 0s
	AstroAngle::setDefaultForPole( dec , ANGLE_DEC_TYPE, ra, ANGLE_RA_TYPE );
}

/*void AstroCCD::calculatePointRightAscension( double HA, double Dec, double& RA )
{
	calculateRightAscension( HA, DEC, RA );

	RA = getSiderealTimeLocal( UT_time, Long ) - HA;
	// RA = getSiderealTimeAtGreenwich( UT_time ) - HA;

	//and if Declination is 90 or -90 degrees there is a default Right Ascension of 0h 0m 0s
	
	AstroAngle::setDefaultForPole( Dec, ANGLE_DEC_TYPE, RA, ANGLE_RA_TYPE );
}*/


void AstroCCD::calculateHourAngleBase( double ra, double dec, double& sid_time, double& ha )
{
	ha = sid_time - ra;

	// but Hour Angle cannot be less then zero and greater then 24h, so...
	AstroAngle::shiftToRange( ha,  ANGLE_HA_TYPE );

	//and if Declination is 90 or -90 degrees there is a default Hour Angle of 0h 0m 0s
	AstroAngle::setDefaultForPole( dec, ANGLE_DEC_TYPE, ha, ANGLE_HA_TYPE );
}




void AstroCCD::calculateEquatorialOrientation()
{
	AstroAngle::cutToRange( Dec_obs , ANGLE_DEC_TYPE );

	if ( ( Lat > PI_2_VALUE ) || ( Lat < -PI_2_VALUE ) )
	{
	  eq_orient = orientation;
	}else{
		eq_orient = orientation
		+ atan2( Cos_Lat*Sin_Azim_obs*Cos_Dec_obs*Cos_Alt_obs,
		         (Sin_Lat - Sin_Dec_obs*Sin_Alt_obs) * Cos_Dec_obs );
	}
}




/*----------------------------   GET FUNCTIONS  -------------------------------------*/

double AstroCCD::JD_to_TJD( double jd )
{
	double tjd = (int)(jd - 2440000.5);
   return tjd;
}

double AstroCCD::getTJD( time_t ut_time )
{
	double jd = getJulianDay( ut_time );
	double tjd = (int)(jd - 2440000.5);
	return tjd;
}

/*  tjd_to_unix is the offset in days between TJD and Unix time */      //kn
static const double tjd_to_unix = -587.0;                               //kn
static const long seconds_per_day = 60*60*24;                           //kn

time_t AstroCCD::unixtime_from_tjdsod( int tjd, double sod )
{
	time_t unixtime = (time_t) (seconds_per_day * (tjd + tjd_to_unix) + (time_t)sod );
	return unixtime;
}

double AstroCCD::jd2hjd( double jd, double ra, double dec )
{
	return ::jd2hjd( jd, ra, dec );
}

double AstroCCD::getHJD( time_t ut_time, double ra, double dec ){
	int y,m,d,h,mi,s;
   get_ymd_hms_ut( ut_time, y,m,d,h,mi,s );

	double ut = h + ((double)mi)/60.0 + ((double)s)/3600.00;
	double hjd = heljuldate( y-2000,m,d, ut,  ra, dec );

	return hjd;
}

time_t AstroCCD::jd2ux( double jd )
{
	int p  = (int)(jd + 0.5);
	int s1 = p + 68569;
   int n  = (int)(4*s1/146097);
   int s2 = s1 - (int)((146097*n + 3)/4);
	int i  = (int)(4000*(s2 + 1)/1461001);
	int s3 = s2 - (int)(1461*i/4) + 31;
   int q  = (int)(80*s3/2447);
   int e  = s3 - (int)(2447*q/80);
   int s4 = (int)(q/11);


	int mm = q + 2 - 12*s4;
   int yy = 100*(n - 49) + i + s4;
   float dd = e + jd - p + 0.5;            // floating point

	// get utc
   float tm = dd;

   dd = (int)(tm);
	tm = 24*(tm - dd);
   int hrs = (int)(tm);
	tm = 60*(tm - hrs);
   int min = (int)(tm);
	tm = 60*(tm - min);
//	int sec = my_round(tm);
	int sec = (int)(tm);


	char szUTC[512];

	sprintf(szUTC,"%d%.2d%.2d_%.2d%.2d%.2d",yy,mm,(int)dd,hrs,min,sec);
	time_t ret = get_gmtime_from_string( szUTC );

	return ret;	
}

double AstroCCD::getJulianDay( time_t ut_time )
{
	struct tm gmtm;
	gmtime_r(&ut_time,&gmtm);	
	int year = gmtm.tm_year+1900;
	int month = gmtm.tm_mon + 1;

	double h = ( gmtm.tm_hour + double(gmtm.tm_min)/60.00 + double(gmtm.tm_sec)/3600.00 );
	double ret = getJulianDay( year, month, gmtm.tm_mday,h );
	return ret;
}

double AstroCCD::getJulianDay( int year, int month, int day, double hour )			//Reference: "Astronomical Algorithms" by Jean Meeus
{
	//valid from 1900/3/1 to 2100/2/28
	if (month<=2) {month=month+12; year=year-1;}
	return (int)(365.25*year) + (int)(30.6001*(month+1)) - 15 + 1720996.5 + day + hour/24.0;
}

double AstroCCD::getJulianDay1( int year, int month, int day, double hour )			//Reference: "Astronomical Algorithms" by Jean Meeus
{
	//valid whenever
	int temp = 2 - year/100 + year/400;
	return int( 365.25 * (year + 4716) ) + int( 30.6001 * (month + 1) )
	       + day + temp - 1524.5 + hour/24.0;
}


double AstroCCD::getSiderealTimeLocal( time_t ut_time, double geo_long )
{
	double gm_sid = getSiderealTimeAtGreenwich( ut_time );	
	double local_sid = gm_sid +  geo_long; // geo-long already in radians !!!
	return local_sid;
}

double AstroCCD::getSiderealTimeAtGreenwich( time_t ut_time )
{
	struct tm _tm;
	gmtime_r( &ut_time, &_tm );	
	double d_hour = AstroAngle::time2hours( _tm.tm_hour, _tm.tm_min, _tm.tm_sec );

	double JD = getJulianDay( _tm.tm_year+1900, _tm.tm_mon+1, _tm.tm_mday, d_hour );
	double T = (JD - 2451545.0)/36525.0;
	double deg_angle =  280.46061837 + 360.98564736629*(JD-2451545.0) + 0.000387933*T*T - T*T*T/38710000.0;
	double in_rad = (deg_angle*PI_VALUE/180.00);

	// printf("test = %.4f\n",in_rad);
	AstroAngle::shiftToRange( in_rad, ANGLE_TIME_TYPE );

	return in_rad;
}

double AstroCCD::getSiderealTimeAtGreenwichFromJD( double JD_date )
{
	double T = (JD_date - 2451545.0)/36525.0;
	double deg_angle =  280.46061837 + 360.98564736629*(JD_date-2451545.0) + 0.000387933*T*T - T*T*T/38710000.0;
	double in_rad = (deg_angle*PI_VALUE/180.00);
	AstroAngle::shiftToRange( in_rad, ANGLE_TIME_TYPE );

	return in_rad;
}



/******************************************************************************************/
/*******************************   POINT FUNCTIONS   **************************************/
/******************************************************************************************/



void AstroCCD::XYAfterTime (double x_start, double y_start, 
									 double& end_x, double& end_y,
									 double sec) // same as above, but you can set X and Y at start
{
	double ha,dec,ra;
	calculatePointEquatorial( x_start, y_start, ra, dec  );
	UT_time += (int)sec;
	calculateRightAscension( ha, dec, ra );
	calculatePointXY( ra, dec, end_x, end_y );
}



/**************************** X,Y to RA,DEC transformation ********************************/




void AstroCCD::calculatePointEquatorial( double x, double y, 
													  double& RA, double& Dec )
{

double rx, ry, distance;
double angular_distance, angle_with_local_meridian, frame_rotation_angle;
double ad, awlm;		//shortcuts


	//coordinates of the center of an image :

	//distance from center of image :
	rx = ( x - x_center );
	ry = ( y - y_center );

	distance = sqrt( rx*rx + ry*ry );

	//angular distance in radians:
	angular_distance = pix2ang( distance );


	//there are two great cirlces that goe straight across the center
	//of the picture horizontaly and verticaly;
	//we rotate these circles clockwise until the "vertical" circle
	//eats our point of interest
	frame_rotation_angle = atan2(rx, ry);


	//this angle together with 'orientation' (if differ from 0)
	//makes an angle between "local meridian" and  our "vertical" circle
	angle_with_local_meridian = frame_rotation_angle + eq_orient;


	//shortcuts instead of these looong names
	ad = angular_distance;
	awlm = angle_with_local_meridian;

	//there was a problem with orientation while Dec_obs was a bit more than 90 degrees, so...
	AstroAngle::cutToRange( Dec_obs, ANGLE_DEC_TYPE ); //not necessary, problem was while 'asin()' was used (atan2() is used now)

	//finally, Declination and Right Ascension are as follows...

	Dec = asin ( cos(ad)*Sin_Dec_obs + sin(ad)*Cos_Dec_obs*cos(awlm) );
	RA = RA_obs - atan2 ( sin(ad)*sin(awlm)*cos(Dec)*Cos_Dec_obs,
	                       ( cos(ad) - sin(Dec)*Sin_Dec_obs ) *  cos(Dec) );

	//but Right Ascension cannot be less then zero and greater then 24h, so...
	AstroAngle::shiftToRange( RA, ANGLE_RA_TYPE );

	//and for Declination 90 or -90 degrees there is a default Right Ascension of 0h 0m 0s
	AstroAngle::setDefaultForPole(Dec, ANGLE_DEC_TYPE, RA, ANGLE_RA_TYPE );
}

void AstroCCD::calculatePointEquatorialNEW( double x, double y, 
													  double& RA, double& Dec )
{

double rx, ry, distance;
double angular_distance, angle_with_local_meridian, frame_rotation_angle;
double ad, awlm;		//shortcuts


	//coordinates of the center of an image :

	//distance from center of image :
	// double eta = in_focal( x - x_center );
	// double teta = in_focal( y - y_center );
	double teta = 0.01772*x + 0.00008*y -2.56216;
	double eta = 0.00377*x + -0.00031*y -0.59796;
	printf("eta=%.5f test=%.5f\n",eta,teta);

	double tg_a_min_RAobs = ( teta/(Cos_Dec_obs-eta*Sin_Dec_obs) );
	RA = RA_obs+atan2( teta, (Cos_Dec_obs-eta*Sin_Dec_obs) );

	
	double tg_dec = ((eta*Cos_Dec_obs+Sin_Dec_obs)*sin(RA-RA_obs))/eta;
	Dec = atan2( ((eta*Cos_Dec_obs+Sin_Dec_obs)*sin(RA-RA_obs)), teta );
}


/*void AstroCCD::XYAfterTimeNEW(double x_start, double y_start, double& end_x, double& end_y, double sec)
{
	double RA,Dec,HA,alt,azim;
	calculatePointEquatorial2( x_start, y_start, RA, Dec, HA, alt, azim );
	
	// changing HA :
	HA += sec*PI_VALUE_TO_HOUR_ANGLE;

	// get new value of RA :
	// calculateRightAscension( HA, Dec, RA );

	// now get new horizontal coordinates :
	double alt_new,azim_new;
	calculateHorizontalCoordinatesBase( Dec, HA,  alt_new, azim_new );

	// now get new position on CCD matrix :
	calculatePointFromHorizonatal( alt_new, azim_new, end_x, end_y );		
}*/


void AstroCCD::XYAfterTimeTEST(double x_start, double y_start, double& end_x, double& end_y, double sec)
{
	double RA,Dec;
	calculatePointEquatorialTest( x_start, y_start, RA, Dec );	

	// changing HA :
   // HA += sec*PI_VALUE_TO_HOUR_ANGLE;

	// no changing of HA instead we change RA :
	RA -= sec*PI_VALUE_TO_HOUR_ANGLE;

	// and now calculate position of "rotated" star on x,y plane :
	calculatePointFromEquatorial( Dec, RA, end_x, end_y );
}

void AstroCCD::XYAfterTimeTEST( double x, double y, int* prevX, int* prevY, int backStart, int backTo,
               			        const double* PrevFramesTime,int time_sign )
{
	double RA,Dec;
	calculatePointEquatorialTest( x, y, RA, Dec );

	for(register int i=backStart;i<=backTo;i++){
		double end_x,end_y;

		// just to skip one call : XYAfterTimeFromEq( Dec, RA, end_x, end_y, PrevFramesTime[i]*time_sign );
		XYAfterTimeFromEq( Dec, RA, end_x, end_y, PrevFramesTime[i]*time_sign );

		// double sec = PrevFramesTime[i]*time_sign;
		// double ra_tmp = (RA - sec*PI_VALUE_TO_HOUR_ANGLE);
		// calculatePointFromEquatorial( Dec, ra_tmp, end_x, end_y );


		prevX[i] = my_round( end_x );
		prevY[i] = my_round( end_y );
	}
}


void AstroCCD::XYAfterTimeFromEq( double dec, double ra, double& end_x, double& end_y, double sec)
{
	ra -= sec*PI_VALUE_TO_HOUR_ANGLE;

	// and now calculate position of "rotated" star on x,y plane :
   calculatePointFromEquatorial( dec, ra, end_x, end_y );
}


void AstroCCD::calculatePointFromHorizonatal( double h, double az, double& x, double& y )
{
	// mozna sprobowac z dA i dh tylko skorzystac i wyliczyc dx , dy ?

	double cos_a = sin(h)*Sin_Alt_obs + cos(h)*Cos_Alt_obs*cos(az-Azim_obs);	
	double a = acos( cos_a );
	double cos_b = ( sin( az-Azim_obs )*cos(h) )/sin(a);
	double beta = acos( cos_b );
	if( h-Alt_obs < 0 )	
		beta = TWO_PI_VALUE - beta;

	double tg_beta = tan(beta);
	double rx = fabs( (focus/pixel_size*tan(a) )/(sqrt(1+tg_beta*tg_beta)) );
	double ry = fabs( rx*tg_beta );
	
	if(beta>PI_2_VALUE && beta<PI_VALUE){
		rx = - rx;
	}else{
		if(beta>PI_VALUE && beta<1.5*PI_VALUE){
			rx = -rx;
			ry = -ry;
		}else{
			if(beta>1.5*PI_VALUE && beta<TWO_PI_VALUE)
				ry = -ry;
		}	
	}

	if(ccd_orient<0){
		rx = -rx;
		ry = -ry;
	}
	/*if( orientation!=0 ){
		rx = ( rx*cos(orientation) - ry*sin(orientation) );
		ry = ( rx*sin(orientation) + ry*cos(orientation) );
	}*/


	// double rot_angle=0.38831872;
	// rx = ( rx*cos(rot_angle) - ry*sin(rot_angle) );
	// ry = ( rx*sin(rot_angle) + ry*cos(rot_angle) );


	x = ( rx + x_center );
   y = ( ry + y_center );
}

void AstroCCD::calculatePointHorizonatal( double x, double y, 
													  double& h, double& az )
{
double rx, ry, distance;
double frame_rotation_angle;
	
	// ccd_orient = -1;
	//coordinates of the center of an image :

	//distance from center of image :
	rx = ( x - x_center );
	ry = ( y - y_center );
	if( ccd_orient<0 ){
		rx = -rx;
		ry = -ry;
	}
	/*if( orientation!=0 ){
		rx = ( rx*cos(orientation) - ry*sin(orientation) );
   	ry = ( rx*sin(orientation) + ry*cos(orientation) );
	}*/

	distance = sqrt( rx*rx + ry*ry );

	//angular distance in radians:
	double ad = pix2ang( distance );


	//there are two great cirlces that goe straight across the center
	//of the picture horizontaly and verticaly;
	//we rotate these circles clockwise until the "vertical" circle
	//eats our point of interest
	frame_rotation_angle = atan2(rx, ry);


	//this angle together with 'orientation' (if differ from 0)
	//makes an angle between "local meridian" and  our "vertical" circle
	// double awlm = frame_rotation_angle + eq_orient;
	double awlm = frame_rotation_angle;


	//there was a problem with orientation while Dec_obs was a bit more than 90 degrees, so...
	// AstroAngle::cutToRange( Dec_obs, ANGLE_DEC_TYPE ); //not necessary, problem was while 'asin()' was used (atan2() is used now)

	//finally, Declination and Right Ascension are as follows...

	// h = asin ( cos(ad)*sin(Alt_obs) + sin(ad)*cos(Alt_obs)*cos(awlm) );	
	// az  = (Azim_obs + atan2( (rx*pixel_size), focus ) );
	// RA = RA_obs - atan2 ( sin(ad)*sin(awlm)*cos(Dec)*cos(Dec_obs),
	//                       ( cos(ad) - sin(Dec)*sin(Dec_obs) ) *  cos(Dec) );

	double beta = atan2( ry, rx ) + orientation;
	h = asin ( cos(ad)*sin(Alt_obs) + sin(ad)*cos(Alt_obs)*sin(beta) );	


 	double dA = asin( ( sin(ad)*cos(beta) )/cos(h) );
	az = Azim_obs + dA;
	

	//but Right Ascension cannot be less then zero and greater then 24h, so...
	//AstroAngle::shiftToRange( h, ANGLE_ALT_TYPE );

	//and for Declination 90 or -90 degrees there is a default Right Ascension of 0h 0m 0s
	AstroAngle::setDefaultForPole(h, ANGLE_ALT_TYPE, az, ANGLE_AZIM1_TYPE );
	AstroAngle::cutToRange( az, ANGLE_AZIM1_TYPE );
}


/**************************** RA,DEC to X,Y transformation ********************************/


void AstroCCD::calculatePointXY( double RA, double Dec, double& x, double& y )
{

double rx, ry, distance;
double angular_distance, angle_with_local_meridian, frame_rotation_angle;
double dRA, ad;		// (dRA = RA - RA_obs) and (ad - shortcut to angular distance)


	//coordinates of the center of an image


	//there was a problem with orientation while Dec_obs was a bit more than 90 degrees, so...
	AstroAngle::cutToRange( Dec_obs, ANGLE_DEC_TYPE ); //not necessary, problem was while 'asin()' was used (atan2() is used now)

	AstroAngle::cutToRange( Dec, ANGLE_DEC_TYPE );


	//first step

	dRA = RA_obs - RA;	//difference between RA of the center of the picture and RA of an object
	angular_distance = acos ( sin(Dec)*Sin_Dec_obs
	                          + cos(Dec)*Cos_Dec_obs*cos(dRA) );

	ad = angular_distance;
	angle_with_local_meridian = atan2 ( cos(Dec)*sin(dRA)*sin(ad)*Cos_Dec_obs,
	                                    ( sin(Dec) - cos(ad)*Sin_Dec_obs ) * sin(ad) );

	//second step

	frame_rotation_angle = angle_with_local_meridian - eq_orient;

	distance = ang2pix( angular_distance );


	//third step

	rx = distance*sin(frame_rotation_angle);
	ry = distance*cos(frame_rotation_angle);


	//final step

	x = ( rx + x_center );
	y = ( ry + y_center );

}

/**************************** RA(HA),Dec to Alt,Azim transformation ********************************/

/*void AstroCCD::calculatePointHorizontalCoordinates( double Dec, double HA,
																	 double& Alt, double& Azim)
{
 AstroAngle::cutToRange( Dec, ANGLE_DEC_TYPE );

 if ( Lat > PI_2_VALUE )
	{
	  Alt = Dec;
	  Azim = HA;
	}
 else
 if ( Lat < -PI_2_VALUE )
	{
	  Alt = - Dec;
	  Azim =  HA;
	}
 else
 {
	// Altitude is...
	Alt = asin ( Sin_Lat*sin(Dec) + Cos_Lat*cos(Dec)*cos(HA) );

	// Azimuth is...
	Azim = atan2 ( cos(Dec)*sin(HA)*cos(Alt)*Cos_Lat,
	               ( sin(Alt)*Sin_Lat - sin(Dec) )*cos(Alt) );
 }

	//and for Altitude 90 or -90 degrees there is a default Azimuth of 0d 0m 0s
	AstroAngle::setDefaultForPole( Alt, ANGLE_ALT_TYPE, Azim, ANGLE_AZIM1_TYPE );

}*/



/*void AstroCCD::calculatePointHourAngle( double RA, double Dec, double& HA )
{
	HA = getSiderealTimeLocal( UT_time, Long ) - RA;
	// HA = getSiderealTimeAtGreenwich(UT_time) - RA;

	//and if Declination is 90 or -90 degrees there is a default Hour Angle of 0h 0m 0s
	
	AstroAngle::setDefaultForPole( Dec, ANGLE_DEC_TYPE, HA, ANGLE_HA_TYPE );
}*/

void AstroCCD::correctAltitudeForRefraction( double& Alt )
{
	double zd = getPointZenithalDistance( Alt );

	double correction = (1.0 - 1.29125e-8 * pow(zd,3.5) ) * tan(zd) * 59.6;
	double correction_in_rad = AstroAngle::arcsec2rad( correction );
	Alt = Alt - correction_in_rad;
	
	// Alt.setInArcSeconds( Alt.inArcSec() - correction );
}


/*void AstroCCD::correctForRefraction()
{
	calculateHourAngle();
	calculatePointHorizontalCoordinates();
	correctAltitudeForRefraction();
	calculatePointEquatorialCoordinates ();
	calculatePointRightAscension ();
}*/


/*-------------------------   angle <-> pixels transformations  --------------------*/



double AstroCCD::calcStandardCoord( double ra, double dec, double& eta, double& teta )
{
	teta = sin(ra-RA_obs)/( sin(Dec_obs)*tan(dec) + cos(Dec_obs)*cos(ra-RA_obs) );
	eta = ( tan(dec) - tan(Dec_obs)*cos(ra-RA_obs) )/(tan(Dec_obs)*tan(dec) +cos(ra-RA_obs) );
	printf("RA_obs=%s Dec_obs=%.2f\n",AstroAngle::toString(RA_obs,ANGLE_RA_TYPE).c_str(),
				AstroAngle::rad2deg(Dec_obs));
	printf("Azim_obs=%.2f Alt_obs=%.2f\n",AstroAngle::rad2deg(Azim_obs),AstroAngle::rad2deg(Alt_obs));
	return 0.00;
}


void AstroCCD::calcHorToEqHour( double h, double az, double phi, double& t, double& dec )
{
	double tan_t = ( cos(h)*sin(az) )/( cos(phi)*sin(h) + sin(phi)*cos(h)*cos(az) );
	double sin_dec   = (sin(phi)*sin(h) - cos(phi)*cos(h)*cos(az) );

	t = atan( tan_t );
	dec = asin( sin_dec );
}


void AstroCCD::calcEqHourToHoriz( double t, double dec, double phi, double& h, double& az )
{
	double tan_A = ( cos(dec)*sin(t) )/( sin(phi)*cos(dec)*cos(t) - cos(phi)*sin(dec) );
	double sin_h = ( sin(phi)*sin(dec) + cos(phi)*cos(dec)*cos(t) );

	az = atan( tan_A );
	h  = asin( sin_h );
}

double AstroCCD::CalcDistRADEC( double rah1, double dec_deg1, double rah2, double dec_deg2 )
{
	return AstroCCD::CalcDistInDeg( rah1*15.0, dec_deg1, rah2*15.00, dec_deg2 );
}

double AstroCCD::CalcDistInDeg( double ra1, double dec1, double ra2, double dec2 )
{
	double ra1_in_rad = AstroAngle::deg2rad( ra1 );
	double ra2_in_rad = AstroAngle::deg2rad( ra2 );
	double dec1_in_rad = AstroAngle::deg2rad( dec1 );
	double dec2_in_rad = AstroAngle::deg2rad( dec2 );

	double cos_x = sin(dec1_in_rad)*sin(dec2_in_rad) + cos(dec1_in_rad)*cos(dec2_in_rad)*cos( ra1_in_rad - ra2_in_rad );
	double x = acos( cos_x );
	return AstroAngle::rad2deg( x );
}

double AstroCCD::CalcDistInRad( double ra1, double dec1, double ra2, double dec2 )
{
	double cos_x = sin(dec1)*sin(dec2) + cos(dec1)*cos(dec2)*cos( ra1 - ra2 );
	double x = acos( cos_x );
	return AstroAngle::rad2deg( x );
}


int AstroCCD::calc_fov_range( double ra0, double dec0, double fov,
                              double& ra_min, double& ra_max,
                              double& dec_min, double& dec_max )
{
	int ret=1;

	double half_fov = ( fov/2.00 );
	double radius_in_deg = ( fov/2.0 )*CMyMathFunc::mysqrt(2.00) + 3.00;
	dec_min = dec0 - half_fov - 3;
	dec_max = dec0 + half_fov + 3;
	double radius_in_h = ( radius_in_deg / 15.00 );

	// value further from equator is put here :
	double dec_eq_futher = dec_min;
	if( fabs(dec_max) > fabs(dec_min) ){
		dec_eq_futher = dec_max;
	}
	
	double ra=ra0;
	double step_ra=( 1.0/60.0 )/15.00; // 1 arcmin 
	int bOK=1;
	double distance = CalcDistRADEC( ra0, dec0, ra, dec_eq_futher );
	while( distance < radius_in_deg && ra > -13.00 ){
		ra = ra - step_ra;
		distance = CalcDistRADEC( ra0, dec0, ra, dec_eq_futher );
	}
	double half_size_ra = fabs( ra0 - ra );
	ra_min = ra;
	if( ra_min < 0 && ra_min>=-12 ){
		ra_min = 24 + ra_min;
	}
	ra_max = ra0 + half_size_ra;	
	
	return ret;
}

int AstroCCD::c__1 = 1;
double AstroCCD::c_b11 = 2e3;
double AstroCCD::c_b13 = 360.;
int AstroCCD::c__5 = 5;
double AstroCCD::c_b35 = 6.2831853071779999;


static double my_d_mod( double* a, double* b )
{
	if( (*b)!=0 ){
		int d = (int)( (*a) / (*b) );
		double ret = (*a) - d*(*b);
		return ret;
	}
	return 0;
}

int AstroCCD::gal2eq(double *gl, double *b, double *raepo,double *depo, double *epoch)
{
    /* Initialized data */

    static double b1950 = 2433282.42345905;
    static double pi = 3.141592653589;
    static int ifirst = 0;
    static double rap = 192.25;
    static double dp = 27.4;
    static double gle = 33.;

    /* Local variables */
    static double gl0, d1950, dpr, ra1950, sind, rapr, pi_180__, djepo, 
	    cosdp, sindp, cdcdra, cdsdra;

/* Converts galactic coordinates (gl,b) given in radians */
/* into the equatorial coordinates (RAepo,Depo) [radians] */
/* at the Julian epoch "Epoch". */
/* Standard galactic pole coordinates and center position at B1950 [degrees]: */
    if (ifirst != 0) {
	goto L1;
    }
    ++ifirst;
    pi_180__ = pi / 180;
    rapr = rap * pi_180__;
    dpr = dp * pi_180__;
    gl0 = gle * pi_180__;
    sindp = sin(dpr);
    cosdp = cos(dpr);
L1:
    sind = sindp * sin(*b) + cosdp * cos(*b) * sin(*gl - gl0);
    d1950 = asin(sind);
    cdsdra = cos(*b) * cos(*gl - gl0) * cosdp;
    cdcdra = sin(*b) - sind * sindp;
    ra1950 = atan2(cdsdra, cdcdra) + rapr;
    djepo = (*epoch - 2e3) * 365.25 + 2451545;
    prenew(&b1950, &djepo, &ra1950, &d1950, raepo, depo);
    return 0;
} /* gal2eq_ */

int AstroCCD::eq2gal2k(double *ra2000, double *d2000, double *gl, double *b)
{
    /* Initialized data */

    static double pi = 3.141592653589;
    static int ifirst = 0;
    static double rap = 192.8593357295;
    static double dp = 27.1282510272;
    static double gl0 = 32.9319193453;

    /* System generated locals */
    double d__1, d__2;

    /* Local variables */
    static double dpr, gl0r, sinb, rapr, pi_180__, cbcdl, cbsdl, cosdp, 
	    sindp;

/* Procedure to convert the equatorial coordinates (ra2000,d2000) [radians] */
/* given for the J2000 epoch to the galactic coordinates (gl,b) [radians] */
/* Here above are galactic pole coordinates and location of the galactic */
/* center precessed to J2000 from standard values at B1950 (in degrees): */
/* RAp=192.25, Dp=27.40 and gl0=33.00. The last parameter, gl0, is */
/* the galactic longitude of the equator. */
    if (ifirst != 0) {
	goto L1;
    }
    ++ifirst;
    pi_180__ = pi / 180;
    rapr = rap * pi_180__;
    dpr = dp * pi_180__;
    sindp = sin(dpr);
    cosdp = cos(dpr);
    gl0r = gl0 * pi / 180;
L1:
    sinb = sindp * sin(*d2000) + cosdp * cos(*d2000) * cos(*ra2000 - rapr);
    *b = asin(sinb);
    cbsdl = sin(*d2000) * cosdp - cos(*d2000) * sindp * cos(*ra2000 - rapr);
    cbcdl = cos(*d2000) * sin(*ra2000 - rapr);
    *gl = atan2(cbsdl, cbcdl) + gl0r;
    d__1 = *gl + pi + pi;
    d__2 = pi + pi;
    *gl = my_d_mod(&d__1, &d__2);
    return 0;
} /* eq2gal2k_ */

int AstroCCD::prenew(double *dje1, double *dje2, double *
	ra1, double *d1, double *ra2, double *d2)
{
    /* Initialized data */

    static double dcsar = 4.848136812e-6;
    static int j2000 = 2451545;
    static double dc1 = 2306.2181;
    static double dc2 = 1.39656;
    static double dc3 = .30188;
    static double dc4 = .017998;
    static double dc5 = 1.09468;
    static double dc6 = 2004.3109;
    static double dc7 = -.8533;
    static double dc8 = -.42665;
    static double dc9 = -.041833;

    /* System generated locals */
    double d__1;

    /* Local variables */
    static double dc, dt, dt0, dtc, dth, dts, dzet, dzeta;

/* Calculates general precession from DJE1 to DJE2 using new IAU theory */
/* Constants according to J.H.Lieske 1979, Astron. Astrophys. 73, 282. */
    dt0 = (*dje1 - j2000) / 36525;
    dt = (*dje2 - *dje1) / 36525;
    dts = dt * dt;
    dtc = dts * dt;
    dzeta = ((dc1 + (dc2 - dt0 * 1.39e-4) * dt0) * dt + (dc3 - dt0 * 3.44e-4) 
	    * dts + dc4 * dtc) * dcsar;
    dzet = dzeta + ((dc5 - dc3 + dt0 * 4.1e-4) * dts + dtc * 2.05e-4f) * 
	    dcsar;
    dth = ((dc6 + (dc7 - dt0 * 2.17e-4) * dt0) * dt + (dc8 - dt0 * 2.17e-4) * 
	    dts + dc9 * dtc) * dcsar;
/* convert Dj1 position to Dj2 using new IAU precession angles */
    dc = cos(*ra1 + dzeta);
    *ra2 = dzet + atan2(sin(*ra1 + dzeta), -tan(*d1) * sin(dth) + cos(dth) * 
	    dc);
    *d2 = asin(sin(*d1) * cos(dth) + cos(*d1) * sin(dth) * dc);
    d__1 = *ra2 + 6.2831853071779999;
    *ra2 = my_d_mod(&d__1, &c_b35);

//	printf("%.8f,%.8f,%.8f,%.8f -> %.8f,%.8f\n",*dje1,*dje2,*ra1,*d1,*ra2,*d2);

    return 0;
} /* prenew */

int AstroCCD::eq2gal0(double *raepo, double *depo, double 
	*epoch, double *gl, double *b)
{
    /* Initialized data */

    static double b1950 = 2433282.42345905;
    static double pi = 3.141592653589;
    static int ifirst = 0;
    static double rap = 192.25;
    static double dp = 27.4;
    static double gl0 = 33.;

    /* System generated locals */
    double d__1, d__2;


    /* Local variables */
    static double d1950, gl0r, ra1950, sinb, rapr, pi_180__, cbcdl, cbsdl,
	     djepo, cosdp, sindp;

/* Converts equatorial coordinates (RAepo,Depo) given in radians */
/* for any Julian epoch "Epoch" (e.g. 2001.526) into the galactic */
/* coordinates (gl,b) in radians. */
/* Standard galactic center and pole coordinates at B1950 [degrees]: */
    djepo = 2451545 + 365.25 * (*epoch - 2e3);
    prenew(&djepo, &b1950, raepo, depo, &ra1950, &d1950);
    if (ifirst != 0) {
	goto L1;
    }
    ++ifirst;
    pi_180__ = pi / 180;
    rapr = rap * pi_180__;
    sindp = sin(dp * pi_180__);
    cosdp = cos(dp * pi_180__);
    gl0r = gl0 * pi_180__;
L1:
    sinb = sindp * sin(d1950) + cosdp * cos(d1950) * cos(ra1950 - rapr);
    *b = asin(sinb);
    cbsdl = sin(d1950) * cosdp - cos(d1950) * sindp * cos(ra1950 - rapr);
    cbcdl = cos(d1950) * sin(ra1950 - rapr);
    *gl = atan2(cbsdl, cbcdl) + gl0r;
    d__1 = *gl + pi + pi;
    d__2 = pi + pi;
    *gl = my_d_mod(&d__1, &d__2);
    return 0;
} /* eq2gal0_ */

int AstroCCD::eq2gal(double *raepo, double *depo, double *
	epoch, double *gl, double *b)
{
    /* Initialized data */

    static double b1950 = 2433282.42345905;
    static double pi = 3.141592653589;
    static int ifirst = 0;
    static double rap = 192.25;
    static double dp = 27.4;
    static double gle = 33.;

    /* System generated locals */
    double d__1, d__2;


    /* Local variables */
    static double gl0, d1950, dpr, ra1950, sinb, rapr, pi_180__, cbcdl, 
	    cbsdl, djepo, cosdp, sindp;

/* Converts equatorial coordinates (RAepo,Depo) given in radians */
/* from any Julian epoch "Epoch" (e.g. 2001.526) into the galactic */
/* coordinates (gl,b) [radians]. */
/* Standard galactic center and pole coordinates at B1950 [degrees]: */
    djepo = 2451545 + 365.25 * (*epoch - 2e3);
    prenew(&djepo, &b1950, raepo, depo, &ra1950, &d1950);
    if (ifirst != 0) {
	goto L1;
    }
    ++ifirst;
    pi_180__ = pi / 180;
    rapr = rap * pi_180__;
    dpr = dp * pi_180__;
/* The angular distance: galactic center to equator */
    gl0 = gle * pi_180__;
    sindp = sin(dpr);
    cosdp = cos(dpr);
L1:
    sinb = sindp * sin(d1950) + cosdp * cos(d1950) * cos(ra1950 - rapr);
    *b = asin(sinb);
    cbsdl = sin(d1950) - sinb * sindp;
    cbcdl = cos(d1950) * sin(ra1950 - rapr) * cosdp;
    *gl = atan2(cbsdl, cbcdl) + gl0;
    d__1 = *gl + pi + pi;
    d__2 = pi + pi;
    *gl = my_d_mod(&d__1, &d__2);
    return 0;
} /* eq2gal_ */

