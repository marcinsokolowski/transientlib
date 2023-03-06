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


#include "Astroangle.h"
#include "Header.h"
#include <stdio.h>
#include "mathfunc.h"


/************************************************************/


AstroAngle::AstroAngle(angle_type type_to_set, double angle_to_set)
{
	setInRadians ( angle_to_set );
	setAngleType ( type_to_set );
} //End of constructor 1


/************************************************************/


AstroAngle::AstroAngle (angle_type type_to_set, double main, double min, double sec)
{
	type = type_to_set;

	setAngle (main, min, sec);

} //End of constructor 2



/************************************************************/


void AstroAngle::setAngle (double main, double min, double sec )	// remember to set correct angle type before setting angle !!!
{

  switch (type)
	{

	// angles that can be positve or negative (negative angle has only one negative component
	// - first, that does not equal zero from left) and are measured in degrees, minutes and seconds
	  case ANGLE_ALT_TYPE:
	  case ANGLE_AZIM2_TYPE:
	  case ANGLE_DEC_TYPE:
	  case ANGLE_LAT_TYPE:
	  case ANGLE_TYPE_1:
	  case ANGLE_TYPE_2:
	  case ANGLE_NOT_RANGED:
	  case ANGLE_LONG2_TYPE:
	  {
		if ( (main < 0.0) || (min < 0.0) || (sec < 0.0) )
			angle = - ( ( abs(main)*60.0 + abs(min) )*60.0 + abs(sec) ) / ARCSEC_IN_RAD;
		else
			angle = ( ( main*60.0 + min )*60.0 + sec ) / ARCSEC_IN_RAD;
		break;
	  }

	// angles that are measured like time (in hours, minutes and seconds) and have no negative values
	  case ANGLE_TIME_TYPE:
	  case ANGLE_HA_TYPE:
	  case ANGLE_RA_TYPE:
	  {
		if ( (main < 0.0) || (min < 0.0) || (sec < 0.0) )
		{
			cout<<"Warning! You have just set negative time. Please make sure you really mean this. "
			    <<"You can also shift it to range. If you set time to -3h after shifting "
			    <<"( your_angle.shiftToRange() ) it will be 21h."<<endl;
			cout<<endl;
			angle = - ( ( abs(main)*60.0 + abs(min) )*60.0 + abs(sec) ) / SEC_IN_RAD;
		}
		else
			angle = ( ( main*60.0 + min)*60.0 + sec ) /  SEC_IN_RAD;
		break;
	  }

	// all other angles that are measured in degrees, minutes and seconds and can not be negative
	  default:
	  {
		if ( (main < 0.0) || (min < 0.0) || (sec < 0.0) )
		{
			cout<<"Warning! You have just set to negative angle that cannot be negative. "
			    <<"Please make sure you really mean this. You can also shift it to range. "
			    <<"If you set -12deg to angle that has range '0 to 360'deg than after "
			    <<"shifting ( your_angle.shiftToRange() ) it will have 348deg."<<endl;
			cout<<endl;
			angle = - ( ( abs(main)*60.0 + abs(min) )*60.0 + abs(sec) ) / ARCSEC_IN_RAD;
		}
		else
			angle = ( ( main*60.0 + min)*60.0 + sec ) / ARCSEC_IN_RAD;
	  }

	}
} //End of setAngle(...)



/************************************************************/


void AstroAngle::setAngleType (angle_type type_to_set)
{
	type = type_to_set;
} //End of setAngleType(...)



/************************************************************/
// conversions :
double AstroAngle::timeangle2rad( double time_angle )
{
	double ret = (time_angle/24.00)*PI_VALUE*2.00;
	return ret;
}

double AstroAngle::timeangle2rad( int hour, int min, double sec  )
{
	double hours = time2hours( hour, min ,sec );
	return (hours/24.00)*TWO_PI_VALUE;
}


double AstroAngle::timeangle2deg( int hour, int min, double sec  )
{
	double hours = time2hours( hour, min ,sec );
	return (hours/24.00)*360.00;
}

void AstroAngle::rad2timeangle( double in_rad, int& hour, int& min, int& sec )
{
	// double in_sec = fabs( in_rad*((3600.000*24.00)/(2.00*PI_VALUE)) );
	// hour = mysign(in_rad)*( ((int)in_sec)/3600 );
	double in_sec =  in_rad*((3600.000*24.00)/(2.00*PI_VALUE));
	hour = ( ((int)in_sec)/3600 );
	min = ( ((int)in_sec)%3600 )/60;
	sec = ( ( ((int)in_sec)%3600 )%60 );
}

void AstroAngle::deg2timeangle( double in_deg, int& hour, int& min, double& sec )
{
	double in_hours = in_deg/15.00;
	hour = (int)in_hours;
	double rest = (fabs(in_hours)-fabs(double(hour)));
	min=(int)(rest*60.000);
   sec = (3600.00*rest-min*60.00);
}

void AstroAngle::rad2timeangle_new( double in_rad, int& hour, int& min, 
												double& sec )
{
	double in_hours = in_rad*(RAD_TO_HOURS);
	hour = (int)in_hours;
	double rest = (fabs(in_hours)-fabs(double(hour)));
	min=(int)(rest*60.000);
	sec = (3600.00*rest-min*60.00);
}

void AstroAngle::rad2timeangle( double in_rad, double& hour, double& min, 
                 	              double& sec )
{
	double in_hours = in_rad*(RAD_TO_HOURS);
	hour = (int)in_hours;
	double rest = (in_hours-hour);
	min=(int)(rest*60);
	sec = (3600*rest-min*60);
}

double AstroAngle::hours2deg( double in_hours )
{
	return ( in_hours * 15.000 );
}

double AstroAngle::hours2rad( double in_hours )
{
	double in_rad = in_hours/(RAD_TO_HOURS);
	return in_rad;
}

double AstroAngle::rad2hours( double in_rad )
{
	double in_hours = (in_rad*RAD_TO_HOURS);
	return in_hours;
}

double AstroAngle::rad2arcsec( double in_rad )
{
	double in_deg = rad2deg( in_rad );
	double ret = in_deg*3600.00;
	return ret;
}

double AstroAngle::degtime2deg( int deg, int min, double sec )
{
	double in_deg = abs(deg) + (abs(min)*60.00 + fabs(sec))/3600.00;
	return in_deg*mysign(deg);
}

double AstroAngle::deg2hours( double in_deg )
{
	double ret = (in_deg/15.00);
	return ret;
}

void AstroAngle::deg2deg( double in_deg, int& deg, int& min, int& sec )
{
	deg = (int)in_deg;

	double partial = ( abs(in_deg) - abs(deg) );
	min = (int)(partial*60.00);
	sec = (int)((partial*60.00 - min)*60.00);
}

void AstroAngle::deg2deg( double in_deg, int& deg, int& min, double& sec )
{
	deg = (int)in_deg;

	double partial = ( abs(in_deg) - abs(deg) );
	min = (int)(partial*60.00);
	sec = (partial*60.00 - min)*60.00;
}



double  AstroAngle::degdeg2deg( int deg, int min, double sec )
{
	double in_deg = abs(deg)+double(min*60+sec)/3600.00;
	if( deg<0 )
		in_deg = -in_deg;
	return in_deg;
}

mystring AstroAngle::rad2degstring( double in_rad )
{
	double in_deg = rad2deg( in_rad );
	int deg,min,sec;
	deg2deg( in_deg, deg, min, sec );
	char szOut[128];
	sprintf(szOut,"%d %d %d",deg,min,sec);

	mystring szRet=szOut;
	return szRet;
}

double AstroAngle::ra2deg( const char* szRA )
{
	int h,m;
   double s;
   sscanf( szRA, "%dh%dm%lfs", &h,&m,&s);

	double in_rad = mysign_non_zero( h )*AstroAngle::timeangle2rad( abs(h) ,m ,s );	
	double in_deg = AstroAngle::rad2deg( in_rad );

	return in_deg;
}


double AstroAngle::deg2rad( double in_deg )
{
	double ret = (in_deg)*PI_VALUE/180.00;
	return ret;
}

double AstroAngle::rad2deg( double in_rad )
{
	double ret = ( in_rad )*(RAD_TO_DEG);
	return ret;
}

double AstroAngle::arcsec2rad( double arcsec )
{
	double ret = ( arcsec  / ARCSEC_IN_RAD );
	return ret;
}

double AstroAngle::time2hours( int hour, int min, double sec )
{
	double ret = hour + ( min*60.00 + sec )/3600.00;
	return ret;
}

void AstroAngle::setInArcDegrees (int deg, int min, double sec)
{
	if ( (deg < 0.0) || (min < 0.0) || (sec < 0.0) )
	{
	    angle = - ( ( abs(deg)*60.0 + abs(min) )*60.0 + abs(sec) ) / ARCSEC_IN_RAD;
	}
	else
	{
	    angle = ( ( deg*60.0 - min )*60.0 - sec ) / ARCSEC_IN_RAD;
	}
}


void AstroAngle::setInArcDegrees (double deg)
{
	angle = ( deg*3600.0 ) / ARCSEC_IN_RAD;
}


void AstroAngle::setInArcMinutes (double min)
{
	angle = ( min*60.0 ) / ARCSEC_IN_RAD;
}


void AstroAngle::setInArcSeconds (double sec)
{
	angle = ( sec ) / ARCSEC_IN_RAD;
}



/************************************************************/


void AstroAngle::setInHours (int hour, int min, double sec)
{
	angle = ( ( hour*60.0 + min)*60.0 + sec ) /  SEC_IN_RAD;
}


void AstroAngle::setInHours (double hour)
{
	angle = ( hour*3600.0 ) /  SEC_IN_RAD;
}


void AstroAngle::setInMinutes (double min)
{
	angle = ( min*60.0 ) /  SEC_IN_RAD;
}


void AstroAngle::setInSeconds (double sec)
{
	angle = ( sec ) /  SEC_IN_RAD;
}



/************************************************************/


void AstroAngle::setInRadians( double rad )
{
	angle = rad;
}



/************************************************************/


string AstroAngle::toString ()
{

	return toString (type);

} //End of giveAngle()

string AstroAngle::toString(double _angle, angle_type in_type)
{
  ostringstream angleString;
  int precision = 4;
  string szSign;
  getSign( _angle, in_type, szSign );	

  if( in_type!=ANGLE_TIME_TYPE && in_type!=ANGLE_RA_TYPE && in_type!=ANGLE_HA_TYPE ){
	  angleString<<szSign;
  }
  switch (in_type)
  	{

	// angles that are measured like time (in hours, minutes and seconds)
	  case ANGLE_TIME_TYPE:
	  case ANGLE_HA_TYPE:
	  case ANGLE_RA_TYPE:
	  {
	  	/*if ( sec( _angle ) < 10.0 )   precision = 3;

		angleString<<setfill('0')<<setw(2)<<hour( _angle)<<'h'
		           <<setfill('0')<<setw(2)<<min( _angle )<<'m';

		if ( sec(_angle) - floor(sec(_angle)) < 0.01 )
			angleString<<setw(2)<<sec(_angle)<<".00s";
		else 
			angleString<<setw(5)<<setprecision(precision)<<sec(_angle)<<'s';*/
		// NEW change on 20041021 due to strange : +00h37m49.43s,15.6323
		char szTmp[64];
		int h,m;
		double s;
		if( _angle<0 ){
			_angle = TWO_PI_VALUE + _angle;
		}
		// printf("_angle = %.2f\n",_angle);
		rad2timeangle_new( _angle, h, m, s );
		sprintf(szTmp,"%.2dh%.2dm%#05.2fs",h,m,s);
		// printf("szSign = %s\n",szSign.c_str());
		// printf("h = %.2d\n",h);
		// printf("szTmp = %s %d %d %d\n",szTmp,h,m,s);
		angleString << szTmp;

		break;
	  }

	// angles that are measured in degrees, minutes and seconds and cannot have a 3 digit value
	  case ANGLE_TYPE_1:
	  case ANGLE_DEC_TYPE:
	  case ANGLE_LAT_TYPE:
	  case ANGLE_ALT_TYPE:
	  {
	  /*	if ( arcsec(_angle) < 10.0 )   precision = 3;

		angleString<<setfill('0')<<setw(2)<<arcdeg(_angle)<<'d'
		           <<setfill('0')<<setw(2)<<arcmin(_angle)<<'m';

		if ( arcsec(_angle) - floor(arcsec(_angle)) < 0.01 )
			angleString<<setw(2)<<arcsec(_angle)<<".00s";
		else angleString<<setw(5)<<setprecision(precision)<<arcsec(_angle)<<'s';
		*/

		char szTmp[64];
		int deg,min;
		double s;
	
		double in_deg = rad2deg( _angle );
		deg2deg( in_deg, deg, min, s );
		sprintf(szTmp,"%.2dd%.2dm%#05.2fs",deg,min,s);
		angleString << szTmp;

		break;
	  }
	  case ANGLE_DEG_TYPE:
		{
			char tmp[20];
			sprintf(tmp,"%.2f deg",rad2deg(_angle));
			angleString<<tmp;
			break;
		}	

	// all other angles that are measured in degrees, minutes and seconds and can have 3 digits
	  default:
	  {
	  	if ( arcsec(_angle) < 10.0 )   precision = 3;

		angleString<<setfill(' ')<<setw(3)<<arcdeg(_angle)<<'d'
		           <<setfill('0')<<setw(2)<<arcmin(_angle)<<'m';

		if ( arcsec(_angle) - floor(arcsec(_angle)) < 0.01 )
			angleString<<setw(2)<<arcsec(_angle)<<".00s";
		else	angleString<<setw(5)<<setprecision(precision)<<arcsec(_angle)<<'s';

		break;
	  }
		
	}

  return angleString.str();


}

string AstroAngle::toString( angle_type in_type)
{

  ostringstream angleString;
  int precision = 4;


  angleString<<getSign()<<' ';


  switch (in_type)
  	{

	// angles that are measured like time (in hours, minutes and seconds)
	  case ANGLE_TIME_TYPE:
	  case ANGLE_HA_TYPE:
	  case ANGLE_RA_TYPE:
	  {
	  	if ( sec() < 10.0 )   precision = 3;

		angleString<<setfill('0')<<setw(2)<<hour()<<'h'
		           <<setfill('0')<<setw(2)<<min()<<'m';

		if ( sec() - floor(sec()) < 0.01 )
			angleString<<setw(2)<<sec()<<".00s";
		else	angleString<<setw(5)<<setprecision(precision)<<sec()<<'s';

		break;
	  }

	// angles that are measured in degrees, minutes and seconds and cannot have a 3 digit value
	  case ANGLE_TYPE_1:
	  case ANGLE_DEC_TYPE:
	  case ANGLE_LAT_TYPE:
	  case ANGLE_ALT_TYPE:
	  {
	  	if ( arcsec() < 10.0 )   precision = 3;

		angleString<<setfill('0')<<setw(2)<<arcdeg()<<'d'
		           <<setfill('0')<<setw(2)<<arcmin()<<'m';

		if ( arcsec() - floor(arcsec()) < 0.01 )
			angleString<<setw(2)<<arcsec()<<".00s";
		else	angleString<<setw(5)<<setprecision(precision)<<arcsec()<<'s';

		break;
	  }

	// all other angles that are measured in degrees, minutes and seconds and can have 3 digits
	  default:
	  {
	  	if ( arcsec() < 10.0 )   precision = 3;

		angleString<<setfill(' ')<<setw(3)<<arcdeg()<<'d'
		           <<setfill('0')<<setw(2)<<arcmin()<<'m';

		if ( arcsec() - floor(arcsec()) < 0.01 )
			angleString<<setw(2)<<arcsec()<<".00s";
		else	angleString<<setw(5)<<setprecision(precision)<<arcsec()<<'s';

		break;
	  }

	}

  return angleString.str();

} //End of giveAngleInType(...)



/************************************************************/


double AstroAngle::inArcDeg ()
{
	return ( angle * ARCSEC_IN_RAD ) / 3600.0;
}


double AstroAngle::inArcMin ()
{
	return ( angle * ARCSEC_IN_RAD ) / 60.0;
}


double AstroAngle::inArcSec ()
{
	return ( angle * ARCSEC_IN_RAD );
}



/************************************************************/


double AstroAngle::inHours ()
{
	return ( angle * SEC_IN_RAD ) / 3600.0;
}


double AstroAngle::inMin ()
{
	return ( angle * SEC_IN_RAD ) / 60.0;
}


double AstroAngle::inSec ()
{
	return ( angle * SEC_IN_RAD );
}



/************************************************************/


double AstroAngle::inRad ()
{
	return angle;
}



/************************************************************/


int AstroAngle::arcdeg( double _angle )
{
	double abs_arcsec = abs ( _angle * ARCSEC_IN_RAD );

	int arcdeg = int ( abs_arcsec / 3600.0 );

	if ( fmod( abs_arcsec , 3600.0) > 3599.9949 ) arcdeg++;

	return arcdeg;
}


int AstroAngle::arcmin( double _angle )
{
	double abs_arcsec = abs ( _angle * ARCSEC_IN_RAD );

	int arcmin = int ( fmod( abs_arcsec ,3600.0) / 60.0 );

	if ( fmod( abs_arcsec , 60.0) > 59.9949 ) arcmin++;

	if ( arcmin > 59 ) arcmin = 0;

	return arcmin;
}


double AstroAngle::arcsec( double _angle)
{
	double abs_arcsec = abs ( _angle * ARCSEC_IN_RAD );
	double arcsec = floor( fmod( abs_arcsec ,60.0) * 100.0 + 0.5 ) / 100.0;

	if ( arcsec > 59.98 ) arcsec = 0.0;		// ??? dlaczego nie dziala 59.99 ???
	return arcsec;
}



/************************************************************/


int AstroAngle::hour( double _angle )
{
	double abs_sec = abs ( _angle * SEC_IN_RAD );

	int hour = int ( abs_sec / 3600.0 );

	if ( fmod( abs_sec , 3600.0) > 3599.9949 ) hour++;

	return hour;

}

int AstroAngle::min( double _angle )
{
	double abs_sec = abs ( _angle * SEC_IN_RAD );

	int min = int ( fmod( abs_sec ,3600.0) / 60.0 );

	if ( fmod( abs_sec ,60.0) > 59.9949 ) min++;

	if ( min > 59 ) min = 0;

	return min;
	
}

double AstroAngle::sec( double _angle )
{
	double abs_sec = abs ( _angle * SEC_IN_RAD );
	double sec = floor( fmod( abs_sec ,60.0) * 100.0 + 0.5 ) / 100.0;

	if ( sec > 59.98 ) sec = 0;
	return sec;
}




/************************************************************/

string& AstroAngle::getSign( double _angle, angle_type _type, string& szSign )
{
	szSign="";
	switch (_type)
	{

	  case ANGLE_LONG1_TYPE:
	  case ANGLE_LONG2_TYPE:
	  {
		if ( _angle < 0.0 )
			szSign = "W";		// 'E' for positive
		else 
			szSign = "E";
		break;
	  }

	  case ANGLE_AZIM1_TYPE:
	  case ANGLE_AZIM2_TYPE:
	  {
		if ( _angle < 0.0 )
			szSign = "E";		// 'W' for positive
		else
			szSign = "W";
		break;
	  }

	//case ANGLE_DEC_TYPE:
	  case ANGLE_LAT_TYPE:
	  case ANGLE_ALT_TYPE:
	  {
		if ( _angle < 0.0 )
			szSign = "S";		// 'N' for positive
		else
			szSign = "N";
		break;
	  }

	  default:
	  {
		if ( _angle < 0.0 )
			szSign = "-";		// '+' for positive
		else
			szSign = "+";
	  }

	}

	return szSign;
}

char AstroAngle::getSign( double _angle, angle_type _type )
{
  switch (_type)
	{

	  case ANGLE_LONG1_TYPE:
	  case ANGLE_LONG2_TYPE:
	  {
		if ( _angle < 0.0 ) return 'W';		// 'E' for positive
		else return 'E';
	  }

	  case ANGLE_AZIM1_TYPE:
	  case ANGLE_AZIM2_TYPE:
	  {
		if ( _angle < 0.0 ) return 'E';		// 'W' for positive
		else return 'W';
	  }

	//case ANGLE_DEC_TYPE:
	  case ANGLE_LAT_TYPE:
	  case ANGLE_ALT_TYPE:
	  {
		if ( _angle < 0.0 ) return 'S';		// 'N' for positive
		else return 'N';
	  }

	  case ANGLE_TIME_TYPE:
	  case ANGLE_HA_TYPE:
	  case ANGLE_RA_TYPE:
	  {
		return ' ';				// no sign for timelike values
	  }

	  default:
	  {
		if ( _angle < 0.0 ) return '-';		// '+' for positive
		else return '+';
	  }

	}
} //End of getSign()


char AstroAngle::getSign ()	// sign of angle ('+' or '-', 'E' or 'W', 'N' or 'S', depending on angle type)
{
	return getSign( angle, type );
} //End of getSign()


/************************************************************/


void AstroAngle::cutToRange( double& _angle, angle_type _type )		// eg. for 0-2PI range cuts to 0 if less then 0 and to 2PI when greater
{
  switch (_type)
	{

	// angles that have '-90 - +90' degrees range
	  case ANGLE_TYPE_1:
	  case ANGLE_DEC_TYPE:
	  case ANGLE_LAT_TYPE:
	  case ANGLE_ALT_TYPE:
	  {
		if ( _angle < - PI_VALUE/2.0 ) { _angle = - PI_VALUE/2.0;   break; }
		if ( _angle >   PI_VALUE/2.0 )   _angle =   PI_VALUE/2.0;
		break;
	  }

	// angles that have '0 - +180' degrees range
	  case ANGLE_TYPE_3:
	  case ANGLE_ZD_TYPE:
	  {
		if ( _angle < 0.0 ) { _angle = 0.0;   break; }
		if ( _angle >  PI_VALUE )   _angle =  PI_VALUE;
		break;
	  }

		case ANGLE_AZIM1_TYPE:
			{
				if(_angle<0)
					_angle += TWO_PI_VALUE;
				if(_angle>TWO_PI_VALUE)
					_angle -= TWO_PI_VALUE;
				break;
			}

	  default:
	  {
		cout<<"You have just tried to cut to range angle that is periodic."<<endl
		    <<"Nothing happend but it would be pleasant if you check why this occured."<<endl;
	  }

	}
} //End of cutToRange()


/************************************************************/

void AstroAngle::shiftToRange( double& _angle, angle_type _type )
{
  switch (_type)
	{

	// angles that have '0 - +360' degrees range
	  case ANGLE_0_2PI:
	  case ANGLE_TIME_TYPE:
	  case ANGLE_HA_TYPE:
	  case ANGLE_RA_TYPE:
	  case ANGLE_LONG1_TYPE:
	  case ANGLE_AZIM1_TYPE:
	  {
		shiftToRangeRAType( _angle );
		break;
	  }

	// angles that have '-180 - +180' degrees range
	  case ANGLE_TYPE_2:
	  case ANGLE_LONG2_TYPE:
	  case ANGLE_AZIM2_TYPE:
	  {
		if ( ( _angle <= - PI_VALUE ) || ( _angle > PI_VALUE ) )
		   _angle = _angle - 2*PI_VALUE * ( floor( (_angle + PI_VALUE)/(2*PI_VALUE) + 0.000000000001 ) );
		break;
	  }

	  default:
	  {
		cout<<"You have just tried to shift to range angle that isn't periodic."<<endl
		    <<"Nothing happend but it would be pleasant if you check why this occured."<<endl;
	  }

	}

}


/************************************************************/


void AstroAngle::setDefaultForPole( double& _angle, angle_type _type, double& out_angle, angle_type _out_type )
{


  switch (_out_type)
	{

	// angles that have a 'zero' value at Poles 
	  case ANGLE_TIME_TYPE:
	  case ANGLE_HA_TYPE:
	  case ANGLE_RA_TYPE:
	  case ANGLE_LONG1_TYPE:
	  case ANGLE_LONG2_TYPE:
	  case ANGLE_AZIM1_TYPE:
	  case ANGLE_AZIM2_TYPE:
	  {

  		switch (_type)
			{

			// angles that have Poles at -90 and +90 deg
			case ANGLE_TYPE_1:
			case ANGLE_DEC_TYPE:
			case ANGLE_LAT_TYPE:
			case ANGLE_ALT_TYPE:
			{
				if ( (_angle > PI_2_VALUE - 0.00000242)
				     || (_angle < -PI_2_VALUE + 0.00000242) )
				out_angle = 0.0;
				break;
			}
			case ANGLE_TYPE_3:
			case ANGLE_ZD_TYPE:
			{
				if ( (_angle > PI_VALUE - 0.00000242)
				     || (_angle < 0.00000242) )
				out_angle = 0.0;
				break;
			}
			default:
			{
				cout<<"You have specified to check for Pole angle that has not any Pole."<<endl
				<<"Nothing happend but it would be pleasant if you check why this occured."<<endl;
			}

			}
		break;
	  }

	  default:
	  {
		cout<<"You have just tried to set default value for angle that has not any Pole."<<endl
		    <<"Nothing happend but it would be pleasant if you check why this occured."<<endl;
	  }

	}

}

double AstroAngle::getDist( double angle1 , double angle2 )
{
	double dist_ang = fabs(angle1-angle2);
   if ( dist_ang > PI_VALUE ){
      dist_ang = TWO_PI_VALUE - dist_ang;
   }

	return dist_ang;
}

double AstroAngle::getDistDeg( double angle1 , double angle2 )
{
	double dist_ang = fabs(angle1-angle2);
   if ( dist_ang > 180. ){
      dist_ang = 360. - dist_ang;
   }

	return dist_ang;
}


double AstroAngle::getDist( double ra1, double dec1, double ra2, double dec2 )
{
	double dist_ra = getDist(ra1,ra2);
	double dist_dec = getDist(dec1,dec2);

	double dist = sqrt( dist_ra*dist_ra + dist_dec*dist_dec );
	return dist;
}



int AstroAngle::CompareAngle( double angle1, double angle2 )
{
	if(angle1>angle2)
		return -1;
	if(angle1<angle2)
		return 1;	
	return 0;
}
