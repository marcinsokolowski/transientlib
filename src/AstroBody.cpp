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
#include "AstroBody.h"
#include <mathfunc.h>
#include "AstroDefines.h"
#include "Astroangle.h"
#include <mydate.h>
#include <myutil.h>

int CAstroBody::m_CalcPosStepInSec=5;

CBodyPosition::CBodyPosition( time_t t, double _ra, double _dec, double _az, 
										double _h, double _illum )
: unix_time(t),ra(_ra),dec(_dec),h(_h),az(_az),illum(_illum)
{
}


CBodyPosition* CBodyPositions::FindRise( int& pos )
{
	BOOL_T bNegat=FALSE;
	for(int i=0;i<size();i++){
		if( (*this)[i].h < 0 ){
			bNegat = TRUE;
		}
		if( (*this)[i].h > 0 && bNegat ){
			pos = i;
			return &((*this)[i]);
		}
	}
	return NULL;
}

CBodyPosition* CBodyPositions::FindSet( int& pos )
{
	BOOL_T bPosit=FALSE;
	for(int i=0;i<size();i++){
		if( (*this)[i].h > 0 ){
			bPosit = TRUE;
		}
		if( (*this)[i].h < 0 && bPosit ){
			pos = i;
			return &((*this)[i]);
		}
	}
	return NULL;
}



void CBodyPositions::Add( time_t unix_time, double ra, double dec, double az, 
								  double h, double illum )
{
	CBodyPosition tmp( unix_time, ra, dec, az, h, illum );
	push_back( tmp );
}

CAstroBody::CAstroBody()
: illum(0),body_ra(0),body_dec(0),body_ra_in_deg(0),body_dec_in_deg(0),
  m_pSatInfo( NULL )
{
}

CAstroBody::CAstroBody( const char* szName, double ra, double dec, double fov )
: body_ra( ra ), body_dec( dec ), body_ra_in_deg(0), body_dec_in_deg(0),
  m_FOV(fov), m_pSatInfo(NULL)
{
	if( szName && szName[0] )
		m_szName = szName;

	body_ra_in_deg = AstroAngle::rad2deg( ra );
	body_dec_in_deg = AstroAngle::rad2deg( dec );
}


double CAstroBody::calc_height( double JD_date , double lat, double lon,
										 double body_ra, double body_dec )
{
	double sidT_in_rad = getSiderealTimeAtGreenwichFromJD( JD_date );
	double sidT = AstroAngle::rad2deg( sidT_in_rad );
	printf("sid_time = %.2f deg\n",sidT);

	double t_h = sidT - body_ra + lon;
	double h = asind( sind(lat)*sind(body_dec) + cosd(lat)*cosd(body_dec)*cosd(t_h) );

	return h;
}

int CAstroBody::calc_sun2( double j_date, double &ra, double &dec )
{
	double k = PI_VALUE/180.00;

	// double T = getSiderealTimeAtGreenwichFromJD( j_date );
	// T = AstroAngle::rad2deg( T );	
	double T = ( j_date - 2451545.0 ) / 36525;
	double M = 357.52910 + 35999.05030*T - 0.0001559*T*T - 0.00000048*T*T*T; // mean anomaly, degree

printf("M=%.5f\n",M);

	double L0 = 280.46645 + 36000.76983*T + 0.0003032*T*T; // mean longitude, degree
	double DL = (1.914600 - 0.004817*T - 0.000014*T*T)*sin(k*M)
				+ (0.019993 - 0.000101*T)*sin(k*2*M) + 0.000290*sin(k*3*M);

	double L = L0 + DL; // true longitude, degree 

	double eps = 23.43999; // obliquity of ecliptic

	double X = cosd(L); 
	double Y = cosd(eps)*sind(L);
	double Z = sind(eps)*sind(L); 
	double R = sqrt(1.0-Z*Z);


   double delta = (180/PI_VALUE)*atan(Z/R); // in degrees

   double RA = (24/PI_VALUE)*atan(Y/(X+R)); // in hours
	
//	if( RA<0 ){
//		RA = (2*PI_VALUE+RA);
//	}

	ra = (RA*15.00);
	dec = delta;

	if(ra<0.00)
		ra = (360 + ra);

	return 1;
}

int CAstroBody::calc_sun( double j_date, double &ra, double &dec )
{

double w_sun, a_sun, e_sun, M_sun, oblecl_sun, L_sun, E_sun,
       x_sun, y_sun, z_sun, r_sun, rs, v_sun, lon_sun,
       xequat_sun, yequat_sun, zequat_sun, RA_sun, Decl_sun,
       xeclip_sun, yeclip_sun, zeclip_sun;

double d = j_date - MILHJD;

double Epoch2, lon_corr;                        // correction for precession

Epoch2 = 2000;
lon_corr = 3.82394E-5 * ( 365.2422 * ( Epoch2 - 2000.0 ) - d );

// Calculate Sun Position

w_sun = 282.9404 + 4.70935E-5 * d;              // longitude of perihelion (deg)
a_sun = 1.000000;                               // mean distance, a.u.
e_sun = 0.016709 - 1.151E-9 * d;                // eccentricity
M_sun = 356.0470 + 0.9856002585 * d;            // mean anomaly (deg)
oblecl_sun = 23.4393 - 3.563E-7 * d;            // obliquity of the ecliptic (deg)
L_sun = w_sun + M_sun;                          // mean longitude (deg)
E_sun = M_sun + (180.0/PI_VALUE) * e_sun              // eccentric anomaly
      * sind(M_sun)*(1 + e_sun*cosd(M_sun));
x_sun = cosd(E_sun) - e_sun;                    // rectangular coordinates
y_sun = sind(E_sun) * sqrt(1 - e_sun*e_sun);
r_sun = sqrt(x_sun*x_sun + y_sun*y_sun);        // distance and true anomaly
rs = r_sun;
v_sun = atan2d( y_sun, x_sun );
lon_sun = v_sun + w_sun + lon_corr;             // longitude
xeclip_sun = r_sun * cosd(lon_sun);             // xs
yeclip_sun = r_sun * sind(lon_sun);             // ys
zeclip_sun = 0.0;
xequat_sun = xeclip_sun;
yequat_sun = yeclip_sun * cosd(oblecl_sun) + zeclip_sun * sind(oblecl_sun);
zequat_sun = yeclip_sun * sind(oblecl_sun) + zeclip_sun * cosd(oblecl_sun);
r_sun = sqrt( x_sun*x_sun + y_sun*y_sun + z_sun*z_sun );

RA_sun = atan2d( yequat_sun, xequat_sun )/ 15.0;	// hours
Decl_sun = atan2d( zequat_sun, sqrt(xequat_sun*xequat_sun+yequat_sun*yequat_sun) );

w_sun = rev(w_sun);
M_sun = rev(M_sun);
L_sun = rev(L_sun);
oblecl_sun = rev(oblecl_sun);
E_sun = rev(E_sun);
v_sun = rev(v_sun);

/*
xe_sun = xeclip_sun;
ye_sun = yeclip_sun;
ze_sun = zeclip_sun;
*/

ra = RA_sun*15.0; // deg
dec = Decl_sun;   // deg 

if( ra<0.00 )
	ra = ( 360 + ra );

return 1;
}

double CAstroBody::when( double j_date, double height, double ra, double dec, double la, double lo, double& ah )
{
double hh, T;                            // ALL in deg
double h_beg, jd_beg;
int bOK=0;

for( double xx=j_date; xx<j_date+1; xx+=0.01 )
{
   ah = calc_sidt( xx ) - ra + lo;
   hh = asind( sind(la)*sind(dec) + cosd(la)*cosd(dec)*cosd(ah) );
   if( fabs(hh-height)<=5.0 )
   {
      for( xx; xx<j_date+1; xx+=0.0001 )
      {
         ah = calc_sidt( xx ) - ra + lo;
         hh = asind( sind(la)*sind(dec) + cosd(la)*cosd(dec)*cosd(ah) );
         if( fabs(hh-height)<=0.1 )
         {
            h_beg = hh;
				jd_beg = xx;
				bOK=1;
				break;
         }
      }
   }
   if( fabs(hh-height)<=0.1 )   break;

}

if( !bOK )
	return -1;

return jd_beg;                               // julian date
}

double CAstroBody::when_new( double j_date, double height, double ra, double dec, 
								 double la, double lo, double& ah, double prec )
{
double hh, T;                            // ALL in deg
double h_beg, jd_beg;
int bOK=0;

for( double xx=j_date; xx<j_date+1; xx+=0.01 )
{
   ah = calc_sidt( xx ) - ra + lo;
   hh = asind( sind(la)*sind(dec) + cosd(la)*cosd(dec)*cosd(ah) );
   if( fabs(hh-height)<=5.0 )
   {
      for( xx; xx<j_date+1; xx+=0.0001 )
      {
         ah = calc_sidt( xx ) - ra + lo;
         hh = asind( sind(la)*sind(dec) + cosd(la)*cosd(dec)*cosd(ah) );
         if( fabs(hh-height)<=prec )
         {
            h_beg = hh;
				jd_beg = xx;
				bOK=1;
				break;
         }
      }
   }
   if( fabs(hh-height)<=0.1 )   break;

}

if( !bOK )
	return -1;

return jd_beg;                               // julian date
}


// ------------------------------------------------------------
double CAstroBody::when_sun( double j_date, double height, double ra, double dec, double la, double lo )
{
double ah, hh, T;                            // ALL in deg
double h_beg, jd_beg;

calc_sun(j_date,ra,dec);

for( double xx=j_date; xx<j_date+1; xx+=0.01 )
{
   ah = calc_sidt( xx ) - ra + lo;
   hh = asind( sind(la)*sind(dec) + cosd(la)*cosd(dec)*cosd(ah) );
   if( fabs(hh-height)<=5.0 )
   {
      for( xx; xx<j_date+1; xx+=0.0001 )
      {
         ah = calc_sidt( xx ) - ra + lo;
         hh = asind( sind(la)*sind(dec) + cosd(la)*cosd(dec)*cosd(ah) );
         if( fabs(hh-height)<=0.1 )
         {
            h_beg = hh;   jd_beg = xx;   break;
         }
         calc_sun(xx,ra,dec);
      }
   }
   calc_sun(xx,ra,dec);
   if( fabs(hh-height)<=0.1 )   break;

}

return jd_beg;                               // julian date
}

double CAstroBody::calc_sidt( double JD_date )
{
	double sidT_in_rad = getSiderealTimeAtGreenwichFromJD( JD_date );
   double sidT = AstroAngle::rad2deg( sidT_in_rad );
	
	return sidT;		
}


void CAstroBody::ShowSkyPath()
{
	printf("######### %s POSITION FOR CURRENT NIGHT #############\n",m_szName.c_str());
   for(int i=0;i<m_PosList.size();i++){
      CBodyPosition& pos = m_PosList[i];
      printf("-------------------------------------------------------\n");
      mystring szDTM = get_date_time_string( pos.unix_time );
      double alt,azim;
      double moon_ra, moon_dec, moon_illum;
                                                                                
      printf("%s %s : ( %s,%.2f) position (%.2f,%.2f)\n", m_szName.c_str(),
						szDTM.c_str(),
                  AstroAngle::toString( pos.ra, ANGLE_RA_TYPE ).c_str(),
                  AstroAngle::rad2deg( pos.dec),
                  AstroAngle::rad2deg( pos.az ),
                  AstroAngle::rad2deg( pos.h ) );
   }
   printf("######################################################\n");
}

CBodyPosition* CAstroBody::SetMidnightRADEC( time_t ut_time )
{
	time_t mid_ut = get_midnight_ut( ut_time ) + (24*3600);

	CBodyPositions::iterator i;
	BOOL_T bNegative=FALSE;
	int min_time_dist=(30*60);
	CBodyPosition* pRet=NULL;
	for(i=m_PosList.begin();i!=m_PosList.end();i++){
		if( abs( i->unix_time - mid_ut )<min_time_dist ){
			pRet = &(*i);
			min_time_dist = abs( i->unix_time - mid_ut );
			body_ra  = i->ra;
			body_dec = i->dec; 
			illum    = i->illum;
		}
	}	
	return pRet;
}

CBodyPosition* CAstroBody::GetPosition( time_t ut_time )
{
	CBodyPositions::iterator i;
	BOOL_T bNegative=FALSE;
	CBodyPosition* pRet=NULL;

	if( m_PosList.size()>0 ){	
		if( ut_time <= m_PosList[0].unix_time )
			return &( m_PosList[0] );

		for(int i=1;i<m_PosList.size();i++){
			if( ut_time>m_PosList[i-1].unix_time && ut_time<=m_PosList[i].unix_time ){
				return &( m_PosList[i] );
			}
		}
	}

	return NULL;
}



CBodyPosition* CAstroBody::FindRiseAbove( double h )
{
	CBodyPositions::iterator i;
	BOOL_T bBelow=FALSE;
	for(i=m_PosList.begin();i!=m_PosList.end();i++){
		if( i->h < h ){
			bBelow = TRUE;
		}
		if( i->h>=h && bBelow ){
			return &(*i);
		}
	}		
	return NULL;
}

CBodyPosition* CAstroBody::FindSetBelow( double h )
{
	CBodyPositions::iterator i;
	BOOL_T bPositive=FALSE;
	for(i=m_PosList.begin();i!=m_PosList.end();i++){
		if( i->h > 0 ){
			bPositive = TRUE;
		}
		if( i->h<0 && bPositive )
			return &(*i);
	}		
	return NULL;
}

double CAstroBody::calcDistDeg( double ra1, double dec1, double ra2, double dec2 )
{
	//  double dist = sqrt( (ra1-ra2)*(ra1-ra2) + (dec1-dec2)*(dec1-dec2) );
	double dist_ra = AstroAngle::getDistDeg( ra1, ra2 );
	double dist_dec = AstroAngle::getDistDeg( dec1, dec2 );		
	double dist = sqrt( dist_ra*dist_ra + dist_dec*dist_dec );
	return dist;
}


BOOL_T CAstroBody::CalcPositions( time_t dtm_start, time_t dtm_end,
                   				    double long_in_deg, double lat_in_deg )
{

	double lat_rad = AstroAngle::deg2rad( lat_in_deg );
   double longit_rad = AstroAngle::deg2rad( long_in_deg );
	
	m_PosList.clear();

	mystring szTimeStart = get_gmtime_string( dtm_start );
	mystring szTimeEnd   = get_gmtime_string( dtm_end );
   printf_now4("calculating %s positions for UT time range ( %s-%s )...",m_szName.c_str(),
					 szTimeStart.c_str(), szTimeEnd.c_str()	);
	for(int t=dtm_start;t<dtm_end;t+=m_CalcPosStepInSec){
		// mystring szDTM = get_gmtime_string( t );
		mystring szDTM = get_date_time_string( t );
		double alt,azim;
		AstroCCD::calculateHorizontalCoordinatesFromEq( body_dec,
																		body_ra,
																		t,
																		longit_rad, lat_rad,
																		alt,azim );

		m_PosList.Add( t, body_ra, body_dec, azim, alt );
		//msgprintf("%d (%.2f,%.2f) (%.2f,%.2f)\n",t, AstroAngle::deg2rad( hete_ra ),
      //                         AstroAngle::deg2rad( hete_dec ), azim, alt);
	}
	printf("OK\n");

	return TRUE;	
}


BOOL_T CAstroBody::GetSimpleRADEC( time_t unix_time, double& ra_in_rad, double& dec_in_rad,
                           double& ra_in_deg, double& dec_in_deg, time_t& start_time,
									time_t& end_time )
{
	if( m_pSatInfo ){
		time_t at_time,track_time,total_track_time=0,max_track_time=0;
		time_t t = unix_time;
		BOOL_T bOK=FALSE;
		double prev_ra=-1000,prev_dec=-1000;

		time_t obs_start,obs_end;
		CSatInfo* pInfo = m_pSatInfo->GetInfo( t, ra_in_deg, dec_in_deg, at_time, track_time, obs_start, obs_end );
		if( pInfo ){
			ra_in_rad = AstroAngle::deg2rad( ra_in_deg );
	      dec_in_rad = AstroAngle::deg2rad( dec_in_deg );
			start_time = obs_start;
			end_time = obs_end;
			return TRUE;
		}
	}				

	ra_in_rad = body_ra;
	dec_in_rad = body_dec;
	ra_in_deg = body_ra_in_deg;
	dec_in_deg = body_dec_in_deg;

	return FALSE;

}

BOOL_T CAstroBody::GetRADEC( time_t unix_time, double& ra_in_rad, double& dec_in_rad,
                 				double& ra_in_deg, double& dec_in_deg, 
									double lat, double longit, double min_acceptable_alt,
								   time_t min_track_time, double max_dist, 
									time_t max_future_time )
{
	double lat_rad = AstroAngle::deg2rad( lat );
   double longit_rad = AstroAngle::deg2rad( longit );


	if( m_pSatInfo ){
		time_t at_time,track_time,total_track_time=0,max_track_time=0;
		time_t t = unix_time;
		BOOL_T bOK=FALSE;
		double prev_ra=-1000,prev_dec=-1000;
		time_t best_time = 0;
		double max_alt=-1000.00,best_azim=-1000.00;

		while( (t-unix_time) <= max_future_time ){
			double test_ra,test_dec;
			time_t track_time,start_time,end_time;
			if( m_pSatInfo->GetInfo( t, test_ra, test_dec, at_time, track_time, start_time, end_time ) ){
				// old version - uncomment :
				/* ra_in_deg = test_ra;
				dec_in_deg = test_dec;
				ra_in_rad = AstroAngle::deg2rad( ra_in_deg );
	         dec_in_rad = AstroAngle::deg2rad( dec_in_deg );
				return TRUE;*/

				double azim_in_rad, alt_in_rad;

			   AstroCCD::calculateHorizontalCoordinatesFromEq( AstroAngle::deg2rad( test_dec ),
                                                     AstroAngle::deg2rad( test_ra ),
                                                     t,
                                                     longit_rad, lat_rad,
                                                     alt_in_rad, azim_in_rad );
			   double azim_in_deg = AstroAngle::rad2deg( azim_in_rad );
			   double alt_in_deg  = AstroAngle::rad2deg( alt_in_rad );


				// finding best time and altitude :
				if( alt_in_deg > max_alt ){
					max_alt = alt_in_deg;
					best_azim = azim_in_deg;
					best_time = t;
				}


				if( alt_in_deg < min_acceptable_alt ){
					// skip to low SWIFT positions :
					printf("position (ra,dec)=(%.2f,%.2f) at %d to low h=%.2f\n",test_ra,test_dec,t,alt_in_deg);
					t += 60;
					continue;
				}


				if( strcmp( m_szName.c_str(), "SWIFT" ) ){
					ra_in_deg = test_ra;
               dec_in_deg = test_dec;
               ra_in_rad = AstroAngle::deg2rad( ra_in_deg );
               dec_in_rad = AstroAngle::deg2rad( dec_in_deg );
					return TRUE;
				}

				if( AstroCCD::CalcDistInDeg( prev_ra, prev_dec, test_ra, test_dec )>max_dist ){
					total_track_time = 0;
				}
				total_track_time += track_time;

				if( total_track_time > max_track_time ){
					ra_in_deg = test_ra;
					dec_in_deg = test_dec;
					ra_in_rad = AstroAngle::deg2rad( ra_in_deg );
					dec_in_rad = AstroAngle::deg2rad( dec_in_deg );
					bOK=TRUE;
					max_track_time = total_track_time;
				}

				if( track_time >= min_track_time ){
					return TRUE;
				}
				t += (track_time+5);

				prev_ra = test_ra;
		      prev_dec = test_dec;
			}else{
				t += 60;
			}
		}

		if( bOK ){
			return TRUE;
		}
		printf("No good position found, using default position\n");

		
		if ( max_alt >= (min_acceptable_alt-m_FOV/2.00) ){
			printf("max_alt = %.2f deg at unix_time=%d\n",max_alt,best_time);	
			printf("Will look at FOV , but not in center !\n");

			double  _ha, _dec, _ra;
			AstroCCD::calculateEquatorialCoordinates( AstroAngle::deg2rad( best_azim ),
																   AstroAngle::deg2rad( 30.00 ),
																	best_time, longit_rad, lat_rad,
																	_ha, _dec, _ra );
			ra_in_rad = _ra;
			dec_in_rad = _dec;
			ra_in_deg = AstroAngle::rad2deg( _ra );
			dec_in_deg = AstroAngle::rad2deg( _dec );						
			return TRUE;
		}
	}				
	
	ra_in_rad = body_ra;
	dec_in_rad = body_dec;
	ra_in_deg = body_ra_in_deg;
	dec_in_deg = body_dec_in_deg;

	return FALSE;
}
