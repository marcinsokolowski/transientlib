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
#include "satcone.h"
#include "ccd_globals.h"
#include <AstroBody.h>
#include <AstroCCD.h>
#include <mathfunc.h>

double gE = 0.0000000;
double gC = 0.0000000;

int CSatCone::calc_min_dist( double ra_in_deg, double dec_in_deg, 
									  time_t unix_time, double& r_min )
{
	double Rz = EARTH_RADIUS;
	double h = ( SUN_DIST / ( (SUN_RADIUS/Rz) - 1 ) );

	double fi_l = AstroAngle::deg2rad(ra_in_deg);
	double cos_fi_l = cos(fi_l);
   double sin_fi_l = sin(fi_l);
	double teta_l = AstroAngle::deg2rad( 90 - dec_in_deg );

	// Julian date :
	double jd = AstroCCD::getJulianDay(unix_time);

	// calc anty-solar point :
	double sun_ra_in_deg,sun_dec_in_deg;
	CAstroBody::calc_sun( jd, sun_ra_in_deg, sun_dec_in_deg );

	double lambda_as = ( 180 + sun_ra_in_deg );
	if( lambda_as>360 )
	{
		lambda_as = lambda_as - 360;
	}	
	lambda_as = AstroAngle::deg2rad(lambda_as);
	double delta_as = ( -AstroAngle::deg2rad(sun_dec_in_deg) );

// TEMPORARY !!!!!!!!!!!
//	lambda_as = 0;
//	delta_as = 0;
// TEMPORARY !!!!!!!!!!!
	

	double sin_teta_l = sin( teta_l );
	double cos_teta_l = cos( teta_l );

	double cos_lambda_as = cos( lambda_as );
	double sin_lambda_as = sin( lambda_as );			

	double cos_delta_as = cos( delta_as );
   double sin_delta_as = sin( delta_as );			

	double cos_phi = cos( gCCDParams.m_GeoLatitude );
	double sin_phi = sin( gCCDParams.m_GeoLatitude );


	// siderial time :
	double T_sid = AstroCCD::getSiderealTimeLocal( unix_time, gCCDParams.m_GeoLongitude );
	printf("T_sid = %.2f [rad] = %.2f [deg]\n",T_sid,AstroAngle::rad2deg(T_sid));
	double sun_ha_in_rad;
	AstroCCD::calculateHourAngleBase( AstroAngle::deg2rad( sun_ra_in_deg ),
												 AstroAngle::deg2rad( sun_dec_in_deg ),
												 T_sid, sun_ha_in_rad );
	double sun_ha_in_deg= AstroAngle::rad2deg( sun_ha_in_rad );
//	double T = ( PI_VALUE - sun_ha_in_rad );
	double T = T_sid;
	double cos_T = cos(T);
	double sin_T = sin(T);


	printf("####################################################\n");
	printf("calculating minimume distatnce to be outside of EARTH shadow\n");
	printf("jd = %.8f\n",jd);
	printf("Sun position     (RA,DEC) = (%.2f,%.2f) [deg]\n",sun_ra_in_deg,sun_dec_in_deg);
	printf("Solar HA = %.2f [rad] = %.2f [deg]\n",sun_ha_in_rad,sun_ha_in_deg);
	printf("Observatory (Long,Lat) = (%.2f,%.2f) [deg]\n",AstroAngle::rad2deg(gCCDParams.m_GeoLongitude), AstroAngle::rad2deg(gCCDParams.m_GeoLatitude));
	printf("cos_phi=%.2f, sin_phi=%.2f\n",cos_phi,sin_phi);
	printf("Anty-solar point (RA,DEC) = (%.2f,%.2f) [deg]\n",AstroAngle::rad2deg(lambda_as),AstroAngle::rad2deg(delta_as));	
	printf("T = %.2f [rad] = %.2f [deg]\n",T,AstroAngle::rad2deg(T));
	printf("Cone h = %.2f [km]\n",h);


	double vec_x = Rz*cos_phi*cos_T;
	double vec_y = Rz*cos_phi*sin_T;
	double vec_z = Rz*sin_phi;
	double A = ( sin_teta_l*cos(fi_l-lambda_as)*cos_delta_as + cos_teta_l*sin_delta_as );
	double B = Rz*( cos_phi*cos(  T - lambda_as )*cos_delta_as  + 
				       sin_phi*sin_delta_as );
		
	double C = sin_teta_l*sin(fi_l-lambda_as);
	double D = Rz*cos_phi*sin(  T - lambda_as );

	double E = -sin_teta_l*cos(fi_l-lambda_as)*sin_delta_as + cos_teta_l*cos_delta_as;
	double F = Rz*( -cos_phi*cos(  T - lambda_as )*sin_delta_as  + sin_phi*cos_delta_as );
	double f = (F - h);
	double p = (B - h );


	double gamma = (Rz/h);
	double gamma2 = gamma*gamma;

	printf("A=%.2f B=%.2f\n",A,B);
	printf("C=%.2f D=%.2f\n",C,D);
	printf("E=%.2f F=%.2f\n",E,F);
	printf("f=%.2f gamma=%.4f\n",f,gamma);
	printf("line in cone frame : ( %.2f t + %2.f, %.2f t + %2.f, %.2f t + %2.f )\n",
				A,B,C,D,E,F);

	// double delta = CMyMathFunc::mysqr(A*B+D*C-gamma*gamma*E*f) - (A*A+C*C-E*E*gamma*gamma)*(B*B+D*D-gamma*gamma*f*f);
	
//	double a = (A*A+C*C-E*E*gamma*gamma);
//	double b = 2*(A*B+D*C-gamma*gamma*E*f);
//	double c = (B*B+D*D-gamma*gamma+f*f);

	double a = (E*E+C*C-gamma2*A*A);
	double b = 2.00*(E*F+C*D-gamma2*A*p);
	double c = (F*F + D*D - gamma2*p*p);

	printf("a=%.2f b=%.2f c=%.2f\n",a,b,c);
	double delta = b*b-4*a*c;
	printf("delta = %.2f\n",delta);	

	if( delta >= 0 ){
		double x1 = (-b+sqrt(delta))/(2*a);
		double x2 = (-b-sqrt(delta))/(2*a);

		printf("x1 = %.2f\n",x1);
		printf("x2 = %.2f\n",x2);

		double x_bis = x1*A+B;
		double y_bis = x1*C+D;
		double z_bis = x1*E+F;
		double r1 = sqrt(x_bis*x_bis+y_bis*y_bis+z_bis*z_bis);

		double x_bis2 = x2*A+B;
		double y_bis2 = x2*C+D;
		double z_bis2 = x2*E+F;
		double r2 = sqrt(x_bis2*x_bis2+y_bis2*y_bis2+z_bis2*z_bis2);


		printf("r1=%.2f [km] (%.2f,%.2f,%.2f)\n",r1,x_bis,y_bis,z_bis);
		printf("r2=%.2f [km] (%.2f,%.2f,%.2f)\n",r2,x_bis2,y_bis2,z_bis2);

		/*if( x_bis>0 && x_bis<h ){
			printf("Out soltion is r1\n");
		}
		if( x_bis2>0 && x_bis2<h ){
			printf("Out soltion is r2\n");
		}*/

/*		double x11,y11,z11,x22,y22,z22,x3,y3,z3;
		CMyMathFunc::calc_rot_y_declin(  x_bis2, y_bis2, z_bis2, delta_as, x11, y11, z11 );
		CMyMathFunc::calc_rot_z(  x11, y11, z11, lambda_as, x22, y22, z22 );
		CMyMathFunc::shift_vec( x22, y22, z22, vec_x, vec_y, vec_z, x3 ,y3, z3 );
		printf("(x3,y3,z3) = (%.2f,%.2f,%.2f)\n",x3 ,y3, z3 );
		double r_prim = sqrt( x3*x3+y3*3+z3*z3 );
		double teta_prim = acos( z3/r_prim );
		double fi_prim = atan2( y3 , x3 );
		printf("(teta,fi) = (%.2f,%.2f) [deg]\n",AstroAngle::rad2deg(teta_prim),
						AstroAngle::rad2deg(fi_prim));
*/
		double vec_x = sin_teta_l*cos_fi_l;
		double vec_y = sin_teta_l*sin_fi_l;
		double vec_z = cos_teta_l;

		double x_start1 = x1*sin_teta_l*cos_fi_l;
		double y_start1 = x1*sin_teta_l*sin_fi_l;
		double z_start1 = x1*cos_teta_l;

		double x_start2 = x2*sin_teta_l*cos_fi_l;
		double y_start2 = x2*sin_teta_l*sin_fi_l;
		double z_start2 = x2*cos_teta_l;
		printf("\n\nIn LCO frame :\n");
		printf("Start pos1 = (%.2f,%.2f,%.2f)\n",x_start1,y_start1,z_start1);
		printf("Start pos2 = (%.2f,%.2f,%.2f)\n",x_start2,y_start2,z_start2);
		printf("Direction of vector (%.2f,%.2f,%.2f)\n",vec_x,vec_y,vec_z);
		
		double il1 = vec_x*x_start1+vec_y*y_start1+vec_z*z_start1;
		double il2 = vec_x*x_start2+vec_y*y_start2+vec_z*z_start2;
		printf("il1=%.2f, il2=%.2f\n",il1,il2);
		if( il1 > 0 ){
			printf("Out solution is r1 = %.2f [km]\n",r1);
			r_min = r1;
		}
		if( il2 > 0 ){
			printf("Out solution is r2 = %.2f [km]\n",r2);
			r_min = r2;
		}
	}else{
		printf("delta < 0 - no solution !\n");
	}
}


int CSatCone::get_earth_centered_line( double ra_in_deg, double dec_in_deg, time_t unix_time, double d )
{
	double Rz = EARTH_RADIUS;
	double h = ( SUN_DIST / ( (SUN_RADIUS/Rz) - 1 ) );

	double fi_l = AstroAngle::deg2rad(ra_in_deg);
	double cos_fi_l = cos(fi_l);
	double sin_fi_l = sin(fi_l);
	double teta_l = AstroAngle::deg2rad( 90 - dec_in_deg );

	// Julian date :
	double jd = AstroCCD::getJulianDay(unix_time);

	double sin_teta_l = sin( teta_l );
	double cos_teta_l = cos( teta_l );

	double cos_phi = cos( gCCDParams.m_GeoLatitude );
	double sin_phi = sin( gCCDParams.m_GeoLatitude );

	// calc anty-solar point :
   double sun_ra_in_deg,sun_dec_in_deg,sun_ha_in_rad;
   CAstroBody::calc_sun( jd, sun_ra_in_deg, sun_dec_in_deg );

	// siderial time :
	double T_sid = AstroCCD::getSiderealTimeLocal( unix_time, gCCDParams.m_GeoLongitude );
	printf("T_sid = %.2f [rad] = %.2f [deg]\n",T_sid,AstroAngle::rad2deg(T_sid));
	AstroCCD::calculateHourAngleBase( AstroAngle::deg2rad( sun_ra_in_deg ),
												 AstroAngle::deg2rad( sun_dec_in_deg ),
												 T_sid, sun_ha_in_rad );
	double sun_ha_in_deg= AstroAngle::rad2deg( sun_ha_in_rad );
	printf("Solar HA = %.2f [rad] = %.2f [deg]\n",sun_ha_in_rad,sun_ha_in_deg);
	printf("Observatory (Long,Lat) = (%.2f,%.2f) [deg]\n",AstroAngle::rad2deg(gCCDParams.m_GeoLongitude), AstroAngle::rad2deg(gCCDParams.m_GeoLatitude));
	printf("cos_phi=%.2f, sin_phi=%.2f\n",cos_phi,sin_phi);
	double T = ( PI_VALUE - sun_ha_in_rad );
	double cos_T = cos(T);
	double sin_T = sin(T);
	
	
	double x_lco = sin_teta_l*cos_fi_l*d;
	double y_lco = sin_teta_l*sin_fi_l*d;
	double z_lco = cos_teta_l*d;

	double ra_test = atan2( y_lco , x_lco );
	double ra2 = AstroAngle::rad2deg( ra_test );

	printf("LCO : (x_lco,y_lco,z_lco) = (%.2f,%.2f,%.2f) [km]\n",x_lco,y_lco,z_lco);
	printf("LCO : (RA,DEC)            = (%.4f,%.4f)      [deg] = (%.4f,%.4f)\n",
				ra_in_deg,dec_in_deg,ra2,dec_in_deg);

	double x_earth = x_lco + Rz*cos_phi*cos_T;
	double y_earth = y_lco + Rz*cos_phi*sin_T;
	double z_earth = z_lco + Rz*sin_phi;

	printf("EARTH : (x_e,y_e,z_e) = (%.2f,%.2f,%.2f) [km]\n",x_earth,y_earth,z_earth);

	double d_earth = sqrt( x_earth*x_earth + y_earth*y_earth + z_earth*z_earth );
	double teta_e = acos( z_earth/d );
	double phi_e = atan2( y_earth, x_earth );
	if( phi_e < 0 )
		phi_e = ( 2*PI_VALUE + phi_e );
	double delta_e = ( PI_VALUE/2.00 - teta_e );	

	printf("EARTH : (RA_e,DEC_e) = (%.4f,%.4f)\n",AstroAngle::rad2deg(phi_e),AstroAngle::rad2deg(delta_e));
}

int CSatCone::get_local_line( double ra_in_deg, double dec_in_deg, time_t unix_time, double d )
{
	double Rz = EARTH_RADIUS;
	double h = ( SUN_DIST / ( (SUN_RADIUS/Rz) - 1 ) );

	double fi_l = AstroAngle::deg2rad(ra_in_deg);
	double cos_fi_l = cos(fi_l);
	double sin_fi_l = sin(fi_l);
	double teta_l = AstroAngle::deg2rad( 90 - dec_in_deg );

	// Julian date :
	double jd = AstroCCD::getJulianDay(unix_time);

	double sin_teta_l = sin( teta_l );
	double cos_teta_l = cos( teta_l );

	double cos_phi = cos( gCCDParams.m_GeoLatitude );
	double sin_phi = sin( gCCDParams.m_GeoLatitude );

	// calc anty-solar point :
   double sun_ra_in_deg,sun_dec_in_deg,sun_ha_in_rad;
   CAstroBody::calc_sun( jd, sun_ra_in_deg, sun_dec_in_deg );

	// siderial time :
	double T_sid = AstroCCD::getSiderealTimeLocal( unix_time, gCCDParams.m_GeoLongitude );
	printf("T_sid = %.2f [rad] = %.2f [deg]\n",T_sid,AstroAngle::rad2deg(T_sid));
	AstroCCD::calculateHourAngleBase( AstroAngle::deg2rad( sun_ra_in_deg ),
												 AstroAngle::deg2rad( sun_dec_in_deg ),
												 T_sid, sun_ha_in_rad );
	double sun_ha_in_deg= AstroAngle::rad2deg( sun_ha_in_rad );
	printf("Solar HA = %.2f [rad] = %.2f [deg]\n",sun_ha_in_rad,sun_ha_in_deg);
	printf("Observatory (Long,Lat) = (%.2f,%.2f) [deg]\n",AstroAngle::rad2deg(gCCDParams.m_GeoLongitude), AstroAngle::rad2deg(gCCDParams.m_GeoLatitude));
	printf("cos_phi=%.2f, sin_phi=%.2f\n",cos_phi,sin_phi);
	double T = ( PI_VALUE - sun_ha_in_rad );
	double cos_T = cos(T);
	double sin_T = sin(T);
	
	
	double x_earth = sin_teta_l*cos_fi_l*d;
   double y_earth = sin_teta_l*sin_fi_l*d;
   double z_earth = cos_teta_l*d;

	double ra_test = atan2( y_earth , x_earth );
	double ra2 = AstroAngle::rad2deg( ra_test );

	printf("EARTH: (x_earth,y_earth,z_earth) = (%.2f,%.2f,%.2f) [km]\n",x_earth,y_earth,z_earth);
	printf("EARTH: (RA,DEC)            = (%.4f,%.4f)      [deg] = (%.4f,%.4f)\n",
				ra_in_deg,dec_in_deg,ra2,dec_in_deg);

	double x_local = x_earth - Rz*cos_phi*cos_T;
	double y_local = y_earth - Rz*cos_phi*sin_T;
	double z_local = z_earth - Rz*sin_phi;

	printf("LOCAL: (x_local,y_local,z_local) = (%.2f,%.2f,%.2f) [km]\n",x_local,y_local,z_local);

	double d_local  = sqrt( x_local*x_local + y_local*y_local + z_local*z_local );
	double teta_local = acos( z_local/d_local );
	double phi_local = atan2( y_local, x_local );
	if( phi_local < 0 )
		phi_local = ( 2*PI_VALUE + phi_local );
	double delta_local = ( PI_VALUE/2.00 - teta_local );	

	int phi_deg, phi_min, delta_deg, delta_min;
	double phi_sec,delta_sec;
	AstroAngle::deg2deg( AstroAngle::rad2deg(phi_local), phi_deg, phi_min, phi_sec );
	AstroAngle::deg2deg( AstroAngle::rad2deg(delta_local), delta_deg, delta_min, delta_sec );

	printf("LOCAL: (RA,DEC) = (%.4f,%.4f) = ( %.2d:%.2d:%.2f, %.2d:%.2d:%.2f )\n",
				AstroAngle::rad2deg(phi_local),AstroAngle::rad2deg(delta_local),
				phi_deg, phi_min, phi_sec, delta_deg, delta_min, delta_sec);
}

int CSatCone::calc_cone_line_crossing()
{
	double A=1.000000000;
	double E=gE;
	double C=gC;
	printf("E=%.2f, C=%.2f\n",E,C);

	double Rz = EARTH_RADIUS;
   double h = ( SUN_DIST / ( (SUN_RADIUS/Rz) - 1 ) );
	double gamma  = (Rz/h);
	double gamma2 = (gamma*gamma);

	double a = (E*E + C*C - gamma2);
	double b = 2.00*gamma2*h;
	double c = -gamma2*h*h;		

	double delta,x1,x2;
	int ret = CMyMathFunc::calc_sqr_eq( a, b, c, delta, x1, x2 );
	if ( ret>0 ){
		printf("Number of crossings = %d\n",ret);
		if( ret==1 ){
			printf("d = %.4f\n",x1);
			double x = A*x1;
			double y = C*x1;
			double z = E*x1;
			printf("(x,y,z) = (%.2f,%.2f,%.2f)\n",x,y,z);
		}else{
			printf("d1 = %.4f , d2 = %.4f\n",x1,x2);
			double x = A*x1;
			double y = C*x1;
			double z = E*x1;
			printf("(x1,y1,z1) = (%.2f,%.2f,%.2f)\n",x,y,z);
			x = A*x2;
			y = C*x2;
			z = E*x2;
			printf("(x2,y2,z2) = (%.2f,%.2f,%.2f)\n",x,y,z);
		}
	}else{
		printf("line does not cross a cone delta=%.20f\n",delta);
	}

	return ret;
}

