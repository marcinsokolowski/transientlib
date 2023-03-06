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
#include "ccd_procstate.h"
#include "ccd_globals.h"
#include "ccd_pipeline.h"
#include <AstroCCD.h>
#include <cfg.h>
#include <mymatrix.h>
#include "ccd_asastransform.h"
#include <mymacros.h>

CCDProcState::CCDProcState( CCDPipeline* pPipeline )
: m_eProcState( eNormalState ), m_pAstroForms(NULL), m_pAsasTransform(NULL),
  m_pPipeline( pPipeline )
{

}

CCDProcState::CCDProcState( CCDAsasTransform* pAsasTransform )
{
	InitAstroForms( pAsasTransform );
}

CCDProcState::~CCDProcState()
{
	if(m_pAstroForms)
		delete m_pAstroForms;
}

void CCDProcState::SetProcState( ePROCSTATE_T eProcState )
{
	m_eProcState = eProcState;
}

void CCDProcState::ResetProcState()
{
	m_eProcState = eNormalState;
}

double AlfaMultPosTab[MAX_PIPELINE_SIZE];
   double AlfaMultNegTab[MAX_PIPELINE_SIZE];

void CCDProcState::ReCalcAlfaMultTab( double alfa, int max_times )
{
	for(register int i=0;i<max_times;i++){
		double mult  = alfa*max_times;
		AlfaMultPosTab[i] = mult;
		AlfaMultNegTab[i] = -mult;
	}
}


CCDProcState::CCDProcState( CCDAsasTransform* pAsasTransform, eObservationMode_T obsMode )
{
	InitAstroForms( pAsasTransform, obsMode );
}

void CCDProcState::InitAstroForms( CCDAsasTransform* pAsasTransform )
{
	m_pAstroForms = new AstroCCD( gCCDParams.m_TransformCCDFocus, 
											gCCDParams.m_TransformCCDPixelSize,
											gCCDParams.m_SizeX,gCCDParams.m_SizeY,
											gCCDParams.m_GeoLatitude,
											gCCDParams.m_GeoLongitude,																						
											gCCDParams.m_HorAzimuth,
											gCCDParams.m_HorAltitude,
											gCCDParams.m_TransformCCDOrientation,
											gCCDParams.m_CCDAxixXOrientation,
											gCCDParams.m_CCDAxixYOrientation,
											gCCDParams.m_ObsMode, gCCDParams.m_DecObs,
											gCCDParams.m_RAObs );

	m_pAsasTransform = pAsasTransform;	
}

void CCDProcState::InitAstroForms( CCDAsasTransform* pAsasTransform, eObservationMode_T obsMode )
{
	m_pAstroForms = new AstroCCD( gCCDParams.m_TransformCCDFocus, 
											gCCDParams.m_TransformCCDPixelSize,
											gCCDParams.m_SizeX,gCCDParams.m_SizeY,
											gCCDParams.m_GeoLatitude,
											gCCDParams.m_GeoLongitude,																						
											gCCDParams.m_HorAzimuth,
											gCCDParams.m_HorAltitude,
											gCCDParams.m_TransformCCDOrientation,
											gCCDParams.m_CCDAxixXOrientation,
											gCCDParams.m_CCDAxixYOrientation,
											obsMode, gCCDParams.m_DecObs,
											gCCDParams.m_RAObs );

	m_pAsasTransform = pAsasTransform;	
}


void CCDProcState::InitAstroForms( double focus, double pixel_size_to_set, int x_size, int y_size,
											  double geoLat, double geoLong,
											  double geoAzim, double geoAlt, double orientation, 
											  int time_zone_to_set, int _X_orient, int _Y_orient,
											  eObservationMode_T obsMode, double dec_obs, double ra_obs,
											  double ut_to_set, int day_to_set, int month_to_set, 
											  int year_to_set, CCDAsasTransform* pAsasTransform )
{
	m_pAstroForms = new AstroCCD( focus, pixel_size_to_set,
											x_size, y_size, 
											geoLat, geoLong,
											geoAzim, geoAlt ,
											orientation, _X_orient, _Y_orient,
											obsMode, dec_obs, ra_obs );

//											gCCDParams.m_TimeZone, ut_angle,
//											day_to_set, month_to_set, year_to_set );
	m_pAsasTransform = pAsasTransform;
}

void CCDProcState::GetObsCoo( time_t ut_time, double& ra, double& dec, 
										double& alt, double& azim, double& ha,
									   eObservationMode_T& obsMode )
{
	obsMode = m_pAstroForms->GetObsMode();
	printf("CCDProcState::GetObsCoo obsMode = %d\n",obsMode);
	
	if( obsMode==eNoMovingMode ){
		m_pAstroForms->GetObsAzim( azim, alt );
		m_pAstroForms->calculateEquatorialCoordinates( azim, alt, ut_time, ha, dec, ra );
	}else{
		m_pAstroForms->calculateHourAngleObs();
		m_pAstroForms->GetObsEq( ra, dec, ha );
		m_pAstroForms->calculateHorizontalCoordinatesFromEq( dec, ra, ut_time, alt, azim );
	}
}

void CCDProcState::SetFrameDateTime( time_t ut_time, CCDPipeline* pPipeline )
{
	// AstroAngle ut(ANGLE_TIME_TYPE);
   // ut.setInHours( _tm->tm_hour, _tm->tm_min,_tm->tm_sec );
   // m_pAstroForms->setDate( _tm->tm_mday, _tm->tm_mon, _tm->tm_year );
	m_pAstroForms->SetUT( ut_time );

	if(m_pAstroForms->GetObsMode()==eNoMovingMode){
		// constant - no earth compensation movement :
		m_pAstroForms->calculateEquatorialCoordinatesObs();
		if( pPipeline ){
			m_pAstroForms->GetObsEq( (pPipeline->m_PipelineCfg).m_RAObs, 
											 (pPipeline->m_PipelineCfg).m_DecObs, 
											 (pPipeline->m_PipelineCfg).m_HAObs );
		}
	}else{			
		// earth rotation comensation :
		m_pAstroForms->calculateHorizontalCoordinatesObs();
		if( pPipeline ){
			m_pAstroForms->GetObsAzim( (pPipeline->m_PipelineCfg).m_HorAzimuth, (pPipeline->m_PipelineCfg).m_HorAltitude );
		}
	}
	// printf("RA_obs=%s , Dec_Obs=%.2f\n",AstroAngle::toString( m_pAstroForms->RA_obs, ANGLE_RA_TYPE ).c_str(), AstroAngle::rad2deg(m_pAstroForms->Dec_obs) );
}


/*void CCDProcState::xyAfterTime( double x, double y, double time_in_sec, double& x_new, double& y_new )
{
//	m_pAstroForms->XYAfterTime( x, y, time_in_sec, x_new, y_new );	
//	m_pAstroForms->setPointXY(x, y);

	m_pAstroForms->XYAfterTimeTEST(x, y, x_new, y_new, time_in_sec);
}

void CCDProcState::xyAfterTimeNEW( double x, double y, double time_in_sec, double& x_new, double& y_new )
{
	// m_pAstroForms->XYAfterTimeNEW( x, y, x_new,  y_new, time_in_sec );
	// m_pAstroForms->XYAfterTime(x, y, x_new, y_new, time_in_sec);
	m_pAstroForms->XYAfterTimeTEST( x, y, x_new,  y_new, time_in_sec );
}*/



BOOL_T CCDProcState::CalcAzimutalCoord( int x, int y, time_t ut_time, 
							  					    double& azim, double& altit )
{
	double ra,dec,ha;
	m_pAstroForms->SetUT( ut_time );
	// m_pAstroForms->calculatePointEquatorial( x, y, dec, ra );
	// m_pAstroForms->calculateHorizontalCoordinatesFromEq( dec, ra, altit, azim );

// test 2
	//m_pAstroForms->calculatePointEquatorial2( x, y, ra, dec, ha, altit, azim );	
//	printf("ra=%s , dec=%.2f\n",AstroAngle::toString( ra, ANGLE_RA_TYPE ).c_str(), AstroAngle::rad2deg(dec) );

	// test 3
	// m_pAstroForms->calculatePointEquatorialNEW( x, y ,ra ,dec );
	// m_pAstroForms->calculateHourAngle( ra, dec, ha );
	// m_pAstroForms->calculateHorizontalCoordinatesFromEq( dec, ha, altit, azim );

	// NEW pretty working version :
	m_pAstroForms->calculatePointEquatorialTest( x,y, ra, dec, ha, altit, azim );	

	return TRUE;				
}

void CCDProcState::azh2ad( double az_in_rad, double h_in_rad, time_t unix_time,
				               double& ra_in_rad, double& dec_in_rad )
{
	double ha_in_rad;
	m_pAstroForms->calculateEquatorialCoordinates( az_in_rad, h_in_rad,
									unix_time, ha_in_rad, dec_in_rad, ra_in_rad );
}

void CCDProcState::ad2azh( double ra_in_rad, double dec_in_rad,
									time_t unix_time,
				               double& az_in_rad, double& h_in_rad )
{
	m_pAstroForms->calculateHorizontalCoordinatesFromEq( dec_in_rad, ra_in_rad, 
																		  unix_time,
																			h_in_rad, az_in_rad );
}

BOOL_T CCDProcState::ad2xy_float( double ra, double dec, double& x, double& y )
{
	if( gCCDParams.m_bUseAsasTransform && m_pAsasTransform && m_pAsasTransform->m_bTransformOK ){
		// first calculate horizontal
		/*double alt,azim;
		m_pAstroForms->calculateHorizontalCoordinatesFromEq( dec , ra, 
																			  m_pAsasTransform->m_TransformUtTime,
																			  alt, azim );


		// then calculate equatorial for new time 
		double ha_new,ra_new,dec_new;
		m_pAstroForms->calculateEquatorialCoordinates( azim, alt, ha_new, dec_new, ra_new );

		// and from ra_new,dec_new calculate x,y :
		m_pAsasTransform->ad2xy( ra_new, dec_new, x, y );*/

		if( m_pAstroForms->GetObsMode()==eNoMovingMode ){
			// IDEA :
			// - calculate azim,alt 
			// - calculate ra,dec - at time of transformation 
			// - calculate x,y from ra_tr,dec_tr
			double dec_tr,ra_tr,azim,alt;
			m_pAstroForms->calculateHorizontalCoordinatesFromEq( dec , ra, 
																			  m_pAsasTransform->m_TransformUtTime,
																			  alt, azim );		
			Calc_RaDec_at_T_transform( azim, alt, ra_tr, dec_tr );
			m_pAsasTransform->ad2xy_float( ra_tr, dec_tr, x , y );
		}else{
			// in tracking mode :
			m_pAsasTransform->ad2xy_float( ra, dec, x , y );
		}
		
		return TRUE;
	}
	return FALSE;	
}


BOOL_T CCDProcState::ad2xy( double ra, double dec,
									 int& x, int& y )
{
	if( gCCDParams.m_bUseAsasTransform && m_pAsasTransform && m_pAsasTransform->m_bTransformOK ){
		// first calculate horizontal
		/*double alt,azim;
		m_pAstroForms->calculateHorizontalCoordinatesFromEq( dec , ra, 
																			  m_pAsasTransform->m_TransformUtTime,
																			  alt, azim );


		// then calculate equatorial for new time 
		double ha_new,ra_new,dec_new;
		m_pAstroForms->calculateEquatorialCoordinates( azim, alt, ha_new, dec_new, ra_new );

		// and from ra_new,dec_new calculate x,y :
		m_pAsasTransform->ad2xy( ra_new, dec_new, x, y );*/

		if( m_pAstroForms->GetObsMode()==eNoMovingMode ){
			// IDEA :
			// - calculate azim,alt 
			// - calculate ra,dec - at time of transformation 
			// - calculate x,y from ra_tr,dec_tr
			double dec_tr,ra_tr,azim,alt;
			m_pAstroForms->calculateHorizontalCoordinatesFromEq( dec , ra, 
																			  m_pAsasTransform->m_TransformUtTime,
																			  alt, azim );		
			Calc_RaDec_at_T_transform( azim, alt, ra_tr, dec_tr );
			m_pAsasTransform->ad2xy( ra_tr, dec_tr, x , y );
		}else{
			// in tracking mode :
			m_pAsasTransform->ad2xy( ra, dec, x , y );
		}
		
		return TRUE;
	}
	return FALSE;	
}
									 

void CCDProcState::Calc_RaDec_at_T_transform( double azim, double alt,
                   				                double& ra, double& dec )
{
	double ha;
	m_pAstroForms->calculateEquatorialCoordinates( azim, alt, m_pAsasTransform->m_TransformUtTime, 
																  ha, dec, ra );	
}

BOOL_T CCDProcState::CalcAstroCoord( double x, double y,
                     		          double& azim, double& altit,
			                            double& dec, double& ra, double& ha )
{
	if( fabs(gCCDParams.m_BadAstroDX) > 0.1 ){
		double old_x =x;
		x = x + gCCDParams.m_BadAstroDX;		
//		printf("Shifted X=%.2f -> X=%.2f\n",old_x,x);
	}
	if( fabs(gCCDParams.m_BadAstroDY) > 0.1 ){
		double old_y=y;
		y = y + gCCDParams.m_BadAstroDY;
//		printf("Shifted Y=%.2f -> Y=%.2f\n",old_y,y);
	}

	
	// original :
	// m_pAstroForms->calculatePointEquatorial( x, y, dec, ra );
	// m_pAstroForms->calculateHourAngle( ra, dec, ha );
	// m_pAstroForms->calculateHorizontalCoordinatesBase( dec, ha, altit, azim );

	
	// test 3
	// m_pAstroForms->calculatePointEquatorialNEW( x, y, ra, dec );
	// m_pAstroForms->calculateHourAngle( ra, dec, ha );	
	// m_pAstroForms->calculateHorizontalCoordinatesBase( dec, ha, altit, azim );

	// test 2
	// m_pAstroForms->calculatePointEquatorial2( x, y, ra, dec, ha, altit, azim );
	// printf("ra=%s , dec=%.2f\n",AstroAngle::toString( ra, ANGLE_RA_TYPE ).c_str(), AstroAngle::rad2deg(dec) );
	
	// NEW pretty working version :
	if( gCCDParams.m_bUseAsasTransform && m_pAsasTransform && m_pAsasTransform->m_bTransformOK ){
		double ra_old,dec_old;
		m_pAsasTransform->xy2ad( x,y, ra_old, dec_old );		

		// changing ASAS ra - hours ,dec 0 degrees -> radians : 
		ra_old = AstroAngle::deg2rad( ra_old*15.00 ); // hour -> rad
		dec_old = AstroAngle::deg2rad( dec_old );     // deg -> rad

		if( m_pAstroForms->GetObsMode()==eNoMovingMode ){
			// in EARTH moving mode (mount is not moving )
			// we must calculated (azim,alt) in original time and 
			// then calculate ra,dec in current time :
			m_pAstroForms->calculateHorizontalCoordinatesFromEq( dec_old , ra_old, 
																				  m_pAsasTransform->m_TransformUtTime,
																				  altit, azim );
			double ha;
			m_pAstroForms->calculateEquatorialCoordinates( azim, altit, ha, dec, ra );
			//ra = AstroAngle::rad2deg( ra );
			//dec = AstroAngle::rad2deg( dec );
		}else{
			// eEarthMovingMode - mount is compensating for earth rotation :
			// simple case transformation does not change due to 
			// mount moving mode :
			ra = ra_old;
			dec = dec_old;

			// calculate azim and alt :
			m_pAstroForms->calculateHorizontalCoordinatesFromEq( dec , ra, 
																				  altit, azim );
		}
	}else{
		m_pAstroForms->calculatePointEquatorialTest( x,y, ra, dec, ha, altit, azim );
	}
	
	return TRUE;	
}

BOOL_T CCDProcState::CalcAstroCoord( double x, double y, time_t ut_time,
                     		          double& azim, double& altit,
			                            double& dec, double& ra, double& ha )
{
	m_pAstroForms->SetUT( ut_time );
	return CalcAstroCoord( x, y, azim, altit, dec, ra, ha );	
}

BOOL_T CCDProcState::TransformCCD1_to_CCD2( double x1, double y1, 
                                      		  double& x2, double& y2 )
{
	(gCCDParams.m_TransformMatrix).TimesVec( x1, y1, x2, y2 );
	return TRUE;
}


/*BOOL_T CCDProcState::TransformCCD1_to_CCD2( double x1, double y1, 
														  CCDParams* pParamCCD1, CCDParams* pParamCCD2,
														  double ccd2_shift_dx, double ccd2_shift_dy,
                                      		  double& x2, double& y2 )
{
	double sum = x1*pParamCCD1->m_TransformCCDPixelSize;
	double minus = ( 0.5*pParamCCD1->m_TransformCCDPixelSize*pParamCCD1->m_SizeX -
						  0.5*pParamCCD2->m_TransformCCDPixelSize*pParamCCD2->m_SizeX +
						  ccd2_shift_dx*pParamCCD2->m_TransformCCDPixelSize );
	double up = sum - minus;
	x2 = up/pParamCCD2->m_TransformCCDPixelSize;

	sum = y1*pParamCCD1->m_TransformCCDPixelSize;
	minus = ( 0.5*pParamCCD1->m_TransformCCDPixelSize*pParamCCD1->m_SizeY -
				 0.5*pParamCCD2->m_TransformCCDPixelSize*pParamCCD2->m_SizeY +
				 ccd2_shift_dy*pParamCCD2->m_TransformCCDPixelSize );
	up = sum - minus;
	y2 = up/pParamCCD2->m_TransformCCDPixelSize;


	BOOL_T bRet=FALSE;
	if( x2>=0 && x2<pParamCCD2->m_SizeX && y2>=0 && y2<=pParamCCD2->m_SizeY )
		bRet = TRUE;

	return bRet;					  
}
*/


void CCDProcState::xyAfterTime( double x, double y, double time_in_sec, double& x_new, double& y_new ){
	m_pAstroForms->XYAfterTimeTEST(x, y, x_new, y_new, time_in_sec);
}

void CCDProcState::xyAfterTimeNEW( double x, double y, double time_in_sec, double& x_new, double& y_new ){
	m_pAstroForms->XYAfterTimeTEST( x, y, x_new,  y_new, time_in_sec );
}

void CCDProcState::xyAfterTimeNEW( double x, double y, int* prevX, int* prevY, int backStart, int backTo,
                               const double* PrevFramesTime,int time_sign ){
	m_pAstroForms->XYAfterTimeTEST( x, y, prevX, prevY, backStart, backTo, PrevFramesTime, time_sign );
}


time_t CCDProcState::GetUtTime()
{
	return m_pAstroForms->GetUT();	
}