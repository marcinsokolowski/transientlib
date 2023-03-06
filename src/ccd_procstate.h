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
#ifndef _CCD_PROC_STATE_H__
#define _CCD_PROC_STATE_H__


#include "ccd_defines.h"

class AstroCCD;
class CCcdCfg;
class CCDParams;
class CCDPipeline;
class CCDAsasTransform;

class CCDProcState
{
public:
	CCDPipeline* m_pPipeline;
	CCDProcState( CCDPipeline* pPipeline );
	CCDProcState( CCDAsasTransform* pAsasTransform );
	CCDProcState( CCDAsasTransform* pAsasTransform, eObservationMode_T obsMode );		
	~CCDProcState();


	static BOOL_T TransformCCD1_to_CCD2( double x1, double y1,
	         			                   double& x2, double& y2 );	                              
	                                                                                                                                    
	
	void SetProcState( ePROCSTATE_T eProcState );
	void ResetProcState();
	
	ePROCSTATE_T m_eProcState;


	double AlfaMultPosTab[MAX_PIPELINE_SIZE];
	double AlfaMultNegTab[MAX_PIPELINE_SIZE];
	void ReCalcAlfaMultTab( double alfa, int max_times );


	void InitAstroForms();
	void InitAstroForms( CCDAsasTransform* pAsasTransform, eObservationMode_T obsMode );
	void InitAstroForms( CCDAsasTransform* pAsasTransform );
	void InitAstroForms( double focus, double pixel_size_to_set, int x_size, int y_size,
								double geoLat, double geoLong,
			 			      double geoAzim, double geoAlt, double orientation, 
			 			      int time_zone_to_set, 
			 			      int _X_orient, int _Y_orient,
			 			      eObservationMode_T obsMode, double dec_obs, double ra_obs,
			 			      double ut_to_set, int day_to_set, int month_to_set, 
						      int year_to_set, CCDAsasTransform* pAsasTransform=NULL );

	void GetObsCoo( time_t ut_time, double& ra, double& dec, double& alt, 
						 double& azim, double& ha, eObservationMode_T& obsMode );
	void SetFrameDateTime( time_t ut_time, CCDPipeline* pPipeline );

	void xyAfterTime( double x, double y, double time_in_sec, double& x_new, double& y_new );
	
	void xyAfterTimeNEW( double x, double y, double time_in_sec, double& x_new, double& y_new );

	void xyAfterTimeNEW( double x, double y, int* prevX, int* prevY, int backStart, int backTo,
								 const double* PrevFramesTime,int time_sign );

	BOOL_T CalcAzimutalCoord( int x, int y, time_t ut_time,
									  double& azim, double& altit );

	BOOL_T ad2xy( double ra, double dec, int& x, int& y );	                            
	BOOL_T ad2xy_float( double ra, double dec, double& x, double& y );	                            
	
	void ad2azh( double ra_in_rad, double dec_in_rad, time_t unix_time,
					 double& az_in_rad, double& h_in_rad );
					 
	void azh2ad( double az_in_rad, double h_in_rad, time_t unix_time,
					 double& ra_in_rad, double& dec_in_ra );

	void Calc_RaDec_at_T_transform( double azim, double alt,
											  double& ra, double& dec );

	BOOL_T CalcAstroCoord( double x, double y, time_t ut_time,
								  double& azim, double& altit,
								  double& dec, double& ra, double& ha );

	BOOL_T CalcAstroCoord( double x, double y,
								  double& azim, double& altit,
								  double& dec, double& ra, double& ha );
								  
	time_t GetUtTime();								  


	// obiects for calculationing posisitons :
	AstroCCD* m_pAstroForms;		


	// asas transformation for given camera :
	CCDAsasTransform* m_pAsasTransform;
};



#endif
