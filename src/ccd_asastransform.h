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
#ifndef _CCD_ASAS_TRANSFORM_H__
#define _CCD_ASAS_TRANSFORM_H__

#include <mytypes.h>
#include <mystring.h>
#include <asas_gparam_def.h>
#include "ccd_defines.h"

// #define ASAS_TRANSFORM_PARAM_COUNT 21
class CCDConfig;
class CCDPipeline;

class CCDAsasTransform : public GPARAM
{
public :
	time_t m_TransformUtTime;   // it is set before astrometry - MAYBE NOT GOOD idea 
										 // but I am afraid to change it now ...
										 // So indicates time of last ASTRO PERFORMED
										 // but NOT GOOD ONE !!!
										 
	time_t m_LastGoodAstroTime; // indicating last good astrometry time 
										 // this is just for send_pos_to_mount_astro
										 // safety reasons - NOT TO USE !
										 
	int m_SizeX;
	int m_SizeY;
	CCDConfig* m_pCfg;
	BOOL_T m_bTransformOK;
	BOOL_T m_bReadDone;
	double m_AzimInRad;
	double m_AltInRad;
	CCDPipeline* m_pPipeline;
	int m_nLastFailedCount;
	int m_nFailedAstrometryCount; // number of failed astrometries
											// only those realy done are counted
											// if skiped due to too low #stars not increased
	int m_LastAstroRetCode;											
	

	CCDAsasTransform( CCDConfig* pCfg, CCDPipeline* pPipeline=NULL );	

	// function for calling ASAS photometry :
	BOOL_T asas_photometry( const char* fits_file, double threshold, mystring& szMagFile );

	// function for calling ASAS astrometry :
	BOOL_T asas_astrometry( const char* mag_file, const char* ast_file,
								 double ra0, double dec0, double fi, double pixscale,
								int ord,
								const char* star_cat_path,
								int verb=0, int _try=0,
								double err=0.3, double errf=0.6 );
									
	// transformations x,y <-> ra,dec
	// output in ASAS hours,dec
	BOOL_T xy2ad( double x, double y, double& ra, double& dec );

	// takes input in radians :
	BOOL_T ad2xy( double ra, double dec, int& x, int& y );
	BOOL_T ad2xy_float( double ra, double dec, double& x, double& y );
	
	// takes hours,degrees, returns int
	void ad2xy_raw( double ra, double dec, int& x, int& y );

	BOOL_T SaveToFile( const char* szFileName, eObservationMode_T obsMode );

	BOOL_T ReadFromFile( const char* szFileName , BOOL_T bInitTransform=FALSE );

	void ResetAstroOK();

	void Dump();
	
	void InitSimple( double _ra, double _dec, double _fi=-1000.00, double _pixscale=-1000.00 );
	
	
	BOOL_T WasAstroCanceled();

	static void ForceAstrometryBreak();
	static int AstroBreakForced();
};


#endif
