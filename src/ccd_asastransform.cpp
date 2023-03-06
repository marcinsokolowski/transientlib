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

// GLOBAL VARIABLE TO BE SURE THAT ASTROMETRY IS IMPLEMENTED 
// in order to try to run without astrometry set it to 0 
int gExitDueToAstrometryNotImpl = 0;

#include "ccd_asastransform.h"
#include <mylock.h>
#include <mystring.h>
#include <myfile.h>
#include "ccd_globals.h"
#include "ccd_pipeline.h"
#include <mymacros.h>
#include <Astroangle.h>
#include <AstroCCD.h>
#include <Astroangle.h>
#include <myutil.h>
#include <fits_file.h>
                                                                                
extern "C" {
// #include <ccdphot.h>
#include <asas_astrometry.h>
}
                                                                                



// lock to protect global resources in ASASlib :
CMyMutex gAsasAstroLock;


CCDAsasTransform::CCDAsasTransform( CCDConfig* pCfg, CCDPipeline* pPipeline )
: m_bTransformOK(FALSE), m_pCfg( pCfg ), m_SizeX(0), m_SizeY(0),
  m_TransformUtTime(0),m_bReadDone(FALSE), m_AzimInRad(0), m_AltInRad(0),
  m_pPipeline( pPipeline ), m_LastGoodAstroTime(0), m_nLastFailedCount(0),
  m_nFailedAstrometryCount(0),m_LastAstroRetCode(0)
{
	for(register int i=0;i<ASAS_TRANSFORM_PARAM_COUNT;i++){
		px[i] = 0;
		py[i] = 0;
	}
	order = 4;
	pixscale = 36.0;
	ra = 0.00;
	dec = 0.00;
	xc = (2062/2);
	yc = (2048/2);
	ast_err = 0.00;
	fi = 0.00;

	if( pCfg ){
	   flip=pCfg->m_eReverseForTransform;
	}else{
		flip = eReverseImageVert;
	}
}


BOOL_T CCDAsasTransform::xy2ad( double x, double y, double& ra, double& dec )
{
	if( m_bTransformOK ){
		// transform x,y -> a,d :
		double x_ft = x;
		double y_ft = y;	

		/*_TRACE_PRINTF_8("before flip : xy2ad( %.2f , %.2f )\n",x, y);
		if( m_pCfg->m_eReverseForTransform!=eReverseImageNone ){
			if( m_pCfg->m_eReverseForTransform==eReverseImageHor || m_pCfg->m_eReverseForTransform==eReverseImageHorVert )
			{
				x_ft = (m_SizeX-x_ft);
			}
			if( m_pCfg->m_eReverseForTransform==eReverseImageVert || m_pCfg->m_eReverseForTransform==eReverseImageHorVert )
			{
				y_ft = (m_SizeY-y_ft);
			}
		}

		_TRACE_PRINTF_8("after flip : xy2ad2( %.2f , %.2f )\n",x_ft, y_ft);*/

		// IMPORTANT NOTE :
		// ra is returned in HOURS !!!!!
		// dec in DEG !!!!!!!
//ASAS code not available :
//		::xy2ad( (struct GPARAM*)this, x_ft, y_ft, &ra, &dec );
// PUT YOUR CODE HERE :
		printf("\n\nERROR : CCDAsasTransform::xy2ad should be filled with proper body of your astrometric code!!!\n");
		
		if( gExitDueToAstrometryNotImpl ){
			printf("ERROR : astrometry not implemented gExitDueToAstrometryNotImpl>0 -> exiting now !\n");
			exit(-1);
		}

		return FALSE;


		return TRUE;
	}
	return FALSE;		
}

void CCDAsasTransform::ad2xy_raw( double ra, double dec, int& x, int& y )
{
	double x_ft,y_ft;

// PUT YOUR CODE BELOW : 
//	::ad2xy( (struct GPARAM*)this, ra, dec, &x_ft, &y_ft );	
// PUT YOUR CODE ABOVE

	x = my_round(x_ft);
	y = my_round(y_ft);
}

BOOL_T CCDAsasTransform::ad2xy( double ra, double dec, int& x, int& y )
{
	if( m_bTransformOK ){
		// transform x,y -> a,d :
		double x_ft,y_ft;

		ra = AstroAngle::rad2hours( ra );
		dec = AstroAngle::rad2deg( dec );


// PUT YOUR CODE BELOW :
//ASAS code not available :	
//		::ad2xy( (struct GPARAM*)this, ra, dec, &x_ft, &y_ft );
// PUT YOUR CODE HERE : 
		printf("\n\nERROR : CCDAsasTransform::ad2xy should be filled with proper body of your astrometric code!!!\n");

		if( gExitDueToAstrometryNotImpl ){
			printf("ERROR : astrometry not implemented gExitDueToAstrometryNotImpl>0 -> exiting now !\n");
			exit(-1);
		}
		return FALSE;
// PUT YOUR CODE ABOVE
		
		/*x = my_round(x_ft);
		y = my_round(y_ft);
		if( m_pCfg->m_eReverseForTransform!=eReverseImageNone ){
			if( m_pCfg->m_eReverseForTransform==eReverseImageHor || m_pCfg->m_eReverseForTransform==eReverseImageHorVert )
			{
				x = my_round(m_SizeX-x_ft);
			}
			if( m_pCfg->m_eReverseForTransform==eReverseImageVert || m_pCfg->m_eReverseForTransform==eReverseImageHorVert )
			{
				y = my_round(m_SizeY-y_ft);
			}
		}*/
		// ra is returned in HOURS !!!!!
		// dec in DEG !!!!!!!
		return TRUE;
	}
	return FALSE;		
}

BOOL_T CCDAsasTransform::ad2xy_float( double ra, double dec, double& x, double& y )
{
	if( m_bTransformOK ){
		// transform x,y -> a,d :

		ra = AstroAngle::rad2hours( ra );
		dec = AstroAngle::rad2deg( dec );


// PUT YOUR CODE BELOW :
//ASAS code not available :	
		printf("\n\nERROR : CCDAsasTransform::ad2xy_float should be filled with proper body of your astrometric code!!!\n");
		if( gExitDueToAstrometryNotImpl ){
			printf("ERROR : astrometry not implemented gExitDueToAstrometryNotImpl>0 -> exiting now !\n");
			exit(-1);
		}
		return FALSE;
// PUT YOUR CODE ABOVE
		
		if( m_pCfg->m_eReverseForTransform!=eReverseImageNone ){
			if( m_pCfg->m_eReverseForTransform==eReverseImageHor || m_pCfg->m_eReverseForTransform==eReverseImageHorVert )
			{
				x = m_SizeX-x;
			}
			if( m_pCfg->m_eReverseForTransform==eReverseImageVert || m_pCfg->m_eReverseForTransform==eReverseImageHorVert )
			{
				y = m_SizeY-y;
			}
		}
		// ra is returned in HOURS !!!!!
		// dec in DEG !!!!!!!
		return TRUE;
	}
	return FALSE;		
}


BOOL_T CCDAsasTransform::asas_photometry( const char* fits_file, double threshold, mystring& szMagFile )
{
	int flags=0;
	mystring szFile = fits_file;
	szMagFile = "";
	szMagFile << getfname( szFile ) << ".mag";
	char sMagFile[1024],szFitsFile[1024];
	strcpy(sMagFile,szMagFile.c_str());
	strcpy(szFitsFile,fits_file);

	gAsasAstroLock.Lock();
	printf("ERROR : function ::phot ( asaslib ) is not available in this version\n");
//ASAS not available 
	printf("ERROR : ASAS version of ::phot is not available in this version\n");
//	int ret = ::phot( szFitsFile, sMagFile, NULL, threshold, NULL, 0, flags, NULL );
	int ret = 0;
	gAsasAstroLock.UnLock();
	return (ret==0);
}


BOOL_T CCDAsasTransform::asas_astrometry( const char* mag_file, 
													const char* ast_file,
												   double ra0, double dec0, 
													double fi, double pixscale,
												   int ord,
											      const char* star_cat_path, 
                                       int verb, int _try,
												   double err, double errf )
{
	char sMagFile[1024],sAstFile[1024],sStarCat[1024];
	strcpy(sMagFile,mag_file);
	strcpy(sAstFile,ast_file);
	strcpy(sStarCat,star_cat_path);

	gAsasAstroLock.Lock();
	double xoff=0,yoff=0;
	this->order = ord;
	MAX_AST_ERR_FATAL = errf;
   MAX_AST_ERR = err;
	

	// reseting OK flag before astrometry :
	// IT MUST ALWAYS BE LIKE THIS !!!!!!!!!!!!!!!!
	// 
	// due to fact that when we go to new field and astrometry is performed
	// flag must be FALSE - so that potential send_pos_to_mount 
	// will not be performed !!!!!!!!!
	// it would make big mess - would sent OLD position to DAQ :
	m_bTransformOK = FALSE;

	// experimental - WORKAROUND of bug inside astrometry 
	// it does not work correctly for RA ~= 0, instead 24 must be used :
	// workaround was moved to function ::astrometry()
/*	if( ra0>=0.0 && ra0<0.1 ){
		printf("WARNING : ra<0.1 , changing to ra=24 !!!!! - check if ok !!!\n");fflush(0);
		ra0 = 24.00;
	}*/
	
	time_t t1=get_dttm();
//	ASAS code not available :
	int ret = 0;
	printf("ERROR : asas astrometry cannot be used in this version !!!\n");
/*	int ret = ::astrometry( sMagFile, sAstFile, ra0, dec0, fi, pixscale,
									xoff, yoff,
									0, 0, 1, 1, (struct GPARAM*)this, 
									sStarCat, verb, _try, 0 );*/
	
// PUT YOUR CODE BELOW :
	mystring szMagFile = mag_file,szExt;
	mystring szBaseName = getfname( szMagFile );
	char szFitsFile[512];
	int len = strlen(szBaseName.c_str())-4;
	strncpy(szFitsFile,szBaseName.c_str(),len);
	szFitsFile[len] = '\0';
	mystring szFITS;
	szFITS << szFitsFile << ".fit";

	printf("WARNING : put hear reading of astrometric paramters for file %s / %s\n",szFitsFile,szFITS.c_str());	
	printf("\n\nERROR : CCDAsasTransform::asas_astrometry should be filled with proper body of your astrometric code!!!\n");
	printf("ERROR : put reading of astrometric transformation of file %s here !!!\n",mag_file);
	if( gExitDueToAstrometryNotImpl ){
		printf("ERROR : astrometry not implemented gExitDueToAstrometryNotImpl>0 -> exiting now !\n");fflush(stdout);
		printf("EXIT !!!\n");
		exit(-1);
	}
	ret = -1000;
// PUT YOUR CODE ABOVE

	time_t t2=get_dttm();
	printf("ASTROMETRY_TOOK : %d sec (ret=%d)\n",(t2-t1),ret);

	if( ret == 0 ){
		printf("Setting final coord from GOOD astrometry (ra,dec)=(%.2f,%.2f)\n",gOutRA2000,gOutDEC2000);
		ra = gOutRA2000;
		dec = gOutDEC2000;
	}
	m_LastAstroRetCode = ret;
	gAsasAstroLock.UnLock();	
	return ( ret == 0 );
}

BOOL_T CCDAsasTransform::SaveToFile( const char* szFileName, eObservationMode_T obsMode )
{
	MyOFile out( szFileName );
	out.Printf("POSANGLE=%.8f\n",fi);
	out.Printf("AST_ORD=%d\n",order);
	for(int i=0;i<14;i++){
		out.Printf("PAR_X_%d=%.8f\n",i,px[i]);
	}
	for(int i=0;i<14;i++){
		out.Printf("PAR_Y_%d=%.8f\n",i,py[i]);
	}
	out.Printf("PIXSCALE=%.8f\n",pixscale);
	out.Printf("RA2000=%.8f\n",gOutRA2000);
	out.Printf("DEC2000=%.8f\n",gOutDEC2000);
	out.Printf("CCD_TRANSFORM_UT_TIME=%d\n",m_TransformUtTime);

	/*AstroCCD astro( gCCDParams.m_TransformCCDFocus, gCCDParams.m_TransformCCDPixelSize,
						 gCCDParams.m_SizeX, gCCDParams.m_SizeY,
						 gCCDParams.m_GeoLatitude, gCCDParams.m_GeoLongitude,
						 0, 0, 0 );*/
	double dec_in_rad, ra_in_rad;
	dec_in_rad = AstroAngle::deg2rad( gOutDEC2000 );
	ra_in_rad = AstroAngle::hours2rad( gOutRA2000 );

	if( m_pPipeline ){
		// using astro object in pipeline :
		m_pPipeline->GetAstroCalcObj()->calculateHorizontalCoordinatesFromEq(
			dec_in_rad, ra_in_rad, m_TransformUtTime, m_AltInRad, m_AzimInRad );
	}else{
		// printf("\n\n\nERROR IN CODE !!!!!!!!!!!!!!!!!!\n");
		_TRACE_PRINTF_3("pointer to Pipeline not set using static AstroCCD::calculateHorizontalCoordinatesFromEq\n");
		AstroCCD::calculateHorizontalCoordinatesFromEq( dec_in_rad, ra_in_rad,
															 m_TransformUtTime, 
															 gCCDParams.m_GeoLongitude,
															 gCCDParams.m_GeoLatitude,
															 m_AltInRad, m_AzimInRad );
	}
	out.Printf("AZIM_CENTER=%.8f\n",AstroAngle::rad2deg( m_AzimInRad ));
	out.Printf("ALT_CENTER=%.8f\n",AstroAngle::rad2deg( m_AltInRad ));
	out.Printf("CCD_OBS_MODE=%d\n",obsMode);

	return TRUE;
}

BOOL_T CCDAsasTransform::ReadFromFile( const char* szFileName , BOOL_T bInitTransform )
{
	BOOL_T bRet=TRUE;
	CSafeKeyTab keyTab;
   CFITSFile<ELEM_TYPE> in;
   if(!in.ReadFITSHeader( keyTab, szFileName )){
      printf("could not read FITS-Header file : %s\n", szFileName );
		return FALSE;
   }

	m_TransformUtTime = (int)CCDMatrix::getObsTime( keyTab, TRUE );
	if( m_TransformUtTime>0 )
		bRet = TRUE;

	if( bInitTransform ){
		const char* szVal = keyTab.getKeyVal( POSANGLE );
		if( szVal && szVal[0] ){
			fi = atof( szVal );
		}else{
			return FALSE;
		}

		szVal = keyTab.getKeyVal( AST_ORD );
		if( szVal && szVal[0] ){
			order = atol( szVal );
		}else{
			return FALSE;
		}
		
		char szTmpKey[64];
		for(int i=0;i<14;i++){
			sprintf(szTmpKey,"PAR_X_%d",i);
			szVal = keyTab.getKeyVal( szTmpKey );	
			if ( szVal && szVal[0] ){
				px[i] = atof( szVal );
			}else{
				px[i] = 0.00;
				return FALSE;
			}
		}
		for(int i=0;i<14;i++){
			sprintf(szTmpKey,"PAR_Y_%d",i);
			szVal = keyTab.getKeyVal( szTmpKey );	
			if ( szVal && szVal[0] ){
				py[i] = atof( szVal );	
			}else{
				py[i] = 0.00;
				return FALSE;
			}
		}

		szVal = keyTab.getKeyVal( PIXSCALE );
		if( szVal && szVal[0] ){
			pixscale = atof( szVal );
		}else{
			return FALSE;
		}

		szVal = keyTab.getKeyVal( RA_OBS );
		if( szVal && szVal[0] ){
         ra = atof( szVal );
      }else{
			return FALSE;
		}

		szVal = keyTab.getKeyVal( DEC_OBS );
		if( szVal && szVal[0] ){
         dec = atof( szVal );
      }else{
			return FALSE;
		}
	}

	return bRet;
}

void CCDAsasTransform::ResetAstroOK()
{
	printf("Reseting astrometry OK flag\n");fflush(0);
	m_bTransformOK = FALSE;
//	m_nFailedAstrometryCount = 0;
}

void CCDAsasTransform::Dump()
{
	printf("ut_time = %d fi=%.2f pixscale=%.2f ra=%.4f dec=%.4f order=%d\n",
		m_TransformUtTime,fi,pixscale,ra,dec,order);
	for(int i=0;i<14;i++){
		printf("%.2f ",px[i]);
	}
	printf("\n");
	for(int i=0;i<14;i++){
		printf("%.2f ",py[i]);
	}
	printf("\n");	
}

void CCDAsasTransform::InitSimple( double _ra, double _dec, double _fi, double _pixscale )
{
	ra = _ra;
	dec = _dec;

	int i=0;
	for(i=0; i<21; i++){
		px[i]=0.00;
		py[i]=0.00;
	}
	if( _pixscale > 0 ){
		pixscale = _pixscale;
	}else{
		pixscale = gCCDParams.m_fPixScale;
	}
	fi = 0.00;
	if( _fi > -1000 ){
		fi = _fi;
	}
	double rd = fi*(PI_VALUE/180.00);

	px[0] = py[0] = 0.;
	px[1] = pixscale * cos(rd);
	px[2] = -pixscale * sin(rd);
	py[1] = pixscale * sin(rd);
	py[2] = pixscale * cos(rd);
	for(i=3; i<21; i++){
		px[i]=0.;
		py[i]=0.;
	}
	
}

void CCDAsasTransform::ForceAstrometryBreak()
{
	gForceAstrometryBreak = 1;
}

int CCDAsasTransform::AstroBreakForced()
{
	return gForceAstrometryBreak;
}

BOOL_T CCDAsasTransform::WasAstroCanceled()
{ 
	return (m_LastAstroRetCode == E_stop_forced_externally);
}


