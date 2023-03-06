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
#ifndef _CCD_STAR_SPY_H__
#define _CCD_STAR_SPY_H__

#include <mytypes.h>
#include <myhisto.h>
#include "ccd_defines.h"
#include "mypoints.h"
#include <vector>

class CCDParams;
class CCDPipeline;
class CCDConfig;

struct CStarDesc
{
	CStarDesc( int low_x=0, int low_y=0, int up_x=0, int up_y=0 );

	CStarDesc& operator=(const CStarDesc& right);

	int m_LowX;
	int m_LowY;
	int m_UpX;
	int m_UpY;

	double m_StarRA;
	double m_StarDec;
	
	double m_StartX;
	double m_StartY;
	
	double m_CurrX;
	double m_CurrY;
	
	double m_PredictedX;
	double m_PredictedY;


	double m_PrevStepX;
	double m_PrevStepY;
	double m_AverageStepX;	
	double m_AverageStepY;

	double m_PrevMAX;
	
	
	double m_MeanValue;
	double m_RMS;

	double m_MaxDX;
	double m_MaxDY;
		
	CPointList m_PosTab;
	
	BOOL_T m_bNotEnd;
	
   double m_starDX;
   double m_starDY;
   double m_starDXPerSec;
   double m_starDYPerSec;
                     


	double m_OrtoLineA;
	double m_OrtoLineB;		
	double m_OrtoLineC;		
	
	BOOL_T m_bGood;
	
	int m_nCurrentStep;
	
	CMyHisto m_StarValues;
	
	int m_nDrasticChangeCount;
	BOOL_T m_bIdentOK;

	CPoint* FindFrame( int frame_index );
	int GetAverageShifts( double& aver_dx, double& aver_dy, int back );
	BOOL_T InitMAXStar( ELEM_TYPE** p_data, double startTime, int frame_index, int SizeX,
							  int SizeY, CStarDesc* pOtherStarsTab, int otherCount, 
							  CCDPipeline* pPipeline );
	BOOL_T IsDistFromOthersOK( int x, int y, CStarDesc* pOtherStarsTab, int otherCount );							  		

	double CalcAverageDX( int new_x , int nPrev );
	double CalcAverageDY( int new_y , int nPrev );		
	
	double UpdateRMS();
	
	int GetFramesCount(){ return m_PosTab.size(); }
};

class CCDStarSpy
{
public:
	CCDPipeline* m_pPipeline;
	CStarDesc* m_StarDesc;
	int m_StarsCount;
	int m_nSize;

	int m_SizeX;
	int m_SizeY;

	double m_StartTime;
   double m_EndTime;
   int m_nFramesCount;
   

	BOOL_T m_bFromFile;
	BOOL_T m_bRotation;
	double m_RotCenterX;
	double m_RotCenterY;
	double m_dAlfa;
	double m_dAlfaPerSec;
	double m_DX;
	double m_DY;
	double m_DXPerSec;
	double m_DYPerSec;		


	CCDStarSpy( int sizeX, int sizeY, int nSize=1, CCDPipeline* pPipeline=NULL );
	CCDStarSpy( int sizeX, int sizeY, double rotX, double rotY, double dAlfa, double dx, double dy, 
					BOOL_T bRot, double dx_per_sec, double dy_per_sec, double dAlfaPerSec,
					CCDPipeline* pPipeline=NULL );
	CCDStarSpy( const CCDStarSpy& right );
	CCDStarSpy& operator=(const CCDStarSpy& right );
	~CCDStarSpy();

	CStarDesc& GetStarDesc( int x , int y );

	void ReInitMAXStarIfNeeded( ELEM_TYPE** p_data, double startTime, int frame_index, BOOL_T bForceReInit=FALSE );	
	void InitWithSpecificStar( ELEM_TYPE** p_data, double startTime, int starX0, int starY0, int frame_index );
	void InitWithMAXStar( ELEM_TYPE** p_data, double startTime=0, int frame_index=0 );

	void LogMAXStarsPositions( CCDPipeline* pPipeline );
	BOOL_T ReCalcNewPosition( ELEM_TYPE** p_data, int frame_time, int frame_index,
									  int searchSize=10, double sizeTolerance=1.0, BOOL_T bUseHint=FALSE,
									  int reCalcSearchSize=5,									  
									  BOOL_T bUseSteps=FALSE, double dx=0, double dy=0,
									  BOOL_T* reInit=NULL );	
									  
	BOOL_T ReCalcNewPositionByRaDec( int frame_time, CCDPipeline* pPipeline );									  

	BOOL_T CheckIfRotation();

	BOOL_T CalcRotation( double& angle_per_frame,  double& x_cross, double& y_cross,
								double& frameDX, double& frameDY, 
							   double& frameDXPerSec, double& frameDYPerSec, BOOL_T& bRotation,
							   double& angle_per_sec,
							   double endTime=0, int min_steps=5, int min_steps_for_single=20 );
								
	BOOL_T SaveToFile( const char* szFileName, int CameraIndex=0, 
							 const char* load_file=NULL, CCDPipeline* pPipeline=NULL );

	void SaveStarsTracks( const char* fname="starstrack.trace", int cam_no=0, 
								 BOOL_T bSaveTrack=FALSE, CCDPipeline* pPipeline=NULL );
								 

	BOOL_T CheckIfHaveTwoWithFrame( int frame1, int frame2 );								 
	BOOL_T FindTransform( int frame1, int frame2, CCDConfig& params, int ccd_index, CCDPipeline* pPipeline );

	void GetOutFileName( mystring& szOutName, const char* szBaseName, int ccd_index, CCDPipeline* pPipeline );	
	void DumpToFile( const char* szReport, const char* szFileName, int ccd_index, CCDPipeline* pPipeline );

	// 
	void SaveCurrentPosition( const char* fname );		
	void LogCurrentPosition(  const char* szBaseName, int ccd_index, CCDPipeline* pPipeline );

	// prediction :
	BOOL_T CalcPredicted( CCDPipeline* pPipeline, double curr_time, int ccd_index );

	// returns number of frames if best star ( the one with 
	// biggest number of frames on which was traced ) : 
	int GetFramesCount();
	
	void ResetTracedStar();
	
	// initializes borders of traced star areas :
	void InitAreas( int BorderSize=50 );
};


class CCDStarSpyTab : public vector<CCDStarSpy>
{
public :
	BOOL_T m_bFromFile;
	mystring m_szName;
	CCDPipeline* m_pPipeline;

	CCDStarSpyTab( const char* szName="spystars");
	
	void SetPipelinePtr( CCDPipeline* pPipeline ){ m_pPipeline=pPipeline; }
	
	BOOL_T SaveToFile( const char* szFileName, const char* load_file=NULL, 
							 CCDPipeline* pPipeline=NULL );
	BOOL_T ReadFromFile(	const char* szFileName );
	void DumpReportToFile();

	void MarkToReInit();

	// clears all traced stars so that they are re-initialized in next step :
	void ResetTracedStars();
};

#endif
