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
#ifndef _IMAGE_CREATOR_H__
#define _IMAGE_CREATOR_H__

#define NUM_POINTS 1000

#include "ccd_matrix.h"
#include <myfile.h>
#include <mycyclictab.h>
#include <vector>

using namespace std;

class mystring;
class cCCD;
class CccdSamples;
class CCDPipeline;
class CccdReport;

struct CSampleInfo 
{
	mystring fname;
	mystring magnitude;
};


class CSamplePutInfo
{
public:
	CSamplePutInfo( int x=0, int y=0, int frame_index=0, int cam_idx=0, const char* szMag="", 
						 double prevTime=0, int last_pos=0, int last_sample_idx=0 );
	CSamplePutInfo( const CSamplePutInfo& right );

		
	int m_X;
	int m_Y;
	int m_FrameIndex;
	int m_CamIdx;
	mystring m_szMag;		
	double m_prevTime;
	int m_LastPos;
	int m_LastSampleIdx;
};

class CSamplePutInfoTab : public vector<CSamplePutInfo>
{
public:
	CSamplePutInfoTab(){}
	CSamplePutInfoTab( const CSamplePutInfoTab& right ){
		(*this) = right;
	}
	CSamplePutInfoTab& operator=( const CSamplePutInfoTab& right ){
		clear();
		for(int i=0;i<right.size();i++){
			push_back( right[i] );
		}
		return (*this);
	}
};

class CImageCreator 
{
protected :
	double 	m_ObjProb;
	BOOL_T m_bFromSampleFile;
	mystring m_MinBright;
	mystring m_MaxBright;
	BOOL_T m_bInit;
	mystring m_szListname; 


	CListFile* FindListFile( const char* fname );

	LONG_T m_SampleIndex;
	
	
	CCycleTab2D<CSamplePutInfo>* m_pPrevPutObj;	
	CCDPipeline* m_pPipeline;
public :
	// where last object was put :
	int m_LastPutObjectX;
	int m_LastPutObjectY;
	double m_LastPutObjectRA;
	double m_LastPutObjectDEC;	
	int m_PrevValue;
	int m_LastPos;
	int m_LastSampleIdx;
	CccdReport* m_pGenEvent;

	vector<CListFile*> m_ListFiles;
	CccdSamples* m_pObjectSamples;
	
	void InitParams();
	void RefreshParams();	

	BOOL_T InitSampleCache();
	
	CImageCreator( CCDPipeline* pPipeline );
	~CImageCreator();

	// getting and setting params :
	const char* GetMinBright(){ return m_MinBright.c_str(); }
	const char* GetMaxBright(){ return m_MaxBright.c_str(); }
	void SetMinBright( const char* MinBright ){ m_MinBright=MinBright; }
	void SetMaxBright( const char* MaxBright ){ m_MaxBright=MaxBright; }
	void SetBright( const char* Bright ){
		SetMinBright( Bright );
		SetMaxBright( Bright );								
	}
		
	void RewindSamples(){
		m_SampleIndex = 0;			
	}

	void ModifyImage(CCDMatrix& Image);	
	
	// randomly put Mion signal on the image
	void AddMuon(CCDMatrix& Image);
	
	// adds some noise to image
	void AddNoise(CCDMatrix& Image);
	
	// adds new star to image
	void AddStar(CCDMatrix& Image,Table2D<ELEM_TYPE>& StarImage);
	
	// adds GRB 
	void AddGRBBurst(CCDMatrix& Image,Table2D<ELEM_TYPE>& GRBImage);

	// add objects : randomly adds objects, depends on ccd.cfg 
	// parameter CCD_USE_SAMPLE_OBJECTS
	void AddNextObject( cCCD& CCDFrame );

	void AddObjects( cCCD& CCDFrame ); 
	
	BOOL_T SimulateSuperNew( CCDMatrix& frame, int frame_index,
          				        const char* szMinMag, const char* szMaxMag,
									  double increase_ratio,
									  int& star_x, int& star_y,
									  double& initial_mag, double& res_mag );
	
	void AddObject( cCCD& CCDFrame, CCDPipeline& ccd_pipeline,
						 const char* fname, 
	                const char* szMag=NULL, BOOL_T bOnlyInAnalysied=FALSE );

	void AddObjectWithCoord( cCCD& CCDFrame, const char* fname, 
						 long& x, long& y,
	                const char* szMag=NULL, BOOL_T bOnlyInAnalysied=FALSE );

	
	void AddObject( CCDMatrix& Image );

	void ScaleSample( const char* szMagnitude, const Table2D<ELEM_SAMPLE_TYPE>& SampleObject, 
							Table2D<ELEM_SAMPLE_TYPE>& ScaledSample, const char* szInSampleMag );
	
	BOOL_T GetSampleObject( Table2D<ELEM_SAMPLE_TYPE>& SampleObject );
	BOOL_T GetSampleObject( Table2D<ELEM_SAMPLE_TYPE>& SampleObject, const char* fname, const char* szMagnitude );
	
	// new function - requires that sample was cut-out of same region 
	BOOL_T GetSampleObject_CheckXY( Table2D<ELEM_SAMPLE_TYPE>& SampleObject,
											  const char* fname, 
											  const char* szMagnitude, 
											  int x, int y, int nFrame=-1 );

	BOOL_T GetSampleObject_ByName(  Table2D<ELEM_SAMPLE_TYPE>& SampleIn,
											  Table2D<ELEM_SAMPLE_TYPE>& SampleOut,
											  const char* fname, 
											  const char* szMagnitude, 
											  int x, int y, int nFrame=-1 );
	
	void AddSampleObject( CCDMatrix& Image, Table2D<ELEM_SAMPLE_TYPE>& SampleObject, const char* szMag);

	void AddSampleObject( CCDMatrix& Image, LONG_T SampleIdx=-1 );
	
	void AddSampleObject( CCDMatrix& Image, Table2D<ELEM_SAMPLE_TYPE>& sample_img,
								 long& x, long& y, const char* szMag, BOOL_T bOnlyInAnalysiedRange,
								 int frame_index=0 );
								 
	void GetRandomPosition( CCDMatrix& Image, Table2D<ELEM_SAMPLE_TYPE>& sample_img,
				LONG_T& x, LONG_T& y, BOOL_T bOnlyInAnalysiedRange );

	void GetRandomPosition( CCDMatrix& Image, int sample_x_size, int sample_y_size,
				LONG_T& x, LONG_T& y, BOOL_T bOnlyInAnalysiedRange );
										 

	void AddSampleObjectAt( CCDMatrix& Image, Table2D<ELEM_SAMPLE_TYPE>& sample_img, long x, long y, 
											const char* szMag, int frame_index,
									BOOL_T bRePut=FALSE, int cam_idx=0,
									double ra_in_rad=0, double dec_in_rad=0 );
	

	void AddSampleObject( CCDMatrix& Image,
								 const char* fname, 
								 long& x, long& y,
								 const char* szMag, BOOL_T bOnlyInAnalysiedRange=FALSE,
								 int frame_index=0 );

	// creates simple simulation star, according to Gauss shape
	void InitSimpleStar(Table2D<ELEM_TYPE>& Image,LONG_T nPoints=NUM_POINTS,
	                    LONG_T max_value=255,LONG_T x=-1,LONG_T y=-1,LONG_T r=-1,
	                    LONG_T step=1);		
	

	// creates simple simulation star, according to Gauss shape
	BOOL_T AddSimpleStar(Table2D<ELEM_TYPE>& Image,double brightness,
			  				   LONG_T width,LONG_T x0,LONG_T y0);		

	BOOL_T AddSimpleStar(Table2D<ELEM_TYPE>& Image,double brightness,LONG_T width);


	void RePutPrevImages( CCDPipeline& ccd_pipeline );

	static double get_gauss( double s0, double sigma );
	
	BOOL_T IsGenEvent( CPoint& evt );


	static void LogSamplePut( Table2D<ELEM_SAMPLE_TYPE>* sample, int cam, int run,
	                          int frame, int x, int y, const char* szError=NULL );
protected :
	// returns file name of sample of given brightness :
	mystring GetRandomSampleFile( const char* szMag, LONG_T SampleIdx=-1  );

	

	// not very usful :
	LONG_T GetGauss(LONG_T x0,LONG_T y0,LONG_T r0,LONG_T x,LONG_T y,LONG_T max_value=1);	
};



#endif
