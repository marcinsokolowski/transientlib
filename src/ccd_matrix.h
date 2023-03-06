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
#ifndef _CCD_MATRIX_H__
#define _CCD_MATRIX_H__

#include <tab2D.h>
#include <mypoints.h>
#include <ccdcfg.h>
#include <basedefines.h>
#include "ccd_defines.h"
#include "ccd_report.h"
#include <ccd_fits_header_defs.h>
#include <ccd_hardware_defines.h>

class CKeyTab;
class cCCD;
class InfoTable2D;
class CCcdCfg;
class CCDPipeline;
class CMyStrTable;
class CSafeKeyTab;

// HINTS :
// in order to auto-alloc , do :
// CCDMatrix ccd_image( 0 , 0 , FALSE, NOT_DEFINED, NULL, 0, FALSE, NULL, TRUE)
// ccd_image.ReadFITSFile( szFile.c_str() , FALSE );
//
//
// In order to force automatic usage of global x_size,y_size parameter call
// construktor with CCDMatrix matrix(0,0);
class CCDMatrix : public Table2D<ELEM_TYPE>
{
protected:
	// list of coordinates on the frame containing potentially
	// intersting events :
	// found - real events are stored here :
	CCDEventList m_FoundEventList;

	// all generated events - also found 
	CCDEventList m_GeneratedEventList;
	
	// list of all identified events - no matter generated or not :
	// all found events - no metter generated or real here :
	// CCDEventList m_IdentifiedEvents;

	BOOL_T m_bIntersting;
	
	
	char* m_pHeaderPtr;
	
	// functions :	
public:
	double m_dx; // dx shift from previous image to current 	
	double m_dy; // dy shift from previous image to current

	const char* GetFileName(){ return m_szFileName.c_str(); }

	BIG_ELEM_TYPE* get_laplace_data()
	{
		if( m_pFrameLaplace ){
			return m_pFrameLaplace->get_data_buffer();
		}
		return NULL;
	}

	BIG_ELEM_TYPE** get_laplace_data_fast() 
	{
		if( m_pFrameLaplace ){
			return m_pFrameLaplace->get_data_buffer_fast();
		}
		return NULL;
	}

	void SetHeaderPtr(  char* pHeaderPtr ) { m_pHeaderPtr=pHeaderPtr; }
	char* GetHeaderPtr() { return m_pHeaderPtr; }

	static double getObsTime( CSafeKeyTab& keys, BOOL_T bIgnoreError );
	double getObsTime( BOOL_T bIgnoreError=FALSE );
	BOOL_T getObsTime( struct tm* _tm );
	time_t getMinTime();
	time_t getMaxTime();		
	BOOL_T IsDark();

	virtual eTableType GetType(){ return eCCDMatrix; }
	
	CCDMatrix();
	CCDMatrix(long x_size,long y_size,BOOL_T bAutoInit=TRUE,
				 LONG_T idx=NOT_DEFINED,cCCD* pFrame=NULL,int nAllocEvents=0,
				 BOOL_T bKeepLaplaceOfCurrent=FALSE, CCcdCfg* pParamsSet=NULL,
				 BOOL_T bAllocHere=TRUE);
	~CCDMatrix();
	
	BOOL_T CCDMatrix_InitConstructor(long x_size=0,long y_size=0,BOOL_T bAutoInit=TRUE,
	             LONG_T idx=NOT_DEFINED,cCCD* pFrame=NULL,int nAllocEvents=0,
	             BOOL_T bKeepLaplaceOfCurrent=FALSE, CCcdCfg* pParamsSet=NULL,
	             BOOL_T bAllocHere=TRUE);

	//-----------------------------------------------------------------------
	// ANALYSING FUNCTIONS 
	//-----------------------------------------------------------------------
	int Subtract(CCDMatrix& right,CCDMatrix& result,BOOL_T bZero=FALSE);
	int Subtract(CCDMatrix& right,BOOL_T bZero=FALSE);
	void Divide(CCDMatrix& right,CCDMatrix& result,int max_val);
	void CatBorder( int border, BOOL_T bPutMean );
	void CatBorderFast( int border );		
	
	void CutBorder( int left, int up, int right, int bottom, 
						 CCDMatrix& out_image );
	
	
	void AdjustFrame();
	void Normalize();
	
	// checking of frame errors :
	int CalcColumnStat( int col, int y_start, int y_end,
	               double& mean, double& sigma );
	BOOL_T CheckFrame( int& bad_line, BOOL_T bForce=FALSE );
	BOOL_T RepairFrame( int start_bad_line );
	void CalcColumnStat( int sigma_col, double& median, double& sigma34_14,
	                     double& mean, double& sigma );
	BOOL_T FindShift( int col , int step_size_adu, int sigma_adu, int& bad_line );        
	//-----------------------------------------------------------------------
	
	
	CCDMatrix& operator=(const CCDMatrix& right);
	CCDMatrix& operator+=(const CCDMatrix& right);

	// accessors :
	inline void SetParamSet( CCcdCfg* pParamsSet ){ m_pMatrixParamsSet = pParamsSet; }

	//
	void ClearState();
	inline void SetInteresting(BOOL_T bFlag=TRUE){ m_bIntersting=bFlag; }
	inline BOOL_T GetInteresting(){ return m_bIntersting; }
	inline long GetFoundEventCount(){ return m_FoundEventList.size(); }

	inline CCDEventList& GetGenEvents() { return m_GeneratedEventList; }
	inline CCDEventList& GetFoundEvents() { return m_FoundEventList; }
	int CountFoundEvents( int x, int y, int radius );
	// CCDEventList& GetIdentifiedEvents() { return m_IdentifiedEvents; }
	


	// calculation of sum :	
	
	// around given point :
	LONGLONG_T CalcSumAround( LONG_T x, LONG_T y , LONG_T r );
	
	// of group of pixels :
	LONGLONG_T CalcSum( CLongList& pixels );
	
	LONGLONG_T CalcSum( LONG_T* pixels, LONG_T cnt );


	// for creating images :
	void PutImage( long x0, long y0,
	               Table2D<ELEM_TYPE>& OtherImage,
	               long start_x=0,long start_y=0,long end_x=-1,long end_y=-1,
 	               BOOL_T bAdd=FALSE,BOOL_T bGen=FALSE,const char* mag=NULL,int frame_index=0,
 	               int* maxX=NULL, int* maxY=NULL );
	// BOOL_T GetImage( long x0, long y0,long len_x, long len_y,
	//                  Table2D<ELEM_TYPE>& Image );

	void AddImage( long x0, long y0,
						Table2D<ELEM_TYPE>& OtherImage,
						long start_x=0,long start_y=0,long end_x=-1,long end_y=-1,
						BOOL_T bGen=FALSE,const char* mag=NULL,int frame_index=0,
						int* maxX=NULL, int* maxY=NULL);

	CccdReport* PutSample( long x0, long y0,
	               Table2D<ELEM_SAMPLE_TYPE>& OtherImage,
	               long start_x=0,long start_y=0,long end_x=-1,long end_y=-1,
 	               BOOL_T bAdd=FALSE,BOOL_T bGen=FALSE,const char* mag=NULL,int frame_index=0,
 	               int* maxX=NULL, int* maxY=NULL );
	// BOOL_T GetImage( long x0, long y0,long len_x, long len_y,
	//                  Table2D<ELEM_TYPE>& Image );

	CccdReport* AddSample( long x0, long y0,
						Table2D<ELEM_SAMPLE_TYPE>& OtherImage,
						long start_x=0,long start_y=0,long end_x=-1,long end_y=-1,
						BOOL_T bGen=FALSE,const char* mag=NULL,int frame_index=0,
						int* maxX=NULL, int* maxY=NULL);

	BOOL_T IsGenerated( int x, int y );
	static void UpdateGenEvent( CccdReport& genEvent, CccdReport& foundEvent );		
	void UpdateGenEventsByArea( int low_x, int up_x, int low_y, int up_y );

	static void UpdateFoundEvent( CccdReport& foundEvent, CccdReport& genEvent );
	static void MoveGeneratedToBegin( CCDEventList& ccdEvents );
	
	static BOOL_T IsGenEventIdentified( const CccdReport& genEvent, const CccdReport& foundEvent );

	static BOOL_T IsGenEventIdentified( int gen_x, int gen_y, int found_x, int found_y );
	

	LONG_T CompileEventReport( CCDEventList& ccdEvents,
										LONG_T idx=-1, BOOL_T bAdd=FALSE, 
										BOOL_T bEraseFoundFromGen=TRUE,
										BOOL_T bProfile=TRUE );

	LONG_T CompileEventReportPtr( CCDEventList* ccdEvents,
										LONG_T idx=-1, BOOL_T bAdd=FALSE, 
										BOOL_T bEraseFoundFromGen=TRUE,
										BOOL_T bProfile=TRUE );
	
	static LONG_T CompileEventReport( CCDEventList& ccdEvents, 
										CCDEventList& foundEvents,
	                           CCDEventList& genEvents,
	                           CCDMatrix* pMatrix = NULL,		
										LONG_T idx=-1, BOOL_T bAdd=FALSE,
										BOOL_T bEraseFoundFromGen=TRUE,
										BOOL_T bProfile=TRUE );

	static LONG_T CompileEventReportPtr( CCDEventList* ccdEvents, 
										CCDEventList& foundEvents,
	                           CCDEventList& genEvents,
	                           CCDMatrix* pMatrix = NULL,		
										LONG_T idx=-1, BOOL_T bAdd=FALSE,
										BOOL_T bEraseFoundFromGen=TRUE,
										BOOL_T bProfile=TRUE );


	int RejectIfMoreThen( CCDPipeline* pPipeline );
	static int RejectIfMoreThen( CCDEventList& foundEventList, 
										  CCDPipeline* pPipeline,
										  BOOL_T bCheckConfirmedOnly=FALSE,
										  BOOL_T bAfterCoicCheck=FALSE );
	
	BOOL_T CheckIfOverlaps( LONG_T* cluster, LONG_T cluster_cnt,
	                        LONG_T& prev_evt_x, LONG_T& prev_evt_y );
	BOOL_T CheckIfOverlapsFast( long x, long y, 
	                            LONG_T* cluster, LONG_T cluster_cnt,
	                            LONG_T& prev_evt_x, LONG_T& prev_evt_y );
	BOOL_T CheckIfOverlapsFastMaxPoint( CCDEventList& found_events,
										 long x, long y, 
	                            LONG_T* cluster, LONG_T cluster_cnt,
	                            LONG_T& prev_evt_x, LONG_T& prev_evt_y );
	BOOL_T CheckIfOverlapsFast( CCDEventList& found_events,
	                               long x, long y,
	                               LONG_T* cluster, LONG_T cluster_cnt,
	                               LONG_T& prev_evt_x, LONG_T& prev_evt_y );
	                                                                                             	                            

	int CountOverlaps( long x, long y, 
	                      LONG_T* cluster, LONG_T cluster_cnt,
	                      LONG_T& prev_evt_x, LONG_T& prev_evt_y );


	void AddFoundEvent( CccdReport& tmp );
	void SetEventTime();

	void AddFoundEvent( long x, long y, CLongList& cluster );	

	CccdReport& AddFoundEvent( CPixelAnalyseIn& in, CPixelAnalyseOut& out );
	
	CccdReport& AddFoundEvent( CCDEventList& event_list,
	                           CPixelAnalyseIn& in, CPixelAnalyseOut& out,
	                          	BOOL_T bCheckUnique=FALSE );

	void CalcAdditionalInfo( ELEM_TYPE* p_data, long x, long y, long pos,
									 CPixelAnalyseOut& out,
                            CccdReport& event );
	               
	CccdReport* AddGeneratedEvent( long start_x, long start_y, long end_x, long end_y,
	                        const char* mag, int otherMaxX, int OtherMaxY, 
	                       	int frame_index,
	                        int& OutMaxX, int& OutMaxY,
	                        eEventType event_type=eFlash, int max_lap=0 );

	void AddGeneratedEvent( long x_s, long y_s, 
	                        const char* mag, int frame_index,
	                        eEventType event_type=eFlash );

	
	Table2D<BIG_ELEM_TYPE>* InitLaplaceFrame();
	void Laplace( CCDMatrix& out, long x0=-1, long y0=-1, long x1=-1, long y1=-1 );
	void Laplace( BOOL_T bEstimateOnEdge=TRUE );
	void LaplaceEdgeOnly( long ignore_edge );		
	void LaplaceEdgeOnly( long ignore_edge,
        		             CLongPoint* plus_list, LONG_T plus_count,
	                      CLongPoint* minus_list, LONG_T minus_count );
	void Laplace( long x0, long y0, long x1, long y1 );
	void LaplaceSafe( long x0, long y0, long x1, long y1 );
	void Laplace( eLaplaceType_T laplace_type, BOOL_T bEstimateOnEdge=TRUE );
	void Transform( eConfShape_T shapeType, double Redial, CCDMatrix& out );
	
	// matrix transformations :
	long CalcGaussFilter( Table2D<double>& gaussMatrix, CCDMatrix& out, 
								 long threshold, long x0=-1, long y0=-1, long x1=-1, long y1=-1 );

	long CalcGaussFilter( Table2D<double>& gaussMatrix, Table2D<double>& out, 
								 long threshold, long x0=-1, long y0=-1, long x1=-1, long y1=-1 );

	/*void UpdateGenEvents( long x, long y, LONGLONG_T newSum, LONGLONG_T otherSum,
           	             long x0, long y0, LONG_T* cluster, LONG_T cluster_cnt,
                         LONGLONG_T newClusterSum, LONGLONG_T* _prevClusterSum,
                         LONG_T prev_frames_cnt );*/
	void UpdateGenEvents( CPixelAnalyseIn& in, CPixelAnalyseOut& out );                         	

	// shifting/rotation/translation etc ...
	BOOL_T CalcAstroCoordinates( time_t ut_time, int x, int y, double& azim, double& altit, double& dec,
	                             double& ra, double& ha, CCDProcState& ccdInfo );

	BOOL_T CalcAstroCoordinates( int x, int y, double& azim, double& altit, double& dec,
										  double& ra, double& ha, CCDProcState& ccdInfo );

	BOOL_T CalcAzimutalCoord( int x, int y, double& azim, double& altit, CCDProcState& ccdInfo );


	// I/O
	static void GetEventDirAndList( CccdReport& event, mystring& szEvtDir, mystring& szEvtList );

	void WriteEventToFITSFile( CccdReport& event, int CurrentFrameIndex, 
										int EventNo, CCDPipeline* pPipeline,
										mystring& szEventDir, 
										const char* szMainSubDir="Events" );
	
	BOOL_T WriteToFITSFile( const char* fname, 
								 int low_x, int low_y, int up_x, int up_y,
								 CKeyTab* pHduTab=NULL, BOOL_T bNoCompress=FALSE );	
	
	BOOL_T WriteToFITSFile( const char* fname="$(DATADIR)/image.fit", BOOL_T bDumpEvents=FALSE,
	                      CKeyTab* pHduTab=NULL, BOOL_T bNeverCompress=FALSE );
	BOOL_T ReadFITSFile(  const char* fname="$(DATADIR)/image.fit", BOOL_T bMustBeAllocated=TRUE );
	BOOL_T ReadFITSHeader(  const char* fname );

	// static member functions:
	static void SaveEventFrameSum( CccdReport& event, const char* dir_path,
                           int startFrame, int CurrentFrameIndex,
                           int EventNo, CCDPipeline* pPipeline );
	static void GetFrameList( const char* dir_path, CccdReport& event,
										int EventNo,
										CMyStrTable& out_list, mystring& szOutFile,
										int first, int last );

	// add keys when saving part of image :
	static void AddPartKeys( CSafeKeyTab& eventHDUInfo,
									  int low_x, int low_y, int up_x, int up_y,
									  int evt_x, int evt_y,
									  int big_size_x, int big_size_y,
									  int evt_frame_index,
									  int curr_frame_index );	


	// getpicture params :
	BOOL_T GetPictureParams( const char* fname, int& sizeX, int& sizeY, int& x, int& y  );

	static BOOL_T ParseSaveArea( const char* szVal,
	                      int& x0, int& y0, int& x1, int& y1 );
	                                 	
	const char* getKeyValue( const char* keyname );
	const char* getKeyValue( const char* keyname, int start_pos );
	void setKeyValue( const char* keyname, const char* value );
	void setKeyValue( const char* keyname, double value );
	void setKeyValue( const char* keyname, int value );


	// reporting events :
	void DumpEventReport( const char* fname );	
	
	void DumpEventReport( const char* fname, CCDEventList& ccdReport );

	BIG_ELEM_TYPE** get_frame_laplace_fast();
	BIG_ELEM_TYPE* get_frame_laplace();
	
		
	double m_Average;
	cCCD* m_pFrame;


	Table2D<BIG_ELEM_TYPE>* m_pFrameLaplace;	
//	CCDMatrix* m_pFrameLaplace;

	// members :
	mystring m_szFileName;
	
	// parameters :
	CCcdCfg* m_pMatrixParamsSet;
	
	// from Skypic class by B&D :
	int createBackgroundLuminosityMap( int hibad_range, int x_start, int y_start,
												  int x_end, int y_end);

	// generation of images for simulator :
	int GenSkyImage( double sigma, double mean,
							double ra=0, double dec=0, 
							double fi=-1000, double pixscale=-1000,
							BOOL_T bCheckAlt=FALSE );

// opened shutter correction :
	void Data_Correction(double Tczyt=0.001, double Teksp=10.00);

};

#endif
