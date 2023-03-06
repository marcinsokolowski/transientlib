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
#ifndef _CCD_CONTROLLER_H__
#define _CCD_CONTROLLER_H__


#include <mytypes.h>
#include <basedefines.h>
#include <paramtab.h>
#include <pthread.h>
#include <ccd_common_defines.h>
#include <mypipe.h>
#include <mylock.h>
#include <mysafekeytab.h>
#include "ccd_defines.h"

class CCDPipeline;
class DeviceBase;
class cCCD;
class CCDMatrix;
class CCDConfig;

class CCDController
{
protected:
	CPipeLock ccdReq;
	CPipeLock ccdAns;

public:
	// for ccddouble - to wait for both cameras to be fired in same moment :
	static int m_ToBeStarted;
	
	mystring m_szBaseFileName;

	CMyMutex ccdLock;		
	CMyMutex getDataLock;

	// 
	int m_nPipelineTotalSize;

	// no more frames in CCD :
	BOOL_T m_bNoMoreFrames;
	
	// 
	BOOL_T m_bCompleted;
	BOOL_T m_bContinue;
	BOOL_T m_bStarted;
	BOOL_T m_bPreGetCalled;
	eCCDTYPE_T m_CamType;
	BOOL_T m_bWaitMode;
	BOOL_T m_bReqForNextCalled; // for CCDPipeline to now if read last 
										 // frame before sleep 
										 
	BOOL_T m_bFirstInAsynchroMode; // should be set from outside to indicate first 
											 // frame in asynchro mode 
	
	// cooling status - TRUE if temperature is OK after waiting :
	BOOL_T m_bTempOK;											 
	
	// pipeline :
	CCDPipeline* m_pPipeline;
	
	cCCD* m_pCCD;

	// ccd-device
	DeviceBase*  m_pCCDDevice;
	
	// header of last collected frame is stored here 
	CSafeKeyTab m_Header;
	
	int m_NextFramePos;
	int m_CurrentFrameIndex;
	int m_FrameCounter;
	

	pthread_t m_ThreadID;
	int id;
	
	// 
	mystring m_szLastSavedFrame;
	
	// data buffer :
	double CCDTemp;
	double CaseTemp;
	double AmbientTemp;
public:
	// accessors :
	int InitDevice( BOOL_T bDoWaitForTemp=TRUE );
	int WaitForTemp( int max_wait_time=-1, int* waited_time_in_sec=NULL );
	DeviceBase* GetDriver(){ return m_pCCDDevice; }
	const char* GetDriverTypeName();
	BOOL_T GetDeviceName( mystring& szName, mystring& szID, int& idx );
	const char* GetName();
		
	inline BOOL_T IsStarted(){ return m_bStarted; }
	void SetCCDType( eCCDTYPE_T ccdType ){ m_CamType = ccdType; }
	
	// set force setting of shutter mode flag :
	void SetForceSetShutterMode( BOOL_T bFlag=TRUE );

	
	CCDController( CCDPipeline* pPipeline, eCCDTYPE_T camType, const char* szListFile="", const char* szBaseFileName="aaa_" );
	CCDController( eCCDTYPE_T camType, const char* szListFile="", const char* szBaseFileName="" );
	~CCDController(); 


	// runs colecting pictures in Second Thread :
	BOOL_T StartCollectingPictures( BOOL_T bDoInit=TRUE );
	void StopAnal();
	void WaitForStop();
	

	// runs collecting pictures in current thread :
	int Run( int nFrames=-1, time_t upto=-1 );


	// 
	void WriteToListFile( const char* fname, const char* fullpath, 
								 BOOL_T bDark, CCDMatrix& matrix, int frameno  );
		
	void HandleNewPicture( CCDMatrix& matrix, int FrameCounter, 
								  double mid_time, mystring* szName );

	ELEM_TYPE* GetNewFrame();

	ELEM_TYPE* GetNextFrame(int cam_idx,int frame_index, CSafeKeyTab& keyTab,
									mystring* szName=NULL, BOOL_T bClaimForNext=TRUE );
									
	ELEM_TYPE* GetNextFrameSynchro( CSafeKeyTab& keyTab,mystring* szName );
									
	int GetCurrentFrameIndex();
	
	BOOL_T IsWaitMode(){ return m_bWaitMode; }
	
	
	// waiting requests to driver :
	void SendRepeatFrameReq();
	void SendNextFrameReq();			
	char WaitForFrameRequest();
	

	// additional :
	void SendFirstFrameReq();
	
	// waiting for frames to be done :		
	void WaitForFrameToBeDone();
	void SendFrameDone();
	
	
	// ccd-operating :
	void UpdateDeviceParams();
	int SetDesiredTemp(double temp);
	int GetTempInfo( double& ccdtemp, double& settemp,
	                 double& casetemp, double& outtertemp );	                                
	int GetRealTemp( double& ccdtemp, double& casetemp, double& outtertemp );
	int SetCooling( int OnOff );		
	double GetElecGainValue();
	int GetCooling();
	void Close();
	void SetCamCfg( CCDConfig* pCfg );
	void SetAsynchroModeOnOff( BOOL_T on_off );
	

	// camera info buffer :
	void UpdateCameraInfo();

	// FITS :
	void GetFITSFileName( char* szFullPathFITS, int frameNo, 
								 BOOL_T bNoCompress=FALSE );
	BOOL_T InitDayFrameCounter();
	int GetDayFrameCounter();
	void SetDayFrameCounter( int nDayFrameCounter );

	// send request for other cameras to repeat theier pictures :
	int SendRepeatPicuturesRequest();
	
	BOOL_T RetryGetData( ELEM_TYPE* p_data_ptr , int size, int retry_count );


	// error loging :
	void LogError( eCamState_T cam_state );

	// simulator purposes :
	static void SetRADEC( double ra, double dec );
};


#endif
