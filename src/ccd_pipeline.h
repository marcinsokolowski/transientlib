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
#ifndef _CCD_PIPELINE_H__
#define _CCD_PIPELINE_H__

#include "ccd.h"
#include <cfgtab.h>
#include <mystring.h>
#include <many_tab2D.h>
#include <pthread.h>
#include "ccd_eventlog.h"
#include "ccd_procstate.h"
#include "ccd_starspy.h"
#include "ccd_globals.h"
#include "ccd_state.h"
#include <mypipe.h>
// #include <myframeslist.h>
#include <asas_interface.h>
#include <ccd_common_struct.h>
#include <mystrtable.h>

#include <queue>
#include <vector>

using namespace std;

class AstroCCD;
class CDAQStatus;
class CListFile;
class CCD_Analyser;
class CRunStat;
class CccdReport;
class InfoTable2D;
class CCDPiSystemManagerInterface;
class CCDPiSystemManagerInterfaceCORBA;
class CCDRequestDef;
class CTriggerList;
class CTriggerAction;
class CImageCreator;
class CCDAsasTransform;
class CPimanInterface;
class CPiDBInterface;

enum eSavePictureType_T { eSaveNormalPicture=0, eSaveCalibPicture };

#define CCDPIPELINE_ACCEPT_NEW_AUTOPOST_DEFAULT FALSE

struct cDAQStateInfo
{
	cDAQStateInfo():m_bWaitingMode(FALSE),m_WorkingMode(eDAQNormalMode),
	m_MountRA_InRad(0),m_MountDec_InRad(0),m_MountAzim_InRad(0),
	m_MountAlt_InRad(0),m_MountObsMode(0),m_MountInfoTimeStamp(0),
	m_bDAQStarted(TRUE),m_eDomeStatus(eDomeOpened),m_LastDomeCheckTime(0),
	m_bScanMode(FALSE),m_LastCommand(eDAQReq_NotDef)
	{
		m_sObject[0] = '\0';
	}

	// sleep / wotking :
	BOOL_T m_bWaitingMode;
	
	// if DAQ started 
	BOOL_T m_bDAQStarted;
	
	// working mode :
	eDAQMode_T m_WorkingMode;
	BOOL_T     m_bScanMode;
	char       m_sObject[64]; // this is currently only for saving current OBJECT after trigger is obtained
	
	// dome status :
	eDomeStatus_T m_eDomeStatus;
	time_t        m_LastDomeCheckTime;

	// current mount info 
	double m_MountRA_InRad;
	double m_MountDec_InRad;
	double m_MountAzim_InRad;
	double m_MountAlt_InRad;
	int m_MountObsMode;	            
	time_t m_MountInfoTimeStamp;
		
	// global system error 
	mystring m_szSystemError;
	
	// last command :
	eCCDRequestType_T m_LastCommand;
};

class CCDPipeline 
{
public:
	void trace(const char* desc);

	// BEGIN OF FUNCTIONS DERIVED FROM CCDRequestExecutor
	// virtual functions derived from request executor class :
	virtual BOOL_T GetSynchroMode(){ return gCCDParams.m_bTakeInSynchroMode; }
	virtual void  GetSynchroMode( BOOL_T& bSynchroMode, BOOL_T& bAsynchroMode );
	virtual int GetXSize( int ccd_index );
	virtual int GetYSize( int ccd_index );
	virtual void GetPosition( double& ra, double& dec, double& azim,
	                          double& alt, int& obsMode );
	virtual int GetCCDCount() { return m_PipelineList.size(); }
	inline int GetPureGetNextTime() { return m_PureGetNextTime; }
	virtual void GetDAQParams( double& fov, double& alert_acceptance, double& min_alt_to_observe );
	virtual BOOL_T IsWaitingMode(){ return m_WorkingMode.m_bWaitingMode; }
	virtual BOOL_T IsNormalMode(){ return (m_WorkingMode.m_WorkingMode==eDAQNormalMode); }
	virtual BOOL_T GetWaitForFrame(){ return gCCDParams.m_bWaitForFrame; }
	virtual int GetWaitForAstroTime(){ return gCCDParams.m_nWaitForAstrometryInSec; }
	virtual BOOL_T DoAstrometryInTakeMode(){ return gCCDParams.m_bDoAstrometryInTakeNMode; }
	virtual BOOL_T IsDomeOpened();
	virtual int GetDarkCount();
	virtual int IsCooling();
	virtual void SetAstroFlagNow();
	virtual BOOL_T GetAstroOK();
	virtual time_t GetAstrometryTime();
	virtual BOOL_T IsAstroRunning();	   
	virtual time_t GetCurrentFrameTime();

	virtual int EmergencyExit( int bDoExit );										 
	virtual int IsScan();
	// END OF FUNCTIONS DERIVED FROM CCDRequestExecutor
	
	
	// device control and consistancy checks :
	static BOOL_T CheckDevices();
	static BOOL_T CheckConsistancy();
	

	// parameter set for all cameras :
	CcfgTable m_CamCfgTab;
	   
	// parameter set specific for pipeline object:
	CCDConfig m_PipelineCfg;
	
	// some custom header definitions :
	CMyMutex customKeyMutex;
	CSafeKeyTab m_CustomFITSKeys;
	void SetCustomKey( const char* keyname, const char* keyvalue );
	         	
	// eDAQNormalMode=0, eDAQSatTriggerMode, eDAQTriggerMovingMode	         		
	static BOOL_T gInStopAll;
	static BOOL_T m_bMountMove;
	static BOOL_T m_bStateRead;		
	static cDAQStateInfo m_WorkingMode;
	static cDAQStateInfo m_SaveStatusInfo;
	CCDState m_CamStateInfo;
	static int m_nAlertsCount;
	static int m_RestartDumpCounter;
	static mystring m_szPrevObjectValue;
	static int m_nStartAnalCounter; // counting StartAnalysis commands 
	
	static void SaveCurrStatus( CCDRequestDef* pReq=NULL, CCDPipeline* pPipeline=NULL ){}
	static void RestorePrevStatus();
	static void DumpCurrState();
	static void CheckIfExitOnError();
	static void DumpDAQStatus( BOOL_T bAll=TRUE, const char* msg=NULL, const char* mode=NULL, BOOL_T bForce=FALSE );
	static BOOL_T ReadCurrState();
	static void ReadCustomKeys();
	static BOOL_T IsCollectionMode() { return ::IsCollectionMode( m_WorkingMode.m_WorkingMode ); }

	// enabling/disabling fast mode :
	static void SetAsynchroModeOnOff( BOOL_T on_off );
	static void ShowAsynchroMode();
	
	// ASAS transformation :	
	CCDAsasTransform* m_pAsasTransform;
	CCDAsasTransform* m_pAsasTransformSav; // object to which transformation
														// from previous frame is saved 
//	CPipeLock runAstrometryLock; // lock for waking up astrometry thread 
	BOOL_T m_bPrevAstrometryOK;  // result of last performed astrometry
	
	CCDAsasTransform* m_pAsasTransformAsynchro;
	
	// dark frames list :
	CMyStrTable m_DarkFramesList;

	// ast files list :
	CMyStrTable m_AstFilesList;
	int m_nGoodAstrometryCount;
	
	// CFramesList m_FramesList;
	
	// astrometry :	
	int m_nAstrometryRetry;
	BOOL_T m_bForceAstroNow;
	BOOL_T m_bForceSynchroMode;
	BOOL_T m_bAstrometryRunning;
	BOOL_T m_bAstrometryExecuted;
	BOOL_T m_bAstrometryIsReady;
	int    m_nPhotometryStarsCount;
	double m_FwhmAver;
	static time_t m_NextAstrometryTime;
	mystring m_szLastAsynchroAstroFrame;
	
	// frames quality and errors :
	int m_LastBadFrame;
	int m_LastBadLineY;


	// finding image shift :
	vector<cStarCat> m_PrevStarList;
	vector<cStarCat> m_CurrStarList;			
/* ------------------------- PROTECTED ----------------------- */	
protected:

	static int m_PipelineCounter;
	static vector<CCDPipeline*> m_PipelineList;
	
	int m_PipelineIndex;
	mystring m_szCameraName;

	
	// all data specific to single camera will be kept in this structure :
	vector<CCDProcState> m_CameraStates;
	
	// star spy for each camera :
	CCDStarSpyTab m_StarSpyTab;
	
	// star spy - for transformation matrix :
	CCDStarSpyTab m_StarsForMatrix;
	
	// controll spy star - to control if calculated position 
	// agrees with real star position :
	CCDStarSpy* m_pControlSpyStar;
	

	const char* GetParam( const char* cfgname );
	

	long m_SizeX;
	long m_SizeY;
	long m_nCCD;	
	long m_nPipelineSize;
	long m_Count;
	int  m_FrameCounter;		
	int  m_DayFrameCounter; // set after read from device to day frame 
	int  m_FirstAnalysiedFrame;
	int  m_StartFrameIndex;
	int  m_SumedFrameIndex;
		
	// celestial coordinates controll :		
	BOOL_T m_bCoordInitialized;

	// restart pipeline on next frame 
	static BOOL_T m_bFinalDumped;
	
	// static 	
	BOOL_T m_bInitialized;
	static BOOL_T m_bPreActionsDoneOK;
	static BOOL_T m_bPostInitActionsDoneOK;

	
	// writing :
	mystring m_szOutputDir;

	void InitAnalyserObj();
	CRunStat* m_pAnalyser;
	
	// queue of frames coming from aparature
	cCCD* m_pCCD;

	// dark frame - commented out to save memory in tests 
	// stors values of uniform picture
	cCCD* m_pDarkFrame;
	
	// result of image subtraction 		
	Table2D<float>* m_pFlatFrame;
	
	// sum of several consecutive frames :
	Table2D<BIG_ELEM_TYPE>* m_pPrevSumOfSeveral;
	Table2D<BIG_ELEM_TYPE>* m_pCurrSumOfSeveral;
	Table2D<BIG_ELEM_TYPE>* m_pPrevSumOfSeveralLaplace;
	Table2D<BIG_ELEM_TYPE>* m_pCurrSumOfSeveralLaplace;				
	
	// old frame in memory :
	Table2D<BIG_ELEM_TYPE>* m_pOldSumOfSeveral;
	Table2D<BIG_ELEM_TYPE>* m_pOldSumOfSeveralLaplace;
	
	// maximal values from previous frames :
	cCCD* m_pPixelMaxFrame;
	
	// homeopatic sum ( or any other thing ) frame 
	CManyTab2D<BIG_ELEM_TYPE>* m_pHomeopaticFrame;

	CManyTab2D<BIG_ELEM_TYPE>* m_pLaplaceFrame;

	// working frames :
	CCDMatrix* m_pWrkMatrix;
	Table2D<BIG_ELEM_TYPE>* m_pWrkMatrixBigElem;
	Table2D<BIG_ELEM_TYPE>* m_pWrkMatrixBigElem2;
	void InitWorkingFrames();
	void CleanWorkingFrames();
	void LockWorkingFrames();
	void UnLockWorkingFrames();
	void LockWorkingTable( void* pFrame ){}
	void UnLockWorkingTable( void* pFrame ){}
	
	// background map :
	InfoTable2D* m_pMatrixInfoMap;

	// optymalization cache ONLY for MonteCarlo :
	cCCDCache* m_pFrameCache;	
	
	// Monte Carlo members :
	//  generator :
	CImageCreator* m_pGenObj;
	
	// for reading frames from log file :
	int m_LastFrameFromLogFile;
	CCDEventList m_EventsFromLogFile;
	time_t m_FrameUnixTime_FromLogFile;
	
	

	long tail; // newest frame 
	long head; // oldest frame 
	long iter;
	

	double m_CurrTime;
	double m_PrevTimeDiff;
	int m_PureGetNextTime;
	
	// lock for astrometry starting :
	BOOL_T m_bAstroThreadRunning;
	pthread_t m_AstroThreadID;
   int AstroThreadID;

	// functions :	
	void ShiftFrame( CManyTab2D<BIG_ELEM_TYPE>* pFrame );
	void DumpFrame( CManyTab2D<BIG_ELEM_TYPE>* pFrame );
	

	// Update homeopatic frame :
	void ShiftHomeopaticFrame();
	void UpdateHomeopaticFrame(){}
	void DumpHomeopaticFrame(int i);
	
	// laplace frame :
	void UpdateLaplaceFrame();		
	
	// Update maximum frame
	void UpdateMaximumsFrame(long iFrame=-1);

	// optimized version :
	void UpdateMaximumsFrameOpt(long iFrame=-1);		
	
	
	void ReCalcMaxFrame();
	
	// MC functions and members :
	CListFile* m_pList;
	LONG_T m_FrameIdx;
	mystring m_szLastFITSFile;
public :	
	// restart pipeline on next frame
   BOOL_T m_bRestartOnNext;

	// public frame info :
	time_t m_FrameUnixTime;

	// comparison to old frame :
	mystring m_szCurrSumOfNFileName;
	mystring m_szPrevSumOfNFileName;
	time_t   m_PrevSumOfFramesTime;    // time of first frame in the sum 
	int      m_nFramesCollectedForSum;
	int      m_nTotalCollectedSumFrames;

	// monte carlo :
	CCDEventList& GetEventFromLogFile(){ return m_EventsFromLogFile; }

	inline Table2D<BIG_ELEM_TYPE>* GetPrevSum(){ return m_pPrevSumOfSeveral; }
	inline Table2D<BIG_ELEM_TYPE>* GetCurrSum(){ return m_pCurrSumOfSeveral; }
	inline Table2D<BIG_ELEM_TYPE>* GetPrevSumLaplace(){ return m_pPrevSumOfSeveralLaplace; }
	inline Table2D<BIG_ELEM_TYPE>* GetCurrSumLaplace(){ return m_pCurrSumOfSeveralLaplace; }
	Table2D<BIG_ELEM_TYPE>* GetPrevSum( const char* szFileName );
	Table2D<BIG_ELEM_TYPE>* GetPrevLaplace( const char* szFileName );
	   
	cCCD* GetDark(){ return m_pDarkFrame; }
	Table2D<float>* GetFlat(){ return m_pFlatFrame; }
	void GetLastFrame( mystring& szFrame );
	virtual void GetLastFrameList( CMyStrTable& framesList );

	// PARAMTERS : 
	// bSaveFlip - if save flip accoring to global paramter
	// flip - overwrites global paramter , if >=0 bSaveFlip does not have to be TRUE 
	static void AddAstrometry( CSafeKeyTab& keyTab, CCDAsasTransform* pTransform, BOOL_T bSaveFlip=TRUE, int flip=-1 );
	
	void AddCelestialCoo( CCDMatrix& matrix, CCDAsasTransform* pAsasTransform=NULL );
	// void AddCelestialCoo( CSafeKeyTab& keyTab );
	void AddCelestialCoo( CSafeKeyTab& keyTab, double ra, double dec, 
								 double alt, double azim, double ut_start, 
								 double ha, eObservationMode_T obsMode,
								 CCDAsasTransform* pAsasTransform=NULL );
	
	inline InfoTable2D* GetBackgroundMap() { return m_pMatrixInfoMap; }
	
	inline eDAQMode_T GetWorkingMode(){ return m_WorkingMode.m_WorkingMode; }
	static void SetNewMode( CCDRequestDef& req );
	static void UpdateMode( BOOL_T bSetNormalMode=FALSE );
	static void InitMode();
	static void SetWaitingMode();
	static void SetWaitingModeOff( BOOL_T bWakeUpCollector=TRUE );	
	static BOOL_T CheckIfDarksExist();
	BOOL_T CheckIfDoFreshDarks();
	static void CheckIfDoFreshDarksBoth( BOOL_T bForce=FALSE );
	
	inline int GetXSize(){ return m_SizeX; }
	inline int GetYSize(){ return m_SizeY; }
	inline int GetCCDNo(){ return m_nCCD; }

	inline vector<CCDProcState>& GetCCDInfoTab(){ return m_CameraStates; }
	inline CCDProcState& GetCamStat(){ return m_CameraStates[0]; }
	AstroCCD* GetAstroCalcObj(){
		if( m_CameraStates.size()>0 ){
			return m_CameraStates[0].m_pAstroForms;
		}
		return NULL;
	}

	inline cCCDCache* GetCache() { return m_pFrameCache; }

	// searches for max star in current frame - in specific area (FIXED)
	void MarkTracedAllToReInit();
	void InitMAXStars( CCDStarSpyTab& starsSpyTab );
	BOOL_T CalcRotations();
	BOOL_T UpdateRotations();
	BOOL_T UpdateMAXStarsPositions( CCDStarSpyTab& starsSpyTab );
	void LogMAXStarsPositions( CCDStarSpyTab& starsSpyTab );
	void CheckIfReInitNeeded( CCDStarSpyTab& starsSpyTab );
	BOOL_T UpdateTransformMatrix( CCDStarSpyTab& starsSpyTab );
	void SetShiftParams( double angle, double x, double y, double dx, double dy,
								double dxpersec, double dypersec, BOOL_T bRot,
								BOOL_T bSetGlobal, int camIdx, double angle_per_sec );
	void PrintShiftParams();								
	void SetShiftParams();		
	void LogShiftChange();
	static void CheckIfShiftNeeded();
	
	// astro formulae usage :
	void GetCurrentFrameCoo( CCDMatrix& matrix,
	                         time_t& ut_time,
	                         double& ra, double& dec, double& alt, double& azim,
	                         double& ha );
	void GetCurrentFrameCoo( CCDMatrix& matrix,
												  time_t& ut_time,
												  double& ra, double& dec, double& alt, 
												  double& azim, double& ha,
												  eObservationMode_T& obsMode	 );
	                         
	void SetCurrentFrameTimes( CCDMatrix& matrix, int i );		
	void SetCurrentFrameTimes( time_t ut_time );		
	void SetCurrentFrameTimes();		
	
	
	CcfgTable& GetCamCfgTab(){ return m_CamCfgTab; }
	inline static vector<CCDPipeline*>& GetPipelineList() { return m_PipelineList; }
	static CCDPipeline* GetPipeline( int pipelineIndex );
	const char* GetCameraName();
	
	inline CCDConfig&  GetPipelineCfg() { return m_PipelineCfg; }

	void UpdateHomeoFrameBase(Table2D<BIG_ELEM_TYPE>& homeoMatrixIn, Table2D<BIG_ELEM_TYPE>& newMatrix, Table2D<BIG_ELEM_TYPE>& homeoMatrixOut ){}
	void UpdateHomeoFrameBase(Table2D<BIG_ELEM_TYPE>& homeoMatrixIn, CCDMatrix& newMatrix, Table2D<BIG_ELEM_TYPE>& homeoMatrixOut ){}
	void UpdateHomeoFrame(Table2D<BIG_ELEM_TYPE>& homeoMatrix, CCDMatrix& newMatrix, BOOL_T bUseMaxOfHomeo=TRUE ){}
	void UpdateHomeoFrame(Table2D<BIG_ELEM_TYPE>& homeoMatrix, Table2D<BIG_ELEM_TYPE>& newMatrix, BOOL_T bUseMaxOfHomeo=TRUE ){}

	// shifts :
	void InitShiftsInfo( InfoTable2D& shiftInfoTab );
	InfoTable2D* GetRotationMap();


	// astrometry functions :
	void UpdateAstrometry( mystring* szFileNames, BOOL_T bForceInAnyMode=FALSE );
	BOOL_T IsAstrometryOK();
	static void ResetAstro();
	BOOL_T SaveReducedFrame( mystring& szReductFile );
	BOOL_T ExecAstrometry( const char* szName  );
	BOOL_T ExecAstrometryAsynchro( const char* szName  );
	int    RunFastPhotometry( CCDMatrix& matrix, double tresh, const char* szMagFile, BOOL_T bSaveMag=TRUE );
	BOOL_T RunAsasAstrometry( const char* szReductFile, 
									  CCDMatrix* pMatrix=NULL,
									  CCDAsasTransform* pAsasTransform=NULL );
	BOOL_T RunPhotometry( const char* szFile, CCDMatrix& matrix, BOOL_T bSaveMag=TRUE );  
	BOOL_T SaveAstrometry( mystring* szFileNames );

	BOOL_T RemoveTemporaryFiles( const char* szMagFile,
	                             const char* szFlipedName,
	                             const char* szAstFile );
	                                                                                    									  
	void StartAstrometryThread();
	BOOL_T IsAsasTransformOK();
	BOOL_T UpdateCurrentFITSFile( const char* szFileName, CCDAsasTransform* pAsasTransform=NULL );
	BOOL_T UpdateCurrentFITSFile_PhotoOnly( const char* szFileName );

	// updating bacground info :	
	void UpdateBackgroundMap();
	void LogBackgroundMap();
	
	// suming frames :
	void UpdateSumOfFrames( mystring* szFileNames );

	double GetTreshold( int x, int y, int CameraIdx, eLaplaceType_T laplaceType, 
							double sig_above, 
							BOOL_T bAver=FALSE, int nAver=1 );

	void GetBackground( int CameraIdx, double& mS1, double& sS1, eLaplaceType_T laplaceType=eRawS );

	void GetBackground( int x, int y, int CameraIdx, double& mS1, double& sS1, eLaplaceType_T laplaceType=eRawS );
	
	static void GetBackgroundStat( CCDPipeline* pPipeline,
                           int CameraIdx, double& mS1, double& sS1,
                           eLaplaceType_T laplaceType=eRawS );

	static void GetBackgroundStat( CCDPipeline* pPipeline, int x, int y,
                           int CameraIdx, double& mS1, double& sS1,
                           eLaplaceType_T laplaceType=eRawS );
                                                                          	
	void SubtractBackground();

	// laplace
	void CalcLaplaceOfNew();


	// calc average of N previous frames :
	BOOL_T CalcAverageOfN( Table2D<BIG_ELEM_TYPE>& out, long prevN, long x0=-1, long y0=-1, long x1=-1, long y1=-1,
								  BOOL_T bLaplace=FALSE, CLongPoint* shapePoints=NULL, 
								  long pointCnt=0, BOOL_T bReCalcLaplace=FALSE);

	BOOL_T CalcAverageOfN( CCDMatrix& out, long prevN, long x0=-1, long y0=-1, long x1=-1, long y1=-1,
								  BOOL_T bLaplace=FALSE, CLongPoint* shapePoints=NULL, 
								  long pointCnt=0, BOOL_T bReCalcLaplace=FALSE );
								  
	long CalcAverageOfPrev( long prevN, long x, long y );

	BOOL_T CalcWeightedAverageOfN( CCDMatrix& out, long prevN, long x0=-1, long y0=-1, long x1=-1, long y1=-1,
			BOOL_T bLaplace=FALSE, CLongPoint* shapePoints=NULL, long pointCnt=0 );
								  

	CCDPipeline(int size=0,int nCCD=0,int xSize=0,int ySize=0);
	CCDPipeline( const CCDPipeline& right );
	~CCDPipeline();
	void Clean();

	CListFile* GetFramesList(){ return m_pList; }
	
	cCCD* GetVerifiedFrame();

	static BOOL_T m_bGlobalOptionsDone;
	static void InitSystemGlobalOptions( CCDPipeline* pPipeline ){}
	
	// initialization :
	CCcdCfg* GetCamParamSet( int CamIdx );	
	void ReadCamerasParams();
	void InitParams();			
	void RefreshParams();
	void SetParam( const char* cfgcode, const char* cfgval );
	void InitPipelineParams();
	void InitSizes( int xSize, int ySize, int nCCD, int size );
	void InitGenObj();
	
	// clean before new frame :
	void CleanBeforeNewFramePre();
	void CleanBeforeNewFramePost();
	
	// logs :
	void ClearConfirmedList();

	// some global actions to be executed at the begining :
	//  copy ccd.cfg file into output directory :
	static void ExecPreActions();
	static void ExecPostInitActions();
	static void DumpParams();
	
	// ACCESSORS :
	
	inline CImageCreator* GetGenObj(){ return m_pGenObj; }
	
	int InitDayFrameCounter();
	int GetDayFrameCounter();
	inline int GetFrameFromStart(){ return (m_FrameCounter-m_StartFrameIndex); }

	// same :
	inline int GetFrameIndex(){ return m_FrameCounter; }
	virtual int GetFrameCounter(){ return m_FrameCounter; }
	virtual int GetFrameCounter1();
	virtual int GetFrameCounter2();
	
	inline void SetDayFrameCounter( int DayFrameCounter ){ m_DayFrameCounter=DayFrameCounter; }
	inline int GetPipelineIndex() { return m_PipelineIndex; }
	static inline int GetPipelineCount() { return m_PipelineCounter; }

	CManyTab2D<BIG_ELEM_TYPE>* GetHomeopaticFrame(){ return m_pHomeopaticFrame; }	
	CManyTab2D<BIG_ELEM_TYPE>* GetLaplaceFrame(){ return m_pLaplaceFrame; }
			
	int SkipIfObjectMatches();
	long GetProcessedFrames(){ return m_FrameIdx;}
	void SkipFrameMC(){ m_FrameIdx++; }
	BOOL_T IsLastFrame();
	long GetCount();
	void Add(cCCD& cNew);
	long GetSizeX(){ return m_SizeX; }
	long GetSizeY(){ return m_SizeY; }
	long GetCCDcount(){ return m_nCCD; }
	long GetPipelineSize(){ return m_nPipelineSize; }
	cCCD* GetMaxFrame(){ return m_pPixelMaxFrame; }
	cCCD& GetNext();	
	cCCD& GetCurrent();
	cCCD* GetPrevFrame();
	cCCD* GetHead();
	cCCD* GetPrev();	
	cCCD* GetCurr();
	cCCD* GetFrame(LONG_T idx);
	cCCD* GetFrameObject(LONG_T idx);
	cCCD& AcceptNew(	BOOL_T bAutoPost=CCDPIPELINE_ACCEPT_NEW_AUTOPOST_DEFAULT,
							int bIncCounter=TRUE);	
	inline CRunStat& GetRunStatObj(){ return (*m_pAnalyser); }
	inline CCD_Analyser& GetAnalObj(){ return (CCD_Analyser&)(*m_pAnalyser); }		
	inline CCD_Analyser* GetAnalPtr(){ return ((CCD_Analyser*)m_pAnalyser); }
	void PostNewFrame();
	cCCD* FindFrame( LONG_T FrameIdx, LONG_T& pos );

	void GetAllMatrixPtrs( LONG_T index, CCDMatrix** pMatrixPtrTab );
	int GetAllMatrixPtrsChronological( LONG_T index, Table2D<ELEM_TYPE>** pMatrixPtrTab, BOOL_T bAddCurrent=TRUE );		

	int GetAllMatrixPtrsInt( LONG_T index, CCDMatrix** pMatrixPtrTab );
	int GetAllMatrixPtrsChronologicalInt( LONG_T index, CPixelAnalyseIn& in, 
													  BOOL_T bAddCurrent=TRUE, int frames_to_add=-1 );
	
	void GetNextMatrixPtrChronological( LONG_T index, LONG_T pos, CPixelAnalyseIn& in );

	int GetNPreviousFrames( LONG_T index, LONG_T nNextToTake, CCDMatrix** NextMatrixPtr );												   	

	// indexing frames :
	LONG_T GetCurrentIndex(){ return tail; }
 	
	LONG_T begin();
	LONG_T end();
	
	// clearing :
	
	// in order to begin MC analysis from the first frame again
	void ResetMCPipeline();	
	
	void ClearState();
	void ClearEvents( LONG_T FrameIndex );	
	void InitPipeline( BOOL_T bResetCounter=TRUE );
	void ResetPipelineFlags( BOOL_T bResetCounter=TRUE, BOOL_T bResetList=TRUE );
	void ResetPipelineLogs( BOOL_T bClearLog=TRUE);
	void ResetFinalLog( BOOL_T bCatToALL=TRUE );
	

	// use this function after changing position of mount 
	// so that shifts are re-determined :
	void RestartPipeline( BOOL_T bClearLog=TRUE );		
	
	void RestartPipelineKeepCurrent( BOOL_T bForce=FALSE );
	
	void InitCamStates();
	
	void InitAsasTransform();
	

	// analysing functions :
	// function prepares new frame to be analysed - subtraction of dark frame etc
	BOOL_T PrepareNewFrame();
	

	
	//-----------------------------------------------------------------
	//
	//        PI-Man communication section :
	//
	// Handling pre-analysis requests - check messages from 
	// System Manager Program and also re-calculate tresholds
	// if auto-calculation is reqired
	//-----------------------------------------------------------------
	static void AddNewAlert( CCDRequestDef& req );
	int ExecuteTriggerAction( CTriggerAction& action );
	int ExecuteTriggerActions();
	int ExecCheckAction( CTriggerAction& action );
	static void DumpToTriggerLog( CCDPipeline* pPipeline, CCDEventList& events );
	static void CheckExternalTriggers( CCDEventList& finallist1, CCDEventList& finallist2 );
	
	

	// call this function to update stored coordinates of frame ceneter
	// - AFTER MOUNT MOVEMENT !!!
	// also sets that now astrometry is required to recalculate transformation :
	void ChangeObsCoordinates( double ra, double dec, double azim, double alt, 
										eObservationMode_T obsMode, BOOL_T bAstroReq=TRUE );
	
	static void SetMountCoord( double ra, double dec, double azim, double alt,
 	       			            eObservationMode_T obsMode );
	
	// compares new coordinates to current - returns TRUE is different
	// and FALSE is same :
	BOOL_T CompareObsCoordinates( double ra, double dec, double azim, double alt, 
											eObservationMode_T obsMode );

	// call this function to update stored coordinates of frame ceneter
	// - AFTER SUCCESSFULL ASTROMETRY !!!
	void UpdateObsCoordinates( double ra, double dec, double azim, double alt );
	void UpdateObsCoordinates( double ra, double dec,
	                           double azim, double alt,
	                           eObservationMode_T& obsMode );
	                                                                                

	eObservationMode_T GetCamObsMode();
	static eObservationMode_T GetObsMode();
	eObservationMode_T m_PrevFrameObsMode;
	
	void LogError( const char* szErrorDesc );
	static void LogRequest( const CCDRequestDef& req, CCDPipeline* pPipeline );
	static void LogRequest( const char* szReqName, const char* szParamName,
	                        const char* szParamValue, CCDPipeline* pPipeline );
	static void LogRequest( const char* szReqName, const char* szParamName,
	                        int dParamValue, CCDPipeline* pPipeline );
	
	void UseAstrometryFromOther();
	static void ForceAstrometryOnNext();
	static BOOL_T ExecResetRequest();
	static BOOL_T ExecExitRequest( BOOL_T bExitFile=TRUE );
	static BOOL_T ExecLoadParamFile( CCDRequestDef& req );
	static BOOL_T ExecParamChange( CCDRequestDef& req );		
	static BOOL_T ExecTakeNPictures( CCDRequestDef& req );
	static BOOL_T ExecDoDarksRequest( CCDRequestDef& req, BOOL_T bWaitModeOff=TRUE );
	static BOOL_T ChangeParameterInAll( const char* szParamName, const char* szParamValue );
	static BOOL_T ChangeParameter( const char* szParamName, const char* szParamValue, int idx );
	static void CreateDaqExitFile();
	static void CleanAstList();
	
	
	BOOL_T WasAnyFrameTaken();
	static BOOL_T HandleRequests( BOOL_T bDoHandle ){ return TRUE; }
	eCCDRequestType_T HandleRequest( CCDRequestDef& req ){}
	static eCCDRequestType_T HandleSystemRequest( CCDRequestDef& req ){}
	static eCCDRequestType_T HandleDAQRequests( CCDPipeline* pPipeline ){}
	static eCCDRequestType_T HandleDAQRequests(){}
	static void WaitForReady( int sleep_time_msec=100);
	static void WaitForStart( int sleep_time_msec=500 );
	static void CheckAndHandleDomeStatus();
	static void WaitForCounter( int& counter, int final_value );
	
	// changing daq mode from SCAN -> NORMAL :
	static void EndScanMode();
	static void SetScanMode();

	// handling problems :
	static int EmergencyEnd( int bDoExit=0 );

	
	// every camera has PI-MAN communication object :
//	CCDController* m_pCCDDeviceInterface;
//	inline CCDController* GetCCDInterface() { return m_pCCDDeviceInterface; }
	
	// 
	mystring m_szCameraID;
	int GetCameraID();
	int GetNight();
	static int GetNightGlobal();
	int GetMinMaxAver( int& min_frame, int& max_frame );
	
	
	// system manager connection :
	
	
	
	void PreAnalysis();
	
	void AutoCalcTresholds();
	
	void LogSamplesInfo();		
	BOOL_T AnalyzeNewFrame(BOOL_T bAutoDump=FALSE,BOOL_T bReport=FALSE,LONG_T idx=0,BOOL_T bAutoSaveEvents=TRUE);
	BOOL_T VerifyOutFrame();	
	
	LONG_T VerifyPreviousFrame( LONG_T FrametoVerify, LONG_T nextMatrixCount );
	void SetNextFrameValues( LONG_T FrametoVerify, LONG_T nextMatrixCount=1 );
	void SetNextFrameValues();
			

	// functions for writing :
	void DumpPipeline(long num=1,long start=-1); 	

	//--------------------------------------------------------------------------------------------
	// functions for reporting :
	//--------------------------------------------------------------------------------------------
	static void GetOutFileName( mystring& szOutName, const char* szSubDir,
										 const char* szBaseName, CCDPipeline* pPipeline, 
										 int ccd_index=-1, BOOL_T bTXT=TRUE );
	void get_out_file_name( mystring& szOutName, const char* szSubDir,
	                        const char* szBaseName );
	                          
	void GetEventsLogFileName( mystring& szFile, const char* szBaseName, BOOL_T bFullPath=TRUE );
	
	CCDEventList& GetFoundEvents( int idx=0 );

	// currently found
	// CCDEventList& GetCurrentFoundEvents( int idx=0);

	CCDEventList& GetNewEvents( int idx=0 );
	
	// currently found 
	// CCDEventList& GetCurrentEvents();
	
	void SaveAllFrames();		
	void SaveCalibFrame();
	mystring SaveCurrentFrame( const char* szSubDirName="/Events/Triggers", eSavePictureType_T saveType=eSaveNormalPicture );
	void SaveEventsToFITS();		

	void SaveSumFrameEventsOnNext();
	void SaveSNSumFrameEventsOnNext();
	
	void SaveSumFrameEvents( const char* szFileName, CCDEventList& events );
	void SaveSNSumFrameEvents( const char* szFileName, CCDEventList& events );
	
	void SaveParts( Table2D<BIG_ELEM_TYPE>* pFrame, CCDEventList& evt_list,
			  	CCDMatrix& currFrame,
				int CamIdx, const char* szFileName, int saveFrameIndex,
				const char* szEvtSubDir="SumEvents" );
				
	int SaveParts( const char* list_file, CCDEventList& evt_list,
	               int CamIdx );
	                                                        				


	void WriteEventDetails( CCDEventList& eventList );
	int SaveEventsToFITS( CCDEventList& newEvents, BOOL_T bFirstDump, 
								 BOOL_T bAll=FALSE, BOOL_T bAutoDetermineIfFirst=FALSE,
								 BOOL_T bAutoCheckRange=FALSE, BOOL_T bSN=FALSE );		
	int SaveEventsToFITS( int lastEventFrameNo );
	
	// internal triggers identification and check :
	// on single events list :
	static int CheckForInternalTriggers( CCDEventList& verEvents, CCDEventList& triggerList );
	
	// on 2 coiciding sets of cameras :
	static int CheckCoicInternalTriggers( CCDPipeline& pipeline1, CCDPipeline& pipeline2 );
		

	void DumpAntyCoicEvents();
	void DumpNewEvents( LONG_T Idx, const char* msg=NULL, BOOL_T bToFile=TRUE );
	void DumpNewEvents( const char* szIndex, const char* msg=NULL, BOOL_T bToFile=TRUE );
	void DumpNewEvents( CCDEventLog* pEventLog, CCDEventList& newEvents, const char* szIndex, 
							  const char* msg, BOOL_T bToFile=TRUE );		
	void DumpNewEvents( CCDEventList* pInternalTriggerList=NULL , BOOL_T bForceAll=FALSE );
	
	void DumpSNCoic();
	static int DumpFinalSN( CCDPipeline* pPipeline1,
	       		            CCDPipeline* pPipeline2 );

	// dump final events and characteristic values :
	void DumpFinalEvents( CCDEventLog* pEventLog, CCDEventList& newEvents, 
								 BOOL_T bWriteToDB=FALSE );

	void AddFinalEvents( CCDEventList& confirmedEvents );
	
	void DumpCurrentEvents();		
	void UpdateSingleCamEventsList();
	
	void DumpCurrentCoicEvents();
	
	void RemoveOldRejectIfMoreTracks();
	
	// dumps events confirmed - index for frame to be printed
	// depends on number of next frames required for confirmation 
	void DumpConfirmedEvents( LONG_T FrameIndex );
	
	void InitLogfileName( CCDEventLog& eventLog, const char* szFileName );
	
	void ForceShiftsDetermination();
	void MoveCurrentShiftsFile();
	void GetShiftFileName( mystring& out );
		
	// 
	int m_ShiftFilesCount;			
	
	// members :
	CCDEventLog m_AllEventsLog;
	CCDEventLog m_VerifiedEventsLog;			
	CCDEventLog m_GenEventsLog;
	CCDEventLog m_ReGenEventsLog;		
	CCDEventLog m_FinalEventsLog;

	CCDEventLog m_SNEventsLog;
	CCDEventLog m_SNCoicEventsLog; // final
	CCDEventLog m_SNAllCoicLog;    // coic 

	CCDEventLog m_SNSumEventsLog;
	CCDEventLog m_SNSumCoicEventsLog; // final
	CCDEventLog m_SNSumAllCoicLog;    // coic
	   
	
	
	// re-generated events list - on current frame :
	CPointList m_ReGenEventsList;
	
	BOOL_T CheckIfReGenEvent( int x, int y );
	   
	
	//----------------------------------------------------------------------------------------------
	
	
	void ConfirmEventsFromPreviousFrames();
	
	void PrepareFoundEventsList();		
	int AddFoundEventsToList();
	void AddGenEventsToList();		
	void UpdateGenEventsList();
	
	
	// compiling events report generated vs identified :		
	int CompileEventsReport( LONG_T& genCount,LONG_T& foundCount,LONG_T& genIdent );
	static int CompileEventsReportFromFile( mystring& szFoundLog,
                                    mystring& szGenLog, mystring& szReGenLog,
												LONG_T& genCount,
											   LONG_T& foundCount, 
												LONG_T& genIdent, 
												CCDEventList& genEventsList );
	int CompileEventsReportFromFile( mystring& szFoundLog,
															 LONG_T& genCount,LONG_T& foundCount, 
														    LONG_T& genIdent, CCDEventList& genEventsList );
												
	int CompileEventsReportFromFileNew( mystring& szFoundLog,
	                                    LONG_T& genCount,
	                                    LONG_T& foundCount,
	                                    LONG_T& genIdent,
	                                    CCDEventList& genEventsList );
	                                                                                                                                                                                                                                      


	int CompileEventsReportFromFile( LONG_T& genCount,LONG_T& foundCount,LONG_T& genIdent, CCDEventList& genEventsList );
	
	
	// dumping all events :			
	int DumpFoundEventsLog( BOOL_T bForceAll=FALSE, CCDEventList* pInternalTriggerList=NULL );
	static int GetFinalCoicReport(  CCDEventList& eventCCD1, CCDEventList& eventCCD2,
												  CCDEventList& finallist1, CCDEventList& finallist2 );

	static void DumpSumEvents( CCDPipeline* pPipeline1, CCDPipeline* pPipeline2 );

	static void DumpEventsOnSumFrame( CCDPipeline* pPipeline1, CCDPipeline* pPipeline2,
					 CCDEventList& events_on1, CCDEventList& events_on2, const char* basename,
					 BOOL_T bSaveParts=TRUE, int limit=100, const char* szEvtSubDir="SumEvents" );

	static void DumpEventsOnOldFrame( CCDPipeline* pPipeline1, CCDPipeline* pPipeline2 );

	static void DumpEventsOnSumFrame( CCDPipeline* pPipeline1, CCDPipeline* pPipeline2 );

	static void DumpSNEventsOnSumFrame( CCDPipeline* pPipeline1, CCDPipeline* pPipeline2 );

	static int DumpFinalCoicReport();
	
	static int DumpFinalCoicReport( CCDPipeline* pPipeline1, CCDPipeline* pPipeline2 );

	static int DumpFinalCoicReport( CCDPipeline* pPipeline1, CCDPipeline* pPipeline2,
												  CCDEventList& finallist1, CCDEventList& finallist2,
												  const char* log1, const char* log2,
												  BOOL_T bAsFinalLog=TRUE );


	static int DumpFinalCoicReport( CCDPipeline* pPipeline1, CCDPipeline* pPipeline2,
									 CCDEventList& eventCCD1, CCDEventList& eventCCD2 );

	static int DumpFinalCoicReport( CCDPipeline* pPipeline1, CCDPipeline* pPipeline2,
												  CCDEventList& eventCCD1, CCDEventList& eventCCD2,
												  CCDEventList& finallist1, CCDEventList& finallist2 );
	static int InsertEventsToDB( CCDEventList& finallist1, CCDEventList& finallist2,
										  CCDPipeline* pPipeline1, CCDPipeline* pPipeline2,
										  BOOL_T bForce=FALSE );
	static int InsertEventsToDB( CCDEventList& finallist1, CCDPipeline* pPipeline1,
										  BOOL_T bForce=FALSE );
	
	int SetFinalRates( CCDEventList& finallist );
	
	void PrintEventDesc( CccdReport& eventReport, int Index,
	                     LONG_T xSize, LONG_T ySize, LONG_T MatrixIdx );
	
	void DumpAllEvents();
	
	// loging tracks :
	void DumpTracks( deque<CFrameEvents>::iterator& start,
	                 deque<CFrameEvents>::iterator& end );		
	                 
	BOOL_T CheckIfEventBelongsToTrack( double x, double y, double radius, CTrackDesc& track );
	BOOL_T CheckIfEventBelongsToTrack( double x, double y, int frame_index,
												  double radius, 
												  CTrackDesc& track, double chi2_limit, 
												  double& chi2, BOOL_T& bBelongs );
	static BOOL_T CheckIfEventBelongsToTrack( double x, double y, double radius, 
												  CTrackDesc& track, CTrackList& oldTracks );
	
	// function for testing purposes only :
	void InitRandom();			
	
	// Simutaltion function
	void Init_MC();		
	BOOL_T EndFileExists();
	BOOL_T ReReadListFile();
	void CheckNewFiles();
	void SkipFramesN( int toSkip );
	const char* GetNextFrameName();		
	int CheckNextDayFrameNo();
	const char* CheckNextFrameFileName();
	BOOL_T GetNextFrame_MC_AutoSum( mystring& szFileNames );
	BOOL_T GetNextFrame_MC(mystring& szFileNames, BOOL_T bUseCache=FALSE);
	BOOL_T GetNextFrame_FromLogFile(mystring& szFileNames, BOOL_T bUseCache=FALSE);
	static int ReadEvents_FromLogFile();
	BOOL_T GetASASAstrometryFromFITS();
	static BOOL_T ReadASASAstrometryFromFITS( CCDMatrix& currFrame, CCDAsasTransform* pAsasTransform, 
															eObservationMode_T& obsMode, int* nStars=NULL );
	static BOOL_T ReadASASAstrometryFromFITS( CSafeKeyTab& keyTab, 
						CCDAsasTransform* pAsasTransform,
						eObservationMode_T& obsMode );

	BOOL_T SetCoordFromFITS();
	
	BOOL_T CheckCoordChange();		
	BOOL_T CheckIfChangeCoordBig( double prev_ra, double prev_dec,											
	                              double prev_azim_in_deg,
	                              double prev_alt_in_deg,
	                              CCDAsasTransform* pAsasTransform );

	// part analysis :
	void UpdateEvents( CCD_Analyser* pAnal );

	// controll stars updates :
	void UpdateControlStar();
	void UpdateTraceOnAll();
	void UpdateAutoShiftStars();
	
	// calculating frame shift :
	BOOL_T CalcFrameShift( int nStarCount=20, double mag_min=-1000, double mag_max=1000, 
								  int min_x=-1, int min_y=-1,int max_x=10000,int max_y=10000, 
								  int min_dist=50, double delta_mag=0.5, BOOL_T bVerb=FALSE);
	
	// dark calculation :
	void UpdateDarkFrame();

	void NewFrameAnalyseEnd();


	void StopAnal();
		
	BOOL_T GetNextFrame_Real( mystring& szFileName, BOOL_T bClaimForNext=TRUE );
	BOOL_T GetNextFrame( mystring* szFileNames=NULL, BOOL_T bUseCache=FALSE, 
								BOOL_T bHandleReq=TRUE, BOOL_T bClaimForNext=TRUE,
								BOOL_T bOnlyCallGet=FALSE, BOOL_T bDoDumpAtLast=TRUE );
	BOOL_T GetNextFrameRaw( mystring* szFileNames, BOOL_T bUseCache,
									BOOL_T bClaimForNext=TRUE );
	BOOL_T GetNextFrameRawAcc( mystring* szFileNames, BOOL_T bUseCache );										
	
	
	void InitAndLockWrkEventList( int size );
	void UnlockWrkEventList(){}; 
	CCDEventList* m_pWrkEventList;

	CCDEventList m_SingleCamEvents;	// events from recent 20 frames are stored here 
	CCDEventList m_ConfirmedEvents;
	CCDEventList m_TriggerEventsList;
	
	// cosmic events list :
	CCDEventList m_AntyCoicEvents;
	
	// brithening events list 
	CCDEventList m_BrightenList;
	CCDEventList m_BrightenVerifList;

	// brithening events list 
	CCDEventList m_BrightenOnSumList;
	CCDEventList m_BrightenOnSumVerifList;

	CCDEventList m_EventsOnSumedFrame; // events on frames which are sums of previous N frames 
	CCDEventList m_OldEventsOnSumedFrame; // list of events to be dumped on future frames 
	CTrackList   m_TracksOnSumedFrames; // tracks found on sumed frames	

	CCDEventList m_EventsFromCompareToOld; // events from comparison to old frame
	
	
	// all found events : 	
	int m_LastWrittenFITSSize;
	int m_LastDumpedIndex;
	deque<CFrameEvents> m_allFoundEvents;
	
	// all generated events :
	deque<CFrameEvents> m_allGenEvents;
		
	// all generated nopt identified and found events :
	vector<CFrameEvents> m_CompiledEventsList;
	
	void ClearNewEvents();		

	CCDEventList m_EmptyEventsList;
	
	// events found in single frames analysis ( to be rejected in 
	// analysis on sumed frames ) :
	CCDEventList m_SingleFramesEvents;		
	int ReadSingleFrameEvents();
	
	CTrackList m_SingleFramesTracks;
	int ReadSingleFrameTracks();

	// found tracks :
	CTrackList m_TrackList;
	CTrackList m_CurrFrameTracks; // tracks from current frame only
	CTrackList m_PlaneTracks;     // tracks from planes, meteors etc - strict chi2 cut
										   // but no verlocity check
	CTrackList m_SingleCamTracks;										   
	
	inline CTrackList&  GetTrackList(){ return m_TrackList; }
	
	void End();
	static void StopAll();
};



class CCDPipelineIterator {
	LONG_T m_Index;
	cCCD* m_pFrame;
	CCDPipeline* m_pPipeline;		

public:
	CCDPipelineIterator(CCDPipeline* pPipeline);
	cCCD* operator->(){ return m_pFrame; }
	cCCD* operator++(int);
	cCCD* begin();
	cCCD* end();						
	cCCD* operator=( cCCD* pFrame){ 
		m_pFrame = pFrame; 
		return m_pFrame;
	};
	BOOL_T operator!=( cCCD* pFrame){
		return (m_pFrame!=pFrame);
	}
	BOOL_T not_end(){
		return (m_pFrame!=NULL);
	}
	BOOL_T last(){
		return (m_Index == m_pPipeline->GetCurrentIndex());
	}
	BOOL_T curr(){
		return (m_Index == m_pPipeline->GetCurrentIndex());
	}
};

// GLOBAL FUNCTIONS :
extern int GlobalCCDPipeline_EmergencyEnd( int bDoExit );

#endif
