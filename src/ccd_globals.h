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
#ifndef _CCD_GLOBALS_H__
#define _CCD_GLOBALS_H__
#include <ccdcfg.h>
#include <mypoints.h>
#include <mytypes.h>
#include <mystring.h>
#include <mymatrix.h>
#include <mylock.h>
#include <ccdrowcollist.h>
#include "ccd_defines.h"
#include <ccd_fits_header_defs.h>
#include <ccd_hardware_defines.h>

// #include <asas_gparam_def.h>

long GetNumberOfPoints( double brightness );

mystring GetEventTypeDesc( eEventType eventType );


class CCDMatrix;
class InfoTable2D;
class CCcdCfg;
class CCfgFile;
class CMatrix3X3;
class CTriggerActionList;
// struct GPARAM;
class CCDAsasTransform;
class CLog4C;

void RefreshStaticParamsInCcdlib();

#define _CCD_DEFAULT_TRACE_FILE_ "ccd.trace"

// default - IP of heplx42 :
#define CCD_DEFAULT_CAMERA_IP   "193.0.84.111"
// in simulator version :
// #define CCD_DEFAULT_CAMERA_PORT 23567
#define CCD_DEFAULT_CAMERA_PORT 1234 
#define CCD_DEFAULT_CAMERA_TIMEOUT 200
#define GLOBAL_PI_LOG "/opt/pi/dev/pisys/log/pi.log"


class CCDConfig
{
protected :
	// for saving current parameter state :
	vector<CEnvVar> m_SaveCfgTab;

	// configuration file parameters :
	BOOL_T m_bInitialized;
	long m_nInitCount;

	void InitBackgrFlagTable();
	void CheckIfBackgroundCalcEnabled();	   

	// if use gobal cfg file or local is provided :
	BOOL_T m_bUseLocalCfgFile;	
	CCcdCfg* m_pCfgFile;		
public:
	// SECTION OF STATIC PARAMETERS :
	static mystring gProgramName;
	
	static void DumpDAQStatus();
	
	CCDConfig( const char* szCfgFile=NULL );
	~CCDConfig();
	void InitDefaultValues();

	CCcdCfg& GetCfg();
	BOOL_T InitLocalCfgFile( const char* szCfgFile, CCDAsasTransform* pTransform=NULL );	
	inline BOOL_T IsInitialized(){ return m_bInitialized; }

	const char* GetParam( const char* name, BOOL_T bAllowNull=FALSE );
	void SetParam( const char* name, const char* value );
	void SetParam( const char* name, double value );
	CCfgFile& GetParamFile();
	CCfgFile* GetParamFilePtr();
	void GetParams( mystring& szParams );

	// indication of how program is used
	BOOL_T m_bCCDDouble;
	BOOL_T m_bParamsTest;

	// reporting formats :
	mystring m_FrameAnalLogFmt;
	mystring m_FrameAnalHeader;

	// tracing 
	int m_CCDTraceLevel;
	mystring m_szCCDTraceFile;
	mystring m_szErrorLogFile;

	// ini files:
	mystring m_szDefaultIniFile;
	
	// ERROR handling 
	int m_nExitOnToManyErrors;
	
	// enable / disable parameter changes logging :
	BOOL_T m_bSaveParamChanges;

	// Monte Carlo parameters :
	BOOL_T m_bMC;
	mystring m_szEndFilePath;
	BOOL_T m_bIgnoreDaqExitFile;
	BOOL_T m_bCheckGenOnly;
	BOOL_T m_bPutSample;
	LONG_T m_bGenEventRedial;
	int m_nPutObjOnNFrames;
	int m_nPutSampleEveryNFrame;
	BOOL_T m_bSumAllMethodsInNparTest;
	int m_nSamplesToPutOnFrame;
	BOOL_T m_bSaveImageWithSample;

	// simulator parameters :
	BOOL_T m_bUseRealFITSFiles;
	
	// samples :
	double m_fPutSampleCutOutInWidth;
	BOOL_T m_bPutSecondSampleByName;
	BOOL_T m_bLogSamplePut;
	int    m_nPutTakenFromFrameInRange;
	mystring m_szSampleXYFromListFile;
	int    m_nSamplePutRetry;
	
	// SuperNovea check parametrs :
	BOOL_T    m_bOnlySuperNewBackgr;
	mystring  m_szMinMagSUPERNEW;
	mystring  m_szMaxMagSUPERNEW;
	double    m_dRatioSUPERNEW;
	double m_fSuperNovaTvInSigma;
	double m_fSuperNovaTnInSigma;
	double  m_fSuperNovaMinPrevValue;
	int m_MaxFinalSN;
	   		

	// OUTPUT :			
	mystring m_szBaseOutDir;
	mystring m_szBaseOutSubDir;		
	mystring m_szOutDir;			

	BOOL_T m_bVisualizeEvents;
	LONG_T m_PartSize;
	LONG_T m_PartScale;

	// $CFGFILE parameters go here :
	
	// matrix description :
	LONG_T m_SizeX;
	LONG_T m_SizeY;
	LONG_T m_nCamNo;
	LONG_T m_nPipelineSize;
	int m_CameraIndex;// index of camera ! - for multi camera system
							// index of ccd_pipelineX.cfg 
		
	
	// if read should strictly check ealier declared size ( DEFAULT TRUE )
	BOOL_T m_bCheckSizeWhenRead;

	// general anaysis section :
	BOOL_T m_bMoveWeighted;
	BOOL_T m_bCorrectForRotation;
	BOOL_T m_bUseShiftTotal;		
	double m_FrameDX;
	double m_FrameDY;
	double m_FrameDXPerSec;
	double m_FrameDYPerSec;
	double m_RotCenterX;
	double m_RotCenterY;
	double m_RotValueDAlfa;
	double m_RotValueDAlfaPerSec;
	double m_MinShiftToCalcRot;
	double m_MinTotalShiftToCalcRot;
	BOOL_T m_bDoNotUseRotation;
	
	// frames shift :
	BOOL_T m_bUseFrameShift;
	double m_MinShiftToUse;
	double m_ForceAstroWhenBigShift;
	double m_ForceAstroWhenBigShiftRMS;

	BOOL_T m_bKeepLocalShift;		
	
	double m_S0;
	double m_S1;
	double m_S2;
	double m_S3;		
	CLongPoint m_overlapedPixels[4];
			
	
	BOOL_T m_bIsCOIC;
	BOOL_T m_bCheckForFlashes;
	BOOL_T m_bAnalyseMaxFromPrevFrames;
	BOOL_T m_bAnalyseSumAround;
	LONG_T m_FramesBack;
	LONG_T m_MaxAllowedVal;
	LONG_T m_MinClusterSize;
	LONG_T m_nPixelsAroundToConfirm;
	BOOL_T m_bUseClusterWithMore;		

	LONG_T m_nIgnoreEdge;
	LONG_T m_nIgnoreEdgeRight;		
	LONG_T m_nIgnoreEdgeLeft;
	LONG_T m_nIgnoreEdgeUp;		
	LONG_T m_nIgnoreEdgeBottom;		
	
	// average of prev N :
	BOOL_T m_bAverageOfPrevRejectMAXandMIN;		
	LONG_T m_nMaxOfAverageOfPrevN;
	BOOL_T m_bUseRotInAverageOfPrevN;
	BOOL_T m_bUseRotPerSec;		
	
	// analysis on sumed frames ( but sumed ealier by pi_red_frame ) :
	BOOL_T m_bOnSumedFrames;
	BOOL_T m_bUseFoundPosition;
	BOOL_T m_bAutoCalcSum; // automatically calculates sum of frames 
	BOOL_T m_bRejectSingleEvents;
	BOOL_T m_bAlwaysReReadSingleFrameEvents;
	BOOL_T m_bRejectSingleTracks;
	BOOL_T m_bAlwaysReReadSingleFrameTracks;
	
	// keep sum of several frames 
	BOOL_T m_bDoSumOfPrevNFrames;
	int m_bKeepSumOfPrevNFrames;
	BOOL_T m_bSaveSumOfNFrames;
	BOOL_T m_bAnalyzeSumOfPrevNFrames;
	BOOL_T m_bCheckNormalTracksOnSumEvents;
	double m_fVetoRadiusOnSumFrame;
	double m_nNewLaplaceInSigmaOnSum;
	double m_nMaxLaplaceOnOtherInSigmaOnSum;
	BOOL_T m_bRejectTracksOnSumedFrames;
	BOOL_T m_nNumBackFramesForTracksOnSum;
	BOOL_T m_bCheckPrevOfMaxPixel;
	
	
	// comparison to old frame, same settings as in m_bKeepSumOfPrevNFrames 
	int m_nCompareToOldFreqInSec;

	// FIRST LEVEL TRIGGER PARAMETERS :
	LONG_T m_TresholdPerPixel;
	LONG_T m_MaxPrevPerPixel;
	
	
	// MAX on all prev frames - possibly difficult to implement
	// due to frame rotations :
	BOOL_T m_bLocalMaxReq;
	BOOL_T m_bConfirmReq;
	BOOL_T m_bCalcClusterReq;
	LONG_T m_ConfRedial;
	LONG_T m_ConfTreshold;
	eConfShape_T m_ConfShape;
	LONG_T m_ConfTresholdPerPixel;
	LONG_T m_ConfMaxPrevPerPixel;
	// LONG_T m_MaxOnPrevAllowed;	


	// VERIFICATION IF NOT MORE THEN N EXCEEDS TRESHOLD :
	LONG_T m_nCheckIfNoneOfNExceedsTresh;
	double m_TreshForCheckIfNoneOfNExceedsInSigma;
	eConfShape_T m_eNotMoreThenNExceedsShape;
	LONG_T m_nNotMoreThenNExceedsRedial;		
	LONG_T m_nNotMoreThenNExceedsBackFramesCount;		
	LONG_T m_nNotMoreThenNExceedsPointsCount;
	CLongPoint mNotMoreThenNExceedsArea[MAX_CLUSTER_SIZE];

	// VETO AREA :
	eConfShape_T m_eVetoShape;
	double m_nVetoRedial;
	LONG_T m_nVetoPointsCount; // recalculated from m_eVetoShape, m_nVetoRedial
	CLongPoint m_VetoArea[MAX_CLUSTER_SIZE];		
	
	// MAX on previous frames (but only those kept in pipeline)
	BOOL_T m_bAnalyseMaxFromPrevInPipeline;
	
	// SUM
	eConfShape_T m_eNeigbShape;
	double m_dNeighbRedial;
	LONG_T m_nNeighbToSumCount;
	// LONG_T m_SumTreshold;		

	// homeopatic :
	BOOL_T m_bUseHomeoSameTreshForNewAndPrev;
	BOOL_T m_bKeepHomeopaticSum;
	LONG_T m_nHomeoAverageOfPrevNFrames;
	BOOL_T m_bCheckHomeoRawCond;

	// homeo alpha parameter :
	double m_HomeopaticFactor;	
	BOOL_T m_bCalcMaxNeighbHomeo;
	LONG_T m_nCalcMaxNieghbRedial;
	double m_nCalcMaxForAboveNSigma;
	BOOL_T m_bCalcMaxForAboveNSigmaOnHomeo;		
	BOOL_T m_bStartHomeoWithFirstFrame;		
	
	// laplacjan :
	BOOL_T m_bCalcLaplaceOfNEW;		
	BOOL_T m_bKeepLaplaceFrame;
	BOOL_T m_bCheckLaplaceCondition;
	BOOL_T m_bCheckLaplaceOnHomeoCondition;
	eLaplaceType_T m_eLaplaceType;		
	LONG_T m_nMaxLaplaceOnOther;
	LONG_T m_nMaxLaplaceOfPrevAverage;
	BOOL_T m_bMinLaplaceOnOtherEnabled;
	LONG_T m_nMinLaplaceOnOther;	
	LONG_T m_nNewLaplace;	
	double m_nNewLaplaceInSigma;
	double m_nMaxLaplaceOnOtherInSigma;
	BOOL_T m_bLaplaceUseSameTresh;
	BOOL_T m_bConfirmLaplaceMedianMinus;

	CLongPoint m_LaplacePlusList[MAX_CLUSTER_SIZE];
	LONG_T m_LaplacePlusCount;
	CLongPoint m_LaplaceMinusList[MAX_CLUSTER_SIZE];
	LONG_T m_LaplaceMinusCount;
	
	// variable objects and light increase :
	BOOL_T m_bCheckForSUPERNEW;
	LONG_T m_MinPrevPerPixel;
	LONG_T m_MinPrevTotal;		
	eBRIGHT_CHECK_TYPE m_eBrigtheningCheckType;
	LONG_T m_IncreaseTreshADUPerPixel;
	LONG_T m_IncreaseTreshADU;		
	double m_IncreaseTreshPercent;
			

	// noise :
	LONG_T m_MaxNoiseLevel;			
	double m_ClusterIfNSigmaAboveBackgr;		

	// conditions definitions :
	BOOL_T m_bCheckDifference;
	BOOL_T m_bCheckTreshAndMaxPrev;
	
	// confirmation conditions :
	BOOL_T m_bConfCheckDifference;
	BOOL_T m_bConfCheckTreshAndMaxPrev;
	
	// combined check params :
	BOOL_T m_bCombinedCheck;
	LONG_T m_UseDoubleCheckAboveADU;
	LONG_T m_UseDoubleCheckAboveADUPerPixel;
	
	// confirmation on next frames params :
	LONG_T m_ConfirmEventsOnNextNFrames;
	int    m_ConfirmOnNextRadius;
	BOOL_T m_bRejectNotConfirmed;		
	LONG_T m_MaxOnNextPerPixel;
	BOOL_T DoRejectNotConfirmed();

	BOOL_T m_bAlwaysUseSameTreshPerPixel;
	LONG_T m_SumTresholdForNewFrame;
	LONG_T m_SumOnPrevLessThen;
	LONG_T m_DiffTreshold;
	
	// difference check types :
	eDiffCheckType_T m_DiffCheckType;
	LONG_T m_DiffTreshSignal;

	// additional frames to be kept in memory :
	BOOL_T m_bDarkFrame;
	BOOL_T m_bFlatFrame;
	BOOL_T m_bOffsetFrame;
	BOOL_T m_bPixelMaxFrame;
	
	// corrections :
	BOOL_T m_bShutterCorr;
	
	// dark/dark frame file :
	mystring m_szDarkFrameFile;		
	mystring m_szFlatFrameFile;		

	// background :
	BOOL_T m_bSubtractBackground;
	eBackgrSubtrType_T m_eBackgrSubtrType;
	BOOL_T m_CalcBackgrTable[MAX_LAPLACE_DEFINED];
	BOOL_T m_bCalcBackgrOfCurrLaplace;

	mystring m_szRNoiseRowColDefFile;
	CWindowList m_RNoiseRopColDefList;
	
	// initial values :
	double m_AverageS1;
	double m_SigmaS1;
	
	// auto calculation of tresholds :
	BOOL_T m_bAutoCalcTresh;		
	BOOL_T m_bCalcTresholdsByBackgrMap;
	
	LONG_T m_BackgrMapXSize;
	LONG_T m_BackgrMapYSize;
	LONG_T m_BackgUpdateFreq;
	BOOL_T m_bBackgrDump;		
	

	// reports dumping :
	int m_nSaveCurrentPicture;	
	BOOL_T m_bSaveEventDescOnly;
	BOOL_T m_bSaveSuperNovaOnly;
	LONG_T m_nSaveFramesBeforeAndAfter;		
	BOOL_T m_bSaveFramesWithEvents;
	BOOL_T m_bSaveFramesWithEventsOnSum;
	LONG_T m_bSaveEventSize;
	BOOL_T m_bSaveOnlyGood;
	BOOL_T m_bSaveAverageParts;
	BOOL_T m_bDumpNewEvents;
	BOOL_T m_bDumpAllEvents;
	LONG_T m_DumpEventsFreq;		
	BOOL_T m_bDumpHomeoFrame;
	BOOL_T m_bStdoutOn;
	BOOL_T m_bGenNotFoundReport;
	int m_MaxStoredEventsFromSingleCam;

	// log files :	
	mystring m_szGenEventsLog;
	mystring m_szReGenEventsLog;
	mystring m_szVerifiedEventsLog;	
	mystring m_szRunEventsLog;
	BOOL_T m_bSameLogFiles; // if re-use log files 
	BOOL_T m_bDumpAllLogToStdout;

	// global pi.log using log4c library and syslog daemon 
	mystring m_szDAQLogFile;	
	BOOL_T   m_bPiLogEnabled;
	

	// performance optimalizations parametres:
	BOOL_T m_bKeepNeighbMap;

	// optimalizations :
	LONG_T m_EventsBufferSize;		


	// flying object cuts :-) :
	BOOL_T m_bSkipOverlaps;		
	LONG_T m_OverlapRedial;
	int m_MaxNumberOfEventsOnFrameAfterTv;
	int m_MaxNumberOfEventsOnFrame;		
	int m_MaxNumberOfEventsOnFrameAfterCoic;
	BOOL_T m_bRejectIfMoreVerb;		
	BOOL_T m_bUseOriginalXY;
	LONG_T m_bSkipIfMoreThen;
	LONG_T m_bSkipIfMoreThenRedial;
	LONG_T m_bSkipIfMoreThenMinDist;		
	eCheckIfMorePoint_T m_eCheckIfMorePoint;
	double m_CheckEventShape;
	BOOL_T m_bRejectBlackPixels;
	double m_fBlackPixelsRatio;
	double m_bBlackPixelsIfNSigmaBelow;
	BOOL_T m_bCheckCenterInCluster;
	int m_MaxPixelsInClusterAllowed;
	int m_MinStarsToAccEvents; // anti-cloud cut rejecting events on cloudy images
	
	// hot pixels rejection - based on list ;
	CPointList m_HotPixels;
	mystring m_szHotList;
	BOOL_T m_bRejectHotByList;
	BOOL_T m_bRejectHotPixelsByAverage;
	double m_nRejectHotPixelsTresholdInSigma;
	
	// list of bad areas of the chip :
	CWindowList m_BadPartsOfChip;

	// TRACKS :
	// verification if not track is found : TRACK :
	BOOL_T m_bCheckTracks;
	int    m_MaxEventRate;
	double m_MaxChi2InTrack;
	LONG_T m_nMinEventNoInTrack;
	LONG_T m_nNumBackFramesForTracks;
	LONG_T m_nNumBackFramesHasEventsForTracks;
	LONG_T m_nMaxEventsToFitTrack;	
	double m_nMinDistOfEventsInTrack;
	BOOL_T m_bCheckVelocity;
	double m_fVelocityError;
	// if more tracks :
	BOOL_T m_bCheckRejectIfMoreTracks;
	double m_MaxChi2ForRejectIfMoreTrack;
	int m_nKeepRejectIfMoreTracksOnN;
	
	//
	int m_nCheckFramesIfNotOlderThenNFrames;		
	
	// PLANE TRACKS - strict chi2 cut, but no verlocity cut :
	BOOL_T m_bCheckPlaneTracks;
	double m_MaxChi2InPlaneTrack;
	double  m_MaxChi2InPlaneTrackToOld;
	int    m_nNumBackFramesForPlaneTrack;
	int m_nMinPointNoOnPlaneTrack;
	int m_MinFramesOnPlaneTrack;
	
	// tracks on single frame :
	BOOL_T m_bFitLineToSingleFrameToAll;
	BOOL_T m_bFitLineToSingleFrame;
	BOOL_T m_fChi2ForLineOnSingleFrame;
	
	// tracks on single cam :
	BOOL_T m_bCheckTracksOnSingleCam;
	   
	

	// for adding new points (less strict criteria) :
	double m_MaxChi2ForPointToMatchLine;

	// checking for satelites :
	BOOL_T m_bCheckIfSatelite;
	BOOL_T m_bRejectSatelite;
	mystring m_szTleFile; // satelite database file 
	mystring m_szQthFile; // ground station location description
	double m_nSatRejRadius;
	double m_nNotVisibleSatRejRadius;
	
	// checking if flash is not just a star coming from outside the clouds :
	BOOL_T m_bCheckIfStar;
	BOOL_T m_bRejectStars;
	double m_fStarRejectRadius;
	double m_fCountBrighterThen;
	BOOL_T m_bCheckStarsInTycho; // if use tycho catalog to reject stars 
	double m_fStarCatMaxMag; // maximum star magnitude to be rejected by catalog
									 // in case of deep catalogs ( 15mag ) it would spoil flash recognition algorithm 
									 // if rejecting to faint objects !!!
	
	// big stars 
	BOOL_T m_bRejectIfBigStarNearBy;
	double m_fBigStarRejectRadiusInArcSec;
	double m_fBigStarMaxMagnitudo; // stars brighter then this value will be used

	// image errors :
	BOOL_T m_bCheckFrame;
	BOOL_T m_bRepairFrame;
	int m_ShiftCol;
	int m_ShiftVal;	
	BOOL_T m_bRejectEventsNearShift;
	
	
	// using DB :
	BOOL_T m_bUseDB;
	BOOL_T m_bUseODBC;
	mystring m_szDBName;
	mystring m_szDBUser;
	mystring m_szDBPass;
	mystring m_szDBHost;
	mystring m_szDB2Schema;
	mystring m_szScanDB;
	BOOL_T m_bSaveEventsToDB;
	BOOL_T m_bSaveFramesToDB;
	BOOL_T m_bSaveVerifEventsToDB;
	BOOL_T m_bSaveAllEventsToDB;
	BOOL_T m_bSaveTracksToDB;
	BOOL_T m_bSaveSumEventsToDB;
	int    m_eRunType;
	
	// queries for past GRB :
	BOOL_T m_bCheckForExternalsInDB;
	double m_fGrbRadiusInDeg;
	int    m_TimePastToCheckInSec;
	double m_SNRadiusInArcSec;
	

	// logs :
	BOOL_T m_bLogTracks;		
	BOOL_T m_bLogFrameStat;
	
	// image id :
	int m_DayFramesCounter;


	// checking if event a fluctuation on edge of big star :
	BOOL_T m_bCheckIfNotEdgeOfBigStar;
	BOOL_T m_bEdgeOfBigByRawData;
	double m_nSigmaAboveMeanInRawCluster;
	int m_nPrevFramesToCheckEdge;
	
	// OPT-3:
	BOOL_T m_bSuperOptimized;		
	BOOL_T m_bAverageOfPrev;
	BOOL_T m_bHomeopatic;
	
	// MonteCarlo optymalization :
	BOOL_T m_bKeepAllInMemory;
	BOOL_T m_bKeepSamplesInMemory;				

	// MC samples directories, filenames etc
	mystring m_szSampleDir;
	BOOL_T m_bPutScaledSamples;
	mystring m_szListName;
	
	// input :
	mystring m_szSampleFramesDir;
	mystring m_szFramesListFile;
	mystring m_szSkipIfObjectMatches;
	BOOL_T   m_bSkipBadFrames;   // MC only
	double   m_SkipWarmImages;
	BOOL_T   m_bCheckFrameOrder; // MC only
	BOOL_T   m_bCheckFramesListCount;

	
	// auto-update list 
	BOOL_T m_bAutoUpdateList;
	int m_nWaitTime;
	int m_nFrameListTimeout;
	
	// real analysis :
	// BOOL_T m_bReadFromDriver;
	BOOL_T m_bIgnoreCamera;
	
	BOOL_T m_bDoNotAnalyse;

	
	// two cameras in COICYDENCE :
	// transformation matrix :
	CMatrix3X3 m_TransformMatrix;
	mystring   m_szTransformFile;

	// instead of shifts transformation is found to convert frame X into frame Y :
	CTransformMatrixInTime m_FrameTransformMatrix;
	BOOL_T m_bTransformMatrixOK;
	int m_nMatrixNotFoundCount;

	// cosmic log - anty coic events :
	BOOL_T m_bLogAntyCoic;

	// coicydence redial :
	double m_nCoicRedial;
	BOOL_T m_bCoicByRaDec;
	double m_nCoicRadiusInRad;
	BOOL_T m_bPutSampleByRaDec;
	
	// matching stars to catalog :
	double m_fMatchStarToCatRadiusInArcSec;
	double m_fDiffMagToClaim;
	double m_fMinMagToClaimAlert;
	int    m_nNextFramesToConfirm;
	
	// cataloging 
	BOOL_T m_bUseAverageMagnitudo;
	double m_fFrameSize;

	// FITS - usful parameters :
	mystring m_DateObs;		
	mystring m_TimeObs;
	mystring m_DateTimeObsFormat;		

	// only for DEBUG version - all FALSE in release - by default !
	BOOL_T m_bDumpClusters;


	// Transforms :
	BOOL_T m_bIgnoreMissingTime;		
	BOOL_T m_bUseFrameTimeForShifts;
	BOOL_T m_bShiftUsesAstroFormulas;
	double m_TransformCCDOrientation;
	int m_CCDAxixXOrientation;
	int m_CCDAxixYOrientation;
	double m_TransformCCDFocus;
	double m_TransformCCDPixelSize;

	// Trigger actions :
	BOOL_T m_bSkipAstroOnAlert;
	int m_PostponeAstrometryTime;
	BOOL_T m_bHandleGCN;
	mystring m_szTriggerActionsDefFile;
	CTriggerActionList* m_DefaultTriggerActions;
	mystring m_szInternalTriggerActionsDefFile;
	CTriggerActionList* m_InternalTriggerActions;		
	int m_TriggerFollowTime;
	int m_InternalTriggerFollowTime;
	int m_TriggerCheckActionPeriod;
	BOOL_T m_bSaveFullFramesOnTrigger;
	int GetFollowTime( int triggerType );
	
	// significance of events and tresholds for internal triggers
	BOOL_T m_bSendInternalTriggers;
	
	// GEO - params
	double m_GeoLatitude;
	double m_GeoLongitude;
	double m_GeoAltitude;
	double m_SinOfGeoLatitude;
	int m_TimeZone;
	double m_HorAzimuth;
	double m_HorAltitude;	
	double m_HorAzimCorr;
	double m_HorAltCorr;
	double m_RACorr;
	double m_DecCorr;
	eObservationMode_T m_ObsMode;
	double m_DecObs;
	double m_RAObs;
	double m_HAObs;
	
	// asas photometry / astrometry :
	BOOL_T m_bUseFastPhotoInAstro;
	int    m_nMinStarCountToRunAstro;
	BOOL_T m_bUseGoodFromOther;		
	BOOL_T m_bDoASASPhotAstr;
	int    m_nSaveReductFramesFreq;
	BOOL_T m_bExecAstroOnFirst;
	BOOL_T m_bDoASASAstroOnStartup;
	mystring m_szLastReducedFrameName;
	double m_fASASPhotoThres;
	// ASAS astrometry :
	double m_fASASAstrometryFi;
	int m_nASASAstrometryOrd;
	mystring m_szASASStarCatalog;
	mystring m_szHIPStarCatalog;
	int m_bASASAstrometryVerb;
	int m_nASASAstrometryTry;
	int m_nASASAstrometryReTry;
	int m_nChangeToSilentAfterNGood;
	int m_nForceSynchroMaxCount;
	double m_fPixScale;
	BOOL_T m_bUseAsasTransform;
	eDriverReverseImage_T m_eReverseForTransform;
	int m_nAsasBorderSize;
	BOOL_T m_bAsasSubtrDark;
	BOOL_T m_bAsasDivideByFlat;
	double m_fAsasError;
	double m_fAsasFatalError;
	int m_nAutoAstrometryFreq;
	int m_nAutoAstrometryFreqInSec;
	BOOL_T m_bExecAstroSynchroMode;
	BOOL_T m_bUseAsasAstrometryFromFITS;
	BOOL_T m_bUseCoordFromFITS;
	BOOL_T m_bKeepMagAndAstFromAstro;
	BOOL_T m_bDoAstrometryInTakeNMode;
	BOOL_T m_bDoPhotoInTakeNMode;
	int    m_nWaitForAstrometryInSec;
	int    m_IgnoreCoordChangeOnCamera;
	BOOL_T m_bSaveGoodAstro;
	BOOL_T m_bSaveAstroInTakeNMode;
	double    m_BadAstroDX;
	double    m_BadAstroDY;

	// auto guiding :
	int m_CamIdxForAG_OnOff;
	
	BOOL_T m_bFixDeviceOrder;
	mystring m_szDeviceID;
	
	// fast-PI photometry
	mystring m_szFastPhotoCorrFile;
		
	
	// working modes :
	// BOOL_T m_bWaitingMode;
	BOOL_T m_bDoChangeShutterTime;
	int m_nDoFlatModeModulo;
	int m_nDoFlatModeResidue;
	int m_nDarksToBeTaken;
	int m_nDarksToBeTakenSav;
	BOOL_T m_bOverwriteOldDark;
	int m_nFramesToBeTaken;
	BOOL_T m_bReadAndSaveStatFile;
	BOOL_T m_bWaitForFrame;
	
	// dome status checking :
	int m_CheckDomeStatusFreqInSec;
	BOOL_T m_bStopOnDomeClose;
	
	//
	BOOL_T m_bDumpDAQStatus;
	mystring m_szDAQStatusFile;
	int m_nDumpStatusFreqInSec;


	// DRIVER PARAMETERS :
	mystring m_szKernelModulesPath;
	BOOL_T m_bReloadModuleOnStartup;
	BOOL_T m_bAsynchroMode;
	BOOL_T m_bParallelMode;
	mystring m_CCDIdentNo;		
	mystring m_CCDName;
	int m_SearchDeviceStartNo;
	eFITSCompressionType m_eCompressFITS;
	BOOL_T m_bDriverWriteAllToFITS;
	mystring m_szFramesOutputSubDir;
	BOOL_T m_bBuildFramesList;
	mystring m_szBaseFileNameFITS;

	// SHUTTER SETTINGS :
	double m_DriverShutterTimeInSec;
	double m_DriverShutterTimeInSecSaved;
	int m_ShutterMode;
	int m_ShutterModeOriginal;
	int  m_bHasBreakVoltages; // if camera has break voltages implemented ?
	char m_OpenBreakDelay;  // shutter breaking times see : http://grb.fuw.edu.pl/piwiki/MarcinSokolowski/soft/camera/driver
	char m_OpenBreakLength; // THE SAME	 
	char m_CloseBreakDelay; // THE SAME
	char m_CloseBreakLength;// THE SAME
	
	int m_DriverAnalBinning;
	eDriverReverseImage_T m_bDriverReverseImage;
	BOOL_T m_bUseCamID;
	int m_ReadFullStatusFreq;
	int m_nPipelineRestartTimeout;
	double m_CoordChangeToRestart;
	int m_MinStarCountToCloseShutter;
	
	// new options for gigabit ethernet camera :
	mystring m_szCameraIP;
	int      m_CameraPortNo;
	int      m_CameraTimeoutMiliSec;
	mystring m_szEthCamLogFile;
	mystring m_szEthCamErrFile;
	int      m_nEthCamLogLevel;
	int      m_nEthCamErrLogLevel;
	int      m_EthCamDataRetries;
	int      m_EthCamCmdTimeout;
	int      m_EthCamCmdRetries;
	int      m_LocalPortNo;
	
	// taking frames in synchro mode - no next frame claim :
	BOOL_T m_bTakeInSynchroMode;
	
	
	
	BOOL_T m_bDAQCommunicationON;
	int m_PortNo;
	mystring m_CorbaOptions;
	BOOL_T m_bCorbaWithNameService;
	int m_SharedMemKey;
	int m_PipelineSafeBufferSize;
 	int m_IntervalBetweenFrames;
 	int m_CameraInterfaceNo;

	// cooling, temperature :
	int m_bCoolingOnOff;
 	double m_CCDTemp;
 	BOOL_T m_bWaitForDesiredTemp;
 	int m_WaitForDesiredTempInSec;
 	eActionOnTempWaitFail m_eActionOnTempWaitFail; 	
 	int m_CoolingToleranceDiff;
 	
 	eCCDTYPE_T m_eCAMType;
 	eMPP_T m_eMPPMode;
 	int m_ReadoutSpeedHorizontal;
 	int m_ReadoutSpeedVertical;
 	unsigned char m_ADCConfReg0;
 	unsigned char m_ADCConfReg1;
 	unsigned char m_ADCConfReg2;
 	unsigned char m_ADCConfReg3;
 	unsigned char m_ADCConfReg4;
 	unsigned char m_ADCConfReg5;
 	unsigned char m_ADCConfReg6;
 	unsigned char m_ADCConfReg7;
 	int m_bForceRegistersUsage;
 	int m_ADCOffset;
 	int m_ADCOffsetCH2; 	
 	int m_ADCGain;
 	int m_ADCGainCH2;
 	BOOL_T m_ADCClampingBitOnOff;
 	eLNAGainValue_T m_eLNAGain;
 	eADCRange m_eADCRange; 	
 	int m_DriverMaxIterTimeout;
 	int m_RetryFrameCount;
 	int m_DriverGetDataRetryCount;
 	int m_MaxAllowedFramesLost;
 	BOOL_T m_bLensHitOnOff;
 	int m_MaxCommErrorCountToEXIT;
 	int m_SendBytesRetryCount;
 	int m_ExitOnCommError;

	// daq-simulator :
	BOOL_T m_bRepeatSameImages;

 	
	// COMMUNICATION WITH PIMAN :
 	BOOL_T m_bPISysManagerON;
 	int    m_PISysManPortNo;
 	BOOL_T m_bAskMountForCoord;
 	
 	// MOUNT INFORMATION :
 	int m_MountID;


	// AUTO Shifts calc :
	mystring m_szShiftsValuesFile;
	int m_nAutoShiftsMinStepsToIgnore;
	int m_nAutoShiftsCalc;
	int m_nAutoShiftMatrixNo;
	double m_MAXStarShapeLimit;
	double m_MAXStarNSigmaAbove;

	double m_MaxDiffSignalInSigma;
	int m_MinDistBetweenTracedStars;
	BOOL_T m_bTraceOnAllFrames;
	BOOL_T m_bUseTransformMatrix;
	BOOL_T m_bDumpTraceStarToFile;
	int    m_nMaxAllowedDrasticChng;
	BOOL_T m_UseNBackFrameForMatrix;
	int m_MaxFramesWithOldMatrix;
	int m_nSkipNFramesAfterChange;
	BOOL_T m_bUseControllStar;


	// AUTO calc size :
	BOOL_T m_bAutoCheckSize;
	mystring m_szMcInputFile;

	
	// SECOND LEVEL TRIGGER ( SLT )
	BOOL_T m_CheckClouds;
	int    m_TypicalStarsCount;
	double m_RejectFrameIfLess;
	double m_HoughDistrMaxLimit;
	double m_HoughTransformTresh;
	double m_HoughDistrTresh;
	double m_LapDiffMinRatio;
	BOOL_T m_bCheckHoughOnSmall;
	int m_nSmallHoughSize;
	double m_SmallHoughDistrTresh;
	int m_MinStarCountOnParts;
	BOOL_T m_bCheckCloudsOnPrev;
	int m_nMinBelowLimitToReject;
	BOOL_T m_bCheckTrackRADEC;
	
	BOOL_T m_bCheckHoughOnRaw;
	double m_HoughTransformOnRawTresh;
	double m_HoughDistrOnRawTresh;		

	// re-using old log files for tests of higher trigger levels :
	BOOL_T m_bReadFirstLevelInfoFromLog;
	BOOL_T m_bFromLogFile_WithFrames;
	mystring m_szFirstLevelTriggerDir;

	// auto turn off
	time_t m_RunUpTo;
	mystring m_szRunUpTo;

	// Observation / Instrument date  :
	mystring m_szOrigin;
	mystring m_szObserver;
	mystring m_szSite;
	mystring m_szInstrume;
	mystring m_szCamOptic;
	mystring m_szFilter;
	mystring m_szObject;
	double   m_FOV; // size on FOV in degrees 
	double   m_FOVBorder; // size of FOV border used in cataloging , 
								 // when selecting stars from DB 
	

	// Pointing paramteres :
	BOOL_T m_bPointIO;
	BOOL_T m_bAutoPointSWIFT;
	BOOL_T m_bAutoPointINTEGRAL;
	mystring m_szPointingPrior;
	mystring m_szSpecialTarget; // target of opportunity 
	
			
	// saving :
	void SaveParams();
	void RestoreParams();		
			
	void ModifyParams();
	BOOL_T RefreshParams();
	BOOL_T InitParams();
	BOOL_T InitParam(const char* szParamName);
	void InitDefaults();
	BOOL_T LoadFile( const char* szCfgFile );

	BOOL_T UpdateParamCommon( const char* szParamName, const char* szName );

	BOOL_T UpdateParam( const char* szParamName, const char* szName,
							  double& new_value, BOOL_T& bRetVal );
	BOOL_T UpdateParam( const char* szParamName, const char* szName,
							  int& new_value, BOOL_T& bRetVal );
	BOOL_T UpdateParam( const char* szParamName, const char* szName,
							  BOOL_T& new_value, BOOL_T& bRetVal );
	BOOL_T UpdateParam( const char* szParamName, const char* szName,
							  mystring& new_value, BOOL_T& bRetVal );

	BOOL_T GetASASTransform( void* TrParam );
	
	// END OF $CFGFILE aPARAMETERS SECTION 

	void SetIgnoreEdges( double left, double right, double bottom, double up );

	LONG_T GetEventsBufferSize();

	inline void SetMC(BOOL_T bMC=TRUE) { m_bMC = bMC; }
	inline BOOL_T GetMC(){ return m_bMC; }
	
//	void SetFromSamples(BOOL_T bFromSamples=TRUE){ m_bPutFromSamples=bFromSamples; }
//	BOOL_T GetFromSamples(){ return m_bPutFromSamples; }	

	void SetOutputDir( const char* szOutDir );
	const char* GetOutputDir();

	BOOL_T GetVisualizeEvents();
	void SetVisualizeEvents( BOOL_T bVisualizeEvents ){ m_bVisualizeEvents=bVisualizeEvents;}

	BOOL_T GetKeepMapFlag();
	BOOL_T GetCalcBackgrFlag();

	void RecalculateParams( CCDMatrix* pMatrix=NULL, InfoTable2D* pMatrixInfoTable=NULL );
	void ReCalcTresholds( CCDMatrix* pMatrix=NULL, InfoTable2D* pMatrixInfoTable=NULL );

	// parameters dependencies and consistancy checks go here :
	BOOL_T CheckConsistancy();

	void CheckEdges();

	void CalcShapes();
	
	void SetShiftsToZero();

	void Dump();

	long GetEdgeSizeInLaplace( eLaplaceType_T type );
	long GetEdgeSizeInLaplace();

	double GetPixScale();
	double GetPixToRad( double nPixels );
	
	int GetTrLevel();
	mystring GetTrFile();

	LONG_T GetAnalNeighbCount(){ 
		InitParams();
		return m_nNeighbToSumCount; 
	}
	
	// qth location file for satelites :
	void InitSatLib();
	void InitQthFile( const char* qthfile );

	// pi log
	BOOL_T InitPiLog( BOOL_T bDoInit );

	// returns OK if value>0 and FAILED when value<=0
	static const char* GetResultFlag( int value );
	
};

extern CCDConfig gCCDParams;
extern time_t gStartUTTime;
extern BOOL_T gCameraError;
extern mystring gRunDate;
extern CLog4C gPiLog;

// lock to protect global "C" functions in asaslib :
extern CMyMutex gAsasAstroLock;

#endif
