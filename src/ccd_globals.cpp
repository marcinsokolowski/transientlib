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
#include "ccd_globals.h"
#include "ccd_analyse.h"
#include "ccd_pipeline.h"
#include "ccd_datares.h"
#include "ccd_photometry.h"
#include "cfg.h"
#include "ccd_trace.h"
#include "ccd_starcat.h"
#include <cexcp.h>
#include <mycmnglobals.h>
#include <mymacros.h>
#include <myhisto.h>
// #include <unistd.h>
#include <ccd_fits_header_defs.h>
#include <ccd_piman_interf_defs.h>
#include <ccd_hardware_defines.h>
#include <math.h>
#include <myutil.h>
#include <mydate.h>
#include <myparser.h>
#include <sat_info.h>

#include <AstroHeader.h>
#include <asas_gparam_def.h>
#include "ccd_asastransform.h"
#include <AstroCCD.h>
#include <myfits.h>
#include <mathfunc.h>


extern "C" {
#include <asas_astrometry.h>
}


// global params object :
CCDConfig gCCDParams;

// time set when program is loaded - to know global start date :
time_t gStartUTTime=get_dttm();
mystring gRunDate=get_night_date_local();

// flag error in camera :
BOOL_T gCameraError=FALSE;

mystring CCDConfig::gProgramName;

// this to compile in values of shared memory keys :
// #include <ccddriver_interface.cpp>

long GetNumberOfPoints( double brightness )
{
	if(brightness>=1.00 && brightness<=10.0)
		return 100;
	if(brightness<1.00 && brightness>0.2)
		return 10000;
	if(brightness<=0.2 && brightness>0.1)
		return 100000;
	if(brightness<=0.1)
		return 1000000;
	return 200;
}

void RefreshStaticParamsInCcdlib()
{
	// CCDPipeline::RefreshParams();
	CCD_Analyser::RefreshParams();
	gCCDParams.RefreshParams();
}



BOOL_T CCDConfig::LoadFile( const char* szCfgFile )
{
	BOOL_T bOverSav = FALSE;
	BOOL_T bLocalSav = FALSE;
	if( m_pCfgFile ){
		bOverSav = m_pCfgFile->m_bOverwriteMode;
		bLocalSav = m_bUseLocalCfgFile;			
	}
	
	CCcdCfg* pSavCfg = m_pCfgFile;
	m_pCfgFile = new CCcdCfg( szCfgFile );
	m_pCfgFile->m_bOverwriteMode = TRUE;
	m_bUseLocalCfgFile = TRUE;
	BOOL_T bRet = RefreshParams();
	m_pCfgFile->m_bOverwriteMode = bOverSav;
	m_bUseLocalCfgFile = bLocalSav;
	delete m_pCfgFile;
	m_pCfgFile = pSavCfg;

	return bRet;
}

BOOL_T CCDConfig::InitLocalCfgFile( const char* szCfgFile, CCDAsasTransform* pTransform )
{
	if( m_pCfgFile )
		delete m_pCfgFile;

	// NEW - 20041020 , change due to viewlist! - when 
	// set CFGFILE=$SRCDIR/ccd/cfg/WZORY_CCD_CFG/ccd_to_view_events.cfg
	// then due to m_bOverwriteMode=TRUE in line : cfg.cpp:329
	// program does not try to use it - so here I try reading same 
	// file as Global cfg was set before - instead of auto-guasing it again :
	CCfgFile* pGlobalParamFile = GetGlobalParamFilePtr();
	mystring szGlobalCfgFileName = DEFAULT_CFG_FILE;
	if( pGlobalParamFile ){
		const char* szName = pGlobalParamFile->GetFileName();
		if( szName && szName[0] ){
			szGlobalCfgFileName = pGlobalParamFile->GetFileName();
		}
	}

	// first read default params from ccd.cfg :
	m_bUseLocalCfgFile = TRUE;
//	m_pCfgFile = new CCcdCfg( szGlobalCfgFileName.c_str() );
	if( gCCDCfg ){
		_TRACE_PRINTF_2("Initializing pipeline config by copy constructor\n");
		m_pCfgFile = new CCcdCfg( (*gCCDCfg) );
	}else{
		_TRACE_PRINTF_2("Initializing pipeline config from file : %s\n",szGlobalCfgFileName.c_str() );
		m_pCfgFile = new CCcdCfg( szGlobalCfgFileName.c_str() );		
	}

	m_pCfgFile->m_bOverwriteMode = TRUE;
	RefreshParams();
	if( pTransform ){
		pTransform->m_bTransformOK = GetASASTransform( pTransform );		
		if( pTransform->m_bTransformOK ){
			pTransform->m_bReadDone = TRUE;
		}
	}	
	delete m_pCfgFile;

	m_pCfgFile = new CCcdCfg( szCfgFile );
	m_pCfgFile->m_bOverwriteMode = TRUE;


	return RefreshParams();
}

const char* CCDConfig::GetParam( const char* name, BOOL_T bAllowNull )
{
	// printf("CCDConfig::GetParam start\n");
	if(!m_bUseLocalCfgFile){
		return GetGlobalParam( name, bAllowNull );
	}else{
		// printf("here 2 ?\n");
		Assert(m_pCfgFile!=NULL,"Local cfg file not initialized");
		return m_pCfgFile->GetParam( name, bAllowNull );
	}
}

void CCDConfig::SetParam( const char* name, const char* value )
{	
	SetGlobalParam( name, value );	
	if(m_bUseLocalCfgFile){
		Assert(m_pCfgFile!=NULL,"Local cfg file not initialized");
		return m_pCfgFile->SetParam( name, value );
	}
}

void CCDConfig::SetParam( const char* name, double value )
{
	char szValue[100];
	sprintf(szValue,"%.10f",value);
	SetParam( name, szValue );
}

CCfgFile& CCDConfig::GetParamFile()
{ 
	if(m_bUseLocalCfgFile){
		Assert(m_pCfgFile!=NULL,"Local cfg file not initialized");
		return m_pCfgFile->GetParamFile();
	}else{
		return GetGlobalParamFile(); 
	}
}

CCfgFile* CCDConfig::GetParamFilePtr()
{
	CCfgFile* ret=NULL;
	if(m_bUseLocalCfgFile){
		// printf("here 1 ??\n");
      Assert(m_pCfgFile!=NULL,"Local cfg file not initialized");
      ret = m_pCfgFile->GetParamFilePtr();
   }else{
		// printf("here 2 ??\n");
      ret = GetGlobalParamFilePtr();
   }
	return ret;
}

void CCDConfig::GetParams( mystring& szParams )
{
	if(m_bUseLocalCfgFile){
		Assert(m_pCfgFile!=NULL,"Local cfg file not initialized");
		m_pCfgFile->GetParams( szParams );	
	}else{
		GetGlobalParams( szParams );
	}
}

CCcdCfg& CCDConfig::GetCfg(){
	Assert(m_pCfgFile!=NULL,"Local cfg file not initialized");
   return (*m_pCfgFile);
}


CCDConfig::CCDConfig( const char* szCfgFile )
: m_bUseLocalCfgFile(FALSE), m_pCfgFile(NULL)
{

	InitDefaultValues();
	if( szCfgFile && szCfgFile[0] ){
		InitLocalCfgFile( szCfgFile );
	}
}

CCDConfig::~CCDConfig()
{
	if(m_pCfgFile)
		delete m_pCfgFile;
}

void CCDConfig::InitDefaultValues()
{
	m_bCCDDouble = FALSE;
	m_bParamsTest = FALSE;

	// enable / disable parameter changes logging :
   m_bSaveParamChanges = TRUE;

	mystring szNight=get_night_date_local();
	m_szErrorLogFile = "";
	m_szErrorLogFile << "/opt/pi/dev/pisys/log/daq_" << szNight << ".err";

	m_bSkipAstroOnAlert = FALSE;
	m_bHandleGCN = FALSE;
	m_PostponeAstrometryTime = 300;
	m_szTriggerActionsDefFile="$(SRCDIR)/ccd/cfg/WZORY_CCD_CFG/trigger_actions.cfg";
	m_szInternalTriggerActionsDefFile="$(SRCDIR)/ccd/cfg/WZORY_CCD_CFG/trigger_actions_internal.cfg";
	m_TriggerFollowTime=0;
	m_InternalTriggerFollowTime=0;
	m_TriggerCheckActionPeriod = 12*3600; // in seconds (default 12h )
	m_bSaveFullFramesOnTrigger = FALSE;

// event significance :
	m_bSendInternalTriggers=FALSE;

	m_CCDTraceLevel=0;
	m_szCCDTraceFile=_CCD_DEFAULT_TRACE_FILE_;

	m_FrameAnalLogFmt="%d %d %d %d %d %d %d %d %d %d %d %d\n";
	m_FrameAnalHeader="Frame# Tn Tv MeanG SigmaG MeanS SigmaS Event# T_hot T_cluster FitG FitS\n";

	m_nInitCount=0;

	// ini file for ccdview
	m_szDefaultIniFile = DEFAULT_INI_FILE;
	m_szDefaultIniFile.env2str();


// defaults should be zero - in case it is zero programs can read CCDMatrix 
// and reading procedure will set size as defined in FITS - otherwise 
// if deifferent sizes then default here - assertion fails :
	m_SizeX=0;
	m_SizeY=0;

	m_nCamNo=1;
	m_nPipelineSize=6;
	m_CameraIndex = -1;

	m_bCheckSizeWhenRead = TRUE;

	m_szEndFilePath=DAQ_EXIT_FILE_STAMP;
	m_bIgnoreDaqExitFile=FALSE;
	m_bCheckGenOnly=FALSE;
	m_bMC=FALSE;
	m_bGenEventRedial=1;
	m_bPutSample=FALSE;
	m_bInitialized=FALSE;
	m_nPutObjOnNFrames=0;
	m_nPutSampleEveryNFrame=-1;
	m_bSumAllMethodsInNparTest=FALSE;
	m_nSamplesToPutOnFrame = 1;
	m_bSaveImageWithSample = FALSE;

	// simulator :
	m_bUseRealFITSFiles = TRUE;

	// analysis on sumed frames :
	m_bOnSumedFrames = FALSE;
	m_bUseFoundPosition = FALSE;
	m_bAutoCalcSum = FALSE;
	m_bRejectSingleEvents = TRUE;
	m_bAlwaysReReadSingleFrameEvents = FALSE;
	m_bRejectSingleTracks = TRUE;
	m_bAlwaysReReadSingleFrameTracks = FALSE;

	m_fPutSampleCutOutInWidth=-1;
   m_bPutSecondSampleByName=FALSE;
	m_bLogSamplePut=FALSE;
	m_nPutTakenFromFrameInRange=-1;
	m_nSamplePutRetry = 5;


	m_szMinMagSUPERNEW="8.5"; 
	m_szMaxMagSUPERNEW="9.0";
	m_dRatioSUPERNEW=1.5;
	m_bOnlySuperNewBackgr=FALSE;
	m_fSuperNovaTvInSigma = 4.0;
	m_fSuperNovaTnInSigma = 8.00;
	m_fSuperNovaMinPrevValue = 2.5;
	m_MaxFinalSN = 10;


// OUTPUT :
	m_szBaseOutDir="RESULTS";
	m_szBaseOutSubDir="";
	m_szOutDir="";


	// run up to:
	m_RunUpTo = 0;

// visualization section :
	m_bVisualizeEvents=FALSE;
	m_PartSize=40;
	m_PartScale=4;

// FIRST LEVEL TRIGGER PARAMS :
	m_TresholdPerPixel=400;
	m_MaxPrevPerPixel=40;


// averarage of previous frames :
	m_bAverageOfPrevRejectMAXandMIN=FALSE;
	m_nMaxOfAverageOfPrevN=0;
	m_bUseRotInAverageOfPrevN=FALSE;
	m_bUseRotPerSec=FALSE;

	// sum of previous N frames :
	m_bDoSumOfPrevNFrames = FALSE;
	m_bKeepSumOfPrevNFrames = -1;	
	m_bSaveSumOfNFrames=FALSE;
	m_bAnalyzeSumOfPrevNFrames = FALSE;
	m_bCheckNormalTracksOnSumEvents = TRUE;
	m_fVetoRadiusOnSumFrame = 3;
	m_nNewLaplaceInSigmaOnSum = 3.5;
	m_nMaxLaplaceOnOtherInSigmaOnSum = 1.5;
	m_bRejectTracksOnSumedFrames = FALSE;
	m_nNumBackFramesForTracksOnSum = 50;

	// 20050407 - change to DEFAULT false - seems that rejects to many
	m_bCheckPrevOfMaxPixel = FALSE;

	// comparison to old frame :
	m_nCompareToOldFreqInSec = -1;
	


// general analysis section :
	m_bMoveWeighted=FALSE;
	m_bCorrectForRotation = TRUE;
	m_FrameDX = 0;
	m_FrameDY = 0;
	m_FrameDXPerSec = 0;
	m_FrameDYPerSec = 0;
	m_bUseShiftTotal=TRUE;
	m_RotCenterX = 0;
	m_RotCenterY = 0;
	m_RotValueDAlfa = 0;
	m_RotValueDAlfaPerSec = 0;

// frames shift :
	m_bUseFrameShift = FALSE;
	m_MinShiftToUse = 0.5;
	m_ForceAstroWhenBigShift = -1.00;
	m_ForceAstroWhenBigShiftRMS = -0.5;

// minimal requirements for calculating rotation
	m_MinShiftToCalcRot=0.5;
	m_MinTotalShiftToCalcRot=100.00;
	m_bDoNotUseRotation=FALSE;

	// dome status :
	m_CheckDomeStatusFreqInSec=-1;
   m_bStopOnDomeClose=FALSE;

//	m_bWaitingMode = FALSE;


	// good value is 10 ( sec ) , mount makes moves with modulo=10s 
	// and residue=5 sec, this means   :  5  , 15 , 25 etc ...
	// and daq will do it every 10 sec :    10 ,	 20 , 30 etc ..
	m_nDoFlatModeModulo = -1;
	m_nDoFlatModeResidue = 0;

	m_ObsMode=eNoMovingMode;
	m_DecObs=0;
	m_RAObs=0;
	m_HAObs=0;
	m_HorAzimuth = 0;
	m_HorAltitude = PI_2_VALUE;
	m_HorAzimCorr=0.00;
	m_HorAltCorr=0.00;
	m_RACorr=0.00;
	m_DecCorr=0.00;
	m_GeoLatitude = AstroAngle::deg2rad(52.25); // WARSAW - for the begining 
	m_GeoLongitude = AstroAngle::deg2rad(21.00); // WARSAW - for the begining
	m_GeoAltitude = AstroAngle::deg2rad(0.00); // WARSAW - for the begining
	m_SinOfGeoLatitude = sin( AstroAngle::deg2rad(52.25) );
	m_TimeZone = 2;
	m_nDarksToBeTaken=0;
	m_nDarksToBeTakenSav=0;
	m_bOverwriteOldDark=TRUE;
	m_nFramesToBeTaken=0;
	m_bWaitForFrame=FALSE;
	m_bReadAndSaveStatFile=FALSE;
	m_bDoChangeShutterTime = TRUE;
	

	// status file WWW :
	m_bDumpDAQStatus = FALSE;
	m_szDAQStatusFile = "/opt/pi/dev/pisys/status/daq.status";
	m_nDumpStatusFreqInSec = 60;

	// asas :
	m_bUseFastPhotoInAstro = FALSE;
	m_nMinStarCountToRunAstro = -1;
	m_bUseGoodFromOther = FALSE;
	m_bDoASASAstroOnStartup = FALSE;
	m_bDoASASPhotAstr = FALSE;
   m_nSaveReductFramesFreq = -1;
	m_bExecAstroOnFirst = FALSE;
	m_fASASPhotoThres = 5.00;
	m_fASASAstrometryFi = 0.00;
	m_nASASAstrometryOrd = 4;
   m_szASASStarCatalog = "$(DATADIR)/cat/act";
	m_szHIPStarCatalog  = "$(DATADIR)/cat/hip";
	m_szASASStarCatalog.env2str();
   m_bASASAstrometryVerb = 0;
   m_nASASAstrometryTry = 20;
	m_nASASAstrometryReTry = 5;
	m_nChangeToSilentAfterNGood = -1;
	m_nForceSynchroMaxCount = 5;
	m_fPixScale = 36.00;
	m_bUseAsasTransform = FALSE;
	m_eReverseForTransform = eReverseImageNone;
	m_nAsasBorderSize = 0;
	m_bAsasSubtrDark=TRUE;
	m_bAsasDivideByFlat=FALSE;
	m_fAsasError =0.15;
   m_fAsasFatalError=0.3;
   m_nAutoAstrometryFreq=200;
	m_nAutoAstrometryFreqInSec=-1;
	m_bExecAstroSynchroMode=FALSE;
	m_bKeepMagAndAstFromAstro = FALSE;
	m_bDoAstrometryInTakeNMode = FALSE;
	m_bDoPhotoInTakeNMode = FALSE;
	m_nWaitForAstrometryInSec = 0;
	m_bSaveGoodAstro = FALSE;
	m_bSaveAstroInTakeNMode = FALSE;
	m_BadAstroDX = 0;
	m_BadAstroDY = 0;

	// gigabit camera :
	m_szCameraIP = "100.100.100.1";
	m_CameraPortNo = 0;
	m_LocalPortNo = 0;
	m_CameraTimeoutMiliSec = 0;
	m_szEthCamLogFile = "cam_driver.log";
   m_szEthCamErrFile = "cam_driver.err";
	m_nEthCamLogLevel = 5;
	m_nEthCamErrLogLevel = 5;
	m_EthCamDataRetries = 0;
   m_EthCamCmdTimeout = 0;
   m_EthCamCmdRetries = 0;

	m_szDAQLogFile = GLOBAL_PI_LOG;
	m_bPiLogEnabled = FALSE;



	// AG - auto guide :
	m_CamIdxForAG_OnOff = 0;

	m_bUseAsasAstrometryFromFITS=FALSE;
	m_bUseCoordFromFITS=TRUE;	

	// ignore coord change on this camera 
	m_IgnoreCoordChangeOnCamera=-1;

	// fast-PI photometry :
	m_szFastPhotoCorrFile = "$(NDIR)/cfg/fast_photo.corr";

	// cataloging :
	m_bUseAverageMagnitudo = TRUE;

	m_bKeepLocalShift=FALSE;

	m_S0=0;
	m_S1=0;
	m_S2=0;
	m_S3=0;


// veto area:
	m_eVetoShape = shapeSquare;
	m_nVetoRedial=2;
	m_nVetoPointsCount=0;

	m_bIsCOIC=FALSE;
	m_bCheckForFlashes=TRUE;
	m_bAnalyseMaxFromPrevFrames=FALSE;
	m_bAnalyseSumAround=FALSE;
	m_FramesBack=0;
	m_MaxAllowedVal=30000;

	m_nIgnoreEdge=5;
	m_nIgnoreEdgeRight=5;
	m_nIgnoreEdgeLeft=5; 
	m_nIgnoreEdgeUp=5;   
	m_nIgnoreEdgeBottom=5;

	m_MinClusterSize=1;
	m_eNeigbShape = shapeStar;
	m_dNeighbRedial=0;
	m_nNeighbToSumCount=5;
// 	m_SumTreshold=0;		
	m_MaxNoiseLevel=200;
	m_ClusterIfNSigmaAboveBackgr=-10000;
	m_ConfTresholdPerPixel=200;
	m_ConfMaxPrevPerPixel=40;
	m_nPixelsAroundToConfirm=0;
	m_bUseClusterWithMore=FALSE;
// 	m_MaxOnPrevAllowed=1000;

// looking for variable objects :
	m_bCheckForSUPERNEW=FALSE;
	m_MinPrevPerPixel=2000;
	m_MinPrevTotal=10000;
	m_eBrigtheningCheckType=ePercentIncrease;
	m_IncreaseTreshADUPerPixel=500;
	m_IncreaseTreshADU=2500;
	m_IncreaseTreshPercent=0.1;

// homeopatic :
	m_bUseHomeoSameTreshForNewAndPrev=FALSE;
	m_bKeepHomeopaticSum=FALSE;
	m_nHomeoAverageOfPrevNFrames=0;
	m_bCheckHomeoRawCond=FALSE;
	m_HomeopaticFactor=(1.000/2.000);
	m_bCalcMaxNeighbHomeo=FALSE;
	m_nCalcMaxNieghbRedial=1;
	m_nCalcMaxForAboveNSigma=-1000; // means always use max of 3x3
	m_bCalcMaxForAboveNSigmaOnHomeo=FALSE;
	m_bStartHomeoWithFirstFrame=FALSE;

// laplacjan :
	m_bKeepLaplaceFrame=FALSE;
	m_bCheckLaplaceCondition=FALSE;
	m_bCalcLaplaceOfNEW=TRUE;
	m_bCheckLaplaceOnHomeoCondition=FALSE;
	m_eLaplaceType = eSinglePoint;
	m_nMaxLaplaceOnOther=20;
	m_nMaxLaplaceOfPrevAverage=70;
	m_nMinLaplaceOnOther=-1000;
	m_bMinLaplaceOnOtherEnabled=FALSE;
	m_nNewLaplace=100;
	m_nNewLaplaceInSigma=-1;
	m_nMaxLaplaceOnOtherInSigma=-1;
	m_bLaplaceUseSameTresh=FALSE;
	m_bConfirmLaplaceMedianMinus=FALSE;

	m_LaplacePlusCount=1;
	m_LaplaceMinusCount=4;

	m_eNotMoreThenNExceedsShape=shapeSquare;
	m_nNotMoreThenNExceedsRedial=2;
	m_nNotMoreThenNExceedsPointsCount=1;
	m_nCheckIfNoneOfNExceedsTresh=-1; // max allowed to exceed
	m_TreshForCheckIfNoneOfNExceedsInSigma=2.5;
	m_nNotMoreThenNExceedsBackFramesCount=0;


	m_bDarkFrame=FALSE;
	m_bFlatFrame=FALSE;
	m_bOffsetFrame=FALSE;
	m_bPixelMaxFrame=FALSE;

	// corrections :
	m_bShutterCorr = FALSE;

	 // image errors :
   m_bCheckFrame = FALSE;
	m_bRepairFrame = FALSE;
   m_ShiftCol = 15;
   m_ShiftVal = 5; // in sigma, in case >200 then in ADU - value of difference when shift is detected
	m_bRejectEventsNearShift = FALSE;


	m_szDarkFrameFile = "";
	m_szFlatFrameFile = "";

// confirmation :
	m_bLocalMaxReq=FALSE;
	m_bConfirmReq=FALSE;
	m_bCalcClusterReq=FALSE;
	m_ConfRedial=4;
	m_ConfTreshold=300;
	m_ConfShape=shapeSquare;

// max on prev in pipeline :
	m_bAnalyseMaxFromPrevInPipeline=FALSE;

// conditions types :
	m_bCheckDifference=FALSE;
	m_bCheckTreshAndMaxPrev=TRUE;

// conditions for confirmation :
	m_bConfCheckDifference=FALSE;
	m_bConfCheckTreshAndMaxPrev=FALSE;

// using same for confirmation and first check :
	m_bAlwaysUseSameTreshPerPixel=TRUE;
	m_SumTresholdForNewFrame=8545;
	m_SumOnPrevLessThen=8440;
	m_DiffTreshold=1000;

// difference check params :
	m_DiffCheckType=DiffPerPixel;
	m_DiffTreshSignal=2000;

// combined check params :
	m_bCombinedCheck=FALSE;
	m_UseDoubleCheckAboveADU=8000;
	m_UseDoubleCheckAboveADUPerPixel=1500;

// confirmation on next frames params :
	m_ConfirmEventsOnNextNFrames=0;
	m_ConfirmOnNextRadius=1;
	m_bRejectNotConfirmed=FALSE;
	m_MaxOnNextPerPixel=40;


// background :
	m_szRNoiseRowColDefFile="";
	m_bSubtractBackground=FALSE;
	m_BackgrMapXSize=10;
	m_BackgrMapYSize=10;
	m_BackgUpdateFreq=10;
	m_bBackgrDump=FALSE;
	m_eBackgrSubtrType=eNoneBackgrSubtr;
	m_CalcBackgrTable[MAX_LAPLACE_DEFINED];
	m_bCalcBackgrOfCurrLaplace=TRUE;

	m_AverageS1=1680;
	m_SigmaS1=40;


// TRACKS :
// verification if not track is found : TRACK :	
	m_bCheckTracks = FALSE;
		m_MaxEventRate=50;
	m_MaxChi2InTrack = 2.00;
	m_nMinEventNoInTrack = 5;
	m_nNumBackFramesForTracks=20;
	m_nNumBackFramesHasEventsForTracks=10;
	m_nMaxEventsToFitTrack=100;
	m_MaxChi2ForPointToMatchLine=36; // corresponding to 5 pixels distance 
	m_nMinDistOfEventsInTrack=-1; // originaly 3 
	m_bCheckVelocity = FALSE;
	m_fVelocityError = 0.20; // allowed error in pixels / sec 
	m_bLogTracks = TRUE;
	m_bLogFrameStat = FALSE;
	m_MaxChi2ForRejectIfMoreTrack=2.00;
	m_bCheckRejectIfMoreTracks=FALSE;
	m_nKeepRejectIfMoreTracksOnN=5;

	// old tracks
	m_nCheckFramesIfNotOlderThenNFrames=200;

	// PLANE tracks
	m_bCheckPlaneTracks = FALSE;
   m_MaxChi2InPlaneTrack = 2.00;
	m_MaxChi2InPlaneTrackToOld = 2.00;
	m_nNumBackFramesForPlaneTrack = 3;
	m_nMinPointNoOnPlaneTrack = 5;
	m_MinFramesOnPlaneTrack = 1;


	// line on single frame - add to all lines :
	m_bFitLineToSingleFrameToAll=FALSE;
	m_bFitLineToSingleFrame=FALSE;
	m_fChi2ForLineOnSingleFrame=2.00;

// 	m_MaxChi2InTrackToMatchLine=10;

	// SATELITES :
	m_bCheckIfSatelite = FALSE;
   m_bRejectSatelite = FALSE;
   m_szTleFile = "satelitesdb.tle";
   m_szQthFile = "satelitesdb.qth";
	m_nSatRejRadius = AstroAngle::deg2rad( 1.00 );
	m_nNotVisibleSatRejRadius = AstroAngle::arcsec2rad( 300.00 );

	// STAR :
	m_bCheckIfStar = FALSE;
   m_bRejectStars = FALSE;
   m_fStarRejectRadius = AstroAngle::arcsec2rad( 120.00 ); // 120 sec ~ 2 pixels -> 
	m_fCountBrighterThen = -1;
	m_bCheckStarsInTycho = FALSE;
	m_fStarCatMaxMag = 13.00;

	m_bRejectIfBigStarNearBy = FALSE;
   m_fBigStarRejectRadiusInArcSec = (60*10); // 10 pixels 
   m_fBigStarMaxMagnitudo = 4;

	m_bUseDB=FALSE;
	m_bUseODBC=FALSE;
   m_szDBName="pidb";
   m_szDBUser="pidb_user";
   m_szDBPass="pidb_user";
   m_szDBHost="localhost";
   m_szDB2Schema="MAIN";
	m_szScanDB="scan";
	m_bSaveEventsToDB=FALSE;
	m_bSaveFramesToDB=TRUE;
	m_bSaveVerifEventsToDB=FALSE;
	m_bSaveAllEventsToDB=FALSE;
	m_bSaveTracksToDB=FALSE;
	m_bSaveSumEventsToDB=FALSE;

	m_eRunType = eRunOnlineCoic;

	// queries for past GRB :
	m_bCheckForExternalsInDB = FALSE;
   m_fGrbRadiusInDeg = 1.00;
   m_TimePastToCheckInSec = 3600*24*10; // ten days TODO - check if enough
	m_SNRadiusInArcSec = 600; // 10 arcmin coic with exterenal SN 



// checking if not EDGE of BIG STAR :
	m_bCheckIfNotEdgeOfBigStar=FALSE;
	m_bEdgeOfBigByRawData=FALSE;
	m_nSigmaAboveMeanInRawCluster=3.00;
	m_nPrevFramesToCheckEdge=5;




// shift RA-DEC -> x,y transformations parameters :
	m_bIgnoreMissingTime=FALSE;
	m_bUseFrameTimeForShifts=FALSE;
	m_bShiftUsesAstroFormulas=FALSE;
	m_TransformCCDOrientation = 0.00;
		m_CCDAxixXOrientation=1;
		m_CCDAxixYOrientation=-1;
	m_TransformCCDFocus = 0.058;
	m_TransformCCDPixelSize = 0.000018;

	//
	m_szOrigin=DEFAULT_ORIGIN;
   m_szSite=DEFAULT_SITE;
   m_szInstrume=DEFAULT_INSTRUME;
   m_szCamOptic=DEFAULT_CAMOPTIC;
   m_szFilter=DEFAULT_FILTER;
   m_szObject=DEFAULT_OBJECT;
	m_szObserver=DEFAULT_OBSERVER;

	// 
	m_FOV=CANON_85mm_FOV_SIZE;
	m_FOVBorder=3.00;
	m_fFrameSize = ( m_FOV/2.0 )*CMyMathFunc::mysqrt(2) + 3.00; // with additional 2 degrees
																// - just to be sure



// auto calculations of tresholds :
	m_bAutoCalcTresh=FALSE;
	m_bCalcTresholdsByBackgrMap=FALSE;



// ERROR  handling :
	m_nExitOnToManyErrors = -1;

// reporting and tracing :
	m_nSaveCurrentPicture=0;
	m_bSaveEventDescOnly=FALSE;
	m_bSaveSuperNovaOnly=FALSE;
	m_nSaveFramesBeforeAndAfter=10;
	m_bSaveAverageParts = FALSE;
	m_bSaveOnlyGood = TRUE;
	m_bSaveFramesWithEvents=FALSE;
	m_bSaveFramesWithEventsOnSum=FALSE;
	m_bSaveEventSize=100;
	m_bDumpNewEvents=FALSE;
	m_bDumpAllEvents=FALSE;
	m_DumpEventsFreq=100;
	m_bDumpHomeoFrame=FALSE;
	m_bStdoutOn=TRUE;
	m_bGenNotFoundReport=TRUE;

	m_MaxStoredEventsFromSingleCam=100;
	m_bCheckTracksOnSingleCam = FALSE;

	m_bDumpAllLogToStdout = TRUE;
	m_bSameLogFiles = FALSE;
	m_szGenEventsLog="genevents_%d.log";
	m_szReGenEventsLog="regenevents_%d.log";
	m_szVerifiedEventsLog="verifiedevents_%d.log";
	m_szRunEventsLog="allevents_%d.log";

	m_EventsBufferSize=1000;

// flying objects cuts :
	m_bUseOriginalXY=FALSE;
	m_bSkipOverlaps=TRUE;
	m_OverlapRedial=10;
	m_eCheckIfMorePoint=eAfterCurrentFrameOnly;
	m_bSkipIfMoreThen=-1;
	m_bRejectIfMoreVerb=FALSE;
	m_bSkipIfMoreThenMinDist=2;
	m_bSkipIfMoreThenRedial=50;
	m_CheckEventShape=-1.00;
	m_bRejectBlackPixels=FALSE;
	m_bBlackPixelsIfNSigmaBelow=3;
	m_fBlackPixelsRatio=0.800;
	m_bCheckCenterInCluster=FALSE;
		m_MaxPixelsInClusterAllowed=10000; 
		m_MaxNumberOfEventsOnFrame=1000;
		m_MaxNumberOfEventsOnFrameAfterCoic=100;
	m_MaxNumberOfEventsOnFrameAfterTv = 3000;
	m_MinStarsToAccEvents = 500;

	// hot pixels :
	// m_HotPixels;
   m_bRejectHotByList=FALSE;
	m_bRejectHotPixelsByAverage = FALSE;
	m_nRejectHotPixelsTresholdInSigma = m_nNewLaplaceInSigma;


	m_szKernelModulesPath = "$(NDIR)/modules/";
	m_bReloadModuleOnStartup = FALSE;
	m_CCDIdentNo="SINGLE_PIPELINE_CCD_0";
	m_SearchDeviceStartNo=0;
	m_ReadFullStatusFreq=100;
	m_bDriverWriteAllToFITS=FALSE;
	m_bBuildFramesList=FALSE;
	m_eCompressFITS=eFITSComprNone;
	m_bAsynchroMode=FALSE;
	m_bParallelMode=FALSE;
	m_bTakeInSynchroMode=FALSE;
	m_szBaseFileNameFITS="aaa";


	// SHUTTER SETTINGS :
	m_DriverShutterTimeInSec=10;
	m_DriverShutterTimeInSecSaved=10;
	m_ShutterMode=2; 
	m_ShutterModeOriginal=m_ShutterMode;
	m_bHasBreakVoltages = 0;
	m_OpenBreakDelay   = 0x80;
	m_OpenBreakLength  = 0x70;
	m_CloseBreakDelay  = 0x60;
	m_CloseBreakLength = 0x55;

	m_bDriverReverseImage=eReverseImageNone;
	m_bUseCamID=FALSE;
		m_DriverAnalBinning=1;
	m_bDAQCommunicationON=FALSE;
		m_PortNo=DAQ_PORT;
	m_bFixDeviceOrder=FALSE;
	m_nPipelineRestartTimeout = 1200; // 10 minutes 
	m_CoordChangeToRestart = 0.2;
	m_MinStarCountToCloseShutter = -1;



// CORBA OPTIONS :
// in case standelone server - no CORBA deamons :
// 	m_CorbaOptions="corba_srv -ORBIIOPAddr inet:localhost.localdomain:12123 -ORBNoCodeSets -ORBIIOPBlocking";
	m_bCorbaWithNameService=TRUE;
// with Naming service :
	m_CorbaOptions=CORBA_NAME_SERVER_OPTIONS_SRV;
	m_CameraInterfaceNo = -1;


		m_SharedMemKey=1122;
		m_PipelineSafeBufferSize=1;
		m_IntervalBetweenFrames=12;
		m_bCoolingOnOff=DEFAULT_COOLING_ON_OFF;
		m_CCDTemp=DEFAULT_TEMPERATURE;
		m_bWaitForDesiredTemp=FALSE;
		m_WaitForDesiredTempInSec=900;
		m_eActionOnTempWaitFail = eTempWaitFail_IGNORE;
		m_CoolingToleranceDiff = 2;

// readout speeds :
		m_ReadoutSpeedHorizontal=DEFAULT_READOUT_SPEED_HORIZONTAL;
		m_ReadoutSpeedVertical=DEFAULT_READOUT_SPEED_VERTICAL;
		m_eMPPMode = DEFAULT_MPP_BC;

// gain and offset :
		m_ADCConfReg0=DEFAULT_ADC_CONF_REG0_SETTING;
		m_ADCConfReg1=DEFAULT_ADC_CONF_REG1_SETTING;
		m_ADCConfReg2=0;
		m_ADCConfReg3=0;
		m_ADCConfReg4=0;
		m_ADCConfReg5=0;
		m_ADCConfReg6=0;
		m_ADCConfReg7=0;
		m_bForceRegistersUsage=0;
		m_ADCOffset=DEFAULT_ADC_OFFSET;
		m_ADCOffsetCH2=m_ADCOffset;
		m_ADCGain=DEFAULT_ADC_GAIN;
		m_ADCGainCH2=m_ADCGain;
		m_ADCClampingBitOnOff = DEFAULT_ADC_CLAMPING;
	   m_eADCRange = DEFAULT_ADC_RANGE;
		m_eLNAGain = DEFAULT_LNA_GAIN;
	   m_bLensHitOnOff = FALSE;

		m_DriverMaxIterTimeout=DEFAULT_DRIVER_ITER_TIMEOUT;
		m_RetryFrameCount=5;
		m_DriverGetDataRetryCount=DEFAULT_RETRY_COUNT;
		m_MaxAllowedFramesLost=100;
		m_MaxCommErrorCountToEXIT = -1;
		m_SendBytesRetryCount = 1;
		m_ExitOnCommError = 0;

	m_eCAMType=eFileSimulator;

// simulator paremters :
	m_bRepeatSameImages = TRUE;


// comunnications with other parts of system :
	m_bPISysManagerON=FALSE;
	m_PISysManPortNo=DAQ_PI_SYS_PORT_NO;
	m_bAskMountForCoord=FALSE;

// MOUNT INFORMATION :
   m_MountID = -1;



// MonteCarlo optimalization :
	m_bKeepAllInMemory=FALSE;
	m_bKeepSamplesInMemory=TRUE;
	m_szSampleDir="";
	m_szListName="list";
	m_bPutScaledSamples=FALSE;

// frames dir/list
	m_szSampleFramesDir="aaa";
	m_szFramesListFile="frames_list";
	m_bCheckFramesListCount = TRUE;

	m_bCheckFrameOrder = FALSE;
	m_bSkipBadFrames = FALSE;
	m_SkipWarmImages = 1000;
	m_bAutoUpdateList = FALSE;
   m_nWaitTime = 60;
   m_nFrameListTimeout = 1800; // if no new frame in 1/2 hour then stop exit 



// performance optimalizations parametres:
	m_bKeepNeighbMap=FALSE;

// only for debuging reasons :
	m_bDumpClusters=FALSE;

// optimizations :
	m_bSuperOptimized=FALSE;
	m_bAverageOfPrev=FALSE;
	m_bHomeopatic=FALSE;

	// cosmic log - anty coic events :
   BOOL_T m_bLogAntyCoic;

// coicydence verification redial :
	m_nCoicRedial = 10.00;
	m_bCoicByRaDec = FALSE;
	m_bPutSampleByRaDec=FALSE;



// real analysis :
//	m_bReadFromDriver=FALSE;
	m_bIgnoreCamera = FALSE;


// FITS :
	m_DateObs = DATE_OBS;
	m_TimeObs = UT_START;
	m_DateTimeObsFormat = "'%d-%d-%dT%d:%d:%d.%d'";


// shifts calculation :
	m_szShiftsValuesFile="ccd_shifts_pipeline%d.cfg";
	m_nAutoShiftsMinStepsToIgnore=500;
		m_nAutoShiftsCalc=-1;
		m_nAutoShiftMatrixNo=2;
	m_MAXStarShapeLimit=0.2;
	m_MAXStarNSigmaAbove=3.00;

	m_MaxDiffSignalInSigma=0.70;
	m_MinDistBetweenTracedStars=20;
	m_UseNBackFrameForMatrix=50;
	m_bTraceOnAllFrames=FALSE;
	m_bUseTransformMatrix=FALSE;
	m_bTransformMatrixOK=FALSE;
	m_nMatrixNotFoundCount=0;
	m_MaxFramesWithOldMatrix=50;
	m_bDumpTraceStarToFile=TRUE;
	m_nMaxAllowedDrasticChng=5;
	m_nSkipNFramesAfterChange=0;
	m_bUseControllStar=FALSE;

	m_bDoNotAnalyse=FALSE;

	// SECOND LEVEL TRIGGER ( SLT )
	m_CheckClouds = FALSE;
   m_TypicalStarsCount = 40000;
	m_RejectFrameIfLess = 0.2;
	m_HoughTransformTresh = 4.00;
	m_HoughDistrMaxLimit = 2.00;
	m_HoughDistrTresh = 5.00;
	m_SmallHoughDistrTresh = 5.00;
	m_LapDiffMinRatio = -0.2;
	m_bCheckHoughOnSmall = FALSE;
	m_nSmallHoughSize = 20;
	m_MinStarCountOnParts = -1;
	m_bCheckCloudsOnPrev = FALSE;
	m_bCheckHoughOnRaw = FALSE;
	m_HoughTransformOnRawTresh = 4.00;
	m_HoughDistrOnRawTresh = 3.00;
	m_nMinBelowLimitToReject = 3;
	m_bCheckTrackRADEC = FALSE;
	

	// re-using old log files :
	m_bReadFirstLevelInfoFromLog = FALSE;
	m_bFromLogFile_WithFrames = FALSE;
	


// visu :
	m_bAutoCheckSize=FALSE;

	m_DayFramesCounter = 0;

// transformation :

	// coic radius in radians :
	m_nCoicRadiusInRad = GetPixToRad( m_nCoicRedial );	

	// matching stars to catalog :
   m_fMatchStarToCatRadiusInArcSec = 120; // or 120 
	m_fDiffMagToClaim = 1.00;
	m_fMinMagToClaimAlert = 10.00;
	m_nNextFramesToConfirm = 0;

	// Pointing paramteres : 
	m_bPointIO = TRUE;
	m_bAutoPointSWIFT = TRUE;
   m_bAutoPointINTEGRAL = TRUE;
	m_szPointingPrior = "SWIFT,INTEGRAL";


	// before was 22 and 24 - but something strange on night 20050309
	// this is too much - can be changed to 20 deg - select will be
	// faster :
	// m_fFrameSize=28.00; // image size in degrees  ~ sqrt(2)*15
	// for new canon 85 mm camera 10deg * sqrt(2) ~= 14.14 -> r ~= 20 deg moze byc
}


BOOL_T CCDConfig::RefreshParams()
{
	m_bInitialized = FALSE;
	return InitParams();
}

void CCDConfig::InitDefaults()
{
	InitBackgrFlagTable();
	mNotMoreThenNExceedsArea[0].x = 0;
	mNotMoreThenNExceedsArea[0].y = 0;
}

void CCDConfig::SaveParams()
{
	m_SaveCfgTab.clear();
	for(register int i=0;i<(GetParamFile()).GetParamTable().size();i++){
		m_SaveCfgTab.push_back( (GetParamFile()).GetParamTable()[i] );
	}
}

void CCDConfig::RestoreParams()
{
	GetParamFile().GetParamTable().clear();
	CCfgFile::CopyParamsTab( GetParamFile().GetParamTable(), m_SaveCfgTab );
	RefreshParams();
}

int CCDConfig::GetTrLevel()
{
	InitParams();
	return m_CCDTraceLevel;
}

mystring CCDConfig::GetTrFile()
{
	InitParams();
	return m_szCCDTraceFile;
}


BOOL_T CCDConfig::InitParams()
{
	if(!m_bInitialized){	
		_TRACE_PRINTF_0("CCDConfig::InitParams : checking for exit file stamp ...\n");fflush(0);

		// printf("In CCDConfig::InitParams\n");		

		const char* szTmp;		
		InitDefaults();

		// Pointing paramteres : 
		szTmp = GetParam("CCD_POINT_IO",TRUE);
		if(szTmp && szTmp[0]){
			m_bPointIO = ( atol( szTmp ) > 0 );
		}

		szTmp = GetParam("CCD_AUTOPOINT_SWIFT",TRUE);
		if(szTmp && szTmp[0]){
			m_bAutoPointSWIFT = ( atol( szTmp ) > 0 );
		}
		szTmp = GetParam("CCD_AUTOPOINT_INTEGRAL",TRUE);
		if(szTmp && szTmp[0]){
		   m_bAutoPointINTEGRAL = ( atol( szTmp ) > 0 );
		}
		szTmp = GetParam("CCD_POINTING_PRIOR",TRUE);
		if(szTmp && szTmp[0]){
			m_szPointingPrior = szTmp;
		}
		szTmp = GetParam("CCD_SPECIAL_TARGET",TRUE);
		if(szTmp && szTmp[0]){
			m_szSpecialTarget = szTmp;
		}

		szTmp = GetParam("CCD_CCDVIEW_INI_FILE",TRUE);
		if(szTmp && szTmp[0]){
			m_szDefaultIniFile = szTmp;
			m_szDefaultIniFile.env2str();
			_TRACE_PRINTF_4("ccdview ini file : %s\n",m_szDefaultIniFile.c_str());
		}
	
		szTmp = GetParam("CCD_ERROR_LOG_FILE",TRUE);
		if(szTmp && szTmp[0]){
			m_szErrorLogFile = szTmp;
		}

		szTmp = GetParam("CCD_SAVE_PARAM_CHANGES",TRUE);
		if(szTmp && szTmp[0]){
			m_bSaveParamChanges = ( atol( szTmp )>0 );
		}

		// run up to :
		szTmp = GetParam("CCD_RUN_UP_TO_DTM",TRUE);
		if(szTmp && szTmp[0]){
			m_szRunUpTo = szTmp;
			m_RunUpTo = get_runupto_dtm( m_szRunUpTo.c_str() );
		}

		// retry fit or no :
		szTmp = GetParam("CCD_DO_RETRY_FIT",TRUE);
		if(szTmp && szTmp[0]){
			CMyHisto::m_DoRetryFit = ( atol(szTmp)>0 );
		}


		// first param if not cfg file return :
		if(!GetParamFilePtr()){
			printf("CRITICAL ERROR : could not obtain config file pointer !!!\n");
			m_bInitialized=TRUE;
	      m_nInitCount++;
			return FALSE;
		}


		szTmp = GetParam("CCD_DAY_FRAME_COUNTER",TRUE);
		if(szTmp && szTmp[0]){			
			m_DayFramesCounter = atol( szTmp );
		}

		szTmp = GetParam("CCD_TRACE_LEVEL",TRUE);                                                                               
		if(szTmp && szTmp[0]){
			m_CCDTraceLevel = atol( szTmp );
			gCCDTrace.SetTrLevel( m_CCDTraceLevel );
		}

		szTmp = GetParam("CCD_TRACE",TRUE);
		if(szTmp && szTmp[0])
			m_szCCDTraceFile = szTmp;

		szTmp = GetParam("CCD_VISUALIZE_EVENTS",TRUE);
		if(szTmp && szTmp[0])
			m_bVisualizeEvents=( atol(szTmp)>0 );

		szTmp = GetParam("CCD_PART_SIZE",TRUE);
		if(szTmp && szTmp[0])
			m_PartSize = atol( szTmp );
		szTmp = GetParam("CCD_PART_SCALE",TRUE);
		if(szTmp && szTmp[0])
			m_PartScale = atol( szTmp );

		// basic :
		szTmp = GetParam("CCD_SIZE_X", TRUE );
		if( szTmp && szTmp[0] )
			m_SizeX = atol( szTmp );
		szTmp = GetParam( "CCD_SIZE_Y",TRUE );
		if( szTmp && szTmp[0] )
			m_SizeY = atol( szTmp );

		szTmp = GetParam( "CCD_NUMBER", TRUE );
		if(szTmp && szTmp[0] )
			m_nCamNo = atol( szTmp );

		szTmp = GetParam( "CCD_PIPELINE_SIZE", TRUE );
		if(szTmp && szTmp[0] )
			m_nPipelineSize = atol( szTmp );

		szTmp = GetParam( "CCD_CAMERA_INDEX", TRUE );
		if(szTmp && szTmp[0] ){
			m_CameraIndex = atol(szTmp);
		}

		// general analysis section :
		szTmp = GetParam("CCD_MOVE_WEIGHTED",TRUE);
		if(szTmp && szTmp[0])
			m_bMoveWeighted = ( atol(szTmp)>0 );

		szTmp = GetParam("CCD_CORRECT_FOR_ROTATION",TRUE);
		if(szTmp && szTmp[0])
			m_bCorrectForRotation = ( atol( szTmp ) >0 );

		szTmp = GetParam("CCD_MIN_SHIFT_TO_CALC_ROT",TRUE);
		if(szTmp && szTmp[0]){
			m_MinShiftToCalcRot = atof( szTmp );
		}

		szTmp = GetParam("CCD_MIN_TOTAL_SHIFT_TO_CALC_ROT",TRUE);
		if(szTmp && szTmp[0]){
			m_MinTotalShiftToCalcRot = atof( szTmp );
		}

		szTmp = GetParam("CCD_DO_NOT_USE_ROT",TRUE);
		if(szTmp && szTmp[0]){
			m_bDoNotUseRotation = ( atol( szTmp )>0 );
		}

// MS - correction on 2008-04-02 , not all paramters are refreshed in m_PipelineCfg
// when RefreshParams is called ! 
// in case it causes problems if can be restored , but only rotation paramteres should be there :
// This was rather obsolate, to read this paremters only when shifting is required, however
// if they are set to zeros, they should not make any troubles :
//		if(m_bCorrectForRotation){
			szTmp = GetParam("CCD_SINGLE_FRAME_DX", TRUE );
			if(szTmp && szTmp[0])			
				m_FrameDX = atof( szTmp );

			szTmp = GetParam("CCD_SINGLE_FRAME_DY", TRUE );
			if(szTmp && szTmp[0])			
				m_FrameDY = atof( szTmp );
			
			szTmp = GetParam("CCD_SINGLE_FRAME_DX_PER_SEC",TRUE);
			if(szTmp && szTmp[0])
				m_FrameDXPerSec = atof( szTmp );

			szTmp = GetParam("CCD_SINGLE_FRAME_DY_PER_SEC",TRUE);
			if(szTmp && szTmp[0])
				m_FrameDYPerSec = atof( szTmp );

			szTmp = GetParam("CCD_USE_SHIFT_TOTAL",TRUE);
			if(szTmp && szTmp[0])
				m_bUseShiftTotal = (atol(szTmp)>0);


			szTmp = GetParam("CCD_USE_FRAME_SHIFT",TRUE);
			if(szTmp && szTmp[0]){
				m_bUseFrameShift = (atol(szTmp)>0);
			}
			szTmp = GetParam("CCD_MIN_SHIFT_TO_USE",TRUE);
         if(szTmp && szTmp[0]){
				m_MinShiftToUse = atof( szTmp );
			}
			szTmp = GetParam("CCD_FORCE_ASTRO_WHEN_BIG_SHIFT",TRUE);
			if(szTmp && szTmp[0]){
				m_ForceAstroWhenBigShift = atof( szTmp );
			}
			szTmp = GetParam("CCD_FORCE_ASTRO_WHEN_BIG_SHIFT_RMS",TRUE);
			if(szTmp && szTmp[0]){
				m_ForceAstroWhenBigShiftRMS = atof( szTmp );
			}


			szTmp = GetParam("CCD_ROT_CENTER_X",TRUE);
			if(szTmp && szTmp[0])	
				m_RotCenterX = atof( szTmp );

			szTmp = GetParam("CCD_ROT_CENTER_Y",TRUE);
			if(szTmp && szTmp[0])
				m_RotCenterY = atof( szTmp );

			szTmp = GetParam("CCD_SINGLE_FRAME_D_ALFA",TRUE);
			if(szTmp && szTmp[0])
				m_RotValueDAlfa = atof( szTmp );

			szTmp = GetParam("CCD_SINGLE_FRAME_D_ALFA_PER_SEC",TRUE);
			if(szTmp && szTmp[0])
				m_RotValueDAlfaPerSec = atof( szTmp );

			// geograpfic coordinates :
			szTmp = GetParam("CCD_GEO_LATITUDE",TRUE);
			if(szTmp && szTmp[0])
				m_GeoLatitude = atof( szTmp );
			szTmp = GetParam("CCD_GEO_LONGITUDE",TRUE);
			if(szTmp && szTmp[0])
				m_GeoLongitude = atof( szTmp );

			szTmp = GetParam("CCD_GEO_LATITUDE_IN_DEG",TRUE);
			if(szTmp && szTmp[0]){
				m_GeoLatitude = AstroAngle::deg2rad( atof( szTmp ) );
			}
			szTmp = GetParam("CCD_GEO_LONGITUDE_IN_DEG",TRUE);
			if(szTmp && szTmp[0]){
				m_GeoLongitude = AstroAngle::deg2rad( atof( szTmp ) );
			}

			szTmp = GetParam("CCD_ALTITUDE",TRUE);
			if(szTmp && szTmp[0]){
				m_GeoAltitude = atof( szTmp );
			}

			// horizontal coordinates of observation :
			szTmp = GetParam("CCD_HOR_AZIMUTH",TRUE);
			if(szTmp && szTmp[0])
				m_HorAzimuth = atof( szTmp );
			szTmp = GetParam("CCD_HOR_ALTITUDE",TRUE);
			if(szTmp && szTmp[0])
				m_HorAltitude = atof( szTmp );

			szTmp = GetParam("CCD_HOR_AZIMUTH_IN_DEG",TRUE);
			if(szTmp && szTmp[0])
				m_HorAzimuth = AstroAngle::deg2rad( atof( szTmp ) );
			szTmp = GetParam("CCD_HOR_ALTITUDE_IN_DEG",TRUE);
			if(szTmp && szTmp[0])
				m_HorAltitude = AstroAngle::deg2rad( atof( szTmp ) );

			szTmp = GetParam("CCD_HOR_ALT_CORR",TRUE);
		   if(szTmp && szTmp[0])
		      m_HorAltCorr = atof( szTmp );
		   szTmp = GetParam("CCD_HOR_AZIM_CORR",TRUE);
		   if(szTmp && szTmp[0])
		      m_HorAzimCorr = atof( szTmp );

			szTmp = GetParam("CCD_RA_CORR",TRUE);
		   if(szTmp && szTmp[0])
		      m_RACorr = atof( szTmp );
		   szTmp = GetParam("CCD_DEC_CORR",TRUE);
		   if(szTmp && szTmp[0])
		      m_DecCorr = atof( szTmp );


			szTmp = GetParam("CCD_OBS_MODE",TRUE);
		   if(szTmp && szTmp[0])
      		m_ObsMode = (eObservationMode_T)(atol(szTmp));

		   szTmp = GetParam("CCD_DEC_OBS",TRUE);
		   if(szTmp && szTmp[0])
      		m_DecObs = atof(szTmp);

		   szTmp = GetParam("CCD_DEC_OBS_IN_DEG",TRUE);
		   if(szTmp && szTmp[0]){
      		m_DecObs = AstroAngle::deg2rad( atof(szTmp) );
			}

		   szTmp = GetParam("CCD_RA_OBS",TRUE);
		   if(szTmp && szTmp[0])
      		m_RAObs = atof(szTmp);

			szTmp = GetParam("CCD_RA_OBS_IN_HOURS",TRUE);
			if(szTmp && szTmp[0]){
				m_RAObs = AstroAngle::hours2rad( atof(szTmp) );
			}

			szTmp = GetParam("CCD_USE_ASAS_ASTROMETRY_FROM_FITS",TRUE);
			if(szTmp && szTmp[0]){
				m_bUseAsasAstrometryFromFITS = (atol(szTmp)>0);
			}

			szTmp = GetParam("CCD_USE_COORD_FROM_FITS",TRUE);
         if(szTmp && szTmp[0]){
				m_bUseCoordFromFITS = (atol(szTmp)>0);
			}

			szTmp = GetParam("CCD_USE_GOOD_FROM_OTHER",TRUE);
			if(szTmp && szTmp[0]){
				m_bUseGoodFromOther = (atol(szTmp)>0);
			}

			// pi cataloging :
			szTmp = GetParam("CCD_USE_AVERAGE_MAG",TRUE);
			if(szTmp && szTmp[0]){
				m_bUseAverageMagnitudo = (atol(szTmp)>0);
			}

			// fast photometry :				
			szTmp = GetParam("CCD_FAST_PHOTO_CORR_FILE",TRUE);
			if(szTmp && szTmp[0]){
				m_szFastPhotoCorrFile = szTmp;
			}

			szTmp = GetParam("CCD_USE_FAST_PHOTO_IN_ASTRO",TRUE);
         if(szTmp && szTmp[0]){
				m_bUseFastPhotoInAstro = (atol(szTmp)>0);
			}

//			szTmp = GetParam("CCD_FAST_PHOTO_TIMEOUT",TRUE);
//			if(szTmp && szTmp[0]){
//				CPiPhotometry::m_TimeLimit = atol(szTmp);
//			}

			szTmp = GetParam("CCD_MIN_STAR_COUNT_TO_RUN_ASTRO",TRUE);
         if(szTmp && szTmp[0]){
				m_nMinStarCountToRunAstro = atol( szTmp );
			}
	
			szTmp = GetParam("CCD_DO_ASAS_ASTRO_ON_STARTUP",TRUE);
			if(szTmp && szTmp[0]){
				m_bDoASASAstroOnStartup = (atol(szTmp)>0);
			}

			szTmp = GetParam("CCD_DO_ASAS_PHOT_ASTR",TRUE);
			if(szTmp && szTmp[0]){
				m_bDoASASPhotAstr = (atol(szTmp)>0);
			}
			szTmp = GetParam("CCD_SAVE_REDUCT_FRAMES_FREQ",TRUE);
			if(szTmp && szTmp[0]){
			   m_nSaveReductFramesFreq = atol(szTmp);
			}
			szTmp = GetParam("CCD_EXEC_ASTRO_ON_FIRST",TRUE);
			if(szTmp && szTmp[0]){
				m_bExecAstroOnFirst = (atol(szTmp)>0);
			}
			szTmp = GetParam("CCD_ASAS_PHOTO_THRES",TRUE);
			if(szTmp && szTmp[0]){				
				m_fASASPhotoThres = atof( szTmp );
			}
			szTmp = GetParam("CCD_ASAS_ASTROMETRY_FI",TRUE);
			if(szTmp && szTmp[0]){			
				m_fASASAstrometryFi = atof( szTmp );
			}
			szTmp = GetParam("CCD_ASAS_ASTROMETRY_ORD",TRUE);
			if(szTmp && szTmp[0]){
				m_nASASAstrometryOrd = atol( szTmp );
			}
			szTmp = GetParam("CCD_ASAS_STAR_CAT",TRUE);
			if(szTmp && szTmp[0]){
				m_szASASStarCatalog = szTmp;
			}
			szTmp = GetParam("CCD_HIP_STAR_CAT",TRUE);
			if(szTmp && szTmp[0]){
				m_szHIPStarCatalog = szTmp;
			}
			m_szHIPStarCatalog.env2str();
		   m_szASASStarCatalog.env2str();

			szTmp = GetParam("CCD_TYCHO_STAR_CAT",TRUE);
			if(szTmp && szTmp[0]){
				gStarCatTYCHO.m_StarCatCache.m_szCatBase = szTmp;
				gStarCatTYCHO.m_StarCatCache.m_szCatBase.env2str();
			}

			szTmp = GetParam("CCD_ASAS_ASTROMETRY_VERB",TRUE);
			if(szTmp && szTmp[0]){
			   m_bASASAstrometryVerb = atol(szTmp);
			}
			szTmp = GetParam("CCD_ASAS_ASTROMETRY_TRY",TRUE);
			if(szTmp && szTmp[0]){
			   m_nASASAstrometryTry = atol(szTmp);
			}
			szTmp = GetParam("CCD_ASAS_ASTROMETRY_RETRY",TRUE);
			if(szTmp && szTmp[0]){
				m_nASASAstrometryReTry = atol(szTmp);
			}
//			szTmp = GetParam("CCD_FIX_RA_CLOSE_ZERO",TRUE);
//			if(szTmp && szTmp[0]){
//				gFixRACloseZero = ( atol(szTmp)>0 );
//			}
			szTmp = GetParam("CCD_MIN_NUMBER_OF_MATCH_STARS",TRUE);
			if(szTmp && szTmp[0]){
				gMinNumberOfMatchesToAccept = atol(szTmp);
			}

			szTmp = GetParam("CCD_ASAS_ASTROMETRY_TIMEOUT",TRUE);
			if(szTmp && szTmp[0]){
				gMaxTimeForAstrometryInSec = atol(szTmp);
			}
			szTmp = GetParam("CCD_ASTRO_TO_SILENT_AFTER_N_GOOD",TRUE);
			if(szTmp && szTmp[0]){
				m_nChangeToSilentAfterNGood = atol(szTmp);
			}

			szTmp = GetParam("CCD_FORCE_SYNCHRO_MAX_COUNT",TRUE);
			if(szTmp && szTmp[0]){
				m_nForceSynchroMaxCount = atol(szTmp);
			}
			szTmp = GetParam("PIXSCALE",TRUE);
			if(szTmp && szTmp[0]){
				m_fPixScale = atof( szTmp );
			}
			szTmp = GetParam("CCD_USE_ASAS_TRANSFORM",TRUE);
			if(szTmp && szTmp[0]){
				m_bUseAsasTransform = ( atol(szTmp)>0 );
			}
			szTmp = GetParam("CCD_ASAS_REVERSE_FOR_TRANSFORM",TRUE);
			if(szTmp && szTmp[0]){
				m_eReverseForTransform = (eDriverReverseImage_T)atol(szTmp);
			}
			szTmp = GetParam("CCD_ASAS_BORDER_SIZE",TRUE);
			if(szTmp && szTmp[0]){							
				m_nAsasBorderSize = atol( szTmp );
			}
			szTmp = GetParam("CCD_KEEP_MAG_AND_AST_FROM_ASTRO",TRUE);
			if(szTmp && szTmp[0]){
				m_bKeepMagAndAstFromAstro = ( atol(szTmp)>0 );
			}
			
			szTmp = GetParam("CCD_ASAS_SUBTR_DARK",TRUE);
			if(szTmp && szTmp[0]){
				m_bAsasSubtrDark = ( atol(szTmp)>0 );
			}

			szTmp = GetParam("CCD_ASAS_DIVIDE_BY_FLAT",TRUE);
			if(szTmp && szTmp[0]){
				m_bAsasDivideByFlat = ( atol(szTmp)>0 );
			}

			szTmp = GetParam("MAX_AST_ERR",TRUE);
			if(szTmp && szTmp[0]){
				m_fAsasError = atof( szTmp );
			}
			szTmp = GetParam("MAX_AST_ERR_FATAL",TRUE);
			if(szTmp && szTmp[0]){
			   m_fAsasFatalError = atof( szTmp );
			}
			szTmp = GetParam("CCD_AUTO_ASTROMETRY_FREQ",TRUE);
			if(szTmp && szTmp[0]){
			   m_nAutoAstrometryFreq = atol( szTmp );
			}
			szTmp = GetParam("CCD_AUTO_ASTROMETRY_FREQ_IN_SEC",TRUE);
			if(szTmp && szTmp[0]){
			   m_nAutoAstrometryFreqInSec = atol( szTmp );
			}
			szTmp = GetParam("CCD_ASTROMETRY_SYNCHRO_MODE",TRUE);
			if(szTmp && szTmp[0]){
				m_bExecAstroSynchroMode = ( atol(szTmp)>0 );
			}

			szTmp = GetParam("CCD_SAVE_GOOD_ASTRO",TRUE);
         if(szTmp && szTmp[0]){
				m_bSaveGoodAstro = ( atol(szTmp)>0 );
			}

			szTmp = GetParam("CCD_SAVE_ASTRO_IN_TAKE_N_MODE",TRUE);
			if(szTmp && szTmp[0]){
				m_bSaveAstroInTakeNMode = ( atol(szTmp)>0 );
			}

			szTmp = GetParam("CCD_BAD_ASTRO_DX",TRUE);
			if(szTmp && szTmp[0]){
				m_BadAstroDX = atof(szTmp);
			}
			szTmp = GetParam("CCD_BAD_ASTRO_DY",TRUE);
			if(szTmp && szTmp[0]){
				m_BadAstroDY = atof(szTmp);
			}

			// AG - auto guide :
			szTmp = GetParam("CCD_CAM_IDX_FOR_AG_ON_OFF",TRUE);
			if(szTmp && szTmp[0]){
				m_CamIdxForAG_OnOff = atol(szTmp);
			}

			// for TAKE_N :
			szTmp = GetParam("CCD_ASTROMETRY_IN_TAKE_N_MODE",TRUE);
			if(szTmp && szTmp[0]){
				m_bDoAstrometryInTakeNMode = ( atol(szTmp)>0 );
			}
			szTmp = GetParam("CCD_PHOTO_IN_TAKE_N_MODE",TRUE);
         if(szTmp && szTmp[0]){
				m_bDoPhotoInTakeNMode = ( atol(szTmp)>0 );
			}
			szTmp = GetParam("CCD_WAIT_FOR_ASTROMETRY_IN_SEC",TRUE);
			if(szTmp && szTmp[0]){
				m_nWaitForAstrometryInSec = atol( szTmp );
			}

			szTmp = GetParam("CCD_IGNORE_COORD_CHANGE_ON_CAMERA",TRUE);
			if(szTmp && szTmp[0]){
				m_IgnoreCoordChangeOnCamera = atol( szTmp );
			}

			// time zone :
			szTmp = GetParam("CCD_TIME_ZONE",TRUE);
			if(szTmp && szTmp[0])
				m_TimeZone = atol( szTmp );


			szTmp = GetParam("CCD_KEEP_LOCAL_SHIFT",TRUE);
			if(szTmp && szTmp[0])
				m_bKeepLocalShift = (atol(szTmp)>0);
// MS : m_bCorrectForRotation - 2008-04-02 - if commented out
//		}


		szTmp = GetParam("CCD_CHECK_DOME_STATUS_FREQ_IN_SEC",TRUE);
		if(szTmp && szTmp[0]){
			m_CheckDomeStatusFreqInSec = atol(szTmp);
		}
		szTmp = GetParam("CCD_STOP_ON_DOME_CLOSE",TRUE);
		if(szTmp && szTmp[0]){
	   	m_bStopOnDomeClose = ( atol(szTmp)>0 );
		}


		szTmp = GetParam("CCD_WAITING_MODE",TRUE);
		if(szTmp && szTmp[0]){
			(CCDPipeline::m_WorkingMode).m_bWaitingMode = ( atol(szTmp)>0 );
		}

		szTmp = GetParam("CCD_DAQ_STARTED",TRUE);
      if(szTmp && szTmp[0]){
			(CCDPipeline::m_WorkingMode).m_bDAQStarted = ( atol(szTmp)>0 );
		}

		szTmp = GetParam("CCD_DO_FLAT_MODE_MODULO",TRUE);
		if(szTmp && szTmp[0]){
			m_nDoFlatModeModulo = atol( szTmp );			
		}
		szTmp = GetParam("CCD_DO_FLAT_MODE_RESIDUE",TRUE);
		if(szTmp && szTmp[0]){
			m_nDoFlatModeResidue = atol( szTmp );
		}
			

		szTmp = GetParam("CCD_READ_AND_SAVE_STAT_FILE",TRUE);
		if(szTmp && szTmp[0]){
			m_bReadAndSaveStatFile = ( atol(szTmp)>0 );
		}

		szTmp = GetParam("CCD_DUMP_DAQ_STATUS",TRUE);
		if(szTmp && szTmp[0]){
			m_bDumpDAQStatus = ( atol(szTmp)>0 );
		}
		szTmp = GetParam("CCD_DAQ_STATUS_FILE",TRUE);
		if(szTmp && szTmp[0]){
		   m_szDAQStatusFile = szTmp;
		}
		szTmp = GetParam("CCD_DUMP_STATUS_FREQ_IN_SEC",TRUE);
		if(szTmp && szTmp[0]){
			m_nDumpStatusFreqInSec = atol(szTmp);
		}			

	
		szTmp = GetParam("CCD_DO_CHANGE_SHUTTER_TIME",TRUE);
		if(szTmp && szTmp[0]){
			m_bDoChangeShutterTime = ( atol(szTmp)>0 );
		}

		szTmp = GetParam("CCD_OVERWRITE_OLD_DARKS",TRUE);
		if(szTmp && szTmp[0]){
			m_bOverwriteOldDark = ( atol(szTmp)>0 );
		}

		szTmp = GetParam("CCD_DARKS_TO_BE_TAKEN",TRUE);
		if(szTmp && szTmp[0]){
			m_nDarksToBeTaken	= atol( szTmp );
			m_nDarksToBeTakenSav = m_nDarksToBeTaken;
		}		

		szTmp = GetParam("CCD_WAIT_FOR_FRAME",TRUE);
		if(szTmp && szTmp[0]){
			m_bWaitForFrame = ( atol(szTmp)>0 );
		}
		
		szTmp = GetParam("CCD_FRAMES_TO_BE_TAKEN",TRUE);
		if(szTmp && szTmp[0]){
			m_nFramesToBeTaken = atol( szTmp );
		}

		szTmp = GetParam("CCD_HOMEO_USE_SAME_TRESH_FOR_NEW_AND_PREV",TRUE);
		if(szTmp && szTmp[0])
			m_bUseHomeoSameTreshForNewAndPrev = ( atol(szTmp)>0 );

		// always read :
		szTmp = GetParam("CCD_CHECK_SUM_AROUND",TRUE);
		if(szTmp && szTmp[0])
			m_bAnalyseSumAround = atol( szTmp );

		szTmp = GetParam("CCD_CHECK_MAX_FROM_PREV",TRUE);
		if(szTmp && szTmp[0])		
			m_bAnalyseMaxFromPrevFrames = atol( szTmp );

		szTmp = GetParam("CCD_CHECK_FOR_FLASHES",TRUE);
		if(szTmp && szTmp[0])
			m_bCheckForFlashes = ( atol(szTmp)>0 );

		szTmp = GetParam("CCD_BACK_FRAMES_USED",TRUE);
		if(szTmp && szTmp[0])			
			m_FramesBack = atol( szTmp );

		szTmp = GetParam("CCD_MAX_VALUE",TRUE);
		if(szTmp && szTmp[0])
			m_MaxAllowedVal = atol( szTmp );

		szTmp = GetParam("CCD_IGNORE_EDGE",TRUE);
		if(szTmp && szTmp[0]){
			m_nIgnoreEdge = atol( szTmp );
			m_nIgnoreEdgeRight = m_nIgnoreEdge;
			m_nIgnoreEdgeLeft = m_nIgnoreEdge;
			m_nIgnoreEdgeUp = m_nIgnoreEdge;
			m_nIgnoreEdgeBottom = m_nIgnoreEdge;
		}
		
		szTmp = GetParam("CCD_IGNORE_EDGE_RIGHT",TRUE);
		if(szTmp && szTmp[0])
			m_nIgnoreEdgeRight = atol(szTmp);
		szTmp = GetParam("CCD_IGNORE_EDGE_LEFT",TRUE);
		if(szTmp && szTmp[0])
			m_nIgnoreEdgeLeft = atol(szTmp);
		szTmp = GetParam("CCD_IGNORE_EDGE_UP",TRUE);
		if(szTmp && szTmp[0])
			m_nIgnoreEdgeUp = atol(szTmp);
		szTmp = GetParam("CCD_IGNORE_EDGE_BOTTOM",TRUE);
		if(szTmp && szTmp[0])
			m_nIgnoreEdgeBottom = atol(szTmp);


		// cluster calculations :
		szTmp = GetParam("CCD_CLUSTER_IF_N_SIGMA_ABOVE_BACKGR",TRUE);
		if(szTmp && szTmp[0])
			m_ClusterIfNSigmaAboveBackgr = atof( szTmp );
		
		szTmp = GetParam("CCD_CLUSTER_MIN_SIZE",TRUE);
		if(szTmp && szTmp[0])
			m_MinClusterSize = atol( szTmp );

		szTmp = GetParam("CCD_MAX_NOISE_LEVEL",TRUE);
		if(szTmp && szTmp[0])
			m_MaxNoiseLevel = atol( szTmp );

		szTmp = GetParam("CCD_CONFIRM_TRESHOLD_BY_PIXEL",TRUE);
		if(szTmp && szTmp[0])
			m_ConfTresholdPerPixel = atol( szTmp );

		
		szTmp = GetParam("CCD_CONFIRM_MAX_NOISE_PER_PIXEL",TRUE);
		if(szTmp && szTmp[0])
			m_ConfMaxPrevPerPixel = atol( szTmp );


		if(m_bUseHomeoSameTreshForNewAndPrev){
			m_ConfMaxPrevPerPixel = m_ConfTresholdPerPixel;
		}
			

		szTmp = GetParam("CCD_PIXELS_AROUND_TO_CONFIRM",TRUE );
		if(szTmp && szTmp[0])
			m_nPixelsAroundToConfirm = atol( szTmp );
		// m_MaxOnPrevAllowed = atol( GetParam("CCD_MAX_PREV") );
		szTmp = GetParam("CCD_DIFFERENCE_CHECK",TRUE);
		if(szTmp && szTmp[0])
			m_bCheckDifference = (atol(szTmp)>0);

		szTmp = GetParam("CCD_SEPERATE_TRESH_AND_MAX_PREV",TRUE);
		if(szTmp && szTmp[0])
			m_bCheckTreshAndMaxPrev = (atol(szTmp)>0);

		szTmp = GetParam("CCD_CONF_DIFFERENCE_CHECK",TRUE);
		if(szTmp && szTmp[0])
			m_bConfCheckDifference = atol( szTmp );

		szTmp = GetParam("CCD_CONF_SEPERATE_TRESH_AND_MAX_PREV",TRUE);
		if(szTmp && szTmp[0])
			m_bConfCheckTreshAndMaxPrev = atol( szTmp );

		// not required :
		szTmp = GetParam("CCD_USE_CLUSTER_WITH_MORE",TRUE);
		if(szTmp && szTmp[0])
	      m_bUseClusterWithMore = ( atol( szTmp )>0 );

		// combined check params :
		szTmp = GetParam("CCD_COMBINED_CHECK",TRUE);
		if(szTmp && szTmp[0])
			m_bCombinedCheck = (atol(szTmp) > 0);

		szTmp = GetParam("CCD_USE_DIFF_UP_TO_PER_PIXEL",TRUE);
		if(szTmp && szTmp[0])
			m_UseDoubleCheckAboveADUPerPixel =  atol( szTmp );

		// confirmation of events on next frames :
		szTmp = GetParam("CCD_CONFIRM_ON_N_NEXT_FRAMES",TRUE);
		if(szTmp && szTmp[0])
			m_ConfirmEventsOnNextNFrames = atol( szTmp );

		szTmp = GetParam("CCD_CONFIRM_ON_NEXT_RADIUS",TRUE);
      if(szTmp && szTmp[0])
			m_ConfirmOnNextRadius = atol( szTmp );


		szTmp = GetParam("CCD_REJECT_NOT_CONFIRMED", TRUE );
		if(szTmp && szTmp[0])
			m_bRejectNotConfirmed = ( atol( szTmp )>0 );

		szTmp = GetParam("CCD_MAX_ON_NEXT_PER_PIXEL",TRUE);
		if(szTmp && szTmp[0])
			m_MaxOnNextPerPixel = atol( szTmp );


		// difference check params :
		if(m_bCheckDifference){
			m_DiffCheckType = (eDiffCheckType_T) atol( GetParam("CCD_DIFF_TRESH_TYPE") );
			if(m_DiffCheckType == DiffTotal){
				m_DiffTreshSignal = atol( GetParam("CCD_DIFF_TRESH_SIGNAL") );
			}
		}

		// 
		szTmp = GetParam("CCD_USE_FOUND_POSITION",TRUE);
      if( szTmp && szTmp[0] ){
			m_bUseFoundPosition = ( atol( szTmp ) > 0 );
		}

		// analysis on sumed frames flag ( number of sumed frames )
		szTmp = GetParam("CCD_ON_SUMED_FRAMES",TRUE);
		if( szTmp && szTmp[0] ){
			m_bOnSumedFrames = atol( szTmp );
		}
		szTmp = GetParam("CCD_AUTO_CALC_SUM",TRUE);
		if( szTmp && szTmp[0] ){
			m_bAutoCalcSum = ( atol( szTmp ) > 0 );
		}
		szTmp = GetParam("CCD_REJECT_SINGLE_EVENTS",TRUE);
		if( szTmp && szTmp[0] ){
			m_bRejectSingleEvents = ( atol( szTmp ) > 0 );
		}

		szTmp = GetParam("CCD_ALWAYS_REREAD_SINGLE_FRAME_EVENTS",TRUE);
		if( szTmp && szTmp[0] ){
			m_bAlwaysReReadSingleFrameEvents	= ( atol( szTmp ) > 0 );
		}

		szTmp = GetParam("CCD_REJECT_SINGLE_TRACKS",TRUE);
		if( szTmp && szTmp[0] ){
			m_bRejectSingleTracks = ( atol( szTmp ) > 0 );
		}

		szTmp = GetParam("CCD_ALWAYS_REREAD_SINGLE_FRAME_TRACKS",TRUE);
		if( szTmp && szTmp[0] ){
			m_bAlwaysReReadSingleFrameTracks	= ( atol( szTmp ) > 0 );
		}


		// m_SumTreshold = atol( GetParam("CCD_SUM_ALG_TRESHOLD") );
		szTmp = GetParam("CCD_SIG_TRESH_PER_PIXEL_LEVEL_1",TRUE);
		if(szTmp && szTmp[0])
			m_TresholdPerPixel = atol( szTmp );

		szTmp = GetParam("CCD_MAX_PREV_PER_PIXEL_LEVEL_1",TRUE);
		if(szTmp && szTmp[0])
			m_MaxPrevPerPixel = atol( szTmp);


		if(m_bUseHomeoSameTreshForNewAndPrev){
			m_MaxPrevPerPixel = m_TresholdPerPixel;
		}

		szTmp = GetParam("CCD_LOCAL_MAX_REQ",TRUE);
		if(szTmp && szTmp[0])
			m_bLocalMaxReq = ( atol( szTmp )>0 );


		szTmp = GetParam("CCD_CONFIRM_REQ",TRUE);
		if(szTmp && szTmp[0])
			m_bConfirmReq = atol( szTmp );

		szTmp = GetParam("CCD_CALC_CLUSTER_REQ",TRUE);
		if(szTmp && szTmp[0])
			m_bCalcClusterReq = ( atol( szTmp )>0 );

		
		szTmp = GetParam("CCD_CONFIRM_REDIAL",TRUE);
		if(szTmp && szTmp[0])
			m_ConfRedial = atol( szTmp );

		szTmp = GetParam("CCD_CONFIRM_TRESHOLD",TRUE);
		if(szTmp && szTmp[0])
			m_ConfTreshold = atol( szTmp );

		szTmp = GetParam("CCD_CONFIRM_SHAPE",TRUE);
		if(szTmp && szTmp[0])			
			m_ConfShape = (eConfShape_T)atol( szTmp );

		szTmp = GetParam("CCD_VETO_SHAPE",TRUE);
		if(szTmp && szTmp[0])
			m_eVetoShape = (eConfShape_T)atol( szTmp );

		szTmp = GetParam("CCD_VETO_REDIAL",TRUE);
		if(szTmp && szTmp[0])
			m_nVetoRedial = atol( szTmp );
	
		szTmp = GetParam("CCD_NEIGHB_SHAPE",TRUE);
		if(szTmp && szTmp[0])
		   m_eNeigbShape = (eConfShape_T)( atol(szTmp) );

		szTmp = GetParam("CCD_NEIGHB_REDIAL",TRUE);
		if(szTmp && szTmp[0])
			m_dNeighbRedial = atof( szTmp );

		szTmp = GetParam("CCD_CHECK_MAX_FROM_PREV_N",TRUE);
		if(szTmp && szTmp[0])
			m_bAnalyseMaxFromPrevInPipeline = ( atol( szTmp )>0 );

		
		// average of previous frames :
		szTmp = GetParam("CCD_AVERAGE_OF_PREV_REJECT_MAX_AND_MIN",TRUE);
		if(szTmp && szTmp[0])
			m_bAverageOfPrevRejectMAXandMIN = (atol(szTmp)>0);

		szTmp = GetParam("CCD_AVERAGE_OF_PREV_N",TRUE);
		if(szTmp && szTmp[0])
			m_nMaxOfAverageOfPrevN = atol( szTmp );

		szTmp = GetParam("CCD_USE_ROT_IN_AVER_OF_PREV",TRUE);
		if(szTmp && szTmp[0])
			m_bUseRotInAverageOfPrevN = atol( szTmp );

		szTmp = GetParam("CCD_USE_ROT_PER_SEC",TRUE);
		if(szTmp && szTmp[0])
			m_bUseRotPerSec = (atol(szTmp)>0);

		// homeopatic section :
		szTmp = GetParam("CCD_KEEP_HOMEOPATIC_SUM",TRUE);
		if(szTmp && szTmp[0])
			m_bKeepHomeopaticSum = (atol(szTmp)>0);

		szTmp = GetParam("CCD_HOMEO_AVERAGE_OF_N",TRUE);
		if(szTmp && szTmp[0])
			m_nHomeoAverageOfPrevNFrames = atol(szTmp);

		szTmp = GetParam("CCD_KEEP_SUM_OF_PREV_N_FRAMES",TRUE);
		if(szTmp && szTmp[0])
			m_bKeepSumOfPrevNFrames = atol(szTmp);

		szTmp = GetParam("CCD_DO_SUM_OF_PREV_N_FRAMES",TRUE);
		if(szTmp && szTmp[0])
			m_bDoSumOfPrevNFrames = (atol(szTmp)>0);

		szTmp = GetParam("CCD_ANALYZE_SUM_OF_PREV_N_FRAMES",TRUE);
		if(szTmp && szTmp[0])
			m_bAnalyzeSumOfPrevNFrames = (atol(szTmp)>0);
		
		szTmp = GetParam("CCD_SAVE_SUM_OF_N_FRAMES",TRUE);
		if(szTmp && szTmp[0]){
			m_bSaveSumOfNFrames = (atol(szTmp)>0);
		}

		szTmp = GetParam("CCD_CHECK_NORMAL_TRACKS_ON_SUM_EVENTS",TRUE);
      if(szTmp && szTmp[0]){
         m_bCheckNormalTracksOnSumEvents = (atol(szTmp)>0);
      }

		szTmp = GetParam("CCD_VETO_RADIUS_ON_SUM_FRAME",TRUE);
		if(szTmp && szTmp[0]){
			m_fVetoRadiusOnSumFrame = atof(szTmp);
		}
		szTmp = GetParam("CCD_REJECT_TRACKS_ON_SUMED_FRAMES",TRUE);
		if(szTmp && szTmp[0]){
			m_bRejectTracksOnSumedFrames = (atol(szTmp)>0);
		}

		szTmp = GetParam("CCD_NUM_BACK_FRAMES_FOR_TRACK_ON_SUM",TRUE);
      if(szTmp && szTmp[0]){		
			m_nNumBackFramesForTracksOnSum = atol(szTmp);
		}

		szTmp = GetParam("CCD_COMPARE_TO_OLD_FREQ_IN_SEC",TRUE);
		if(szTmp && szTmp[0]){		
			m_nCompareToOldFreqInSec = atol(szTmp);
		}
			
		if(m_bKeepHomeopaticSum){
			m_bCheckHomeoRawCond = (atol(GetParam("CCD_CHECK_HOMEO_RAW_COND"))>0);
		}

		szTmp = GetParam("CCD_HOMEOPATIC_FACTOR",TRUE);
		if(szTmp && szTmp[0])
			m_HomeopaticFactor = atof( szTmp );
		Assert(m_HomeopaticFactor>0 && m_HomeopaticFactor<=1,"Homeopatic factor %f not handled",m_HomeopaticFactor);

		szTmp = GetParam("CCD_CALC_MAX_NEIGHB_HOMEO",TRUE);
		if(szTmp && szTmp[0])
			m_bCalcMaxNeighbHomeo = (atol(szTmp)>0);
		if(m_bCalcMaxNeighbHomeo){
			szTmp = GetParam("CCD_CALC_MAX_FOR_ABOVE_N_SIGMA",TRUE);
			if(szTmp && szTmp[0])
				m_nCalcMaxForAboveNSigma = atof(szTmp);
			szTmp = GetParam("CCD_CALC_MAX_FOR_ABOVE_N_SIGMA_ON_HOMEO",TRUE);
			if(szTmp && szTmp[0])
				m_bCalcMaxForAboveNSigmaOnHomeo = (atol(szTmp)>0);
			
			szTmp = GetParam("CCD_CALC_MAX_NEIGHB_REDIAL",TRUE);
			if(szTmp && szTmp[0])
				m_nCalcMaxNieghbRedial = atol(szTmp);
		}

		szTmp = GetParam("CCD_START_HOMEO_WITH_FIRST_FRAME",TRUE);
		if(szTmp && szTmp[0])
			m_bStartHomeoWithFirstFrame = (atol(szTmp)>0);

		// laplacjan frame :
		szTmp = GetParam("CCD_KEEP_LAPLACE_FRAME",TRUE);
		if(szTmp && szTmp[0])
			m_bKeepLaplaceFrame = (atol(szTmp)>0 );
		szTmp = GetParam("CCD_CHECK_LAPLACE_CONDITION",TRUE);
		if(szTmp && szTmp[0])
			m_bCheckLaplaceCondition = (atol(szTmp) >0 );


		szTmp = GetParam("CCD_CALC_LAPLACE_OF_NEW",TRUE);
		if(szTmp && szTmp[0])
			m_bCalcLaplaceOfNEW = (atol(szTmp)>0);



		szTmp = GetParam("CCD_MAX_LAPLACE_OF_PREV_AVERAGE",TRUE);
		if(szTmp && szTmp[0])
			m_nMaxLaplaceOfPrevAverage = atol( szTmp );

		szTmp = GetParam("CCD_CHECK_LAPLACE_ON_HOMEO_CONDITION",TRUE);
		if(szTmp && szTmp[0])
			m_bCheckLaplaceOnHomeoCondition = atol( szTmp );


		szTmp = GetParam("CCD_LAPLACE_TYPE",TRUE);
		if(szTmp && szTmp[0])
			m_eLaplaceType = (eLaplaceType_T)atol(szTmp);

		szTmp = GetParam("CCD_MAX_LAPLACE_ON_OTHER",TRUE);
		if(szTmp && szTmp[0])
			m_nMaxLaplaceOnOther = atol( szTmp );
						
		// checking if prev laplace is not less then min allowed :
		szTmp = GetParam("CCD_MIN_LAPLACE_ON_OTHER_ENABLED",TRUE);
		if(szTmp && szTmp[0])
			m_bMinLaplaceOnOtherEnabled = (atol(szTmp)>0);

		szTmp = GetParam("CCD_MIN_LAPLACE_ON_OTHER",TRUE);
		if(szTmp && szTmp[0])
			m_nMinLaplaceOnOther = atol( szTmp );

		szTmp = GetParam("CCD_NEW_LAPLACE_TRESHOLD",TRUE);
		if(szTmp && szTmp[0])	
			m_nNewLaplace = atol( szTmp );			

		szTmp = GetParam("CCD_MAX_LAPLACE_ON_OTHER_IN_SIGMA",TRUE);
      if(szTmp && szTmp[0]){
			// mystring szTmpStr=szTmp;
			// szTmpStr.env2str();
         m_nMaxLaplaceOnOtherInSigma = atof( szTmp );
		}
		szTmp = GetParam("CCD_NEW_LAPLACE_TRESHOLD_IN_SIGMA",TRUE);
      if(szTmp && szTmp[0]){
			// mystring szTmpStr=szTmp;
         // szTmpStr.env2str();
         m_nNewLaplaceInSigma = atof( szTmp );
		}
		m_nRejectHotPixelsTresholdInSigma = m_nNewLaplaceInSigma;

		szTmp = GetParam("CCD_CHECK_PREV_OF_MAX_PIXEL",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckPrevOfMaxPixel = (atol(szTmp)>0);
		}

		m_nNewLaplaceInSigmaOnSum  = m_nNewLaplaceInSigma;
		szTmp = GetParam("CCD_NEW_LAPLACE_TRESHOLD_IN_SIGMA_ON_SUM",TRUE);
		if(szTmp && szTmp[0]){
			m_nNewLaplaceInSigmaOnSum = atof(szTmp);
		}
		m_nMaxLaplaceOnOtherInSigmaOnSum = m_nMaxLaplaceOnOtherInSigma;
		szTmp = GetParam("CCD_MAX_LAPLACE_ON_OTHER_IN_SIGMA_ON_SUM",TRUE);
		if(szTmp && szTmp[0]){
			m_nMaxLaplaceOnOtherInSigmaOnSum = atof(szTmp);
		}



		szTmp = GetParam("CCD_LAPLACE_USE_SAME_TRESHOLD",TRUE);
		if(szTmp && szTmp[0])
			m_bLaplaceUseSameTresh = (atol(szTmp)>0);
		if(m_bLaplaceUseSameTresh)
			m_nMaxLaplaceOnOther = m_nNewLaplace;

		szTmp = GetParam("CCD_CONFIRM_LAPLACE_MEDIAN_MINUS",TRUE);
		if(szTmp && szTmp[0])
			m_bConfirmLaplaceMedianMinus = (atol(szTmp)>0);


		// parameters for checking if not more then N exceeds threshold :
		szTmp = GetParam("CCD_CHECK_IF_NONE_OF_N_EXCEEDS_TRESH",TRUE);
		if(szTmp && szTmp[0])
			m_nCheckIfNoneOfNExceedsTresh = atol( szTmp );
		szTmp = GetParam("CCD_TRESH_FOR_CHECK_IF_NONE_OF_N_EXCEEDS_IN_SIGMA",TRUE);
		if(szTmp && szTmp[0])
			m_TreshForCheckIfNoneOfNExceedsInSigma = atof( szTmp );
		szTmp = GetParam("CCD_NOT_MORE_THEN_N_EXEEDS_SHAPE",TRUE);
		if(szTmp && szTmp[0])
			m_eNotMoreThenNExceedsShape=(eConfShape_T)atol(szTmp);
		szTmp = GetParam("CCD_NOT_MORE_THEN_N_EXEEDS_REDIAL",TRUE);
		if(szTmp && szTmp[0])
			m_nNotMoreThenNExceedsRedial = atol(szTmp);			
		szTmp = GetParam("CCD_NOT_MORE_THEN_N_EXEEDS_BACK_FRAMES_COUNT",TRUE);
		if(szTmp && szTmp[0])
			m_nNotMoreThenNExceedsBackFramesCount = atol(szTmp);

		//	m_nNotMoreThenNExceedsPointsCount=1;
		//CLongPoint CCDConfig::mNotMoreThenNExceedsArea


		// variable objects section :
		szTmp = GetParam("CCD_MAX_SN_ON_FINAL",TRUE);
		if( szTmp && szTmp[0] ){
			m_MaxFinalSN = atol(szTmp);
		}

		szTmp = GetParam("CCD_SUPER_NOVA_TV_IN_SIGMA",TRUE);
		if( szTmp && szTmp[0] ){
			m_fSuperNovaTvInSigma = atof( szTmp );
		}
		szTmp = GetParam("CCD_SUPER_NOVA_TN_IN_SIGMA",TRUE);
		if( szTmp && szTmp[0] ){
   		m_fSuperNovaTnInSigma = atof( szTmp );
		}
		szTmp = GetParam("CCD_SUPER_NOVA_MIN_PREV_IN_SIGMA",TRUE);
		if( szTmp && szTmp[0] ){
			m_fSuperNovaMinPrevValue = atof( szTmp );
		}
		szTmp = GetParam("CCD_CHECK_FOR_SUPER_NEW",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckForSUPERNEW = ( atol(szTmp)>0 );
		}

		// supernovae simulation :
		szTmp = GetParam("CCD_ONLY_SUPERNEW_BACKGR",TRUE);
		if(szTmp && szTmp[0]){
			m_bOnlySuperNewBackgr = ( atol(szTmp)>0 );
		}
		szTmp = GetParam("CCD_SIMUL_SUPERNOVA_MIN_MAG",TRUE);
		if(szTmp && szTmp[0]){
			m_szMinMagSUPERNEW = szTmp;
		}
		szTmp = GetParam("CCD_SIMUL_SUPERNOVA_MAX_MAG",TRUE);
		if(szTmp && szTmp[0]){
			m_szMaxMagSUPERNEW = szTmp;
		}
		szTmp = GetParam("CCD_SIMUL_SUPERNOVA_INCREASE_RATIO",TRUE);
		if(szTmp && szTmp[0]){
			m_dRatioSUPERNEW = atof( szTmp );
		}

		szTmp = GetParam("CCD_DARK_FRAME",TRUE );
		if(szTmp && szTmp[0]){
			m_bDarkFrame = atol( szTmp );
		}

		szTmp = GetParam("CCD_FLAT_FRAME",TRUE);
		if(szTmp && szTmp[0]){		
	      m_bFlatFrame = ( atol(szTmp)>0 );
		}

		szTmp = GetParam("CCD_KEEP_PIXEL_MAX",TRUE);
		if(szTmp && szTmp[0]){
			m_bPixelMaxFrame = (atol( szTmp )>0);
		}

		szTmp = GetParam("CCD_SHUTTER_CORR",TRUE);
		if(szTmp && szTmp[0]){
			m_bShutterCorr = (atol( szTmp )>0);
		}

		szTmp = GetParam("CCD_CHECK_FRAME",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckFrame = (atol( szTmp )>0);
		}
		szTmp = GetParam("CCD_REPAIR_FRAME",TRUE);
      if(szTmp && szTmp[0]){
			m_bRepairFrame = (atol( szTmp )>0);
		}
		szTmp = GetParam("CCD_SHIFT_COL",TRUE);
		if(szTmp && szTmp[0]){
		   m_ShiftCol = atol( szTmp );
		}
		szTmp = GetParam("CCD_SHIFT_VAL",TRUE);
		if(szTmp && szTmp[0]){
		   m_ShiftVal = atol( szTmp );
		}
		szTmp = GetParam("CCD_REJECT_EVENTS_NEAR_SHIFT",TRUE);
		if(szTmp && szTmp[0]){
			m_bRejectEventsNearShift = (atol( szTmp )>0);
		}

		szTmp = GetParam("CCD_DARK_FRAME_FILE",TRUE);
		if(szTmp && szTmp[0])
			m_szDarkFrameFile = szTmp;

		szTmp = GetParam("CCD_FLAT_FRAME_FILE",TRUE);
		if(szTmp && szTmp[0]){
			m_szFlatFrameFile = szTmp;
		}
		

		// report dumping :
		szTmp = GetParam("CCD_DUMP_NEW_EVENTS",TRUE);
		if(szTmp && szTmp[0])
			m_bDumpNewEvents = atol( szTmp );

		szTmp = GetParam("CCD_DUMP_ALL_EVENTS",TRUE);
		if(szTmp && szTmp[0])
			m_bDumpAllEvents = atol( szTmp );
		szTmp = GetParam("CCD_DUMP_EVENT_FREQ",TRUE);
		if(szTmp && szTmp[0])
			m_DumpEventsFreq = atol( szTmp );
		szTmp = GetParam("CCD_SAVE_N_FRAMES_BEFORE_AND_AFTER",TRUE);
		if(szTmp && szTmp[0])
			m_nSaveFramesBeforeAndAfter = atol( szTmp );

		if( szTmp && szTmp[0] )
			m_bSaveOnlyGood = ( atol( szTmp )>0 );

		szTmp = GetParam("CCD_SAVE_FRAMES_WITH_EVENT",TRUE);
		if(szTmp && szTmp[0])
			m_bSaveFramesWithEvents = ( atol( szTmp )>0 );

		szTmp = GetParam("CCD_SAVE_FRAMES_WITH_EVENT_ON_SUM",TRUE);
		if(szTmp && szTmp[0])
			m_bSaveFramesWithEventsOnSum = ( atol( szTmp )>0 );

		szTmp = GetParam("CCD_SAVE_AVERAGE_PARTS",TRUE);
		if(szTmp && szTmp[0])
			m_bSaveAverageParts = ( atol( szTmp )>0 );
		
		szTmp = GetParam("CCD_SAVE_EVENT_SIZE",TRUE);
		if(szTmp && szTmp[0])
			m_bSaveEventSize = atol( szTmp );
		szTmp = GetParam("CCD_SAVE_SUPER_NOVA_ONLY",TRUE);
		if(szTmp && szTmp[0])
			m_bSaveSuperNovaOnly = ( atol( szTmp )>0 );
		szTmp = GetParam("CCD_SAVE_EVENT_DESC_ONLY",TRUE);
		if(szTmp && szTmp[0])
			m_bSaveEventDescOnly = ( atol( szTmp )>0 );

		szTmp = GetParam("CCD_SAME_LOG_FILES",TRUE);
		if(szTmp && szTmp[0])
			m_bSameLogFiles = ( atol( szTmp )>0 );

		szTmp = GetParam("CCD_DUMP_ALL_LOG_TO_STDOUT",TRUE);
		if(szTmp && szTmp[0]){
			m_bDumpAllLogToStdout = ( atol( szTmp )>0 );
		}

		szTmp = GetParam("CCD_VERIFIED_EVENTS_LOG",TRUE);
		if(szTmp && szTmp[0])
			m_szVerifiedEventsLog = szTmp;

		szTmp = GetParam("CCD_RUN_EVENT_LOG",TRUE);
		if(szTmp && szTmp[0])
			m_szRunEventsLog = szTmp;

		szTmp = GetParam("CCD_GEN_EVENT_LOG",TRUE);
		if(szTmp && szTmp[0])
			m_szGenEventsLog = szTmp;

		szTmp = GetParam("CCD_MAX_STORED_EVENTS_FROM_SINGLE_CAM",TRUE);
		if(szTmp && szTmp[0]){
			m_MaxStoredEventsFromSingleCam = atol( szTmp );
		}
		szTmp = GetParam("CCD_CHECK_TRACKS_ON_SINGLE_CAM",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckTracksOnSingleCam = ( atol( szTmp )>0 );
		}


		szTmp = GetParam("CCD_DUMP_HOMEO_FRAME",TRUE);
		if(szTmp && szTmp[0])
			m_bDumpHomeoFrame = ( atol( szTmp )>0 );

		szTmp = GetParam("CCD_GEN_NOT_FOUND_REPORT",TRUE);
		if(szTmp && szTmp[0])
			m_bGenNotFoundReport = ( atol( szTmp )>0 );
		
		
		// background section :
		szTmp = GetParam("CCD_CALC_S_DISTR",TRUE);
		if(szTmp && szTmp[0])
			m_CalcBackgrTable[(int)eRawS] = ( atol( szTmp )>0 );
		szTmp = GetParam("CCD_CALC_G1_DISTR",TRUE);
		if(szTmp && szTmp[0])
			m_CalcBackgrTable[(int)eSinglePoint] = ( atol( szTmp )>0 );
		szTmp = GetParam("CCD_CALC_G2_DISTR",TRUE);
		if(szTmp && szTmp[0])
			m_CalcBackgrTable[(int)eTwoPoints] = ( atol( szTmp )>0 );
		szTmp = GetParam("CCD_CALC_G4_DISTR",TRUE);
		if(szTmp && szTmp[0])
			m_CalcBackgrTable[(int)eFourPoints] = ( atol( szTmp )>0 );
		szTmp = GetParam("CCD_CALC_G5_DISTR",TRUE);
		if(szTmp && szTmp[0])
			m_CalcBackgrTable[(int)eFivePoints] = ( atol( szTmp )>0 );
		szTmp = GetParam("CCD_CALC_G54_DISTR",TRUE);
		if(szTmp && szTmp[0])
			m_CalcBackgrTable[(int)eFivePlusFourMin] = ( atol( szTmp )>0 );
		szTmp = GetParam("CCD_CALC_G985_DISTR",TRUE);
		if(szTmp && szTmp[0])
			m_CalcBackgrTable[(int)eNineEightFive] = ( atol( szTmp )>0 );
		szTmp = GetParam("CCD_CALC_G987_DISTR",TRUE);
		if(szTmp && szTmp[0])
			m_CalcBackgrTable[(int)eNineEightSevenVeryBig] = ( atol( szTmp )>0 );
		szTmp = GetParam("CCD_CALC_G84_DISTR",TRUE);
		if(szTmp && szTmp[0])
			m_CalcBackgrTable[(int)eEightFour] = ( atol( szTmp )>0 );
		szTmp = GetParam("CCD_CALC_G810_DISTR",TRUE);
		if(szTmp && szTmp[0])
			m_CalcBackgrTable[(int)eEightTen] = ( atol( szTmp )>0 );
		szTmp = GetParam("CCD_CALC_G58_DISTR",TRUE);
		if(szTmp && szTmp[0])
			m_CalcBackgrTable[(int)eFiveEight] = ( atol( szTmp )>0 );

		szTmp = GetParam("CCD_CALC_G412_DISTR",TRUE);
		if(szTmp && szTmp[0])
			m_CalcBackgrTable[(int)eFourTwelve] = ( atol( szTmp )>0 );
		szTmp = GetParam("CCD_CALC_G412F_DISTR",TRUE);
		if(szTmp && szTmp[0])
			m_CalcBackgrTable[(int)eFourTwelveFar] = ( atol( szTmp )>0 );

		szTmp = GetParam("CCD_CALC_BACKGR_OF_CURR_LAPLACE",TRUE);
		if(szTmp && szTmp[0])
			m_bCalcBackgrOfCurrLaplace = ( atol( szTmp )>0 );
		if(m_bCalcBackgrOfCurrLaplace){
			m_CalcBackgrTable[(int)m_eLaplaceType] = TRUE;
		}
		

		szTmp = GetParam("CCD_AVERAGE_S1",TRUE);
		if(szTmp && szTmp[0]){
			m_AverageS1 = atof( szTmp );	
			gAverageBacgroundG[eRawS] = m_AverageS1;
		}
		szTmp = GetParam("CCD_SIGMA_S1",TRUE);
		if(szTmp && szTmp[0]){
			m_SigmaS1 = atof( szTmp );
			gSigmaBacgroundG[eRawS] = m_SigmaS1;
		}

		szTmp = GetParam("CCD_AUTO_CALCULATE_TRESHOLDS",TRUE );
		if(szTmp && szTmp[0]){
			m_bAutoCalcTresh = ( atol( szTmp )>0 );
		}

		szTmp = GetParam("CCD_CALC_TRESHOLDS_BY_MAP",TRUE);
		if(szTmp && szTmp[0]){
			m_bCalcTresholdsByBackgrMap = ( atol( szTmp )>0 );
		}

		szTmp = GetParam("CCD_RNOISE_DEF_FILE",TRUE);
		if(szTmp && szTmp[0]){
			m_szRNoiseRowColDefFile = szTmp;
			if(!m_RNoiseRopColDefList.ReadFromFile( m_szRNoiseRowColDefFile.c_str() )){
				_TRACE_PRINTF_0("could not rnoise row-col def file : %s\n",m_szRNoiseRowColDefFile.c_str());
				// log to error log 
			}
		}

		if(GetKeepMapFlag()){
			szTmp = GetParam("CCD_BACKGR_MAP_X_SIZE",TRUE );
			if(szTmp && szTmp[0])
				m_BackgrMapXSize = atol( szTmp );
			
			szTmp = GetParam("CCD_BACKGR_MAP_Y_SIZE",TRUE );
			if(szTmp && szTmp[0])
				m_BackgrMapYSize = atol( szTmp );
		
			szTmp = GetParam("CCD_BACKGR_UPDATE_FREQ",TRUE);
			if(szTmp && szTmp[0])
				m_BackgUpdateFreq = atol( szTmp );

			szTmp = GetParam("CCD_BACKGR_DUMP",TRUE);
			if(szTmp && szTmp[0])
				m_bBackgrDump = ( atol( szTmp )>0 );
		}

		szTmp = GetParam("CCD_SUBTRACT_BACKGROUND",TRUE);
		if(szTmp && szTmp[0])
			m_bSubtractBackground = ( atol( szTmp )>0 );
		if(m_bSubtractBackground){
			szTmp = GetParam("CCD_BACKGROUND_SUBTR_TYPE",TRUE);
			if(szTmp && szTmp[0])
				m_eBackgrSubtrType = (eBackgrSubtrType_T)( atol(szTmp));
		}


		szTmp = GetParam("CCD_STDOUT_ON",TRUE);
		if(szTmp && szTmp[0]){
			m_bStdoutOn      = ( atol(szTmp)>0 );
		}
		

		szTmp = GetParam("CCD_EVENTS_BUFFER_SIZE",TRUE);
		if(szTmp && szTmp[0]){
			m_EventsBufferSize = atol(szTmp);
		}



		//-------------------------------------------------------------------------------
		// flying objects cuts :
		szTmp = GetParam("CCD_USE_ORIGINAL_XY",TRUE);
      if(szTmp && szTmp[0]){
			m_bUseOriginalXY = ( atol( szTmp )>0 );
		}

		szTmp = GetParam("CCD_SKIP_OVERLAPS",TRUE);
		if(szTmp && szTmp[0]){
			m_bSkipOverlaps = ( atol( szTmp )>0 );
		}
		szTmp = GetParam("CCD_OVERLAP_REDIAL",TRUE);
		if(szTmp && szTmp[0]){
			m_OverlapRedial = atol( szTmp );
		}

		szTmp = GetParam("CCD_MAX_NUMBER_OF_EVENTS",TRUE);
		if(szTmp && szTmp[0]){
			m_MaxNumberOfEventsOnFrame = atol( szTmp );
		}

		szTmp = GetParam("CCD_MAX_NUMBER_OF_EVENTS_AFTER_TV",TRUE);
		if(szTmp && szTmp[0]){
			m_MaxNumberOfEventsOnFrameAfterTv = atol( szTmp );
		}

		szTmp = GetParam("CCD_MAX_NUMBER_OF_EVENTS_AFTER_COIC",TRUE);
		if(szTmp && szTmp[0]){
			m_MaxNumberOfEventsOnFrameAfterCoic = atol( szTmp );
		}

		szTmp = GetParam("CCD_SKIP_IF_MORE_THEN",TRUE);
		if(szTmp && szTmp[0]){
			m_bSkipIfMoreThen = atol( szTmp );
		}
		szTmp = GetParam("CCD_SKIP_IF_MORE_THEN_REDIAL",TRUE);
		if(szTmp && szTmp[0]){
			m_bSkipIfMoreThenRedial = atol( szTmp );
		}

		szTmp = GetParam("CCD_SKIP_IF_MORE_THEN_MIN_DIST",TRUE);
		if(szTmp && szTmp[0]){
			m_bSkipIfMoreThenMinDist = atol( szTmp );
		}

		szTmp = GetParam("CCD_CHECK_IF_MORE_POINT",TRUE);
		if(szTmp && szTmp[0]){		
			m_eCheckIfMorePoint = (eCheckIfMorePoint_T)atol(szTmp);
		}

		szTmp = GetParam("CCD_CHECK_EVENT_SHAPE",TRUE);
		if(szTmp && szTmp[0]){
			m_CheckEventShape = atof( szTmp );
		}

		szTmp = GetParam("CCD_CHECK_CENTER_MASS_IN_CLUSTER",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckCenterInCluster = ( atol( szTmp )>0 );
		}

		// hot pixels :
		// CPointList m_HotPixels;
		szTmp = GetParam("CCD_REJECT_HOT_PIXELS_BY_LIST",TRUE);
		if(szTmp && szTmp[0]){
			m_bRejectHotByList = ( atol( szTmp )>0 );
		}		
		szTmp = GetParam("CCD_HOT_PIXELS_LIST",TRUE);
		if(szTmp && szTmp[0]){
		  	m_szHotList = szTmp;
			int nHotPixels = m_HotPixels.ReadFromFile( m_szHotList.c_str() );
			printf("Read %d hot pixels from cfg file %s\n",nHotPixels,m_szHotList.c_str() );
		}
		szTmp = GetParam("CCD_REJECT_HOT_PIXELS_BY_AVERAGE",TRUE);
      if(szTmp && szTmp[0]){
         m_bRejectHotPixelsByAverage = ( atol( szTmp )>0 );
      }
		szTmp = GetParam("CCD_REJECT_HOT_PIXELS_BY_AVERAGE_TRESH",TRUE);
		if(szTmp && szTmp[0]){
			m_nRejectHotPixelsTresholdInSigma = atof( szTmp );
		}


/*		szTmp = GetParam("CCD_LIST_BAD_COLUMNS",TRUE);
		if(szTmp && szTmp[0]){
			MyParser pars = szTmp;
			CMyStrTable bad_columns;
			pars.GetItems( bad_columns );
			for(int i=0;i<bad_columns.size();i++){
				int bad_col = atol( bad_columns[i].c_str() );
				m_BadPartsOfChip.Add( bad_col-1, 0, bad_col+1, 5000 );
			}
		}*/

		szTmp = GetParam("CCD_REJECT_BLACK_PIXELS",TRUE);
		if(szTmp && szTmp[0]){
			m_bRejectBlackPixels = ( atol( szTmp )>0 );

			// in case Black Pixels are checked also RAW distribution 
			// calculation is required :
			// m_CalcBackgrTable[(int)eRawS] = TRUE;
		}

		szTmp = GetParam("CCD_BLACK_PIXELS_IF_N_SIGMA_BELOW",TRUE);
		if(szTmp && szTmp[0]){
			m_bBlackPixelsIfNSigmaBelow = atof( szTmp );

			// in case Black Pixels are checked also RAW distribution
         // calculation is required :
         m_CalcBackgrTable[(int)eRawS] = TRUE;
		}
		
		szTmp = GetParam("CCD_BLACK_PIXELS_RATIO",TRUE);
		if(szTmp && szTmp[0]){
			m_fBlackPixelsRatio = atof( szTmp );
		}

		szTmp = GetParam("CCD_MAX_PIXELS_IN_CLUSTER_ALLOWED",TRUE);
		if(szTmp && szTmp[0])
			m_MaxPixelsInClusterAllowed = atol( szTmp );


		szTmp = GetParam("CCD_MIN_STARS_TO_ACC_EVENTS",TRUE);
		if(szTmp && szTmp[0]){
			m_MinStarsToAccEvents = atol( szTmp );
		}
		//--------------------------------------------------------------------------------
		szTmp = GetParam("CCD_CHANGE_SEED",TRUE );
		if(szTmp && szTmp[0]){
			gChangeSeed = ( atol( szTmp )>0 );
		}

		szTmp = GetParam("CCD_EXTERNAL_SEED",TRUE );
		if(szTmp && szTmp[0]){
			gExternalSeed = atol( szTmp );
		}

		// monte carlo parameters :
		szTmp = GetParam("CCD_MONTE_CARLO",TRUE );
		if(szTmp && szTmp[0]){
			m_bMC = ( atol( szTmp )>0 );
		}

		szTmp = GetParam("CCD_END_FILE_NAME",TRUE);
		if(szTmp && szTmp[0]){
			m_szEndFilePath = szTmp;
		}

		szTmp = GetParam("CCD_N_SAMPLES_TO_PUT_ON_FRAME",TRUE);
      if(szTmp && szTmp[0]){
			m_nSamplesToPutOnFrame = atol(szTmp);
		}


		szTmp = GetParam("CCD_SAVE_IMAGE_WITH_SAMPLE",TRUE);
		if(szTmp && szTmp[0]){
			m_bSaveImageWithSample = (atol(szTmp)>0);
		}

		szTmp = GetParam("CCD_PUT_SAMPLE_EVERY_N_FRAME",TRUE);
		if(szTmp && szTmp[0]){
			m_nPutSampleEveryNFrame = atol(szTmp);
		}

		szTmp = GetParam("CCD_CHECK_GEN_ONLY",TRUE);
		if(szTmp && szTmp[0])
			m_bCheckGenOnly = (atol(szTmp)>0);

		szTmp = GetParam("CCD_MC_GEN_IDENT_REDIAL",TRUE);
		if(szTmp && szTmp[0])
			m_bGenEventRedial = atol( szTmp );


		// testing with confirmation on next frames :
		// paramters to keep puting object in same place on subsequent
		// N frames :
		szTmp = GetParam("CCD_PUT_OBJ_ON_N_FRAMES",TRUE);
		if(szTmp && szTmp[0])
			m_nPutObjOnNFrames = atol( szTmp );

		szTmp = GetParam("CCD_SUM_ALL_METHODS_IN_NPAR_TEST",TRUE);
		if(szTmp && szTmp[0])
			m_bSumAllMethodsInNparTest = ( atol(szTmp)>0 );

		//--------------------------------------------------------------------------------------
		// SIMULATOR :
		szTmp = GetParam("CCD_SIMUL_USE_REAL_FITS_FILES",TRUE);
		if( szTmp && szTmp[0] ){
			m_bUseRealFITSFiles = ( atol(szTmp)>0 );
		}

		//--------------------------------------------------------------------------------------

		szTmp = GetParam("CCD_KEEP_ALL_IN_MEMORY_MC",TRUE);
		if(szTmp && szTmp[0])
			m_bKeepAllInMemory = ( atol(szTmp)>0 );

		szTmp = GetParam("CCD_KEEP_SAMPLES_IN_MEMORY",TRUE);
		if(szTmp && szTmp[0])
			m_bKeepSamplesInMemory = (atol(szTmp)>0);


		szTmp = GetParam("CCD_IMAGES_SAMPLE_DIR",TRUE);
		if(szTmp){
			m_szSampleDir = szTmp;
		}

		szTmp = GetParam("CCD_PUT_SAMPLE_CUT_OUT_IN_WIDTH",TRUE);
		if(szTmp && szTmp[0]){
			m_fPutSampleCutOutInWidth = atof(szTmp);
		}
		szTmp = GetParam("CCD_PUT_SECOND_SAMPLE_BY_NAME",TRUE);
		if(szTmp && szTmp[0]){
			m_bPutSecondSampleByName = ( atol(szTmp)>0 );
		}
		szTmp = GetParam("CCD_PUT_TAKEN_FROM_FRAME_IN_RANGE",TRUE);
		if(szTmp && szTmp[0]){		
			m_nPutTakenFromFrameInRange = atol( szTmp );
		}


		szTmp = GetParam("CCD_PUT_SCALED_SAMPLES",TRUE);
		if(szTmp && szTmp[0]){
			m_bPutScaledSamples = (atol(szTmp)>0);
		}

		szTmp = GetParam("CCD_LOG_SAMPLE_PUT",TRUE);
		if(szTmp && szTmp[0]){
			m_bLogSamplePut = (atol(szTmp)>0);
		}
		
		szTmp = GetParam("CCD_SAMPLE_LIST_BASENAME",TRUE);
		if(szTmp && szTmp[0]){
			m_szListName = szTmp;
		}

		// input dir/list :
		szTmp = GetParam( "CCD_FULL_FRAMES_SAMPLE_DIR",TRUE);
		if(szTmp){
			m_szSampleFramesDir = szTmp;
			// printf("Frames dir := %s\n",m_szSampleFramesDir.c_str());
		}

		szTmp = GetParam( "CCD_CHECK_FRAMES_LIST_COUNT",TRUE );
		if(szTmp && szTmp[0]){
			m_bCheckFramesListCount = (atol(szTmp)>0);
		}

		szTmp = GetParam( "CCD_FULL_FRAMES_LIST", TRUE );
		if(szTmp && szTmp[0]){
			m_szFramesListFile = szTmp;
		}else{
			SetParam( "CCD_FULL_FRAMES_LIST", m_szFramesListFile);
		}
		
		szTmp = GetParam( "CCD_SKIP_IF_OBJECT_MATCHES", TRUE );
		if(szTmp && szTmp[0]){
			m_szSkipIfObjectMatches = szTmp;
		}


		szTmp = GetParam( "CCD_SKIP_BAD_FRAMES",TRUE );
      if(szTmp && szTmp[0]){
			m_bSkipBadFrames = ( atol(szTmp)>0 );
		}

		szTmp = GetParam( "CCD_SKIP_WARM_IMAGES",TRUE );
		if(szTmp && szTmp[0]){
			m_SkipWarmImages = atof(szTmp);
		}

		szTmp = GetParam( "CCD_CHECK_FRAMES_ORDER",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckFrameOrder = ( atol(szTmp)>0 );
		}

		szTmp = GetParam( "CCD_AUTO_UPDATE_LIST",TRUE);
		if(szTmp && szTmp[0]){
			m_bAutoUpdateList = (atol(szTmp)>0);
		}
		szTmp = GetParam( "CCD_FRAME_WAIT_TIME",TRUE );
		if(szTmp && szTmp[0]){
		   m_nWaitTime = atol(szTmp);
		}
		szTmp = GetParam( "CCD_FRAME_LIST_TIMEOUT", TRUE );
		if(szTmp && szTmp[0]){
	   	m_nFrameListTimeout = atol(szTmp);
		}
			
		
		// performance optimalizations parametres:
		szTmp = GetParam( "CCD_KEEP_NEIGHB_MAP", TRUE );
		if(szTmp && szTmp[0]){
			m_bKeepNeighbMap = ( atol(szTmp)>0 );
		}

		szTmp = GetParam( "CCD_ALWAYS_USE_SAME_TRESH_PER_PIXEL", TRUE );
		if(szTmp && szTmp[0]){
			m_bAlwaysUseSameTreshPerPixel = ( atol(szTmp)>0 );
		}


		// new parameters for excluding defintions of algorithms :
		szTmp = GetParam("CCD_SUPER_OPTIMIZED",TRUE);
		if(szTmp && szTmp[0])
			m_bSuperOptimized = (atol(szTmp)>0);
		if(m_bSuperOptimized){
			szTmp = GetParam("CCD_AVERAGE_OF_PREV_ALG",TRUE);
			if(szTmp && szTmp[0])
				m_bAverageOfPrev = (atol(szTmp)>0);
			szTmp = GetParam("CCD_HOMEOPATIC_ALG",TRUE);
			if(szTmp && szTmp[0])
				m_bHomeopatic = (atol(szTmp)>0);
		}

// ERROR handling
		szTmp = GetParam("CCD_EXIT_ON_TOO_MANY_ERRORS",TRUE);
		if(szTmp && szTmp[0])
			m_nExitOnToManyErrors = atol(szTmp);		

		// common params :
		// tracing, printfs :
		szTmp = GetParam("CCD_DEBUG_TRACKS",TRUE);
		if(szTmp && szTmp[0])
			gDebugTracks = atol(szTmp);

		szTmp = GetParam("CCD_PRINTF_LEVEL",TRUE);
		if(szTmp && szTmp[0]){
			gPrintfLevel = atol(szTmp);
			gSatlibTraceLevel = gPrintfLevel;
		}

		szTmp = GetParam("CCD_DO_DUMP_BAD_FIT",TRUE);
		if(szTmp && szTmp[0])
			gDoDumpBadFit = ( atol(szTmp)>0 );

		szTmp = GetParam("CCD_DO_DUMP_ALL_HISTO",TRUE);
		if(szTmp && szTmp[0])
			gDoDumpAllHisto = ( atol(szTmp)>0 );

		CalcShapes();
		RecalculateParams();

		szTmp = GetParam("CCD_DEBUG_DUMP_CLUSTERS",TRUE);
		if(szTmp && szTmp[0]){
			m_bDumpClusters = ( atol( szTmp )>0 );
		}

		//szTmp = GetParam("CCD_READ_FRAMES_FROM_DRIVER",TRUE);
		//if(szTmp && szTmp[0]){
		//	m_bReadFromDriver = ( atol( szTmp )>0 );
		//}

		// real analysis :
		szTmp = GetParam("CCD_IGNORE_CAMERA",TRUE);
		if(szTmp && szTmp[0]){
			m_bIgnoreCamera = ( atol( szTmp )>0 );
		}

		szTmp = GetParam("CCD_LOG_ANTY_COIC",TRUE);
		if(szTmp && szTmp[0]){
			m_bLogAntyCoic = ( atol( szTmp )>0 );
		}

		szTmp = GetParam("CCD_COICYDENCE_REDIAL",TRUE);
		if(szTmp && szTmp[0]){
			m_nCoicRedial = atof(szTmp);
		}

		szTmp = GetParam("CCD_COIC_BY_RA_DEC",TRUE);
		if(szTmp && szTmp[0]){
			m_bCoicByRaDec = ( atol( szTmp )>0 );
		}

		szTmp = GetParam("CCD_PUT_SAMPLE_BY_RA_DEC",TRUE);
		if(szTmp && szTmp[0]){
			m_bPutSampleByRaDec = ( atol( szTmp )>0 );
		}

		BOOL_T bCoicRadiusInRad=FALSE;
		szTmp = GetParam("CCD_COIC_RADIUS_IN_RAD",TRUE);
		if(szTmp && szTmp[0]){
			m_nCoicRadiusInRad = atof( szTmp );
			bCoicRadiusInRad = TRUE;
		}		
		szTmp = GetParam("CCD_COIC_RADIUS_IN_SEC",TRUE);
		if(szTmp && szTmp[0]){
			m_nCoicRadiusInRad = AstroAngle::arcsec2rad( atof( szTmp ) );
			bCoicRadiusInRad = TRUE;			
		}	

		szTmp = GetParam("CCD_SAMPLE_XY_FROM_LIST",TRUE);
		if(szTmp && szTmp[0]){
			m_szSampleXYFromListFile = szTmp;
		}	

		
		szTmp = GetParam("CCD_SAMPLE_PUT_RETRY",TRUE);
		if(szTmp && szTmp[0]){
			m_nSamplePutRetry = atol( szTmp );
		}
		

		// matching stars to catalog :
		szTmp = GetParam("CCD_NEXT_FRAMES_TO_CONFIRM",TRUE);
      if(szTmp && szTmp[0]){
			m_nNextFramesToConfirm = atol( szTmp );
		}

/*		szTmp = GetParam("CCD_FRAME_SIZE_IN_DEG",TRUE);
      if(szTmp && szTmp[0]){
			m_fFrameSize = atof( szTmp );
		}*/

		szTmp = GetParam("CCD_MATCH_STAR_TO_CAT_RADIUS_IN_ARCSEC",TRUE);
		if(szTmp && szTmp[0]){
		   m_fMatchStarToCatRadiusInArcSec = atof( szTmp );
		}

		szTmp = GetParam("CCD_DIFF_MAG_TO_CLAIM",TRUE);
		if(szTmp && szTmp[0]){
			m_fDiffMagToClaim = atof( szTmp );
		}
	
		szTmp = GetParam("CCD_MIN_MAG_TO_CLAIM_ALERT",TRUE);
		if(szTmp && szTmp[0]){
			m_fMinMagToClaimAlert = atof( szTmp );
		}


		//----------------------------------------------------------------
		// verification if not track is found : TRACK :
		szTmp = GetParam("CCD_CHECK_TRACKS",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckTracks = ( atol( szTmp )>0 );
		}
		szTmp = GetParam("CCD_MIN_DIST_OF_EVENTS_IN_TRACK",TRUE);
		if(szTmp && szTmp[0]){
			m_nMinDistOfEventsInTrack = atof( szTmp );
		}
		szTmp = GetParam("CCD_MAX_EVENT_RATE_FOR_TRACKS",TRUE);
		if(szTmp && szTmp[0]){	
			m_MaxEventRate = atol( szTmp );
		}
		szTmp = GetParam("CCD_MAX_CHI2_IN_TRACK",TRUE);
		if(szTmp && szTmp[0]){
			m_MaxChi2InTrack = atof( szTmp );
		}
		szTmp = GetParam("CCD_MAX_CHI2_FOR_POINT_TO_MATCH_LINE",TRUE);
		if(szTmp && szTmp[0]){
			m_MaxChi2ForPointToMatchLine = atof( szTmp );
		}
		szTmp = GetParam("CCD_MIN_EVENTS_NO_IN_TRACK",TRUE);
		if(szTmp && szTmp[0]){
			m_nMinEventNoInTrack = atol( szTmp );
		}
		szTmp = GetParam("CCD_NUM_BACK_FRAMES_FOR_TRACK",TRUE);
		if(szTmp && szTmp[0]){
			m_nNumBackFramesForTracks = atol( szTmp );
		}
		szTmp = GetParam("CCD_NUM_BACK_FRAMES_HAS_EVENTS_FOR_TRACKS",TRUE);
		if(szTmp && szTmp[0]){
			m_nNumBackFramesHasEventsForTracks = atol( szTmp );
		}
		szTmp = GetParam("CCD_MAX_EVENTS_TO_FIT_TRACK",TRUE);
		if(szTmp && szTmp[0]){
			m_nMaxEventsToFitTrack = atol( szTmp );
		}
		szTmp = GetParam("CCD_LOG_TRACKS",TRUE);
		if(szTmp && szTmp[0]){
			m_bLogTracks = ( atol( szTmp )>0 );
		}
		szTmp = GetParam("CCD_LOG_FRAME_STAT",TRUE);
		if(szTmp && szTmp[0]){
			m_bLogFrameStat = ( atol( szTmp )>0 );
		}
		szTmp = GetParam("CCD_CHECK_VELOCITY_FOR_TRACK",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckVelocity = ( atol( szTmp )>0 );
		}
		szTmp = GetParam("CCD_CHECK_VELOCITY_ERROR_FOR_TRACK",TRUE);
		if(szTmp && szTmp[0]){
			m_fVelocityError = atof( szTmp );		
		}
		szTmp = GetParam("CCD_MIN_VELO_TO_CHECK",TRUE);
		if(szTmp && szTmp[0]){
			gMinVeloToCheck = atof( szTmp );
		}	

		szTmp = GetParam("CCD_CHECK_REJECT_IF_MORE_TRACKS",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckRejectIfMoreTracks = ( atol( szTmp )>0 );
		}
		szTmp = GetParam("CCD_REJECT_IF_MORE_TRACKS_MAX_CHI2",TRUE);
		if(szTmp && szTmp[0]){
			m_MaxChi2ForRejectIfMoreTrack = atof( szTmp );
		}
		szTmp = GetParam("CCD_KEEP_REJECT_IF_MORE_TRACKS_ON_N_FRAMES",TRUE);
		if(szTmp && szTmp[0]){
			m_nKeepRejectIfMoreTracksOnN = atol( szTmp );
		}

		// track on single cam :
		szTmp = GetParam("CCD_FIT_LINE_TO_SINGLE_FRAME",TRUE);
		if(szTmp && szTmp[0]){
			m_bFitLineToSingleFrame = ( atol( szTmp )>0 );
		}
		szTmp = GetParam("CCD_FIT_LINE_TO_SINGLE_FRAME_TO_ALL",TRUE);
      if(szTmp && szTmp[0]){
			m_bFitLineToSingleFrameToAll = ( atol( szTmp )>0 );
		}
		szTmp = GetParam("CCD_CHI2_FOR_LINE_ON_SINGLE_FRAME",TRUE);
		if(szTmp && szTmp[0]){
		   m_fChi2ForLineOnSingleFrame = atof( szTmp );
		}

		szTmp = GetParam("CCD_CHECK_OLD_TRACKS_NOT_OLDER_THEN",TRUE);
		if(szTmp && szTmp[0]){
			m_nCheckFramesIfNotOlderThenNFrames = atol( szTmp );
		}

		szTmp = GetParam("CCD_SHOW_CHI2_3POINTS",TRUE);
		if(szTmp && szTmp[0]){
			gShowChi2Points3 = ( atol( szTmp )>0 );
		}
		szTmp = GetParam("CCD_SHOW_NPOINTS_VS_CHI2",TRUE);
      if(szTmp && szTmp[0]){
         gShowChi2_NPoints = ( atol( szTmp )>0 );
      }

		szTmp = GetParam("CCD_SHOW_CHI2_OF_ADDED",TRUE);
		if(szTmp && szTmp[0]){
			gShowChi2OfAdded = ( atol( szTmp )>0 );
		}

		szTmp = GetParam("CCD_CHECK_PLANE_TRACKS",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckPlaneTracks = ( atol( szTmp )>0 );
		}
		szTmp = GetParam("CCD_MAX_CHI2_IN_PLANE_TRACK",TRUE);
		if(szTmp && szTmp[0]){
	   	m_MaxChi2InPlaneTrack = atof( szTmp );
		}
		szTmp = GetParam("CCD_MAX_CHI2_IN_PLANE_TRACK_TO_OLD",TRUE);
		if(szTmp && szTmp[0]){
			m_MaxChi2InPlaneTrackToOld = atof( szTmp );
		}
		szTmp = GetParam("CCD_NUM_BACK_FRAMES_FOR_PLANE_TRACK",TRUE);
		if(szTmp && szTmp[0]){
			m_nNumBackFramesForPlaneTrack = atol( szTmp );
		}
		szTmp = GetParam("CCD_MIN_POINT_NO_ON_PLANE_TRACK",TRUE);
		if(szTmp && szTmp[0]){
			m_nMinPointNoOnPlaneTrack = atol( szTmp );
		}
		szTmp = GetParam("CCD_MIN_FRAMES_ON_PLANE_TRACK",TRUE);
		if(szTmp && szTmp[0]){
			m_MinFramesOnPlaneTrack = atol( szTmp );
		}
		//----------------------------------------------------------------		

		//----------------------------------------------------------------
		// check if not event on EDGE of BIG STAR :
		szTmp = GetParam("CCD_CHECK_IF_NOT_EDGE_OF_BIG_STAR",TRUE);
		if(szTmp && szTmp[0]){			
			m_bCheckIfNotEdgeOfBigStar = ( atol( szTmp )>0 );
		}

		szTmp = GetParam("CCD_EDGE_OF_BIG_BY_RAW_DATA",TRUE);
		if(szTmp && szTmp[0]){
			m_bEdgeOfBigByRawData = ( atol( szTmp )>0 );
		}

		szTmp = GetParam("CCD_SIGMA_ABOVE_MEAN_IN_RAW_CLUSTER",TRUE);
		if(szTmp && szTmp[0]){
			m_nSigmaAboveMeanInRawCluster = atof( szTmp );
		}
	
		szTmp = GetParam("CCD_PREV_FRAMES_CHECK_EDGE",TRUE);
		if(szTmp && szTmp[0]){
			m_nPrevFramesToCheckEdge = atol( szTmp );
		}

		//--- FITS header usful values :
		szTmp = GetParam("CCD_DATE_OBS",TRUE);
		if(szTmp){
			m_DateObs = szTmp;
		}
		szTmp = GetParam("CCD_TIME_OBS",TRUE);
		if(szTmp){
			m_TimeObs = szTmp;
		}
		szTmp = GetParam("CCD_DATE_TIME_OBS_FMT",TRUE);
		if(szTmp && szTmp[0]){
			m_DateTimeObsFormat = szTmp;
		}


		//----------------------------------------------------------------
		szTmp = GetParam("CCD_DO_NOT_ANALYSE",TRUE);
		if(szTmp && szTmp[0]){
			m_bDoNotAnalyse = ( atol( szTmp )>0 );
		}
				

		//----------------------------------------------------------------
		// OUTPUT :
		if(m_nInitCount==0){
			// this parametrs are initialized only once , later calls do not 
			// re-read them :

			szTmp = GetParam("CCD_OUTPUT_DIR",TRUE);			
			if(szTmp && szTmp[0]){
				m_szBaseOutDir = szTmp;
			}else{
				const char* envoutdir = getenv("OUTPUT_DIR");
				if(envoutdir && strlen(envoutdir))
					m_szBaseOutDir = envoutdir;
			}

			szTmp = GetParam("CCD_OUTPUT_SUBDIR",TRUE);
			if(szTmp && szTmp[0]){
				m_szBaseOutSubDir = szTmp;
			}else{
				const char* envoutdir = getenv("OUTPUT_SUBDIR");
				if(envoutdir && strlen(envoutdir))
					m_szBaseOutSubDir = envoutdir;
			}

			m_szOutDir="";
	
			m_szOutDir << m_szBaseOutDir << "/";
				
			if(strlen(m_szBaseOutSubDir.c_str()))
				 m_szOutDir << m_szBaseOutSubDir;
		}

		szTmp = GetParam("CCD_USE_FRAME_TIME_FOR_SHIFTS",TRUE);		
		if(szTmp && szTmp[0]){
			m_bUseFrameTimeForShifts = (atol(szTmp)>0);
		}

		szTmp = GetParam("CCD_IGNORE_MISSING_TIME",TRUE);
		if(szTmp && szTmp[0]){
			m_bIgnoreMissingTime = (atol(szTmp)>0);
		}

		szTmp = GetParam("CCD_SHIFT_USES_ASTRO_FORMULAS",TRUE);
		if(szTmp && szTmp[0]){
			m_bShiftUsesAstroFormulas = (atol(szTmp)>0);
		}
		szTmp = GetParam("CCD_TRANSFORM_ORIENTATION",TRUE);
		if(szTmp && szTmp[0]){
			m_TransformCCDOrientation = atof(szTmp);
		}

		szTmp = GetParam("CCD_AXIS_X_ORIENTATION",TRUE);
	   if(szTmp && szTmp[0]){
   		m_CCDAxixXOrientation = atol(szTmp);
	   }

   	szTmp = GetParam("CCD_AXIS_Y_ORIENTATION",TRUE);
	   if(szTmp && szTmp[0]){
   		m_CCDAxixYOrientation = atol(szTmp);
	   }
			
		szTmp = GetParam("CCD_TRANSFORM_CCD_FOCUS",TRUE);
		if(szTmp && szTmp[0]){
			m_TransformCCDFocus = atof( szTmp );
		}

		szTmp = GetParam("CCD_TRANSFORM_CCD_PIXEL_SIZE",TRUE);
		if(szTmp && szTmp[0]){
			m_TransformCCDPixelSize = atof( szTmp );
		}

		szTmp = GetParam("CCD_ORIGIN",TRUE);
		if(szTmp && szTmp[0]){
			m_szOrigin = szTmp;
		}
		szTmp = GetParam("CCD_SITE",TRUE);
		if(szTmp && szTmp[0]){
		   m_szSite = szTmp;
		}
		szTmp = GetParam("CCD_INSTRUME",TRUE);
		if(szTmp && szTmp[0]){
	   	m_szInstrume = szTmp;
		}
		szTmp = GetParam("CCD_CAMOPTIC",TRUE);
		if(szTmp && szTmp[0]){
		   m_szCamOptic = szTmp;
		}
		szTmp = GetParam("CCD_FILTER",TRUE);
		if(szTmp && szTmp[0]){
	   	m_szFilter = szTmp;
		}
		szTmp = GetParam("CCD_OBSERVER",TRUE);
		if(szTmp && szTmp[0]){
			m_szObserver = szTmp;
		}
		szTmp = GetParam("CCD_OBJECT",TRUE);
		if(szTmp && szTmp[0]){
			m_szObject = szTmp;
		}
		szTmp = GetParam("CCD_FOV",TRUE);
		if(szTmp && szTmp[0]){
			m_FOV = atof( szTmp );
			m_fFrameSize = ( m_FOV/2.0 )*CMyMathFunc::mysqrt((double)2.00) + 3.00;
		}
		szTmp = GetParam("CCD_FOV_BORDER",TRUE);
		if(szTmp && szTmp[0]){
			m_FOVBorder = atof( szTmp );
		}




		// DRIVER parameters :
		szTmp = GetParam("CCD_KERNEL_MODULES_PATH",TRUE);
		if(szTmp && szTmp[0]){
			m_szKernelModulesPath = szTmp;
		}
		m_szKernelModulesPath.env2str();

		szTmp = GetParam("CCD_RELOAD_MODULE_ON_STARTUP",TRUE);
		if(szTmp && szTmp[0]){
			m_bReloadModuleOnStartup = ( atol( szTmp )>0 );
		}

		szTmp = GetParam("CCD_DRIVER_ASYNCHRO_MODE",TRUE);
		if(szTmp && szTmp[0]){
			m_bAsynchroMode = ( atol( szTmp )>0 );
		}

		szTmp = GetParam("CCD_TAKE_IN_SYNCHRO_MODE",TRUE);
		if(szTmp && szTmp[0]){
			m_bTakeInSynchroMode = ( atol( szTmp )>0 );
		}

		szTmp = GetParam("CCD_PARALEL_MODE",TRUE);
		if(szTmp && szTmp[0]){
			m_bParallelMode = ( atol( szTmp )>0 );			
		}

		szTmp = GetParam("CCD_FIX_DEVICE_ORDER",TRUE);
		if(szTmp && szTmp[0]){
			m_bFixDeviceOrder = ( atol( szTmp )>0 );
		}

		szTmp = GetParam("CCD_DEVICE_ID",TRUE);
		if(szTmp && szTmp[0]){
			m_szDeviceID = szTmp;
		}

		szTmp = GetParam("CCD_DRIVER_MAX_ITER_TIMEOUT",TRUE);
		if(szTmp && szTmp[0]){
			m_DriverMaxIterTimeout = atol(szTmp);
		}

		szTmp = GetParam("CCD_RETRY_FRAME_COUNT",TRUE);
		if(szTmp && szTmp[0]){
			m_RetryFrameCount = atol(szTmp);
		}

		szTmp = GetParam("CCD_GETDATA_RETRY_COUNT",TRUE);
		if(szTmp && szTmp[0]){
			m_DriverGetDataRetryCount = atol(szTmp);
		}
	
		szTmp = GetParam("CCD_MAX_ALLOWED_FRAMES_LOST",TRUE);
		if(szTmp && szTmp[0]){
			m_MaxAllowedFramesLost = atol( szTmp );
		}

		szTmp = GetParam("CCD_FORCE_REGISTERS_USAGE",TRUE);
		if(szTmp && szTmp[0]){
			m_bForceRegistersUsage = ( atol(szTmp) > 0 );
		}

		szTmp = GetParam("CCD_ADC_CONF_REG0",TRUE);
		if(szTmp && szTmp[0]){
			m_ADCConfReg0 = atol(szTmp);
		}
		szTmp = GetParam("CCD_ADC_CONF_REG1",TRUE);
		if(szTmp && szTmp[0]){
			m_ADCConfReg1 = atol(szTmp);
		}

		szTmp = GetParam("CCD_ADC_CONF_REG2",TRUE);
		if(szTmp && szTmp[0]){
			m_ADCConfReg2 = atol(szTmp);
		}
		szTmp = GetParam("CCD_ADC_CONF_REG3",TRUE);
		if(szTmp && szTmp[0]){
			m_ADCConfReg3 = atol(szTmp);
		}
		szTmp = GetParam("CCD_ADC_CONF_REG4",TRUE);
		if(szTmp && szTmp[0]){
			m_ADCConfReg4 = atol(szTmp);
		}
		szTmp = GetParam("CCD_ADC_CONF_REG5",TRUE);
		if(szTmp && szTmp[0]){
			m_ADCConfReg5 = atol(szTmp);
		}
		szTmp = GetParam("CCD_ADC_CONF_REG6",TRUE);
		if(szTmp && szTmp[0]){
			m_ADCConfReg6 = atol(szTmp);
		}
		szTmp = GetParam("CCD_ADC_CONF_REG7",TRUE);
		if(szTmp && szTmp[0]){
			m_ADCConfReg7 = atol(szTmp);
		}














		szTmp = GetParam("CCD_ADC_OFFSET",TRUE);
		if(szTmp && szTmp[0]){		
			m_ADCOffset = atol(szTmp);
			m_ADCOffsetCH2 = m_ADCOffset;
		}
		szTmp = GetParam("CCD_ADC_OFFSET_CH2",TRUE);
		if(szTmp && szTmp[0]){
			m_ADCOffsetCH2 = atol(szTmp);
		}

		szTmp = GetParam("CCD_ADC_GAIN",TRUE);
		if(szTmp && szTmp[0]){
			m_ADCGain = (int)(atof(szTmp));
			m_ADCGainCH2 = m_ADCGain;
		}

		szTmp = GetParam("CCD_ADC_GAIN_CH2",TRUE);
		if(szTmp && szTmp[0]){
			m_ADCGainCH2 = (int)(atof(szTmp));
		}


		szTmp = GetParam("CCD_ADC_CLAMPING",TRUE);
		if(szTmp && szTmp[0]){
			m_ADCClampingBitOnOff = ( atol( szTmp )>0 );
		}
		szTmp = GetParam("CCD_ADC_RANGE",TRUE);
		if(szTmp && szTmp[0]){
		   m_eADCRange = (eADCRange)atol( szTmp ); 
		}

		szTmp = GetParam("CCD_LNA_GAIN",TRUE);
		if(szTmp && szTmp[0]){
			m_eLNAGain = (eLNAGainValue_T)atol(szTmp);
		}

		szTmp = GetParam("CCD_LENS_HIT_ON",TRUE);
		if(szTmp && szTmp[0]){
	      m_bLensHitOnOff = (atol(szTmp)>0);
		}

		szTmp = GetParam("CCD_MAX_COMM_ERROR_COUNT_TO_EXIT",TRUE);
      if(szTmp && szTmp[0]){
			m_MaxCommErrorCountToEXIT = atol( szTmp );
		}

		szTmp = GetParam("CCD_EXIT_ON_COMM_ERROR",TRUE);
		if(szTmp && szTmp[0]){
			m_ExitOnCommError = atol(szTmp);
		}

		szTmp = GetParam("CCD_SENDBYTES_RETRY_COUNT",TRUE);
		if(szTmp && szTmp[0]){
			m_SendBytesRetryCount = atol( szTmp );
		}


		szTmp = GetParam("CCD_MPP_MODE",TRUE);
		if( szTmp && szTmp[0]){
			m_eMPPMode = (eMPP_T)(atol(szTmp));
		}
		szTmp = GetParam("CCD_READOUT_SPEED_HORIZONTAL",TRUE);
		if(szTmp && szTmp[0]){
			m_ReadoutSpeedHorizontal = atol(szTmp);
		}
		szTmp = GetParam("CCD_READOUT_SPEED_VERTICAL",TRUE);
		if(szTmp && szTmp[0]){
			m_ReadoutSpeedVertical = atol(szTmp);
		}

		szTmp = GetParam("CCD_TEMPERATURE",TRUE);	
		if(szTmp && szTmp[0]){	
			m_CCDTemp = atof( szTmp );
		}

		szTmp = GetParam("CCD_COOLING_ON_OFF",TRUE);
		if(szTmp && szTmp[0]){
			m_bCoolingOnOff = ( atol( szTmp )>0 );
		}

		szTmp = GetParam("CCD_WAIT_FOR_DESIRED_TEMP",TRUE);
		if(szTmp && szTmp[0]){
			m_bWaitForDesiredTemp = ( atol( szTmp )>0 );
		}

		szTmp = GetParam("CCD_ACTION_ON_TEMP_WAIT_FAIL",TRUE);
		if(szTmp && szTmp[0]){
			m_eActionOnTempWaitFail = (eActionOnTempWaitFail)atol( szTmp );
		}

		szTmp = GetParam("CCD_COOLING_TOLERANCE_DIFF",TRUE);
		if(szTmp && szTmp[0]){
			m_CoolingToleranceDiff = atol( szTmp );
		}

		szTmp = GetParam("CCD_WAIT_FOR_DESIRED_TEMP_IN_SEC",TRUE);
		if(szTmp && szTmp[0]){
			m_WaitForDesiredTempInSec = atol( szTmp );
		}

		szTmp = GetParam("CCD_CAM_TYPE",TRUE);
		if(szTmp && szTmp[0]){
			m_eCAMType = (eCCDTYPE_T)atol(szTmp);
		}

		szTmp = GetParam("CCD_CAM_IDENT_NO",TRUE);
		if(szTmp && szTmp[0]){
			m_CCDIdentNo = szTmp;
		}

		// szTmp = GetParam("CCD_CAM_NO",TRUE); m_CCDNo
		szTmp = GetParam("CCD_SEARCH_DEVICE_START_NO",TRUE);
		if(szTmp && szTmp[0]){
			m_SearchDeviceStartNo = atol(szTmp);
		}		

		szTmp = GetParam("CCD_CAM_NAME",TRUE);
		if(szTmp && szTmp[0]){
			m_CCDName = szTmp;
		}

		szTmp = GetParam("CCD_INTERVAL_BETWEEN_IMAGES",TRUE);
		if(szTmp && szTmp[0]){
			m_IntervalBetweenFrames = atol( szTmp );
		}

		szTmp = GetParam("CCD_BASE_FILE_NAME_FITS",TRUE);
		if(szTmp && szTmp[0]){
			m_szBaseFileNameFITS = szTmp;
		}

		szTmp = GetParam("CCD_FRAMES_OUTPUT_SUBDIR",TRUE);
		if(szTmp && szTmp[0]){
			m_szFramesOutputSubDir = szTmp;
		}

		szTmp = GetParam("CCD_DRIVER_WRITE_ALL_TO_FITS",TRUE);
		if(szTmp && szTmp[0]){
			m_bDriverWriteAllToFITS = ( atol(szTmp)>0 );
		}

		szTmp = GetParam("CCD_WRITE_FRAMES_LIST",TRUE);
		if(szTmp && szTmp[0]){
			m_bBuildFramesList = ( atol(szTmp)>0 );
		}

		szTmp = GetParam("CCD_READ_FULL_STATUS_FREQ",TRUE);
		if(szTmp && szTmp[0]){
			m_ReadFullStatusFreq = atol(szTmp);
		}
		
		szTmp = GetParam("CCD_FITS_COMPRESS_TYPE",TRUE);
		if(szTmp && szTmp[0]){		
			m_eCompressFITS = (eFITSCompressionType)( atol(szTmp) );
		}
		szTmp = GetParam("CCD_DRIVER_SHUTTER_TIME_IN_SEC",TRUE);
		if(szTmp && szTmp[0]){
			m_DriverShutterTimeInSec = atof( szTmp );
			m_DriverShutterTimeInSecSaved = m_DriverShutterTimeInSec;
		}

		// shutter mode with save field : m_ShutterModeOriginal
		szTmp = GetParam("CCD_SHUTTER_MODE",TRUE);
		if(szTmp && szTmp[0]){
			m_ShutterMode = atol( szTmp );
		}
		m_ShutterModeOriginal = m_ShutterMode;

		// SHUTTER BREAKING VOLTAGES TIMES :
		szTmp = GetParam("CCD_CAMERA_HAS_BRAKE_VOLTAGES",TRUE);
		if(szTmp && szTmp[0]){
			m_bHasBreakVoltages = ( atol( szTmp ) > 0 );
		}
	

		szTmp = GetParam("CCD_DRIVER_REVERSE_IMAGE",TRUE);
		if(szTmp && szTmp[0]){
			m_bDriverReverseImage = (eDriverReverseImage_T)(atol(szTmp));
		}
		szTmp = GetParam("CCD_USE_CAM_ID",TRUE);
		if(szTmp && szTmp[0]){
			m_bUseCamID = ( atol(szTmp)>0 );
		}

		szTmp = GetParam("CCD_OPEN_SHUTTER_WHEN_STARCOUNT_BELOW",TRUE);
		if(szTmp && szTmp[0]){
			m_MinStarCountToCloseShutter = atol( szTmp );
		}

// gigabit camera :
		szTmp = GetParam("CCD_CAMERA_IP",TRUE);
		if(szTmp && szTmp[0]){
			m_szCameraIP = szTmp;
		}

		szTmp = GetParam("CCD_CAMERA_PORT_NO",TRUE);
		if(szTmp && szTmp[0]){
			m_CameraPortNo = atol(szTmp);
		}

		szTmp = GetParam("CCD_LOCAL_PORT_NO",TRUE);
		if(szTmp && szTmp[0]){
			m_LocalPortNo = atol(szTmp);
		}

		szTmp = GetParam("CCD_CAMERA_TIMEOUT_MILI_SEC", TRUE);
		if(szTmp && szTmp[0]){
			m_CameraTimeoutMiliSec = atol(szTmp);
		}
	
		szTmp = GetParam("CCD_ETHCAM_DATA_RETRIES", TRUE );
		if(szTmp && szTmp[0]){
			m_EthCamDataRetries = atol(szTmp);
		}
		szTmp = GetParam("CCD_ETHCAM_COMMAND_TIMEOUT", TRUE );
		if(szTmp && szTmp[0]){
		   m_EthCamCmdTimeout = atol(szTmp);
		}
		szTmp = GetParam("CCD_ETHCAM_COMMAND_RETRIES", TRUE );
		if(szTmp && szTmp[0]){
		   m_EthCamCmdRetries = atol(szTmp);
		}

		szTmp = GetParam("CCD_ETH_CAM_LOGFILE",TRUE);
		if(szTmp){
			m_szEthCamLogFile = szTmp;
		}

		szTmp = GetParam("CCD_ETH_CAM_LOG_LEVEL", TRUE);
		if(szTmp && szTmp[0]){
			m_nEthCamLogLevel = atol( szTmp );
		}

		szTmp = GetParam("CCD_ETH_CAM_ERRLOG_LEVEL", TRUE);
		if(szTmp && szTmp[0]){
			m_nEthCamErrLogLevel = atol( szTmp );
		}

		szTmp = GetParam("CCD_ETH_CAM_ERRFILE",TRUE);
      if(szTmp){
			m_szEthCamErrFile = szTmp;
		}

		szTmp = GetParam("CCD_DAQ_LOG_FILE",TRUE);
		if(szTmp && szTmp[0]){
			m_szDAQLogFile = szTmp;
		}

		szTmp = GetParam("CCD_GLOBAL_LOG_ENABLED",TRUE);
		if(szTmp && szTmp[0]){
			m_bPiLogEnabled = ( atol(szTmp) > 0 );
/*	temporary changed not to mix messages from other programs then ccdsingle
   and ccddouble in same log file, there should be different one : 
		pi_offline.log for example
			if( m_bPiLogEnabled ){
				InitPiLog( m_bPiLogEnabled );
			}*/
		}

// simulator parameters :
		szTmp = GetParam("CCD_REPEAT_SAME_IMAGES", TRUE);
		if(szTmp && szTmp[0]){
			m_bRepeatSameImages = ( atol(szTmp) > 0 );
		}

		szTmp = GetParam("CCD_DAQ_COMMUNICATION_ON",TRUE);
		if(szTmp && szTmp[0]){
			m_bDAQCommunicationON = ( atol(szTmp)>0 );
		}

		szTmp = GetParam("CCD_PI_SYS_MANAGER_ON",TRUE);
		if(szTmp && szTmp[0]){
			m_bPISysManagerON = ( atol(szTmp)>0 );
		}
		
		szTmp = GetParam("CCD_PI_SYS_MANAGER_PORT_NO",TRUE);
		if(szTmp && szTmp[0]){
			m_PISysManPortNo = atol( szTmp );
		}
		
		szTmp = GetParam("CCD_ASK_MOUNT_FOR_COORD",TRUE);
		if(szTmp && szTmp[0]){
			m_bAskMountForCoord = ( atol(szTmp)>0 );
		}

		// MOUNT INFORMATION :
		szTmp = GetParam("CCD_MOUNT_ID",TRUE);
		if(szTmp && szTmp[0]){
		   m_MountID = ( atol(szTmp)>0 );
		}

		szTmp = GetParam("CCD_DRIVER_BINNING",TRUE);
		if(szTmp && szTmp[0]){
			m_DriverAnalBinning = atol( szTmp );
		}
		szTmp = GetParam("CCD_PORT_NO",TRUE);
		if(szTmp && szTmp[0]){
			m_PortNo = atol( szTmp );
		}

		szTmp = GetParam("CCD_PIPELINE_RESTART_TIMEOUT",TRUE);
		if(szTmp && szTmp[0]){
			m_nPipelineRestartTimeout = atol( szTmp );
		}

		szTmp = GetParam("CCD_COORD_CHANGE_TO_RESTART",TRUE);
      if(szTmp && szTmp[0]){
			m_CoordChangeToRestart = atof( szTmp );
		}

		szTmp = GetParam("CCD_CAMERA_INTERFACE_NO",TRUE);
		if(szTmp && szTmp[0]){
			m_CameraInterfaceNo = atol( szTmp );
		}

		szTmp = GetParam("CCD_CORBA_OPTIONS",TRUE);
		if(szTmp && szTmp[0]){
			m_CorbaOptions = szTmp;
		}
		szTmp = GetParam("CCD_CORBA_WITH_NAME_SERVICE",TRUE);
		if(szTmp && szTmp[0]){
			m_bCorbaWithNameService = ( atol(szTmp)>0 );
		}

		szTmp = GetParam("CCD_SHARED_MEM_KEY",TRUE);
		if(szTmp && szTmp[0]){
			m_SharedMemKey = atol( szTmp );
		}

		szTmp = GetParam("CCD_PIPELINE_SAFE_BUFFER_SIZE",TRUE);
		if(szTmp && szTmp[0]){
			m_PipelineSafeBufferSize = atol( szTmp );
		}

		// automatic determination of shifts :
		szTmp = GetParam("CCD_SHIFTS_VALUES_FILE",TRUE);
		if(szTmp && szTmp[0]){
			m_szShiftsValuesFile = szTmp;
		}
		
		szTmp = GetParam("CCD_TRACE_ON_ALL_FRAMES",TRUE);
		if(szTmp && szTmp[0]){
			m_bTraceOnAllFrames = ( atol(szTmp)>0 );
		}

		szTmp = GetParam("CCD_USE_TRANSFORM_MATRIX",TRUE);
		if(szTmp && szTmp[0]){
			m_bUseTransformMatrix = ( atol(szTmp)>0 );
		}

		szTmp = GetParam("CCD_MAX_DIFF_SIGNAL_IN_SIGMA",TRUE);
		if(szTmp && szTmp[0]){
			m_MaxDiffSignalInSigma = atof(szTmp);
		}

		szTmp = GetParam("CCD_DUMP_TRACE_STAR_TO_FILE",TRUE);
		if(szTmp && szTmp[0]){
			m_bDumpTraceStarToFile = ( atol(szTmp)>0 );
		}

		szTmp = GetParam("CCD_USE_CONTROL_STAR",TRUE);
		if(szTmp && szTmp[0]){
			m_bUseControllStar = ( atol(szTmp)>0 );
		}

		szTmp = GetParam("CCD_MAX_ALLOWED_DRASTIC_CHNG",TRUE);
		if(szTmp && szTmp[0]){
			m_nMaxAllowedDrasticChng = atol(szTmp);
		}

		szTmp = GetParam("CCD_AUTO_SHIFTS_MIN_STEPS_TO_IGNORE",TRUE);
		if(szTmp && szTmp[0]){
			m_nAutoShiftsMinStepsToIgnore = atol(szTmp);
		}

		szTmp = GetParam("CCD_AUTO_SHIFTS_CALC",TRUE);
		if(szTmp && szTmp[0]){
			m_nAutoShiftsCalc = atol(szTmp);
		}
		szTmp = GetParam("CCD_AUTO_SHIFTS_CALC_MATRIX_NO",TRUE);
		if(szTmp && szTmp[0]){
			m_nAutoShiftMatrixNo = atol(szTmp);
		}

		szTmp = GetParam("CCD_MAX_STAR_SHAPE_LIMIT",TRUE);
		if(szTmp && szTmp[0]){
			m_MAXStarShapeLimit = atof( szTmp );
		}

		szTmp = GetParam("CCD_MAX_STAR_N_SIGMA_ABOVE",TRUE);
		if(szTmp && szTmp[0]){
			m_MAXStarNSigmaAbove = atof( szTmp );
		}

		szTmp = GetParam("CCD_TRANSFORM_MATRIX_FILE",TRUE);
		if(szTmp && szTmp[0]){
			m_szTransformFile = szTmp;
			if(!m_TransformMatrix.ReadFromFile( m_szTransformFile.c_str() )){
				Assert(FALSE,"Could not read transformation matrix from file : %s\n",m_szTransformFile.c_str());
			}
			if( gPrintfLevel>=2 ){
				printf("Transformation matrix :\n");
				m_TransformMatrix.Dump();
			}
		}


		szTmp = GetParam("CCD_AUTO_CHECK_SIZE",TRUE);
		if( szTmp && szTmp[0] ){
			m_bAutoCheckSize = ( atol(szTmp)>0 );
		}

		// trigger actions :

		szTmp = GetParam("CCD_SKIP_ASTRO_ON_ALERT",TRUE);
		if( szTmp && szTmp[0] ){
			m_bSkipAstroOnAlert = ( atol(szTmp)>0 );
		}

		szTmp = GetParam("CCD_POSTPONE_ASTROMETRY_TIME",TRUE);
		if( szTmp && szTmp[0] ){
			m_PostponeAstrometryTime = atol(szTmp);
		}

		szTmp = GetParam("CCD_HANDLE_GCN",TRUE);
		if( szTmp && szTmp[0] ){
			m_bHandleGCN = ( atol(szTmp)>0 );
		}

		szTmp = GetParam("CCD_TRIGGER_ACTIONS_DEF_FILE",TRUE);
		if( szTmp && szTmp[0] ){
			m_szTriggerActionsDefFile = szTmp;
		}
	   m_szTriggerActionsDefFile.env2str();
		szTmp = GetParam("CCD_INTERNAL_TRIGGER_ACTIONS_DEF_FILE",TRUE);
		if( szTmp && szTmp[0] ){
			m_szInternalTriggerActionsDefFile = szTmp;
		}
		m_szInternalTriggerActionsDefFile.env2str();

		szTmp = GetParam("CCD_SAVE_FULL_FRAMES_ON_TRIGGER",TRUE);
		if( szTmp && szTmp[0] ){
			m_bSaveFullFramesOnTrigger = atol( szTmp );
		}

		szTmp = GetParam("CCD_TRIGGER_CHECK_ACTION_PERIOD",TRUE);
		if( szTmp && szTmp[0] ){
			m_TriggerCheckActionPeriod = atol(szTmp);
		}


		// Event significance :
		szTmp = GetParam("CCD_SEND_INTERNAL_TRIGGERS",TRUE);
		if( szTmp && szTmp[0] ){
			m_bSendInternalTriggers = ( atol(szTmp)>0 );
		}

		//-------- SATELITES 
		szTmp = GetParam("CCD_CHECK_IF_SATELITE",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckIfSatelite = ( atol( szTmp )>0 );
		}
		szTmp = GetParam("CCD_REJECT_SATELITES",TRUE);
		if(szTmp && szTmp[0]){
		   m_bRejectSatelite = ( atol( szTmp )>0 );
		}
		szTmp = GetParam("CCD_TLE_FILE",TRUE);
		if(szTmp && szTmp[0]){
		   m_szTleFile = szTmp;
			m_szTleFile.env2str();
		}
		szTmp = GetParam("CCD_QTH_FILE",TRUE);
		if(szTmp && szTmp[0]){
	   	m_szQthFile = szTmp;
			m_szQthFile.env2str();
		}
		szTmp = GetParam("CCD_SAT_REJ_RADIUS",TRUE);
		if(szTmp && szTmp[0]){
			m_nSatRejRadius = atof( szTmp );
		}
		szTmp = GetParam("CCD_SAT_REJ_RADIUS_IN_SEC",TRUE);
		if(szTmp && szTmp[0]){
			m_nSatRejRadius = AstroAngle::arcsec2rad( atof( szTmp ) );
		}

		szTmp = GetParam("CCD_NOT_VISIBLE_SAT_REJ_RADIUS_IN_SEC",TRUE);
		if(szTmp && szTmp[0]){
			m_nNotVisibleSatRejRadius = AstroAngle::arcsec2rad( atof( szTmp ) );
		}

		if( m_bCheckIfSatelite ){
			InitSatLib();
		}
		//-------- END SATELITES 

		//--------- STAR REJECTION :
		szTmp = GetParam("CCD_CHECK_IF_STAR",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckIfStar = ( atol( szTmp )>0 );
		}
		szTmp = GetParam("CCD_REJECT_STARS",TRUE);
		if(szTmp && szTmp[0]){
		   m_bRejectStars = ( atol( szTmp )>0 );
		}
		szTmp = GetParam("CCD_STAR_REJ_RADIUS",TRUE);
		if(szTmp && szTmp[0]){
	   	m_fStarRejectRadius = atof( szTmp );
		}
		szTmp = GetParam("CCD_STAR_REJ_RADIUS_IN_SEC",TRUE);
		if(szTmp && szTmp[0]){
	   	m_fStarRejectRadius = AstroAngle::arcsec2rad( atof(szTmp) );
		}
		szTmp = GetParam("CCD_COUNT_BRIGHTER_THEN",TRUE);
		if(szTmp && szTmp[0]){
			m_fCountBrighterThen = atof( szTmp );
		}
		szTmp = GetParam("CCD_REJECT_IF_BIG_STAR_NEAR_BY",TRUE);
		if(szTmp && szTmp[0]){
			m_bRejectIfBigStarNearBy = ( atol( szTmp )>0 );
		}
		szTmp = GetParam("CCD_BIG_STAR_REJ_RADIUS_IN_ARCSEC",TRUE);
		if(szTmp && szTmp[0]){
		   m_fBigStarRejectRadiusInArcSec = atof( szTmp );
		}
		szTmp = GetParam("CCD_BIG_STAR_REJ_RADIUS_IN_PIXEL",TRUE);
		if(szTmp && szTmp[0]){
		   m_fBigStarRejectRadiusInArcSec = atof( szTmp )*m_fPixScale;
		}
		szTmp = GetParam("CCD_BIG_STAR_MAX_MAGNITUDO",TRUE);		
		if(szTmp && szTmp[0]){
		   m_fBigStarMaxMagnitudo = atof( szTmp );
		}

		szTmp = GetParam("CCD_CHECK_STARS_IN_TYCHO",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckStarsInTycho = ( atol(szTmp)>0 );
		}
		szTmp = GetParam("CCD_STARCAT_MAX_MAG",TRUE);
		if(szTmp && szTmp[0]){
			m_fStarCatMaxMag = atof( szTmp );
		}

		//--------- END OF STAR REJECTION

		//--------- PI DB SECTION :
		szTmp = GetParam("CCD_USE_DB",TRUE);
		if(szTmp && szTmp[0]){
			m_bUseDB = ( atol(szTmp)>0 );
		}
		szTmp = GetParam("CCD_USE_ODBC",TRUE);
		if(szTmp && szTmp[0]){
			m_bUseODBC = ( atol(szTmp)>0 );
		}
		szTmp = GetParam("CCD_DBNAME",TRUE);
		if(szTmp && szTmp[0]){
	   	m_szDBName = szTmp;
		}
		szTmp = GetParam("CCD_DBUSER",TRUE);
		if(szTmp && szTmp[0]){
   		m_szDBUser = szTmp;
		}
		szTmp = GetParam("CCD_DBPASS",TRUE);
		if(szTmp && szTmp[0]){
	   	m_szDBPass = szTmp;
		}
		szTmp = GetParam("CCD_DBHOST",TRUE);
		if(szTmp && szTmp[0]){
	   	m_szDBHost = szTmp;
		}
		szTmp = GetParam("CCD_DB2_SCHEMA",TRUE);
		if(szTmp && szTmp[0]){
			m_szDB2Schema = szTmp;
		}

/*		szTmp = GetParam("CCD_DB_LOG_LEVEL",TRUE);
		if(szTmp && szTmp[0]){
			CPiDBInterface::SetLogLevel( atol( szTmp ) );
		}*/

		szTmp = GetParam("CCD_SCAN_DBNAME",TRUE);
		if(szTmp && szTmp[0]){
			m_szScanDB = szTmp;
		}
		szTmp = GetParam("CCD_SAVE_EVENTS_TO_DB",TRUE);
		if(szTmp && szTmp[0]){
			m_bSaveEventsToDB = ( atol(szTmp)>0 );
		}
		szTmp = GetParam("CCD_SAVE_SUM_EVENTS_TO_DB",TRUE);
      if(szTmp && szTmp[0]){
			m_bSaveSumEventsToDB = ( atol(szTmp)>0 );
		}
		szTmp = GetParam("CCD_SAVE_VERIF_EVENTS_TO_DB",TRUE);
		if(szTmp && szTmp[0]){
			m_bSaveVerifEventsToDB = ( atol(szTmp)>0 );
		}
		szTmp = GetParam("CCD_SAVE_ALL_EVENTS_TO_DB",TRUE);
		if(szTmp && szTmp[0]){
			m_bSaveAllEventsToDB = ( atol(szTmp)>0 );
		}

		szTmp = GetParam("CCD_SAVE_FRAMES_TO_DB",TRUE);
		if(szTmp && szTmp[0]){
			m_bSaveFramesToDB = ( atol(szTmp)>0 );
		}
		szTmp = GetParam("CCD_TRACKS_TO_DB",TRUE);
		if(szTmp && szTmp[0]){
			m_bSaveTracksToDB = ( atol(szTmp)>0 );
		}

		szTmp = GetParam("CCD_RUN_TYPE",TRUE);
		if(szTmp && szTmp[0]){
			m_eRunType = atol( szTmp );
		}
		//---------- END OF PI DB SECTION


		//---------- PAST GRB QUERY SECTION :

		szTmp = GetParam("CCD_SN_RADIUS_IN_ARCSEC",TRUE);
      if(szTmp && szTmp[0]){
			m_SNRadiusInArcSec = atof(szTmp);
		}
		szTmp = GetParam("CCD_CHECK_FOR_EXTERNALS_IN_DB",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckForExternalsInDB = ( atol(szTmp)>0 );
		}
		szTmp = GetParam("CCD_GRB_RADIUS_IN_DEG",TRUE);
		if(szTmp && szTmp[0]){
			m_fGrbRadiusInDeg = atof( szTmp );
		}
		szTmp = GetParam("CCD_TIME_PAST_TO_CHECK_GRB_IN_SEC",TRUE);
		if(szTmp && szTmp[0]){
		   m_TimePastToCheckInSec = atof( szTmp );
		}
		//----------- END OF PAST GRB QUERY SECTION 


		// recalculation of parameters :
		if( !bCoicRadiusInRad ){
			// radius in rad was not passed - calculating from 
			// pixel coic radius
			m_nCoicRadiusInRad = GetPixToRad( m_nCoicRedial ); 
		}
	
		szTmp = GetParam("CCD_READ_FIRST_LEVEL_INFO_FROM_LOG", TRUE );
		if(szTmp && szTmp[0]){
			m_bReadFirstLevelInfoFromLog = ( atol(szTmp)>0 );
		}

		szTmp = GetParam("CCD_FROM_LOG_FILE_WITH_FRAMES", TRUE );
      if(szTmp && szTmp[0]){
			m_bFromLogFile_WithFrames = ( atol(szTmp)>0 );
		}

		szTmp = GetParam("CCD_FIRST_LEVEL_TRIGGER_DIR",TRUE);
		if(szTmp && szTmp[0]){
			m_szFirstLevelTriggerDir = szTmp;
		}

		// SECOND LEVEL TRIGGER ( SLT )
		szTmp = GetParam("CCD_SLT_CHECK_CLOUDS",TRUE);
		if(szTmp && szTmp[0]){
			m_CheckClouds = ( atol(szTmp)>0 );
		}
		szTmp = GetParam("CCD_SLT_TYPICAL_STARS_COUNT",TRUE);
		if(szTmp && szTmp[0]){
	   	m_TypicalStarsCount = atol(szTmp);
		}
		szTmp = GetParam("CCD_SLT_REJECT_FRAME_IF_LESS",TRUE);
		if(szTmp && szTmp[0]){
			m_RejectFrameIfLess = atof( szTmp );
		}
		szTmp = GetParam("CCD_HOUGH_DISTR_MAX_LIMIT",TRUE);
		if(szTmp && szTmp[0]){
			m_HoughDistrMaxLimit = atof( szTmp );
		}
		szTmp = GetParam("CCD_HOUGH_TRANSFORM_TRESH",TRUE);
		if(szTmp && szTmp[0]){
			m_HoughTransformTresh = atof( szTmp );
		}
		szTmp = GetParam("CCD_HOUGH_DISTR_TRESH",TRUE);
      if(szTmp && szTmp[0]){
			m_HoughDistrTresh = atof( szTmp );
		}
		szTmp = GetParam("CCD_SMALL_HOUGH_DISTR_TRESH",TRUE);
		if(szTmp && szTmp[0]){
         m_SmallHoughDistrTresh = atof( szTmp );
      }
		szTmp = GetParam("CCD_LAP_DIFF_MIN_RATIO",TRUE);
		if(szTmp && szTmp[0]){
			m_LapDiffMinRatio = atof( szTmp );
		}
		szTmp = GetParam("CCD_CHECK_HOUGH_ON_SMALL",TRUE);		
		if(szTmp && szTmp[0]){
			m_bCheckHoughOnSmall = ( atol(szTmp)>0 );
		}
		szTmp = GetParam("CCD_SMALL_HOUGH_SIZE",TRUE);
		if(szTmp && szTmp[0]){
		   m_nSmallHoughSize = atol(szTmp);
		}
		szTmp = GetParam("CCD_MIN_STAR_COUNT_ON_PART",TRUE);
		if(szTmp && szTmp[0]){
			m_MinStarCountOnParts = atol(szTmp);
		}
		szTmp = GetParam("CCD_CHECK_CLOUDS_ON_PREV",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckCloudsOnPrev = ( atol(szTmp)>0 );
		}
		szTmp = GetParam("CCD_CHECK_HOUGH_ON_RAW",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckHoughOnRaw = ( atol(szTmp)>0 );
		}
		szTmp = GetParam("CCD_HOUGH_TRANSFORM_ON_RAW_TRESH",TRUE);
		if(szTmp && szTmp[0]){
			m_HoughTransformOnRawTresh = atof( szTmp );
		}
		szTmp = GetParam("CCD_HOUGH_DISTR_ON_RAW_TRESH",TRUE);
      if(szTmp && szTmp[0]){
			m_HoughDistrOnRawTresh = atof( szTmp );
		}
		szTmp = GetParam("CCD_MIN_BACK_FRAMES_BELOW_TO_REJECT",TRUE);
		if(szTmp && szTmp[0]){
			m_nMinBelowLimitToReject = atol(szTmp);
		}
		szTmp = GetParam("CCD_CHECK_TRACK_RADEC",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckTrackRADEC = ( atol(szTmp)>0 );
		}

		CheckConsistancy();
	
		if(gCCDTrace.GetTrLevel()>=5)
			Dump();
		m_bInitialized=TRUE;
		m_nInitCount++;


		ModifyParams();

		// printf("End of CCDConfig::InitParams\n");
	}
	return TRUE;
} // end of CCDConfig::InitParams


void CCDConfig::ModifyParams()
{
	if( m_eRunType==eRunOnlineCoic ){
		if( GetMC() ){
			m_eRunType=eOfflineCoic;
			if( m_bDoSumOfPrevNFrames ){
				m_eRunType=eRunOfflineSum;
			}
			if( m_bCheckForSUPERNEW ){
				m_eRunType=eRunOfflineSN;
			}
		}		
	}

}

BOOL_T CCDConfig::InitParam( const char* szParamName )
{
	const char* szTmp = NULL;
	BOOL_T bRetVal;
	double tmp_val;
	int tmp_val_int;

	// run up to :
	if(strcmp( szParamName, "CCD_RUN_UP_TO_DTM" )==0){
		szTmp = GetParam("CCD_RUN_UP_TO_DTM",TRUE);
		if(szTmp && szTmp[0]){
			m_szRunUpTo = szTmp;
			m_RunUpTo = get_runupto_dtm( m_szRunUpTo.c_str() );
		}
		return TRUE;
	}


	if(strcmp( szParamName, "CCD_WAITING_MODE" )==0){
      szTmp = GetParam("CCD_WAITING_MODE",TRUE);
      if(szTmp && szTmp[0]){
         (CCDPipeline::m_WorkingMode).m_bWaitingMode = ( atol( szTmp )>0 );
      }
		return TRUE;
   }

//	if( UpdateParam( szParamName, "CCD_SIMUL_ROTATE_EVERY_N_IMAGES", CDeviceSimulator::m_RotateEveryNImages, bRetVal ) ){
//		return bRetVal;
//	}

	if( UpdateParam( szParamName, "CCD_SAVE_ASTRO_IN_TAKE_N_MODE" , m_bSaveAstroInTakeNMode, bRetVal ) ){
		return bRetVal;
	}


/*	if( UpdateParam( szParamName, "CCD_GLOBAL_LOG_ENABLED", CLog4C::m_bEnabled, bRetVal ) ){
		return bRetVal;
	}*/


	if( UpdateParam( szParamName, "CCD_CAM_IDX_FOR_AG_ON_OFF", m_CamIdxForAG_OnOff, bRetVal ) ){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_HANDLE_GCN", m_bHandleGCN, bRetVal ) ){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_MIN_STAR_COUNT_TO_RUN_ASTRO", m_nMinStarCountToRunAstro, bRetVal ) ){
		return bRetVal;
	}

// real analysis / cameras :
	if( UpdateParam( szParamName, "CCD_IGNORE_CAMERA", m_bIgnoreCamera, bRetVal ) ){
		return bRetVal;
	}	

	if( UpdateParam( szParamName, "CCD_SAVE_PARAM_CHANGES", m_bSaveParamChanges, bRetVal ) ){
		return bRetVal;
	}	

// ERROR handling
	if( UpdateParam( szParamName, "CCD_EXIT_ON_TOO_MANY_ERRORS", m_nExitOnToManyErrors, bRetVal ) ){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_CHECK_PLANE_TRACKS", m_bCheckPlaneTracks, bRetVal ) ){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_CHECK_DOME_STATUS_FREQ_IN_SEC", m_CheckDomeStatusFreqInSec , bRetVal ) ){
		return bRetVal;
	}
	if( UpdateParam( szParamName, "CCD_STOP_ON_DOME_CLOSE", m_bStopOnDomeClose, bRetVal ) ){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_DO_FLAT_MODE_MODULO", m_nDoFlatModeModulo, bRetVal ) ){
		return bRetVal;	
	}

	if( UpdateParam( szParamName, "CCD_USE_DB", m_bUseDB, bRetVal ) ){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_WAIT_FOR_FRAME", m_bWaitForFrame, bRetVal ) ){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_SAVE_AVERAGE_PARTS", m_bSaveAverageParts , bRetVal ) ){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_ASAS_ASTROMETRY_TIMEOUT", gMaxTimeForAstrometryInSec, bRetVal ) ){
		return bRetVal;
	}
	
	if( UpdateParam( szParamName, "CCD_ASTRO_TO_SILENT_AFTER_N_GOOD", m_nChangeToSilentAfterNGood, bRetVal ) ){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_KEEP_MAG_AND_AST_FROM_ASTRO", m_bKeepMagAndAstFromAstro, bRetVal ) ){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_ASTROMETRY_SYNCHRO_MODE", m_bExecAstroSynchroMode , bRetVal ) ){
		return bRetVal;
	}	

	if( UpdateParam( szParamName, "CCD_PARALEL_MODE", m_bParallelMode , bRetVal ) ){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_MIN_LAPLACE_ON_OTHER_ENABLED",
					  m_bMinLaplaceOnOtherEnabled, bRetVal ) ){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_MIN_DIST_OF_EVENTS_IN_TRACK", m_nMinDistOfEventsInTrack, bRetVal ) ){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_CHECK_VELOCITY_FOR_TRACK", m_bCheckVelocity, bRetVal )){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_CHECK_VELOCITY_ERROR_FOR_TRACK", m_fVelocityError, bRetVal )){
		return bRetVal;
	}	


	if(strcmp( szParamName, "CCD_LOG_TRACKS" )==0){
		szTmp = GetParam("CCD_LOG_TRACKS",TRUE);
		if(szTmp && szTmp[0]){
			m_bLogTracks = ( atol( szTmp )>0 );
		}
		return TRUE;
	}

	if(strcmp( szParamName, "CCD_LOG_FRAME_STAT" )==0){
		szTmp = GetParam("CCD_LOG_FRAME_STAT",TRUE);
		if(szTmp && szTmp[0]){
			m_bLogFrameStat = ( atol( szTmp )>0 );
		}
		return TRUE;
	}
	
	if(strcmp( szParamName, "CCD_TRACE_LEVEL" )==0){
		szTmp = GetParam("CCD_TRACE_LEVEL",TRUE);
      if(szTmp && szTmp[0]){
         m_CCDTraceLevel = atol( szTmp );
         gCCDTrace.SetTrLevel( m_CCDTraceLevel );
      }
		return TRUE;
	}

	if(strcmp( szParamName, "CCD_DO_NOT_USE_ROT" )==0){
		szTmp = GetParam("CCD_DO_NOT_USE_ROT",TRUE);
      if(szTmp && szTmp[0]){
         m_bDoNotUseRotation = ( atol( szTmp )>0 );
      }
		return TRUE;
	}
	
	if(strcmp( szParamName, "CCD_SINGLE_FRAME_DX" )==0){
		szTmp = GetParam("CCD_SINGLE_FRAME_DX", TRUE );
		if(szTmp && szTmp[0]){			
			m_FrameDX = atof( szTmp );			
			return TRUE;
		}
	}

	if(strcmp( szParamName, "CCD_SINGLE_FRAME_DY" )==0){
		szTmp = GetParam("CCD_SINGLE_FRAME_DY", TRUE );
		if(szTmp && szTmp[0]){			
			m_FrameDY = atof( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName, "CCD_SINGLE_FRAME_DX_PER_SEC")==0){			
		szTmp = GetParam("CCD_SINGLE_FRAME_DX_PER_SEC",TRUE);
		if(szTmp && szTmp[0]){
			m_FrameDXPerSec = atof( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName, "CCD_SINGLE_FRAME_DY_PER_SEC")==0){
		szTmp = GetParam("CCD_SINGLE_FRAME_DY_PER_SEC",TRUE);
		if(szTmp && szTmp[0]){
			m_FrameDYPerSec = atof( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName, "CCD_USE_SHIFT_TOTAL" )==0){
			szTmp = GetParam("CCD_USE_SHIFT_TOTAL",TRUE);
			if(szTmp && szTmp[0]){
				m_bUseShiftTotal = (atol(szTmp)>0);
			return TRUE;
		}
	}


	if(strcmp( szParamName, "CCD_ROT_CENTER_X" )==0){
			szTmp = GetParam("CCD_ROT_CENTER_X",TRUE);
			if(szTmp && szTmp[0]){	
				m_RotCenterX = atof( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName, "CCD_ROT_CENTER_Y" )==0){
			szTmp = GetParam("CCD_ROT_CENTER_Y",TRUE);
			if(szTmp && szTmp[0]){
				m_RotCenterY = atof( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName, "CCD_SINGLE_FRAME_D_ALFA" )==0){
			szTmp = GetParam("CCD_SINGLE_FRAME_D_ALFA",TRUE);
			if(szTmp && szTmp[0]){
				m_RotValueDAlfa = atof( szTmp );
			return TRUE;
		}
	}


	if(strcmp( szParamName, "CCD_SINGLE_FRAME_D_ALFA_PER_SEC" )==0){
			szTmp = GetParam("CCD_SINGLE_FRAME_D_ALFA_PER_SEC",TRUE);
			if(szTmp && szTmp[0]){
				m_RotValueDAlfaPerSec = atof( szTmp );
			return TRUE;
		}
	}


	if(strcmp( szParamName, "CCD_HOMEO_USE_SAME_TRESH_FOR_NEW_AND_PREV")==0){
		szTmp = GetParam("CCD_HOMEO_USE_SAME_TRESH_FOR_NEW_AND_PREV",TRUE);
		if(szTmp && szTmp[0]){
			m_bUseHomeoSameTreshForNewAndPrev = ( atol(szTmp)>0 );
			return TRUE;
		}
	}

	if(strcmp( szParamName, "CCD_CHECK_SUM_AROUND" )==0){
		// always read :
		szTmp = GetParam("CCD_CHECK_SUM_AROUND",TRUE);
		if(szTmp && szTmp[0]){
			m_bAnalyseSumAround = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName, "CCD_CHECK_MAX_FROM_PREV" )==0){
		szTmp = GetParam("CCD_CHECK_MAX_FROM_PREV",TRUE);
		if(szTmp && szTmp[0]){
			m_bAnalyseMaxFromPrevFrames = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName, "CCD_CHECK_FOR_FLASHES" )==0){
		szTmp = GetParam("CCD_CHECK_FOR_FLASHES",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckForFlashes = ( atol(szTmp)>0 );
			return TRUE;
		}
	}

	if(strcmp( szParamName, "CCD_BACK_FRAMES_USED" )==0){
		szTmp = GetParam("CCD_BACK_FRAMES_USED",TRUE);
		if(szTmp && szTmp[0]){
			m_FramesBack = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName, "CCD_MAX_VALUE")==0){
		szTmp = GetParam("CCD_MAX_VALUE",TRUE);
		if(szTmp && szTmp[0]){
			m_MaxAllowedVal = atol( szTmp );
			return TRUE;
		}
	}


	if(strcmp( szParamName, "CCD_IGNORE_EDGE" )==0){
		szTmp = GetParam("CCD_IGNORE_EDGE", TRUE );
		if(szTmp && szTmp[0]){
			m_nIgnoreEdge = atol( szTmp );
			m_nIgnoreEdgeRight = m_nIgnoreEdge;
			m_nIgnoreEdgeLeft = m_nIgnoreEdge;
			m_nIgnoreEdgeUp = m_nIgnoreEdge;
			m_nIgnoreEdgeBottom = m_nIgnoreEdge;
			return TRUE;
		}
	}

	if(strcmp( szParamName, "CCD_IGNORE_EDGE_RIGHT" )==0){		
		szTmp = GetParam("CCD_IGNORE_EDGE_RIGHT",TRUE);
		if(szTmp && szTmp[0]){
			m_nIgnoreEdgeRight = atol(szTmp);
			return TRUE;
		}
	}
	if(strcmp( szParamName, "CCD_IGNORE_EDGE_LEFT" )==0){
		szTmp = GetParam("CCD_IGNORE_EDGE_LEFT",TRUE);
		if(szTmp && szTmp[0]){
			m_nIgnoreEdgeLeft = atol(szTmp);
			return TRUE;
		}
	}
	if(strcmp( szParamName, "CCD_IGNORE_EDGE_UP" )==0){
		szTmp = GetParam("CCD_IGNORE_EDGE_UP",TRUE);
		if(szTmp && szTmp[0]){
			m_nIgnoreEdgeUp = atol(szTmp);
			return TRUE;
		}
	}
	if(strcmp( szParamName, "CCD_IGNORE_EDGE_BOTTOM" )==0){
		szTmp = GetParam("CCD_IGNORE_EDGE_BOTTOM",TRUE);
		if(szTmp && szTmp[0]){
			m_nIgnoreEdgeBottom = atol(szTmp);
			return TRUE;
		}
	}


	// cluster calculations :
	if(strcmp( szParamName, "CCD_CLUSTER_IF_N_SIGMA_ABOVE_BACKGR" )==0){
		szTmp = GetParam("CCD_CLUSTER_IF_N_SIGMA_ABOVE_BACKGR",TRUE);
		if(szTmp && szTmp[0]){
			m_ClusterIfNSigmaAboveBackgr = atof( szTmp );
			return TRUE;
		}
	}

	
	if(strcmp( szParamName, "CCD_CLUSTER_MIN_SIZE" )==0){	
		szTmp = GetParam("CCD_CLUSTER_MIN_SIZE",TRUE);
		if(szTmp && szTmp[0]){
			m_MinClusterSize = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName, "CCD_MAX_NOISE_LEVEL" )==0){
		szTmp = GetParam("CCD_MAX_NOISE_LEVEL",TRUE);
		if(szTmp && szTmp[0]){
			m_MaxNoiseLevel = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName, "CCD_CONFIRM_TRESHOLD_BY_PIXEL" )==0){
		szTmp = GetParam("CCD_CONFIRM_TRESHOLD_BY_PIXEL",TRUE);
		if(szTmp && szTmp[0]){
			m_ConfTresholdPerPixel = atol( szTmp );
			return TRUE;
		}
	}


	if(strcmp( szParamName, "CCD_CONFIRM_MAX_NOISE_PER_PIXEL" )==0 ){
		szTmp = GetParam("CCD_CONFIRM_MAX_NOISE_PER_PIXEL",TRUE);
		if(szTmp && szTmp[0]){
			m_ConfMaxPrevPerPixel = atol( szTmp );
			if(m_bUseHomeoSameTreshForNewAndPrev){
				m_ConfMaxPrevPerPixel = m_ConfTresholdPerPixel;
			}
			return TRUE;
		}
	}				

	if(strcmp( szParamName, "CCD_PIXELS_AROUND_TO_CONFIRM" )==0 ){
		szTmp = GetParam("CCD_PIXELS_AROUND_TO_CONFIRM",TRUE );
		if(szTmp && szTmp[0]){
			m_nPixelsAroundToConfirm = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_DIFFERENCE_CHECK" )==0 ){
		szTmp = GetParam("CCD_DIFFERENCE_CHECK",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckDifference = (atol(szTmp)>0);
			return TRUE;
		}
	}


	if(strcmp( szParamName,"CCD_SEPERATE_TRESH_AND_MAX_PREV" )==0 ){
		szTmp = GetParam("CCD_SEPERATE_TRESH_AND_MAX_PREV",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckTreshAndMaxPrev = (atol(szTmp)>0);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_CONF_DIFFERENCE_CHECK" )==0 ){
		szTmp = GetParam("CCD_CONF_DIFFERENCE_CHECK",TRUE);
		if(szTmp && szTmp[0]){
			m_bConfCheckDifference = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_CONF_SEPERATE_TRESH_AND_MAX_PREV" )==0 ){
		szTmp = GetParam("CCD_CONF_SEPERATE_TRESH_AND_MAX_PREV",TRUE);
		if(szTmp && szTmp[0]){
			m_bConfCheckTreshAndMaxPrev = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_USE_CLUSTER_WITH_MORE" )==0 ){
		// not required :
		szTmp = GetParam("CCD_USE_CLUSTER_WITH_MORE",TRUE);
		if(szTmp && szTmp[0]){
	      m_bUseClusterWithMore = ( atol( szTmp )>0 );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_COMBINED_CHECK" )==0 ){
		// combined check params :
		szTmp = GetParam("CCD_COMBINED_CHECK",TRUE);
		if(szTmp && szTmp[0]){
			m_bCombinedCheck = (atol(szTmp) > 0);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_USE_DIFF_UP_TO_PER_PIXEL" )==0 ){
		szTmp = GetParam("CCD_USE_DIFF_UP_TO_PER_PIXEL",TRUE);
		if(szTmp && szTmp[0]){
			m_UseDoubleCheckAboveADUPerPixel =  atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_CONFIRM_ON_N_NEXT_FRAMES" )==0 ){
		// confirmation of events on next frames :
		szTmp = GetParam("CCD_CONFIRM_ON_N_NEXT_FRAMES",TRUE);
		if(szTmp && szTmp[0]){
			m_ConfirmEventsOnNextNFrames = atol( szTmp );
			return TRUE;
		}
	}


	if(strcmp( szParamName,"CCD_REJECT_NOT_CONFIRMED" )==0 ){
		szTmp = GetParam("CCD_REJECT_NOT_CONFIRMED", TRUE );
		if(szTmp && szTmp[0]){
			m_bRejectNotConfirmed = ( atol( szTmp )>0 );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_MAX_ON_NEXT_PER_PIXEL" )==0 ){
		szTmp = GetParam("CCD_MAX_ON_NEXT_PER_PIXEL",TRUE);
		if(szTmp && szTmp[0]){
			m_MaxOnNextPerPixel = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_SIG_TRESH_PER_PIXEL_LEVEL_1" )==0 ){
		szTmp = GetParam("CCD_SIG_TRESH_PER_PIXEL_LEVEL_1",TRUE);
		if(szTmp && szTmp[0]){
			m_TresholdPerPixel = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_MAX_PREV_PER_PIXEL_LEVEL_1" )==0 ){
		szTmp = GetParam("CCD_MAX_PREV_PER_PIXEL_LEVEL_1",TRUE);
		if(szTmp && szTmp[0]){
			m_MaxPrevPerPixel = atol( szTmp);
			if(m_bUseHomeoSameTreshForNewAndPrev){
				m_MaxPrevPerPixel = m_TresholdPerPixel;
			}
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_LOCAL_MAX_REQ" )==0 ){
		szTmp = GetParam("CCD_LOCAL_MAX_REQ",TRUE);
		if(szTmp && szTmp[0]){
			m_bLocalMaxReq = ( atol( szTmp )>0 );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_CONFIRM_REQ" )==0 ){
		szTmp = GetParam("CCD_CONFIRM_REQ",TRUE);
		if(szTmp && szTmp[0]){
			m_bConfirmReq = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_CALC_CLUSTER_REQ" )==0 ){
		szTmp = GetParam("CCD_CALC_CLUSTER_REQ",TRUE);
		if(szTmp && szTmp[0]){
			m_bCalcClusterReq = ( atol( szTmp )>0 );
			return TRUE;
		}
	}
		
	if(strcmp( szParamName,"CCD_CONFIRM_REDIAL" )==0 ){
		szTmp = GetParam("CCD_CONFIRM_REDIAL",TRUE);
		if(szTmp && szTmp[0]){
			m_ConfRedial = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_CONFIRM_TRESHOLD" )==0 ){
		szTmp = GetParam("CCD_CONFIRM_TRESHOLD",TRUE);
		if(szTmp && szTmp[0]){
			m_ConfTreshold = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_CONFIRM_SHAPE" )==0 ){
		szTmp = GetParam("CCD_CONFIRM_SHAPE",TRUE);
		if(szTmp && szTmp[0]){		
			m_ConfShape = (eConfShape_T)atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_VETO_SHAPE" )==0 ){
		szTmp = GetParam("CCD_VETO_SHAPE",TRUE);
		if(szTmp && szTmp[0]){
			m_eVetoShape = (eConfShape_T)atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_VETO_REDIAL" )==0 ){
		szTmp = GetParam("CCD_VETO_REDIAL",TRUE);
		if(szTmp && szTmp[0]){
			m_nVetoRedial = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_NEIGHB_SHAPE" )==0){
		szTmp = GetParam("CCD_NEIGHB_SHAPE",TRUE);
		if(szTmp && szTmp[0]){
			m_eNeigbShape = (eConfShape_T)(atol( szTmp ));
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_NEIGHB_REDIAL" )==0){
		szTmp = GetParam("CCD_NEIGHB_REDIAL",TRUE);
		if(szTmp && szTmp[0]){
			m_dNeighbRedial = atof( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_CHECK_MAX_FROM_PREV_N" )==0 ){
		szTmp = GetParam("CCD_CHECK_MAX_FROM_PREV_N",TRUE);
		if(szTmp && szTmp[0]){
			m_bAnalyseMaxFromPrevInPipeline = ( atol( szTmp )>0 );
			return TRUE;
		}
	}
		
		// average of previous frames :
	if(strcmp( szParamName,"CCD_AVERAGE_OF_PREV_REJECT_MAX_AND_MIN" )==0 ){
		szTmp = GetParam("CCD_AVERAGE_OF_PREV_REJECT_MAX_AND_MIN",TRUE);
		if(szTmp && szTmp[0]){
			m_bAverageOfPrevRejectMAXandMIN = (atol(szTmp)>0);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_AVERAGE_OF_PREV_N" )==0 ){
		szTmp = GetParam("CCD_AVERAGE_OF_PREV_N",TRUE);
		if(szTmp && szTmp[0]){
			m_nMaxOfAverageOfPrevN = atol( szTmp );
			return TRUE;
		}
	}


	if(strcmp( szParamName,"CCD_USE_ROT_IN_AVER_OF_PREV" )==0 ){
		szTmp = GetParam("CCD_USE_ROT_IN_AVER_OF_PREV",TRUE);
		if(szTmp && szTmp[0]){
			m_bUseRotInAverageOfPrevN = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_USE_ROT_PER_SEC" )==0 ){
		szTmp = GetParam("CCD_USE_ROT_PER_SEC",TRUE);
		if(szTmp && szTmp[0]){
			m_bUseRotPerSec = (atol(szTmp)>0);
			return TRUE;
		}
	}

		// homeopatic section :
	if(strcmp( szParamName,"CCD_KEEP_HOMEOPATIC_SUM" )==0 ){
		szTmp = GetParam("CCD_KEEP_HOMEOPATIC_SUM",TRUE);
		if(szTmp && szTmp[0]){
			m_bKeepHomeopaticSum = (atol(szTmp)>0);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_HOMEO_AVERAGE_OF_N" )==0 ){
		szTmp = GetParam("CCD_HOMEO_AVERAGE_OF_N",TRUE);
		if(szTmp && szTmp[0]){
			m_nHomeoAverageOfPrevNFrames = atol(szTmp);
			return TRUE;
		}
	}
	
	if(strcmp( szParamName,"CCD_CHECK_HOMEO_RAW_COND" )==0){
		szTmp = GetParam("CCD_CHECK_HOMEO_RAW_COND",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckHomeoRawCond = ( atol(szTmp)>0 );
			return TRUE;
		}
	}	

	if(strcmp( szParamName,"CCD_HOMEOPATIC_FACTOR" )==0 ){
		szTmp = GetParam("CCD_HOMEOPATIC_FACTOR",TRUE);
		if(szTmp && szTmp[0]){
			if( m_HomeopaticFactor>0 && m_HomeopaticFactor<=1 ){
				m_HomeopaticFactor = atof( szTmp );
				return TRUE;
			}
		}
		return FALSE;
	}

	if(strcmp( szParamName,"CCD_CALC_MAX_NEIGHB_HOMEO" )==0 ){
		szTmp = GetParam("CCD_CALC_MAX_NEIGHB_HOMEO",TRUE);
		if(szTmp && szTmp[0]){
			m_bCalcMaxNeighbHomeo = (atol(szTmp)>0);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_CALC_MAX_FOR_ABOVE_N_SIGMA")==0){
		szTmp = GetParam("CCD_CALC_MAX_FOR_ABOVE_N_SIGMA",TRUE);
		if(szTmp && szTmp[0]){
			m_nCalcMaxForAboveNSigma = atof(szTmp);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_CALC_MAX_FOR_ABOVE_N_SIGMA_ON_HOMEO")==0){
		szTmp = GetParam("CCD_CALC_MAX_FOR_ABOVE_N_SIGMA_ON_HOMEO",TRUE);
		if(szTmp && szTmp[0]){
			m_bCalcMaxForAboveNSigmaOnHomeo = (atol(szTmp)>0);
			return TRUE;
		}
	}
			
	if(strcmp( szParamName,"CCD_CALC_MAX_NEIGHB_REDIAL" )==0 ){
		szTmp = GetParam("CCD_CALC_MAX_NEIGHB_REDIAL",TRUE);
		if(szTmp && szTmp[0]){
			m_nCalcMaxNieghbRedial = atol(szTmp);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_START_HOMEO_WITH_FIRST_FRAME" )==0 ){
		szTmp = GetParam("CCD_START_HOMEO_WITH_FIRST_FRAME",TRUE);
		if(szTmp && szTmp[0]){
			m_bStartHomeoWithFirstFrame = (atol(szTmp)>0);
			return TRUE;
		}
	}

		// laplacjan frame :
	if(strcmp( szParamName,"CCD_KEEP_LAPLACE_FRAME" )==0 ){
		szTmp = GetParam("CCD_KEEP_LAPLACE_FRAME",TRUE);
		if(szTmp && szTmp[0]){
			m_bKeepLaplaceFrame = (atol(szTmp)>0 );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_CHECK_LAPLACE_CONDITION" )==0 ){
		szTmp = GetParam("CCD_CHECK_LAPLACE_CONDITION",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckLaplaceCondition = (atol(szTmp) >0 );
			return TRUE;
		}
	}


	if(strcmp( szParamName,"CCD_CALC_LAPLACE_OF_NEW" )==0 ){
		szTmp = GetParam("CCD_CALC_LAPLACE_OF_NEW",TRUE);
		if(szTmp && szTmp[0]){
			m_bCalcLaplaceOfNEW = (atol(szTmp)>0);
			return TRUE;
		}
	}



	if(strcmp( szParamName,"CCD_MAX_LAPLACE_OF_PREV_AVERAGE" )==0 ){
		szTmp = GetParam("CCD_MAX_LAPLACE_OF_PREV_AVERAGE",TRUE);
		if(szTmp && szTmp[0]){
			m_nMaxLaplaceOfPrevAverage = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_CHECK_LAPLACE_ON_HOMEO_CONDITION" )==0 ){
		szTmp = GetParam("CCD_CHECK_LAPLACE_ON_HOMEO_CONDITION",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckLaplaceOnHomeoCondition = atol( szTmp );
			return TRUE;
		}
	}


	if(strcmp( szParamName,"CCD_LAPLACE_TYPE" )==0 ){
		szTmp = GetParam("CCD_LAPLACE_TYPE",TRUE);
		if(szTmp && szTmp[0]){
			m_eLaplaceType = (eLaplaceType_T)atol(szTmp);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_MAX_LAPLACE_ON_OTHER" )==0 ){
		szTmp = GetParam("CCD_MAX_LAPLACE_ON_OTHER",TRUE);
		if(szTmp && szTmp[0]){
			m_nMaxLaplaceOnOther = atol( szTmp );
			return TRUE;
		}
	}
						
	if(strcmp( szParamName,"CCD_MIN_LAPLACE_ON_OTHER" )==0 ){
		szTmp = GetParam("CCD_MIN_LAPLACE_ON_OTHER",TRUE);
		if(szTmp && szTmp[0]){
			m_nMinLaplaceOnOther = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_NEW_LAPLACE_TRESHOLD" )==0 ){
		szTmp = GetParam("CCD_NEW_LAPLACE_TRESHOLD",TRUE);
		if(szTmp && szTmp[0]){
			m_nNewLaplace = atol( szTmp );			
			return TRUE;
		}
	}


	if( UpdateParam( szParamName, "CCD_MAX_LAPLACE_ON_OTHER_IN_SIGMA", m_nMaxLaplaceOnOtherInSigma, bRetVal )){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_NEW_LAPLACE_TRESHOLD_IN_SIGMA", m_nNewLaplaceInSigma, bRetVal )){
		return bRetVal;
	}


	if(strcmp( szParamName,"CCD_LAPLACE_USE_SAME_TRESHOLD" )==0 ){
		szTmp = GetParam("CCD_LAPLACE_USE_SAME_TRESHOLD",TRUE);
		if(szTmp && szTmp[0]){
			m_bLaplaceUseSameTresh = (atol(szTmp)>0);
		if(m_bLaplaceUseSameTresh)
			m_nMaxLaplaceOnOther = m_nNewLaplace;
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_CONFIRM_LAPLACE_MEDIAN_MINUS" )==0 ){
		szTmp = GetParam("CCD_CONFIRM_LAPLACE_MEDIAN_MINUS",TRUE);
		if(szTmp && szTmp[0]){
			m_bConfirmLaplaceMedianMinus = (atol(szTmp)>0);
			return TRUE;
		}
	}


		// parameters for checking if not more then N exceeds threshold :
	if(strcmp( szParamName,"CCD_CHECK_IF_NONE_OF_N_EXCEEDS_TRESH" )==0 ){
		szTmp = GetParam("CCD_CHECK_IF_NONE_OF_N_EXCEEDS_TRESH",TRUE);
		if(szTmp && szTmp[0]){
			m_nCheckIfNoneOfNExceedsTresh = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp(szParamName,"CCD_TRESH_FOR_CHECK_IF_NONE_OF_N_EXCEEDS_IN_SIGMA" )==0 ){
		szTmp = GetParam("CCD_TRESH_FOR_CHECK_IF_NONE_OF_N_EXCEEDS_IN_SIGMA",TRUE);
		if(szTmp && szTmp[0]){
			m_TreshForCheckIfNoneOfNExceedsInSigma = atof( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_NOT_MORE_THEN_N_EXEEDS_SHAPE" )==0 ){
		szTmp = GetParam("CCD_NOT_MORE_THEN_N_EXEEDS_SHAPE",TRUE);
		if(szTmp && szTmp[0]){
			m_eNotMoreThenNExceedsShape=(eConfShape_T)atol(szTmp);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_NOT_MORE_THEN_N_EXEEDS_REDIAL" )==0 ){
		szTmp = GetParam("CCD_NOT_MORE_THEN_N_EXEEDS_REDIAL",TRUE);
		if(szTmp && szTmp[0]){
			m_nNotMoreThenNExceedsRedial = atol(szTmp);			
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_NOT_MORE_THEN_N_EXEEDS_BACK_FRAMES_COUNT" )==0 ){
		szTmp = GetParam("CCD_NOT_MORE_THEN_N_EXEEDS_BACK_FRAMES_COUNT",TRUE);
		if(szTmp && szTmp[0]){
			m_nNotMoreThenNExceedsBackFramesCount = atol(szTmp);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_CHECK_FOR_SUPER_NEW" )==0 ){
		szTmp = GetParam("CCD_CHECK_FOR_SUPER_NEW",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckForSUPERNEW = ( atol(szTmp)>0 );
			return TRUE;
		}
	}	

	
	if(strcmp( szParamName,"CCD_MIN_PREV_ADU_PER_PIXEL")==0){
		szTmp = GetParam("CCD_MIN_PREV_ADU_PER_PIXEL",TRUE);
		if(szTmp && szTmp[0]){
			m_MinPrevPerPixel = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_CHECK_TYPE" )==0){
		szTmp = GetParam("CCD_CHECK_TYPE",TRUE);
		if(szTmp && szTmp[0]){
			m_eBrigtheningCheckType = (eBRIGHT_CHECK_TYPE)atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_INCREASE_TRESHOLD_ADU_PER_PIXEL" )==0){
		szTmp = GetParam("CCD_INCREASE_TRESHOLD_ADU_PER_PIXEL",TRUE);
		if(szTmp && szTmp[0]){
			m_IncreaseTreshADUPerPixel = atol(szTmp);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_INCREASE_TRESHOLD_PERCENT")==0){
		szTmp = GetParam("CCD_INCREASE_TRESHOLD_PERCENT",TRUE);
		if(szTmp && szTmp[0]){
			m_IncreaseTreshPercent = atof( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_DARK_FRAME" )==0){
		szTmp = GetParam("CCD_DARK_FRAME",TRUE);
		if(szTmp && szTmp[0]){
			m_bDarkFrame = (atol( szTmp )>0);
			return TRUE;
		}
	}


	if(strcmp( szParamName,"CCD_FLAT_FRAME" )==0 ){
		szTmp = GetParam("CCD_FLAT_FRAME",TRUE);
		if(szTmp && szTmp[0]){		
	      m_bFlatFrame = ( atol(szTmp)>0 );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_KEEP_PIXEL_MAX" )==0 ){
		szTmp = GetParam("CCD_KEEP_PIXEL_MAX",TRUE);
		if(szTmp && szTmp[0]){
			m_bPixelMaxFrame = (atol( szTmp )>0);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_DARK_FRAME_FILE" )==0 ){
		szTmp = GetParam("CCD_DARK_FRAME_FILE",TRUE);
		if(szTmp && szTmp[0]){
			m_szDarkFrameFile = szTmp;
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_FLAT_FRAME_FILE" )==0 ){
		szTmp = GetParam("CCD_FLAT_FRAME_FILE",TRUE);
		if(szTmp && szTmp[0]){
			m_szFlatFrameFile = szTmp;
			return TRUE;
		}
	}
		

		// report dumping :
	
	if(strcmp( szParamName,"CCD_DUMP_NEW_EVENTS" )==0){
		szTmp = GetParam("CCD_DUMP_NEW_EVENTS",TRUE);
		if(szTmp && szTmp[0]){
			m_bDumpNewEvents = atol( szTmp );
			return TRUE;
		}
	}
	if(strcmp( szParamName,"CCD_DUMP_ALL_EVENTS" )==0){
		szTmp = GetParam("CCD_DUMP_ALL_EVENTS",TRUE);
		if(szTmp && szTmp[0]){
			m_bDumpAllEvents = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_DUMP_EVENT_FREQ")==0){
		szTmp = GetParam("CCD_DUMP_EVENT_FREQ",TRUE);
		if(szTmp && szTmp[0]){
			m_DumpEventsFreq = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_SAVE_N_FRAMES_BEFORE_AND_AFTER" )==0 ){
		szTmp = GetParam("CCD_SAVE_N_FRAMES_BEFORE_AND_AFTER",TRUE);
		if(szTmp && szTmp[0]){
			m_nSaveFramesBeforeAndAfter = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_SAVE_FRAMES_WITH_EVENT" )==0 ){
		szTmp = GetParam("CCD_SAVE_FRAMES_WITH_EVENT",TRUE);
		if(szTmp && szTmp[0]){
			m_bSaveFramesWithEvents = ( atol( szTmp )>0 );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_SAVE_EVENT_SIZE" )==0 ){
		szTmp = GetParam("CCD_SAVE_EVENT_SIZE",TRUE);
		if(szTmp && szTmp[0]){
			m_bSaveEventSize = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_VERIFIED_EVENTS_LOG" )==0 ){
		szTmp = GetParam("CCD_VERIFIED_EVENTS_LOG",TRUE);
		if(szTmp && szTmp[0]){
			m_szVerifiedEventsLog = szTmp;
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_RUN_EVENT_LOG" )==0 ){
		szTmp = GetParam("CCD_RUN_EVENT_LOG",TRUE);
		if(szTmp && szTmp[0]){
			m_szRunEventsLog = szTmp;
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_GEN_EVENT_LOG" )==0 ){
		szTmp = GetParam("CCD_GEN_EVENT_LOG",TRUE);
		if(szTmp && szTmp[0]){
			m_szGenEventsLog = szTmp;
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_DUMP_HOMEO_FRAME" )==0 ){
		szTmp = GetParam("CCD_DUMP_HOMEO_FRAME",TRUE);
		if(szTmp && szTmp[0]){
			m_bDumpHomeoFrame = ( atol( szTmp )>0 );
		m_bGenNotFoundReport = ( atol( GetParam("CCD_GEN_NOT_FOUND_REPORT") )>0 );
			return TRUE;
		}
	}			
		
		// background section :
	if(strcmp( szParamName,"CCD_CALC_S_DISTR" )==0 ){
		szTmp = GetParam("CCD_CALC_S_DISTR",TRUE);
		if(szTmp && szTmp[0]){
			m_CalcBackgrTable[(int)eRawS] = ( atol( szTmp )>0 );
			return TRUE;
		}
	}


	if(strcmp( szParamName,"CCD_CALC_G1_DISTR" )==0 ){
		szTmp = GetParam("CCD_CALC_G1_DISTR",TRUE);
		if(szTmp && szTmp[0]){
			m_CalcBackgrTable[(int)eSinglePoint] = ( atol( szTmp )>0 );
			return TRUE;
		}
	}
	if(strcmp( szParamName,"CCD_CALC_G2_DISTR" )==0 ){
		szTmp = GetParam("CCD_CALC_G2_DISTR",TRUE);
		if(szTmp && szTmp[0]){
			m_CalcBackgrTable[(int)eTwoPoints] = ( atol( szTmp )>0 );
			return TRUE;
		}
	}
	if(strcmp( szParamName,"CCD_CALC_G4_DISTR" )==0 ){
		szTmp = GetParam("CCD_CALC_G4_DISTR",TRUE);
		if(szTmp && szTmp[0]){
			m_CalcBackgrTable[(int)eFourPoints] = ( atol( szTmp )>0 );
			return TRUE;
		}
	}
	if(strcmp( szParamName,"CCD_CALC_G5_DISTR" )==0 ){
		szTmp = GetParam("CCD_CALC_G5_DISTR",TRUE);
		if(szTmp && szTmp[0]){
			m_CalcBackgrTable[(int)eFivePoints] = ( atol( szTmp )>0 );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_CALC_G54_DISTR" )==0 ){
		szTmp = GetParam("CCD_CALC_G54_DISTR",TRUE);
		if(szTmp && szTmp[0]){
			m_CalcBackgrTable[(int)eFivePlusFourMin] = ( atol( szTmp )>0 );
			return TRUE;
		}
	}
	if(strcmp( szParamName,"CCD_CALC_G985_DISTR" )==0 ){
		szTmp = GetParam("CCD_CALC_G985_DISTR",TRUE);
		if(szTmp && szTmp[0]){
			m_CalcBackgrTable[(int)eNineEightFive] = ( atol( szTmp )>0 );
			return TRUE;
		}
	}
	if(strcmp( szParamName,"CCD_CALC_G987_DISTR" )==0 ){
		szTmp = GetParam("CCD_CALC_G987_DISTR",TRUE);
		if(szTmp && szTmp[0]){
			m_CalcBackgrTable[(int)eNineEightSevenVeryBig] = ( atol( szTmp )>0 );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_CALC_G84_DISTR" )==0 ){
		szTmp = GetParam("CCD_CALC_G84_DISTR",TRUE);
		if(szTmp && szTmp[0]){
			m_CalcBackgrTable[(int)eEightFour] = ( atol( szTmp )>0 );
			return TRUE;
		}
	}
	if(strcmp( szParamName,"CCD_CALC_G810_DISTR" )==0 ){
		szTmp = GetParam("CCD_CALC_G810_DISTR",TRUE);
		if(szTmp && szTmp[0]){
			m_CalcBackgrTable[(int)eEightTen] = ( atol( szTmp )>0 );
			return TRUE;
		}
	}
	if(strcmp( szParamName,"CCD_CALC_G58_DISTR" )==0 ){
		szTmp = GetParam("CCD_CALC_G58_DISTR",TRUE);
		if(szTmp && szTmp[0]){
			m_CalcBackgrTable[(int)eFiveEight] = ( atol( szTmp )>0 );
			return TRUE;
		}
	}
	if(strcmp( szParamName,"CCD_CALC_G412_DISTR" )==0 ){
		szTmp = GetParam("CCD_CALC_G412_DISTR",TRUE);
		if(szTmp && szTmp[0]){
			m_CalcBackgrTable[(int)eFourTwelve] = ( atol( szTmp )>0 );
			return TRUE;
		}
	}
	if(strcmp( szParamName,"CCD_CALC_G412F_DISTR" )==0 ){
		szTmp = GetParam("CCD_CALC_G412F_DISTR",TRUE);
		if(szTmp && szTmp[0]){
			m_CalcBackgrTable[(int)eFourTwelveFar] = ( atol( szTmp )>0 );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_AVERAGE_S1" )==0 ){
		szTmp = GetParam("CCD_AVERAGE_S1",TRUE);
		if(szTmp && szTmp[0]){
			m_AverageS1 = atof( szTmp );	
			gAverageBacgroundG[eRawS] = m_AverageS1;
			return TRUE;
		}
	}
	if(strcmp( szParamName,"CCD_SIGMA_S1" )==0 ){
		szTmp = GetParam("CCD_SIGMA_S1",TRUE);
		if(szTmp && szTmp[0]){
			m_SigmaS1 = atof( szTmp );
			gSigmaBacgroundG[eRawS] = m_SigmaS1;
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_AUTO_CALCULATE_TRESHOLDS" )==0 ){
		szTmp = GetParam("CCD_AUTO_CALCULATE_TRESHOLDS",TRUE );
		if(szTmp && szTmp[0]){
			m_bAutoCalcTresh = ( atol( szTmp )>0 );
			return TRUE;
		}
	}

	if( UpdateParam( szParamName, "CCD_DUMP_ALL_LOG_TO_STDOUT", m_bDumpAllLogToStdout, bRetVal )){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_SAVE_VERIF_EVENTS_TO_DB", m_bSaveVerifEventsToDB, bRetVal )){
		return bRetVal;
	}
	if( UpdateParam( szParamName, "CCD_SAVE_ALL_EVENTS_TO_DB", m_bSaveAllEventsToDB, bRetVal )){
		return bRetVal;
	}
	
	if( UpdateParam( szParamName, "CCD_TRACKS_TO_DB", m_bSaveTracksToDB , bRetVal )){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_CALC_TRESHOLDS_BY_MAP", tmp_val_int, bRetVal )){
      if( bRetVal ){
			m_bCalcTresholdsByBackgrMap =  ( tmp_val_int > 0 );
      }
      return bRetVal;
   }

	if( UpdateParam( szParamName, "CCD_OBS_MODE", tmp_val_int, bRetVal )){
		if( bRetVal ){
   		m_ObsMode = (eObservationMode_T)tmp_val_int;
		}
		return bRetVal;
	}


	if(strcmp( szParamName,"CCD_BACKGR_UPDATE_FREQ" )==0){
		szTmp = GetParam("CCD_BACKGR_UPDATE_FREQ",TRUE);
		if(szTmp && szTmp[0]){
			m_BackgUpdateFreq = atol( szTmp );
			return TRUE;
		}
		return FALSE;
	}

	if(strcmp( szParamName,"CCD_BACKGR_DUMP" )==0){
		szTmp =  GetParam("CCD_BACKGR_DUMP", TRUE);
		if(szTmp && szTmp[0]){
			m_bBackgrDump = ( atol(szTmp)>0 );
			return TRUE;
		}		
		return FALSE;
	}


	if(strcmp( szParamName,"CCD_SUBTRACT_BACKGROUND" )==0 ){
		szTmp = GetParam("CCD_SUBTRACT_BACKGROUND",TRUE);
		if(szTmp && szTmp[0]){
			m_bSubtractBackground = ( atol( szTmp )>0 );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_BACKGROUND_SUBTR_TYPE" )==0 ){
		szTmp = GetParam("CCD_BACKGROUND_SUBTR_TYPE",TRUE);
		if(szTmp && szTmp[0]){
			m_eBackgrSubtrType = (eBackgrSubtrType_T)( atol(szTmp));
			return TRUE;
		}
	}


	if(strcmp( szParamName,"CCD_STDOUT_ON" )==0 ){
		szTmp = GetParam("CCD_STDOUT_ON",TRUE);
		if(szTmp && szTmp[0]){
			m_bStdoutOn      = ( atol(szTmp)>0 );
			return TRUE;
		}
	}
		

	if(strcmp( szParamName,"CCD_EVENTS_BUFFER_SIZE" )==0 ){
		szTmp = GetParam("CCD_EVENTS_BUFFER_SIZE",TRUE);
		if(szTmp && szTmp[0]){
			m_EventsBufferSize = atol(szTmp);
			return TRUE;
		}
	}



	if(strcmp( szParamName,"CCD_SKIP_OVERLAPS" )==0 ){
		szTmp = GetParam("CCD_SKIP_OVERLAPS",TRUE);
		if(szTmp && szTmp[0]){
			m_bSkipOverlaps = ( atol( szTmp )>0 );
			return TRUE;
		}
	}
	if(strcmp( szParamName,"CCD_OVERLAP_REDIAL" )==0 ){
		szTmp = GetParam("CCD_OVERLAP_REDIAL",TRUE);
		if(szTmp && szTmp[0]){
			m_OverlapRedial = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_CHECK_IF_MORE_POINT" )==0){
		szTmp = GetParam("CCD_CHECK_IF_MORE_POINT",TRUE);
		if(szTmp && szTmp[0]){		
			m_eCheckIfMorePoint = (eCheckIfMorePoint_T)atol(szTmp);
			return TRUE;
		}
	}
	if(strcmp( szParamName,"CCD_MAX_NUMBER_OF_EVENTS")==0){
		szTmp = GetParam("CCD_MAX_NUMBER_OF_EVENTS",TRUE);
		if(szTmp && szTmp[0]){
			m_MaxNumberOfEventsOnFrame = atol( szTmp );
			return TRUE;
		}
	}
	if(strcmp( szParamName,"CCD_MAX_NUMBER_OF_EVENTS_AFTER_COIC")==0){
		szTmp = GetParam("CCD_MAX_NUMBER_OF_EVENTS_AFTER_COIC",TRUE);
   	if(szTmp && szTmp[0]){
   		m_MaxNumberOfEventsOnFrameAfterCoic = atol( szTmp );
			return TRUE;
		}
	}
      

	if(strcmp( szParamName,"CCD_SKIP_IF_MORE_THEN" )==0 ){
		szTmp = GetParam("CCD_SKIP_IF_MORE_THEN",TRUE);
		if(szTmp && szTmp[0]){
			m_bSkipIfMoreThen = atol( szTmp );
			return TRUE;
		}
	}
	if(strcmp( szParamName,"CCD_SKIP_IF_MORE_THEN_REDIAL" )==0 ){
		szTmp = GetParam("CCD_SKIP_IF_MORE_THEN_REDIAL",TRUE);
		if(szTmp && szTmp[0]){
			m_bSkipIfMoreThenRedial = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_SKIP_IF_MORE_THEN_MIN_DIST" )==0 ){
		szTmp = GetParam("CCD_SKIP_IF_MORE_THEN_MIN_DIST",TRUE);
		if(szTmp && szTmp[0]){
			m_bSkipIfMoreThenMinDist = atol( szTmp );
			return TRUE;
		}
	}		

	if(strcmp( szParamName,"CCD_CHECK_EVENT_SHAPE" )==0 ){
		szTmp = GetParam("CCD_CHECK_EVENT_SHAPE",TRUE);
		if(szTmp && szTmp[0]){
			m_CheckEventShape = atof( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_CHECK_CENTER_MASS_IN_CLUSTER" )==0 ){
		szTmp = GetParam("CCD_CHECK_CENTER_MASS_IN_CLUSTER",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckCenterInCluster = ( atol( szTmp )>0 );
			return TRUE;
		}
	}

		
	if(strcmp( szParamName,"CCD_MAX_PIXELS_IN_CLUSTER_ALLOWED" )==0 ){
		szTmp = GetParam("CCD_MAX_PIXELS_IN_CLUSTER_ALLOWED",TRUE);
		if(szTmp && szTmp[0]){
			m_MaxPixelsInClusterAllowed = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_CHECK_GEN_ONLY" )==0 ){
		szTmp = GetParam("CCD_CHECK_GEN_ONLY",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckGenOnly = (atol(szTmp)>0);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_MC_GEN_IDENT_REDIAL" )==0 ){
		szTmp = GetParam("CCD_MC_GEN_IDENT_REDIAL",TRUE);
		if(szTmp && szTmp[0]){
			m_bGenEventRedial = atol( szTmp );
			return TRUE;
		}
	}


	if(strcmp( szParamName,"CCD_PUT_OBJ_ON_N_FRAMES" )==0 ){
		szTmp = GetParam("CCD_PUT_OBJ_ON_N_FRAMES",TRUE);
		if(szTmp && szTmp[0]){
			m_nPutObjOnNFrames = atol( szTmp );
			return TRUE;
		}
	}

		//--------------------------------------------------------------------------------------

	if(strcmp( szParamName,"CCD_KEEP_ALL_IN_MEMORY_MC" )==0 ){
		szTmp = GetParam("CCD_KEEP_ALL_IN_MEMORY_MC",TRUE);
		if(szTmp && szTmp[0]){
			m_bKeepAllInMemory = ( atol(szTmp)>0 );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_KEEP_SAMPLES_IN_MEMORY" )==0 ){
		szTmp = GetParam("CCD_KEEP_SAMPLES_IN_MEMORY",TRUE);
		if(szTmp && szTmp[0]){
			m_bKeepSamplesInMemory = (atol(szTmp)>0);
			return TRUE;
		}
	}


	if(strcmp( szParamName,"CCD_IMAGES_SAMPLE_DIR" )==0 ){
		szTmp = GetParam("CCD_IMAGES_SAMPLE_DIR",TRUE);
		if(szTmp){
			m_szSampleDir = szTmp;
			return TRUE;
		}
	}


	if(strcmp( szParamName,"CCD_PUT_SCALED_SAMPLES" )==0 ){
		szTmp = GetParam("CCD_PUT_SCALED_SAMPLES",TRUE);
		if(szTmp && szTmp[0]){
			m_bPutScaledSamples = (atol(szTmp)>0);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_ALWAYS_USE_SAME_TRESH_PER_PIXEL" )==0 ){
		szTmp = GetParam( "CCD_ALWAYS_USE_SAME_TRESH_PER_PIXEL", TRUE );
		if(szTmp && szTmp[0]){
			m_bAlwaysUseSameTreshPerPixel = ( atol(szTmp)>0 );
			return TRUE;
		}
	}


		// new parameters for excluding defintions of algorithms :
	if(strcmp( szParamName,"CCD_SUPER_OPTIMIZED" )==0 ){
		szTmp = GetParam("CCD_SUPER_OPTIMIZED",TRUE);
		if(szTmp && szTmp[0]){
			m_bSuperOptimized = (atol(szTmp)>0);
			return TRUE;
		}
	}
	
	if(strcmp( szParamName,"CCD_AVERAGE_OF_PREV_ALG" )==0 ){
		szTmp = GetParam("CCD_AVERAGE_OF_PREV_ALG",TRUE);
			if(szTmp && szTmp[0]){
				m_bAverageOfPrev = (atol(szTmp)>0);
			return TRUE;
		}
	}
	if(strcmp( szParamName,"CCD_HOMEOPATIC_ALG" )==0 ){
		szTmp = GetParam("CCD_HOMEOPATIC_ALG",TRUE);
		if(szTmp && szTmp[0]){
			m_bHomeopatic = (atol(szTmp)>0);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_PRINTF_LEVEL" )==0 ){
		szTmp = GetParam("CCD_PRINTF_LEVEL",TRUE);
		if(szTmp && szTmp[0]){
			gPrintfLevel = atol(szTmp);
			return TRUE;
		}
	}


	if(strcmp( szParamName,"CCD_DEBUG_DUMP_CLUSTERS" )==0 ){
		szTmp = GetParam("CCD_DEBUG_DUMP_CLUSTERS",TRUE);
		if(szTmp && szTmp[0]){
			m_bDumpClusters = ( atol( szTmp )>0 );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_COICYDENCE_REDIAL" )==0 ){
		szTmp = GetParam("CCD_COICYDENCE_REDIAL",TRUE);
		if(szTmp && szTmp[0]){
			m_nCoicRedial = atol(szTmp);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_CHECK_TRACKS" )==0 ){
		szTmp = GetParam("CCD_CHECK_TRACKS",TRUE);
		if(szTmp && szTmp[0]){
			m_bCheckTracks = ( atol( szTmp )>0 );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_MAX_EVENT_RATE_FOR_TRACKS" )==0 ){
		szTmp = GetParam("CCD_MAX_EVENT_RATE_FOR_TRACKS",TRUE);
		if(szTmp && szTmp[0]){	
			m_MaxEventRate = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_MAX_CHI2_IN_TRACK" )==0 ){
		szTmp = GetParam("CCD_MAX_CHI2_IN_TRACK",TRUE);
		if(szTmp && szTmp[0]){
			m_MaxChi2InTrack = atof( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_MAX_CHI2_FOR_POINT_TO_MATCH_LINE" )==0 ){
		szTmp = GetParam("CCD_MAX_CHI2_FOR_POINT_TO_MATCH_LINE",TRUE);
		if(szTmp && szTmp[0]){
			m_MaxChi2ForPointToMatchLine = atof( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_MIN_EVENTS_NO_IN_TRACK" )==0 ){
		szTmp = GetParam("CCD_MIN_EVENTS_NO_IN_TRACK",TRUE);
		if(szTmp && szTmp[0]){
			m_nMinEventNoInTrack = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_NUM_BACK_FRAMES_FOR_TRACK" )==0 ){
		szTmp = GetParam("CCD_NUM_BACK_FRAMES_FOR_TRACK",TRUE);
		if(szTmp && szTmp[0]){
			m_nNumBackFramesForTracks = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_NUM_BACK_FRAMES_HAS_EVENTS_FOR_TRACKS" )==0 ){
		szTmp = GetParam("CCD_NUM_BACK_FRAMES_HAS_EVENTS_FOR_TRACKS",TRUE);
		if(szTmp && szTmp[0]){
			m_nNumBackFramesHasEventsForTracks = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_MAX_EVENTS_TO_FIT_TRACK" )==0 ){
		szTmp = GetParam("CCD_MAX_EVENTS_TO_FIT_TRACK",TRUE);
		if(szTmp && szTmp[0]){
			m_nMaxEventsToFitTrack = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_DATE_OBS" )==0 ){
		szTmp = GetParam("CCD_DATE_OBS",TRUE);
		if(szTmp){
			m_DateObs = szTmp;
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_TIME_OBS" )==0 ){
		szTmp = GetParam("CCD_TIME_OBS",TRUE);
		if(szTmp){
			m_TimeObs = szTmp;
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_DATE_TIME_OBS_FMT" )==0 ){
		szTmp = GetParam("CCD_DATE_TIME_OBS_FMT",TRUE);
		if(szTmp && szTmp[0]){
			m_DateTimeObsFormat = szTmp;
			return TRUE;
		}
	}


	if(strcmp( szParamName,"CCD_DO_NOT_ANALYSE" )==0 ){
		szTmp = GetParam("CCD_DO_NOT_ANALYSE",TRUE);
		if(szTmp && szTmp[0]){
			m_bDoNotAnalyse = ( atol( szTmp )>0 );
			return TRUE;
		}
	}			
	
	if(strcmp( szParamName,"CCD_USE_FRAME_TIME_FOR_SHIFTS" )==0 ){
		szTmp = GetParam("CCD_USE_FRAME_TIME_FOR_SHIFTS",TRUE);		
		if(szTmp && szTmp[0]){
			m_bUseFrameTimeForShifts = (atol(szTmp)>0);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_IGNORE_MISSING_TIME" )==0 ){
		szTmp = GetParam("CCD_IGNORE_MISSING_TIME",TRUE);
		if(szTmp && szTmp[0]){
			m_bIgnoreMissingTime = (atol(szTmp)>0);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_SHIFT_USES_ASTRO_FORMULAS" )==0 ){
		szTmp = GetParam("CCD_SHIFT_USES_ASTRO_FORMULAS",TRUE);
		if(szTmp && szTmp[0]){
			m_bShiftUsesAstroFormulas = (atol(szTmp)>0);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_TRANSFORM_ORIENTATION" )==0 ){
		szTmp = GetParam("CCD_TRANSFORM_ORIENTATION",TRUE);
		if(szTmp && szTmp[0]){
			m_TransformCCDOrientation = atof(szTmp);
			return TRUE;
		}
	}
			
	if(strcmp( szParamName,"CCD_TRANSFORM_CCD_FOCUS" )==0 ){
		szTmp = GetParam("CCD_TRANSFORM_CCD_FOCUS",TRUE);
		if(szTmp && szTmp[0]){
			m_TransformCCDFocus = atof( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_TRANSFORM_CCD_PIXEL_SIZE" )==0 ){
		szTmp = GetParam("CCD_TRANSFORM_CCD_PIXEL_SIZE",TRUE);
		if(szTmp && szTmp[0]){
			m_TransformCCDPixelSize = atof( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_RA_OBS" )==0 ){
		szTmp = GetParam("CCD_RA_OBS",TRUE);
		if(szTmp && szTmp[0]){
			m_RAObs = atof( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_DEC_OBS" )==0 ){
		szTmp = GetParam("CCD_DEC_OBS",TRUE);
		if(szTmp && szTmp[0]){
			m_DecObs = atof( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName, "CCD_DRIVER_ASYNCHRO_MODE" )==0){	
		szTmp = GetParam("CCD_DRIVER_ASYNCHRO_MODE",TRUE);
		if(szTmp && szTmp[0]){
			m_bAsynchroMode = ( atol( szTmp )>0 );
//			CCDPipeline::SetAsynchroModeOnOff( m_bAsynchroMode );
			return TRUE;
		}else{
			return FALSE;
		}
	}

	if(strcmp( szParamName, "CCD_TAKE_IN_SYNCHRO_MODE" )==0){
		szTmp = GetParam("CCD_TAKE_IN_SYNCHRO_MODE",TRUE);
		if(szTmp && szTmp[0]){
			m_bTakeInSynchroMode = ( atol( szTmp )>0 );
			return TRUE;
		}else{
			return FALSE;
		}
	}


		// DRIVER parameters :
	if(strcmp( szParamName,"CCD_DRIVER_MAX_ITER_TIMEOUT" )==0 ){
		szTmp = GetParam("CCD_DRIVER_MAX_ITER_TIMEOUT",TRUE);
		if(szTmp && szTmp[0]){
			m_DriverMaxIterTimeout = atol(szTmp);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_ADC_OFFSET" )==0 ){
		szTmp = GetParam("CCD_ADC_OFFSET",TRUE);
		if(szTmp && szTmp[0]){		
			m_ADCOffset = atol(szTmp);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_ADC_CLAMPING" )==0 ){
		szTmp = GetParam("CCD_ADC_CLAMPING",TRUE);
		if(szTmp && szTmp[0]){		
			m_ADCClampingBitOnOff = ( atol( szTmp )>0 );
			return TRUE;
		}
		return FALSE;
	}

	if( strcmp( szParamName,"CCD_ADC_RANGE" )==0 ){
		szTmp = GetParam("CCD_ADC_RANGE",TRUE);
      if(szTmp && szTmp[0]){
         m_eADCRange = (eADCRange)atol( szTmp );
			return TRUE;
      }
		return FALSE;
	}

	if( UpdateParam( szParamName, "CCD_LNA_GAIN", tmp_val_int, bRetVal )){
		if( bRetVal ){
			m_eLNAGain = (eLNAGainValue_T)tmp_val_int;
		}
		return bRetVal;
	}

	if(strcmp( szParamName,"CCD_ADC_GAIN" )==0 ){
		szTmp = GetParam("CCD_ADC_GAIN",TRUE);
		if(szTmp && szTmp[0]){
			m_ADCGain = atol(szTmp);
			return TRUE;
		}
		return FALSE;
	}

	if( UpdateParam( szParamName, "CCD_MPP_MODE", tmp_val_int, bRetVal )){
		if( bRetVal ){
			m_eMPPMode = (eMPP_T)tmp_val_int;
		}
		return bRetVal;
	}

	if(strcmp( szParamName,"CCD_READOUT_SPEED_HORIZONTAL" )==0 ){
		szTmp = GetParam("CCD_READOUT_SPEED_HORIZONTAL",TRUE);
		if(szTmp && szTmp[0]){
			m_ReadoutSpeedHorizontal = atol(szTmp);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_READOUT_SPEED_VERTICAL" )==0 ){
		szTmp = GetParam("CCD_READOUT_SPEED_VERTICAL",TRUE);
		if(szTmp && szTmp[0]){
			m_ReadoutSpeedVertical = atol(szTmp);
			return TRUE;
		}
	}


	if( UpdateParam( szParamName, "CCD_LENS_HIT_ON", m_bLensHitOnOff, bRetVal )){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_FRAMES_OUTPUT_SUBDIR", m_szFramesOutputSubDir, bRetVal )){
		return bRetVal;
	}

		
	if( UpdateParam( szParamName, "CCD_TEMPERATURE", m_CCDTemp, bRetVal )){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_COOLING_ON_OFF", m_bCoolingOnOff, bRetVal )){
		return bRetVal;
	}

	if(strcmp( szParamName,"CCD_CAM_TYPE" )==0 ){
		szTmp = GetParam("CCD_CAM_TYPE",TRUE);
		if(szTmp && szTmp[0]){
			m_eCAMType = (eCCDTYPE_T)atol(szTmp);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_CAM_IDENT_NO" )==0 ){
		szTmp = GetParam("CCD_CAM_IDENT_NO",TRUE);
		if(szTmp && szTmp[0]){
			m_CCDIdentNo = szTmp;
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_INTERVAL_BETWEEN_IMAGES" )==0 ){
		szTmp = GetParam("CCD_INTERVAL_BETWEEN_IMAGES",TRUE);
		if(szTmp && szTmp[0]){
			m_IntervalBetweenFrames = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_BASE_FILE_NAME_FITS" )==0 ){
		szTmp = GetParam("CCD_BASE_FILE_NAME_FITS",TRUE);
		if(szTmp && szTmp[0]){
			m_szBaseFileNameFITS = szTmp;
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_DRIVER_WRITE_ALL_TO_FITS" )==0 ){
		szTmp = GetParam("CCD_DRIVER_WRITE_ALL_TO_FITS",TRUE);
		if(szTmp && szTmp[0]){
			m_bDriverWriteAllToFITS = ( atol(szTmp)>0 );
			return TRUE;
		}
	}

	if( UpdateParam( szParamName, "CCD_WRITE_FRAMES_LIST", m_bBuildFramesList, bRetVal )){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_GETDATA_RETRY_COUNT", m_DriverGetDataRetryCount, bRetVal )){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_FITS_COMPRESS_TYPE", tmp_val_int, bRetVal )){
		if( bRetVal ){
			m_eCompressFITS = (eFITSCompressionType)tmp_val_int;
		}
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_HOR_ALT_CORR", m_HorAltCorr, bRetVal )){
		return bRetVal;
	}
	if( UpdateParam( szParamName, "CCD_HOR_AZIM_CORR", m_HorAzimCorr, bRetVal )){
		return bRetVal;
	}
	if( UpdateParam( szParamName, "CCD_RA_CORR", m_RACorr, bRetVal )){
		return bRetVal;
	}
	if( UpdateParam( szParamName, "CCD_DEC_CORR", m_DecCorr, bRetVal )){
		return bRetVal;
	}

	int tmp_shutter_mode;
	if( UpdateParam( szParamName, "CCD_SHUTTER_MODE", tmp_shutter_mode, bRetVal )){
		m_ShutterModeOriginal = m_ShutterMode;
		m_ShutterMode = tmp_shutter_mode;
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_MAX_COMM_ERROR_COUNT_TO_EXIT", m_MaxCommErrorCountToEXIT , bRetVal )){
		return bRetVal;
	}	

	if( UpdateParam( szParamName, "CCD_EXIT_ON_COMM_ERROR", m_ExitOnCommError , bRetVal )){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_OPEN_SHUTTER_WHEN_STARCOUNT_BELOW", m_MinStarCountToCloseShutter, bRetVal )){
		return bRetVal;
   }

	if( UpdateParam( szParamName, "CCD_DRIVER_SHUTTER_TIME_IN_SEC", m_DriverShutterTimeInSec, bRetVal )){
		if( bRetVal ){
			m_DriverShutterTimeInSecSaved = m_DriverShutterTimeInSec;
		}
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_MOUNT_ID", m_MountID, bRetVal )){
		return bRetVal;
	}

	if(strcmp( szParamName,"CCD_DAQ_COMMUNICATION_ON" )==0 ){
		szTmp = GetParam("CCD_DAQ_COMMUNICATION_ON",TRUE);
		if(szTmp && szTmp[0]){
			m_bDAQCommunicationON = ( atol(szTmp)>0 );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_PI_SYS_MANAGER_ON" )==0 ){
		szTmp = GetParam("CCD_PI_SYS_MANAGER_ON",TRUE);
		if(szTmp && szTmp[0]){
			m_bPISysManagerON = ( atol(szTmp)>0 );
			return TRUE;
		}
	}
		
	if(strcmp( szParamName,"CCD_PI_SYS_MANAGER_PORT_NO" )==0 ){
		szTmp = GetParam("CCD_PI_SYS_MANAGER_PORT_NO",TRUE);
		if(szTmp && szTmp[0]){
			m_PISysManPortNo = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_REJECT_BLACK_PIXELS" )==0 ){
		szTmp = GetParam("CCD_REJECT_BLACK_PIXELS",TRUE);
		if(szTmp && szTmp[0]){
			m_bRejectBlackPixels = ( atol( szTmp )>0 );
			// in case Black Pixels are checked also RAW distribution 
			// calculation is required :
			// m_CalcBackgrTable[(int)eRawS] = TRUE;
			return TRUE;
		}
		return FALSE;
	}

	if(strcmp( szParamName,"CCD_BLACK_PIXELS_IF_N_SIGMA_BELOW" )==0){
		szTmp = GetParam("CCD_BLACK_PIXELS_IF_N_SIGMA_BELOW",TRUE);
		if(szTmp && szTmp[0]){
			m_bBlackPixelsIfNSigmaBelow = atof( szTmp );

			// in case Black Pixels are checked also RAW distribution
         // calculation is required :
         m_CalcBackgrTable[(int)eRawS] = TRUE;
			return TRUE;
		}
		return FALSE;
	}
		
	if(strcmp( szParamName,"CCD_BLACK_PIXELS_RATIO" )==0){
		szTmp = GetParam("CCD_BLACK_PIXELS_RATIO",TRUE);
		if(szTmp && szTmp[0]){
			m_fBlackPixelsRatio = atof( szTmp );
			return TRUE;
		}
		return FALSE;
	}

	// Internal triggers :
	if( UpdateParam( szParamName, "CCD_SEND_INTERNAL_TRIGGERS", m_bSendInternalTriggers, bRetVal )){
		return bRetVal;
	}

	//-------- SATELITES 
	if( UpdateParam( szParamName, "CCD_CHECK_IF_SATELITE", m_bCheckIfSatelite, bRetVal )){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_REJECT_SATELITES", m_bRejectSatelite , bRetVal )){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_SAT_REJ_RADIUS" , m_nSatRejRadius, bRetVal )){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_SAT_REJ_RADIUS_IN_SEC" , tmp_val , bRetVal )){
		if( bRetVal ){
			m_nSatRejRadius = AstroAngle::arcsec2rad( tmp_val );
		}
      return bRetVal;
   }

	//--------- STAR REJECTION :
	if( UpdateParam( szParamName, "CCD_CHECK_IF_STAR" , m_bCheckIfStar, bRetVal )){
      return bRetVal;
   }

	if( UpdateParam( szParamName, "CCD_REJECT_STARS" , m_bRejectStars, bRetVal )){
      return bRetVal;
   }
	if( UpdateParam( szParamName, "CCD_STAR_REJ_RADIUS" , m_fStarRejectRadius, bRetVal )){
      return bRetVal;
   }
	if( UpdateParam( szParamName, "CCD_CHECK_IF_STAR" , m_bCheckIfStar, bRetVal ) ){
      return bRetVal;
   }
	if( UpdateParam( szParamName, "CCD_STAR_REJ_RADIUS_IN_SEC" , m_fStarRejectRadius, bRetVal ) ){
      return bRetVal;
   }	
	//--------- END OF STAR REJECTION

	if( UpdateParam( szParamName, "CCD_USE_FAST_PHOTO_IN_ASTRO", m_bUseFastPhotoInAstro, bRetVal ) ){
		return bRetVal;
	}

	if( UpdateParam( szParamName, "CCD_ASTROMETRY_IN_TAKE_N_MODE" , m_bDoAstrometryInTakeNMode, bRetVal ) ){
		return bRetVal;
	}
	if( UpdateParam( szParamName, "CCD_PHOTO_IN_TAKE_N_MODE", m_bDoPhotoInTakeNMode, bRetVal ) ){
		return bRetVal;
	}
	if( UpdateParam( szParamName, "CCD_WAIT_FOR_ASTROMETRY_IN_SEC" , m_nWaitForAstrometryInSec, bRetVal ) ){
		return bRetVal;
	}

	if(strcmp( szParamName,"CCD_AUTO_ASTROMETRY_FREQ" )==0 ){
		szTmp = GetParam("CCD_AUTO_ASTROMETRY_FREQ",TRUE);
      if(szTmp && szTmp[0]){
      	m_nAutoAstrometryFreq = atol( szTmp );
			return TRUE;
      }
		return FALSE;	
	}
	if(strcmp( szParamName,"CCD_AUTO_ASTROMETRY_FREQ_IN_SEC" )==0 ){
		szTmp = GetParam("CCD_AUTO_ASTROMETRY_FREQ_IN_SEC",TRUE);
      if(szTmp && szTmp[0]){
      	m_nAutoAstrometryFreqInSec = atol( szTmp );
			return TRUE;
      }
		return FALSE;	
	}

	if( UpdateParam( szParamName, "CCD_DO_ASAS_PHOT_ASTR", m_bDoASASPhotAstr, bRetVal ) ){
		return bRetVal;
	}	

	if(strcmp( szParamName,"CCD_ASAS_PHOTO_THRES" )==0 ){
		szTmp = GetParam("CCD_ASAS_PHOTO_THRES",TRUE);
		if(szTmp && szTmp[0]){				
			m_fASASPhotoThres = atof( szTmp );
			return TRUE;
		}
		return FALSE;
	}
	if(strcmp( szParamName,"CCD_ASAS_ASTROMETRY_FI" )==0 ){
		szTmp = GetParam("CCD_ASAS_ASTROMETRY_FI",TRUE);
		if(szTmp && szTmp[0]){			
			m_fASASAstrometryFi = atof( szTmp );
			return TRUE;
		}
		return FALSE;
	}
	if(strcmp( szParamName,"CCD_ASAS_ASTROMETRY_ORD" )==0 ){
		szTmp = GetParam("CCD_ASAS_ASTROMETRY_ORD",TRUE);
		if(szTmp && szTmp[0]){
			m_nASASAstrometryOrd = atol( szTmp );
			return TRUE;
		}
		return FALSE;
	}
	if(strcmp( szParamName,"CCD_ASAS_ASTROMETRY_VERB" )==0 ){
		szTmp = GetParam("CCD_ASAS_ASTROMETRY_VERB",TRUE);
		if(szTmp && szTmp[0]){
		   m_bASASAstrometryVerb = atol(szTmp);
			return TRUE;
		}
		return FALSE;
	}
	if(strcmp( szParamName,"CCD_ASAS_ASTROMETRY_TRY" )==0 ){
		szTmp = GetParam("CCD_ASAS_ASTROMETRY_TRY",TRUE);
		if(szTmp && szTmp[0]){
		   m_nASASAstrometryTry = atol(szTmp);
			return TRUE;
		}
		return FALSE;
	}
	if(strcmp( szParamName,"CCD_ASAS_ASTROMETRY_RETRY" )==0 ){
		szTmp = GetParam("CCD_ASAS_ASTROMETRY_RETRY",TRUE);
		if(szTmp && szTmp[0]){
			m_nASASAstrometryReTry = atol(szTmp);
			return TRUE;
		}
		return FALSE;
	}
	if(strcmp( szParamName,"PIXSCALE" )==0 ){
		szTmp = GetParam("PIXSCALE",TRUE);
		if(szTmp && szTmp[0]){
			m_fPixScale = atof( szTmp );
			return TRUE;
		}
		return FALSE;
	}
	if(strcmp( szParamName,"CCD_USE_ASAS_TRANSFORM" )==0 ){	
		szTmp = GetParam("CCD_USE_ASAS_TRANSFORM",TRUE);
		if(szTmp && szTmp[0]){
			m_bUseAsasTransform = ( atol(szTmp)>0 );
			return TRUE;
		}
		return FALSE;
	}
	if(strcmp( szParamName,"CCD_ASAS_REVERSE_FOR_TRANSFORM" )==0){
		szTmp = GetParam("CCD_ASAS_REVERSE_FOR_TRANSFORM",TRUE);
		if(szTmp && szTmp[0]){
			m_eReverseForTransform = (eDriverReverseImage_T)atol(szTmp);
			return TRUE;
		}
		return FALSE;
	}
	if(strcmp( szParamName,"CCD_ASAS_BORDER_SIZE" )==0){
		szTmp = GetParam("CCD_ASAS_BORDER_SIZE",TRUE);
		if(szTmp && szTmp[0]){							
			m_nAsasBorderSize = atol( szTmp );
			return TRUE;
		}
		return FALSE;
	}
	if(strcmp( szParamName,"CCD_ASAS_SUBTR_DARK" )==0){
		szTmp = GetParam("CCD_ASAS_SUBTR_DARK",TRUE);
		if(szTmp && szTmp[0]){
			m_bAsasSubtrDark = ( atol(szTmp)>0 );
			return TRUE;
		}
		return FALSE;
	}
	if(strcmp( szParamName,"CCD_ASAS_DIVIDE_BY_FLAT" )==0){
		szTmp = GetParam("CCD_ASAS_DIVIDE_BY_FLAT",TRUE);
		if(szTmp && szTmp[0]){
			m_bAsasDivideByFlat = ( atol(szTmp)>0 );
			return TRUE;
		}
		return FALSE;
	}
	if(strcmp( szParamName,"MAX_AST_ERR" )==0){
		szTmp = GetParam("MAX_AST_ERR",TRUE);
		if(szTmp && szTmp[0]){
			m_fAsasError = atof( szTmp );
			return TRUE;
		}
		return FALSE;
	}
	if(strcmp( szParamName,"MAX_AST_ERR_FATAL" )==0){	
		szTmp = GetParam("MAX_AST_ERR_FATAL",TRUE);
		if(szTmp && szTmp[0]){
		   m_fAsasFatalError = atof( szTmp );
			return TRUE;
		}
		return FALSE;
	}


	if(strcmp( szParamName,"CCD_DRIVER_BINNING" )==0 ){
		szTmp = GetParam("CCD_DRIVER_BINNING",TRUE);
		if(szTmp && szTmp[0]){
			m_DriverAnalBinning = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_PORT_NO" )==0 ){
		szTmp = GetParam("CCD_PORT_NO",TRUE);
		if(szTmp && szTmp[0]){
			m_PortNo = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_CORBA_OPTIONS" )==0 ){
		szTmp = GetParam("CCD_CORBA_OPTIONS",TRUE);
		if(szTmp && szTmp[0]){
			m_CorbaOptions = szTmp;
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_SHARED_MEM_KEY" )==0 ){
		szTmp = GetParam("CCD_SHARED_MEM_KEY",TRUE);
		if(szTmp && szTmp[0]){
			m_SharedMemKey = atol( szTmp );
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_PIPELINE_SAFE_BUFFER_SIZE" )==0 ){
		szTmp = GetParam("CCD_PIPELINE_SAFE_BUFFER_SIZE",TRUE);
		if(szTmp && szTmp[0]){
			m_PipelineSafeBufferSize = atol( szTmp );
			return TRUE;
		}
	}

		// automatic determination of shifts :
	if(strcmp( szParamName,"CCD_SHIFTS_VALUES_FILE" )==0 ){
		szTmp = GetParam("CCD_SHIFTS_VALUES_FILE",TRUE);
		if(szTmp && szTmp[0]){
			m_szShiftsValuesFile = szTmp;
			return TRUE;
		}
	}		

	if(strcmp( szParamName,"CCD_AUTO_SHIFTS_CALC" )==0 ){
		szTmp = GetParam("CCD_AUTO_SHIFTS_CALC",TRUE);
		if(szTmp && szTmp[0]){
			m_nAutoShiftsCalc = atol(szTmp);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_AUTO_SHIFTS_CALC_MATRIX_NO" )==0 ){
		szTmp = GetParam("CCD_AUTO_SHIFTS_CALC_MATRIX_NO",TRUE);
		if(szTmp && szTmp[0]){
			m_nAutoShiftMatrixNo = atol(szTmp);
			return TRUE;
		}
	}

	if(strcmp( szParamName,"CCD_TRANSFORM_MATRIX_FILE" )==0 ){
		szTmp = GetParam("CCD_TRANSFORM_MATRIX_FILE",TRUE);
		if(szTmp && szTmp[0]){
			m_szTransformFile = szTmp;
			if(!m_TransformMatrix.ReadFromFile( m_szTransformFile.c_str() )){
				Assert(FALSE,"Could not read transformation matrix from file : %s\n",m_szTransformFile.c_str());
			}
			return TRUE;
		}
	}

	return FALSE;
}


LONG_T CCDConfig::GetEventsBufferSize()
{
	InitParams();
	return m_EventsBufferSize;
}

void CCDConfig::CheckEdges()
{
	return;
	if(m_FrameDX>0){
		Assert((m_nIgnoreEdgeLeft>=m_FramesBack*m_FrameDX && m_nIgnoreEdgeLeft>=m_nMaxOfAverageOfPrevN*m_FrameDX ),"Left ignore edge to small");
	}	
	if(m_FrameDX<0){
		Assert((m_nIgnoreEdgeRight>=fabs(m_FramesBack*m_FrameDX) && m_nIgnoreEdgeRight>=fabs(m_nMaxOfAverageOfPrevN*m_FrameDX) ),"Right ignore edge to small");
	}	
	if(m_FrameDY>0){
		Assert((m_nIgnoreEdgeBottom>=m_FramesBack*m_FrameDY && m_nIgnoreEdgeBottom>=m_nMaxOfAverageOfPrevN*m_FrameDY ),"Bottom ignore edge to small");
	}	
	if(m_FrameDY<0){
		Assert((m_nIgnoreEdgeUp>=fabs(m_FramesBack*m_FrameDY) && m_nIgnoreEdgeUp>=fabs(m_nMaxOfAverageOfPrevN*m_FrameDY) ),"Right ignore edge to small");
	}	
}

BOOL_T CCDConfig::CheckConsistancy()
{
//	if(m_bCheckForSUPERNEW){
//		Assert( m_bAnalyseMaxFromPrevInPipeline || m_bKeepHomeopaticSum,
//				  "In case of exisiting object flashes search - one of parameters : CCD_CHECK_MAX_FROM_PREV_N or CCD_KEEP_HOMEOPATIC_SUM, must also be set !");
//	}

	CheckEdges();

	if(m_nPipelineSize>MAX_PIPELINE_SIZE){
		Assert(FALSE,"Pipeline size %d exceeds maximum size of %d",m_nPipelineSize,MAX_PIPELINE_SIZE);
	}
	Assert(m_nMaxOfAverageOfPrevN<m_nPipelineSize,"Cannot analyse %d previous frames in pipeline of size %d",m_nMaxOfAverageOfPrevN,m_nPipelineSize);
	Assert(m_FramesBack<m_nPipelineSize,"Cannot analyse %d previous frames in pipeline of size %d",m_FramesBack,m_nPipelineSize);
	Assert(m_nNotMoreThenNExceedsBackFramesCount<m_nPipelineSize,"NOT MORE THEN ... - Cannot analyse %d previous frames in pipeline of size %d",m_nNotMoreThenNExceedsBackFramesCount,m_nPipelineSize);

	if(m_bSuperOptimized){
		if(!m_bAverageOfPrev && !m_bHomeopatic && !m_bCheckGenOnly){
			Assert(FALSE,"In optimized version algorimth must be defined");
		}
		int sum_gt_zero = 0;
		if(m_bAverageOfPrev)
			sum_gt_zero++;
		if(m_bHomeopatic)
			sum_gt_zero++;
		if(m_bCheckGenOnly)
			sum_gt_zero++;
		Assert(sum_gt_zero==1,"For optimized version only one algorithm can be enabled !");
	}
	
	if(m_ConfirmEventsOnNextNFrames>0){
		Assert(m_ConfirmEventsOnNextNFrames<m_nPipelineSize,"Cannot confirm events on more frames then in pipeline");
	}
	
	if(m_bAutoCalcTresh){
		Assert(m_nNewLaplaceInSigma>=0 && m_nMaxLaplaceOnOtherInSigma>=0,"Auto-recalation of tresholds requires that parameters CCD_MAX_LAPLACE_ON_OTHER_IN_SIGMA and CCD_NEW_LAPLACE_TRESHOLD_IN_SIGMA are defined");
		CheckIfBackgroundCalcEnabled();
	}

	/*if(m_bSkipIfMoreThen>0){
		// overlaps rejecting must be turned-off 
		if(m_bSkipOverlaps){
			m_bSkipOverlaps = FALSE;
			MYTRACE3(gCCDTrace,"Changing CCD_SKIP_OVERLAPS=FALSE, due to CCD_SKIP_MORE_THEN>0");
		}		
	}*/

	//if(m_bDarkFrame){
	//	if(strlen(m_szDarkFrameFile.c_str())==0){
	//		Assert(FALSE,"Parameter CCD_DARK_FRAME_FILE must be defined for CCD_DARK_FRAME=1");
	//	}
	//}

	if(m_bUseFrameTimeForShifts){
		if(m_bIgnoreMissingTime){
			Assert(FALSE,"Cannot ignore missing time of frame in case shifts are calculated with time differences");
		}
	}

	if(m_bCheckIfNotEdgeOfBigStar){
		if(!m_CalcBackgrTable[eRawS])
			Assert(FALSE,"When check if not edge of big star raw backgr must be determined ,set CCD_CALC_S_DISTR=1");
	}

	if( m_bMC ){
		m_bParallelMode = FALSE;
		m_bReadAndSaveStatFile = FALSE;
	}

	if( m_nCompareToOldFreqInSec>0 ){
		if ( m_bKeepSumOfPrevNFrames<=0 ){
			Assert(FALSE,"When old frame compare enabled parameter CCD_KEEP_SUM_OF_PREV_N_FRAMES must be >0");
		}
	}

	if( m_bReadFirstLevelInfoFromLog ){
		if( !m_bMC ){
			printf("Reading of events from old log files can only be done in monte carlo mode !!!\n");
			printf("Verify parameters : CCD_MONTE_CARLO / CCD_READ_FIRST_LEVEL_INFO_FROM_LOG\n");
			printf("Cannot continue, exiting now\n");
			exit(0);
		}
		if( m_bDoASASPhotAstr ){
			printf("Astrometry must be disabled when running re-analysis of FLT log files !!!\n");
			printf("Verify parameters : CCD_READ_FIRST_LEVEL_INFO_FROM_LOG / CCD_DO_ASAS_PHOT_ASTR \n");
			printf("Cannot continue, exiting now\n");
         exit(0);
		}
	}

	return TRUE;
}

void CCDConfig::InitBackgrFlagTable()
{
	for(int i=0;i<MAX_LAPLACE_DEFINED;i++){
		m_CalcBackgrTable[i] = FALSE;
	}
}

void CCDConfig::CheckIfBackgroundCalcEnabled()
{
	if(!m_CalcBackgrTable[(int)m_eLaplaceType]){
		Assert( FALSE, "Background calculation of current laplace must be enabled");
	}
}

void CCDConfig::Dump()
{
	if(m_bUseLocalCfgFile){
		m_pCfgFile->Dump();
	}else{
		DumpGlobalParams();	
	}
	printf("m_SumTresholdForNewFrame=%d\n",m_SumTresholdForNewFrame);
	printf("m_SumOnPrevLessThen=%d\n",m_SumOnPrevLessThen);
	printf("m_DiffTreshold=%d\n",m_DiffTreshold);
	printf("m_TresholdPerPixel=%d\n",m_TresholdPerPixel);
	printf("m_MaxPrevPerPixel=%d\n",m_MaxPrevPerPixel);	
}

void CCDConfig::CalcShapes()
{
		m_nVetoPointsCount = CCD_Analyser::CalcShapePoints( m_VetoArea, m_eVetoShape, m_nVetoRedial );

		_TRACE_PRINTF_6("VETO SHAPE :\n");
		for(int i=0;i<m_nVetoPointsCount;i++){
			_TRACE_PRINTF_6("(%d,%d)\n",m_VetoArea[i].x,m_VetoArea[i].y);
		}
		// sleep(5);

		m_nNeighbToSumCount = CCD_Analyser::CalcShapePoints( NULL, m_eNeigbShape, m_dNeighbRedial );
      CCD_Analyser::CalcShapeToConfirm();
}

void CCDConfig::RecalculateParams( CCDMatrix* pMatrix/*=NULL*/, InfoTable2D* pMatrixInfoTable/*=NULL*/ )
{
		if(m_bAlwaysUseSameTreshPerPixel){
			m_SumTresholdForNewFrame = m_nNeighbToSumCount*m_ConfTresholdPerPixel;
			m_SumOnPrevLessThen = m_nNeighbToSumCount*m_ConfMaxPrevPerPixel;
	
			if(m_DiffCheckType==DiffPerPixel)
				m_DiffTreshold = m_nNeighbToSumCount*m_ConfTresholdPerPixel;
			else
				m_DiffTreshold = m_DiffTreshSignal;

			m_TresholdPerPixel = m_ConfTresholdPerPixel;
			m_MaxPrevPerPixel = m_ConfMaxPrevPerPixel;
		}else{
			m_SumTresholdForNewFrame = m_nNeighbToSumCount*m_TresholdPerPixel;
			m_SumOnPrevLessThen = m_nNeighbToSumCount*m_MaxPrevPerPixel;

			if(m_DiffCheckType==DiffPerPixel)
				m_DiffTreshold = m_nNeighbToSumCount*m_TresholdPerPixel;
			else
				m_DiffTreshold = m_DiffTreshSignal;
		}
		m_UseDoubleCheckAboveADU = m_nNeighbToSumCount*m_UseDoubleCheckAboveADUPerPixel;

		// variable objects parameters :
		m_MinPrevTotal = m_nNeighbToSumCount*m_MinPrevPerPixel;
		m_IncreaseTreshADU = m_nNeighbToSumCount*m_IncreaseTreshADUPerPixel;


		m_nVetoPointsCount = CCD_Analyser::CalcShapePoints( m_VetoArea, m_eVetoShape, m_nVetoRedial );

		m_nNotMoreThenNExceedsPointsCount =  CCD_Analyser::CalcShapePoints( mNotMoreThenNExceedsArea, m_eNotMoreThenNExceedsShape, m_nNotMoreThenNExceedsRedial );


		Table2D<ELEM_TYPE>::CalcOverlapParts( m_FrameDX, m_FrameDY, m_S0, m_S1, m_S2, m_S3, m_overlapedPixels );
				

		Table2D<ELEM_TYPE>::GetPlusMinusList( m_eLaplaceType, m_LaplacePlusList, m_LaplacePlusCount, 
											  m_LaplaceMinusList, m_LaplaceMinusCount );

		m_SinOfGeoLatitude = sin( m_GeoLatitude );

		ReCalcTresholds( pMatrix, pMatrixInfoTable );


		if(m_ConfirmEventsOnNextNFrames<0)
			m_ConfirmEventsOnNextNFrames = 0;

		CCD_Analyser::AutoCalculateIgnoreEdges();
}

BOOL_T CCDConfig::GetVisualizeEvents()
{
	InitParams();
	return m_bVisualizeEvents;
}


void CCDConfig::SetOutputDir( const char* szOutDir )
{ 
	m_szOutDir = szOutDir;
}

const char* CCDConfig::GetOutputDir() 
{ 
//	printf("GetOutputDir = %s\n",m_szBaseOutDir.c_str());
	return m_szOutDir.c_str(); 
}


mystring GetEventTypeDesc( eEventType eventType )
{
   mystring szRet = "unkonwn";
   if(eventType == eFlash)
      szRet = "flash";
   if(eventType == eBrighten)
      szRet = "brightening";

   return szRet;
}

long CCDConfig::GetEdgeSizeInLaplace()
{
	return GetEdgeSizeInLaplace( m_eLaplaceType );
}

long CCDConfig::GetEdgeSizeInLaplace( eLaplaceType_T type )
{
	//if( type == eNineEightSevenVeryBig)
	//	return 4;
	return 3;
	if( type==eSinglePoint || type==eFivePoints){
		return 1;
	}
	if( type==eFourPoints || type==eFivePlusFourMin)	
		return 2;
	if( type==eFourPoints )
		return 3;
	return 3;
}

BOOL_T CCDConfig::GetKeepMapFlag()
{
	return (m_bKeepLocalShift || GetCalcBackgrFlag());
}

BOOL_T CCDConfig::GetCalcBackgrFlag()
{
	for(int i=0;i<MAX_LAPLACE_DEFINED;i++){
		if(m_CalcBackgrTable[i]){
			return TRUE;
		}
	}
	return FALSE;
}


void CCDConfig::ReCalcTresholds( CCDMatrix* pMatrix/*=NULL*/, InfoTable2D* pMatrixInfoTable/*=NULL*/ )
{
	if(m_nMaxLaplaceOnOtherInSigma>=0 && (!m_bAutoCalcTresh || m_nInitCount==0)){
		if(m_bSuperOptimized){
			//
			if(m_bAverageOfPrev){
				m_nMaxLaplaceOnOther = (long)(CCDDataResults::GetSigmaBackgroundAverageOfPrevN( m_eLaplaceType, m_nMaxOfAverageOfPrevN, pMatrix, pMatrixInfoTable )*m_nMaxLaplaceOnOtherInSigma);
			}else{
				if(m_bHomeopatic){
					m_nMaxLaplaceOnOther = (long)( CCDDataResults::GetSigmaBackgroundHomeo( m_eLaplaceType, m_HomeopaticFactor, pMatrix, pMatrixInfoTable )*m_nMaxLaplaceOnOtherInSigma );
				}else{
					m_nMaxLaplaceOnOther = (long)(CCDDataResults::GetSigmaBackground(m_eLaplaceType,pMatrix,pMatrixInfoTable)*m_nMaxLaplaceOnOtherInSigma);
				}
			}
		}else{
			if(m_nMaxOfAverageOfPrevN>0){
				m_nMaxLaplaceOnOther = (long)(CCDDataResults::GetSigmaBackgroundAverageOfPrevN( m_eLaplaceType, m_nMaxOfAverageOfPrevN, pMatrix, pMatrixInfoTable )*m_nMaxLaplaceOnOtherInSigma);
			}else{
				// assuming homeopatic in other case :
				m_nMaxLaplaceOnOther = (long)( CCDDataResults::GetSigmaBackgroundHomeo( m_eLaplaceType, m_HomeopaticFactor, pMatrix, pMatrixInfoTable )*m_nMaxLaplaceOnOtherInSigma );
			}
		}
	}
	if(m_nNewLaplaceInSigma>=0 && ( !m_bAutoCalcTresh || m_nInitCount==0) ){
		m_nNewLaplace =(long)(CCDDataResults::GetSigmaBackground(m_eLaplaceType,pMatrix,pMatrixInfoTable)*m_nNewLaplaceInSigma);
	}
}

BOOL_T CCDConfig::DoRejectNotConfirmed()
{
	return (m_ConfirmEventsOnNextNFrames>0 && m_bRejectNotConfirmed);
}


void CCDConfig::SetIgnoreEdges( double left, double right, double bottom, double up )
{
	m_nIgnoreEdgeRight = my_round(right);
   m_nIgnoreEdgeLeft = my_round(left);
   m_nIgnoreEdgeUp = my_round(up);
   m_nIgnoreEdgeBottom = my_round(bottom);	

	SetParam( "CCD_IGNORE_EDGE_RIGHT", right );
	SetParam( "CCD_IGNORE_EDGE_LEFT", left );
	SetParam( "CCD_IGNORE_EDGE_UP", up );
	SetParam( "CCD_IGNORE_EDGE_BOTTOM", bottom );
}

void CCDConfig::InitSatLib()
{
	InitQthFile( m_szQthFile.c_str() );
	CSatInfo::InitSateliteLibrary( m_szQthFile.c_str(), m_szTleFile.c_str(), NULL );
}

void CCDConfig::InitQthFile( const char* qthfile )
{
	if(!MyFile::DoesFileExist( qthfile )){
		MyOFile out( qthfile );
		out.Printf("%s\n",m_szSite.c_str());
		out.Printf(" %.4f\n",AstroAngle::rad2deg(m_GeoLatitude));
		out.Printf(" %.4f\n",AstroAngle::rad2deg(m_GeoLongitude));
		out.Printf(" %d\n",m_GeoAltitude);
	}	
}

double CCDConfig::GetPixScale()
{
	if(m_fPixScale>0)
		return m_fPixScale;
	double pixscale=AstroAngle::rad2arcsec( atan2( m_TransformCCDPixelSize, m_TransformCCDFocus ) );
	return pixscale;
}

double CCDConfig::GetPixToRad( double nPixels )
{
	double pixscale = GetPixScale();
	double in_rad = AstroAngle::arcsec2rad( pixscale*nPixels );
	return in_rad;
}
	
BOOL_T CCDConfig::UpdateParamCommon( const char* szParamName, const char* szValue )
{
	int cam_no=-1;
	int i=0;
	vector<CCDPipeline*>& cam_list = (CCDPipeline::GetPipelineList());
	for( i=0;i<cam_list.size();i++){
		CCDConfig* pCfg = &( (cam_list[i])->m_PipelineCfg );
		if( pCfg == this ){
			cam_no = i;
			break;
		}
	}

	if( gCCDParams.m_bSaveParamChanges &&
		 strcmp( szParamName, "CCD_SAVE_PARAM_CHANGES" )  ){
		// enablig / disabling of saving parameters is not written 
	
		if( cam_no < 0 ){
			CCfgFile dynamic_cfg( DYNAMIC_CFG_FILE );
			dynamic_cfg.SetValue( szParamName, szValue );
			dynamic_cfg.SaveToFile( DYNAMIC_CFG_FILE );
		}else{
			mystring szCfgName;
			szCfgName << "dynamic" << cam_no << ".cfg";
			CCfgFile dynamic_cfg( szCfgName.c_str() );
      	dynamic_cfg.SetValue( szParamName, szValue );
	      dynamic_cfg.SaveToFile( szCfgName.c_str() );
		}
	}

	return TRUE;
}

BOOL_T CCDConfig::UpdateParam( const char* szParamName, const char* szName,
										 double& new_value, BOOL_T& bRetVal )
{
	bRetVal=FALSE;		
	if( strcmp( szParamName, szName )==0 ){
		const char* szTmp = GetParam( szName , TRUE );
		if(szTmp && szTmp[0]){	
			UpdateParamCommon( szParamName, szTmp );
			new_value = atof( szTmp );
			bRetVal = TRUE;

//			gPiLog.TraceFmt("Parameter %s changed new value = %s",szParamName, szTmp );
		}else{
			bRetVal = FALSE;
		}
		return TRUE;
	}
	return FALSE;
}

BOOL_T CCDConfig::UpdateParam( const char* szParamName, const char* szName,
										 int& new_value, BOOL_T& bRetVal )
{
	bRetVal=FALSE;		
	if( strcmp( szParamName, szName )==0 ){
		const char* szTmp = GetParam( szName , TRUE );
		if(szTmp && szTmp[0]){	
			UpdateParamCommon( szParamName, szTmp );
			new_value = atol( szTmp );
			bRetVal = TRUE;
		}else{
			bRetVal = FALSE;
		}
		return TRUE;
	}
	return FALSE;
}

BOOL_T CCDConfig::UpdateParam( const char* szParamName, const char* szName,
										 mystring& new_value, BOOL_T& bRetVal )
{
	bRetVal=FALSE;		
	if( strcmp( szParamName, szName )==0 ){
		const char* szTmp = GetParam( szName , TRUE );
		if(szTmp){	
			UpdateParamCommon( szParamName, szTmp );
			new_value = szTmp;
			bRetVal = TRUE;
		}else{
			bRetVal = FALSE;
		}
		return TRUE;
	}
	return FALSE;
}



BOOL_T CCDConfig::UpdateParam( const char* szParamName, const char* szName,
										 BOOL_T& new_value, BOOL_T& bRetVal )
{
	bRetVal=FALSE;		
	if( strcmp( szParamName, szName )==0 ){
		const char* szTmp = GetParam( szName , TRUE );
		if(szTmp && szTmp[0]){	
			UpdateParamCommon( szParamName, szTmp );
			new_value = ( atol( szTmp )>0 );
			bRetVal = TRUE;
		}else{
			bRetVal = FALSE;
		}
		return TRUE;
	}
	return FALSE;
}


BOOL_T CCDConfig::GetASASTransform( void* TrParam )
{
	BOOL_T bRet=TRUE;
	CCDAsasTransform* pAsasTransform = (CCDAsasTransform*)TrParam;

	const char* szTmp = GetParam("POSANGLE",TRUE);
	if(szTmp && szTmp[0]){
		pAsasTransform->fi = atof(szTmp);
	}else{
		bRet=FALSE;
	}
	szTmp = GetParam("AST_ORD",TRUE);
	if(szTmp && szTmp[0]){
		pAsasTransform->order = atol( szTmp );
	}else{
		bRet=FALSE;
	}

	int xCount=0,yCount=0;
	for(register int i=0;i<ASAS_TRANSFORM_PARAM_COUNT;i++){
		mystring szParamName;
		szParamName << "PAR_X_" << i;
		szTmp = GetParam( szParamName.c_str(), TRUE );
		if(szTmp && szTmp[0]){
			pAsasTransform->px[i] = atof( szTmp );
			xCount++;
		}

		szParamName="PAR_Y_";
		szParamName << i;
		szTmp = GetParam( szParamName.c_str(), TRUE );
		if(szTmp && szTmp[0]){
			pAsasTransform->py[i] = atof( szTmp );				
			yCount++;
		}
	}
	if(xCount<14 || yCount<14){
		bRet = FALSE;
	}

	szTmp = GetParam("PIXSCALE",TRUE);
	if(szTmp && szTmp[0]){
		pAsasTransform->pixscale = atof( szTmp );
	}else{
		bRet = FALSE;
	}
	szTmp = GetParam("RA2000",TRUE);
	if(szTmp && szTmp[0]){
		pAsasTransform->ra = atof( szTmp );
	}else{
      bRet = FALSE;
   }
	szTmp = GetParam("DEC2000",TRUE);
	if(szTmp && szTmp[0]){
		pAsasTransform->dec = atof( szTmp );
	}else{
      bRet = FALSE;
   }
	
	szTmp = GetParam("CCD_TRANSFORM_UT_TIME",TRUE);
	if(szTmp && szTmp[0]){
		pAsasTransform->m_TransformUtTime = atol( szTmp );
	}else{
		bRet = FALSE;
	}

	szTmp = GetParam("AZIM_CENTER",TRUE);
	if(szTmp && szTmp[0]){
		pAsasTransform->m_AzimInRad = AstroAngle::deg2rad( atof(szTmp) );
	}else{
		// bRet = FALSE;
		AstroCCD::calculateHorizontalCoordinatesFromEq( AstroAngle::hours2rad( pAsasTransform->ra ),	
																		AstroAngle::deg2rad( pAsasTransform->dec ),
																		pAsasTransform->m_TransformUtTime,
																		m_GeoLongitude, m_GeoLatitude,
																		pAsasTransform->m_AltInRad,
																		pAsasTransform->m_AzimInRad );
	}
	szTmp = GetParam("ALT_CENTER",TRUE);
	if(szTmp && szTmp[0]){
		pAsasTransform->m_AltInRad = AstroAngle::deg2rad( atof(szTmp) );
	}else{
		// bRet = FALSE;
		AstroCCD::calculateHorizontalCoordinatesFromEq( AstroAngle::hours2rad( pAsasTransform->ra ),	
																		AstroAngle::deg2rad( pAsasTransform->dec ),
																		pAsasTransform->m_TransformUtTime,
																		m_GeoLongitude, m_GeoLatitude,
																		pAsasTransform->m_AltInRad,
																		pAsasTransform->m_AzimInRad );
	}

	szTmp = GetParam("MOUNT_RA",TRUE);
	if(szTmp && szTmp[0]){
		(CCDPipeline::m_WorkingMode).m_MountRA_InRad = atof( szTmp );
	}		
	szTmp = GetParam("MOUNT_DEC",TRUE);
	if(szTmp && szTmp[0]){
		(CCDPipeline::m_WorkingMode).m_MountDec_InRad = atof( szTmp );
	}		
	szTmp = GetParam("MOUNT_AZIM",TRUE);
	if(szTmp && szTmp[0]){
		(CCDPipeline::m_WorkingMode).m_MountAzim_InRad = atof( szTmp );
	}		
	szTmp = GetParam("MOUNT_ALT",TRUE);
	if(szTmp && szTmp[0]){
		(CCDPipeline::m_WorkingMode).m_MountAlt_InRad = atof( szTmp );
	}		
	szTmp = GetParam("MOUNT_OBSMODE",TRUE);
	if(szTmp && szTmp[0]){
		(CCDPipeline::m_WorkingMode).m_MountObsMode = atol( szTmp );
	}		


	szTmp = GetParam("CCD_OBS_MODE",TRUE);
	if( szTmp && szTmp[0] ){
		m_ObsMode = (eObservationMode_T)(atol(szTmp));
		gCCDParams.m_ObsMode = m_ObsMode;
	}

	pAsasTransform->xc = (m_SizeX/2);
	pAsasTransform->yc = (m_SizeY/2);
	pAsasTransform->m_SizeX = m_SizeX;
	pAsasTransform->m_SizeY = m_SizeY;
	
	return bRet;
}

void CCDConfig::SetShiftsToZero()
{
	m_bCorrectForRotation = FALSE;
   m_FrameDX = 0;
   m_FrameDY = 0;
   m_FrameDXPerSec = 0;
   m_FrameDYPerSec = 0;
   m_RotCenterX = 0;
   m_RotCenterY = 0;
   m_RotValueDAlfa = 0;
   m_RotValueDAlfaPerSec = 0;
}

const char* CCDConfig::GetResultFlag( int value )
{
	if( value>0 )
		return "OK";
	return "FAILED";	
}

BOOL_T CCDConfig::InitPiLog( BOOL_T bDoInit )
{
/*	if( bDoInit ){
		CLog4C::m_bEnabled = bDoInit;

		if( !gPiLog.WasInitDone() && CLog4C::m_bEnabled ){
			// Initialize global log , because it was not done before ! :
			gPiLog.Init( "DAQ", TRUE );
		}
	}*/
	return TRUE;
}


// good exit indicator :
class CCDFindExit
{
public :
	CCDFindExit();
	~CCDFindExit();
};


CCDFindExit::CCDFindExit()
{
	// this constructor checks if program can continue
	// if exitOK.txt file exist in current dir it means 
	// that exit was requested and no further analysis 
	// will be performed until this file is removed 

	// [NEW] change - check moved to InitParams function 
 	// CCDPipeline::CheckExitFile();
}

CCDFindExit::~CCDFindExit()
{
	// endOK, exitOK, End
	if( strlen(gCCDParams.GetOutputDir())){
		mystring szFName;
		szFName << gCCDParams.GetOutputDir() << "/exitOK.txt";
		MyOFile out(  szFName.c_str() );
		mystring szOut = get_date_time_string();
		out.Printf("Program exited CORRECTLY at : %s\n",szOut.c_str());
		out.Close();

		// copy ccd.cfg file to output directory - to have it stored :
		// mystring szCfgName = GetGlobalParamFile().GetFileName();
		// mystring cpCmd;
		// cpCmd << "cp " << szCfgName.c_str() << " " << CCDConfig::GetOutputDir() << "/ccd.cfg.sav";
		// system(cpCmd.c_str());
	}
}

static CCDFindExit gExitIndicator;

