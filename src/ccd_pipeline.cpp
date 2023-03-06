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
#include "ccd.h"
#include <mypixellist.h>
#include <cfg.h>
#include <myutil.h>
#include "ccd_pipeline.h"
#include "ccd_trace.h"
#include "ccd_analyse.h"
#include "ccd_runstat.h"
#include "ccd_util.h"
#include <cexcp.h>
#include "ccd_globals.h"
#include <mystrtable.h>
#include <myparser.h>
#include <mykeytab.h>
#include <tab2Ddesc.h>
#include <mymacros.h>
#include <fits_file.h>
#include <math.h>
#include "ccd_datares.h"
#include "ccd_eventfile.h"
#include "ccd_controller.h"
#include "ccd_log.h"
#include "ccd_report.h"
#include <AstroCCD.h>
#include "ccd_image_creator.h"
#include <ccd_fits_header_defs.h>
#include <limits.h>
#include "ccd_asastransform.h"
#include "ccd_errors.h"
#include "ccd_photometry.h"
#include "ccd_file_formats.h"
#include <algorithm> // lower_bound on Ubuntu 

// asas astrometry:
#include <asas_errors.h>


mystring szStartDateTime = get_date_time_string();
BOOL_T gRunNameSaved=FALSE;

int CCDPipeline::m_nAlertsCount=0;

// eDAQNormalMode=0, eDAQSatTriggerMode, eDAQTriggerMovingMode
cDAQStateInfo CCDPipeline::m_WorkingMode;

time_t CCDPipeline::m_NextAstrometryTime=0;
mystring CCDPipeline::m_szPrevObjectValue;
BOOL_T CCDPipeline::m_bMountMove=FALSE;
int CCDPipeline::m_RestartDumpCounter=0;
BOOL_T CCDPipeline::m_bStateRead=FALSE;
int CCDPipeline::m_PipelineCounter=0;
vector<CCDPipeline*> CCDPipeline::m_PipelineList;
BOOL_T CCDPipeline::m_bGlobalOptionsDone=FALSE;
BOOL_T CCDPipeline::m_bPreActionsDoneOK=FALSE;
BOOL_T CCDPipeline::m_bPostInitActionsDoneOK=FALSE;
cDAQStateInfo CCDPipeline::m_SaveStatusInfo;
BOOL_T CCDPipeline::m_bFinalDumped=FALSE;
int CCDPipeline::m_nStartAnalCounter=0;


LONG_T gPipelineNeighbList[MAX_CLUSTER_SIZE];

void CCDPipeline::trace(const char* desc)
{
	printf("Pipeline:%d %s\n",m_PipelineIndex,desc);
}

void CCDPipeline::InitGenObj()
{
	m_pGenObj = new CImageCreator( this );
	m_pGenObj->InitParams();
}

int CCDPipeline::GetFrameCounter1()
{
	if( m_PipelineList.size()>0 ){
		int ret = (m_PipelineList[0])->m_FrameCounter;
		printf("CCDPipeline::GetFrameCounter1 ret=%d\n",ret);
	}
	return 0;	
}

int CCDPipeline::GetFrameCounter2()
{
	if( m_PipelineList.size()>=2 ){
		int ret = (m_PipelineList[1])->m_FrameCounter;
		printf("CCDPipeline::GetFrameCounter2, pipeline_count>=2, ret=%d\n",ret);
	}else{
		if( m_PipelineList.size()>=1 ){
			int ret = (m_PipelineList[0])->m_FrameCounter;
			printf("CCDPipeline::GetFrameCounter2, pipeline_count>=1, ret=%d\n",ret);
		}
	}
	return 0;	

}

void CCDPipeline::LogError( const char* szErrorDesc )
{
	if( strlen( gCCDParams.m_szErrorLogFile.c_str() ) ){
		MyOFile out( gCCDParams.m_szErrorLogFile.c_str(), "a+" );


		mystring szDTM = CMyDate::getdate();
		out.Printf("[ %s ] cam=%d frame=%d %s\n",szDTM.c_str(),m_PipelineIndex,
								m_DayFrameCounter,szErrorDesc);
		printf("[ %s ] cam=%d frame=%d %s\n",szDTM.c_str(),m_PipelineIndex,
                        m_DayFrameCounter,szErrorDesc);fflush(0);

		MyOFile out2( "daq.err" , "a+" );
		out2.Printf("[ %s ] cam=%d frame=%d %s\n",szDTM.c_str(),m_PipelineIndex,
                        m_DayFrameCounter,szErrorDesc);
	}
}

int CCDPipeline::EmergencyExit( int bDoExit )
{
	return CCDPipeline::EmergencyEnd( bDoExit );
}

int CCDPipeline::GetXSize( int ccd_index )
{
	CCDPipeline* pPipeline = GetPipeline( ccd_index );
	if( pPipeline )
		return pPipeline->GetXSize();
	return 0;
}

int CCDPipeline::GetYSize( int ccd_index )
{	
	CCDPipeline* pPipeline = GetPipeline( ccd_index );
	if( pPipeline )
		return pPipeline->GetYSize();
	return 0;
}


CCDPipeline* CCDPipeline::GetPipeline( int pipelineIndex )
{
	vector<CCDPipeline*>::iterator i;
	for(i=m_PipelineList.begin();i!=m_PipelineList.end();i++){
		if((*i)->m_PipelineIndex == pipelineIndex)
			return (*i);
	}
	return NULL;
}


void CCDPipeline::InitAndLockWrkEventList( int size )
{
	if(!m_pWrkEventList){
		m_pWrkEventList = new CCDEventList( MAX(1000,size ));
	}
	m_pWrkEventList->clear();
}


void CCDPipeline::SetParam( const char* cfgcode, const char* cfgval )
{
	gCCDParams.SetParam( cfgcode, cfgval );
	m_bInitialized = FALSE;
	InitParams();
}


CCcdCfg* CCDPipeline::GetCamParamSet( int CamIdx )
{
	if(m_CamCfgTab.size()>CamIdx){
		return ((CCcdCfg*)(m_CamCfgTab[CamIdx]));
	}	
	return NULL;
}

const char* CCDPipeline::GetParam( const char* cfgname )
{
	const char* ret = m_PipelineCfg.GetParam( cfgname, TRUE );
	if(!ret){
		ret = gCCDParams.GetParam( cfgname, TRUE );
	}
	return ret;
}

void CCDPipeline::ReadCamerasParams()
{
	if(m_CamCfgTab.size()==0){
		// initialize cfg-filenames :
		for(register int i=0;i<m_nCCD;i++){
			char cfg_file_name[100];
 			// sprintf(cfg_file_name,"ccd_pipeline%d_cam%d.cfg",m_PipelineIndex,i);
	
			// currently using same file as pipeline :
			sprintf(cfg_file_name,"ccd_pipeline%d.cfg",m_PipelineCfg.m_CameraIndex);
			// sprintf(cfg_file_name,"ccd_pipeline%d.cfg",m_PipelineIndex);

			CCcdCfg* pNewCfg = new CCcdCfg( cfg_file_name );
			m_CamCfgTab.Add( pNewCfg );
		}
	}
	for(register int i=0;i<m_nCCD;i++){
		CCcdCfg* pCfg = (CCcdCfg*)(m_CamCfgTab[i]);
		
		if(!pCfg->IsInitialized()){

			// initializing some parameters with global values - they can be overwritten
			// later - if defined in ccd_pipeline%d_cam%d.cfg file :
			strcpy( (pCfg->m_CCDParams).m_szFramesListFile, gCCDParams.m_szFramesListFile.c_str() );
			(pCfg->m_CCDParams).m_eCAMType = gCCDParams.m_eCAMType;
			(pCfg->m_CCDParams).m_SizeX = gCCDParams.m_SizeX;
			(pCfg->m_CCDParams).m_SizeY = gCCDParams.m_SizeY;
			(pCfg->m_CCDParams).m_nAutoShiftsCalc = gCCDParams.m_nAutoShiftsCalc;
			(pCfg->m_CCDParams).m_FrameDX = gCCDParams.m_FrameDX;
			(pCfg->m_CCDParams).m_FrameDY = gCCDParams.m_FrameDY;
			(pCfg->m_CCDParams).m_FrameDXPerSec = gCCDParams.m_FrameDXPerSec;
			(pCfg->m_CCDParams).m_FrameDYPerSec = gCCDParams.m_FrameDYPerSec;
			(pCfg->m_CCDParams).m_RotCenterX = gCCDParams.m_RotCenterX;
			(pCfg->m_CCDParams).m_RotCenterY = gCCDParams.m_RotCenterY;
			(pCfg->m_CCDParams).m_RotValueDAlfa = gCCDParams.m_RotValueDAlfa;
			(pCfg->m_CCDParams).m_RotValueDAlfaPerSec = gCCDParams.m_RotValueDAlfaPerSec;
			(pCfg->m_CCDParams).m_bUseRotInAverageOfPrevN = gCCDParams.m_bUseRotInAverageOfPrevN;
			(pCfg->m_CCDParams).m_TransformCCDOrientation = gCCDParams.m_TransformCCDOrientation;

			(pCfg->m_CCDParams).m_CCDAxixXOrientation = gCCDParams.m_CCDAxixXOrientation;
			(pCfg->m_CCDParams).m_CCDAxixYOrientation = gCCDParams.m_CCDAxixYOrientation;

			(pCfg->m_CCDParams).m_TransformCCDFocus = gCCDParams.m_TransformCCDFocus;
			(pCfg->m_CCDParams).m_TransformCCDPixelSize = gCCDParams.m_TransformCCDPixelSize;
			(pCfg->m_CCDParams).m_HorAzimuth = gCCDParams.m_HorAzimuth;
			(pCfg->m_CCDParams).m_HorAltitude = gCCDParams.m_HorAltitude; 
			(pCfg->m_CCDParams).m_HorAzimCorr = gCCDParams.m_HorAzimCorr;
			(pCfg->m_CCDParams).m_HorAltCorr = gCCDParams.m_HorAltCorr; 
			(pCfg->m_CCDParams).m_RACorr = gCCDParams.m_RACorr;
			(pCfg->m_CCDParams).m_DecCorr = gCCDParams.m_DecCorr; 
			(pCfg->m_CCDParams).m_ObsMode = gCCDParams.m_ObsMode;
			(pCfg->m_CCDParams).m_DecObs = gCCDParams.m_DecObs;
			(pCfg->m_CCDParams).m_RAObs = gCCDParams.m_RAObs;

			// algo parameters :
			(pCfg->m_CCDParams).m_nMaxLaplaceOnOther = gCCDParams.m_nMaxLaplaceOnOther;
			(pCfg->m_CCDParams).m_nNewLaplace        = gCCDParams.m_nNewLaplace;
			(pCfg->m_CCDParams).m_nNewLaplaceInSigma = gCCDParams.m_nNewLaplaceInSigma;
	      (pCfg->m_CCDParams).m_nMaxLaplaceOnOtherInSigma = gCCDParams.m_nMaxLaplaceOnOtherInSigma;

			(pCfg->m_CCDParams).m_nIgnoreEdge = gCCDParams.m_nIgnoreEdge;
			(pCfg->m_CCDParams).m_nIgnoreEdgeLeft = gCCDParams.m_nIgnoreEdgeLeft;
			(pCfg->m_CCDParams).m_nIgnoreEdgeRight = gCCDParams.m_nIgnoreEdgeRight;
			(pCfg->m_CCDParams).m_nIgnoreEdgeUp = gCCDParams.m_nIgnoreEdgeUp;
			(pCfg->m_CCDParams).m_nIgnoreEdgeBottom = gCCDParams.m_nIgnoreEdgeBottom;


			if(!pCfg->Init(FALSE,FALSE)){
				MYTRACE1(gCCDTrace,"Could not read camera cfg file : " << pCfg->GetName());
			}	
		}
	}	

	InitPipelineParams();	
}

void CCDPipeline::InitSizes( int xSize, int ySize, int nCCD, int size )
{
   if(!xSize){
		if(m_PipelineCfg.m_SizeX>0)
			m_SizeX = m_PipelineCfg.m_SizeX;
		else
	      m_SizeX = gCCDParams.m_SizeX;
	}

   if(!ySize){
		if(m_PipelineCfg.m_SizeY>0)
			m_SizeY = m_PipelineCfg.m_SizeY;
		else
	      m_SizeY = gCCDParams.m_SizeY;
	}

   if(!nCCD)
      m_nCCD = gCCDParams.m_nCamNo;
   if(!size)
      m_nPipelineSize = gCCDParams.m_nPipelineSize;
}

void CCDPipeline::InitPipelineParams()
{
	if( !m_PipelineCfg.IsInitialized() ){
		char cfg_file_name[100];
		if( gCCDParams.m_CameraIndex >= 0 ){
			m_PipelineCfg.m_CameraIndex = gCCDParams.m_CameraIndex;
		}
		

		sprintf(cfg_file_name,"ccd_pipeline%d.cfg",m_PipelineCfg.m_CameraIndex);
		gCCDParams.InitParams();
		// m_PipelineCfg.InitParams();

		// m_PipelineCfg.m_RNoiseRopColDefList = gCCDParams.m_RNoiseRopColDefList;
		if( !m_PipelineCfg.InitLocalCfgFile( cfg_file_name, m_pAsasTransform ) ){
			printf("could not read cfg file : %s\n",cfg_file_name);
			exit(-1);
		}

		CCD_Analyser::AutoCalculateIgnoreEdges( this );

		mystring szCustomFITS;
		szCustomFITS << "header" << m_PipelineIndex << ".cfg";
		customKeyMutex.Lock();
		m_CustomFITSKeys.ReadFromFile( szCustomFITS.c_str() );
		customKeyMutex.UnLock();
	}
}

void CCDPipeline::ExecPreActions()
{
	if(!m_bPreActionsDoneOK){
		// copy ccd.cfg file into OUTPUT DIR :
		mystring szCfgName = GetGlobalParamFile().GetFileName();
      mystring cpCmd;
		mystring szCfgSavDir;
		szCfgSavDir << gCCDParams.GetOutputDir() << "/cfg/";
		MyFile::CreateDir( szCfgSavDir.c_str() );
      cpCmd << "cp " << "*.cfg" << " " << szCfgSavDir.c_str();
      system(cpCmd.c_str());

		if(strlen(gCCDParams.m_szMcInputFile.c_str())){
			cpCmd ="";
			cpCmd << "cp " << gCCDParams.m_szMcInputFile.c_str() << " " << gCCDParams.GetOutputDir() << "mcinput.txt.sav";
			system(cpCmd.c_str());
		}



		m_bPreActionsDoneOK = TRUE;
	}
}


void CCDPipeline::ExecPostInitActions()
{
	if( !m_bPostInitActionsDoneOK ){
		DumpParams();
		m_bPostInitActionsDoneOK = TRUE;
	}
}

void CCDPipeline::DumpParams()
{
	vector<CCDPipeline*>::iterator i;
	mystring szCfg;
   for(i=m_PipelineList.begin();i!=m_PipelineList.end();i++){
		if( ((*i)->m_PipelineCfg).GetParamFilePtr() ){
			mystring szCfg1;
			szCfg1 << gCCDParams.GetOutputDir() << "/cfg/" << "ccd_values" 
					 << (*i)->m_PipelineIndex << "_dump.cfg";
			((*i)->m_PipelineCfg).GetParamFile().SaveToFile( szCfg1.c_str() );
		}
	}
	
	if( gCCDParams.GetParamFilePtr() ){
		szCfg << gCCDParams.GetOutputDir() << "/cfg/" << "ccd_values_dump.cfg";
		gCCDParams.GetParamFile().SaveToFile( szCfg.c_str() );
	}
}


void CCDPipeline::InitParams()
{
	if(!m_bInitialized){
		gCCDParams.InitParams();
	   mystring szOutDir;
	
		if( strlen(gCCDParams.m_szBaseOutSubDir.c_str()) == 0 )
			szOutDir << gCCDParams.m_szBaseOutDir << "/" << szStartDateTime << "/"; 
		else
			szOutDir << gCCDParams.m_szBaseOutDir << "/" << gCCDParams.m_szBaseOutSubDir << "/";


		// m_szOutputDir = szOutDir;
	   gCCDParams.SetOutputDir( szOutDir.c_str() );
		
		if( !gRunNameSaved ){
			MyOFile out("ccd_run_list","a+");
			out.Printf("%s\n",szOutDir.c_str() );
			gRunNameSaved = TRUE;
		}

		InitPipelineParams();

		ExecPreActions();
		m_bInitialized=TRUE;		
	}	
}

void CCDPipeline::RefreshParams()
{
	m_bInitialized = FALSE;
	InitParams();
}

void CCDPipeline::InitAnalyserObj()
{
	if(m_pAnalyser)
		delete m_pAnalyser;
	m_pAnalyser = new CRunStat( this, m_SizeX, m_SizeY );
	m_pAnalyser->Init();	
	m_pAnalyser->SetPipelineObj( this );
}

void CCDPipeline::ResetMCPipeline()
{
	InitAnalyserObj();

	InitPipeline();
	
	if(m_pHomeopaticFrame)
		m_pHomeopaticFrame->SetData(0);

	if(m_pLaplaceFrame)
		m_pLaplaceFrame->SetData(0);

}

void CCDPipeline::ResetPipelineFlags( BOOL_T bResetCounter, BOOL_T bResetList )
{
	my_printf_now("CCDPipeline::ResetPipelineFlags called\n");
	tail = m_nPipelineSize;
   head = m_nPipelineSize;
   iter = -1;

	if( bResetCounter ){
		m_FrameCounter = 0;
		m_StartFrameIndex = 0;
		m_LastFrameFromLogFile = -1;
	}else{
		m_StartFrameIndex = m_FrameCounter;
	}
	if( bResetList ){
		m_FrameIdx = 0;
	}
	m_Count = 0;
	if(m_pHomeopaticFrame){
		m_pHomeopaticFrame->ClearState();
		m_pHomeopaticFrame->InitFrame();
	}
	if(m_pLaplaceFrame){
		m_pLaplaceFrame->ClearState();
		m_pLaplaceFrame->InitFrame();
	}
	// ClearState();
	if(m_pFrameCache){
		m_pFrameCache->Reset();
	}
	m_FirstAnalysiedFrame = -1;

	// [NEW 20040727 - commented : ]
	// m_nAstrometryRetry = 0;

	m_PrevSumOfFramesTime = 0;
	m_nFramesCollectedForSum = 0;
	m_nTotalCollectedSumFrames = 0;
}

void CCDPipeline::ResetPipelineLogs( BOOL_T bClearLog /*=TRUE*/ )
{
	// first dumping events kept in memory :
	//DumpFoundEventsLog( TRUE );
	
	//if( gCCDParams.m_bCCDDouble && m_PipelineList.size()==2 ){
	//	DumpFinalCoicReport(); // safed by flag and lock - only in first call 
	//}

	// m_FramesList.clear(); // cleaning list of collected frames
	m_allFoundEvents.clear();
   m_allGenEvents.clear();
	m_CompiledEventsList.clear();
	m_LastDumpedIndex = -1;
	
	if( bClearLog && !gCCDParams.m_bSameLogFiles ){
		// clean log files , but in case m_bSameLogFiles=TRUE - NO CLEAN -
		// same log files are re-used

		gCCDErrLog.DumpToFile1( this, "WRN_CLEANING_LOGS", "Cleaning event logs" );

		if( gCCDParams.m_bCheckForSUPERNEW ){
			m_SNAllCoicLog.Init();
			m_SNCoicEventsLog.Init();
			m_SNEventsLog.Init();

			m_SNSumAllCoicLog.Init( TRUE );
			m_SNSumCoicEventsLog.Init( TRUE );
			m_SNSumEventsLog.Init( TRUE );
		}

		m_AllEventsLog.Init();
		m_VerifiedEventsLog.Init();
		m_GenEventsLog.Init();

		m_ReGenEventsLog.Init();

		if( gCCDParams.m_bAnalyzeSumOfPrevNFrames ){
			// if analysis of sum of several frames is enabled :
			char szLog[512];
		   sprintf(szLog,AVER_FRAME_EVENTS,m_PipelineIndex);
			CCDEventLog::Init( szLog, TRUE, eVerif );

			sprintf(szLog,AVER_VERIF_EVENTS,m_PipelineIndex);
			CCDEventLog::Init( szLog, TRUE, eVerif );
		}

		/*if( gCCDParams.m_bAnalyzeSumOfPrevNFrames && gCCDParams.m_nCompareToOldFreqInSec>0 ){
			char szLog[512];
			sprintf(szLog,COMPARE_TO_OLD_EVENTS,m_PipelineIndex);
			CCDEventLog::Init( szLog );
		}*/
	}
	// m_FinalEventsLog.Init();

	m_TrackList.clear();

	m_DarkFramesList.clear();

	m_szCurrSumOfNFileName = "";
   m_szPrevSumOfNFileName = "";
}

void CCDPipeline::ResetFinalLog( BOOL_T bCatToALL )
{
	if( bCatToALL ){
		mystring szCmd;
		szCmd << "cat " << gCCDParams.GetOutputDir() << "/" << m_FinalEventsLog.m_szFileName.c_str() 
 			   << " >> " << gCCDParams.GetOutputDir() << "/all_" << m_FinalEventsLog.m_szFileName.c_str();
		system(szCmd.c_str());
		
		m_FinalEventsLog.Init( TRUE , eFinal );
	}
}

void CCDPipeline::InitPipeline( BOOL_T bResetCounter/*=TRUE*/)
{
	// initilizing request handling thread
	InitSystemGlobalOptions( this );
	
	// reseting flags and counters :
	ResetPipelineFlags( bResetCounter );

	// reading parameters specific for each camera :
   ReadCamerasParams();

	// must be after ReadCamerasParams to have local pipeline parameters initialized
	Init_MC();


	// reseting logs :
	ResetPipelineLogs();

	// initializing states of cameras :
	InitCamStates(); 


	InitAsasTransform();

	SaveCurrStatus();
}

void CCDPipeline::RestartPipelineKeepCurrent( BOOL_T bForce )
{
	if( m_bRestartOnNext || bForce ){
		printf("CCDPipeline restart of pipeline called at frame : %d\n",m_DayFrameCounter);

		int curr_pos = tail;
		RestartPipeline( FALSE ); // not cleaning log files 

		if( gCCDParams.m_bReadFirstLevelInfoFromLog ){
			printf("re-analysis of first level trigger log files, restart of pipeline ignored !\n");
			return;
		}

		(GetNext()) = m_pCCD[curr_pos];
		m_FrameUnixTime = 0;

		// ne but do not increment frame counter :
		AcceptNew( CCDPIPELINE_ACCEPT_NEW_AUTOPOST_DEFAULT, FALSE );
		m_bRestartOnNext = FALSE;
	}
}

void CCDPipeline::WaitForCounter( int& counter, int final_value )
{
	int i=0;
	while( counter<final_value ){
		printf_now3("waiting for counter to be %d, now =%d\n",final_value,counter);
		sleep(5);
		i++;
	}
}

void CCDPipeline::RestartPipeline( BOOL_T bClearLog )
{
	printf("restart pipeline called, pipeline : %d\n",m_PipelineIndex);
	gCCDErrLog.DumpToFile1( this, "WRN_RESTART_PIPELINE", "restart of pipeline called" );

	// dumping of logs must be before Flags reset !
	// first dumping events kept in memory :
	DumpFoundEventsLog( TRUE );

	if( gCCDParams.m_bReadFirstLevelInfoFromLog ){
		printf("re-analysis of first level trigger log files, restart of pipeline ignored !\n");


		// CHANGES BELOW ARE NEW 20050422 : 
		// correction due to fact that after restart m_confirmedList was 
		// not cleared and dumped again - causing duplication of final events
		// after restart 
		// cleaning event lists in memory :
		ResetPipelineLogs( bClearLog );
	
		// 
	   m_StarSpyTab.ResetTracedStars();
		m_TrackList.clear();
		m_CurrFrameTracks.clear();
		m_PlaneTracks.clear();
		m_SingleCamTracks.clear();

		// 
		CleanBeforeNewFramePre();

		return;
	}


	// currently commented - DumpFinalCoicReport shoud be exected
	// in CCDDouble 
	/*if( gCCDParams.m_bCCDDouble && m_PipelineList.size()==2 && !m_bFinalDumped ){
		// now must wait until both pipelines already called DumpFoundEventsLog( TRUE );

		// first make sure both are dumped !
		WaitForCounter( m_RestartDumpCounter , m_PipelineList.size() );

		DumpFinalCoicReport(); // safed by flag and lock - only in first call 
	}*/
	
	ResetPipelineFlags( FALSE, FALSE );
	ResetPipelineLogs( bClearLog );
	
	// 
   m_StarSpyTab.ResetTracedStars();
	m_TrackList.clear();
	m_CurrFrameTracks.clear();
	m_PlaneTracks.clear();
	m_SingleCamTracks.clear();
}


void CCDPipeline::InitCamStates()
{
	CCDProcState tmp( this );

	m_CameraStates.clear();
	for(int i=0;i<m_nCCD;i++){
		m_CameraStates.push_back( tmp );		
	}

	for(int i=0;i<m_nCCD;i++){		
		CCcdCfg* pCfg = (CCcdCfg*)(m_CamCfgTab[i]);
		CCDConfig& params = m_PipelineCfg;

		double azim = params.m_HorAzimuth;
		double alt  = params.m_HorAltitude;
		double ra   = params.m_RAObs;
		double dec  = params.m_DecObs;

		/*if(params.m_ObsMode==eNoMovingMode){
			azim += params.m_HorAzimCorr;
			alt  += params.m_HorAltCorr;
		}else{
			ra += params.m_RACorr;
			dec += params.m_DecCorr;
		}*/

		m_CameraStates[i].InitAstroForms( params.m_TransformCCDFocus,
												 params.m_TransformCCDPixelSize,
												 m_SizeX, m_SizeY,
												 gCCDParams.m_GeoLatitude, gCCDParams.m_GeoLongitude,
												 azim, alt,
												 params.m_TransformCCDOrientation, gCCDParams.m_TimeZone, 
												 params.m_CCDAxixXOrientation, params.m_CCDAxixYOrientation,
												 params.m_ObsMode, dec, ra,
												 0, 0, 0, 0, m_pAsasTransform );												 
												 
	}
}

CCDPipeline::CCDPipeline( const CCDPipeline& right )
: m_AllEventsLog(this,"",eAllEventsLog),m_VerifiedEventsLog(this,"",eVerif),
	m_GenEventsLog(this,"",eGenEventsLog),
	m_ReGenEventsLog(this,"",eGenEventsLog),m_FinalEventsLog(this,"",eFinal),
	m_SNEventsLog(this),
   m_SNCoicEventsLog(this),m_SNAllCoicLog(this),m_SNSumEventsLog(this),
   m_SNSumCoicEventsLog(this),m_SNSumAllCoicLog(this)
{
	Assert(FALSE,"Do not use this CCDPipeline::CCDPipeline function - copy constructor");
}

CCDPipeline::CCDPipeline(int size,int nCCD,int xSize,int ySize)
:m_Count(0),m_nPipelineSize(size),m_SizeX(xSize),m_SizeY(ySize),m_nCCD(nCCD),
 m_pDarkFrame(NULL),m_pFlatFrame(NULL),
 m_pPixelMaxFrame(NULL),m_pList(NULL),m_pHomeopaticFrame(NULL),
 m_pFrameCache(NULL),
 m_pMatrixInfoMap(NULL),m_pLaplaceFrame(NULL),
 m_pWrkMatrix(NULL),m_pWrkMatrixBigElem(NULL),m_pWrkMatrixBigElem2(NULL),
 m_pWrkEventList(NULL),m_LastDumpedIndex(-1),
 m_AllEventsLog(this,"",eAllEventsLog),m_VerifiedEventsLog(this,"",eVerif),
 m_GenEventsLog(this,"",eGenEventsLog),
 m_ReGenEventsLog(this,"",eGenEventsLog),
 m_FinalEventsLog(this,"",eFinal),
 m_FirstAnalysiedFrame(-1), m_bInitialized(FALSE),
 m_StarSpyTab("starspy"), m_StarsForMatrix("starsfortransform"),
 m_StartFrameIndex(0), m_bCoordInitialized(FALSE), m_pGenObj(NULL),
 m_pAnalyser(NULL),m_pAsasTransform(NULL), m_bAstroThreadRunning(FALSE),
 m_bForceAstroNow(FALSE),m_DayFrameCounter(0),
 m_bAstrometryRunning(FALSE),m_nAstrometryRetry(0),m_bForceSynchroMode(FALSE),
 m_bAstrometryExecuted(FALSE),m_pPrevSumOfSeveral(NULL),m_pCurrSumOfSeveral(NULL),
 m_pPrevSumOfSeveralLaplace(NULL),m_pCurrSumOfSeveralLaplace(NULL),
 m_FrameUnixTime(0),m_bRestartOnNext(FALSE),m_pOldSumOfSeveral(NULL),
 m_pOldSumOfSeveralLaplace(NULL),m_SNEventsLog(this),m_SNCoicEventsLog(this),
 m_SNAllCoicLog(this),m_SNSumEventsLog(this),m_SNSumCoicEventsLog(this),
 m_SNSumAllCoicLog(this),m_LastWrittenFITSSize(0),m_LastFrameFromLogFile(-1),
 m_bAstrometryIsReady(FALSE), m_pAsasTransformAsynchro(NULL),
 m_nPhotometryStarsCount(0),m_FwhmAver(0),m_SumedFrameIndex(0),m_nGoodAstrometryCount(0),
 m_LastBadFrame(-1),m_LastBadLineY(-1),m_pAsasTransformSav(NULL),
 m_bPrevAstrometryOK(FALSE)

{	
	m_pAsasTransform = new CCDAsasTransform( &m_PipelineCfg, this );
	m_pAsasTransformAsynchro = new CCDAsasTransform( &m_PipelineCfg, this );
	m_pAsasTransformSav = new CCDAsasTransform( NULL, NULL );

	m_PipelineIndex = m_PipelineCounter;
	m_PipelineCfg.m_CameraIndex = m_PipelineIndex;

	m_StarSpyTab.SetPipelinePtr( this );
	m_StarsForMatrix.SetPipelinePtr( this );

	m_PipelineCounter++;
	m_PipelineList.push_back( this );

	// initializing parameters :
	InitParams();

	// sizes are requies by call to  new CRunStat
	InitSizes( xSize, ySize, nCCD, size );


	// initialization of analysing object :
	InitAnalyserObj();	


	InitLogfileName( m_AllEventsLog, (gCCDParams.m_szRunEventsLog).c_str() );
	InitLogfileName( m_VerifiedEventsLog, (gCCDParams.m_szVerifiedEventsLog).c_str() );
	InitLogfileName( m_GenEventsLog, (gCCDParams.m_szGenEventsLog).c_str() );
	InitLogfileName( m_ReGenEventsLog, (gCCDParams.m_szReGenEventsLog).c_str() );
	InitLogfileName( m_FinalEventsLog, FINAL_EVENTS_LOG );

	// SN :
	InitLogfileName( m_SNEventsLog, SN_ALLEVENTS_LOG );
	InitLogfileName( m_SNCoicEventsLog, SN_EVENTS_LOG );
	InitLogfileName( m_SNAllCoicLog, SN_COIC_LOG );

	// SN on sum :
	InitLogfileName( m_SNSumEventsLog, SN_SUM_ALLEVENTS_LOG );
	InitLogfileName( m_SNSumCoicEventsLog, SN_SUM_EVENTS_LOG );
	InitLogfileName( m_SNSumAllCoicLog, SN_SUM_COIC_LOG );
	


	InitPipeline();

	// after paramters initialization update flip info :
	m_pAsasTransform->flip = m_PipelineCfg.m_eReverseForTransform;
	m_pAsasTransformAsynchro->flip = m_PipelineCfg.m_eReverseForTransform;
	m_pAsasTransform->fi = m_PipelineCfg.m_fASASAstrometryFi;
	m_pAsasTransformAsynchro->fi = m_PipelineCfg.m_fASASAstrometryFi;
	
	time_t t1=get_dttm();	
	MYTRACE3(gCCDTrace,"Initializing data pipeline, SizeX=" << m_SizeX
                      << ", SizeY=" << m_SizeY 
                      << ", nCCD=" << m_nCCD << ", nPipelineSize=" << m_nPipelineSize);

	BOOL_T bLaplaceOfCurrent = (gCCDParams.m_bCheckLaplaceCondition || gCCDParams.m_bKeepLaplaceFrame );

//	m_pCCD = new cCCD[m_nPipelineSize](m_SizeX,m_SizeY,m_nCCD,NOT_DEFINED,
//								bLaplaceOfCurrent,this,gCCDParams.m_bMC);		
// version for gcc4.0 - NEW :
	 m_pCCD = new cCCD[m_nPipelineSize];
	for(int ii=0;ii<m_nPipelineSize;ii++){
		m_pCCD[ii].cCCD_InitConstructor( m_SizeX,m_SizeY,m_nCCD,NOT_DEFINED,
                        bLaplaceOfCurrent,this,gCCDParams.m_bMC);
	}



	/*for(int i=0;i<m_nPipelineSize;i++){
		 m_pCCD[i].Init(m_SizeX,m_SizeY,m_nCCD,NOT_DEFINED,bLaplaceOfCurrent);
		 MYTRACE3(gCCDTrace,"CCD-Frame " << i << " allocated");		 
	}*/	
	time_t t2=get_dttm();
	MYTRACE3(gCCDTrace,"Pipeline initialized, time=" << (t2-t1) << " sec");


	BOOL_T bOffsetFrame = gCCDParams.m_bOffsetFrame;
	BOOL_T bPixelMaxFrame = gCCDParams.m_bPixelMaxFrame;

	if(gCCDParams.m_bDarkFrame){
		MYTRACE1(gCCDTrace,"Initializing dark frame");
		t1=get_dttm();	
		m_pDarkFrame = new cCCD( m_SizeX,m_SizeY,m_nCCD);
		t2=get_dttm();

		const char* szDarkFile = gCCDParams.m_szDarkFrameFile.c_str();
		if(strlen( m_PipelineCfg.m_szDarkFrameFile ))
			szDarkFile = m_PipelineCfg.m_szDarkFrameFile.c_str();

		if( MyFile::DoesFileExist( szDarkFile ) ){
			printf("Dark file OK : |%s|\n",szDarkFile);
		}else{
			printf("Dark file NOT FOUND : |%s|\n",szDarkFile);
			mystring szDarkFITC=szDarkFile;
			szDarkFITC << "c";
			if( MyFile::DoesFileExist( szDarkFITC.c_str() ) ){
				printf("Dark FITC file found : |%s|\n",szDarkFITC.c_str());
				printf("USING FITC file\n");

				mystring szOldDark = szDarkFile;
				if( strcmp( szOldDark.c_str() , gCCDParams.m_szDarkFrameFile.c_str() )==0){
					gCCDParams.m_szDarkFrameFile = szDarkFITC.c_str();
				}
				if( strcmp( szOldDark.c_str() , m_PipelineCfg.m_szDarkFrameFile.c_str() )==0 ){ 
					m_PipelineCfg.m_szDarkFrameFile = szDarkFITC.c_str();
				}
				szDarkFile = gCCDParams.m_szDarkFrameFile.c_str();
		      if(strlen( m_PipelineCfg.m_szDarkFrameFile ))
      		   szDarkFile = m_PipelineCfg.m_szDarkFrameFile.c_str();
			}
		}

		printf("reading dark frame : %s (|%s|)\n",szDarkFile,szDarkFile);
		BOOL_T bRet = FALSE;
		if( gCCDParams.m_bReadFirstLevelInfoFromLog ){
			if( gCCDParams.m_bFromLogFile_WithFrames ){
				bRet = m_pDarkFrame->ReadFrame( szDarkFile );
			}else{
				printf("Dark read SKIPPED ! - analysing according to FLT log files !\n");
				bRet = TRUE;
			}
		}else{
			bRet = m_pDarkFrame->ReadFrame( szDarkFile );			
		}
		if(!bRet){
			printf("could not read frame : %s (|%s|)\n",szDarkFile,szDarkFile);
			if( gCCDParams.GetMC() || gCCDParams.m_nDarksToBeTaken<=0 ){
				exit(0);
			}
		}
		((*m_pDarkFrame)[0]).FlipImage( m_PipelineCfg.m_bDriverReverseImage );

		Assert(m_pDarkFrame->GetCount()>0,"Could not initialize dark frame");
		Assert((*m_pDarkFrame)[0].GetXSize()==m_SizeX && (*m_pDarkFrame)[0].GetYSize()==m_SizeY,
             "Dark frame size (%d,%d) does not agree with pipeline frame sizes (%d,%d)",
             (*m_pDarkFrame)[0].GetXSize(),(*m_pDarkFrame)[0].GetYSize(),
             m_SizeX,m_SizeY);

		
		MYTRACE1(gCCDTrace,"Dark frame initialized, time=" << (t2-t1) << " sec");
	}


	if(m_PipelineCfg.m_bFlatFrame || m_PipelineCfg.m_bAsasDivideByFlat){
		MYTRACE1(gCCDTrace,"Initializing flat frame");
		t1=get_dttm();	
		m_pFlatFrame = new Table2D<float>( m_SizeX,m_SizeY );
		t2=get_dttm();
		if(!CCDUtil::ReadFITSFile( *m_pFlatFrame, m_PipelineCfg.m_szFlatFrameFile.c_str() )){
			printf("could not read frame : %s\n",m_PipelineCfg.m_szFlatFrameFile.c_str());
         exit(0);
		}
		m_pFlatFrame->FlipImage( m_PipelineCfg.m_bDriverReverseImage );			

		// Assert(m_pFlatFrame->GetCount()>0,"Could not initialize dark frame");
		Assert((*m_pFlatFrame).GetXSize()==m_SizeX && (*m_pFlatFrame).GetYSize()==m_SizeY,
             "Flat frame size (%d,%d) does not agree with pipeline frame sizes (%d,%d)",
             (*m_pFlatFrame).GetXSize(),(*m_pFlatFrame).GetYSize(),
             m_SizeX,m_SizeY);

		
		MYTRACE1(gCCDTrace,"Flat frame initialized, time=" << (t2-t1) << " sec");
	}

	if( gCCDParams.m_bDoSumOfPrevNFrames || gCCDParams.m_nCompareToOldFreqInSec>0 ){
		m_pPrevSumOfSeveral = new Table2D<BIG_ELEM_TYPE>( m_SizeX, m_SizeY );
		m_pCurrSumOfSeveral = new Table2D<BIG_ELEM_TYPE>( m_SizeX, m_SizeY );
		m_pPrevSumOfSeveralLaplace = new Table2D<BIG_ELEM_TYPE>( m_SizeX, m_SizeY );
		m_pCurrSumOfSeveralLaplace = new Table2D<BIG_ELEM_TYPE>( m_SizeX, m_SizeY );
	}
	if( gCCDParams.m_bOnSumedFrames && gCCDParams.m_bAutoCalcSum ){
		m_pCurrSumOfSeveral = new Table2D<BIG_ELEM_TYPE>( m_SizeX, m_SizeY );
	}
	if( gCCDParams.m_nCompareToOldFreqInSec>0 ){
		m_pOldSumOfSeveral = new Table2D<BIG_ELEM_TYPE>( m_SizeX, m_SizeY );
		m_pOldSumOfSeveralLaplace = new Table2D<BIG_ELEM_TYPE>( m_SizeX, m_SizeY );
	}

	if(bPixelMaxFrame){
		MYTRACE1(gCCDTrace,"Initializing pixel maximums frame");
		t1=get_dttm();	
		m_pPixelMaxFrame = new cCCD( m_SizeX,m_SizeY,m_nCCD);
		m_pPixelMaxFrame->SetData(0);
		t2=get_dttm();
		MYTRACE1(gCCDTrace,"Pixel maximums frame initialized, time=" << (t2-t1) << " sec");
	}	

	if(gCCDParams.m_bKeepHomeopaticSum){
		MYTRACE1(gCCDTrace,"Initializing homeopatic frame");
		t1=get_dttm();
		m_pHomeopaticFrame = new CManyTab2D<BIG_ELEM_TYPE>( m_SizeX,m_SizeY,m_nCCD);
		m_pHomeopaticFrame->SetData(0);							
		t2=get_dttm();
		MYTRACE1(gCCDTrace,"Homeopatic frame initialized, time=" << (t2-t1) << " sec");		
	}
	
	if(gCCDParams.m_bKeepLaplaceFrame){			
		MYTRACE1(gCCDTrace,"Initializing homeopatic frame");
		t1=get_dttm();
		m_pLaplaceFrame = new CManyTab2D<BIG_ELEM_TYPE>( m_SizeX,m_SizeY,m_nCCD);
		m_pLaplaceFrame->SetData(0);
		t2=get_dttm();
		MYTRACE1(gCCDTrace,"Laplace frame initialized, time=" << (t2-t1) << " sec");				
	}

	if(gCCDParams.GetKeepMapFlag()){
		// initializing m_nCCD background maps
/*	   m_pMatrixInfoMap = new InfoTable2D[m_nCCD]( m_SizeX, m_SizeY,
                                          gCCDParams.m_BackgrMapXSize,
                                          gCCDParams.m_BackgrMapYSize, 
														FALSE, gCCDParams.m_nIgnoreEdge );*/
// version for gcc4.0 :
		m_pMatrixInfoMap = new InfoTable2D[m_nCCD];
		for(int ii=0;ii<m_nCCD;ii++){
			m_pMatrixInfoMap[ii].InfoTable2D_InitConstructor( m_SizeX, m_SizeY,
                                          gCCDParams.m_BackgrMapXSize,
                                         gCCDParams.m_BackgrMapYSize,
                                        FALSE, gCCDParams.m_nIgnoreEdge );
		}

		if(gCCDParams.m_bKeepLocalShift){
			InitShiftsInfo( m_pMatrixInfoMap[0] );
		}
	}

	ClearState();

	InitWorkingFrames();	


	BOOL_T bReadFromFile=FALSE;
	if(m_PipelineCfg.m_nAutoShiftsCalc>0){
		// auto calculation of shifts on first frames :
		mystring szShiftFileName;
		GetShiftFileName( szShiftFileName );

		if(m_StarSpyTab.ReadFromFile( szShiftFileName.c_str() )){
			// in case found in file - skip determining on first gCCDParams.m_nAutoShiftsCalc frames 
			printf("\nPipeline : %d , re-using shift data from file : %s\n",m_PipelineIndex,szShiftFileName.c_str() );
			

			m_PipelineCfg.m_nAutoShiftsCalc = -abs( m_PipelineCfg.m_nAutoShiftsCalc );					
			gCCDParams.m_nAutoShiftsCalc = -abs( gCCDParams.m_nAutoShiftsCalc );
			gCCDParams.SetParam("CCD_AUTO_SHIFTS_CALC",-abs( gCCDParams.m_nAutoShiftsCalc ));

			SetShiftParams();
			PrintShiftParams();
			bReadFromFile=TRUE;
		}
		
		if( !bReadFromFile ){
			printf("Pipeline : %d, determining shifts in first %d steps\n",m_PipelineIndex,m_PipelineCfg.m_nAutoShiftsCalc);
		}
	}
	if( !bReadFromFile ){
		// printf("Pipeline : %d, determining shifts in first %d steps\n",m_PipelineIndex,m_PipelineCfg.m_nAutoShiftsCalc);
		// printf("Initialzing m_StarSpyTab ...");
		for(int i=0;i<m_nCCD;i++){
			CCDStarSpy tmp(m_SizeX,m_SizeY,gCCDParams.m_nAutoShiftMatrixNo,this);				
			m_StarSpyTab.push_back( tmp );				
		}
		// printf("OK\n");
	}

	if(gCCDParams.m_bTraceOnAllFrames){
		for(int i=0;i<m_nCCD;i++){
			CCDStarSpy tmp(m_SizeX,m_SizeY,gCCDParams.m_nAutoShiftMatrixNo);
			m_StarsForMatrix.push_back( tmp );
		}
	}

	// controll :
	m_pControlSpyStar = new CCDStarSpy( m_SizeX, m_SizeY );

	InitAsasTransform();

	SaveCurrStatus();
}

void CCDPipeline::InitAsasTransform()
{
	if( gCCDParams.m_bDoASASPhotAstr || gCCDParams.m_bUseAsasTransform ){
		if( !m_pAsasTransform->m_bReadDone ){
			m_pAsasTransform->m_bTransformOK = m_PipelineCfg.GetASASTransform( m_pAsasTransform );
			m_pAsasTransform->m_bReadDone = TRUE;
		}
		if( m_pAsasTransform->m_bTransformOK ){
			UpdateObsCoordinates( AstroAngle::hours2rad( m_pAsasTransform->ra ),
										 AstroAngle::deg2rad( m_pAsasTransform->dec ),
										 m_pAsasTransform->m_AzimInRad,
										 m_pAsasTransform->m_AltInRad );
		}
	}
}


void CCDPipeline::MoveCurrentShiftsFile()
{
	mystring szOldFileName,szSav;
	GetShiftFileName( szOldFileName );
	szSav << szOldFileName << ".sav#" << m_FrameCounter;


	mystring szSaveFile;
	
	GetOutFileName( szSaveFile, "AUTO_SHIFTS", szSav.c_str(), this, -1, FALSE );
	MyFile::CreateDir( szSaveFile.c_str() );
	mystring szCmd;
	szCmd << "mv " << szOldFileName << " " << szSaveFile;
	system( szCmd.c_str() );
}

void CCDPipeline::GetShiftFileName( mystring& out )
{
	out = "";
	out = gCCDParams.m_szShiftsValuesFile;
	if(strstr(out.c_str(),"%d")){
		char tmp[256];
		sprintf(tmp,out.c_str(),m_PipelineIndex);
		out = tmp;
	}
}

void CCDPipeline::InitLogfileName( CCDEventLog& eventLog, const char* szFileName )
{
	mystring szFile = szFileName;
	if(strstr(szFileName,"%d")){
		char* ptmp = new char[strlen(szFileName)+100];
		sprintf(ptmp,szFileName,m_PipelineIndex);
		szFile = ptmp;
		delete [] ptmp;
	}
	eventLog.SetFileName( szFile.c_str() );
}

void CCDPipeline::InitWorkingFrames()
{
	if( gCCDParams.m_bKeepLaplaceFrame || gCCDParams.m_bKeepHomeopaticSum ){
		if(!m_pWrkMatrix)
			m_pWrkMatrix = new CCDMatrix( m_SizeX, m_SizeY );
		if(!m_pWrkMatrixBigElem)
			m_pWrkMatrixBigElem = new Table2D<BIG_ELEM_TYPE>( m_SizeX, m_SizeY );
		if(gCCDParams.m_bCalcMaxNeighbHomeo){
			if(!m_pWrkMatrixBigElem2)
				m_pWrkMatrixBigElem2 = new Table2D<BIG_ELEM_TYPE>( m_SizeX, m_SizeY );
		}
	}
}

void CCDPipeline::CleanWorkingFrames()
{
	if(m_pWrkMatrix){
		delete m_pWrkMatrix;	
		m_pWrkMatrix = NULL;
	}
	if(m_pWrkMatrixBigElem){
		delete m_pWrkMatrixBigElem;
		m_pWrkMatrixBigElem = NULL;
	}
	if(m_pWrkMatrixBigElem2){
		delete m_pWrkMatrixBigElem2;
		m_pWrkMatrixBigElem2 = NULL;
	}
}

void CCDPipeline::LockWorkingFrames(){}

void CCDPipeline::UnLockWorkingFrames(){}


void CCDPipeline::End()
{
	DumpFoundEventsLog( TRUE );
}

void CCDPipeline::StopAnal()
{
   End();
}


void CCDPipeline::Clean()
{
//	if(m_pCCDDeviceInterface){
//		delete m_pCCDDeviceInterface;
//		m_pCCDDeviceInterface = NULL;
//	}


	DumpAllEvents();

	if(m_pAnalyser){
		delete m_pAnalyser;
		m_pAnalyser = NULL;
	}
	
	if(m_pCCD){
		delete [] m_pCCD;
		m_pCCD = NULL;
	}
	if(m_pFlatFrame){
		delete m_pFlatFrame;
		m_pFlatFrame = NULL;
	}
	if(m_pDarkFrame){
		delete m_pDarkFrame;
		m_pDarkFrame = NULL;
	}
	if(m_pPixelMaxFrame){
		delete m_pPixelMaxFrame;
		m_pPixelMaxFrame = NULL;
	}
	if(m_pHomeopaticFrame){
		delete m_pHomeopaticFrame;
		m_pHomeopaticFrame = NULL;
	}
	if(m_pMatrixInfoMap){
      delete [] m_pMatrixInfoMap;
		m_pMatrixInfoMap = NULL;
	}
	if(m_pList){
		delete m_pList;
		m_pList = NULL;
	}
	if(m_pLaplaceFrame){
		delete m_pLaplaceFrame;
		m_pLaplaceFrame = NULL;
	}
	if(m_pWrkEventList){
		delete m_pWrkEventList;
		m_pWrkEventList = NULL;
	}
	CleanWorkingFrames();

	if(m_pGenObj){
		delete m_pGenObj;
		m_pGenObj = NULL;
	}

	if(m_pControlSpyStar){
		delete m_pControlSpyStar;
		m_pControlSpyStar = NULL;
	}

	if(m_pAsasTransform){
		delete m_pAsasTransform;
	}

	if( m_pAsasTransformAsynchro ){
		delete m_pAsasTransformAsynchro;
	}

	if( m_pAsasTransformSav ){
		delete m_pAsasTransformSav;
	}

	if( m_pPrevSumOfSeveral ){
		delete m_pPrevSumOfSeveral;
	}
	if( m_pCurrSumOfSeveral ){
		delete m_pCurrSumOfSeveral;
	}
	if( m_pPrevSumOfSeveralLaplace ){
		delete m_pPrevSumOfSeveralLaplace;
	}
	if( m_pCurrSumOfSeveralLaplace ){
		delete m_pCurrSumOfSeveralLaplace;
	}
	if( m_pOldSumOfSeveral ){
		delete m_pOldSumOfSeveral;
	}
	if( m_pOldSumOfSeveralLaplace ){
		delete m_pOldSumOfSeveralLaplace;
	}
}

CCDPipeline::~CCDPipeline(){	
	m_PipelineCounter--;
	if( m_PipelineCounter<0 )
		m_PipelineCounter = 0;
	_TRACE_PRINTF_3("~CCDPipeline() ????\n");fflush(0);
	Clean();
}

CCDEventList& CCDPipeline::GetFoundEvents(int idx/*=0*/){
	/*if(gCCDParams.m_ConfirmEventsOnNextNFrames){
		cCCD* pFrame = GetVerifiedFrame();
		if(pFrame){
			return (*pFrame)[idx].GetFoundEvents();
		}else{
			return (CCDEventList&)m_EmptyEventsList;
		}
	}*/
	return (GetCurrent()[idx]).GetFoundEvents();
}	


/*CCDEventList& CCDPipeline::GetCurrentEvents()
{
	return m_pAnalyser->GetNewEvents();
}*/

CCDEventList& CCDPipeline::GetNewEvents(int idx/*=0*/)
{ 
	/*if(gCCDParams.m_ConfirmEventsOnNextNFrames){
		cCCD* pFrame = GetVerifiedFrame();
		m_ConfirmedEvents.clear();
		if(pFrame){
			m_ConfirmedEvents = (*pFrame)[idx].GetFoundEvents();
			m_ConfirmedEvents += (*pFrame)[idx].GetGenEvents();
		}
		return m_ConfirmedEvents;
	}*/
	return m_pAnalyser->GetNewEvents();
}


cCCD* CCDPipeline::GetVerifiedFrame()
{
		// confirmation of events on next frames is required -
		// no dumping of new events - checking if can confirm 
		// previous frames, analysing and dumping results :
		if(gCCDParams.m_ConfirmEventsOnNextNFrames>=GetPipelineSize()){
			printf("Pipeline to small to check %d next frames and confirm events, exiting ...\n",gCCDParams.m_ConfirmEventsOnNextNFrames);
			exit(0);
		}			
		if(GetFrameIndex()>=(gCCDParams.m_nPipelineSize + gCCDParams.m_ConfirmEventsOnNextNFrames)){
			// if already possible to analyse 
			LONG_T FrameToVerify = GetFrameIndex()-gCCDParams.m_ConfirmEventsOnNextNFrames;
			LONG_T pos;
         cCCD* pFrame = FindFrame( FrameToVerify, pos );
			if(pFrame){
				return pFrame;				
			}
		}
		return NULL;
}

void CCDPipeline::Add(cCCD& cNew)
{
	m_FrameCounter++;
	long new_pos = tail-1;
	BOOL_T bFirstRound=(m_Count<m_nPipelineSize);
	if(new_pos<0){
		new_pos = m_nPipelineSize-1;
	}
	m_pCCD[new_pos].ClearState();
	m_pCCD[new_pos] = cNew;
	m_pCCD[new_pos].SetFrameIndex( m_FrameCounter );
	tail = new_pos;
	if (bFirstRound){
		head = m_nPipelineSize-1;
	}else{
		head = new_pos-1;
		if(head<0)
			head = m_nPipelineSize-1;
	}
	if(m_Count<m_nPipelineSize)
		m_Count++;
	ClearEvents( tail );
}

cCCD* CCDPipeline::FindFrame( LONG_T FrameIdx, LONG_T& pos )
{
	for(int i=0;i<m_nPipelineSize;i++){
		if(m_pCCD[i].GetFrameIndex()==FrameIdx){
			pos = i;
			return (&(m_pCCD[i]));
		}
	}
	return NULL;
}

long CCDPipeline::GetCount()
{
	return m_Count;
}

cCCD* CCDPipeline::GetCurr()
{
	if(m_Count>0){
		iter = ((tail+1) % m_Count);
		return &(m_pCCD[tail]);
	}
	return NULL;
}

cCCD& CCDPipeline::GetCurrent()
{
	Assert(m_Count>0 ,"No image in pipeline");	
	iter = ((tail+1) % m_Count);
	return m_pCCD[tail];
}

cCCD* CCDPipeline::GetPrevFrame()
{
	if( m_Count>=2 ){
		int prev = ((tail+1) % m_nPipelineSize);
		return &(m_pCCD[prev]);		
	}	
	return NULL;
}

cCCD* CCDPipeline::GetHead()
{
	Assert(m_Count>0 ,"No image in pipeline");
	return &(m_pCCD[head]);
}

cCCD* CCDPipeline::GetPrev()
{
	Assert(m_Count==m_nPipelineSize,"Do not use iterator before pipeline is full !");
	if(iter!=head){
		cCCD* ret = &(m_pCCD[iter]);
		iter++;
		if(iter==m_Count)
			iter = 0;
		return ret;
	}else{
		return NULL;
	}
}

LONG_T CCDPipeline::begin()
{
		Assert(m_Count>0 ,"No image in pipeline");
      LONG_T beg = ((tail+1) % m_Count);
      return beg;	
}


LONG_T CCDPipeline::end()
{
	return head;
}


cCCD* CCDPipeline::GetFrameObject(LONG_T idx)
{
	Assert(idx>=0 && idx<m_nPipelineSize,"Index of frame %d out of range",idx);
	return &(m_pCCD[idx]);
}

cCCD* CCDPipeline::GetFrame(LONG_T idx)
{
	Assert(idx>=0 && idx<m_Count,"Index of frame %d out of range",idx); 
	return &(m_pCCD[idx]);
}

cCCD& CCDPipeline::GetNext()
{
	long new_pos = tail-1;
	if(new_pos<0){
		new_pos = m_nPipelineSize-1;
	}
	m_pCCD[new_pos].ClearState();
	return m_pCCD[new_pos];
}

cCCD& CCDPipeline::AcceptNew(BOOL_T bAutoPost,int bIncCounter /*=TRUE*/ )
{
	if( bIncCounter ){
		if( !gCCDParams.m_bReadFirstLevelInfoFromLog ){
			m_FrameCounter++;
		}
	}
	long new_pos = tail-1;
	BOOL_T bFirstRound=(m_Count<m_nPipelineSize);
	if(new_pos<0){
		new_pos = m_nPipelineSize-1;
		bFirstRound = FALSE;
	}
	tail = new_pos;
	m_pCCD[tail].SetFrameIndex( m_FrameCounter );
	if (bFirstRound){
		head = m_nPipelineSize-1;
	}else{
		head = new_pos-1;
		if(head<0)
			head = m_nPipelineSize-1;
	}
	if(m_Count<m_nPipelineSize)
      m_Count++;

	if(bAutoPost)
		PostNewFrame();
	ClearEvents( tail );
	m_pCCD[new_pos].ClearState();

	time_t newFrameUnixTime = 0;

	if( gCCDParams.m_bReadFirstLevelInfoFromLog ){
		if( gCCDParams.m_bFromLogFile_WithFrames ){
			newFrameUnixTime = (time_t)(GetCurrent()[0].getObsTime());
		}else{
			newFrameUnixTime = m_FrameUnixTime_FromLogFile;
		}
	}else{
		newFrameUnixTime = (time_t)(GetCurrent()[0].getObsTime());
	}

	if( m_FrameCounter>1 && m_FrameUnixTime>0 && 
		( newFrameUnixTime - m_FrameUnixTime ) > gCCDParams.m_nPipelineRestartTimeout ){
		printf("New frame came after : %d sec ( timeout = %d sec ) , restart of pipeline required !\n",( newFrameUnixTime - m_FrameUnixTime ),gCCDParams.m_nPipelineRestartTimeout);
		m_bRestartOnNext = TRUE;
	}
	m_FrameUnixTime = newFrameUnixTime;
	return m_pCCD[new_pos];
}

void CCDPipeline::PostNewFrame()
{
	if(m_pPixelMaxFrame){
		// UpdateMaximumsFrame();
      UpdateMaximumsFrameOpt();
   }
	if(m_pHomeopaticFrame){		
		UpdateHomeopaticFrame();
	}
	if(m_pLaplaceFrame){
		UpdateLaplaceFrame();
	}
}

void CCDPipeline::PreAnalysis()
{
	if(gCCDParams.m_bAutoCalcTresh){
      AutoCalcTresholds();
   }
}

void CCDPipeline::AutoCalcTresholds()
{
	BOOL_T bCallRefresh = FALSE;
	if(gCCDParams.m_nNewLaplaceInSigma>=0){
		double SigmaB;
		SigmaB = gCCDParams.m_nNewLaplaceInSigma*m_pMatrixInfoMap->GetSigmaBackground( 0, 0, gCCDParams.m_eLaplaceType );
		
		mystring szSigma;	
		szSigma = mystring::double2string( SigmaB );
		gCCDParams.SetParam( "CCD_NEW_LAPLACE_TRESHOLD", szSigma.c_str() );
		bCallRefresh = TRUE;
			
		_TRACE_PRINTF_4("CCD_NEW_LAPLACE_TRESHOLD := %s\n",szSigma.c_str() );

      m_PipelineCfg.m_nNewLaplace = (int)SigmaB;
   }
   if(gCCDParams.m_nMaxLaplaceOnOtherInSigma>=0){
		double SigmaB;
      SigmaB = gCCDParams.m_nMaxLaplaceOnOtherInSigma*m_pMatrixInfoMap->GetSigmaBackground( 0, 0, gCCDParams.m_eLaplaceType );
		// SigmaB = CCDDataResults::GetSigmaBackgroundAverageOfPrevN( gCCDParams.m_eLaplaceType, gCCDParams.m_nMaxOfAverageOfPrevN );
		

		mystring szSigma;
      szSigma = mystring::double2string( SigmaB );
      gCCDParams.SetParam( "CCD_MAX_LAPLACE_ON_OTHER", szSigma.c_str() );
		bCallRefresh = TRUE;		

		_TRACE_PRINTF_4("CCD_MAX_LAPLACE_ON_OTHER = %s\n",szSigma.c_str() );

		m_PipelineCfg.m_nMaxLaplaceOnOther = (int)SigmaB;
   }
	if(bCallRefresh){
		gCCDParams.InitParam("CCD_NEW_LAPLACE_TRESHOLD");
		gCCDParams.InitParam("CCD_MAX_LAPLACE_ON_OTHER");
	}
}

void CCDPipeline::LogSamplesInfo()
{
	if( gCCDParams.GetMC() && gCCDParams.m_bLogSamplePut ){
		if( gCCDParams.m_nSamplesToPutOnFrame<=1 ){
			if( m_pGenObj && m_pGenObj->m_LastPutObjectX>0 && m_pGenObj->m_LastPutObjectY>0 ){
				int max_val=-100000,max_lap=-1000000;
				CCDMatrix& frame = GetCurrent()[0];
				ELEM_TYPE** p_data = frame.get_data_buffer_fast();
				BIG_ELEM_TYPE** p_lap = frame.get_laplace_data_fast();
				int low_x = MAX(m_pGenObj->m_LastPutObjectX-5,0);
				int low_y = MAX(m_pGenObj->m_LastPutObjectY-5,0);

				int up_x = MIN(m_pGenObj->m_LastPutObjectX+5,(m_SizeX-1));
	         int up_y = MIN(m_pGenObj->m_LastPutObjectY+5,(m_SizeY-1));

				if( p_data ){			
					for(register int y=low_y;y<=up_y;y++){
						for(register int x=low_x;x<=up_x;x++){
							if( p_data[y][x]>max_val ){
								max_val = p_data[y][x];
							}						
						}
					}
				}
				if( p_lap ){			
					for(register int y=low_y;y<=up_y;y++){
						for(register int x=low_x;x<=up_x;x++){
							if( p_lap[y][x]>max_lap ){
									max_lap = p_lap[y][x];
							}						
						}
					}
				}

				mystring szRA = AstroAngle::toString( m_pGenObj->m_LastPutObjectRA, ANGLE_RA_TYPE ).c_str();
				double dec = AstroAngle::rad2deg( m_pGenObj->m_LastPutObjectDEC );

				CCDLog samplesLog( "%d\t%d\t%s\t%.4f\t%d\t%d\t%d\n","x y RA DEC Laplace LapPrev RawValue\n","Samples", "samples");
         	samplesLog.DumpToFile1( this,m_pGenObj->m_LastPutObjectX,	
											m_pGenObj->m_LastPutObjectY,
											szRA.c_str(),dec,
										 	max_lap,
											(m_pGenObj->m_pGenEvent)->m_PixelAnalResults.maxAverageOfPrev,
										   max_val);
			}
		}else{
			// when multiple samples are put on frame 
			if( m_allGenEvents.size() && m_allGenEvents.back().size() ){
				for(int i=0;i<(m_allGenEvents.back()[0]).size();i++){
					CccdReport& evt = (m_allGenEvents.back())[0][i];

					mystring szRA = AstroAngle::toString( m_pGenObj->m_LastPutObjectRA, ANGLE_RA_TYPE ).c_str();
					double dec = AstroAngle::rad2deg( m_pGenObj->m_LastPutObjectDEC );

					CCDLog samplesLog( "%d\t%d\t%s\t%.4f\t%d\t%d\t%d\n","x y RA DEC Laplace LapPrev RawValue\n","Samples", "samples");
					samplesLog.DumpToFile1( this, (int)evt.m_MaxPoint.x,	
										(int)evt.m_MaxPoint.y,
										szRA.c_str(),dec,
									 	(int)evt.m_PixelAnalResults.laplaceSum,
										evt.m_PixelAnalResults.maxAverageOfPrev,
										(int)evt.m_PixelAnalResults.PixelRawValue);
				}
			}
		}
	}
}

BOOL_T CCDPipeline::AnalyzeNewFrame(BOOL_T bAutoDump,BOOL_T bReport,LONG_T idx,BOOL_T bAutoSaveEvents)
{
	if(!idx)
		idx = m_FrameCounter;

/*	lines commented - currently - after each step rotation is determined 
   more and more precisly - and used - so ne skiping due to tracing star :
	if(!gCCDParams.m_bShiftUsesAstroFormulas){
		// if not using astro forumals frst frames are skiped to
		// calculate rotation values :

		if(m_PipelineCfg.m_nAutoShiftsCalc>0){
			if(m_FrameCounter<m_PipelineCfg.m_nAutoShiftsCalc){
				// skiping frame before deterimnation of shifts :

				return FALSE;
			}
		}
	}*/

	if(m_WorkingMode.m_WorkingMode!=eDAQNormalMode){
		// specific ations for special modes (Satelite-Trigger or moving )


		// no saving of frames on trigger now :
		/*if(m_WorkingMode.m_WorkingMode == eDAQSatTriggerMode){
			// satelite trigger mode - dumping all frames :
			// but performing analysis also :
			SaveCurrentFrame();
		}*/

	}
	
	if(GetCount()==m_nPipelineSize){
		BOOL_T bRet = FALSE;		


		
		// check events - verify next frame values :
		SetNextFrameValues();

		// confirmation of previous events :
		if(gCCDParams.m_ConfirmEventsOnNextNFrames){
			ConfirmEventsFromPreviousFrames();
		}
	
		if(gCCDParams.GetMC()){
			// adding generated events to list :
         AddGenEventsToList();
		}
		PrepareFoundEventsList();



		PROFILER_START
		// some pre-analysis actions - line recalulating automatic tresholds :
		PreAnalysis();


//		if( gCCDParams.GetMC() ){
//			LogSamplesInfo();
//		}

		if(!gCCDParams.m_bDoNotAnalyse){
			if(m_FirstAnalysiedFrame<0)
				m_FirstAnalysiedFrame = m_FrameCounter;
			bRet = m_pAnalyser->AnalyseNewFrameOpt( *this, bReport, idx );		
		}else{
			// instead of adding found events to list we only add empty list of 
			// events to have consistency :
			 cCCD& currFrame = GetCurrent();

			if( m_allFoundEvents.back().size()==0 ){
			   for(register int i=0;i<currFrame.GetCount();i++){
		   	   m_allFoundEvents.back().push_back( m_EmptyEventsList );
   			}
			}
		}

		// logging samples info after it is filled in _pAnalyser->AnalyseNewFrameOpt
		if( gCCDParams.GetMC() ){
			LogSamplesInfo();
		}


		PROFILER_END("Analysis of new frame took : ");


		// add found events to global list :
		// AddFoundEventsToList();

// Anti-cloud check - verify number of stars on image and reject all if to
// small number of stars :
		if( gCCDParams.m_MinStarsToAccEvents > 0 ){
			int star_count = m_nPhotometryStarsCount;
			if( star_count <= 0 && m_CurrStarList.size()>star_count ){
				star_count = m_CurrStarList.size();
			}
			if( star_count <=0 ){
				CCD_Analyser* pAnal1 = GetAnalPtr(); 
				if( pAnal1 && pAnal1->backgrStat.size()>=1 ){
					if( (pAnal1->backgrStat).back().nTnewCut > star_count ){
						star_count = (pAnal1->backgrStat).back().nTnewCut;
					}
				}
			}
			printf("CLOUD_CHECK , stars count = %d\n",star_count);
			if( star_count <= gCCDParams.m_MinStarsToAccEvents ){
				printf("Number of stars on image = %d < %d , probably heavy clouds -> all events rejected !\n",star_count,gCCDParams.m_MinStarsToAccEvents);
				GetNewEvents().clear();
			}
		}


		// saving events to FITS :
		if(bAutoSaveEvents){
			SaveEventsToFITS();
		}


		if(bRet){	
			if(bAutoDump){
				MYTRACE2(gCCDTrace,"Interesting events found on frame# " << m_FrameCounter << ", dumping pipeline");				
				DumpPipeline();
			}else{
				MYTRACE2(gCCDTrace,"Interesting events found on frame# " << m_FrameCounter);							
			}
		}
		return bRet;
	}
	return FALSE;		
}

void  CCDPipeline::ReCalcMaxFrame()
{
	m_pPixelMaxFrame->SetData(0);
	for(int i=0;i<m_Count;i++){
		if(i!=head){
			UpdateMaximumsFrame(i);
		}
	}
}

LONG_T CCDPipeline::VerifyPreviousFrame( LONG_T FrametoVerify, LONG_T nextMatrixCount )
{
	LONG_T pos=-1;
	cCCD* pFrame = FindFrame( FrametoVerify, pos );
	LONG_T nVerified=0;
	if(pFrame){		
		LONG_T MatrixCount = m_pCCD[0].GetCount();
		CPixelAnalyseIn in;
		in.pPipeline = this;

		for(int i=0;i<MatrixCount;i++){
			in.ccd_index = i;
			in.pCamCfg = (CCcdCfg*)(GetCamCfgTab()[in.ccd_index]);
			in.pCCDInfo = &(GetCCDInfoTab()[in.ccd_index]);

			CCDMatrix* pMatrixToVerify = pFrame->GetMatrix(i);			
			GetNextMatrixPtrChronological( i, pos, in );
			// NextCount = GetNPreviousFrames( i, nextMatrixCount, NextMatrixPtr );
			Assert(in.PrevMatrixPtrCnt==(nextMatrixCount+1),"Number of next frames must be equal to : %d, but it is %d - error in code",(nextMatrixCount+1), in.PrevMatrixPtrCnt);

			// analysing frames after :	
			if(i < m_allFoundEvents.back().size())
				nVerified += m_pAnalyser->VerifyEventsOnImage( m_allFoundEvents.back()[i], in  );


			
		}
		if( !gCCDParams.m_bCCDDouble ){
   	   m_pAnalyser->VerifyTracks( *this, gCCDParams.m_ConfirmEventsOnNextNFrames );
   	}
	}
	return nVerified;
}

void CCDPipeline::SetNextFrameValues( LONG_T FrametoVerify, LONG_T nextMatrixCount )
{
	LONG_T pos=-1;
	cCCD* pFrame = FindFrame( FrametoVerify, pos );
	LONG_T nVerified=0;
	if(pFrame){		
		LONG_T MatrixCount = m_pCCD[0].GetCount();
		CPixelAnalyseIn in;
		in.pPipeline = this;

		for(int i=0;i<MatrixCount;i++){
			in.ccd_index = i;
			in.pCamCfg = (CCcdCfg*)(GetCamCfgTab()[in.ccd_index]);
			in.pCCDInfo = &(GetCCDInfoTab()[in.ccd_index]);

			CCDMatrix* pMatrixToVerify = pFrame->GetMatrix(i);			
			GetNextMatrixPtrChronological( i, pos, in );
			// NextCount = GetNPreviousFrames( i, nextMatrixCount, NextMatrixPtr );
			Assert(in.PrevMatrixPtrCnt==(nextMatrixCount+1),"Number of next frames must be equal to : %d, but it is %d - error in code",(nextMatrixCount+1), in.PrevMatrixPtrCnt);

			// analysing frames after :	
			if( m_allFoundEvents.size()>0 ){
				if(i < m_allFoundEvents.back().size()){
					m_pAnalyser->SetNextFrameValues( m_allFoundEvents.back()[i], in  );
				}
			}
		}
	}
}


BOOL_T CCDPipeline::VerifyOutFrame()
{
	if(GetCount()==m_nPipelineSize){
		// analysing outgoing frame if it was first flaged as intersting
		if(m_pCCD[head].GetInteresting()){
			clock_t t1 = clock();		
			ReCalcMaxFrame();
			long start = (head-1>=0 ? head-1 : m_nPipelineSize-1);			
			DumpPipeline( m_Count+1, start );
			clock_t t2 = clock();
	      mystring msg = get_clock_in_sec_string( t2-t1 );

	      MYTRACE1(gCCDTrace,"Verification of event frame took : " << msg );
			ClearState();
			return TRUE;
		}
	}
	return FALSE;
}

void CCDPipeline::DumpHomeopaticFrame(int i)
{
	if(gCCDParams.m_bDumpHomeoFrame){
		CFITSFile<BIG_ELEM_TYPE> out;
		mystring szName,szError;
		szName << gCCDParams.GetOutputDir() << "/homeo_" << m_FrameCounter << ".fit";
		if(!out.WriteToFITSFile( ((*m_pHomeopaticFrame)[i]).get_data_buffer(),
									   m_SizeX, m_SizeY, szError, szName.c_str() )){
			printf("could not write homeopatic frame to file %s\n",szName.c_str());
		}
	}
}

void CCDPipeline::ShiftFrame( CManyTab2D<BIG_ELEM_TYPE>* pFrame )
{
	if(pFrame){
		for(long i=0;i<m_nCCD;i++){
			_TRACE_PRINTF_6("shfting frame : %d\n",m_FrameCounter);
			((*pFrame)[i]).ShiftFrame( m_FrameCounter, gCCDParams.m_FrameDX, 
												gCCDParams.m_FrameDY, gCCDParams.m_bUseShiftTotal );
			//if(gCCDParams.m_bDumpHomeoFrame){
			//	DumpHomeopaticFrame(i);
		   //}
		}		
	}
}

void CCDPipeline::DumpFrame( CManyTab2D<BIG_ELEM_TYPE>* pFrame )
{
	for(int i=0;i<pFrame->GetCount();i++){
		CFITSFile<BIG_ELEM_TYPE> out;
		mystring szName,szError;
		szName << gCCDParams.GetOutputDir() << "/xxx_" << m_FrameCounter << ".fit";
		if(!out.WriteToFITSFile( ((*pFrame)[i]).get_data_buffer(),
										 m_SizeX, m_SizeY, szError, szName.c_str() )){
			printf("could not write frame to file %s\n",szName.c_str());
		}
	}
}

void CCDPipeline::ShiftHomeopaticFrame()
{
	if(m_pHomeopaticFrame){
		ShiftFrame( m_pHomeopaticFrame );
		if(gCCDParams.m_bDumpHomeoFrame){
			for(long i=0;i<m_nCCD;i++){
				DumpHomeopaticFrame(i);
		   }
		}
	}
}


void CCDPipeline::UpdateLaplaceFrame()
{
	long size = m_SizeX*m_SizeY;
	cCCD& newFrame = GetCurrent();
	LONG_T* neighb_list = gPipelineNeighbList;
	LONG_T ncnt;
	long ignore_edge = gCCDParams.GetEdgeSizeInLaplace();

	Assert( gCCDParams.GetAnalNeighbCount() < MAX_CLUSTER_SIZE,"To many neighbours required !");

	PROFILER_START
	LockWorkingFrames();
	for(long i=0;i<m_nCCD;i++){
		CCDMatrix& newMatrix = newFrame[i];
		Table2D<BIG_ELEM_TYPE>& laplaceMatrix = (*m_pLaplaceFrame)[i];
		ELEM_TYPE* pNewMatrix = newMatrix.get_data_buffer();
		ELEM_TYPE** pNewMatrixFast = newMatrix.get_data_buffer_fast();
		BIG_ELEM_TYPE* pHomeoMatrix = NULL;
		BIG_ELEM_TYPE** pHomeoMatrixFast = NULL;		

		// Assert(newMatrix.m_pFrameLaplace!=NULL,"Laplace frame of new frame must be calculated");
	   // BIG_ELEM_TYPE** p_laplace_of_new = (newMatrix.m_pFrameLaplace)->get_data_buffer_fast();
	
		BIG_ELEM_TYPE* pLaplaceMatrix = laplaceMatrix.get_data_buffer();		
		BIG_ELEM_TYPE** pLaplaceMatrixFast = laplaceMatrix.get_data_buffer_fast();
		long j=0;

		if(newMatrix.m_pFrameLaplace){
			UpdateHomeoFrame( laplaceMatrix, *newMatrix.m_pFrameLaplace );
		}else{
			newMatrix.Laplace( *m_pWrkMatrix );
			UpdateHomeoFrame( laplaceMatrix, *m_pWrkMatrix );
		}
	}

	PROFILER_END("Updating laplace frame took : ")
}

BOOL_T CCDPipeline::CalcAverageOfN( Table2D<BIG_ELEM_TYPE>& out, long prevN, long x0, long y0, long x1, long y1,
												BOOL_T bLaplace, CLongPoint* shapePoints, long pointCnt,
											   BOOL_T bReCalcLaplace/*=FALSE*/  )
{
	if(prevN>=m_Count){
		prevN = m_Count-1;
	}
	if(prevN<=0)
		return FALSE;

	CPixelAnalyseIn in;
   GetAllMatrixPtrsChronologicalInt( 0, in, TRUE );
	in.xSize = m_SizeX;
	in.ySize = m_SizeY;
	in.treshold_for_max = (long)(gCCDParams.m_nCalcMaxForAboveNSigma*CCDDataResults::GetSigmaBackground( gCCDParams.m_eLaplaceType ));


	if(x0<0)
		x0 = 0;
	if(y0<0)
		y0 = 0;
	if(x1<0 || x1>=m_SizeX)
		x1 = (m_SizeX-1);
	if(y1<0 || y1>=m_SizeY)
		y1 = (m_SizeY-1);

	BIG_ELEM_TYPE** p_av_data = out.get_data_buffer_fast();
	if(!bReCalcLaplace){
		for(register long y=y0;y<=y1;y++){
			in.y = y;
			for(register long x=x0;x<=x1;x++){
				in.x = x;
				p_av_data[y][x] = CCD_Analyser::CalcAverageOfPrevN( in, prevN, bLaplace, shapePoints, pointCnt );
			}
		}
	}else{
		for(register long y=y0;y<=y1;y++){
			in.y = y;
			for(register long x=x0;x<=x1;x++){
				in.x = x;
				p_av_data[y][x] = CCD_Analyser::CalcAverageOfPrevNRecalc( 
												in.x, in.y , in, prevN, 
												bLaplace, shapePoints, pointCnt );
			}
		}		
	}
	return TRUE;
}


BOOL_T CCDPipeline::CalcWeightedAverageOfN( CCDMatrix& out, long prevN, long x0, long y0, long x1, long y1,
												BOOL_T bLaplace, CLongPoint* shapePoints, long pointCnt )
{
	if(prevN>=m_Count){
		prevN = m_Count-1;
	}
	if(prevN<=0)
		return FALSE;

	
	InfoTable2D* pRotMap = GetRotationMap();
	Assert( pRotMap!=NULL,"Cannot calculate CalcWeightedAverageOfN");

/*	InfoTable2D* pRotMapOfPrev = new InfoTable2D[prevN+1]( pRotMap[0].m_SizeX,
																		  pRotMap[0].m_SizeY,
																		  pRotMap[0].m_dX,
																		  pRotMap[0].m_dY );*/
// version for gcc4.0 - NEW :
	InfoTable2D* pRotMapOfPrev = new InfoTable2D[prevN+1];
	int prevN_1 = prevN+1;
	for(int ii=0;ii<prevN_1;ii++){
		pRotMapOfPrev[ii].InfoTable2D_InitConstructor( pRotMap[0].m_SizeX,
                                                        pRotMap[0].m_SizeY,
                                                        pRotMap[0].m_dX,
                                                        pRotMap[0].m_dY );
	}

	/*ShiftInfo* pShiftInfo = new ShiftInfo[prevN];
	for(int i=0;i<prevN;i++){
		double dx = i*

		pShiftInfo[i]
	}*/

/*	CPixelAnalyseIn in;
   long allCount = GetAllMatrixPtrsChronologicalInt( 0, in, TRUE );
	in.xSize = m_SizeX;
	in.ySize = m_SizeY;

	for(register long f=0;f<=prevN;f++){
		for(long x=0;x<pRotMapOfPrev[f].m_X_count;x++){
			for(long y=0;y<pRotMapOfPrev[f].m_Y_count;y++){

				double dx_from_1_to_f = -(f-1)*(pRotMap[0].GetElem(x,y).m_LocalShiftInfo).m_LocalDX;
				double dy_from_1_to_f = -(f-1)*(pRotMap[0].GetElem(x,y).m_LocalShiftInfo).m_LocalDY;

				Table2D<ELEM_TYPE>::CalcOverlapParts2( dx_from_1_to_f, dy_from_1_to_f,
								pRotMapOfPrev[f].GetElem( x, y ).m_LocalShiftInfo );
				pRotMapOfPrev[f].GetElem( x, y ).m_LocalShiftInfo.m_LocalDX = dx_from_1_to_f;
				pRotMapOfPrev[f].GetElem( x, y ).m_LocalShiftInfo.m_LocalDY = dy_from_1_to_f;
			}
		}
	}

	if(x0<0)
		x0 = 0;
	if(y0<0)
		y0 = 0;
	if(x1<0 || x1>=m_SizeX)
		x1 = (m_SizeX-1);
	if(y1<0 || y1>=m_SizeY)
		y1 = (m_SizeY-1);
	
	ELEM_TYPE** p_av_data = out.get_data_buffer_fast();
	for(register long y=y0;y<=y1;y++){
		in.y = y;
		for(register long x=x0;x<=x1;x++){
			in.x = x;
			// Area2DInfo&  rotInfo = pRotMap[0].GetAreaDesc( x, y );
			
			long avValue = CCD_Analyser::CalcWeightedAverageOfPrevN( in, prevN, 
														pRotMapOfPrev, bLaplace );
			p_av_data[y][x] = avValue;
		}
	}*/
	delete [] pRotMapOfPrev;
	return TRUE;
}


BOOL_T CCDPipeline::CalcAverageOfN( CCDMatrix& out, long prevN, 
												long x0, long y0, long x1, long y1,
												BOOL_T bLaplace, CLongPoint* shapePoints, 
												long pointCnt, BOOL_T bReCalcLaplace /*=FALSE*/ )
{
	if(prevN>=m_Count){
		prevN = m_Count-1;
	}
	if(prevN<=0)
		return FALSE;

	CPixelAnalyseIn in;
   GetAllMatrixPtrsChronologicalInt( 0, in, TRUE );
	in.xSize = m_SizeX;
	in.ySize = m_SizeY;
	in.treshold_for_max = (long)(gCCDParams.m_nCalcMaxForAboveNSigma*CCDDataResults::GetSigmaBackground( gCCDParams.m_eLaplaceType ));

	if(x0<0)
		x0 = 0;
	if(y0<0)
		y0 = 0;
	if(x1<0 || x1>=m_SizeX)
		x1 = (m_SizeX-1);
	if(y1<0 || y1>=m_SizeY)
		y1 = (m_SizeY-1);

	ELEM_TYPE** p_av_data = out.get_data_buffer_fast();
	if(!bReCalcLaplace){	
		for(register long y=y0;y<=y1;y++){
			in.y = y;
			for(register long x=x0;x<=x1;x++){
				in.x = x;
				long avValue = CCD_Analyser::CalcAverageOfPrevN( in, prevN, bLaplace, 
																		shapePoints, pointCnt );
				p_av_data[y][x] = avValue;
			}
		}
	}else{
		for(register long y=y0;y<=y1;y++){
			in.y = y;
			for(register long x=x0;x<=x1;x++){
				in.x = x;
				long avValue = CCD_Analyser::CalcAverageOfPrevNRecalc( 
																		in.x, in.y , in, 
																		prevN, bLaplace, 
																		shapePoints, pointCnt );
				p_av_data[y][x] = avValue;
			}
		}
	}
	return TRUE;
}

long CCDPipeline::CalcAverageOfPrev( long prevN, long x, long y )
{
	if(prevN>=m_Count){
		prevN = m_Count-1;
	}
	if(prevN<=0)
		return FALSE;

	CPixelAnalyseIn in;
   GetAllMatrixPtrsChronologicalInt( 0, in, TRUE );
	in.xSize = m_SizeX;
	in.ySize = m_SizeY;
	in.x = x;
	in.y = y;
	in.treshold_for_max = (long)(gCCDParams.m_nCalcMaxForAboveNSigma*CCDDataResults::GetSigmaBackground( gCCDParams.m_eLaplaceType ));
	if(x<0 || y<0 || x>=in.xSize || y>=in.ySize){
		return 0;		
	}	
	long avValue = CCD_Analyser::CalcAverageOfPrevNRecalc( in.x, in.y , in,
                                                      prevN, TRUE , NULL, 0, TRUE );
	return avValue;
}

void CCDPipeline::CalcLaplaceOfNew()
{
	PROFILER_START
	cCCD& newFrame = GetCurrent();
	for(long i=0;i<m_nCCD;i++){
		CCDMatrix& newMatrix = newFrame[i];
		newMatrix.Laplace( TRUE );			
	}
	PROFILER_END("Calculating of laplace of new frame took : ");
}

void CCDPipeline::SubtractBackground()
{
	PROFILER_START
	cCCD& newFrame = GetCurrent();
	for(long i=0;i<m_nCCD;i++){
  	   CCDMatrix& newMatrix = newFrame[i];
		ELEM_TYPE** p_data = newMatrix.get_data_buffer_fast();
		Area2DInfo** pMap = m_pMatrixInfoMap[i].m_pTable2DMap;
		for(register int j=0;j<m_pMatrixInfoMap[i].m_Y_count;j++){
			for(register int k=0;k<m_pMatrixInfoMap[i].m_X_count;k++){
				LONG_T backgrValue = (LONG_T)(pMap[j][k].m_DataInfo[eRawS].m_Average);
				if( gCCDParams.m_eBackgrSubtrType == eMedianValue ){
					backgrValue = pMap[j][k].m_MedianValue;
				}
				// printf("subtracting baqckground = %d\n",backgrValue);
				for(register int y=(int)pMap[j][k].m_LowLeft.y;y<pMap[j][k].m_UpRight.y;y++){
					for(register int x=(int)pMap[j][k].m_LowLeft.x;x<pMap[j][k].m_UpRight.x;x++){
						p_data[y][x] -= backgrValue;
					}			
				}
			}
		}
	}
	PROFILER_END("PROFILER : Subtracting of background map took : ")
}

InfoTable2D* CCDPipeline::GetRotationMap()
{
	if(m_pMatrixInfoMap)
		return &(m_pMatrixInfoMap[0]);
	else
		return NULL;
}

void CCDPipeline::InitShiftsInfo( InfoTable2D& shiftInfoTab )
{
	Assert( shiftInfoTab.m_X_count==3 && shiftInfoTab.m_Y_count==2, "Currently only 3x2 map is handled for rotation");

	shiftInfoTab.GetElem( 0, 0 ).m_LocalShiftInfo.m_LocalDX = -1.64;
	shiftInfoTab.GetElem( 0, 0 ).m_LocalShiftInfo.m_LocalDY = 0.06;

	shiftInfoTab.GetElem( 1, 0 ).m_LocalShiftInfo.m_LocalDX = -1.71;
	shiftInfoTab.GetElem( 1, 0 ).m_LocalShiftInfo.m_LocalDY = 0.14;

	shiftInfoTab.GetElem( 2, 0 ).m_LocalShiftInfo.m_LocalDX = -1.65;
	shiftInfoTab.GetElem( 2, 0 ).m_LocalShiftInfo.m_LocalDY = 0.2;

	shiftInfoTab.GetElem( 0, 1 ).m_LocalShiftInfo.m_LocalDX = -1.74;
	shiftInfoTab.GetElem( 0, 1 ).m_LocalShiftInfo.m_LocalDY = 0.06;

	shiftInfoTab.GetElem( 1, 1 ).m_LocalShiftInfo.m_LocalDX = -1.77;
	shiftInfoTab.GetElem( 1, 1 ).m_LocalShiftInfo.m_LocalDY = 0.14;

	shiftInfoTab.GetElem( 2, 1 ).m_LocalShiftInfo.m_LocalDX = -1.77;
	shiftInfoTab.GetElem( 2, 1 ).m_LocalShiftInfo.m_LocalDY = 0.18;


	for(int x=0;x<shiftInfoTab.m_X_count;x++){
		for(int y=0;y<shiftInfoTab.m_Y_count;y++){
			Table2D<ELEM_TYPE>::CalcOverlapParts( shiftInfoTab.GetElem( x, y ).m_LocalShiftInfo.m_LocalDX,
															  shiftInfoTab.GetElem( x, y ).m_LocalShiftInfo.m_LocalDY,
															  shiftInfoTab.GetElem( x, y ).m_LocalShiftInfo.m_LocalS0,
															  shiftInfoTab.GetElem( x, y ).m_LocalShiftInfo.m_LocalS1,
															  shiftInfoTab.GetElem( x, y ).m_LocalShiftInfo.m_LocalS2,
															  shiftInfoTab.GetElem( x, y ).m_LocalShiftInfo.m_LocalS3,
															  shiftInfoTab.GetElem( x, y ).m_LocalShiftInfo.overlapledPoints );
		}
	}		
}

void CCDPipeline::GetBackgroundStat( CCDPipeline* pPipeline,
												 int CameraIdx, double& mS1, double& sS1, 
												 eLaplaceType_T laplaceType/*=eRawS*/ )
{
	if(pPipeline){
		pPipeline->GetBackground( CameraIdx, mS1, sS1, laplaceType );
	}else{
		mS1 = gAverageBacgroundG[ laplaceType ];
      sS1 = gSigmaBacgroundG[ laplaceType ];
	}
}

void CCDPipeline::GetBackgroundStat( CCDPipeline* pPipeline, int x, int y,
												 int CameraIdx, double& mS1, double& sS1, 
												 eLaplaceType_T laplaceType/*=eRawS*/ )
{
	if(pPipeline){
		pPipeline->GetBackground( x, y, CameraIdx, mS1, sS1, laplaceType );
	}else{
		mS1 = gAverageBacgroundG[ laplaceType ];
      sS1 = gSigmaBacgroundG[ laplaceType ];
	}
}

void CCDPipeline::GetBackground( int CameraIdx, double& mS1, double& sS1, eLaplaceType_T laplaceType/*=eRawS*/ )
{
	if( m_pMatrixInfoMap && gCCDParams.m_CalcBackgrTable[(int)laplaceType] ){
		mS1 = m_pMatrixInfoMap[CameraIdx].GetElem(0,0).m_DataInfo[(int)laplaceType].m_Average;
		sS1 = m_pMatrixInfoMap[CameraIdx].GetElem(0,0).m_DataInfo[(int)laplaceType].m_Sigma;
	}else{
		mS1 = gAverageBacgroundG[ laplaceType ];
		sS1 = gSigmaBacgroundG[ laplaceType ];
	}
}

void CCDPipeline::GetBackground( int x, int y, int CameraIdx, double& mS1, double& sS1, eLaplaceType_T laplaceType )
{
	Area2DInfo& desc = m_pMatrixInfoMap[CameraIdx].GetAreaDesc( x, y );
	mS1 = desc.m_DataInfo[(int)laplaceType].m_Average;
	sS1 = desc.m_DataInfo[(int)laplaceType].m_Sigma;
}

double CCDPipeline::GetTreshold( int x, int y, int CameraIdx, eLaplaceType_T laplaceType, 
										double sig_above, 
										BOOL_T bAver, int nAver )
{
	Area2DInfo& desc = m_pMatrixInfoMap[CameraIdx].GetAreaDesc( x, y );
	double mS1 = desc.m_DataInfo[(int)laplaceType].m_Average;
	double sS1 = desc.m_DataInfo[(int)laplaceType].m_Sigma;

	double treshold=0.0;
	if(!bAver){
		treshold = mS1 + sig_above*sS1;
	}else{
		treshold = mS1 + sig_above*(sS1/sqrt((double)nAver));
	}
	return treshold;
}

static CMyMutex oneAstrometryLock;
static BOOL_T DoAstrometry( CCDPipeline* pPipeline, BOOL_T bClearFlag=TRUE )
{
		pPipeline->m_bAstrometryRunning = TRUE;

		printf("Running DoAstrometry\n");fflush(0);		
		int level_sav=gPrintfLevel;
		gPrintfLevel=-10;
		mystring szLastFrame;
		pPipeline->GetLastFrame( szLastFrame );
		/*BOOL_T bSaveHere = FALSE;
		if( strlen(szLastFrame.c_str())==0 ){
			CCDMatrix& matrix = (pPipeline->GetCurrent())[0];
			szLastFrame << "for_astro_k" << pPipeline->GetPipelineIndex() << ".fit";
			matrix.WriteToFITSFile( szLastFrame.c_str(), FALSE, NULL, TRUE );
		}*/
		if( strlen(szLastFrame.c_str())==0 ){
			(pPipeline->GetPipelineCfg()).m_nSaveCurrentPicture=1;
		}
		while( strlen(szLastFrame.c_str())==0 ){
			printf_now("Waiting for frame to be saved ...\n");
			sleep(2);
			pPipeline->GetLastFrame( szLastFrame );
		}
		printf("Last frame = %s\n",szLastFrame.c_str());fflush(stdout);
		if( strlen(szLastFrame.c_str())>0 ){
			oneAstrometryLock.Lock();
			if(!pPipeline->ExecAstrometry( szLastFrame.c_str() )){
				mystring szErr;
				szErr << "Asas astrometry failed";
				gCCDErrLog.DumpToFile1( pPipeline, ASTRO_FAILED, szErr.c_str() );
			}
			oneAstrometryLock.UnLock();
		}
		gPrintfLevel = level_sav;
		
		if( bClearFlag ){			
			pPipeline->m_bAstrometryRunning = FALSE;
			printf("Cleared flag m_bAstrometryRunning\n");fflush(0);
		}
		
		return pPipeline->IsAstrometryOK();
}

static BOOL_T DoAstrometryAsynchro( CCDPipeline* pPipeline, BOOL_T bClearFlag=TRUE )
{
		BOOL_T bRet=FALSE;		
		pPipeline->m_bAstrometryRunning = TRUE;

		printf("Running DoAstrometryAsynchro\n");fflush(0);
		int level_sav=gPrintfLevel;
		gPrintfLevel=-10;
		mystring szLastFrame;
		pPipeline->GetLastFrame( szLastFrame );
		if( strlen(szLastFrame.c_str())==0 ){
			(pPipeline->GetPipelineCfg()).m_nSaveCurrentPicture=1;
		}
		while( strlen(szLastFrame.c_str())==0 ){
			printf_now("Waiting for frame to be saved ...\n");
			sleep(2);
			pPipeline->GetLastFrame( szLastFrame );
		}	
		if( strlen(szLastFrame.c_str())>0 ){
			if( strcmp( pPipeline->m_szLastAsynchroAstroFrame.c_str(), szLastFrame.c_str()) ){
				printf("NEW frame %s != %s, running astro/photo in asynchro mode ...\n",
							szLastFrame.c_str(),pPipeline->m_szLastAsynchroAstroFrame.c_str());
			}else{
				printf("NEW frame %s = %s, astrometry skiped\n",
							szLastFrame.c_str(),pPipeline->m_szLastAsynchroAstroFrame.c_str());
				pPipeline->m_bAstrometryRunning = FALSE;
				return FALSE;
			}

			oneAstrometryLock.Lock();
			bRet = pPipeline->ExecAstrometryAsynchro( szLastFrame.c_str() );
			if(!bRet){
				mystring szErr;
				szErr << "Asas astrometry failed";
				gCCDErrLog.DumpToFile1( pPipeline, ASTRO_FAILED, szErr.c_str() );
			}
			pPipeline->m_szLastAsynchroAstroFrame = szLastFrame.c_str();
			oneAstrometryLock.UnLock();
		}
		gPrintfLevel = level_sav;
		
		if( bClearFlag ){
			pPipeline->m_bAstrometryRunning = FALSE;
			printf("Cleared flag m_bAstrometryRunning\n");fflush(0);
		}
		
		return bRet;
}



BOOL_T CCDPipeline::IsAstrometryOK()
{ 
	return m_pAsasTransform->m_bTransformOK;
}


BOOL_T CCDPipeline::ExecAstrometry( const char* szName  )
{
	CCDMatrix image1(0,0);
	if( !image1.ReadFITSFile( szName )){
		printf("could not read file : %s\n",szName );
		return FALSE;
	}
	if( m_PipelineCfg.m_bAsasSubtrDark ){
		cCCD* pDark1 = GetDark();
		if( pDark1 ){
        	image1.Subtract( (*pDark1)[0], image1, TRUE );
			printf("Dark image subtracted\n");fflush(stdout);
		}
	}
	if( m_PipelineCfg.m_bAsasDivideByFlat ){
		Table2D<float>* pFlat1 = GetFlat();
		if( pFlat1 ){
			CCDUtil::Divide( image1, (*pFlat1), 0 );					
		}
	}
	mystring szRedName1;
	szRedName1 << getfname( szName ) << "_red.fit";
	image1.WriteToFITSFile( szRedName1.c_str(), FALSE, NULL, TRUE );

	printf("DEBUG : before RunAsasAstrometry\n");fflush(stdout);
	BOOL_T bRet = RunAsasAstrometry( szRedName1.c_str(), &image1 );

// ms 20070507 - updating also failed astrometry
//	if( bRet ){
		if( !gCCDParams.GetMC() ){
			printf("Updating FITS header in current file : %s ...\n",szName );
			UpdateCurrentFITSFile( szName );
		}	
//	}

	mystring szCmd;
	szCmd << "rm -f " << szRedName1.c_str();
	system( szCmd.c_str() );

	return bRet;
}

BOOL_T CCDPipeline::ExecAstrometryAsynchro( const char* szName  )
{
	CCDMatrix image1(0,0);
	if( !image1.ReadFITSFile( szName )){
		printf("could not read file : %s\n",szName );
		return FALSE;
	}
	if( m_PipelineCfg.m_bAsasSubtrDark ){
		cCCD* pDark1 = GetDark();
		if( pDark1 ){
        	image1.Subtract( (*pDark1)[0], image1, TRUE );
			printf("Dark image subtracted\n");fflush(stdout);
		}
	}
	if( m_PipelineCfg.m_bAsasDivideByFlat ){
		Table2D<float>* pFlat1 = GetFlat();
		if( pFlat1 ){
			CCDUtil::Divide( image1, (*pFlat1), 0 );					
		}
	}
	mystring szRedName1;
	szRedName1 << getfname( szName ) << "_red.fit";
	image1.WriteToFITSFile( szRedName1.c_str(), FALSE, NULL, TRUE );


	// copying current VALID astrometry structure :
	memcpy( m_pAsasTransformAsynchro, m_pAsasTransform, sizeof(CCDAsasTransform));

	// running astrometry using temporary structure , special for asynchro
	// astrometry VALID one remains unchanged :	
	BOOL_T bRet = RunAsasAstrometry( szRedName1.c_str(), &image1, m_pAsasTransformAsynchro );

//ms 20070507 - updating also failed astrometry !
//	if( bRet ){
		if( !gCCDParams.GetMC() ){
			printf("Updating FITS header in current file : %s ...\n",szName );
			UpdateCurrentFITSFile( szName, m_pAsasTransformAsynchro );
		}	
//	}

	mystring szCmd;
	szCmd << "rm -f " << szRedName1.c_str();
	system( szCmd.c_str() );

	return bRet;
}

BOOL_T CCDPipeline::RemoveTemporaryFiles( const char* szMagFile,
														const char* szFlipedName,
														const char* szAstFile )
{
	printf("cleaning FITS/MAG/AST used for astrometry ...\n");
   mystring szCmd;

	if( szFlipedName && szFlipedName[0] ){		
	   szCmd << "rm -f " << szFlipedName;
		printf("Cleaning fliped files : %s\n",szCmd.c_str());
   	system( szCmd.c_str() );
	}

	if( !gCCDParams.m_bKeepMagAndAstFromAstro ){
		if( szMagFile && szMagFile[0] ){	
			szCmd = "";
			szCmd << "rm -f " << szMagFile;
			system( szCmd.c_str() );
		}

		if( szAstFile && szAstFile[0] ){
			szCmd = "";
   		szCmd << "rm -f " << szAstFile;
		   system( szCmd.c_str() );
		}

		const char* ptr_mag = "";
		const char* ptr_ast = "";
		if( szMagFile ) ptr_mag = szMagFile;
		if( szAstFile ) ptr_ast = szAstFile;	
		printf("Cleaning mag/ast files : %s / %s\n",ptr_mag,ptr_ast);
	}

	return TRUE;
}

BOOL_T CCDPipeline::RunPhotometry( const char* szFile, CCDMatrix& matrix, BOOL_T bSaveMag )
{
	BOOL_T bRet=FALSE;
	mystring szMagFile;
   szMagFile << getfname( szFile ) << ".mag";
	
	if( gCCDParams.m_bUseFastPhotoInAstro ){
		PROFILER_START

		printf_now2("CCDPipeline::RunPhotometry : running FAST photometry on file : %s...\n",szFile);fflush(stdout);
		m_nPhotometryStarsCount = RunFastPhotometry( matrix, m_PipelineCfg.m_fASASPhotoThres, szMagFile, bSaveMag );
		if( m_nPhotometryStarsCount>0 )
			bRet = TRUE;
		else
			bRet = FALSE;
		printf("FAST_PHOTOMETRY stars# = %d\n",m_nPhotometryStarsCount);
		if(bRet){
			printf("OK : %d stars detected on frame\n",m_nPhotometryStarsCount);
		}else{
			printf("FAILED\n");
		}
		fflush(stdout);
      PROFILER_END("FAST-photometry took :");
	}else{
		printf("This option is not supported, only woring with fast photometry enabled\n");
	}

	return bRet;	
}

BOOL_T CCDPipeline::RunAsasAstrometry( const char* szReductFile, 
													CCDMatrix* pMatrix,
													CCDAsasTransform* pAsasTransform )
{	
	if( !pAsasTransform ){
		pAsasTransform = m_pAsasTransform;
	}
	printf("Previous astrometry parameters (ra,dec,fi,pixscale)=(%.2f,%.2f,%.2f,%.2f)\n",
				pAsasTransform->ra,pAsasTransform->dec,pAsasTransform->fi,
				pAsasTransform->pixscale);fflush(stdout);

	// m_nPhotometryStarsCount=0;
	BOOL_T bRet=FALSE;
	const char* szFileForAstro = szReductFile;

	time_t t1 = get_dttm();

	if( strlen( szReductFile )){
		// [NEW] 20041019 - in case no flip needed init file names :
		mystring szFlipedName=szReductFile;
		mystring szMagFile;
		BOOL_T bFliped=FALSE;
		szMagFile << getfname( szReductFile ) << ".mag";

		PROFILER_START
		// bRet=TRUE;szMagFile << szReductFile << ".mag";
		CCDMatrix matrix(0,0);
		if( m_PipelineCfg.m_eReverseForTransform!=eReverseImageNone ){
			if(matrix.ReadFITSFile( szReductFile )){
				BOOL_T bFlip=FALSE;
				if( m_PipelineCfg.m_eReverseForTransform==eReverseImageHor || m_PipelineCfg.m_eReverseForTransform==eReverseImageHorVert ){
					matrix.FlipImage( eReverseImageHor );
					bFlip=TRUE;
				}
				if( m_PipelineCfg.m_eReverseForTransform==eReverseImageVert || m_PipelineCfg.m_eReverseForTransform==eReverseImageHorVert ){
					matrix.FlipImage( eReverseImageVert );
					bFlip=TRUE;
				}
								
				szFlipedName <<= getfname(szReductFile);
				if( bFlip ){
					mystring szFlipDesc = GetFlipDesc( m_PipelineCfg.m_eReverseForTransform );
					szFlipedName << "_" << szFlipDesc.getlower();
				}
				szFlipedName << ".fit";
				CKeyTab keys;
				keys.Add( FLIP, GetFlipDesc( m_PipelineCfg.m_eReverseForTransform ) );
				keys.Set( "FLAT_DIV", "FLAT" );
				keys.Set( "DARK_SUB", "DARK" );
				keys.Set( "PIXSCALE", (float)m_PipelineCfg.m_fPixScale );
				if( m_PipelineCfg.m_nAsasBorderSize>0 && !gCCDParams.m_bUseFastPhotoInAstro ){
					// matrix.ClearBorder( m_PipelineCfg.m_nAsasBorderSize );
					printf("cuting off image border ...");fflush(0);

					// old version CatBorder :
					// matrix.CatBorder( m_PipelineCfg.m_nAsasBorderSize , TRUE );
					matrix.CatBorderFast( m_PipelineCfg.m_nAsasBorderSize );

					printf("OK\n");fflush(0);
				}
				matrix.WriteToFITSFile( szFlipedName.c_str(), FALSE, &keys, TRUE );		
				szReductFile = szFlipedName.c_str();
				bFliped = TRUE;
			}else{
				printf("Error reading reduced file : %s\n",szReductFile);
				mystring szErr;
				szErr << "Error reading reduced file : " << szReductFile;
				gCCDErrLog.DumpToFile1( this, ASTRO_READ_REDUCED_FAILED, szErr.c_str() );
				m_bPrevAstrometryOK = FALSE;
				return FALSE;
			}
		}else{
			szFlipedName <<= getfname(szReductFile);
			szFlipedName << "_none.fit";
         if(matrix.ReadFITSFile( szReductFile )){
				CKeyTab keys;
				keys.Add( FLIP, "NONE" );
            keys.Set( "FLAT_DIV", "FLAT" );
            keys.Set( "DARK_SUB", "DARK" );
            keys.Set( "PIXSCALE", (float)m_PipelineCfg.m_fPixScale );
				matrix.WriteToFITSFile( szFlipedName.c_str(), FALSE, &keys, TRUE );
			}
		}

		m_FwhmAver = 0;
		if( gCCDParams.m_bUseFastPhotoInAstro ){
			printf_now2("CCDPipeline::RunAsasAstrometry : running FAST photometry on file : %s...\n",szFlipedName.c_str());fflush(0);
			m_nPhotometryStarsCount = RunFastPhotometry( matrix, m_PipelineCfg.m_fASASPhotoThres, szMagFile );
			if( m_nPhotometryStarsCount>0 )
				bRet = TRUE;
			else
				bRet = FALSE;
			printf("FAST_PHOTOMETRY stars# = %d\n",m_nPhotometryStarsCount);
			if(bRet){
				printf("OK : %d stars detected on frame\n",m_nPhotometryStarsCount);
			}else{
				printf("FAILED\n");
			}
			fflush(0);
         PROFILER_END("FAST-photometry took :");
		}else{
			// cannot check number of stars when using ASAS photometry :
			m_nPhotometryStarsCount = 0;

			printf_now2("running ASAS photometry on file : %s...\n",szFlipedName.c_str());fflush(0);
			bRet = pAsasTransform->asas_photometry( szFlipedName.c_str(),
																m_PipelineCfg.m_fASASPhotoThres,
																szMagFile );			
			if( MyFile::DoesFileExist( szMagFile.c_str() ) ){
				CSafeKeyTab keyTab;
			   CFITSFile<ELEM_TYPE> in;
			   if(!in.ReadFITSHeader( keyTab, szMagFile.c_str() )){
			      printf("could not read FITS file : %s\n",szMagFile.c_str() );
			      exit(-1);
			   }
				const char* szNY = keyTab.getKeyVal( FH_NAXIS2 );
				if( szNY && szNY[0] ){
					m_nPhotometryStarsCount = atol( szNY );
				}
				const char* szFWHM = keyTab.getKeyVal( "FWHM" );
				if( szFWHM && szFWHM[0] ){
					m_FwhmAver = atof( szFWHM );
				}
			}

			if(bRet){
				printf("OK : %d stars detected on frame\n",m_nPhotometryStarsCount);
			}else{
				printf("FAILED\n");
			}
			fflush(0);
			PROFILER_END("ASAS-photometry took :");
		}


		if( m_PipelineCfg.m_nMinStarCountToRunAstro>0 ){
			if( m_nPhotometryStarsCount < m_PipelineCfg.m_nMinStarCountToRunAstro ){
				printf("ASTROMETRY : to few stars %d ( limit is %d ) identified on frame, astrometry skiped\n",m_nPhotometryStarsCount,m_PipelineCfg.m_nMinStarCountToRunAstro);

				// changed 20060801 :
				pAsasTransform->ResetAstroOK();
				// pAsasTransform->m_bTransformOK = FALSE;

				printf("FAILED - on star number check\n");
				RemoveTemporaryFiles( szMagFile.c_str(), szFlipedName.c_str(), NULL );
				m_bPrevAstrometryOK = FALSE;
				return FALSE;
			}
		}

		time_t ut_time=0;
		if(bRet){
			mystring szAstFile;
			szAstFile << szReductFile << ".ast";
	
			double ra0,dec0,alt,azim,ha;
			eObservationMode_T obsMode;

			if( pMatrix ){
				GetCurrentFrameCoo( *pMatrix, ut_time, ra0, dec0, alt, azim, ha, obsMode );				
			}else{
				GetCurrentFrameCoo( GetCurrent()[0], ut_time, ra0, dec0, alt, azim, ha, obsMode );
			}
			ra0 = AstroAngle::rad2hours( ra0 );
			dec0 = AstroAngle::rad2deg( dec0 );
			pAsasTransform->ra = ra0;
			pAsasTransform->dec = dec0;
			pAsasTransform->m_TransformUtTime = ut_time;
			// ra0 = AstroAngle::rad2deg( m_PipelineCfg.m_RAObs + m_PipelineCfg.m_RACorr );
			// dec0 = AstroAngle::rad2deg( m_PipelineCfg.m_DecObs + m_PipelineCfg.m_DecCorr );
			double pixscale = m_PipelineCfg.GetPixScale();
			double fi = m_PipelineCfg.m_fASASAstrometryFi;
			if( pAsasTransform->m_bTransformOK ){
				pixscale = pAsasTransform->pixscale;
// problems with astrometry in LCO - what for was this line ????
// I comment it now - I hope forever !!!
//				fi = (180.00 - pAsasTransform->fi);
			}
			
			printf("running asas astrometry on file (ra,dec)=(%.2f,%.2f) , (ra0,dec0)=(%.2f,%.2f) : %s...\n",
						pAsasTransform->ra,pAsasTransform->dec,
						ra0, dec0,
						szMagFile.c_str());fflush(0);
		
			CCDAsasTransform sav( NULL, NULL );
			memcpy( &sav, pAsasTransform, sizeof(CCDAsasTransform));
			PROFILER_RESTART					

			BOOL_T bAstroVerbose=gCCDParams.m_bASASAstrometryVerb;
			if( gCCDParams.m_nChangeToSilentAfterNGood>0 ){
				if( m_nGoodAstrometryCount >= gCCDParams.m_nChangeToSilentAfterNGood ){
					printf("Subsequent %d astrometries OK, disabling verbose mode now\n",m_nGoodAstrometryCount);fflush(0);
					bAstroVerbose = FALSE;	
				}
			}

			int retry_initial=2;
			int retry=retry_initial;

			BOOL_T bAstroBreak=FALSE;
			if( m_WorkingMode.m_LastCommand == eDAQReq_StopAnalysis ){
				printf("Last Command = STOP ANALYSIS - skipping current astrometry !\n");fflush(stdout);
				retry=0;					
				bAstroBreak = TRUE;
			}

			bRet = FALSE;	// to force at least 1 astrometry 
			
			double ast_err_req = m_PipelineCfg.m_fAsasError;
			double ast_err_fatal_req = m_PipelineCfg.m_fAsasFatalError;

			while( retry>0 && !bRet ){
				if( retry < 2 ){
					printf("RETRYING ASTROMETRY !!! with fi=%.2f, retry=%d\n",fi,retry);
				}
				printf("Starting astrometry on file : %s, start params (ra,dec,fi,pixscale,ord,try,ast_err,ast_err_fatal)=(%.4f,%.4f,%.2f,%.2f,%d,%d,%.2f,%.2f)\n",
						szMagFile.c_str(), ra0,dec0,fi,pixscale,
						gCCDParams.m_nASASAstrometryOrd,m_PipelineCfg.m_nASASAstrometryTry,
						ast_err_req, ast_err_fatal_req);fflush(0);
				bRet = pAsasTransform->asas_astrometry( szMagFile.c_str(), szAstFile.c_str(), 
															   ra0, dec0, 
																fi,
																pixscale,
																gCCDParams.m_nASASAstrometryOrd, 
																gCCDParams.m_szASASStarCatalog.c_str(),
																bAstroVerbose,	
																m_PipelineCfg.m_nASASAstrometryTry,
																ast_err_req, ast_err_fatal_req );
				printf("asas_astrometry bRet = %d , ast_err = %.2f, last_err_code = %d\n",bRet,pAsasTransform->ast_err,pAsasTransform->m_LastAstroRetCode);fflush(stdout);
	
				if( !bRet && retry==2 ){
					if( fabs(fi)<5.00 ){
						fi = 180.00 - fi;
					}else{
						fi = 0.00;
					}
				}
//				if( !bRet && pAsasTransform->AstroBreakForced()>0 ){
				if( !bRet && pAsasTransform->WasAstroCanceled() ){
					printf("INFO : astrometry break was forced, no retry performed\n");
					bAstroBreak = TRUE;
					break;
				}

				// checking astrometry error in case to large and still one retry 
				// is possible force repetition of astrometry with stricter limit 
				// for astrometry error :
//				if( bRet && pAsasTransform->ast_err > 2.00 && retry>1 ){
//					printf("WARNING : astrometry error is large, astrometry will be repeated !\n");fflush(stdout);
//					bAstroVerbose = TRUE;
//					ast_err_req = ast_err_req / 2.00;
//					ast_err_fatal_req = ast_err_fatal_req / 2.00;
//					bRet = FALSE;
//				}
				if( !bRet && pAsasTransform->m_LastAstroRetCode == E_astrometry_timeout && retry>1 ){
					printf("WARNING : astrometry failed, will repeat with larger error tolerance\n");
					ast_err_req = ast_err_req * 2.00;
					ast_err_fatal_req = ast_err_fatal_req * 2.00;
					bAstroVerbose = TRUE;
				}

				retry--;
			}

			if( bAstroBreak ){
				// astrometry was canceled due to StopAnalysis command :
				printf("INFO : astrometry was externally cancelled\n");fflush(stdout);
				return FALSE;
			}

			m_bPrevAstrometryOK = bRet;
			if(bRet){
				printf("ASTRO-OK ( file %s )\n",szMagFile.c_str());
				m_AstFilesList.Add( szAstFile.c_str() );
				m_nGoodAstrometryCount++;
				pAsasTransform->m_nFailedAstrometryCount=0; // cleaning failed counter
			}else{
				printf("ASTRO-FAILED ( file %s ) , cleaning good astro counter\n",szMagFile.c_str());
				m_nGoodAstrometryCount=0;
				pAsasTransform->m_nFailedAstrometryCount++; // increase failed counter
													 // when skiped due to to low # of stars not increased
													 // this is indication of bad coord in DAQ
													 // astrometry is performed but not succeeded
			}
			fflush(0);			
			PROFILER_END("Asas astrometry took :");			
						
			if( bRet ){
				if( CheckIfChangeCoordBig( sav.ra, sav.dec, 
													AstroAngle::rad2deg( sav.m_AzimInRad ),  
													AstroAngle::rad2deg( sav.m_AltInRad ),
													pAsasTransform ) ){
					printf("\n\nINFO : BIG change of coordinates detected after astrometry, will restart daq pipeline NOW\n\n");
				}
			}	

			if( bFliped ){
				printf("cleaning FITS used for astrometry ...\n");
				mystring szCmd;
				szCmd << "rm -f " << szFlipedName.c_str();
				system( szCmd.c_str() );
				printf("OK\n");
			}

			if(!bRet){
				printf("\n\nASAS astrometry FAILED !!!\n");
				// NEW - 20041027 , now when astrometry is started must be done
				// cannot move back to old one :
				//if( !m_pAsasTransform->m_bTransformOK ){
				//	   printf("restoring old one ...\n");
				// 	memcpy( m_pAsasTransform, &sav, sizeof(CCDAsasTransform));
				//}
			}else{
				char szFile[1024];
				sprintf(szFile, CCD_ASTROMETRY_FILE, m_PipelineIndex );
				pAsasTransform->SaveToFile( szFile, GetCamObsMode() );

					
				mystring szSave,szLogFile;
				szLogFile << "frame_" << GetDayFrameCounter() << "_" << szFile;
				CCDLog::GetOutPutFileName( szSave, this, "Astrometry", szLogFile.c_str() );
				pAsasTransform->SaveToFile( szSave.c_str(), GetCamObsMode() );

				// set new coordinates :
				UpdateObsCoordinates( AstroAngle::hours2rad(pAsasTransform->ra), 
											 AstroAngle::deg2rad( pAsasTransform->dec ),
											 pAsasTransform->m_AzimInRad,
											 pAsasTransform->m_AltInRad ); 
			}
//			if( !gCCDParams.m_bKeepMagAndAstFromAstro && bRet ){
			if( !gCCDParams.m_bKeepMagAndAstFromAstro ){ // only in this version transienlib
				// cleaning mag and ast file :
				mystring szCmd;
				szCmd << "rm -f " << szMagFile.c_str();
				system( szCmd.c_str() );
				szCmd = "";
				szCmd << "rm -f " << szAstFile.c_str();
				system( szCmd.c_str() );
			}else{
				printf("Astrometry failed keeping files : %s and %s\n",szMagFile.c_str(),szAstFile.c_str());
			}

			CCDLog astrometryLog( "%d\t%.8f\t%.8f\t%.8f\t%.8f %d %s\n","Result RA(h) DEC(deg) AZIM(deg) ALT(deg) MODE MAG_FILE","Astrometry","astrometry");
			astrometryLog.DumpToFile1( this, (int)bRet, pAsasTransform->ra, pAsasTransform->dec,
												AstroAngle::rad2deg( pAsasTransform->m_AzimInRad ),
												AstroAngle::rad2deg( pAsasTransform->m_AltInRad ),
												obsMode,
												szMagFile.c_str() );
		}else{
			printf("\n\nASAS astrometry FAILED !!!\n");			
		}
		if( bRet ){
			// setting only if good - if not good remain as it was 
			// - if ok , use old one , if not - 
			// then nothing changes still BAD
			printf("Setting AstroOK to %d, last_good_time=%d (astro_time=%d)\n",bRet,ut_time,pAsasTransform->m_TransformUtTime);fflush(0);
			pAsasTransform->m_bTransformOK = bRet;
			pAsasTransform->m_LastGoodAstroTime = ut_time;
		}
	}

	time_t t2 = get_dttm();
	printf("ASTROMETRY_TIME (%d) = %d sec results=%d, m_bTransformOK=%d\n",m_PipelineIndex,(t2-t1),bRet,pAsasTransform->m_bTransformOK);fflush(stdout);

	return bRet;
}


BOOL_T CCDPipeline::IsAsasTransformOK()
{
	if( m_pAsasTransform ){
		return m_pAsasTransform->m_bTransformOK;
	}
	return FALSE;
}

BOOL_T CCDPipeline::SaveReducedFrame( mystring& szReductFile )
{
	if(m_PipelineCfg.m_nSaveReductFramesFreq>0){
		if(((m_FrameCounter-m_StartFrameIndex) % m_PipelineCfg.m_nSaveReductFramesFreq==0) || 
			(m_FrameCounter-m_StartFrameIndex)==1){
			CKeyTab redKeys;
			redKeys.Add( "DARK_SUB", "DARK" );
			redKeys.Add( "FLAT_DIV", "FLAT" );
			

			cCCD& cNew = GetCurrent();
			char szName[1024];
			sprintf( szName, "red_%.5d.fits",m_FrameCounter);		
			cNew[0].WriteToFITSFile( szName, FALSE, &redKeys );
			m_PipelineCfg.m_szLastReducedFrameName = szName;
			szReductFile = m_PipelineCfg.m_szLastReducedFrameName;			
			return TRUE;
		}		
	}
	return FALSE;
}

void CCDPipeline::UpdateBackgroundMap()
{
	if(gCCDParams.GetCalcBackgrFlag() && m_pMatrixInfoMap){
		int rest = ((m_FrameCounter-m_StartFrameIndex) % gCCDParams.m_BackgUpdateFreq);
		_TRACE_PRINTF_3("m_FrameCounter=%d, m_StartFrameIndex=%d\n",m_FrameCounter,m_StartFrameIndex);
		if( ( rest==0 ) || (m_FrameCounter-m_StartFrameIndex)==1){
			PROFILER_START
			cCCD& newFrame = GetCurrent();
			long i;
			for(i=0;i<m_nCCD;i++){
	   	   CCDMatrix& newMatrix = newFrame[i];
				newMatrix.CalcTableMapInfo( m_pMatrixInfoMap[i], 
													 gCCDParams.m_MaxAllowedVal-1,
													 gCCDParams.m_CalcBackgrTable,
													 gCCDParams.m_eLaplaceType,
													 newMatrix.get_frame_laplace_fast(),
													 gCCDParams.m_eBackgrSubtrType==eMostPopularValue,
													 gCCDParams.m_eBackgrSubtrType==eMedianValue,
													 m_FrameCounter );
				// log baground map to log :
				LogBackgroundMap();

				if(gPrintfLevel>=1){
					mystring szBackgr;
					newMatrix.GetBackgrDescStr( m_pMatrixInfoMap[i], szBackgr );
					_TRACE_PRINTF_1("Background map calculated : %s\n",szBackgr.c_str());
				}
			}
			PROFILER_END("Update of background map took : ")
		}
	}
}

void CCDPipeline::LogBackgroundMap()
{
	if(gCCDParams.m_bBackgrDump){
		for(register int lap=0;lap<MAX_LAPLACE_DEFINED;lap++){
			eLaplaceType_T lap_type = (eLaplaceType_T)lap;
			if( gCCDParams.m_CalcBackgrTable[lap] ){
				mystring szMapSigma,szMapAverage;
				char szTmp[128];
				szMapSigma << "Frame# " << (GetDayFrameCounter()+1) << "\n";
				szMapAverage << "Frame# " << (GetDayFrameCounter()+1) << "\n";
				for(int j=(m_pMatrixInfoMap[0].m_Y_count-1);j>=0;j--){
					for(int k=0;k<m_pMatrixInfoMap[0].m_X_count;k++){
						sprintf(szTmp,"%.4d",(int)m_pMatrixInfoMap[0].m_pTable2DMap[j][k].m_DataInfo[lap].m_Sigma);
						szMapSigma << szTmp << " ";

 						sprintf(szTmp,"%.4d",(int)m_pMatrixInfoMap[0].m_pTable2DMap[j][k].m_DataInfo[lap].m_Average);
						szMapAverage << szTmp << " ";							
					}
					szMapSigma   << "\n";
					szMapAverage << "\n";
				}
				mystring szBase,szLogMean,szLogSigma;
				szBase << CCDMatrix::GetLaplaceName( lap_type ) << "_mean_map.log" ;
				GetOutFileName( szLogMean, "FrameAnalyseInfoDetails", szBase.c_str(),
									 this, -1, FALSE );
				szBase = "";
				szBase << CCDMatrix::GetLaplaceName( lap_type ) << "_sigma_map.log" ;
				GetOutFileName( szLogSigma, "FrameAnalyseInfoDetails", szBase.c_str(),
									 this, -1, FALSE );

				MyOFile outM( szLogMean.c_str(), "a"  );
				outM.Printf("%s\n",szMapAverage.c_str());

				MyOFile outS( szLogSigma.c_str(), "a" );
				outS.Printf("%s\n",szMapSigma.c_str());
			}
		}
	}	
}

void CCDPipeline::UpdateMaximumsFrameOpt(long iFrame)
{
	long size = m_SizeX*m_SizeY;
	cCCD& newFrame = GetCurrent();
	if(iFrame>=0)
		newFrame = m_pCCD[iFrame];

	clock_t t1 = clock();
	LONG_T neighb_list[MAX_CLUSTER_SIZE],out_list[MAX_CLUSTER_SIZE];
	LONG_T ncnt,ocnt;
	
	Assert( gCCDParams.GetAnalNeighbCount() < MAX_CLUSTER_SIZE,"To many neighbours required !");

	for(long i=0;i<m_nCCD;i++){
		CCDMatrix& newMatrix = newFrame[i];
		CCDMatrix& maxMatrix = (*m_pPixelMaxFrame)[i];
		ELEM_TYPE* pNewMatrix = newMatrix.get_data_buffer();
		ELEM_TYPE* pMaxMatrix = maxMatrix.get_data_buffer();
		for(long j=0;j<size;j++){
			long x = ( j % m_SizeX );
		   long y =  ( j / m_SizeX );
			long ncnt = CCD_Analyser::GetAnalNeighbOpt( x , y , m_SizeX, m_SizeY,
                                                     neighb_list, ncnt,
                                                     out_list, ocnt );
			LONGLONG_T newSum = CCD_Analyser::CalcSumOpt( neighb_list, ncnt, pNewMatrix );			
			LONGLONG_T maxSum = pMaxMatrix[j]; // already calculated sum of neighbours
			if(newSum>maxSum){
				pMaxMatrix[j] = newSum;
			}
		}
	}
	clock_t t2 = clock();
   mystring msg = get_clock_in_sec_string( t2-t1 );

   MYTRACE2(gCCDTrace,"Updating maximum frame took : " << msg );
}


void CCDPipeline::UpdateMaximumsFrame(long iFrame)
{
	long size = m_SizeX*m_SizeY;
	cCCD& newFrame = GetCurrent();
	if(iFrame>=0)
		newFrame = m_pCCD[iFrame];

	clock_t t1 = clock();
	for(long i=0;i<m_nCCD;i++){
		CCDMatrix& newMatrix = newFrame[i];
		CCDMatrix& maxMatrix = (*m_pPixelMaxFrame)[i];
		ELEM_TYPE* pNewMatrix = newMatrix.get_data_buffer();
		ELEM_TYPE* pMaxMatrix = maxMatrix.get_data_buffer();
		CLongList neighb_list,outer;
		for(long j=0;j<size;j++){
			long x = ( j % m_SizeX );
		   long y =  ( j / m_SizeX );
			long neighb_count = CCD_Analyser::GetAnalNeighb( x , y , m_SizeY, m_SizeY, neighb_list, outer);
			LONGLONG_T newSum = CCD_Analyser::CalcSum( neighb_list, pNewMatrix );			
			LONGLONG_T maxSum = pMaxMatrix[j]; // already calculated sum of neighbours
			if(newSum>maxSum){
				pMaxMatrix[j] = newSum;
			}
		}
	}
	clock_t t2 = clock();
   mystring msg = get_clock_in_sec_string( t2-t1 );

   MYTRACE2(gCCDTrace,"Updating maximum frame took : " << msg );
}

void CCDPipeline::InitRandom()
{
	if(m_pDarkFrame)
		m_pDarkFrame->InitRandomDarkFrame();		
}

void CCDPipeline::ClearEvents( LONG_T FrameIndex )
{
	m_pCCD[FrameIndex].ClearState();
}
void CCDPipeline::ClearState()
{
	m_szOutputDir = "";
	for(long i=0;i<m_nCCD;i++){
		m_pCCD[i].ClearState();
	}
}

void CCDPipeline::DumpPipeline(long num,long start)
{
	if(m_Count==m_nPipelineSize){
		clock_t t1 = clock();
		mystring szFileName,szFullFileName;
		if(m_szOutputDir.empty()){
			m_szOutputDir << gCCDParams.GetOutputDir()  << "/" << get_date_time_string();
		}
		szFileName = m_szOutputDir;
		MyFile::CreateDir(szFileName.c_str());
		if(start==-1)
			start = head;
		long i=start;
		LONGLONG_T idx = m_FrameCounter-m_Count+1;
		while(i!=tail){
			cCCD& Frame = m_pCCD[i];
			szFullFileName = szFileName;
			szFullFileName << "/Frame_" << idx;
			Frame.WriteFrameToFile( szFullFileName.c_str() ); 			
			i--;
			if(i<0)
				i=m_nPipelineSize-1;			
			idx++;
		}
		szFullFileName = szFileName;
		szFullFileName << "/";
		if(m_pCCD[tail].GetInteresting())
			szFullFileName << "Hit_"; 
      szFullFileName << "Frame_" << idx;
		m_pCCD[tail].WriteFrameToFile( szFullFileName.c_str(), TRUE );		
		clock_t t2 = clock();
	   mystring msg = get_clock_in_sec_string( t2-t1 );

	   MYTRACE1(gCCDTrace,"Dumping to file of event frame took : " << msg );
	}
}

void CCDPipeline::Init_MC()
{
	if(gCCDParams.GetMC()){
		const char* szListFileName = GetParam("CCD_FULL_FRAMES_LIST");
		if(m_pList)
			delete m_pList;
		m_pList = new CListFile( szListFileName );
		m_FrameIdx = 0;		
		m_szLastFITSFile = "";
	}
	
}

BOOL_T CCDPipeline::EndFileExists()
{
	if( MyFile::DoesFileExist( gCCDParams.m_szEndFilePath.c_str() ) ){
		printf("End file %s found\n",gCCDParams.m_szEndFilePath.c_str());
		return TRUE;
	}	
	return FALSE;
}

BOOL_T CCDPipeline::ReReadListFile()
{
	if(gCCDParams.GetMC()){
		my_printf_now("re-reading list file ...");

		int start_count=0;

		if(!m_pList){
			const char* szListFileName = GetParam("CCD_FULL_FRAMES_LIST");
			m_pList = new CListFile( szListFileName );
		}else{
			start_count = m_pList->GetCount();
			time_t start_time = get_dttm();
			
			while( get_dttm()<(start_time+gCCDParams.m_nFrameListTimeout) && 
					 m_pList->GetCount()==start_count ){

				// if( EndFileExists() && get_dttm()>(start_time+1800) ){
				if( EndFileExists() ){
					// but break minimume 30 min after file found - to be sure all copied
					printf("night analysis already finished, EndFile found\n");
					break;
				}
			
				m_pList->ReadListFile( m_pList->GetFileName() );
				if( m_pList->GetCount()==start_count ){
					printf_now2("waiting %d sec for new frames to come ...\n",gCCDParams.m_nWaitTime);
					sleep( gCCDParams.m_nWaitTime );
				}
			}
			if( m_pList->GetCount()==start_count ){
				printf_now2("Timeout of %d sec reached, no new frames found\n",gCCDParams.m_nFrameListTimeout);
			}else{
				my_printf_now("waiting for all files to be copied ...\n");
				CheckNewFiles();
			}
		}

		printf("setting current frame index ...\n");
		for(int i=0;i<m_pList->GetCount();i++){
			if( strstr( m_szLastFITSFile.c_str(), (m_pList->GetListTable())[i].c_str() ) ){		
				m_FrameIdx = i+1;
				printf("After re-reading frame idx set to %d\n");
				if( m_FrameIdx < m_pList->GetCount() ){
					printf("After re-reading frame idx set on frame %s\n",(m_pList->GetListTable())[m_FrameIdx].c_str() );
				}
				break;
			}
		}

		return (start_count!=m_pList->GetCount());
	}
	return FALSE;
}

void CCDPipeline::CheckNewFiles()
{
	BOOL_T bOK=FALSE;
	time_t start_time = get_dttm();

	while( !bOK && get_dttm()<(start_time+gCCDParams.m_nFrameListTimeout) ){
		BOOL_T bAll=TRUE;
		for(int i=m_FrameIdx;i<m_pList->GetCount();i++){
			mystring szFile;
         szFile << gCCDParams.m_szSampleFramesDir << "/" << (m_pList->GetListTable())[i]; 

			if( !MyFile::DoesFileExist( szFile.c_str() ) ){
				printf_now2("File %s is missing\n",szFile.c_str() );
				bAll=FALSE;
				break;
			}
		}	
		if(!bAll){
			my_printf_now("Not all new files already copied, waiting 20 sec\n");
			sleep(20);
		}else{
			my_printf_now("All files found, continuing\n");
			bOK=TRUE;
		}
	}
}

BOOL_T CCDPipeline::IsLastFrame()
{
	return ( m_FrameIdx >= m_pList->GetCount() );
}

void CCDPipeline::SkipFramesN( int toSkip )
{
	m_FrameIdx = toSkip;
}

const char* CCDPipeline::GetNextFrameName()
{
	if(m_pList){
		return (m_pList->GetListTable())[m_FrameIdx].c_str();		
	}
	return "";
}



BOOL_T CCDPipeline::GetNextFrame_Real( mystring& szFileName, BOOL_T bClaimForNext )
{
	return FALSE;
}

void CCDPipeline::NewFrameAnalyseEnd()
{
	if(!gCCDParams.GetMC()){
	}
}

void CCDPipeline::ForceShiftsDetermination()
{
	printf("forcing shifts determination\n");
	gCCDErrLog.DumpToFile1( this, "FORCE_NEW_SHIFTS_WRN", "Forcing shifts determination" );
	

	m_PipelineCfg.m_nAutoShiftsCalc = abs( m_PipelineCfg.m_nAutoShiftsCalc );
	MoveCurrentShiftsFile();
	m_StarSpyTab[0].ResetTracedStar();
}

void CCDPipeline::GetEventsLogFileName( mystring& szFile, const char* szBaseName,
													 BOOL_T bFullPath )
{
	char szTmp[128];		
	if( strstr( szBaseName, "%d" )){
		sprintf(szTmp,szBaseName,m_PipelineIndex);
	}else{
		strcpy( szTmp,szBaseName );
	}
	szFile = "";
	if( bFullPath )
		szFile << gCCDParams.GetOutputDir() << "/";
	szFile << szTmp;
}

void CCDPipeline::GetPosition( double& ra, double& dec, double& azim,
                  	          double& alt, int& obsMode )
{
	eObservationMode_T _obsMode;
	CCDProcState* pProcState = &(m_CameraStates[0]);
	pProcState->m_pAstroForms->GetObsCoo( ra, dec, azim, alt, _obsMode );
	obsMode = (int)_obsMode;
	printf("###################################################\n");
	printf("GetPosition returns :\n");
	printf("(ra,dec)=(%.8f,%8f)   (azim,alt)=(%.8f,%8f)   obsMode=%d\n",
			 AstroAngle::rad2deg(ra),AstroAngle::rad2deg(dec),
			 AstroAngle::rad2deg(azim),AstroAngle::rad2deg(alt),obsMode);
	printf("###################################################\n");
}


void CCDPipeline::UpdateObsCoordinates( double ra, double dec, 
													 double azim, double alt,
													 eObservationMode_T& obsMode )
{
	CCDProcState* pProcState = &(m_CameraStates[0]);

	pProcState->m_pAstroForms->ChangeObsCoordinates( ra, dec, azim, alt, obsMode );
}


void CCDPipeline::UpdateObsCoordinates( double ra, double dec, double azim, double alt )
{
	CCDProcState* pProcState = &(m_CameraStates[0]);

	pProcState->m_pAstroForms->ChangeObsCoordinates( ra, dec, azim, alt );


	CCDLog log("%.4f\t%.4f\t%.4f\t%.4f\t%s\n","RA(h)\tDEC(deg)\tAZIM(deg)\tALT(deg)\tMODE\n");

	const char* szMode = GetObsModeDesc( pProcState->m_pAstroForms->GetObsMode() );
	log.DumpToFile( this, "Coord", "coordupdatelog", 
						 AstroAngle::rad2hours(ra), AstroAngle::rad2deg(dec), 
						 AstroAngle::rad2deg(azim), AstroAngle::rad2deg(alt), szMode );
}

BOOL_T CCDPipeline::CompareObsCoordinates( double ra, double dec, double azim, double alt, 
														eObservationMode_T obsMode )
{
	eObservationMode_T _obsMode;
	double _ra,_dec,_azim,_alt;
	CCDProcState* pProcState = &(m_CameraStates[0]);
	// pProcState->m_pAstroForms->GetObsCoordinates( _ra, _dec, _azim, _alt, _obsMode );
	if( !m_pAsasTransform )
		return TRUE;

	if( obsMode!=((eObservationMode_T)m_WorkingMode.m_MountObsMode) )
		return TRUE;
	if( obsMode==eEarthMovingMode ){
		if( !compare_double( ra, m_WorkingMode.m_MountRA_InRad ) )
			return TRUE;
		if( !compare_double( dec, m_WorkingMode.m_MountDec_InRad ) )
			return TRUE;
	}
	if( obsMode==eNoMovingMode ){
		if( !compare_double( azim, m_WorkingMode.m_MountAzim_InRad ) )
			return TRUE;
		if( !compare_double( alt, m_WorkingMode.m_MountAlt_InRad ) )
			return TRUE;
	}

	return FALSE;
}


void CCDPipeline::ChangeObsCoordinates( double ra, double dec, double azim, double alt, 
													 eObservationMode_T obsMode, BOOL_T bAstroReq )
{
	CCDProcState* pProcState = &(m_CameraStates[0]);

	// TODO : change DecObs,AzimObs etc ... :

	if( m_CameraStates[0].GetUtTime()<=0 ){
		SetCurrentFrameTimes( get_dttm() );
	}

	double ra_set = ra;
   double dec_set = dec;
   double azim_set = azim;
   double alt_set = alt;

	if( bAstroReq || !IsAstrometryOK() ){
		if( obsMode == eNoMovingMode ){
			m_PipelineCfg.m_HorAzimuth = azim;
			// m_PipelineCfg.m_HorAzimuth += m_PipelineCfg.m_HorAzimCorr;
			m_PipelineCfg.m_HorAltitude = alt;
			// m_PipelineCfg.m_HorAltitude  += m_PipelineCfg.m_HorAltCorr;

			double ha;
			pProcState->m_pAstroForms->calculateEquatorialCoordinates( m_PipelineCfg.m_HorAzimuth, m_PipelineCfg.m_HorAltitude,
													  ha, m_PipelineCfg.m_DecObs, m_PipelineCfg.m_RAObs );							  
		}else{
			m_PipelineCfg.m_RAObs = ra;
			// m_PipelineCfg.m_RAObs += m_PipelineCfg.m_RACorr;
			m_PipelineCfg.m_DecObs = dec;
			// m_PipelineCfg.m_DecObs += m_PipelineCfg.m_DecCorr;

			pProcState->m_pAstroForms->calculateHorizontalCoordinatesFromEq( m_PipelineCfg.m_DecObs, m_PipelineCfg.m_RAObs, 
															  m_PipelineCfg.m_HorAltitude, m_PipelineCfg.m_HorAzimuth );
		}

		ra_set = m_PipelineCfg.m_RAObs;
		dec_set = m_PipelineCfg.m_DecObs;
		azim_set = m_PipelineCfg.m_HorAzimuth;
		alt_set = m_PipelineCfg.m_HorAltitude;			
		pProcState->m_pAstroForms->ChangeObsCoordinates( ra_set, dec_set, 
																		 azim_set, alt_set, obsMode );
	}
	

	CCDLog log("%.4f\t%.4f\t%.4f\t%.4f\t%s\n","RA(h)\tDEC(deg)\tAZIM(deg)\tALT(deg)\tMODE\n");

	const char* szMode = GetObsModeDesc( obsMode );
	log.DumpToFile( this, "Coord", "coordlog", 
                   AstroAngle::rad2hours(ra_set), 
                   AstroAngle::rad2deg(dec_set),
						 AstroAngle::rad2deg(azim_set), 
						 AstroAngle::rad2deg(alt_set), szMode );

	m_bCoordInitialized = TRUE;

	if( bAstroReq ){
// changed in 20060801 :
		m_pAsasTransform->ResetAstroOK();
//		m_pAsasTransform->m_bTransformOK = FALSE;
	}

	// saving coordinates from mount :
	SetMountCoord( ra_set, dec_set, azim_set, alt_set, obsMode );

	// for simulator purposes only :
	// CCDController::SetRADEC( AstroAngle::rad2hours(ra_set), AstroAngle::rad2deg(dec_set) );
}

void CCDPipeline::SetMountCoord( double ra, double dec, double azim, double alt,
 	                              eObservationMode_T obsMode )
{
	time_t mountTimeStamp = get_dttm();
	printf("Changing mount coordinates to (ra,dec)=(%.2f,%.2f) , (azim,alt)=(%.2f,%.2f), obsMode=%d\n",
				AstroAngle::rad2deg(ra),AstroAngle::rad2deg(dec),
				AstroAngle::rad2deg(azim),AstroAngle::rad2deg(alt),obsMode);fflush(0);


	// saving coordinates from mount :
	m_WorkingMode.m_MountRA_InRad = ra;
   m_WorkingMode.m_MountDec_InRad = dec;
   m_WorkingMode.m_MountAzim_InRad = azim;
   m_WorkingMode.m_MountAlt_InRad = alt;
   m_WorkingMode.m_MountObsMode = obsMode;
	m_WorkingMode.m_MountInfoTimeStamp = mountTimeStamp;
}



void CCDPipeline::ForceAstrometryOnNext()
{
	vector<CCDPipeline*>::iterator i;
   for(i=m_PipelineList.begin();i!=m_PipelineList.end();i++){
		if( !(*i)->m_bAstrometryRunning ){
			(*i)->m_bForceAstroNow = TRUE;
			(*i)->m_bForceSynchroMode = TRUE;
			// (*i)->m_pAsasTransform->m_bTransformOK = FALSE;			
		}
	}
}


BOOL_T CCDPipeline::ExecResetRequest()
{
	vector<CCDPipeline*>::iterator i;
   for(i=m_PipelineList.begin();i!=m_PipelineList.end();i++){
		(*i)->InitPipeline();
	}
	return TRUE;
}

int GlobalCCDPipeline_EmergencyEnd( int bDoExit )
{
	return CCDPipeline::EmergencyEnd( bDoExit );
}

int CCDPipeline::EmergencyEnd( int bDoExit /*=0*/ )
{
	printf("EmergencyEnd : entering function\n");fflush(stdout);

	// closing pipelines - dumping events 
	vector<CCDPipeline*>::iterator i;
	printf("Ending analysis ...");fflush(0);
   for(i=m_PipelineList.begin();i!=m_PipelineList.end();i++){
		// ending analysis 
		(*i)->End();		
	}	

	printf("EmergencyEnd : closing devices ...\n");fflush(stdout);
	// closing cameras :

	printf("EmergencyEnd : do exit ?\n");fflush(stdout);
	if( bDoExit ){
		printf("EmergencyEnd : exiting now\n");fflush(stdout);
		exit(-1);
	}
	
	printf("EmergencyEnd : finished\n");fflush(stdout);

	return 1;
}

void CCDPipeline::SetAstroFlagNow()
{
	vector<CCDPipeline*>::iterator i;
   for(i=m_PipelineList.begin();i!=m_PipelineList.end();i++){
		(*i)->m_bForceAstroNow = TRUE;
	}
}

void CCDPipeline::ResetAstro()
{
	_STDOUT_TRACE_("RESETING ASTROMETRY\n");
	vector<CCDPipeline*>::iterator i;
   for(i=m_PipelineList.begin();i!=m_PipelineList.end();i++){
// changed on 20060801 :	
		(*i)->m_pAsasTransform->ResetAstroOK();
//		(*i)->m_pAsasTransform->m_bTransformOK	= FALSE;

		(*i)->m_bForceSynchroMode = TRUE;
		(*i)->m_nGoodAstrometryCount = 0; // 20060604
		(*i)->m_pAsasTransform->m_nFailedAstrometryCount = 0;
	}
}

BOOL_T CCDPipeline::IsAstroRunning()
{
	if( m_PipelineList.size()>0 ){
		return (m_PipelineList[0])->m_bAstrometryRunning;
	}	
	return FALSE;
}

time_t CCDPipeline::GetCurrentFrameTime()
{
	return m_FrameUnixTime;
}

time_t CCDPipeline::GetAstrometryTime()
{
	if( m_PipelineList.size()>0 ){
		if( (m_PipelineList[0])->m_pAsasTransform ){
			return (m_PipelineList[0])->m_pAsasTransform->m_LastGoodAstroTime;
		}
	}	
	return 0;
}

BOOL_T CCDPipeline::GetAstroOK()
{
	if( m_PipelineList.size()>0 ){
		return (m_PipelineList[0])->IsAstrometryOK();
	}	
	return FALSE;
}

int CCDPipeline::IsScan()
{
	return m_WorkingMode.m_bScanMode;
}

int CCDPipeline::IsCooling()
{
	return FALSE;
}

int CCDPipeline::GetDarkCount()
{
	return gCCDParams.m_nDarksToBeTaken;
}

BOOL_T CCDPipeline::IsDomeOpened()
{
	return ( m_WorkingMode.m_eDomeStatus != eDomeClosed );
}


void CCDPipeline::GetCurrentFrameCoo( CCDMatrix& matrix,
												  time_t& ut_time,
												  double& ra, double& dec, double& alt, 
												  double& azim, double& ha )
{
	ut_time = (int)matrix.getObsTime();
	eObservationMode_T obsMode;	
	m_CameraStates[0].GetObsCoo( ut_time, ra, dec, alt, azim, ha, obsMode );
}

void CCDPipeline::GetCurrentFrameCoo( CCDMatrix& matrix,
												  time_t& ut_time,
												  double& ra, double& dec, double& alt, 
												  double& azim, double& ha,
												  eObservationMode_T& obsMode	 )
{
	ut_time = (int)matrix.getObsTime();
	m_CameraStates[0].GetObsCoo( ut_time, ra, dec, alt, azim, ha, obsMode );
}

void CCDPipeline::SetCurrentFrameTimes( CCDMatrix& matrix, int i )
{
	time_t ut_time = 0;
	if( gCCDParams.m_bReadFirstLevelInfoFromLog ){
		if( gCCDParams.m_bFromLogFile_WithFrames ){
			ut_time = (int)matrix.getObsTime();
		}else{
			ut_time = m_FrameUnixTime;
		}
	}else{
		ut_time = (int)matrix.getObsTime();
	}
	m_CameraStates[i].SetFrameDateTime( ut_time, this );
}

void CCDPipeline::SetCurrentFrameTimes()
{
	cCCD& currentFrame = GetCurrent();

	Assert( m_CameraStates.size()==currentFrame.GetCount(), "Object m_CameraStates not initialized correctly !!!");

	for(int i=0;i<currentFrame.GetCount();i++){
		CCDMatrix& matrix = currentFrame[i];
		SetCurrentFrameTimes( matrix, i );
	}
}

void CCDPipeline::SetCurrentFrameTimes( time_t ut_time )
{
	for(int i=0;i<m_CameraStates.size();i++){
		m_CameraStates[i].SetFrameDateTime( ut_time, this );
	}
}

BOOL_T CCDPipeline::GetNextFrameRaw( mystring* szFileNames, BOOL_T bUseCache, 
												 BOOL_T bClaimForNext )
{
	BOOL_T bRet=FALSE;
	time_t get_start=get_dttm();
	if(gCCDParams.GetMC()){
		if( gCCDParams.m_bReadFirstLevelInfoFromLog ){
			if( gCCDParams.m_bFromLogFile_WithFrames ){
				// in this case must read events from log files now :
				ReadEvents_FromLogFile();
			}
		}

		if( gCCDParams.m_bReadFirstLevelInfoFromLog && !gCCDParams.m_bFromLogFile_WithFrames ){
			bRet = GetNextFrame_FromLogFile( *szFileNames, bUseCache );
		}else{
			BOOL_T bOK=FALSE;
			while( !bOK ){
				bRet = GetNextFrame_MC( *szFileNames, bUseCache );
				bOK=TRUE;
				if( bRet ){
					if( strlen( gCCDParams.m_szSkipIfObjectMatches.c_str() )>0 ){
						const char* szOBJECT = GetNext()[0].getKey( OBJECT );
						if( szOBJECT && szOBJECT[0] ){
							if( strstr( szOBJECT, gCCDParams.m_szSkipIfObjectMatches.c_str() )){
								printf("OBJECT %s required to be skiped - check value of parameter CCD_SKIP_IF_OBJECT_MATCHES ( = %s )\n",szOBJECT,gCCDParams.m_szSkipIfObjectMatches.c_str());
								bOK = FALSE;
							}
						}
					}
				}
			}
		}
	}else{
		// bRet = GetNextFrame_Real( *szFileNames, bClaimForNext );
	}
	time_t get_end=get_dttm();
   m_PureGetNextTime = (get_end-get_start);
	 if(!gCCDParams.m_bCCDDouble){
      _TRACE_PRINTF_2("\n\nPURE GET TOOK : %d\n",(get_end-get_start));
   }
	
	return bRet;
}

BOOL_T CCDPipeline::GetNextFrameRawAcc( mystring* szFileNames, BOOL_T bUseCache )
{
	BOOL_T bRet = GetNextFrameRaw( szFileNames, bUseCache );
	if( bRet ){
		AcceptNew();
                                                                                
	   // set times of current frame to AstroForumla calculating object it is
   	// done one after new frame is obtained from Driver ( file ) :
	   SetCurrentFrameTimes();
	}
	return bRet;
}

void CCDPipeline::CleanBeforeNewFramePre()
{
	// cleaning flag that  claimed for next frame :
//	if( GetCCDInterface() ){
//		GetCCDInterface()->m_bReqForNextCalled = FALSE; // [NEW]
//	}

	// cleaning events lists :
   m_EventsOnSumedFrame.clear();
   m_EventsFromCompareToOld.clear();
	m_BrightenOnSumList.clear();

	// flag that events were dumped 
	m_bFinalDumped = FALSE;

	// cleaning confirmed events list :
	ClearConfirmedList();

	// clear anty-coic ( cosmic ) event list :
	m_AntyCoicEvents.clear();

	// clear brighten list
	m_BrightenList.clear();
}

CMyMutex gAsynchroMutex;
BOOL_T CCDPipeline::GetNextFrame( mystring* szFileNames, BOOL_T bUseCache, 
											 BOOL_T bHandleReq, BOOL_T bClaimForNext,
											 BOOL_T bOnlyCallGet, BOOL_T bDoDumpAtLast )
{	
	if( !gCCDParams.m_bExecAstroSynchroMode ){
		printf("cam%d - Astrometry in asynchro mode , checking if fresh done ...\n",m_PipelineIndex);fflush(0);
		if( m_bAstrometryIsReady ){
			double ra_old=0.00,dec_old=0.00,ra_new=0.00,dec_new=0.00;
			if( m_pAsasTransform ){
				ra_old = m_pAsasTransform->ra;
				dec_old = m_pAsasTransform->dec;
			}
			if( m_pAsasTransformAsynchro ){
				ra_new = m_pAsasTransformAsynchro->ra;
				dec_new = m_pAsasTransformAsynchro->dec;
			}
			printf("cam%d - New astrometry detected, will update now (%.2f,%.2f) -> (%.2f,%.2f)...\n",m_PipelineIndex,ra_old,dec_old,ra_new,dec_new);
			
			double max_acceptable_diff=10.00; // [deg] 
			if( fabs(ra_old-ra_new)>(max_acceptable_diff/15.00) || fabs(dec_old-dec_new)> max_acceptable_diff ){
				printf("WARNING : overwritting old transform with much different position\n");
			}

			memcpy( m_pAsasTransform, m_pAsasTransformAsynchro, sizeof(CCDAsasTransform) );
			printf("cam%d - Astrometry updated according to asynchro result, flag cleaned (astro_time=%d)\n",m_PipelineIndex,m_pAsasTransform->m_TransformUtTime);

			// new astrometry maybe performed in asynchro mode
			m_bAstrometryIsReady = FALSE;
		}else{
			int max_allowed_time=600;
			printf("cam%d - No fresh astrometry detected\n",m_PipelineIndex);
			if( ( m_FrameUnixTime - m_pAsasTransform->m_TransformUtTime ) > max_allowed_time ){
				printf("AstroAsynchroMode : ASAS-Transformation on camera %d is older then %d sec, reseting AstroOK flag\n",m_PipelineIndex,max_allowed_time);  
// changed on 20060801 :
				m_pAsasTransform->ResetAstroOK();
//				m_pAsasTransform->m_bTransformOK = FALSE;
			}else{
				printf("AstroAsynchroMode : astrometry fresher then %d sec, ok, astroOK=%d\n",max_allowed_time,m_pAsasTransform->m_bTransformOK);
			}
		}		
	}

	CleanBeforeNewFramePre();


	if( m_PipelineIndex==0 ){
		// only at first camera :
		if( gCCDParams.m_nExitOnToManyErrors>0 ){
			CheckIfExitOnError();
		}
	}


	if( m_bRestartOnNext ){
		printf_now2("RESTART OF PIPELINE %d required due to big postion change or change in observation mode\n",m_PipelineIndex);
		RestartPipeline( FALSE );
	}

	m_bRestartOnNext = FALSE;
	BOOL_T bRet = FALSE;
	if( szFileNames ){
		(*szFileNames) = "";
	}

	// protect global resources with lock - nessesary since GetNextFrame 
	// is called parrallely for both CCD in CCDDouble :
	gAsynchroMutex.Lock();
	
	// only at first frame :
	ExecPostInitActions();
	

	// before each frame check requests from PI-SYSTEM Manager:
	if( !gCCDParams.m_bCCDDouble ){
		// handle requests but only in single CCD mode , otherwise 
		// it is called from CCDDouble::Run :
		if(!HandleRequests( bHandleReq )){
			// exit request !!!
			gAsynchroMutex.UnLock();
			return FALSE;
		}
	}

	gAsynchroMutex.UnLock();	

	bRet = GetNextFrameRaw( szFileNames, bUseCache, bClaimForNext );

	gAsynchroMutex.Lock();

	// in real analysis flip is done in driver level now !!!
	// but in case of MonteCarlo it is here :
	if( gCCDParams.GetMC() && m_PipelineCfg.m_bDriverReverseImage!=eReverseImageNone ){
		GetNext()[0].FlipImage( m_PipelineCfg.m_bDriverReverseImage );
	}


	if(!bRet){
		// last frame probably , dumping rest of events :
		if( bDoDumpAtLast ){
			DumpFoundEventsLog( TRUE );
		}
		gAsynchroMutex.UnLock();

		return FALSE;
	}


	// sets also m_FrameUnixTime 
	AcceptNew();
	// m_FramesList.AddFrame( szFileNames->c_str(), m_FrameUnixTime );

	// set times of current frame to AstroForumla calculating object it is 
	// done one after new frame is obtained from Driver ( file ) :
	SetCurrentFrameTimes();

	if( IsCollectionMode() || bOnlyCallGet ){
		printf_now("Only get called exiting GetNextFrame\n");
		// in collection mode do not perform any further analysis :

		if( IsCollectionMode() ){
			printf("Collection mode - doing astro/photo\n");fflush(stdout); 
			if( gCCDParams.m_bDoAstrometryInTakeNMode ){
				printf("Executing astrometry in TakeN mode enabled ...\n");
				m_bForceAstroNow = TRUE;
				// m_pAsasTransform->m_bTransformOK = FALSE;
				UpdateAstrometry( szFileNames, TRUE );
			}else{
				printf("Astrometry in TakeN mode disabled ...\n");
				if( gCCDParams.m_bDoPhotoInTakeNMode && szFileNames ){
					if( gCCDParams.m_nDarksToBeTaken <= 0 ){
						RunPhotometry( szFileNames->c_str() , GetCurrent()[0], FALSE );
						UpdateCurrentFITSFile_PhotoOnly( szFileNames->c_str() );
					}else{
						printf("In dark collecting mode, fast photometry not performed in TakeNMode\n");fflush(stdout);
					}
				}
			}
		}else{
			printf("Only-call get - no astro/photo actions\n");
		}
		gAsynchroMutex.UnLock();
		return TRUE;
	}

	// checking frame quality if needed :
	if( gCCDParams.m_bCheckFrame ){
		if( !GetCurrent()[0].CheckFrame( m_LastBadLineY ) ){
			const char* dImage = GetCurrent()[0].getKeyValue( DIMAGE );
			int frame_no = GetDayFrameCounter();
			if( dImage && dImage[0] ){
				frame_no = atol( dImage );
			}
			printf("Error detected in frame %d, bad line at y=%d\n",frame_no,m_LastBadLineY);
			m_LastBadFrame = frame_no;		
			
			if( gCCDParams.m_bRepairFrame && gCCDParams.GetMC() ){
				// repairing frames in MC mode here, in real mode in CCDController::Run
				GetCurrent()[0].RepairFrame( m_LastBadLineY );
			}
		}
	}


	clock_t t0 = clock();
	PrepareNewFrame();
	clock_t t1 = clock();
   clock_t dt = (t1-t0);
	mystring msg = get_clock_in_sec_string( dt );
   _TRACE_PRINTF_5("Preperation of new frame, took %s sec\n",msg.c_str());

	if( gCCDParams.GetMC() ){
		if( !gCCDParams.m_bReadFirstLevelInfoFromLog ){
		   // m_DayFrameCounter = m_FrameCounter+m_PipelineCfg.m_DayFramesCounter;
			const char* szDIMAGE = (GetCurrent()[0]).getKeyValue( DIMAGE );
			if( szDIMAGE && szDIMAGE[0] ){
				m_DayFrameCounter = atol( szDIMAGE );
			}else{
				m_DayFrameCounter = m_FrameCounter+m_PipelineCfg.m_DayFramesCounter;
			}
		}
	}


	GetCurrent().AdjustFrames();

	if( gCCDParams.GetMC() ){
		BOOL_T bOK=FALSE;
		if( gCCDParams.m_bUseAsasAstrometryFromFITS ){
			bOK = GetASASAstrometryFromFITS();
			if( bOK ){				
				_TRACE_PRINTF_0("CAM%d : Astrometry read from FITS\n",m_PipelineIndex);
			}else{
				_TRACE_PRINTF_0("CAM%d : Could not retrive astrometry from FITS header, using from cfg\n",m_PipelineIndex);
			}	
		}
		if( !bOK ){
			if( gCCDParams.m_bUseCoordFromFITS ){
				SetCoordFromFITS();
			}
		}
	}

	// auto shifts determination :
   UpdateAutoShiftStars();

	UpdateAstrometry( szFileNames );

	if(szFileNames )
		GetCurrent().SetFileNames( szFileNames->c_str() );

	// moved before astrometry - so in case of big shift astrometry is performed 
	// auto shifts determination :
	// UpdateAutoShiftStars();

	if(gCCDParams.m_bTraceOnAllFrames){
		UpdateTraceOnAll();
	}

	UpdateControlStar();

	gAsynchroMutex.UnLock();

	RemoveOldRejectIfMoreTracks();
	
	if( !gCCDParams.m_bParamsTest || GetCount()<GetPipelineSize() ){
		// in case it is Nparamstest / Nparamstest2C, this update of 
		// averaged frame is performed after sample is put 
		// or put is not performed yet :
		UpdateSumOfFrames( szFileNames );
	}

	CleanBeforeNewFramePost();

	// check if there was coordinates change :
	CheckCoordChange();


	// saving current astro transformation to Sav structure :
	memcpy( m_pAsasTransformSav, m_pAsasTransform, sizeof(CCDAsasTransform));

	return bRet;
}

BOOL_T CCDPipeline::CheckCoordChange()
{
	if( GetFrameIndex() > 1 ){
		printf("Checking coordinates change ...\n");fflush(0);
		if( CheckIfChangeCoordBig( m_pAsasTransformSav->ra, m_pAsasTransformSav->dec,
                              AstroAngle::rad2deg( m_pAsasTransformSav->m_AzimInRad ),
                              AstroAngle::rad2deg( m_pAsasTransformSav->m_AltInRad ),
                              m_pAsasTransform ) ){
   		printf("\n\n!!! BIG change of coordinates detected after astrometry !!!!");
	      printf("will restart pipeline NOW !\n\n");
			return TRUE;
	   }
		printf("Check of coord change DONE\n");fflush(0);
	}else{
		printf("First frame check of coordinates change ignored\n");
	}

	return FALSE;
}

void CCDPipeline::CleanBeforeNewFramePost()
{
}


Table2D<BIG_ELEM_TYPE>* CCDPipeline::GetPrevSum( const char* szFileName )
{
	return m_pOldSumOfSeveral;
}

Table2D<BIG_ELEM_TYPE>* CCDPipeline::GetPrevLaplace( const char* szFileName )
{
	return m_pOldSumOfSeveral;
}

void CCDPipeline::UpdateSumOfFrames( mystring* szFileNames )
{
	if( m_pCurrSumOfSeveral ){
		time_t ut_time = m_FrameUnixTime;
		int prev_frame=(m_FrameCounter-1);

		if( gCCDParams.m_bDoSumOfPrevNFrames ){
			// in case of sum of N frames analysis or n case comparison to old frame 
			// is enabled and its time to compare :
			// m_nFramesCollectedForSum++;

			// first check if already collected required number of frames 
			// in case it was ( on previous frame ), set prev sum frame and
			// reset current sum frame :
   	   if( m_nFramesCollectedForSum==gCCDParams.m_bKeepSumOfPrevNFrames ){
				// update of stored prev frame - after new sum is analysied :				
				printf("AVG_ANAL : updating prev sum frame, reseting current sum frame\n");

				(*m_pPrevSumOfSeveral) = (*m_pCurrSumOfSeveral);
				(*m_pPrevSumOfSeveralLaplace) = (*m_pCurrSumOfSeveralLaplace);
				if( m_pOldSumOfSeveral ){
					if( gCCDParams.m_nCompareToOldFreqInSec>0 && ( ut_time - m_PrevSumOfFramesTime ) > gCCDParams.m_nCompareToOldFreqInSec ){
						(*m_pOldSumOfSeveral) = (*m_pCurrSumOfSeveral);
						(*m_pOldSumOfSeveralLaplace) = (*m_pCurrSumOfSeveralLaplace);
						m_szPrevSumOfNFileName = m_szCurrSumOfNFileName;
	               m_szCurrSumOfNFileName = "";
					}
				}
				m_pCurrSumOfSeveral->SetData( 0 );	

				
				m_PrevSumOfFramesTime = ut_time;
				m_nFramesCollectedForSum = 0;		

				if( !gCCDParams.m_bDoSumOfPrevNFrames ){
					// in case sum is taken every N sec ( comparison to old frame )
					// no need to add new frame now :
					return;
				}
			}

			if( gCCDParams.m_nCompareToOldFreqInSec>0 ){			
				if( ( ut_time - m_PrevSumOfFramesTime ) > gCCDParams.m_nCompareToOldFreqInSec ){
					if( m_nFramesCollectedForSum==1 ){
						printf("AVG_ANAL : Previous sum frame taken %s sec ago ( interval=%d sec ), starting new average now\n",( ut_time - m_PrevSumOfFramesTime ),gCCDParams.m_nCompareToOldFreqInSec);
					}else{			
						printf("AVG_ANAL : Updating average (%d frame added)\n",m_nFramesCollectedForSum );
					}
				}
			}



			CCDUtil::SumFrames( *m_pCurrSumOfSeveral, GetCurrent()[0] );
			m_nFramesCollectedForSum++;

			if( m_nFramesCollectedForSum == gCCDParams.m_bKeepSumOfPrevNFrames ){
				printf("AVG_ANAL : Calculating average of %d frames ...",gCCDParams.m_bKeepSumOfPrevNFrames);

				printf("Normalizing sum frame by factor : %d\n",gCCDParams.m_bKeepSumOfPrevNFrames);
				double factor = (1.00 / gCCDParams.m_bKeepSumOfPrevNFrames);
				m_pCurrSumOfSeveral->MultiplyByConst( factor );

				// calculating laplace of new average :
				m_pCurrSumOfSeveral->Laplace( m_pCurrSumOfSeveralLaplace->get_data_buffer_fast(),
														gCCDParams.m_eLaplaceType, 
														5, 5, m_pCurrSumOfSeveral->GetXSize()-5,
														m_pCurrSumOfSeveral->GetYSize()-5 );
				m_nTotalCollectedSumFrames++;

				if( gCCDParams.m_bAnalyzeSumOfPrevNFrames ){
   			   SaveSumFrameEventsOnNext();

					if( gCCDParams.m_bCheckForSUPERNEW ){
						SaveSNSumFrameEventsOnNext();
					}
			   }

				printf("OK\n");
				if( gCCDParams.m_bSaveSumOfNFrames ){
					// saving sum from previous frame - OPTIONAL :

					CCDMatrix tmp( m_pCurrSumOfSeveral->GetXSize(), m_pCurrSumOfSeveral->GetYSize() );
					CCDUtil::CopyFromBigElemToElem( *m_pCurrSumOfSeveral, tmp );
					mystring szAverageFrame;

					szAverageFrame  << "aver/" << getfname( *szFileNames ) << "_aver.fit";
					tmp.GetKeyTab() = (GetCurrent()[0]).GetKeyTab();
					printf("AVG_ANAL : saving average to file : %s ",szAverageFrame.c_str() );
					tmp.WriteToFITSFile( szAverageFrame.c_str() );	
					printf("OK\n");


					if( gCCDParams.m_nCompareToOldFreqInSec>0 && ( ut_time - m_PrevSumOfFramesTime ) > gCCDParams.m_nCompareToOldFreqInSec ){				
						m_szCurrSumOfNFileName = szAverageFrame;

						// in this case save also laplace frame :
						//szAverageFrame = "";
						//szAverageFrame << "aver/" << getfname( *szFileNames ) << "_aver_lap.fit";
						//CCDUtil::WriteToFITSFile( *m_pCurrSumOfSeveralLaplace, szAverageFrame.c_str() );
					}
				}
			}
		}
	}
}

BOOL_T CCDPipeline::SaveAstrometry( mystring* szFileNames )
{
	CCDMatrix image( m_SizeX, m_SizeY );
	if( !image.ReadFITSFile( szFileNames->c_str() ) ){
		printf("could not read fits file : %s\n",szFileNames->c_str());
		return FALSE;
	}
	AddAstrometry( image.GetKeyTab(), m_pAsasTransform );

	image.WriteToFITSFile( szFileNames->c_str() );
	printf("Astrometry saved to fits file : %s\n",szFileNames->c_str() );

	return TRUE;
}

void CCDPipeline::UpdateAstrometry( mystring* szFileNames, BOOL_T bForceInAnyMode )
{
	if( gCCDParams.m_PostponeAstrometryTime ){
		if( m_NextAstrometryTime > 0 ){
			time_t current_time = get_dttm();
			if( current_time < m_NextAstrometryTime )
			{
				printf("Taking picture of GCN - astrometry postponed to unix_time=%d\n",m_NextAstrometryTime);
				return;
			}
		}
	}

	if( ( m_PipelineCfg.m_nAutoAstrometryFreq>0 || m_PipelineCfg.m_nAutoAstrometryFreqInSec>0 ) && 
		 m_PipelineCfg.m_bDoASASPhotAstr && 
		 !m_bAstrometryRunning && ( !IsCollectionMode() || bForceInAnyMode ) ){
		int nDayFrame = GetDayFrameCounter();

		// 20050126 added : && m_PipelineCfg.m_nAutoAstrometryFreq>0 
		if( ( (nDayFrame % m_PipelineCfg.m_nAutoAstrometryFreq)==0 && m_PipelineCfg.m_nAutoAstrometryFreq>0 ) // every N frame
			||  ( GetFrameCounter()==1 && gCCDParams.m_bDoASASAstroOnStartup && !IsAstrometryOK() ) // or on first frame 
			|| m_bForceAstroNow || m_nAstrometryRetry>0 || // or forced or retry
			( m_PipelineCfg.m_nAutoAstrometryFreqInSec>0 && 
			 (m_FrameUnixTime-m_pAsasTransform->m_TransformUtTime)>m_PipelineCfg.m_nAutoAstrometryFreqInSec) ){ // too old 

			if( m_PipelineCfg.m_nAutoAstrometryFreqInSec>0 && 
				 (m_FrameUnixTime-m_pAsasTransform->m_TransformUtTime)>m_PipelineCfg.m_nAutoAstrometryFreqInSec )
			{
				printf("ASAS transformation older then limit = %d sec, determining new ...\n",m_PipelineCfg.m_nAutoAstrometryFreqInSec);fflush(0);
				
			}

			// mystring szReductFile;
			// SaveReducedFrame( szReductFile );
			if( gCCDParams.GetMC() && szFileNames ){
				m_PipelineCfg.m_szLastReducedFrameName = (*szFileNames);
			}

			if( gCCDParams.m_bDoASASPhotAstr ){
				if( m_nAstrometryRetry<=0 ){
					m_nAstrometryRetry = gCCDParams.m_nASASAstrometryReTry;
				}
				if( gCCDParams.m_bExecAstroSynchroMode || m_bForceSynchroMode || 
					 gCCDParams.m_bDoAstrometryInTakeNMode ){
					printf("cam%d REQUEST FOR ASTROMETRY in SYNCHRO MODE (%d,%d,%d,%d,%d) :\n",m_PipelineIndex,
										gCCDParams.m_bExecAstroSynchroMode,m_bForceSynchroMode,gCCDParams.m_bDoAstrometryInTakeNMode,
										m_nAstrometryRetry,gCCDParams.m_nASASAstrometryReTry);

					// in some cases synchronous astrometry is FORCED 
					// - after position change
					// - when taking frames in SynchroMode - take_n_pictures_synchro

					// run astrometry in current thread and wait for end :

					if( DoAstrometry( this ) ){
						m_nAstrometryRetry =  0;


						// astrometry ok - save results to original fits file 
						// ( if required )
//ms 20070907 						
						if( (gCCDParams.GetMC() && gCCDParams.m_bSaveGoodAstro) || gCCDParams.m_bSaveAstroInTakeNMode ){
							printf("Astrometry succeeded in scan mode\n");
							printf("saving astrometry to fits header ...\n");
							if( !SaveAstrometry( szFileNames ) ){
								printf("could not save fits file : %s\n",szFileNames->c_str());
							}
						}else{
							printf("saving astrometry to fits header - NOT REQUIRED\n");
						}
					}else{
						m_nAstrometryRetry--;
						// if( gCCDParams.m_bUseGoodFromOther ){
						// 	UseAstrometryFromOther();
						// }
					}
					if( m_nAstrometryRetry<=0 )
						m_bForceSynchroMode = FALSE;
					if( !gCCDParams.m_bExecAstroSynchroMode && m_nAstrometryRetry>0 ){
						int retry_count = (gCCDParams.m_nASASAstrometryReTry-m_nAstrometryRetry);
						printf("ASTROMETRY : generally in asynchro mode retry = %d\n",retry_count);
						if( retry_count >= gCCDParams.m_nForceSynchroMaxCount ){
							printf("ASTROMETRY : disabling ForceSynchroMode in cam%d\n",m_PipelineIndex);
							m_bForceSynchroMode = FALSE;
						}		
					}
				}else{
					printf("cam%d REQUEST FOR ASTROMETRY in ASYNCHRO MODE :\n",m_PipelineIndex);

					if( !m_bAstrometryRunning ){
						// in case previous astrometry is not running start new :
//						if( runAstrometryLock.Release()<1 ){
//							printf("ERROR runAstrometryLock.Release FAILED !!!\n");
//						}
						m_nAstrometryRetry=0;	
//						printf("runAstrometryLock released\n");fflush(0);
					}else{
						printf("Another astrometry request is currently running, cannot execute astrometry now\n");
						printf("In case of problems send param change CCD_ASTROMETRY_SYNCHRO_MODE=1\n");
					}
				}
				m_bAstrometryExecuted = TRUE;
			}
			m_bForceAstroNow = FALSE;
		}
	}
}

BOOL_T CCDPipeline::CheckIfChangeCoordBig( double prev_ra, double prev_dec,
														 double prev_azim_in_deg, 
														 double prev_alt_in_deg,
														 CCDAsasTransform* pAsasTransform )
{
	double diff_to_do_astro=0.2;
	double limit_to_restart = gCCDParams.m_CoordChangeToRestart;
	eObservationMode_T obsMode  = GetObsMode();
	BOOL_T bRestartPipeline=FALSE;

	if( obsMode == eEarthMovingMode ){
		if( fabs( AstroAngle::hours2deg( pAsasTransform->ra ) - AstroAngle::hours2deg( prev_ra ) ) > diff_to_do_astro ){
			if( m_PipelineCfg.m_IgnoreCoordChangeOnCamera != m_PipelineIndex ){
				printf("Big change in RA - forcing astrometry (%.2f -> %.2f)\n",prev_ra,pAsasTransform->ra);fflush(0);
				if( fabs( AstroAngle::hours2deg( pAsasTransform->ra ) - AstroAngle::hours2deg( prev_ra ) ) > limit_to_restart ){
					bRestartPipeline=TRUE;
				}
			}else{
				printf("Big change ignored no pipeline restart\n");fflush(0);
			}
		}

		if( fabs( pAsasTransform->dec - prev_dec ) > diff_to_do_astro ){
			if( m_PipelineCfg.m_IgnoreCoordChangeOnCamera != m_PipelineIndex ){
				printf("Big change in DEC - forcing astrometry (%.2f -> %.2f)\n",prev_dec,pAsasTransform->dec);
				if( fabs( pAsasTransform->dec - prev_dec ) > limit_to_restart ){
					bRestartPipeline=TRUE;
				}
			}else{
				printf("Big change ignored no pipeline restart\n");fflush(0);
			}
		}
	}

	if( obsMode == eNoMovingMode ){
		if( fabs( AstroAngle::rad2deg( pAsasTransform->m_AzimInRad ) - prev_azim_in_deg ) > diff_to_do_astro ){
			if( m_PipelineCfg.m_IgnoreCoordChangeOnCamera != m_PipelineIndex ){
				printf("Big change in AZIM - forcing astrometry (%.2f -> %.2f)\n",prev_azim_in_deg,AstroAngle::rad2deg( pAsasTransform->m_AzimInRad ));
				if( fabs( AstroAngle::rad2deg( pAsasTransform->m_AzimInRad ) - prev_azim_in_deg ) > limit_to_restart ){
					bRestartPipeline=TRUE;
				}else{
					printf("Big change ignored no pipeline restart\n");fflush(0);
				}
			}
		}
		if( fabs( AstroAngle::rad2deg( pAsasTransform->m_AltInRad ) - prev_alt_in_deg ) > diff_to_do_astro ){
			if( m_PipelineCfg.m_IgnoreCoordChangeOnCamera != m_PipelineIndex ){
				printf("Big change in ALT - forcing astrometry (%.2f -> %.2f)\n",prev_alt_in_deg,AstroAngle::rad2deg( pAsasTransform->m_AltInRad ));
				if( fabs( AstroAngle::rad2deg( pAsasTransform->m_AltInRad ) - prev_alt_in_deg ) > limit_to_restart ){
					bRestartPipeline=TRUE;
				}else{
					printf("Big change ignored no pipeline restart\n");fflush(0);
				}
			}
		}		
	}

	if( bRestartPipeline ){
		if( m_FrameCounter>1 && m_FrameUnixTime>0 ){
			printf("Due to big change of coordinates restart of pipeline is required (frame=%d, day_frame=%d)\n",m_FrameCounter,m_DayFrameCounter);
			m_bRestartOnNext = TRUE;
		}else{
			printf("Frame = %d (%d) , time = %d, restart not needed\n",m_FrameCounter,m_DayFrameCounter,m_FrameUnixTime);
		}
	}


	return bRestartPipeline;			
}

BOOL_T CCDPipeline::SetCoordFromFITS()
{
	CCDMatrix& currFrame = (GetCurrent())[0];
   CSafeKeyTab& keyTab = currFrame.GetKeyTab();
	int tmp_int;
	double tmp_dbl;
	double diff_to_do_astro=0.2; // in deg 
	double limit_to_restart = gCCDParams.m_CoordChangeToRestart;
	BOOL_T bRestartPipeline=FALSE;

	printf("Retriving coordinates from FITS header ...\n");

	if(!keyTab.Get( OBSMODE, tmp_int ))
		return FALSE;
	eObservationMode_T obsMode = (eObservationMode_T)tmp_int;

	if( obsMode != GetObsMode() ){
		printf("Change of obs mode ( %d -> %d ) - forcing astrometry\n",GetObsMode(),obsMode);
		m_bForceAstroNow = TRUE;
		m_nAstrometryRetry = gCCDParams.m_nASASAstrometryReTry;
		if( obsMode == eEarthMovingMode ){
			bRestartPipeline=TRUE;
		}
	}


	if(!keyTab.Get( RA_OBS, tmp_dbl ))
		return FALSE;
	if( obsMode == eEarthMovingMode ){
		// change of ra only in tracking mode :
		if( fabs( AstroAngle::hours2deg( m_pAsasTransform->ra ) - AstroAngle::hours2deg( tmp_dbl ) ) > diff_to_do_astro ){
			if( m_PipelineCfg.m_IgnoreCoordChangeOnCamera != m_PipelineIndex ){			
				printf("Big change in RA - forcing astrometry (%.2f -> %.2f)\n",m_pAsasTransform->ra,tmp_dbl);fflush(0);
				m_bForceAstroNow = TRUE;
				m_nAstrometryRetry = gCCDParams.m_nASASAstrometryReTry;
				if( fabs( AstroAngle::hours2deg( m_pAsasTransform->ra ) - AstroAngle::hours2deg( tmp_dbl ) ) > limit_to_restart ){
					bRestartPipeline=TRUE;
				}
			}else{
				printf("Big change ignored no pipeline restart , no astrometry forced \n");fflush(0);
			}
		}
	}
	m_pAsasTransform->ra = tmp_dbl;
		

	if(!keyTab.Get( DEC_OBS, tmp_dbl ))
		return FALSE;
	if( obsMode == eEarthMovingMode ){
		// change of dec only in tracking mode :
		if( fabs( m_pAsasTransform->dec - tmp_dbl ) > diff_to_do_astro ){
			if( m_PipelineCfg.m_IgnoreCoordChangeOnCamera != m_PipelineIndex ){
				printf("Big change in DEC - forcing astrometry (%.2f -> %.2f)\n",m_pAsasTransform->dec,tmp_dbl);
				m_bForceAstroNow = TRUE;
				m_nAstrometryRetry = gCCDParams.m_nASASAstrometryReTry;
				if( fabs( m_pAsasTransform->dec - tmp_dbl ) > limit_to_restart ){
					bRestartPipeline=TRUE;
				}
			}else{
				printf("Big change ignored no pipeline restart , no astrometry forced \n");fflush(0);
			}
		}
	}
	m_pAsasTransform->dec = tmp_dbl;

	double azim_deg,alt_deg;
	if(!keyTab.Get( AZIM_OBS, azim_deg ))
		return FALSE;
	if(!keyTab.Get( ALT_OBS, alt_deg ))
		return FALSE;


	if( obsMode == eNoMovingMode ){
		// change of azim,alt only in tracking mode :

		if( fabs( AstroAngle::rad2deg( m_pAsasTransform->m_AzimInRad ) - azim_deg ) > diff_to_do_astro ){
			if( m_PipelineCfg.m_IgnoreCoordChangeOnCamera != m_PipelineIndex ){
				printf("Big change in AZIM - forcing astrometry (%.2f -> %.2f)\n",AstroAngle::rad2deg( m_pAsasTransform->m_AzimInRad ),azim_deg);
   		   m_bForceAstroNow = TRUE;
				m_nAstrometryRetry = gCCDParams.m_nASASAstrometryReTry;
				if( fabs( AstroAngle::rad2deg( m_pAsasTransform->m_AzimInRad ) - azim_deg ) > limit_to_restart ){
					bRestartPipeline=TRUE;
				}
			}else{
				printf("Big change ignored no pipeline restart , no astrometry forced \n");fflush(0);
			}
		}
		if( fabs( AstroAngle::rad2deg( m_pAsasTransform->m_AltInRad ) - alt_deg ) > diff_to_do_astro ){
			if( m_PipelineCfg.m_IgnoreCoordChangeOnCamera != m_PipelineIndex ){
				printf("Big change in ALT - forcing astrometry (%.2f -> %.2f)\n",AstroAngle::rad2deg( m_pAsasTransform->m_AltInRad ),alt_deg);
  		   	m_bForceAstroNow = TRUE;
				m_nAstrometryRetry = gCCDParams.m_nASASAstrometryReTry;
				if( fabs( AstroAngle::rad2deg( m_pAsasTransform->m_AltInRad ) - alt_deg ) > limit_to_restart ){
					bRestartPipeline=TRUE;
				}	
			}else{
				printf("Big change ignored no pipeline restart , no astrometry forced \n");fflush(0);
			}
		}
	}

	m_pAsasTransform->m_AzimInRad = AstroAngle::deg2rad( azim_deg );
	m_pAsasTransform->m_AltInRad = AstroAngle::deg2rad( alt_deg );
	

	if( bRestartPipeline ){
		if( m_FrameCounter>1 && m_FrameUnixTime>0 ){
			printf("Due to big change of coordinates restart of pipeline is required (frame=%d, day_frame=%d)\n",m_FrameCounter,m_DayFrameCounter);
			m_bRestartOnNext = TRUE;
		}else{
			printf("Frame = %d (%d) , time = %d, restart not needed\n",m_FrameCounter,m_DayFrameCounter,m_FrameUnixTime);
		}
	}

	UpdateObsCoordinates( AstroAngle::hours2rad( m_pAsasTransform->ra ),
								 AstroAngle::deg2rad( m_pAsasTransform->dec ),
								 m_pAsasTransform->m_AzimInRad, m_pAsasTransform->m_AltInRad,
								 obsMode );
	return TRUE;	
}

BOOL_T CCDPipeline::ReadASASAstrometryFromFITS( CCDMatrix& currFrame, CCDAsasTransform* pAsasTransform,
																eObservationMode_T& obsMode, int* nStars )
{
	CSafeKeyTab& keyTab = currFrame.GetKeyTab();
	int tmp_int;	
	if(!keyTab.Get( POSANGLE, pAsasTransform->fi )){
		printf("%s keyword missing !\n",POSANGLE);
		return FALSE;
	}
   if(!keyTab.Get( AST_ORD,  pAsasTransform->order )){
		printf("%s keyword is missing !\n",AST_ORD);
		return FALSE;
	}
   if(!keyTab.Get( PIXSCALE, pAsasTransform->pixscale )){
		printf("%s keyword is missing !\n",PIXSCALE);
		return FALSE;
	}
	if(!keyTab.Get( RA_OBS, pAsasTransform->ra )){
		printf("%s keyword is missing !\n",RA_OBS);
		return FALSE;
	}
	if(!keyTab.Get( DEC_OBS, pAsasTransform->dec )){
		printf("%s keyword is missing !\n",DEC_OBS);
		return FALSE;
	}

	double azim_deg,alt_deg;
	if(!keyTab.Get( AZIM_OBS, azim_deg )){
		printf("%s keyword is missing !\n",AZIM_OBS);
		return FALSE;
	}
	if(!keyTab.Get( ALT_OBS, alt_deg )){
		printf("%s keyword is missing !\n",ALT_OBS);
		return FALSE;
	}
	pAsasTransform->m_AzimInRad = AstroAngle::deg2rad( azim_deg );
	pAsasTransform->m_AltInRad = AstroAngle::deg2rad( alt_deg );
	
	obsMode = eEarthMovingMode;
	if(!keyTab.Get( OBSMODE, tmp_int )){
		printf("WARNING : %s keyword is missing !\n",OBSMODE);
// changed on 20070811 by MS , it is not so important, in fact most of 
// image are taken in tracking mode , so constant mount mode is not used 
// almost at all
//		return FALSE;
	}else{
		obsMode = (eObservationMode_T)tmp_int;
	}

   if(!keyTab.Get( AST_UTTIME, tmp_int )){
		printf("%s keyword is missing !\n",AST_UTTIME);
		return FALSE;
	}
	pAsasTransform->m_TransformUtTime = tmp_int;

	for(int i=0;i<14;i++){
  		mystring szParX;
      szParX << PARX << i;
		if(!keyTab.Get( szParX.c_str(), pAsasTransform->px[i] )){
			printf("%s keyword is missing !\n",szParX.c_str());
			return FALSE;
		}
   }
   for(int i=0;i<14;i++){
      mystring szParY;
      szParY << PARY << i;
      if(!keyTab.Get( szParY.c_str(), pAsasTransform->py[i] )){
			printf("%s keyword is missing !\n",szParY.c_str());
			return FALSE;
		}
   }

	// [ NEW ]
	const char* szSizeX = currFrame.GetKeyTab().getKeyVal( ORGSIZEX );
   const char* szSizeY = currFrame.GetKeyTab().getKeyVal( ORGSIZEY );
	if( szSizeX && szSizeY && strlen( szSizeX ) && strlen( szSizeY ) ){
      if( pAsasTransform ){
         pAsasTransform->m_SizeX = atol( szSizeX );
         pAsasTransform->m_SizeY = atol( szSizeY );
         pAsasTransform->xc = (pAsasTransform->m_SizeX/2);
         pAsasTransform->yc = (pAsasTransform->m_SizeY/2);
      }
   }	

	if( !pAsasTransform->m_SizeX || !pAsasTransform->m_SizeY ){
		pAsasTransform->m_SizeX = currFrame.GetXSize();
		pAsasTransform->m_SizeY = currFrame.GetYSize();
		pAsasTransform->xc = (pAsasTransform->m_SizeX/2);
      pAsasTransform->yc = (pAsasTransform->m_SizeY/2);
	}

	pAsasTransform->m_bTransformOK = TRUE;

	if( nStars ){
		keyTab.Get( NSTARS, (*nStars) );
		printf("Number of stars from fits = %d\n",(*nStars));
	}

	return TRUE;
}

BOOL_T CCDPipeline::ReadASASAstrometryFromFITS( CSafeKeyTab& keyTab, 
																CCDAsasTransform* pAsasTransform,
																eObservationMode_T& obsMode )
{
	int tmp_int;	
	if(!keyTab.Get( POSANGLE, pAsasTransform->fi )){
		printf("%s keyword missing !\n",POSANGLE);
		return FALSE;
	}
   if(!keyTab.Get( AST_ORD,  pAsasTransform->order )){
		printf("%s keyword is missing !\n",AST_ORD);
		return FALSE;
	}
   if(!keyTab.Get( PIXSCALE, pAsasTransform->pixscale )){
		printf("%s keyword is missing !\n",PIXSCALE);
		return FALSE;
	}
	if(!keyTab.Get( RA_OBS, pAsasTransform->ra )){
		printf("%s keyword is missing !\n",RA_OBS);
		return FALSE;
	}
	if(!keyTab.Get( DEC_OBS, pAsasTransform->dec )){
		printf("%s keyword is missing !\n",DEC_OBS);
		return FALSE;
	}

	double azim_deg,alt_deg;
	if(!keyTab.Get( AZIM_OBS, azim_deg )){
		printf("%s keyword is missing !\n",AZIM_OBS);
		return FALSE;
	}
	if(!keyTab.Get( ALT_OBS, alt_deg )){
		printf("%s keyword is missing !\n",ALT_OBS);
		return FALSE;
	}
	pAsasTransform->m_AzimInRad = AstroAngle::deg2rad( azim_deg );
	pAsasTransform->m_AltInRad = AstroAngle::deg2rad( alt_deg );
	
	obsMode = eEarthMovingMode;
	if(!keyTab.Get( OBSMODE, tmp_int )){
		printf("WARNING : %s keyword is missing !\n",OBSMODE);
// changed on 20070811 by MS , it is not so important, in fact most of
// image are taken in tracking mode , so constant mount mode is not used
// almost at all
//		return FALSE;
	}else{
		obsMode = (eObservationMode_T)tmp_int;
	}

   if(!keyTab.Get( AST_UTTIME, tmp_int )){
		printf("%s keyword is missing !\n",AST_UTTIME);
		return FALSE;
	}
	pAsasTransform->m_TransformUtTime = tmp_int;

	for(int i=0;i<14;i++){
  		mystring szParX;
      szParX << PARX << i;
		if(!keyTab.Get( szParX.c_str(), pAsasTransform->px[i] )){
			printf("%s keyword is missing !\n",szParX.c_str());
			return FALSE;
		}
   }
   for(int i=0;i<14;i++){
      mystring szParY;
      szParY << PARY << i;
      if(!keyTab.Get( szParY.c_str(), pAsasTransform->py[i] )){
			printf("%s keyword is missing !\n",szParY.c_str());
			return FALSE;
		}
   }

	// [ NEW ]
	const char* szSizeX = keyTab.getKeyVal( ORGSIZEX );
   const char* szSizeY = keyTab.getKeyVal( ORGSIZEY );
	if( szSizeX && szSizeY && strlen( szSizeX ) && strlen( szSizeY ) ){
      if( pAsasTransform ){
         pAsasTransform->m_SizeX = atol( szSizeX );
         pAsasTransform->m_SizeY = atol( szSizeY );
         pAsasTransform->xc = (pAsasTransform->m_SizeX/2);
         pAsasTransform->yc = (pAsasTransform->m_SizeY/2);
      }
   }	

	if( !pAsasTransform->m_SizeX || !pAsasTransform->m_SizeY ){
		szSizeX = keyTab.getKeyVal( FH_NAXIS1 );
		szSizeY = keyTab.getKeyVal( FH_NAXIS2 );

		if(  szSizeX && szSizeY && strlen( szSizeX ) && strlen( szSizeY ) ){
			pAsasTransform->m_SizeX = atol( szSizeX );
			pAsasTransform->m_SizeY = atol( szSizeY );
			pAsasTransform->xc = (pAsasTransform->m_SizeX/2);
	      pAsasTransform->yc = (pAsasTransform->m_SizeY/2);
		}
	}

	pAsasTransform->m_bTransformOK = TRUE;

	return TRUE;
}



BOOL_T CCDPipeline::GetASASAstrometryFromFITS()
{
	CCDMatrix& currFrame = (GetCurrent())[0];
	eObservationMode_T obsMode;

	if(!ReadASASAstrometryFromFITS( currFrame, m_pAsasTransform, obsMode, &m_nPhotometryStarsCount )){
		printf("Could not retrive astrometry information from fits file : %s\n",currFrame.m_szFileName.c_str());
		return FALSE;
	}

	UpdateObsCoordinates( AstroAngle::hours2rad( m_pAsasTransform->ra ),
								 AstroAngle::deg2rad( m_pAsasTransform->dec ),
								 m_pAsasTransform->m_AzimInRad, m_pAsasTransform->m_AltInRad,
								 obsMode );

	return TRUE;
}


void CCDPipeline::RemoveOldRejectIfMoreTracks()
{
	BOOL_T bDoDel=FALSE;
	CTrackList::iterator i;
	CTrackList::iterator pDelEnd;
	for(i=m_CurrFrameTracks.begin();i!=m_CurrFrameTracks.end();i++){
		if( (m_FrameCounter - i->frame_index) > gCCDParams.m_nKeepRejectIfMoreTracksOnN ){
			pDelEnd = i;
			bDoDel=TRUE;
		}		
	}
	if( bDoDel ){
		m_CurrFrameTracks.erase( m_CurrFrameTracks.begin(), pDelEnd+1 );
	}
}

BOOL_T CCDPipeline::CalcFrameShift( int nStarCount, double mag_min, double mag_max,
												int min_x, int min_y,int max_x,int max_y,
											   int min_dist, double delta_mag, BOOL_T bVerb )
{
	if( m_Count >= 2 ){
		cCCD& pCurr = GetCurrent();
		cCCD* pPrev = GetPrevFrame();

		if( pPrev ){
			CCDMatrix& image1 = (*pPrev)[0];
			CCDMatrix& image2 = pCurr[0];

			Area2DInfo info;

			double lapm2,laps2;
			image2.GetVariableMeanAndSigma( gCCDParams.m_eLaplaceType,
									lapm2, laps2, info, 20, 20,
									(image2.GetXSize()-20), (image2.GetYSize()-20),
									0, image2.get_laplace_data_fast() );
			printf("Laplace=%d mean=%.2f sigma=%.2f\n",gCCDParams.m_eLaplaceType,lapm2, laps2); 				

			time_t start=get_dttm();
			CPixelAnalyseOut out;
			CPixelAnalyseIn in;
			CPixelList pixel_list(image1.GetXSize()*image1.GetYSize());
			CPixelList alreadyUsed(image1.GetXSize()*image1.GetYSize());
			in.pPixelList = &alreadyUsed;

			double n_sigma_min=10,n_sigma_max=15;

			m_PrevStarList = m_CurrStarList;
			CPiPhotometry::GetFastPhotoList( image2, m_CurrStarList, out, in, pixel_list,
					  							   n_sigma_min,n_sigma_max,5,laps2, 20 );
			if( m_PrevStarList.size()==0 ){
				pixel_list.Clear();
				alreadyUsed.Clear();

				double lapm1,laps1;
				image1.GetVariableMeanAndSigma( gCCDParams.m_eLaplaceType, 
									lapm1, laps1, info, 20, 20,
									(image1.GetXSize()-20), (image1.GetYSize()-20),
									0, image1.get_laplace_data_fast() );
				printf("Laplace=%d mean=%.2f sigma=%.2f\n",gCCDParams.m_eLaplaceType,lapm1,laps1);				

				CPiPhotometry::GetFastPhotoList( image1, m_PrevStarList, out, in, pixel_list, 
															n_sigma_min,n_sigma_max,5,laps1, 20 );
			}

/*			vector<cStarCat> star_list1, star_list2;
			CPiPhotometry::GetFastPhotoList( image1, star_list1, out, in, pixel_list,
				  							   n_sigma_min,n_sigma_max,5,laps1, 20 );	
			pixel_list.Clear();
			alreadyUsed.Clear();
			CPiPhotometry::GetFastPhotoList( image2, star_list2, out, in, pixel_list,
					  							   n_sigma_min,n_sigma_max,5,laps2, 20 );*/

			printf("Number of found stars on previous image = %d\n",m_PrevStarList.size());
			printf("Number of found stars on current image = %d\n",m_CurrStarList.size());
							

			double dx,dy,sigma_dx,sigma_dy;
			int match_count = CPiPhotometry::CalcImageShift( m_PrevStarList, m_CurrStarList, 
											 dx, dy, sigma_dx, sigma_dy,
											 nStarCount, mag_min, mag_max, 
											 min_x, max_x, min_y, max_y, min_dist, 
											 delta_mag, bVerb, 2 );
			time_t end=get_dttm();
			printf("Number of matched = %d\n",match_count);
			printf("Average shift is : (%.2f,%.2f)\n",dx,dy);
			printf("RMS of shift     : (%.2f,%.2f)\n",sigma_dx,sigma_dy);
			printf("Time             : %d sec\n",(end-start));			

			image2.m_dx = 0;
			image2.m_dy = 0;
			if( sigma_dx < 2.00 && sigma_dy < 2.00 && match_count>1000 ){	
				printf("Shift RMS accepted\n");
				if( fabs(dx) >= gCCDParams.m_MinShiftToUse ){
					image2.m_dx = dx;
				}else{
					printf("shift dx=%.2f to small - ignored\n",dx);
				}
				if( fabs(dy) >= gCCDParams.m_MinShiftToUse ){
					image2.m_dy = dy;
				}else{
					printf("shift dy=%.2f to small - ignored\n",dy);
				}
			}else{
				printf("Shift RMS = (%.2f,%.2f) to high or matched_count=%d to low , probably mount move to big\n",sigma_dx,sigma_dy,match_count);
				// set flag to restart pipeline 				
			}

			if( gCCDParams.m_ForceAstroWhenBigShiftRMS > 0 ){
				// condition for big shift to make astrometry is >MaxShift 
				// or RMS of shift distribution >Limit - which means that 
				// probably many stars were misidentified on new image
				if( (fabs(dx)>gCCDParams.m_ForceAstroWhenBigShift && sigma_dx<gCCDParams.m_ForceAstroWhenBigShiftRMS ) ||
					 (fabs(dy)>gCCDParams.m_ForceAstroWhenBigShift && sigma_dy<gCCDParams.m_ForceAstroWhenBigShiftRMS ) ||
					 sigma_dx>gCCDParams.m_ForceAstroWhenBigShiftRMS || sigma_dy>gCCDParams.m_ForceAstroWhenBigShiftRMS ){
					printf("Large shift detected (dx,dy)=(%.2f,%.2f) , (sigma_dx,sigma_dy)=(%.2f,%.2f) -> ASTROMETRY FORCED\n",dx,dy,sigma_dx,sigma_dy);fflush(stdout);
					m_bForceAstroNow = TRUE;
					m_bForceSynchroMode = TRUE;
				}
			}
		}
	}

	return FALSE;
}

void CCDPipeline::UpdateAutoShiftStars()
{
	if( gCCDParams.m_bReadFirstLevelInfoFromLog && !gCCDParams.m_bFromLogFile_WithFrames ){
		printf("CCDPipeline::UpdateAutoShiftStars re-analysing of log files without frames, calculation of shifts ignored\n");
		return;
	}

	if( gCCDParams.m_bUseFrameShift ){
		CalcFrameShift();
		return;
	}

	
	eObservationMode_T curr_mode = GetCamObsMode();
	if( curr_mode==eEarthMovingMode ){
		// turn off shifts :
		_TRACE_PRINTF_1("Pipeline%d EARTH ROTATION MODE - turning off shifts !!!\n",m_PipelineIndex);fflush(0);
		m_PipelineCfg.SetShiftsToZero();
		if( m_PipelineCfg.m_nAutoShiftsCalc>0 ){
			m_PipelineCfg.m_nAutoShiftsCalc = -abs( m_PipelineCfg.m_nAutoShiftsCalc );
		}
		if( curr_mode != m_PrevFrameObsMode || m_FrameCounter<=1 ){
			PrintError( this, GetDayFrameCounter(), "WRN_SHFITS_OFF_OBSMODE1",  "Shifts disabled due to tracking mode" );
		}
	}
	m_PrevFrameObsMode = curr_mode;

	if( m_PipelineCfg.m_nAutoShiftsCalc>0 && m_StarSpyTab[0].GetFramesCount()<=m_PipelineCfg.m_nAutoShiftsCalc ){
	//if( m_PipelineCfg.m_nAutoShiftsCalc>0 && (m_FrameCounter-m_StartFrameIndex)<=m_PipelineCfg.m_nAutoShiftsCalc ){
		if( m_StarSpyTab[0].GetFramesCount()==0 ){
			// if((m_FrameCounter-m_StartFrameIndex)==1){
			InitMAXStars( m_StarSpyTab );		
		}else{
			BOOL_T bUpdateOK=TRUE;
			if( m_StarSpyTab[0].GetFramesCount()<=m_PipelineCfg.m_nAutoShiftsCalc ){
				// if ((m_FrameCounter-m_StartFrameIndex)<=m_PipelineCfg.m_nAutoShiftsCalc){
				bUpdateOK = UpdateMAXStarsPositions( m_StarSpyTab );
			}
			if( m_StarSpyTab[0].GetFramesCount()==m_PipelineCfg.m_nAutoShiftsCalc || !bUpdateOK){
				// if( (m_FrameCounter-m_StartFrameIndex)==m_PipelineCfg.m_nAutoShiftsCalc || !bUpdateOK){
				if(!CalcRotations()){
					printf("\n\nERROR : could not determine rotations values, exiting ...!!!\n\n\n");
					MYTRACE0(gCCDTrace,"could not determine rotations values, exiting ...!!!\n");
					/*if(!gCCDParams.m_bShiftUsesAstroFormulas){
						// exiting only in NON-ASTROFORMULAS MODE !
						exit(0);
					}*/
					
					// no exit now - retry again :
					InitMAXStars( m_StarSpyTab );	
				}else{
					// to ignore this parameter in future :
					m_PipelineCfg.m_nAutoShiftsCalc=-abs(m_PipelineCfg.m_nAutoShiftsCalc);
					gCCDParams.SetParam("CCD_AUTO_SHIFTS_CALC","0");
				}
			}else{
				if( m_StarSpyTab[0].GetFramesCount()<m_PipelineCfg.m_nAutoShiftsCalc ){
				// if( (m_FrameCounter-m_StartFrameIndex)<m_PipelineCfg.m_nAutoShiftsCalc){
					UpdateRotations();
				}
			}
		}
	}
}

void CCDPipeline::UpdateTraceOnAll()
{
	if(gCCDParams.m_bTraceOnAllFrames){
		double currTime = (GetCurrent()[0]).getObsTime( gCCDParams.m_bIgnoreMissingTime );
		if((m_FrameCounter-m_StartFrameIndex)==1){
			m_CurrTime = currTime - fabs( gCCDParams.m_DriverShutterTimeInSec );
		}
		double diff = (currTime - m_CurrTime);
		m_CurrTime = currTime;
		if( diff>200*fabs(gCCDParams.m_DriverShutterTimeInSec) && (m_FrameCounter-m_StartFrameIndex)>5){
			// require re-initialization of all traced stars :
			printf("Very long delay between frames : %2f\n",diff);
			printf("Please VERIFY if shutter time agrees with real-shutter TIME !!!!!!!!!\n");
			printf("Re-initializing traced stars\n");
			MarkTracedAllToReInit();			
			gCCDParams.m_nSkipNFramesAfterChange = gCCDParams.m_nMaxOfAverageOfPrevN;
		}

		if((m_FrameCounter-m_StartFrameIndex)==1){
			if( m_StarSpyTab.size()>0 && !m_StarSpyTab.m_bFromFile){
				m_StarsForMatrix = m_StarSpyTab;
			}else{
				InitMAXStars( m_StarsForMatrix );
			}
		}else{
			if(!UpdateMAXStarsPositions( m_StarsForMatrix )){
				// special actions - one of stars must re-initialized :
			}
			CheckIfReInitNeeded( m_StarsForMatrix );
		}
		m_StarsForMatrix.DumpReportToFile();
		if(gCCDParams.m_bUseTransformMatrix){
			UpdateTransformMatrix( m_StarsForMatrix );
		}
	}

}

void CCDPipeline::UpdateControlStar()
{
	if(GetCount()==m_nPipelineSize && gCCDParams.m_bUseControllStar){
		CCDMatrix& current = GetCurrent()[0];
		double frame_time = current.getObsTime( TRUE );
		// update controll star :
		if(m_FrameCounter==m_nPipelineSize){
			trace("Initializing controll star\n");	
			m_pControlSpyStar->InitWithMAXStar( current.get_data_buffer_fast(), frame_time , m_FrameCounter );			
		}else{
			m_pControlSpyStar->ReCalcNewPosition( GetCurrent()[0].get_data_buffer_fast(),
														     (int)frame_time, m_FrameCounter );
		}
		m_pControlSpyStar->CalcPredicted( this, frame_time, 0 );
		m_pControlSpyStar->LogCurrentPosition( "controlStar", 0, this );		
		if(! (m_pControlSpyStar->m_StarDesc)[0].m_bNotEnd ){
			trace("Re-Initializing controll star\n");
			m_pControlSpyStar->InitWithMAXStar( current.get_data_buffer_fast(), frame_time , m_FrameCounter );
		}
	}
}

void CCDPipeline::PrintShiftParams()
{
	printf("dAlfa=%.4f , center = (%.4f,%.4f), dx=%.4f, dy=%.4f,"\
			 "dx_per_sec=%.4f, dx_per_sec=%.4f, rot=%d\n",m_StarSpyTab[0].m_dAlfa,
			 m_StarSpyTab[0].m_RotCenterX,m_StarSpyTab[0].m_RotCenterY,
			 m_StarSpyTab[0].m_DX, m_StarSpyTab[0].m_DY,
          m_StarSpyTab[0].m_DXPerSec, m_StarSpyTab[0].m_DYPerSec,
          m_StarSpyTab[0].m_bRotation );
			
}

void CCDPipeline::SetShiftParams()
{
	for(int i=0;i<m_nCCD;i++){
		SetShiftParams( m_StarSpyTab[i].m_dAlfa	, m_StarSpyTab[i].m_RotCenterX,
							 m_StarSpyTab[i].m_RotCenterY, 
							 m_StarSpyTab[i].m_DX, m_StarSpyTab[i].m_DY, 
							 m_StarSpyTab[i].m_DXPerSec, m_StarSpyTab[i].m_DYPerSec,
      					 m_StarSpyTab[i].m_bRotation, (m_nCCD==1), i, m_StarSpyTab[i].m_dAlfaPerSec );
	}
}

void CCDPipeline::SetShiftParams( double angle, double x, double y, double dx, double dy, 
											 double dxpersec, double dypersec, BOOL_T bRot,
											 BOOL_T bSetGlobal, int camIdx, double angle_per_sec )
{
	if( bSetGlobal ){
		if(bRot && !gCCDParams.m_bDoNotUseRotation){
			gCCDParams.m_bUseRotInAverageOfPrevN=TRUE;
			gCCDParams.SetParam("CCD_USE_ROT_IN_AVER_OF_PREV","1");
		}else{
			gCCDParams.m_bUseRotInAverageOfPrevN=FALSE;
      	gCCDParams.SetParam("CCD_USE_ROT_IN_AVER_OF_PREV","0");
		}			

		gCCDParams.m_RotValueDAlfa = angle;
		gCCDParams.m_RotValueDAlfaPerSec = angle_per_sec;
		gCCDParams.m_RotCenterX = x;
		gCCDParams.m_RotCenterY = y;
		gCCDParams.m_FrameDX = dx;
		gCCDParams.m_FrameDY = dy;
		gCCDParams.m_FrameDXPerSec = dxpersec;
      gCCDParams.m_FrameDYPerSec = dypersec;
			
		gCCDParams.SetParam("CCD_SINGLE_FRAME_DX",dx);
		gCCDParams.SetParam("CCD_SINGLE_FRAME_DY",dy);
		gCCDParams.SetParam("CCD_ROT_CENTER_X",x);
		gCCDParams.SetParam("CCD_ROT_CENTER_Y",y);
		gCCDParams.SetParam("CCD_SINGLE_FRAME_D_ALFA",angle);
		gCCDParams.SetParam("CCD_SINGLE_FRAME_D_ALFA_PER_SEC",angle_per_sec);
		gCCDParams.SetParam("CCD_SINGLE_FRAME_DX_PER_SEC",dxpersec);
		gCCDParams.SetParam("CCD_SINGLE_FRAME_DY_PER_SEC",dypersec);
	}

	if(camIdx<m_CamCfgTab.size()){	
		CCcdCfg* pCfg = (CCcdCfg*)(m_CamCfgTab[camIdx]);

		pCfg->m_CCDParams.m_RotCenterX = x;
		pCfg->m_CCDParams.m_RotCenterY = y;
		if(bRot && !gCCDParams.m_bDoNotUseRotation){
			pCfg->m_CCDParams.m_bUseRotInAverageOfPrevN = TRUE;
		}else{
			pCfg->m_CCDParams.m_bUseRotInAverageOfPrevN = FALSE;
		}
		pCfg->m_CCDParams.m_RotValueDAlfa = angle;
		pCfg->m_CCDParams.m_RotValueDAlfaPerSec = angle_per_sec;
		pCfg->m_CCDParams.m_FrameDX = dx;
		pCfg->m_CCDParams.m_FrameDY = dy;
		pCfg->m_CCDParams.m_FrameDXPerSec = dxpersec;
		pCfg->m_CCDParams.m_FrameDYPerSec = dypersec;
	}

	m_PipelineCfg.m_RotCenterX = x;
	m_PipelineCfg.m_RotCenterY = y;
 	if( gCCDParams.m_bDoNotUseRotation )
		m_PipelineCfg.m_bUseRotInAverageOfPrevN = FALSE;
	else
		m_PipelineCfg.m_bUseRotInAverageOfPrevN = bRot;
	m_PipelineCfg.m_FrameDX = dx;
	m_PipelineCfg.m_FrameDY = dy;
	m_PipelineCfg.m_FrameDXPerSec = dxpersec;
	m_PipelineCfg.m_FrameDYPerSec = dypersec;
	m_PipelineCfg.m_RotValueDAlfa = angle;
	m_PipelineCfg.m_RotValueDAlfaPerSec = angle_per_sec;


	CCD_Analyser::AutoCalculateIgnoreEdges( this );
	
	LogShiftChange();
}

void CCDPipeline::CheckIfShiftNeeded()
{
	eObservationMode_T obsMode = GetObsMode();
	if( obsMode==eEarthMovingMode ){
		// shifts not needed in this case :
		if( gCCDParams.m_nAutoShiftsCalc>0 ){
			gCCDParams.m_nAutoShiftsCalc = -abs( gCCDParams.m_nAutoShiftsCalc );
		}
		vector<CCDPipeline*>::iterator i;
      for(i=m_PipelineList.begin();i!=m_PipelineList.end();i++){
			if( ((*i)->m_PipelineCfg).m_nAutoShiftsCalc>0 ){
				 ((*i)->m_PipelineCfg).m_nAutoShiftsCalc = -abs( ((*i)->m_PipelineCfg).m_nAutoShiftsCalc );				
			}
			(*i)->SetShiftParams( 0 , 0 , 0 , 0 , 0 , 0 ,0 , 0 , TRUE, -1 , 0 );
		}
	}
}

void CCDPipeline::LogShiftChange()
{
	mystring szOut;
	GetOutFileName( szOut, TRANSFORM_SUBDIR, "shiftChange", this );
	if(!MyFile::DoesFileExist( szOut.c_str() )){
		MyOFile out( szOut.c_str() );
		out.Printf(SHIFT_LOG_FMT_DESC);
	}
	MyOFile out( szOut.c_str() , "a" );
	out.Printf(SHIFT_LOG_FMT,GetFrameIndex(),m_PipelineCfg.m_bUseRotInAverageOfPrevN,
			m_PipelineCfg.m_RotCenterX,m_PipelineCfg.m_RotCenterY,
			m_PipelineCfg.m_FrameDX,m_PipelineCfg.m_FrameDY,
			m_PipelineCfg.m_FrameDXPerSec,m_PipelineCfg.m_FrameDYPerSec,
			m_PipelineCfg.m_RotValueDAlfa,m_PipelineCfg.m_RotValueDAlfaPerSec );
}


BOOL_T CCDPipeline::UpdateRotations()
{
	BOOL_T bRet=TRUE;
	cCCD& currentFrame = GetCurrent();
   for(int i=0;i<m_nCCD;i++){
      double angle,x,y,dx,dy,dxpersec,dypersec,angle_per_sec;
		BOOL_T bRotation=FALSE;
		double endTime = currentFrame[i].getObsTime();

		if(!m_StarSpyTab[i].CalcRotation( angle,x,y,dx,dy, dxpersec,dypersec, bRotation,
                                        angle_per_sec,endTime, 1 )){
         return FALSE;
      }else{
         SetShiftParams( angle, x, y, dx, dy, dxpersec,dypersec, bRotation, (m_nCCD==1), i, angle_per_sec );
			mystring szTrace;
			szTrace << "For Pipeline=" << m_PipelineIndex << " ccd=" << i 
					  << " => updated rotations to dAlfa=" << angle << " dAlfaPerSec=" << angle_per_sec
					  << " dx=" << dx << " dy=" << dy << " x=" << x << " y=" << y << " dxpersec=" 
					  << dxpersec << " dypersec=" << dypersec << " bRot=" << bRotation;
			_TRACE_PRINTF_0("%s\n",szTrace.c_str());
      }
	}
}

void CCDPipeline::MarkTracedAllToReInit()
{
	cCCD& currentFrame = GetCurrent();
	m_StarSpyTab.MarkToReInit();
	m_StarsForMatrix.MarkToReInit();		
}

BOOL_T CCDPipeline::CalcRotations()
{
	cCCD& currentFrame = GetCurrent();
   for(int i=0;i<m_nCCD;i++){
		double angle,x,y,dx,dy,dxpersec,dypersec,angle_per_sec;
		BOOL_T bRotation=FALSE;
		double endTime = currentFrame[i].getObsTime( gCCDParams.m_bIgnoreMissingTime );
		
		if(!m_StarSpyTab[i].CalcRotation( angle,x,y,dx,dy, dxpersec,dypersec, bRotation, 
													 angle_per_sec,endTime )){
			return FALSE;
		}else{
			SetShiftParams( angle, x, y, dx, dy, dxpersec,dypersec, bRotation, (m_nCCD==1), i, angle_per_sec );
		}								
	}


	mystring szShiftFileName;
   GetShiftFileName( szShiftFileName );	

	mystring szDateFileName;
	szDateFileName << get_mydate() << "_" << m_PipelineIndex << ".cfg";

	m_StarSpyTab.SaveToFile( szShiftFileName.c_str(), szDateFileName.c_str(), this );

	return TRUE;
}

BOOL_T CCDPipeline::UpdateTransformMatrix( CCDStarSpyTab& starsSpyTab )
{
	cCCD& currentFrame = GetCurrent();
	BOOL_T bRet = TRUE;
   for(int i=0;i<m_nCCD;i++){
		int startFrom=(m_FrameCounter-gCCDParams.m_UseNBackFrameForMatrix);
		if(startFrom<=0){
			startFrom = 1;
		}
		if(m_FrameCounter > startFrom){
			bRet = bRet && starsSpyTab[i].FindTransform( startFrom, m_FrameCounter, m_PipelineCfg, i, this );
			if(!bRet){
				printf("WARNING Pipeline:%d , could not determine transformation matrix\n",m_PipelineIndex);
			}
		}else{
			printf("INFO : Skiping determination of transform matrix on frame :%d\n",m_FrameCounter);
		}
	}
	return bRet;
}

void CCDPipeline::CheckIfReInitNeeded( CCDStarSpyTab& starsSpyTab )
{
	cCCD& currentFrame = GetCurrent();
   for(int i=0;i<m_nCCD;i++){
		ELEM_TYPE** p_data = currentFrame[i].get_data_buffer_fast();
		double frameTime = currentFrame[i].getObsTime( gCCDParams.m_bIgnoreMissingTime );

		starsSpyTab[i].ReInitMAXStarIfNeeded( p_data, (int)frameTime ,m_FrameCounter );
	}	
}

void CCDPipeline::LogMAXStarsPositions( CCDStarSpyTab& starsSpyTab )
{
	if( starsSpyTab.size()>0 ){
		starsSpyTab[0].LogMAXStarsPositions( this );
	}
}

BOOL_T CCDPipeline::UpdateMAXStarsPositions( CCDStarSpyTab& starsSpyTab )
{	
	cCCD& currentFrame = GetCurrent();
   for(int i=0;i<m_nCCD;i++){
		ELEM_TYPE** p_data = currentFrame[i].get_data_buffer_fast();
		double frameTime = currentFrame[i].getObsTime( gCCDParams.m_bIgnoreMissingTime );

		if(!starsSpyTab[i].ReCalcNewPosition( p_data, (int)frameTime ,m_FrameCounter )){
      	printf("Could not continue tracing star\n");
			MYTRACE0(gCCDTrace,"Could not continue tracing star frame = " << (int)m_FrameCounter);
			return FALSE;
		}
	}
	LogMAXStarsPositions( starsSpyTab );
	return TRUE;
}



void CCDPipeline::InitMAXStars( CCDStarSpyTab& starsSpyTab )
{
	cCCD& currentFrame = GetCurrent();
	for(int i=0;i<m_nCCD;i++){
		ELEM_TYPE** p_data = currentFrame[i].get_data_buffer_fast();
		double startTime = currentFrame[i].getObsTime( gCCDParams.m_bIgnoreMissingTime );


		starsSpyTab[i].InitWithMAXStar( p_data, startTime, m_FrameCounter );
	}
	LogMAXStarsPositions( starsSpyTab );
}

const char* CCDPipeline::CheckNextFrameFileName()
{
	if(m_FrameIdx<m_pList->GetCount()){
		return (m_pList->GetListTable())[m_FrameIdx].c_str();
	}
	return "";
}

int CCDPipeline::CheckNextDayFrameNo()
{
	int ret=0;
	if(m_FrameIdx<m_pList->GetCount()){
		mystring szNext = (m_pList->GetListTable())[m_FrameIdx];
		char cam;
		int dt;
		if( sscanf(szNext.c_str(),"k2%c_%d_%d.fitc",&cam,&dt,&ret)!=3 ){
			if( sscanf(szNext.c_str(),"k2%c_%d_%d.fit",&cam,&dt,&ret)!=3 ){
				ret=0;
			}
		}
	}
	return ret;
}

int CCDPipeline::ReadEvents_FromLogFile()
{
	vector<CCDPipeline*>::iterator i;
	int nMinDayFrame=100000;
   for(i=m_PipelineList.begin();i!=m_PipelineList.end();i++){
		
		if( ((*i)->m_EventsFromLogFile).size()==0 ){
		   mystring szAllEvents;
   		szAllEvents << gCCDParams.m_szFirstLevelTriggerDir.c_str() << "/"
     			         << "allevents_" << (*i)->m_PipelineIndex << ".log";
		   printf("reading log file %s\n",szAllEvents.c_str());fflush(0);
		
			((*i)->m_EventsFromLogFile).Read( szAllEvents.c_str() );
			if( ((*i)->m_EventsFromLogFile).size() <= 0 ){
				printf("No frames in log file : %s\n",szAllEvents.c_str());
				printf("cannot continue, exiting now\n");
				exit(0);
			}
		}
		int first_frame = (((*i)->m_EventsFromLogFile)[0]).m_DayFrameIndex;
		if( first_frame < nMinDayFrame ){
			nMinDayFrame = first_frame;
		}
	}	

	int startFrame = nMinDayFrame-gCCDParams.m_nPipelineSize+1;
	if( startFrame<=0 )
		startFrame = 1;
	
	for(i=m_PipelineList.begin();i!=m_PipelineList.end();i++){
		(*i)->m_LastFrameFromLogFile = startFrame;
	}
	return startFrame;
}

BOOL_T CCDPipeline::GetNextFrame_FromLogFile( mystring& szFileNames, BOOL_T bUseCache)
{
	if( m_PipelineIndex==0 ){
		// update new frame index only when called for pipeline 0 :

		if( m_LastFrameFromLogFile	< 0 ){
			printf("CCDPipeline::GetNextFrame_FromLogFile - initializing frame number pipeline %d\n",m_PipelineIndex);

			ReadEvents_FromLogFile();
		}else{
			vector<CCDPipeline*>::iterator i;
			int newFrameIndex = m_LastFrameFromLogFile+1;

			for(i=m_PipelineList.begin();i!=m_PipelineList.end();i++){
				(*i)->m_LastFrameFromLogFile = newFrameIndex;
			}
		}
	}

	int size = m_EventsFromLogFile.size();
	for(int i=0;i<size;i++){
		if( m_EventsFromLogFile[i].m_DayFrameIndex == m_LastFrameFromLogFile ){
			m_FrameUnixTime_FromLogFile = m_EventsFromLogFile[i].m_Time;
			break;
		}
	}
		
	m_DayFrameCounter = m_LastFrameFromLogFile;
	m_FrameCounter    = m_LastFrameFromLogFile;	

	printf("Read frame %d FromLogFile , at unix_time=%d\n",m_DayFrameCounter,m_FrameUnixTime_FromLogFile);

	int last_frame = m_EventsFromLogFile.back().m_DayFrameIndex;
	return (m_DayFrameCounter<=last_frame);
}

BOOL_T CCDPipeline::GetNextFrame_MC_AutoSum( mystring& szFileNames )
{
	BOOL_T bOK=FALSE;
	int sumed_count=0;
	mystring szPrevObject,szProcessedList,szOutDir,szOutList;
   double prev_ra=0,prev_dec=0;
   int first=0;
   int count=1;
   CMyStrTable aver_file_list;
	int n_aver = gCCDParams.m_bOnSumedFrames;
	int x_size = m_SizeX;
	int y_size = m_SizeY;
	time_t minTime=get_dttm()+10000;
   time_t maxTime=0;
	double min_diff_in_deg=0.1;
	double min_diff_in_h = min_diff_in_deg / 15.00;

	while( !bOK ){
		// normal behaviour - read next frame from the disk :

		// looking for next frame - skiping same as previous :
		BOOL_T bNextFrameOK=FALSE;
		while( !bNextFrameOK ){
			// look for frame different then previous :
			while( m_FrameIdx < m_pList->GetCount() ){
				if( strlen( m_szLastFITSFile.c_str() )>0 && 
					 strstr( m_szLastFITSFile.c_str(), (m_pList->GetListTable())[m_FrameIdx].c_str() ) ){
					printf("Frame %d/%d = %s same as previous (%s) one skipped\n",
									m_FrameIdx,(m_pList->GetListTable()).GetCount(),
									(m_pList->GetListTable())[m_FrameIdx].c_str(),
									m_szLastFITSFile.c_str());																		
					m_FrameIdx++;										
				}else{
					bNextFrameOK = TRUE;
					break;
				}
			}

			// in case already at end of list try to update list or 
			// finish analysis :
			if(m_FrameIdx>=m_pList->GetCount()){
				if( !gCCDParams.m_bAutoUpdateList ){
					// this is usual situtation
					return FALSE;
				}else{	
					// this is when updating list on-line :
					my_printf_now("at end of frames list, re-reading list file now\n");
					if( !ReReadListFile() )
						return FALSE;
				}
			}
		}
	
		printf("Current position on list : %d / %d\n",m_FrameIdx,(m_pList->GetListTable()).GetCount());

		MyParser szNewFile = (m_pList->GetListTable())[m_FrameIdx].c_str();
		CMyStrTable CCD_List;
		const char* pFile;	
		szFileNames = "";
	
		pFile = szNewFile.GetNextItem();
		if( !pFile ){
			return FALSE;
		}

		mystring fname;
		fname << m_PipelineCfg.m_szSampleFramesDir << "/" << pFile;
		CCDMatrix& new_image = GetNext()[0];
		if(!new_image.ReadFITSFile( fname.c_str() )){
			printf("ERROR could not read frame : %s\n",fname.c_str());
			printf("PLEASE CHECK fits file : %s\n",fname.c_str());					

			if( !gCCDParams.m_bSkipBadFrames ){
				return FALSE;
			}else{
				printf("trying to read next frame from list ...");
			}
		}
		if(!szFileNames.empty())
			szFileNames << ",";
		szFileNames << fname;

		time_t ut_time = (time_t)(new_image.getObsTime( TRUE ));
		const char* szCurrObject = new_image.getKey( OBJECT );
      const char* szCurrRA     = new_image.getKey( RA_OBS );
      const char* szCurrDEC    = new_image.getKey( DEC_OBS );

		// if summation of N frames is required :
	   double curr_ra=0,curr_dec=0;
	   if( szCurrRA ){
   		curr_ra = atof( szCurrRA );
      }
	   if( szCurrDEC ){
   		curr_dec = atof( szCurrDEC );
      }

		// finding current and previous field names :
	   mystring szPrevField=szPrevObject.c_str(),szCurrField=szCurrObject;
   	if( strlen( szCurrObject ) >= 1 )
     		szCurrField = szCurrObject+1;
	   if( strlen( szPrevObject.c_str() ) >= 1 )
   		szPrevField = szPrevObject.c_str() + 1;

		
		BOOL_T bCoordChange = FALSE;

		if( strlen( szPrevField.c_str() ) && 
				( strcmp( szPrevField.c_str(), szCurrField.c_str() ) || 
				  fabs( curr_ra - prev_ra )>=min_diff_in_h ||
				  fabs( curr_dec - prev_dec )>=min_diff_in_deg
				)
			){
			bCoordChange = TRUE;
		}


		if( bCoordChange || aver_file_list.size()==n_aver ){
			if( bCoordChange ){
				printf("field changed %s -> %s\n",szPrevField.c_str(),szCurrField.c_str());
				printf("coordinates changed : (%.2f,%.2f) -> (%.2f,%.2f)\n",
						prev_ra,prev_dec,curr_ra,curr_dec);
			}

// temporary check - in case change due to coordinates check 
// use last frames sum - no matter how many frames where averaged :
			if( aver_file_list.size() != n_aver ){
				printf("using average of only %d images (expected %d)\n",
							aver_file_list.size(),n_aver);
			}
//			if( aver_file_list.size() == n_aver ){
				printf("%d frames summed - averaging now ...\n",n_aver);

				CCDMatrix sav( new_image.GetXSize(), new_image.GetYSize() );
				CCDUtil::CalcAndSaveAver( aver_file_list, aver_file_list, 0, 
										(aver_file_list.size()-1), m_pCurrSumOfSeveral, 
									  sav, NULL, NULL, minTime, maxTime,
									  szProcessedList, szOutDir, m_SumedFrameIndex,
									  FALSE, FALSE, FALSE, x_size, y_size,
									  FALSE, NULL, FALSE, FALSE, szOutList,
									  m_PipelineCfg.m_szSampleFramesDir.c_str() );
				new_image = sav;
				bOK = TRUE;
//			}

			// cleaning prev object :
			szPrevObject = "";

			m_pCurrSumOfSeveral->SetData(0);
			aver_file_list.clear();
			minTime=get_dttm()+10000;
			maxTime=0;

			if( bOK ){
				// average frame build :
				return TRUE;
			}
		}

		printf("adding frame : %s\n",pFile );
      CCDUtil::AddFrame( *m_pCurrSumOfSeveral, new_image );
      aver_file_list.Add( pFile );
                                                                                
      szPrevObject = szCurrObject;
      prev_ra = curr_ra;
      prev_dec = curr_dec;

		// after saving - so that frame for new series is not taken into account
      if( ut_time < minTime ){
         minTime = ut_time;
      }
      if( ut_time > maxTime ){
         maxTime = ut_time;
      }

		m_FrameIdx++;


		// next frames read now check if already collected required number of frames :
	}
}

BOOL_T CCDPipeline::GetNextFrame_MC(mystring& szFileNames,BOOL_T bUseCache)
{
	if( gCCDParams.m_bOnSumedFrames && gCCDParams.m_bAutoCalcSum ){
		// in case of automatic sumation caculate sum when 
		// proper number of frames collected and return :
		return GetNextFrame_MC_AutoSum( szFileNames );
	}

	mystring szBaseFileNames;

	if(bUseCache){		
		szFileNames = "";
		// in case stored in memory 
		if( !m_pFrameCache ){
			printf("Reading all frames from the disk ...\n");
			// intializing cache :
			m_pFrameCache = new cCCDCache();
			m_pFrameCache->Init( m_SizeX, m_SizeY, m_nCCD, m_pList->GetCount() );

			for(int i=0;i<m_pList->GetCount();i++){
				MyParser szNewFile = (m_pList->GetListTable())[i].c_str();
				CMyStrTable CCD_List;
		      const char* pFile;
				while(pFile = szNewFile.GetNextItem()){
					CCD_List.Add(pFile);	
				}
				for(int j=0;j<CCD_List.GetCount();j++){
					mystring fname;
					fname << gCCDParams.m_szSampleFramesDir << "/" << CCD_List[j];
					printf("DIR : %s\n",gCCDParams.m_szSampleFramesDir.c_str());
					if( ((m_pFrameCache->m_pCache)[i])[j].ReadFITSFile(fname.c_str()) ){
						_CCDLIB_PRINTF_("Frame %d-%d initialized from file : %s\n",i,j,fname.c_str());
						if(gCCDParams.m_bDarkFrame && m_pDarkFrame){
							((m_pFrameCache->m_pCache)[i])[j].Subtract( (*m_pDarkFrame)[j], ((m_pFrameCache->m_pCache)[i])[j], TRUE /* zero if newgative */ );
						}
					}else{
						_CCDLIB_PRINTF_("Could not read frame %s, exiting ...\n",fname.c_str());
					}
				}
			}
			m_pFrameCache->m_CurIdx = 0;

			// this is to turn off dark frame subtraction - this 
			// is done only once - here when loading frames to cache
			gCCDParams.m_bDarkFrame = FALSE;
			gCCDParams.SetParam("CCD_DARK_FRAME","0");
		}
		cCCD* CurrentFrame = m_pFrameCache->GetCurrent();			
		if(CurrentFrame){
			szFileNames = (m_pList->GetListTable())[m_pFrameCache->m_CurIdx].c_str();
			GetNext() = (*CurrentFrame);
			(*m_pFrameCache)++;
		}else{
			m_pFrameCache->Reset();
			return FALSE;
		}
	}else{

		BOOL_T bOK=FALSE;

		while( !bOK ){
			// normal behaviour - read next frame from the disk :

			// looking for next frame - skiping same as previous :
			BOOL_T bNextFrameOK=FALSE;
			while( !bNextFrameOK ){

				// look for frame different then previous :
				while( m_FrameIdx < m_pList->GetCount() ){
					if( strlen( m_szLastFITSFile.c_str() )>0 && 
						 strstr( m_szLastFITSFile.c_str(), (m_pList->GetListTable())[m_FrameIdx].c_str() ) ){
						printf("Frame %d/%d = %s same as previous (%s) one skipped\n",
										m_FrameIdx,(m_pList->GetListTable()).GetCount(),
										(m_pList->GetListTable())[m_FrameIdx].c_str(),
										m_szLastFITSFile.c_str());																		
						m_FrameIdx++;										
					}else{
						bNextFrameOK = TRUE;
						break;
					}
				}

				// in case already at end of list try to update list or 
				// finish analysis :
				if(m_FrameIdx>=m_pList->GetCount()){
					if( !gCCDParams.m_bAutoUpdateList ){
						// this is usual situtation
						return FALSE;
					}else{	
						// this is when updating list on-line :
						my_printf_now("at end of frames list, re-reading list file now\n");
						if( !ReReadListFile() )
							return FALSE;
					}
				}
			}
	
			printf("Current position on list : %d / %d\n",m_FrameIdx,(m_pList->GetListTable()).GetCount());

			MyParser szNewFile = (m_pList->GetListTable())[m_FrameIdx].c_str();
			CMyStrTable CCD_List;
			const char* pFile;	
			szFileNames = "";
	
			while(pFile = szNewFile.GetNextItem()){
				CCD_List.Add(pFile);	
			}
			for(int i=0;i<CCD_List.GetCount();i++){
				mystring fname;
				// fname << gCCDParams.m_szSampleFramesDir << "/" << CCD_List[i];
				// change on 2004-10-17 - for SLT eventanal - reading of event parts :
				fname << m_PipelineCfg.m_szSampleFramesDir << "/" << CCD_List[i];
				if(!GetNext()[i].ReadFITSFile( fname.c_str() )){
					printf("ERROR could not read frame : %s\n",fname.c_str());
					printf("PLEASE CHECK fits file : %s\n",fname.c_str());					

					if( !gCCDParams.m_bSkipBadFrames ){
						return FALSE;
					}else{
						printf("trying to read next frame from list ...");
					}
				}else{
					const char* szTemp = GetNext()[i].getKeyValue( CHIP_TEMP );
					printf("ChipTemp on image %s T = %s\n",fname.c_str(),szTemp);
					double chip_temp = atof(szTemp);
					if( chip_temp >= gCCDParams.m_SkipWarmImages ){
						printf("Image too warm (T=%.2f), skiped\n",chip_temp);						
						bOK = FALSE;
					}else{
						printf("ChipTemp = %.2f - OK ( limit = %.2f )\n",chip_temp,gCCDParams.m_SkipWarmImages);
						bOK=TRUE;
					}
				}
				if(!szFileNames.empty())
					szFileNames << ",";
				szFileNames << fname;
				if(strlen(szBaseFileNames.c_str())){
					szBaseFileNames << ",";
				}
				szBaseFileNames << CCD_List[i].c_str();
			}
			m_FrameIdx++;
		}
	}

	m_szLastFITSFile = szFileNames.c_str();

	// flip frame if needed :
	//if( m_PipelineCfg.m_bDriverReverseImage!=eReverseImageNone ){
	//	GetNext()[0].FlipImage( m_PipelineCfg.m_bDriverReverseImage );
	//}
	

	return TRUE;
}



BOOL_T CCDPipeline::PrepareNewFrame()
{
	cCCD& cNew = GetCurrent();
	if(m_pDarkFrame && m_PipelineCfg.m_bDarkFrame){
		PROFILER_START
		cNew.Subtract(*m_pDarkFrame,cNew,TRUE);
		PROFILER_END("Subtraction of dark took : ");
	}

	if(  m_PipelineCfg.m_bShutterCorr ){
		PROFILER_START
		cNew[0].Data_Correction( 0.001 , 10.00 );
		PROFILER_END("Correction for shutter opened effect took :");
	}

	if(m_pFlatFrame && m_PipelineCfg.m_bFlatFrame){
		// cNew.Divide( *m_pFlatFrame , cNew );
		//for(register int i=0;i<m_nCCD;i++){
		//	cNew[i].Multiply( (*m_pFlatFrame)[i], cNew[i], (*m_pFlatFrame)[i].GetMedianValue() );
		//}				
		// CCDUtil::Divide( cNew[0], (*m_pFlatFrame), m_PipelineCfg.m_nIgnoreEdge );
		CCDUtil::Divide( cNew[0], (*m_pFlatFrame) );
	}
	if(m_pHomeopaticFrame){
		ShiftHomeopaticFrame();
	}
	if(m_pLaplaceFrame){
		ShiftFrame( m_pLaplaceFrame );
	}
	if(gCCDParams.m_bCheckLaplaceCondition || gCCDParams.m_bKeepLaplaceFrame){
		CalcLaplaceOfNew();	
	}
	if(gCCDParams.GetCalcBackgrFlag()){
		UpdateBackgroundMap();
	}
	if(gCCDParams.m_bSubtractBackground){
		SubtractBackground();
	}

	return TRUE;
}



void CCDPipeline::PrintEventDesc( CccdReport& eventReport, int Index,
                                  LONG_T xSize, LONG_T ySize, LONG_T FrameIdx )
{
	mystring szEventDesc(2000,2000);	
	eventReport.GetDetailEventDesc( szEventDesc, Index, xSize, this );
	printf("%s\n",szEventDesc.c_str());fflush(0);
}


// dumps events confirmed - index for frame to be printed
// depends on number of next frames required for confirmation
void CCDPipeline::DumpConfirmedEvents( LONG_T FrameIndex )
{
	LONG_T pos=-1;
   cCCD* pFrame = FindFrame( FrameIndex, pos );	
   if(pFrame){
		LONG_T MatrixCount = m_pCCD[0].GetCount();
		mystring szIndex;

		szIndex << FrameIndex;
      for(int i=0;i<MatrixCount;i++){			
			CCDEventList EventsToShow = pFrame->GetMatrix(i)->GetGenEvents();
			EventsToShow += pFrame->GetMatrix(i)->GetFoundEvents();
		
			DumpNewEvents( NULL, EventsToShow, szIndex.c_str(), "" );
		}
	}
}



void CCDPipeline::PrepareFoundEventsList()
{
	CFrameEvents newFrameEvents( m_FrameCounter );
	m_allFoundEvents.push_back( newFrameEvents );
}

void CCDPipeline::ClearNewEvents()
{
	if(m_allFoundEvents.size()>0){
		CFrameEvents::iterator i;
		for(i=m_allFoundEvents.back().begin();i!=m_allFoundEvents.back().end();i++){
      	i->clear();
	   }
	}
}

int CCDPipeline::AddFoundEventsToList()
{	
	cCCD& currFrame = GetCurrent();

	int nCount=0;
	for(register int i=0;i<currFrame.GetCount();i++){
		// printf("BEFORE ADD : m_allFoundEvents.back().size() = %d\n",m_allFoundEvents.back().size());

		m_allFoundEvents.back().push_back( currFrame[i].GetFoundEvents() );				
		// printf("AFTER ADD : m_allFoundEvents.back().size() = %d\n",m_allFoundEvents.back().size());
		nCount += currFrame[i].GetFoundEvents().size();
	}
	return nCount;
}


void CCDPipeline::UpdateGenEventsList()
{
	cCCD& currFrame = GetCurrent();
	/*m_allGenEvents.back().clear();
	for(register int i=0;i<currFrame.GetCount();i++){		
		m_allGenEvents.back().push_back( currFrame[i].GetGenEvents() );
	}*/

	for(register int i=0;i<currFrame.GetCount();i++){
		CCDEventList& genListUpdt = currFrame[i].GetGenEvents();
		CCDEventList& genListOld = m_allGenEvents.back()[i];


		for(CCDEventList::iterator evt=genListOld.begin();evt!=genListOld.end();evt++){
			CccdReport* pEvtUpdt = genListUpdt.FindEvent( (int)evt->m_Point.x, (int)evt->m_Point.y );
			if( pEvtUpdt ){
				evt->m_PixelAnalResults.laplaceSum = pEvtUpdt->m_PixelAnalResults.laplaceSum;
				evt->m_PixelAnalResults.maxAverageOfPrev = pEvtUpdt->m_PixelAnalResults.maxAverageOfPrev;
			}
		}
	}
}

void CCDPipeline::AddGenEventsToList()
{
	cCCD& currFrame = GetCurrent();

	CFrameEvents newFrameEvents( m_FrameCounter );

	m_allGenEvents.push_back( newFrameEvents );
	for(register int i=0;i<currFrame.GetCount();i++){
		m_allGenEvents.back().push_back( currFrame[i].GetGenEvents() );				
	}
	
}


int CCDPipeline::CompileEventsReportFromFile( LONG_T& genCount,LONG_T& foundCount, 
														    LONG_T& genIdent, CCDEventList& genEventsList )
{
	mystring szFoundLog,szGenLog,szReGenLog;
	// szFoundLog << gCCDParams.GetOutputDir() << "/" << m_VerifiedEventsLog.m_szFileName.c_str();
	szFoundLog << gCCDParams.GetOutputDir() << "/" << m_FinalEventsLog.m_szFileName.c_str();
	szGenLog << gCCDParams.GetOutputDir() << "/" << m_GenEventsLog.m_szFileName.c_str();
	szReGenLog << gCCDParams.GetOutputDir() << "/" << m_ReGenEventsLog.m_szFileName.c_str();

	int ret = CompileEventsReportFromFile( szFoundLog,szGenLog,szReGenLog,
														genCount, foundCount, genIdent,	
														genEventsList );
	return ret;
}

int CCDPipeline::CompileEventsReportFromFile( mystring& szFoundLog,
															 LONG_T& genCount,LONG_T& foundCount, 
														    LONG_T& genIdent, CCDEventList& genEventsList )
{
	mystring szGenLog,szReGenLog;
	szGenLog << gCCDParams.GetOutputDir() << "/" << m_GenEventsLog.m_szFileName.c_str();
	szReGenLog << gCCDParams.GetOutputDir() << "/" << m_ReGenEventsLog.m_szFileName.c_str();

	int ret = CompileEventsReportFromFile( szFoundLog,szGenLog,szReGenLog,
														genCount, foundCount, genIdent,	
														genEventsList );
	return ret;
}


int CCDPipeline::CompileEventsReportFromFile( mystring& szFoundLog,
                                              mystring& szGenLog, mystring& szReGenLog,
															 LONG_T& genCount,
															 LONG_T& foundCount, 
															 LONG_T& genIdent, 
															 CCDEventList& genEventsList )
{ 
	mystring szDTM = get_date_time_string();
	printf("Compiling event eff/bkg report starting at : %s\n",szDTM.c_str());fflush(0);

	int frame_index = 0;
	BOOL_T bContinue=TRUE; 	
	const char* pLine = NULL;

	BOOL_T bRead=TRUE;
	BOOL_T bTakeNew=TRUE;
	CccdReport genEvt,foundEvt;
	CCDEventList reGenEvents,genUsed;	

	genEventsList.ReadEvents( szGenLog.c_str() );
	reGenEvents.ReadEvents( szReGenLog.c_str(), FALSE );

	BOOL_T bSortOK=FALSE;
	if( genEventsList.CheckSortByFrame() ){
		printf("Generated events list sorted by frame - OK\n");
		bSortOK = TRUE;
	}else{
		printf("Generated events list not-sorted by frame - compilation may be very slow !!!\n");
	}


	genCount = genEventsList.size();
   foundCount = 0;
   genIdent = 0;
	
	CCDEventList foundEventsList;
	foundEventsList.Read( szFoundLog.c_str() );

	for(int i=0;i<foundEventsList.size();i++){
		CccdReport& foundEvt = foundEventsList[i];

		CCDEventList reGenList;
		reGenEvents.GetEventsByFrameIndex( reGenList, foundEvt.m_DayFrameIndex );

		BOOL_T bReGenEvt=FALSE;
		CCDEventList::iterator pReGen;
		for(pReGen = reGenList.begin();pReGen!=reGenList.end();pReGen++){
			if( CCDMatrix::IsGenEventIdentified( *pReGen, foundEvt ) )
				bReGenEvt=TRUE;
		}
		
		if(!bReGenEvt){
			foundEvt.m_bGenerated = FALSE;

			if( bSortOK ){
				CCDEventList::iterator frame_it = lower_bound( genEventsList.begin() , 
												genEventsList.end() ,
												foundEvt , 
												event_finder() );
				for(CCDEventList::iterator it=frame_it;
					(it!=genEventsList.end() && (*it).m_DayFrameIndex == foundEvt.m_DayFrameIndex);
					it++){
//					printf("%d %d %d\n",(*it).m_DayFrameIndex,(int)(*it).m_MaxPoint.x,(int)(*it).m_MaxPoint.y);
	
					CccdReport& genEvt = (*it);			
					if( !genUsed.FindEvent( genEvt.m_MaxPoint.x, genEvt.m_MaxPoint.y, genEvt.m_DayFrameIndex, 1 ) && 
						 genEvt.m_DayFrameIndex==foundEvt.m_DayFrameIndex &&
						 CCDMatrix::IsGenEventIdentified( genEvt, foundEvt ) ){
						genEvt.m_bIdentified = TRUE;
						foundEvt.m_bGenerated = TRUE;
						genUsed.push_back( genEvt );
						break;					
					}				
				}
			}else{				
				for(int k=0;k<genEventsList.size();k++){
					CccdReport& genEvt = genEventsList[k];

					if( !genUsed.FindEvent( genEvt.m_MaxPoint.x, genEvt.m_MaxPoint.y, genEvt.m_DayFrameIndex, 1 ) && 
						 genEvt.m_DayFrameIndex==foundEvt.m_DayFrameIndex &&
						 CCDMatrix::IsGenEventIdentified( genEvt, foundEvt ) ){
						genEvt.m_bIdentified = TRUE;
						foundEvt.m_bGenerated = TRUE;
						genUsed.push_back( genEvt );
						break;					
					}
				}
			}

			if(foundEvt.m_bGenerated){
				genIdent++;
			}else{
				if(!gCCDParams.m_bOnlySuperNewBackgr || 
					foundEvt.m_PixelAnalResults.eventType==eBrighten ){

					_TRACE_PRINTF_5("BACGR : %d (%d,%d)\n",(int)foundEvt.m_DayFrameIndex,
							(int)foundEvt.m_MaxPoint.x,(int)foundEvt.m_MaxPoint.y );

					foundCount++;
				}
			}
		}
				
	}

	szDTM = get_date_time_string();
	printf("Compilation of report finished at : %s\n",szDTM.c_str());fflush(0);

	
	return (genIdent+foundCount);	
}


// this works in different way - goes over generated events list 
// and finds identified events from list :
int CCDPipeline::CompileEventsReportFromFileNew( 
															 mystring& szFoundLog,
															 LONG_T& genCount,
															 LONG_T& foundCount, 
															 LONG_T& genIdent, 
															 CCDEventList& genEventsList )
{ 
	mystring szGenLog,szReGenLog;
	szGenLog << gCCDParams.GetOutputDir() << "/" << m_GenEventsLog.m_szFileName.c_str();
	szReGenLog << gCCDParams.GetOutputDir() << "/" << m_ReGenEventsLog.m_szFileName.c_str();


	int frame_index = 0;
	BOOL_T bContinue=TRUE; 	
	const char* pLine = NULL;

	BOOL_T bRead=TRUE;
	BOOL_T bTakeNew=TRUE;
	CccdReport genEvt;
	CCDEventList reGenEvents,foundEvts;

	if( genEventsList.size()==0 ){
		printf("WARNING : re-reading genEventsList !!! in CCDPipeline::CompileEventsReportFromFileNew(\n");
		genEventsList.ReadEvents( szGenLog.c_str() );
	}
	reGenEvents.ReadEvents( szReGenLog.c_str(), FALSE );
	foundEvts.ReadEvents( szFoundLog.c_str() );

	genCount = genEventsList.size();
   foundCount = 0;
   genIdent = 0;

	for(int i=0;i<genEventsList.size();i++){
		if( gCCDParams.m_bSumAllMethodsInNparTest ){
			if( genEventsList[i].m_bIdentified ){
				// in case sum of all methods is determined, add only 
				// those events which are not identified by previous methods :
				continue;
			}
		}


		CCDEventList foundEvtsOnFrame;
		int end_frame = genEventsList[i].m_FrameIndex+gCCDParams.m_nPutObjOnNFrames;
		int rest = ( end_frame % gCCDParams.m_bKeepSumOfPrevNFrames);
		if( rest!=0 ){
			end_frame += ( gCCDParams.m_bKeepSumOfPrevNFrames - rest );
		}

		//foundEvts.GetEventsByFrameIndex( foundEvtsOnFrame, genEventsList[i].m_FrameIndex,
		//									end_frame );

		BOOL_T bFound=FALSE;
		CCDEventList::iterator pFound;

		for( pFound = foundEvts.begin();pFound!=foundEvts.end();pFound++){
			// for(pFound = foundEvtsOnFrame.begin();pFound!=foundEvtsOnFrame.end();pFound++){

			if( pFound->m_FrameIndex>=genEventsList[i].m_FrameIndex && pFound->m_FrameIndex<=end_frame ){
				if( CCDMatrix::IsGenEventIdentified( *pFound, genEventsList[i] ) ){
					bFound = TRUE;
					pFound->m_bGenerated = TRUE;
				}
			}
		}
						
		if( bFound ){
			genIdent++;
		}
	}	

	// foundCount = ( foundEvts.size() - genIdent );
	CCDEventList::iterator pFound;
	for(pFound = foundEvts.begin();pFound!=foundEvts.end();pFound++){
		if( !pFound->m_bGenerated ){
			foundCount++;
		}
	}

	return foundEvts.size();	
}



int CCDPipeline::CompileEventsReport(LONG_T& genCount,LONG_T& foundCount,
												 LONG_T& genIdent)
{
	deque<CFrameEvents>::iterator f;
	deque<CFrameEvents>::iterator g;

	m_CompiledEventsList.clear();
	for(f=m_allFoundEvents.begin(),g=m_allGenEvents.begin();
		f!=m_allFoundEvents.end() && g!=m_allGenEvents.end();f++,g++){
		Assert(f->size()==g->size(),"Number of gen frames differs from found frames (%d!=%d)",g->size(),f->size());
		
		CFrameEvents newFrame( f->m_FrameCounter );
		m_CompiledEventsList.push_back( newFrame );
		for(int j=0;j<f->size();j++){
			(m_CompiledEventsList.back()).push_back( m_EmptyEventsList );
			CCDEventList& compEvents = m_CompiledEventsList.back().back();


			printf("FOUND COUNT=%d, GEN COUNT=%d\n",f->size(),g->size());
			CCDMatrix::CompileEventReport( compEvents, (*f)[j], (*g)[j], NULL, -1 , FALSE, FALSE );
		}
	}

	genCount = 0;
	foundCount = 0;
	genIdent = 0;
	vector<CFrameEvents>::iterator i;
	for(i=m_CompiledEventsList.begin();i!=m_CompiledEventsList.end();i++){
		CFrameEvents::iterator pFrameEvents;
		for(pFrameEvents=i->begin();pFrameEvents!=i->end();pFrameEvents++){
			CCDEventList::iterator pEvt;
			for(pEvt=pFrameEvents->begin();pEvt!=pFrameEvents->end();pEvt++){
				if(pEvt->m_bGenerated){
					genCount++;
					if(pEvt->m_bIdentified){
						genIdent++;
					}
				}else{
					if(!(pEvt->m_PixelAnalResults).m_bRejectedDueToTrack)
						foundCount++;
				}				
			}	
		}		
	}

	return (genIdent+foundCount);
}

void CCDPipeline::SaveAllFrames()
{
	if(m_FrameCounter>0){
		mystring szBaseFileName;
		szBaseFileName << gCCDParams.GetOutputDir() << "/Events/Triggers//pipeline" << m_PipelineIndex << "_frame";

		for(int i=0;i<m_nPipelineSize;i++){
			if( m_pCCD[i].GetFrameIndex()>=0 ){
				mystring szFName;
				szFName << szBaseFileName << m_pCCD[i].GetFrameIndex() << "_ccd";
				for(int j=0;j<m_nCCD;j++){
	      	   CCDMatrix& matrix = m_pCCD[i][j];
					szFName << j << ".fit";
					matrix.WriteToFITSFile( szFName.c_str() );
				}				
			}
   	}
	}
}

void CCDPipeline::SaveCalibFrame()
{
	SaveCurrentFrame( "/Calibration" , eSaveCalibPicture );
}

mystring CCDPipeline::SaveCurrentFrame( const char* szSubDirName /*=/Events/Triggers"*/, eSavePictureType_T saveType )
{
	mystring szRet;
	if(m_FrameCounter>0){
		// saving current frame :
		cCCD& currFrame = GetCurrent();

		if( saveType==eSaveCalibPicture ){
			// get RA and DEC from mount and add to saved header
		}


		CCDMatrix& matrix = currFrame[0];
		mystring szFname,szBase;
		szBase << "frame_save_req" << GetFrameIndex() << ".fit";
		GetOutFileName( szFname, szSubDirName, szBase.c_str(), this, -1, FALSE );
		matrix.WriteToFITSFile( szFname.c_str() );		
		szRet = szFname;
	}
	return szRet;
}

void CCDPipeline::SaveEventsToFITS()
{
if( m_DayFrameCounter>=415 )
	printf("odo");

	
	// now check if save found events to disk - but only those confirmed on
	// next gCCDParams.m_ConfirmEventsOnNextNFrames frames :
	if(gCCDParams.m_bSaveFramesWithEvents){
		// in case no confirmation requierd on next frames m_ConfirmEventsOnNextNFrames=0 - ALWAYS !!!
		Assert(gCCDParams.m_ConfirmEventsOnNextNFrames>=0,"m_ConfirmEventsOnNextNFrames must be >=0 !!!");
		int ConfirmedFrameIndex = (m_allFoundEvents.size()-1)-gCCDParams.m_ConfirmEventsOnNextNFrames;

		if(ConfirmedFrameIndex<0)
			return;

		CFrameEvents& framesEvents = m_allFoundEvents[ConfirmedFrameIndex];
		

		for(int cam=0;cam<framesEvents.size();cam++){
			SaveEventsToFITS( framesEvents[cam], TRUE );
			WriteEventDetails( framesEvents[cam] );
		}


		int dumpPrevStart=MAX( (ConfirmedFrameIndex-gCCDParams.m_nSaveFramesBeforeAndAfter), 0);
		for(int frame=dumpPrevStart;frame<ConfirmedFrameIndex;frame++){		
			CFrameEvents& prevFramesEvents = m_allFoundEvents[frame];

			for(int cam=0;cam<prevFramesEvents.size();cam++){
				SaveEventsToFITS( prevFramesEvents[cam], FALSE );
			}
		}		
	}
}

int CCDPipeline::GetFinalCoicReport(  CCDEventList& eventCCD1, CCDEventList& eventCCD2,
												  CCDEventList& finallist1, CCDEventList& finallist2 )
{
	static int evt_id=1;

	finallist1.clear();
	finallist2.clear();
	for(int i=0;i<eventCCD1.size();i++){
		CccdReport* pEvent2 = eventCCD2.FindEvent( (int)eventCCD1[i].m_PointTransformed.x,
																 (int)eventCCD1[i].m_PointTransformed.y,
																 (int)eventCCD1[i].m_DayFrameIndex,
																 (int)gCCDParams.m_nCoicRedial );									


		if( pEvent2 && pEvent2->IsIdentified() && eventCCD1[i].IsIdentified() && 
			 CCD_Analyser::VerifyEvent(eventCCD1[i]) && CCD_Analyser::VerifyEvent(*pEvent2) ){

			if( !gCCDParams.GetMC() ){
				mystring szID;
				eventCCD1[i].GetEvtID( szID );
				eventCCD1[i].m_EvtID = szID;
				pEvent2->m_EvtID = szID;
			}else{
				char szID[64];
				sprintf(szID,"%.11d",evt_id);
				eventCCD1[i].m_EvtID = szID;
				pEvent2->m_EvtID = szID;
				evt_id++;
			}

			finallist1.push_back( eventCCD1[i] );
			finallist2.push_back( *pEvent2 );
		}
	}
	return finallist1.size();
}

int CCDPipeline::DumpFinalCoicReport( CCDPipeline* pPipeline1, CCDPipeline* pPipeline2,
												  CCDEventList& eventCCD1, CCDEventList& eventCCD2 )
{
	CCDEventList finallist1,finallist2;
	int ret = DumpFinalCoicReport( pPipeline1, pPipeline2, eventCCD1, eventCCD2,
											 finallist1, finallist2 );

	return ret;
}

int CCDPipeline::DumpFinalCoicReport( CCDPipeline* pPipeline1, CCDPipeline* pPipeline2,
												  CCDEventList& finallist1, CCDEventList& finallist2,
												  const char* log1, const char* log2,
												  BOOL_T bAsFinalLog/*=TRUE*/ )
{
	printf("CCDPipeline::DumpFinalCoicReport %d-%d\n",finallist1.size(),finallist2.size());
	if( finallist1.size() > 0 ){
		finallist1.CalcConeDist();

		// here check for previous GRB events and SUPERNOVA and other events :
		CheckExternalTriggers( finallist1, finallist2 );

		// then check SN events :
		// ... 

		// then check Galaxies :
		// ...


		if( bAsFinalLog ){
			CCDEventList::DumpFinalEvents( log1, finallist1 );
			CCDEventList::DumpFinalEvents( log2, finallist2 );
		}else{
			CCDEventList::DumpEventReport( log1, finallist1, -1, FALSE, TRUE, FALSE );
			CCDEventList::DumpEventReport( log2, finallist2, -1, FALSE, TRUE, FALSE );
		}

		/*if( m_WorkingMode.m_WorkingMode == eDAQSatTriggerMode ){
			DumpToTriggerLog( pPipeline1, finallist1 );
			DumpToTriggerLog( pPipeline2, finallist2 );
		}


		if( gCCDParams.m_bUseDB ){
			InsertEventsToDB( finallist1, finallist2, pPipeline1, pPipeline2 );
		}

		if( gCCDParams.m_bLogFrameStat ){
			pPipeline1->SetFinalRates( finallist1 );
			pPipeline2->SetFinalRates( finallist2 );
		}*/
	}
	return finallist1.size();
}


CMyMutex gFinalLogDumpLock;
int CCDPipeline::DumpFinalCoicReport()
{
	int ret=0;
	gFinalLogDumpLock.Lock();
	if( gCCDParams.m_bCCDDouble && !m_bFinalDumped ){
		if( m_PipelineList.size()==2 ){
			my_printf_now("FINAL_EVENTS : dumping final events log\n");
			ret = DumpFinalCoicReport( m_PipelineList[0] ,  m_PipelineList[1] );

			// dumping sumed events :
			if( gCCDParams.m_bAnalyzeSumOfPrevNFrames ){
				DumpSumEvents( m_PipelineList[0] ,  m_PipelineList[1] );
			}
		}
	}
	m_bFinalDumped = TRUE;
	gFinalLogDumpLock.UnLock();	

	return ret;
}

int CCDPipeline::DumpFinalCoicReport( CCDPipeline* pPipeline1, CCDPipeline* pPipeline2 )
{
	int ret = DumpFinalCoicReport( pPipeline1, pPipeline2,
											 pPipeline1->m_ConfirmedEvents,
											 pPipeline2->m_ConfirmedEvents,
											 pPipeline1->m_TriggerEventsList,
											 pPipeline2->m_TriggerEventsList );
	return ret;
}

int CCDPipeline::DumpFinalCoicReport( CCDPipeline* pPipeline1, CCDPipeline* pPipeline2,
												  CCDEventList& eventCCD1, CCDEventList& eventCCD2,
												  CCDEventList& finallist1, CCDEventList& finallist2 )
{
	printf("CCDPipeline::DumpFinalCoicReport2 %d-%d\n",finallist1.size(),finallist2.size());

	GetFinalCoicReport( eventCCD1, eventCCD2, finallist1,finallist2 );

	if( finallist1.size() > 0 ){
		finallist1.CalcConeDist();

		CCDPipeline::m_nAlertsCount += finallist1.size();

		// here check for previous GRB events and SUPERNOVA and other events :
		CheckExternalTriggers( finallist1, finallist2 );

		// then check SN events :
		// ... 

		// then check Galaxies :
		// ...

		pPipeline1->DumpFinalEvents( &(pPipeline1->m_FinalEventsLog), finallist1 );
		pPipeline2->DumpFinalEvents( &(pPipeline2->m_FinalEventsLog), finallist2 );	

/*		if( m_WorkingMode.m_WorkingMode == eDAQSatTriggerMode ){
			DumpToTriggerLog( pPipeline1, finallist1 );
			DumpToTriggerLog( pPipeline2, finallist2 );
		}*/


		if( gCCDParams.m_bUseDB ){
			InsertEventsToDB( finallist1, finallist2, pPipeline1, pPipeline2 );
		}

		if( gCCDParams.m_bLogFrameStat ){
			pPipeline1->SetFinalRates( finallist1 );
			pPipeline2->SetFinalRates( finallist2 );
		}
	}
	return finallist1.size();
}


int CCDPipeline::SetFinalRates( CCDEventList& finallist )
{
	CCDEventList::iterator i;
	CFrameEventStatTab* pBackgrStat=NULL;
	if( GetAnalPtr() ){
		pBackgrStat = &(GetAnalPtr()->backgrStat);
	}
	int final=0;
	if( pBackgrStat ){
		for( i=finallist.begin();i!=finallist.end();i++){
			final += pBackgrStat->IncFinal( i->m_FrameIndex );
		}
	}
	return final;
}

void CCDPipeline::ClearConfirmedList()
{
	m_ConfirmedEvents.clear();
}


int CCDPipeline::DumpFoundEventsLog( BOOL_T bForceAll/*=FALSE*/, CCDEventList* pInternalTriggerList/*=NULL*/ )
{
	int nKeepNBack = MAX(gCCDParams.m_nNumBackFramesForTracks,gCCDParams.m_ConfirmEventsOnNextNFrames);
	deque<CFrameEvents>::iterator f;

	//if(m_allFoundEvents.empty())
	//	return 0;

	int upTo = m_allFoundEvents.size()-nKeepNBack;

	if(bForceAll){
		upTo = m_allFoundEvents.size();
		m_pAnalyser->LogAllEventRates( this );
	}

	if(pInternalTriggerList){
		// clearing list of triggers :
		pInternalTriggerList->clear();
	}
	
	// clear confirmed list 
	// change 20041015 - moved to GetNextFrame due to 
	// fact that this function can be called >1 and list was unwilingly cleaned
	// when RestartPipeline called 
	// ClearConfirmedList(); 

	int count=0;
	for(int toDump=m_LastDumpedIndex+1;toDump<upTo;toDump++){
		count += m_allFoundEvents[toDump].size();
		for(int j=0;j<m_allFoundEvents[toDump].size();j++){
			mystring szIndex;
			szIndex << m_allFoundEvents[toDump].m_FrameCounter;

			// in case confirmation on next required :
			if( gCCDParams.m_ConfirmEventsOnNextNFrames>0 ){
				if(toDump==(m_allFoundEvents.size()-1)){
					// setting that rejected on next frame (no next frame used
					// for confirmation - it means not confirmed ...)
					CCDEventList& lastFrameEvents = m_allFoundEvents[toDump][j];
					for(int r=0;r<lastFrameEvents.size();r++){
						lastFrameEvents[r].RejectByNextFrames();
					}
				}
			}

			// saving real final events to file :
			// DumpFinalEvents( &m_FinalEventsLog, m_allFoundEvents[toDump][j] );

			DumpNewEvents( &m_VerifiedEventsLog, m_allFoundEvents[toDump][j], szIndex.c_str(), "", TRUE );			
			// build dumped events list :
			// m_ConfirmedEvents.clear();			
			AddFinalEvents( m_allFoundEvents[toDump][j] );

			if( pInternalTriggerList ){
				CheckForInternalTriggers( m_allFoundEvents[toDump][j], *pInternalTriggerList );
			}
		}			


		if(toDump<m_allGenEvents.size()){
			for(int j=0;j<m_allGenEvents[toDump].size();j++){
				mystring szIndex;
      	   szIndex << m_allGenEvents[toDump].m_FrameCounter;
				DumpNewEvents( &m_GenEventsLog, m_allGenEvents[toDump][j], szIndex.c_str(), "", TRUE );
			}
		}

		m_LastDumpedIndex = toDump;
	}

	if(count>0){
		// maybe clear tables : m_allFoundEvents and m_allGenEvents here ?
		// in such case set m_LastDumpedIndex=-1
		deque<CFrameEvents>::iterator pToDel;
		deque<CFrameEvents>::iterator pFirstToDel = m_allFoundEvents.begin();
		deque<CFrameEvents>::iterator pLastToDel  = m_allFoundEvents.begin();

		int i=0;
		for(pToDel=m_allFoundEvents.begin();pToDel!=m_allFoundEvents.end();i++,pToDel++){
			if(i==m_LastDumpedIndex){
				pLastToDel = pToDel+1;
				break;
			}		
		}
		// pLastToDel - removing up to one elem before this one 

		deque<CFrameEvents>::iterator pFirstGenToDel = m_allGenEvents.begin();
		deque<CFrameEvents>::iterator pLastGenToDel  = m_allGenEvents.begin();

		i=0;
		for(pToDel=m_allGenEvents.begin();pToDel!=m_allGenEvents.end();i++,pToDel++){
			if(i==m_LastDumpedIndex){
				pLastGenToDel = pToDel+1;
				break;
			}		
		}
		// pLastGenToDel - removing up to one elem before this one

		// before deleting dump tracks :
		//if( gCCDParams.m_bLogTracks && gCCDParams.m_bCheckTracks){
		//	DumpTracks( pFirstToDel, pLastToDel );
		//}

		m_allFoundEvents.erase( pFirstToDel, pLastToDel );

		if(m_allGenEvents.size())
			m_allGenEvents.erase( pFirstGenToDel, pLastGenToDel );

		m_LastDumpedIndex=-1;

		if(gCCDParams.GetMC()){
			Assert(m_allFoundEvents.size()==m_allGenEvents.size(),"Error in cleaing procedure %d!=%d !!!",m_allFoundEvents.size(),m_allGenEvents.size());
		}


		// temporary to monitor memory :
		if( gPrintfLevel>=4){	
			int mem_usage = 0;
			printf("------------- MEMORY REPORT -------------------\n");	
			printf("# of stored CFrameDesc = %d\n",m_allFoundEvents.size());
			for(int i=0;i<m_allFoundEvents.size();i++){
				mem_usage += m_allFoundEvents[i].GetMemSize();
			}
			printf("object allFoundEvents uses %d MB of memory\n",(mem_usage/1000000));

			if((mem_usage/1000000)>20){
				for(int i=0;i<m_allFoundEvents.size();i++){
					printf("Frame %d\n",i);
					m_allFoundEvents[i].GetMemSize( TRUE );
				}
				sleep(60);
			}
			printf("-----------------------------------------------\n");
			if((mem_usage/1000000)>20){
				exit(0);
			}			
		}
	}
	m_RestartDumpCounter++;

	if( gCCDParams.m_bCheckForSUPERNEW ){
		DumpSNCoic();
	}

	return count;
}

void CCDPipeline::AddFinalEvents( CCDEventList& confirmedEvents )
{
	CCDEventList::iterator i;
	for( i=confirmedEvents.begin();i!=confirmedEvents.end();i++){
		// NEW 20050714 CCD_Analyser::VerifyEvent added :		
		if( i->IsIdentified() && CCD_Analyser::VerifyEvent(*i) ){
			m_ConfirmedEvents.push_back( *i );
		}
	}
}

void CCDPipeline::DumpTracks( deque<CFrameEvents>::iterator& start,
										deque<CFrameEvents>::iterator& end )
{
	/*deque<CFrameEvents>::iterator pToDump;
	mystring szTrackDesc( 2048, 2048 );
	for( pToDump=start ; pToDump!=end ; pToDump++){
			CTrackList::iterator t;
			if( pToDump->m_TracksOnFrame.size()>0 ){
				for(t=pToDump->m_TracksOnFrame.begin();t!=pToDump->m_TracksOnFrame.end();t++){
					szTrackDesc << t->a << " " << t->b << "\t";
					vector<CEventBaseInfo>::iterator evt;
					for(evt=t->m_EventsOnTrack.begin();evt!=t->m_EventsOnTrack.end();evt++){
						szTrackDesc << (int)evt->m_FrameIndex << "-(" << (int)evt->m_MaxPoint.x << "," 
						<< (int)evt->m_MaxPoint.y << "),";
					}
				}
				CCDLog trackLog( "%s\n", "a b Events_on_Track", "Tracks", "tracklog" );
				trackLog.DumpToFile1( this, szTrackDesc.c_str() );
			}
	}*/
}

BOOL_T CCDPipeline::CheckIfEventBelongsToTrack( double x, double y, double radius, 
						CTrackDesc& track )
{
	CTrackList& oldTracks = GetTrackList();
	CTrackList::iterator it;
   for(it=oldTracks.begin();it!=oldTracks.end();it++){
		for(int i=0;i<(it->m_EventsOnTrack).size();i++){
			CEventBaseInfo& evt = (it->m_EventsOnTrack)[i];

			if( fabs( x - (evt.m_MaxPoint).x )<=radius && fabs( y - (evt.m_MaxPoint).y )<=radius ){
				track = (*it);
				return TRUE;
			}
		}
	}	
	return FALSE;
}

BOOL_T CCDPipeline::CheckIfEventBelongsToTrack( double x, double y, int frame_index,
						double radius, 
						CTrackDesc& track, double chi2_limit, double& chi2, BOOL_T& bBelongs )
{
	CTrackList& oldTracks = GetTrackList();
	CTrackList::iterator it;
	double min_chi2 = 1000000.00;
	bBelongs = FALSE;

   for(it=oldTracks.begin();it!=oldTracks.end();it++){
		for(int i=0;i<(it->m_EventsOnTrack).size();i++){
			CEventBaseInfo& evt = (it->m_EventsOnTrack)[i];

			if( fabs( x - (evt.m_MaxPoint).x )<=radius && 
				 fabs( y - (evt.m_MaxPoint).y )<=radius && 
				 evt.m_DayFrameIndex == frame_index ){
				track = (*it);
				bBelongs = TRUE;
				chi2 = 0.00;
				return TRUE;
			}
		}

		// checking chi2 :
		it->CheckEvent( x, y, chi2 );

		if( chi2 < min_chi2 ){
			min_chi2 = chi2;			
		}

		if( chi2<chi2_limit ){
			track = (*it);
			return TRUE;
		}
	}	

	chi2 = min_chi2;

	return FALSE;
}



BOOL_T CCDPipeline::CheckIfEventBelongsToTrack( double x, double y, double radius, 
						CTrackDesc& track, CTrackList& oldTracks )
{
	CTrackList::iterator it;
   for(it=oldTracks.begin();it!=oldTracks.end();it++){
		for(int i=0;i<(it->m_EventsOnTrack).size();i++){
			CEventBaseInfo& evt = (it->m_EventsOnTrack)[i];

			if( fabs( x - (evt.m_MaxPoint).x )<=radius && fabs( y - (evt.m_MaxPoint).y )<=radius ){
				track = (*it);
				return TRUE;
			}
		}
	}	
	return FALSE;
}


void CCDPipeline::ConfirmEventsFromPreviousFrames()
{
	if(gCCDParams.m_ConfirmEventsOnNextNFrames){
		if(gCCDParams.m_ConfirmEventsOnNextNFrames>=GetPipelineSize()){
			printf("Pipeline to small to check %d next frames and confirm events, exiting ...\n",gCCDParams.m_ConfirmEventsOnNextNFrames);
			exit(0);
		}			
		


		if(m_FirstAnalysiedFrame>0 && (GetFrameIndex()>=(m_FirstAnalysiedFrame + gCCDParams.m_ConfirmEventsOnNextNFrames))){				
			LONG_T FrameToVerify = GetFrameIndex()-gCCDParams.m_ConfirmEventsOnNextNFrames;
			VerifyPreviousFrame( FrameToVerify, gCCDParams.m_ConfirmEventsOnNextNFrames );
			LONG_T pos;
			cCCD* pFrame = FindFrame( FrameToVerify, pos );
			// InitAndLockWrkEventList( (*pFrame)[0].GetFoundEvents().size()+20 );

			for(register int i=0;i<pFrame->GetCount();i++){
		     	(*pFrame)[i].CompileEventReportPtr( NULL , FrameToVerify );
			}
			// UnlockWrkEventList();
		}
	}
}

void CCDPipeline::SetNextFrameValues()
{
	int nNext=1;
	if(m_FirstAnalysiedFrame>0 && (GetFrameIndex()>=(m_FirstAnalysiedFrame + nNext))){				
		LONG_T FrameToVerify = GetFrameIndex()-1;
		SetNextFrameValues( FrameToVerify, 1 );
	}
}



// functions currently found events - not confirmed with track analysis and coicydence :
void CCDPipeline::DumpCurrentEvents()
{
	DumpNewEvents( m_FrameCounter );

	UpdateSingleCamEventsList();

	//if( m_AntyCoicEvents.size() ){
	//	DumpAntyCoicEvents();
	//}
}

void CCDPipeline::UpdateSingleCamEventsList()
{
	// m_SingleCamEvents += GetNewEvents();
	int total = m_SingleCamEvents.size() + GetNewEvents().size();
	if( total > gCCDParams.m_MaxStoredEventsFromSingleCam ){
		CCDEventList::iterator i;
		int j=0;
		for(i=m_SingleCamEvents.begin();
			 i!=m_SingleCamEvents.end() && j<(total-gCCDParams.m_MaxStoredEventsFromSingleCam);i++,j++){
		}

		// something strange : this makes a core dump :
		// m_SingleCamEvents.erase( m_SingleCamEvents.begin(), i );

		printf_now2("CCDPipeline::UpdateSingleCamEventsList, before clean : %d\n",m_SingleCamEvents.size());
		CCDEventList tmpList;
      tmpList = m_SingleCamEvents;
      m_SingleCamEvents.clear();
		for(int k=(total-gCCDParams.m_MaxStoredEventsFromSingleCam);k<tmpList.size();k++){
			m_SingleCamEvents.push_back( tmpList[k] );
		}
		printf_now2("CCDPipeline::UpdateSingleCamEventsList, after clean : %d\n",m_SingleCamEvents.size());
	}
	m_SingleCamEvents += GetNewEvents();
}

void CCDPipeline::DumpCurrentCoicEvents()
{
	mystring szMsg;
	szMsg << "CCD# " << m_PipelineIndex;


	// _TRACE_PRINTF_("\n\nCOICIDING EVENTS ON FRAME# %d\n",m_FrameCounter);
	for(int i=0;i<m_nCCD;i++){
		CCDEventList& newEvents = m_allFoundEvents.back()[i];
		mystring szIndex;
		szIndex << m_FrameCounter;
		DumpNewEvents( NULL, newEvents, szIndex.c_str(), szMsg.c_str(), FALSE );     
	}
}

void CCDPipeline::DumpAntyCoicEvents()
{
	char szLog[128];
	sprintf( szLog, ANTY_COIC_LOG, m_PipelineIndex );	
	
	mystring szFName;
	szFName << gCCDParams.GetOutputDir() << "/" << szLog;
	
	m_AntyCoicEvents.Dump( szFName.c_str() );
}

void CCDPipeline::DumpNewEvents( CCDEventList* pInternalTriggerList/*=NULL*/, BOOL_T bForceAll )
{
	DumpCurrentEvents();

	// first param - default FALSE, second list of triggers OUT
	DumpFoundEventsLog( bForceAll , pInternalTriggerList );
}

void CCDPipeline::DumpNewEvents( LONG_T Idx, const char* msg, BOOL_T bToFile/*=TRUE*/ )
{
	mystring szIndex;
	szIndex << Idx;
	DumpNewEvents( szIndex.c_str(), msg, bToFile );
}

void CCDPipeline::DumpNewEvents( const char* szIndex, const char* msg, BOOL_T bToFile /*=TRUE*/ )
{
	CCDEventList& newEvents = GetNewEvents();
	DumpNewEvents( NULL, newEvents, szIndex, msg, bToFile );
}

void CCDPipeline::DumpNewEvents( CCDEventLog* pEventLog,
											CCDEventList& newEvents, const char* szIndex, 
											const char* msg, BOOL_T bToFile  )
{
	//mystring szFname;
	//szFname << gCCDParams.GetOutputDir() << "/" << szIndex;
	// printf("\n\nDumpNewEvents=%d\n\n",gCCDParams.m_bDumpNewEvents );

	/*CCDEventList::DumpEventReport( szFileName, newEvents, -1,
											 gCCDParams.m_bStdoutOn, 
											 gCCDParams.m_bDumpNewEvents,	
											 m_bFirstDump );
	if(gCCDParams.m_bSaveFramesWithEvents){
		SaveEventsToFITS( newEvents );
		int last = GetFrameIndex()-gCCDParams.m_nSaveFramesBeforeAndAfter;
		SaveEventsToFITS( last );
	}

	m_bFirstDump = FALSE;
	if(msg && strlen(msg))
		_CCDLIB_PRINTF_("%s\n",msg);		*/

	if(!pEventLog)
		pEventLog = &m_AllEventsLog;
	pEventLog->DumpNewEvents( newEvents, szIndex, msg, bToFile );

	// dumping SN events :
	if( gCCDParams.m_bCheckForSUPERNEW ){
		m_SNEventsLog.DumpNewEvents( m_BrightenList, szIndex, "SN Events:", TRUE );
	}
}

void CCDPipeline::DumpSNCoic()
{
	if( gCCDParams.m_bCheckForSUPERNEW ){
		int max_frames_diff = 50;

		for(int i=0;i<m_BrightenVerifList.size();i++){
			CccdReport& evt = m_BrightenVerifList[i];
			
			if( evt.m_PixelAnalResults.m_bRejectedDueToTrack )
				continue;

			CTrackList bestlist;
			m_TrackList.GetBestTracks( (int)evt.m_MaxPoint.x,
												(int)evt.m_MaxPoint.y,
												evt.m_DayFrameIndex,
												5 , bestlist );

			CTrackList::iterator tr;
			for(tr=bestlist.begin();tr!=bestlist.end();tr++){
				int min_frame,max_frame;
				tr->get_frame_range( min_frame, max_frame );

				if( tr->chi2 < gCCDParams.m_MaxChi2ForPointToMatchLine && 
					 tr->bMoveOK && 
					 ( ( evt.m_DayFrameIndex>=min_frame && evt.m_DayFrameIndex<=max_frame ) || 
						abs(evt.m_DayFrameIndex-min_frame)<=max_frames_diff || 
						abs(evt.m_DayFrameIndex-max_frame)<=max_frames_diff ) ){
					evt.m_PixelAnalResults.m_bRejectedDueToTrack = TRUE;
					break;
				}
			}
		}

		// m_SNCoicEventsLog.DumpNewEvents( m_BrightenList, "0", "SN Coic Events:", TRUE );
		

		// saving parts of new events - up to current frame :
		SaveEventsToFITS( m_BrightenVerifList, FALSE, TRUE, TRUE, TRUE, TRUE );

		// saving parts :
		/*CCDEventList removed,todump;
		m_BrightenVerifList.remove_older( m_DayFrameCounter - 20, &removed );

		for(int i=0;i<removed.size();i++){
			if( !removed[i].m_PixelAnalResults.m_bRejectedDueToTrack &&
				  removed[i].m_EventType!=EVENT_TYPE_SAT ){
				todump.push_back( removed[i] );
			}
		}		
		
		m_SNCoicEventsLog.DumpNewEvents( todump, "0", "SN Coic Events:", TRUE );*/
	}
}

int CCDPipeline::DumpFinalSN( CCDPipeline* pPipeline1, 
									   CCDPipeline* pPipeline2 )
{
	CCDEventList removed1,todump1,removed2,todump2;
	(pPipeline1->m_BrightenVerifList).remove_older( pPipeline1->m_DayFrameCounter - 20, &removed1 );
	(pPipeline2->m_BrightenVerifList).remove_older( pPipeline2->m_DayFrameCounter - 20, &removed2 );

	int count=MIN(removed1.size(),removed2.size());
	for(int i=0;i<count;i++){
		if( !removed1[i].m_PixelAnalResults.m_bRejectedDueToTrack &&
			  removed1[i].m_EventType!=EVENT_TYPE_SAT && 
			 !removed2[i].m_PixelAnalResults.m_bRejectedDueToTrack &&
			  removed2[i].m_EventType!=EVENT_TYPE_SAT ){
			todump1.push_back( removed1[i] );
			todump2.push_back( removed2[i] );
		}
	}		

	// checking if limit not exceeded :	
	if( todump1.size() >= gCCDParams.m_MaxFinalSN ){
		printf("Checking if SN events limit (%d) not exceeded\n",gCCDParams.m_MaxFinalSN);

		CCDEventList final_list1,final_list2;
		CMyStrTable checked_list;

		for(int j=0;j<todump1.size();j++){
			CccdReport& evt = todump1[j];
			int frame_no = evt.m_DayFrameIndex;
			int nSameFrame=0;
			mystring szTmp;
         szTmp << frame_no;

			// skip already checked frames :
			if( checked_list.Find( szTmp.c_str() )){
				continue;
			}


			// count number of events on frame_no
			for(int k=0;k<todump1.size();k++){
				if( todump1[k].m_DayFrameIndex == frame_no ){
					nSameFrame++;
				}
			}


			// add events to out list in case not exceeds limit :
			if( nSameFrame < gCCDParams.m_MaxFinalSN ){
				printf("Frame %d - have %d events all remain on list\n",frame_no,nSameFrame);				

				for(int k=0;k<todump1.size();k++){
					if( todump1[k].m_DayFrameIndex == frame_no ){
						final_list1.push_back( todump1[k] );
					}
				}
				for(int k=0;k<todump2.size();k++){
					if( todump2[k].m_DayFrameIndex == frame_no ){
						final_list2.push_back( todump2[k] );
					}
				}
			}else{
				printf("limit exceeded, all %d events from frame %d skiped !\n",nSameFrame,frame_no);
			}

			// add frame_no to checked list
			checked_list.push_back( szTmp.c_str() );
		}

		// re-write final list :
		todump1 = final_list1;
		todump2 = final_list2;

		printf("After check %d-%d events are left\n",todump1.size(),todump2.size());
	}

	(pPipeline1->m_SNCoicEventsLog).DumpNewEvents( todump1, "0", "SN Coic Events Cam0:", TRUE );	
	(pPipeline2->m_SNCoicEventsLog).DumpNewEvents( todump2, "0", "SN Coic Events Cam1:", TRUE );

	(pPipeline1->m_SNAllCoicLog).DumpNewEvents( removed1, "0", "SN All Coic Events Cam0:", TRUE );
	(pPipeline2->m_SNAllCoicLog).DumpNewEvents( removed2, "0", "SN All Coic Events Cam0:", TRUE );
}

void CCDPipeline::DumpFinalEvents( CCDEventLog* pEventLog, CCDEventList& newEvents, 
											  BOOL_T bWriteToDB )
{
	printf("calculating cone distance for events %d events\n",newEvents.size());fflush(0);
	newEvents.CalcConeDist();
	printf("cone distance calculated\n");fflush(0);

	// event ID :
	static int evtsingle_id=1;
	for(int i=0;i<newEvents.size();i++){
		CccdReport& evt = newEvents[i];

		if( !gCCDParams.GetMC() ){
      	mystring szID;
         evt.GetEvtID( szID );
         evt.m_EvtID = szID;
      }else{
         char szID[64];
         sprintf(szID,"%.11d",evtsingle_id);
         evt.m_EvtID = szID;
         evtsingle_id++;
      }
	}

	pEventLog->DumpFinalEvents( newEvents );

	if( bWriteToDB ){
		InsertEventsToDB( newEvents, this );
	}
}

int CCDPipeline::CheckForInternalTriggers( CCDEventList& verEvents, CCDEventList& triggerList )
{
	CCDEventList::iterator pEvt;
	int count=0;
	for( pEvt=verEvents.begin();pEvt!=verEvents.end();pEvt++){
		if( CCD_Analyser::CheckIfInternalTrigger( *pEvt )){
			count++;
			triggerList.push_back( *pEvt );
		}
	}
	return count;
}

int CCDPipeline::CheckCoicInternalTriggers( CCDPipeline& pipeline1, CCDPipeline& pipeline2 )
{
	pipeline1.m_TriggerEventsList.clear();
	pipeline2.m_TriggerEventsList.clear();
	
	register int CamNo = (pipeline1.GetCurrent()).GetCount();
   for(register int i=0;i<CamNo;i++){
		if(i<pipeline1.m_allFoundEvents.back().size() && i<pipeline2.m_allFoundEvents.back().size() ){		
			CCDEventList& eventsOn1 = pipeline1.m_allFoundEvents.back()[i];
         CCDEventList& eventsOn2 = pipeline2.m_allFoundEvents.back()[i];
			
			Assert( eventsOn1.size()==eventsOn2.size(),"Number of coicyding events must be equal" );
			for(int evt=0;evt<eventsOn1.size();evt++){
				if( CCD_Analyser::CheckIfInternalTrigger( eventsOn1[evt], eventsOn2[evt] )){
					// if criteria for internal trigger are satisifed add to both lists :
					pipeline1.m_TriggerEventsList.push_back( eventsOn1[evt] );
					pipeline2.m_TriggerEventsList.push_back( eventsOn2[evt] );
				}
			}
		}
	}

	return (pipeline1.m_TriggerEventsList.size());
}


void CCDPipeline::WriteEventDetails( CCDEventList& eventList )
{
   LONG_T pos;
   register int TotalEventNo=0;


   CCDEventList::iterator pEvt;
	int i=0;
   for(pEvt=eventList.begin();pEvt!=eventList.end();pEvt++){
		if( CccdReport::IsGoodToSave( *pEvt ) ){
         pEvt->WriteEventDetails( this, i );
      }
		i++;
   }
}


int CCDPipeline::SaveEventsToFITS( CCDEventList& eventList, BOOL_T bFirstDump, 
											  BOOL_T bAll, BOOL_T bAutoDetermineIfFirst,
											  BOOL_T bAutoCheckRange, BOOL_T bSN )
{
	if(gCCDParams.m_bSaveEventDescOnly)
		return 0;

	LONG_T pos;
	register int TotalEventNo=0;
	mystring szEventDir;
	int EventNo=0;
	mystring szMainEvtDir = "Events";
	if( bSN ){
		szMainEvtDir = "SN_Events";
	}


	CCDEventList::iterator pEvt;	
	for(pEvt=eventList.begin();pEvt!=eventList.end();pEvt++){
		if( bAutoCheckRange ){
			if( m_DayFrameCounter - pEvt->m_DayFrameIndex > gCCDParams.m_nSaveFramesBeforeAndAfter ){
				continue;
			}
		}
	

		if( bAll || CccdReport::IsGoodToSave( *pEvt ) ){
			if( bAutoDetermineIfFirst ){
				if( pEvt->m_LastDumpedFrame<=0 ){
					bFirstDump = TRUE;
				}
			}
			int DumpStartPos = pEvt->m_LastDumpedFrame+1;
			if(bFirstDump)
				DumpStartPos = m_FrameCounter-gCCDParams.m_ConfirmEventsOnNextNFrames-gCCDParams.m_nSaveFramesBeforeAndAfter;
			int firstDumpPos=-1;

			int nSaved=0;
			for(register int ToDump=DumpStartPos;ToDump<=m_FrameCounter;ToDump++){
				LONG_T posToDump;
				cCCD* pFrameToDump = FindFrame( ToDump , posToDump );	
				if(pFrameToDump){
					cCCD& frameToDump = (*pFrameToDump);				
					frameToDump[pEvt->m_CameraIdx].WriteEventToFITSFile( *pEvt, ToDump, 
																						  pEvt->EvtIdx, this, 
 																						  szEventDir,
																						  szMainEvtDir.c_str() );
					if( firstDumpPos<0 ){
						firstDumpPos = ToDump;
					}
					nSaved++;
				}
			}
			if( gCCDParams.m_bSaveAverageParts ){
				if( firstDumpPos<0 )
	            firstDumpPos = DumpStartPos;
				_TRACE_PRINTF_3("Saving averages , event dir = %s, saved_count=%d\n",szEventDir.c_str(),nSaved);fflush(0);
				if( strlen( szEventDir.c_str() ) ){
					CCDMatrix::SaveEventFrameSum( *pEvt, szEventDir.c_str(), 
															 firstDumpPos,GetFrameCounter(), 
															 pEvt->EvtIdx, this );
				}else{
					LogError("ERROR : Saving averages, event dir not defined !!!!\n");
				}
			}
			EventNo++;
		}
	}
	TotalEventNo += EventNo;
	return TotalEventNo;
}


int CCDPipeline::SaveEventsToFITS( int lastEventFrameNo )
{
	if(gCCDParams.m_bSaveEventDescOnly)
		return 0;

	mystring szEventDir;
	LONG_T pos;
	register int TotalEventNo=0;
	for(register int f=lastEventFrameNo;f<=m_FrameCounter;f++){
		cCCD* pFrame = FindFrame( f , pos );
		Assert(pFrame!=NULL,"Frame no found SaveEventsToFITS-1");
		cCCD& frame = (*pFrame);	
		for(register int m=0;m<frame.GetCount();m++){
			CCDEventList::iterator pEvt;
			CCDEventList& eventList = frame[m].GetFoundEvents();
			int EventNo=0;
			int nSaved=0;
			for(pEvt=eventList.begin();pEvt!=eventList.end();pEvt++){
				int startPos = -1;								
				if( CccdReport::IsGoodToSave( *pEvt ) ){
					for(register int ToDump=pEvt->m_LastDumpedFrame+1;ToDump<=m_FrameCounter;ToDump++){
						LONG_T posToDump;
						cCCD* pFrameToDump = FindFrame( ToDump , posToDump );
						if(pFrameToDump){
							cCCD& frameToDump = (*pFrameToDump);
							frameToDump[pEvt->m_CameraIdx].WriteEventToFITSFile( *pEvt, ToDump, pEvt->EvtIdx, this, szEventDir );
							if( startPos<0 )
								startPos = ToDump;
						}
						nSaved++;
					}				
				}
				if( gCCDParams.m_bSaveAverageParts ){
					if( startPos<0 )
						startPos = pEvt->m_LastDumpedFrame+1;
					_TRACE_PRINTF_3("Saving averages , event dir = %s, saved_count=%d\n",szEventDir.c_str(),nSaved);fflush(0);
					if( strlen( szEventDir.c_str() ) ){
						CCDMatrix::SaveEventFrameSum( *pEvt, szEventDir.c_str(), 
															   startPos,GetFrameCounter(), pEvt->EvtIdx, this );
					}else{
						LogError("ERROR : Saving averages, event dir not defined !!!!\n");
					}
				}
				EventNo++;					
			}
			TotalEventNo += EventNo;
		}
	}
	return TotalEventNo;
}

void CCDPipeline::DumpAllEvents()
{
	if(gCCDParams.m_bDumpAllEvents){
		CCDEventList& allEvents = CCD_Analyser::GetAllEvents();
		allEvents.DumpAllEvents();
	}
}

CCDPipelineIterator::CCDPipelineIterator(CCDPipeline* pPipeline)
: m_pPipeline(pPipeline)
{
	m_Index = m_pPipeline->begin();	
	m_pFrame = m_pPipeline->GetFrame( m_Index );
}


cCCD* CCDPipelineIterator::operator++(int)
{
	m_Index++;
   if(m_Index==m_pPipeline->GetCount())
   	m_Index = 0;	
	m_pFrame = m_pPipeline->GetFrame( m_Index );
	return m_pFrame;
}

cCCD* CCDPipelineIterator::begin()
{
	LONG_T idx = m_pPipeline->begin();
	return m_pPipeline->GetFrame( idx );
}

cCCD* CCDPipelineIterator::end()
{
	LONG_T idx = m_pPipeline->end();
   return m_pPipeline->GetFrame( idx );
}


void CCDPipeline::GetAllMatrixPtrs( LONG_T index, CCDMatrix** pMatrixPtrTab )
{
	GetAllMatrixPtrsInt( index, pMatrixPtrTab );	
}


// it returns frames from frames at position pos up to current frame
// - ADDS ALSO CURRENT FRAME AND START FRAME !
// so that we have all frames taken >= frame at POS up to CURRENT
void CCDPipeline::GetNextMatrixPtrChronological( LONG_T index, LONG_T pos, CPixelAnalyseIn& in )
{
	int i=pos;
	int cnt=0;	
	BOOL_T bAddNext = TRUE;

	if(i<0)
		i = m_nPipelineSize-1;		
	while(bAddNext){
		// here used to store NEXT frames (not PREV)
		in.PrevMatrixPtr[cnt] = GetFrame(i)->GetMatrix( index );				

		in.PrevFramesShiftTab[cnt].frameDX = cnt*gCCDParams.m_FrameDX;
		in.PrevFramesShiftTab[cnt].frameDY = cnt*gCCDParams.m_FrameDY;



		// this 0.82 - why it is required to have exact rotation (dAlfa) as taken from image ????
		// 0.82 = sin( 56 deg ) - 56 deg = szerkosc geograficzna Warszawy	
		BOOL_T bIgnoreTimeError = gCCDParams.m_bIgnoreMissingTime;

		double t1 = 0;
		double t2 = 0;

		if( !gCCDParams.m_bReadFirstLevelInfoFromLog || !gCCDParams.GetMC() || gCCDParams.m_bFromLogFile_WithFrames ){
			t1 = ((CCDMatrix*)in.PrevMatrixPtr[0])->getObsTime( bIgnoreTimeError );
			t2 = ((CCDMatrix*)in.PrevMatrixPtr[cnt])->getObsTime( bIgnoreTimeError );
		}

		if(t1>0 && t2>0)
			in.PrevFramesTime[cnt] = ( t2 - t1  );
		else
			in.PrevFramesTime[cnt] = 0;

		if(i==tail)
			bAddNext=FALSE;

		i--;
		cnt++;
		if(i<0)
			i = m_nPipelineSize-1;		
	}
	in.PrevMatrixPtrCnt = cnt; 
	
	/*CPixelAnalyseIn in;
   register int allCount = GetAllMatrixPtrsChronologicalInt( index, in, TRUE );
	register int nFramesToAdd = (pos-tail);
	if(nFramesToAdd<0)
		nFramesToAdd = (pos+m_nPipelineSize-tail);

	cnt = 0;
	for(register int i=0;i<nFramesToAdd;i++){
		NextMatrixPtr[cnt] = (CCDMatrix*)(in.PrevMatrixPtr[i]);
		cnt++;
	}*/
}

int CCDPipeline::GetNPreviousFrames( LONG_T index, LONG_T nNextToTake, CCDMatrix** NextMatrixPtr )
{
	CPixelAnalyseIn in;
   int allCount = GetAllMatrixPtrsChronologicalInt( index, in, TRUE );
	Assert(allCount>=nNextToTake,"Not enough frames in pipeline");
	register int cnt=0;
	for(register int i=nNextToTake-1;i>=0;i--){
		NextMatrixPtr[cnt] = (CCDMatrix*)(in.PrevMatrixPtr[i]);
		cnt++;
	}
	return cnt;
}


int CCDPipeline::GetAllMatrixPtrsInt( LONG_T index, CCDMatrix** pMatrixPtrTab )
{
	CCDPipelineIterator i(this);
	
	int cnt=0;
	for(;;i++){
		pMatrixPtrTab[cnt] = i->GetMatrix( index );
		cnt++;
		if(i.curr())
			break;
	}	
	return cnt;
}


// in.PrevMatrixPtr[0] - previous (or current)
// in.PrevMatrixPtr[6] - oldest
int CCDPipeline::GetAllMatrixPtrsChronologicalInt( LONG_T index, 
						CPixelAnalyseIn& in, 
						BOOL_T bAddCurrent, int frames_to_add /* =-1 */ )
{
	if(frames_to_add<0){
		frames_to_add = m_nPipelineSize;
		if(!bAddCurrent)
			frames_to_add--;
	}
	if(frames_to_add>m_Count)
		frames_to_add = m_Count;


	// at tail we have newest frame : ...tail....OLDEST, or ...OLDEST tail...
	// tail - newest frame
	// head - oldest frame 

	// starting from oldest frame :
	// int start_pos=head;

	int start_pos = tail;
	if(!bAddCurrent)
		start_pos = tail+1;

	if(start_pos>=m_nPipelineSize)
		start_pos=0;

	int cnt = 0;
   int i=start_pos;

	
	CCDMatrix& currentFrame = GetCurrent()[0];

	// currentFrame.CalcTotalShift( m_FrameCounter, gCCDParams.m_FrameDX, gCCDParams.m_FrameDY);
	// long shiftCurrentDX = my_round( currentFrame.m_totalShiftX );
	// long shiftCurrentDY = my_round( currentFrame.m_totalShiftY );

	// filling in chronological order 
   while(cnt<frames_to_add && cnt<m_Count){
	   in.PrevMatrixPtr[cnt] = GetFrame(i)->GetMatrix( index );
		in.PrevLaplacePtr[cnt] = ((CCDMatrix*)in.PrevMatrixPtr[cnt])->get_frame_laplace_fast();

		long nSteps = cnt;
		long FrameIndex = m_FrameCounter-cnt;
		if( !bAddCurrent ){
			nSteps = (cnt+1);
			FrameIndex--;
		}

		in.PrevFramesShiftTab[cnt].frameDX = -nSteps*gCCDParams.m_FrameDX;
		in.PrevFramesShiftTab[cnt].frameDY = -nSteps*gCCDParams.m_FrameDY;



		// this 0.82 - why it is required to have exact rotation (dAlfa) as taken from image ????
		// 0.82 = sin( 56 deg ) - 56 deg = szerkosc geograficzna Warszawy	
		BOOL_T bIgnoreTimeError = gCCDParams.m_bIgnoreMissingTime;
		
		double t1 = ((CCDMatrix*)in.PrevMatrixPtr[0])->getObsTime( bIgnoreTimeError );
		double t2 = ((CCDMatrix*)in.PrevMatrixPtr[cnt])->getObsTime( bIgnoreTimeError);

		if(t1>0 && t2>0)
			in.PrevFramesTime[cnt] = ( t1 - t2 );
		else
			in.PrevFramesTime[cnt] = 0;


		// new version of shift calculation - calculates what was the real
      // shift :
		// (in.PrevMatrixPtr[cnt])->CalcTotalShift( FrameIndex, gCCDParams.m_FrameDX, gCCDParams.m_FrameDY );		
		// in.PrevFramesShiftTab[cnt].frameDX = ( my_round((in.PrevMatrixPtr[cnt])->m_totalShiftX) - my_round( currentFrame.m_totalShiftX ) );
		// in.PrevFramesShiftTab[cnt].frameDY = ( my_round((in.PrevMatrixPtr[cnt])->m_totalShiftY) - my_round( currentFrame.m_totalShiftY ) );
		

   	cnt++;
		i++;
      if(i>=m_nPipelineSize)
         i = 0;
   }
	in.PrevMatrixPtrCnt = cnt;
   return cnt;
}


int CCDPipeline::GetAllMatrixPtrsChronological( LONG_T index, Table2D<ELEM_TYPE>** pMatrixPtrTab,BOOL_T bAddCurrent )
{
	int frames_to_add = m_nPipelineSize;
   if(!bAddCurrent)
      frames_to_add--;
	CPixelAnalyseIn in;
	long backFrames = GetAllMatrixPtrsChronologicalInt( index, in, bAddCurrent );
	for(register long i=0;i<backFrames;i++){
		pMatrixPtrTab[i] = in.PrevMatrixPtr[i];
	}
   return backFrames;
}


void CCDPipeline::get_out_file_name( mystring& szOutName, const char* szSubDir,
                                   const char* szBaseName  )
{
	szOutName = "";
	szOutName << gCCDParams.GetOutputDir() << "/" << szSubDir << "/" << szBaseName;
}

void CCDPipeline::GetOutFileName( mystring& szOutName, const char* szSubDir, 
											 const char* szBaseName, CCDPipeline* pPipeline, 
											 int ccd_index, BOOL_T bTXT /*=TRUE*/ )
{
   szOutName = "";
   szOutName << gCCDParams.GetOutputDir() << "/" << szSubDir << "/";
	if( pPipeline ){
	   if( pPipeline->GetPipelineCount()>1 )
   	   szOutName << "Cam" << pPipeline->GetPipelineIndex() << "/";
	}		
   szOutName << szBaseName;
	if(ccd_index>=0)
		szOutName << "_ccd" << ccd_index;
	if(bTXT){
		szOutName << ".txt";
	}
}


BOOL_T CCDPipeline::CheckIfReGenEvent( int x, int y )
{
	CPointList::iterator i;
	for(i=m_ReGenEventsList.begin();i!=m_ReGenEventsList.end();i++){
		if( CCDMatrix::IsGenEventIdentified( (int)i->x, (int)i->y, x, y )){
			//if( fabs(i->x - x)<=1 && fabs(i->y - y)<=1 ){
			return TRUE;
		}
	}
	return FALSE;
}

BOOL_T CCDPipeline::UpdateCurrentFITSFile( const char* szFileName, CCDAsasTransform* pAsasTransform )
{
	CCDAsasTransform* pAsasTransformTmp=pAsasTransform;

// NEW 2008-04-23, in case of problem comment out 3 lines below :
	if( !pAsasTransformTmp ){
		printf("INFO : undefined pointer pAsasTransform using member m_pAsasTransform instead\n");
		pAsasTransformTmp = m_pAsasTransform;
	}

	if( pAsasTransformTmp ){
		if( pAsasTransformTmp->WasAstroCanceled() || pAsasTransformTmp->AstroBreakForced() || (m_WorkingMode.m_LastCommand==eDAQReq_StopAnalysis) ){
			// skip when previous astrometry was canceled or last command was StopAnalysis
			// next succesfull/error astrometry will update all flags and 
			// astrometry will be saved to FITS/DB 
			_STDOUT_TRACE_1("INFO : previous astrometry was cancelled ( e.g. StopAnalysis ), FITS/DB not updated\n");
			return FALSE;
		}
	}else{
		printf("INFO : pAsasTransformTmp is null check for StopAnalysis command skipped\n");
	}

	CCDMatrix tmp( m_SizeX, m_SizeY );
	if( !tmp.ReadFITSFile( szFileName )){
		printf("could not read frame for update of Astrometry : %s\n",szFileName);
		return FALSE;
	}
	AddCelestialCoo( tmp, pAsasTransform );

//ms 20070910 - name of file was scan - not correct !!!
//	const char* szOutName = tmp.getKeyValue( FILENAME );
	const char* szOutName = szFileName;
	if( szOutName && szOutName[0] ){
		tmp.WriteToFITSFile( szOutName );	
		printf("Saved astrometry to file : %s\n",szOutName);

		return TRUE;
	}
	printf("could not update astrometry in file : %s\n",szFileName);


	return FALSE;
}

BOOL_T CCDPipeline::UpdateCurrentFITSFile_PhotoOnly( const char* szFileName )
{
	CCDMatrix tmp( m_SizeX, m_SizeY );
	if( !tmp.ReadFITSFile( szFileName )){
		printf("could not read frame for update of Astrometry : %s\n",szFileName);
		return FALSE;
	}

	CSafeKeyTab& keyTab = tmp.GetKeyTab();
	keyTab.Set( NSTARS, m_nPhotometryStarsCount );

	tmp.WriteToFITSFile( szFileName );	
	printf("Saved astrometry to file : %s\n",szFileName);

	return TRUE;
}


void CCDPipeline::AddAstrometry( CSafeKeyTab& keyTab, CCDAsasTransform* pTransform, BOOL_T bSaveFlip, int flip )
{
	if( pTransform ){
		// saving astrometry to FITS :
		keyTab.Set( ASTROOK , "1" );
		keyTab.Set( POSANGLE, pTransform->fi );
		keyTab.Set( AST_ORD,  pTransform->order );
		keyTab.Set( PIXSCALE, pTransform->pixscale );
		keyTab.Set( AST_UTTIME, (int)pTransform->m_TransformUtTime );
		keyTab.Set( AST_ERR, pTransform->ast_err );
		if( bSaveFlip ){
			keyTab.Set( FLIP, (int)gCCDParams.m_eReverseForTransform );
		}
		if( flip >= 0 ){
			keyTab.Set( FLIP, flip );
		}
		
		for(int i=0;i<14;i++){
			mystring szParX;
			szParX << PARX << i;
			keyTab.Set( szParX.c_str(), pTransform->px[i] );
		} 
		for(int i=0;i<14;i++){
			mystring szParY;
			szParY << PARY << i;
			keyTab.Set( szParY.c_str(), pTransform->py[i] );
		} 

		char szTmp[64];
		sprintf(szTmp,"%.8f",pTransform->ra);
		keyTab.Set( RA_OBS, szTmp );

		sprintf(szTmp,"%.8f",pTransform->dec);
		keyTab.Set( DEC_OBS, szTmp );
	}
	
}

void CCDPipeline::AddCelestialCoo( CSafeKeyTab& keyTab, double ra, double dec, 
											  double alt, double azim, double ut_start, 
											  double ha, eObservationMode_T obsMode,
											  CCDAsasTransform* pAsasTransform )
{
	if( !pAsasTransform ){
		pAsasTransform = m_pAsasTransform;
	}

	_TRACE_PRINTF_3("in CCDPipeline::AddCelestialCoo ...");fflush(0);
	double ra_h = AstroAngle::rad2hours( ra );	
	keyTab.Set( fitsHeaderKeyDefTab[eRA_OBS].szKeyName, ra_h, AstroAngle::toString( ra, ANGLE_RA_TYPE).c_str() );	

	double dec_deg = AstroAngle::rad2deg( dec );
	keyTab.Set( fitsHeaderKeyDefTab[eDEC_OBS].szKeyName, dec_deg, AstroAngle::toString( dec, ANGLE_DEC_TYPE ).c_str() );	

	double ha_h = AstroAngle::rad2hours( ha );
	keyTab.Set( fitsHeaderKeyDefTab[eHA_OBS].szKeyName, ha_h, AstroAngle::toString( ha, ANGLE_HA_TYPE).c_str() );	

	double azim_deg = AstroAngle::rad2deg( azim );
	keyTab.Set( fitsHeaderKeyDefTab[eAZIM_OBS].szKeyName, azim_deg, AstroAngle::toString( azim, ANGLE_AZIM1_TYPE).c_str() );	

	double alt_deg = AstroAngle::rad2deg( alt );
	keyTab.Set( fitsHeaderKeyDefTab[eALT_OBS].szKeyName, alt_deg, AstroAngle::toString( alt, ANGLE_ALT_TYPE ).c_str() );	

	double zenith_d = PI_2_VALUE - alt;
	double zenit_deg = AstroAngle::rad2deg( zenith_d );
	keyTab.Set( fitsHeaderKeyDefTab[eZENITH_D].szKeyName, zenit_deg, AstroAngle::toString( zenith_d, ANGLE_DEG_TYPE ).c_str() );

	keyTab.Set( fitsHeaderKeyDefTab[eOBSMODE].szKeyName, (int)obsMode );


	// Sideral time :
	char szTmp[128];
	double sid = AstroCCD::getSiderealTimeLocal(	(time_t)ut_start, m_PipelineCfg.m_GeoLongitude );
	sprintf(szTmp,"%.8f",sid);
	keyTab.Set( fitsHeaderKeyDefTab[eST].szKeyName, szTmp );

	// julian day :
	double jd = AstroCCD::getJulianDay( (time_t)ut_start );
	sprintf(szTmp,"%.8f",jd);
	keyTab.Set( fitsHeaderKeyDefTab[eJD].szKeyName, szTmp );

	// currently HJD=JD in FITS :
	double hjd = AstroCCD::getHJD( (time_t)ut_start, ra, dec );
	sprintf(szTmp,"%.8f",hjd);
	keyTab.Set( fitsHeaderKeyDefTab[eHJD].szKeyName, szTmp );


	customKeyMutex.Lock();
	if(m_CustomFITSKeys.GetCount()){
		for(int i=0;i<m_CustomFITSKeys.GetCount();i++){
			CEnvVar& key = m_CustomFITSKeys[i];
			keyTab.Set( key.szName.c_str(), key.szValue.c_str(), key.szComment.c_str() );
		}
	}
	customKeyMutex.UnLock();

	// interanal camera index - index in PipelineList table - important only for me :
	keyTab.Set( CAMIIDX, m_PipelineIndex );

	// keyTab.Add( DIMAGE, GetDayFrameCounter() );		


	if( pAsasTransform ){
		if( pAsasTransform->m_bTransformOK ){
			// saving astrometry to FITS :
			keyTab.Set( ASTROOK , "1" );
			keyTab.Set( POSANGLE, pAsasTransform->fi );
			keyTab.Set( AST_ORD,  pAsasTransform->order );
			keyTab.Set( PIXSCALE, pAsasTransform->pixscale );
			keyTab.Set( AST_UTTIME, (int)pAsasTransform->m_TransformUtTime );
			keyTab.Set( AST_ERR, pAsasTransform->ast_err );
			keyTab.Set( FLIP, (int)m_PipelineCfg.m_eReverseForTransform );
		
			for(int i=0;i<14;i++){
				mystring szParX;
				szParX << PARX << i;
				keyTab.Set( szParX.c_str(), pAsasTransform->px[i] );
			} 
			for(int i=0;i<14;i++){
				mystring szParY;
				szParY << PARY << i;
				keyTab.Set( szParY.c_str(), pAsasTransform->py[i] );
			} 
		}else{
			keyTab.Set( ASTROOK , "0" );
		}
	}

	if( !gCCDParams.GetMC() && gCCDParams.m_bPISysManagerON ){
		keyTab.Set( MOUNTRA, AstroAngle::rad2hours( m_WorkingMode.m_MountRA_InRad ) );
		keyTab.Set( MOUNTDEC, AstroAngle::rad2deg( m_WorkingMode.m_MountDec_InRad ) );
		keyTab.Set( MOUNTAZIM, AstroAngle::rad2deg( m_WorkingMode.m_MountAzim_InRad ) );
		keyTab.Set( MOUNTALT, AstroAngle::rad2deg( m_WorkingMode.m_MountAlt_InRad ) );
		keyTab.Set( MOUNTTRK, (int)m_WorkingMode.m_MountObsMode );

		char tempstring[128];
		struct tm* newtime = localtime( &(m_WorkingMode.m_MountInfoTimeStamp) );
 		sprintf(tempstring,"%.2u:%.2u:%.2u",newtime->tm_hour,newtime->tm_min,newtime->tm_sec);
		keyTab.Set( MOUNTTM, tempstring );

		mystring szTmpDTM;
		szTmpDTM << m_WorkingMode.m_MountInfoTimeStamp;
		keyTab.Set( MOUNTDTM , szTmpDTM.c_str() );
	
		double ha;
		CCDProcState* pProcState = &(m_CameraStates[0]);
	   pProcState->m_pAstroForms->calculateHourAngle( m_WorkingMode.m_MountRA_InRad, m_WorkingMode.m_MountDec_InRad, ha );
		keyTab.Set( MOUNTHA, AstroAngle::rad2hours( ha ) );

		if( m_bMountMove ){
			keyTab.Set( "MOUNTMV" , "1" );
		}
	}

//	mystring szDomeStatus = CAsasInterface::GetDomeStatusDesc( m_WorkingMode.m_eDomeStatus );
//	keyTab.Set( DOME, szDomeStatus.c_str() );

	keyTab.Set( NSTARS, m_nPhotometryStarsCount );

	printf("OK\n");fflush(0);
}


void CCDPipeline::AddCelestialCoo( CCDMatrix& matrix, CCDAsasTransform* pAsasTransform )
{
	if( m_CameraStates.size()>0 ){
		time_t ut_time = (time_t)matrix.getObsTime();
		CCDProcState* pProcState = &((m_CameraStates)[0]);
		eObservationMode_T obsMode=eEarthMovingMode;
		double ra=0,dec=0,azim=0,alt=0,ha=0;
   	// pProcState->m_pAstroForms->GetObsCoo( ra, dec, azim, alt, obsMode );
		pProcState->GetObsCoo( ut_time, ra, dec, alt, azim, ha, obsMode );
		printf("CCDPipeline::AddCelestialCoo %d (%.2f,%.2f) (%.2f,%.2f) ha=%.2f %d\n",
					ut_time, AstroAngle::rad2deg(ra), AstroAngle::rad2deg(dec), 
					AstroAngle::rad2deg(alt), AstroAngle::rad2deg(azim), 
					AstroAngle::rad2deg(ha), obsMode );
   	AddCelestialCoo( matrix.GetKeyTab(), ra, dec, alt, azim, ut_time, ha, obsMode, pAsasTransform );
	}
}

/*void CCDPipeline::AddCelestialCoo( CSafeKeyTab& keyTab )
{
//	AddCelestialCoo( keyTab, m_PipelineCfg.m_RAObs, m_PipelineCfg.m_DecObs,
//						  m_PipelineCfg.m_HorAltitude, m_PipelineCfg.m_HorAzimuth, 
//						  0, 0 );
	CCDProcState* pProcState = &((m_CameraStates)[0]);
	eObservationMode_T obsMode;
   double ra,dec,azim,alt;
	pProcState->m_pAstroForms->GetObsCoo( ra, dec, azim, alt, obsMode );
	AddCelestialCoo( keyTab, ra, dec, alt, azim, 0, 0 );

}*/


void CCDPipeline::UpdateDarkFrame()
{
	if( m_DarkFramesList.size()>0 ){
		if(! m_pDarkFrame ){
			m_pDarkFrame = new cCCD( m_SizeX,m_SizeY,m_nCCD);
		}

		if( MyFile::DoesFileExist( m_PipelineCfg.m_szDarkFrameFile.c_str() ) ){
			if( !gCCDParams.m_bOverwriteOldDark ){
				printf("DARK already exists - not overwriten , will use old one - PLEASE VERIFY\n");	
				return;
			}
		}

		mystring szError;
		
		mystring szCmd;
		szCmd << "mv " << m_PipelineCfg.m_szDarkFrameFile.c_str() << " "
				<< m_PipelineCfg.m_szDarkFrameFile.c_str() << ".sav";
		system( szCmd.c_str() );

		if(!CCDUtil::BuildMedian( m_DarkFramesList,(*m_pDarkFrame)[0], szError )){
			printf("could not calculate dark frame (frame#=%d)\n",m_PipelineCounter);
		}else{
			printf("Dark frame calculated (frame#=%d)\n",m_PipelineCounter);
			
			BOOL_T bNoCompr=FALSE;
			if( !strstr( m_PipelineCfg.m_szDarkFrameFile.c_str(), ".fitc" ) ){
				bNoCompr = TRUE;
			}
			((*m_pDarkFrame)[0]).WriteToFITSFile( m_PipelineCfg.m_szDarkFrameFile.c_str(),
															  FALSE, NULL, bNoCompr  );
		}
		RestartPipeline( FALSE );
	}
}

void CCDPipeline::GetLastFrameList( CMyStrTable& framesList )
{
	framesList.clear();

	vector<CCDPipeline*>::iterator i;
   for(i=m_PipelineList.begin();i!=m_PipelineList.end();i++){	
		mystring szFrame;
//		if( (*i)->m_pCCDDeviceInterface ){
//			szFrame = ((*i)->m_pCCDDeviceInterface->m_szLastSavedFrame).c_str();
//		}else{
			szFrame = ((*i)->m_PipelineCfg).m_szLastReducedFrameName;
//		}
		if( strlen( szFrame.c_str() )>0 ){
			mystring szFullPath,szDir;
			MyFile::GetCWD( szDir );
			szFullPath << szDir.c_str() << "/" << szFrame;
			framesList.Add( szFullPath.c_str() );
		}
	}
}

void CCDPipeline::GetLastFrame( mystring& szFrame )
{
//	if( m_pCCDDeviceInterface ){
//		szFrame = (m_pCCDDeviceInterface->m_szLastSavedFrame).c_str();
//	}else{
		szFrame = m_PipelineCfg.m_szLastReducedFrameName;
//	}
}

void CCDPipeline::CheckIfExitOnError()
{
	if( gCCDParams.m_nExitOnToManyErrors>0 ){
		mystring szMsg;
		szMsg << "CHECK : Checking if number of errors does not exceed limit of " << gCCDParams.m_nExitOnToManyErrors;
		printf("%s\n",szMsg.c_str());
		vector<CCDPipeline*>::iterator i;
   	for(i=m_PipelineList.begin();i!=m_PipelineList.end();i++){
			if( (*i)->m_CamStateInfo.m_CamErrorCount > gCCDParams.m_nExitOnToManyErrors ){
				char szTmp[512];
				sprintf(szTmp,"Number of errors in cam%d = %d > %d, exiting daq now !\n",
					(*i)->m_PipelineIndex,(*i)->m_CamStateInfo.m_CamErrorCount,
					gCCDParams.m_nExitOnToManyErrors );
				(*i)->LogError( szTmp );
				ExecExitRequest( FALSE );
			}
		}
	}
}

BOOL_T CCDPipeline::ExecExitRequest( BOOL_T bExitFile )
{
	exit(-1);
	return TRUE;
}

const char* CCDPipeline::GetCameraName()
{
//	if( m_pCCDDeviceInterface ){
//		return m_pCCDDeviceInterface->GetName();
//	}

	m_szCameraName="";
	m_szCameraName << "cam" << m_PipelineIndex;
	return m_szCameraName.c_str();
}

eObservationMode_T CCDPipeline::GetCamObsMode()
{
	CCDProcState* pProcState = &((m_CameraStates)[0]);
	eObservationMode_T obsMode;
   double ra,dec,azim,alt;
	pProcState->m_pAstroForms->GetObsCoo( ra, dec, azim, alt, obsMode );

	return obsMode;
}

eObservationMode_T CCDPipeline::GetObsMode()
{
	vector<CCDPipeline*>::iterator i;
	eObservationMode_T _obsMode,obsMode;
	double ra,dec,azim,alt;

	int cnt=0;
	for(i=m_PipelineList.begin();i!=m_PipelineList.end();i++){
	   CCDProcState* pProcState = &(((*i)->m_CameraStates)[0]);
	   pProcState->m_pAstroForms->GetObsCoo( ra, dec, azim, alt, _obsMode );
		if(cnt){
			if( _obsMode!=obsMode ){
				printf("ERROR !!! Different obsMode in cameras !!!\n");
			}
		}
		obsMode = _obsMode;
		cnt++;
	}

	return obsMode;	
}

void CCDPipeline::GetDAQParams( double& fov, double& alert_acceptance, 
										  double& min_alt_to_observe )
{
	printf("Getting DAQ parameters fov,acc,min_alt\n");fflush(0);
	fov = gCCDParams.m_FOV;
	alert_acceptance = fov/2.00;
	min_alt_to_observe = 20.00;
}

BOOL_T CCDPipeline::CheckConsistancy()
{
	BOOL_T bRet=TRUE;
	vector<CCDPipeline*>::iterator i;
   for(i=m_PipelineList.begin();i!=m_PipelineList.end();i++){
		if( ((*i)->m_PipelineCfg).m_ObsMode != gCCDParams.m_ObsMode ){
			char szTmp[1024];
			sprintf(szTmp,"camera %d have obsmode=%s , global mode is %s, error !!!\n",
						(*i)->GetPipelineIndex(),GetObsModeDesc( ((*i)->m_PipelineCfg).m_ObsMode ),
						GetObsModeDesc( gCCDParams.m_ObsMode ) );
			printf("%s\n",szTmp);
			PrintError( (*i), (*i)->GetDayFrameCounter(), "OBS_MODE_ERR", szTmp );
			bRet=FALSE;
		}
	}
	return bRet;	
}

void CCDPipeline::DumpSumEvents( CCDPipeline* pPipeline1, CCDPipeline* pPipeline2 )
{
/*	DumpEventsOnSumFrame( pPipeline1, pPipeline2, 								 
							 	 pPipeline1->m_OldEventsOnSumedFrame,
								 pPipeline2->m_OldEventsOnSumedFrame,
								 AVER_VERIF_EVENTS );*/
	// now dump events verified against tracks :
	CCDEventList tmp_list1,tmp_list2;
	CCDEventList verif1, verif2;

	int f = pPipeline1->GetDayFrameCounter();

	(pPipeline1->m_OldEventsOnSumedFrame).GetEvents( tmp_list1, -1, PLUS_INF, TRUE );
	(pPipeline2->m_OldEventsOnSumedFrame).GetEvents( tmp_list2, -1, PLUS_INF, TRUE );
	GetFinalCoicReport( tmp_list1, tmp_list2, verif1, verif2 );

	printf("CCDPipeline::DumpSumEvents (%d,%d) -> (%d,%d)\n",tmp_list1.size(),
				tmp_list2.size(),verif1.size(),verif2.size());

	if( verif1.size()!=verif2.size() ){
		printf("ERROR in code CCDPipeline::DumpSumEvents %d!=%d\n",verif1.size(),verif2.size());
		return;
	}

	if( verif1.size()>0 && verif2.size()>0 ){
		DumpEventsOnSumFrame( pPipeline1, pPipeline2, verif1, verif2, AVER_VERIF_EVENTS, FALSE );

		mystring szListFile;
		szListFile << (pPipeline1->m_PipelineCfg).m_szFramesListFile;
		printf("PIPELINE%d saving single frames of sum events %d-%d\n",
				pPipeline1->m_PipelineIndex,verif1.size(),verif2.size());
		pPipeline1->SaveParts( szListFile.c_str(), verif1, pPipeline1->m_PipelineIndex );

		szListFile = "";
		szListFile << (pPipeline2->m_PipelineCfg).m_szFramesListFile;
		pPipeline2->SaveParts( szListFile.c_str(), verif2, pPipeline2->m_PipelineIndex );

/*		if( gCCDParams.m_bSaveEventsToDB ){
			printf("Saving aververif-events to database ...\n");fflush(0);
			CCDPipeline::InsertEventsToDB( verif1, verif2, pPipeline1, pPipeline2 );
			printf("Sum final events ( %d/%d ) saved to DB\n",verif1.size(),verif2.size());
		}*/
	}

	pPipeline1->m_OldEventsOnSumedFrame.clear();
	pPipeline2->m_OldEventsOnSumedFrame.clear();	
}

void CCDPipeline::DumpEventsOnSumFrame( CCDPipeline* pPipeline1, CCDPipeline* pPipeline2 )
{
	if( pPipeline1->m_EventsOnSumedFrame.size()>0 && pPipeline2->m_EventsOnSumedFrame.size()>0 ){
		// dumping current events on sumed frame :
		DumpEventsOnSumFrame( pPipeline1, pPipeline2, 								 
								 	 pPipeline1->m_EventsOnSumedFrame,
         	                pPipeline2->m_EventsOnSumedFrame,
									 AVER_FRAME_EVENTS, TRUE, -1 );
	}

	// now dump events verified against tracks :
	CCDEventList verif1, verif2, tmp_list1, tmp_list2;
	(pPipeline1->m_OldEventsOnSumedFrame).GetEvents( tmp_list1, -1, (pPipeline1->m_DayFrameCounter-gCCDParams.m_nNumBackFramesForTracksOnSum), TRUE );
	(pPipeline2->m_OldEventsOnSumedFrame).GetEvents( tmp_list2, -1, (pPipeline2->m_DayFrameCounter-gCCDParams.m_nNumBackFramesForTracksOnSum), TRUE );
	GetFinalCoicReport( tmp_list1, tmp_list2, verif1, verif2 );

	printf("CCDPipeline::DumpSumEvents (%d,%d) -> (%d,%d)\n",tmp_list1.size(),
				tmp_list2.size(),verif1.size(),verif2.size());

	if( verif1.size()!=verif2.size() ){
		printf("ERROR in code CCDPipeline::DumpSumEvents %d!=%d\n",verif1.size(),verif2.size());
		return;
	}

	if( verif1.size()>0 && verif2.size()>0 ){
		DumpEventsOnSumFrame( pPipeline1, pPipeline2, verif1, verif2, AVER_VERIF_EVENTS, FALSE );
		printf("PIPELINE%d dumping of sum events %d-%d\n",
				pPipeline1->m_PipelineIndex,verif1.size(),verif2.size());

		mystring szListFile;
		szListFile << (pPipeline1->m_PipelineCfg).m_szFramesListFile;
		printf("PIPELINE%d saving single frames of sum events %d-%d\n",
				pPipeline1->m_PipelineIndex,verif1.size(),verif2.size());
		pPipeline1->SaveParts( szListFile.c_str(), verif1, pPipeline1->m_PipelineIndex );

		szListFile = "";
		szListFile << (pPipeline2->m_PipelineCfg).m_szFramesListFile;
		pPipeline2->SaveParts( szListFile.c_str(), verif2, pPipeline2->m_PipelineIndex );

/*		if( gCCDParams.m_bSaveEventsToDB ){
			printf("Saving aververif-events to database ...\n");fflush(0);
			CCDPipeline::InsertEventsToDB( verif1, verif2, pPipeline1, pPipeline2 );
			printf("Sum final events ( %d/%d ) saved to DB\n",verif1.size(),verif2.size());
		}*/
	}


	// remove old events from list :
	(pPipeline1->m_OldEventsOnSumedFrame).remove_older( (pPipeline1->m_DayFrameCounter-gCCDParams.m_nNumBackFramesForTracksOnSum) );
	(pPipeline2->m_OldEventsOnSumedFrame).remove_older( (pPipeline2->m_DayFrameCounter-gCCDParams.m_nNumBackFramesForTracksOnSum) );
	if( (pPipeline1->m_OldEventsOnSumedFrame).size()>1000 ){
		printf("WARNING : to many events in buffer m_OldEventsOnSumedFrame, removing all - please verify !!!\n");
		(pPipeline1->m_OldEventsOnSumedFrame).clear();
		(pPipeline2->m_OldEventsOnSumedFrame).clear();
	}

	// to keep events for next averaged frame :
	pPipeline1->m_OldEventsOnSumedFrame += pPipeline1->m_EventsOnSumedFrame;
	pPipeline2->m_OldEventsOnSumedFrame += pPipeline2->m_EventsOnSumedFrame;

	if( gCCDParams.m_bCheckForSUPERNEW ){
		DumpSNEventsOnSumFrame( pPipeline1, pPipeline2 );
	}
}

void CCDPipeline::DumpSNEventsOnSumFrame( CCDPipeline* pPipeline1, CCDPipeline* pPipeline2 )
{
	if( pPipeline1->m_BrightenOnSumList.size()>0 && pPipeline2->m_BrightenOnSumList.size()>0 ){
		// dumping current events on sumed frame :
		DumpEventsOnSumFrame( pPipeline1, pPipeline2, 								 
								 	 pPipeline1->m_BrightenOnSumList,
         	                pPipeline2->m_BrightenOnSumList,
									 SN_SUM_COIC_LOG, TRUE, 0, "SN_SumEvents" );
	}

	// now dump events verified against tracks :
	/*CCDEventList verif1, verif2, tmp_list1, tmp_list2;
	(pPipeline1->m_OldEventsOnSumedFrame).GetEvents( tmp_list1, -1, (pPipeline1->m_DayFrameCounter-gCCDParams.m_nNumBackFramesForTracksOnSum), TRUE );
	(pPipeline2->m_OldEventsOnSumedFrame).GetEvents( tmp_list2, -1, (pPipeline2->m_DayFrameCounter-gCCDParams.m_nNumBackFramesForTracksOnSum), TRUE );
	GetFinalCoicReport( tmp_list1, tmp_list2, verif1, verif2 );

	printf("CCDPipeline::DumpSumEvents (%d,%d) -> (%d,%d)\n",tmp_list1.size(),
				tmp_list2.size(),verif1.size(),verif2.size());

	if( verif1.size()!=verif2.size() ){
		printf("ERROR in code CCDPipeline::DumpSumEvents %d!=%d\n",verif1.size(),verif2.size());
		return;
	}

	if( verif1.size()>0 && verif2.size()>0 ){
		DumpEventsOnSumFrame( pPipeline1, pPipeline2, verif1, verif2, AVER_VERIF_EVENTS, FALSE );
		printf("PIPELINE%d dumping of sum events %d-%d\n",
				pPipeline1->m_PipelineIndex,verif1.size(),verif2.size());

		mystring szListFile;
		szListFile << (pPipeline1->m_PipelineCfg).m_szFramesListFile;
		printf("PIPELINE%d saving single frames of sum events %d-%d\n",
				pPipeline1->m_PipelineIndex,verif1.size(),verif2.size());
		pPipeline1->SaveParts( szListFile.c_str(), verif1, pPipeline1->m_PipelineIndex );

		szListFile = "";
		szListFile << (pPipeline2->m_PipelineCfg).m_szFramesListFile;
		pPipeline2->SaveParts( szListFile.c_str(), verif2, pPipeline2->m_PipelineIndex );
	}


	// remove old events from list :
	(pPipeline1->m_OldEventsOnSumedFrame).remove_older( (pPipeline1->m_DayFrameCounter-gCCDParams.m_nNumBackFramesForTracksOnSum) );
	(pPipeline2->m_OldEventsOnSumedFrame).remove_older( (pPipeline2->m_DayFrameCounter-gCCDParams.m_nNumBackFramesForTracksOnSum) );
	if( (pPipeline1->m_OldEventsOnSumedFrame).size()>1000 ){
		printf("WARNING : to many events in buffer m_OldEventsOnSumedFrame, removing all - please verify !!!\n");
		(pPipeline1->m_OldEventsOnSumedFrame).clear();
		(pPipeline2->m_OldEventsOnSumedFrame).clear();
	}

	// to keep events for next averaged frame :
	pPipeline1->m_OldEventsOnSumedFrame += pPipeline1->m_EventsOnSumedFrame;
	pPipeline2->m_OldEventsOnSumedFrame += pPipeline2->m_EventsOnSumedFrame;
	*/
}



void CCDPipeline::DumpEventsOnOldFrame( CCDPipeline* pPipeline1, CCDPipeline* pPipeline2 )
{
	DumpEventsOnSumFrame( pPipeline1, pPipeline2, 								 
							 	 pPipeline1->m_EventsFromCompareToOld,
                         pPipeline2->m_EventsFromCompareToOld,
								 COMPARE_TO_OLD_EVENTS );
}


void CCDPipeline::DumpEventsOnSumFrame( CCDPipeline* pPipeline1, CCDPipeline* pPipeline2,
													 CCDEventList& events_on1, CCDEventList& events_on2,
													 const char* basename, BOOL_T bSaveParts, 
													 int limit, const char* szEvtSubDir/*=SumEvents*/ )
{
	char szLog1[512],szLog2[512];

	sprintf(szLog1,basename,0);
	sprintf(szLog2,basename,1);

	if( limit>0 ){
		if( events_on1.size()>limit ){
			printf("CCDDouble::DumpEventsOnSumFrame - all rejected : %d\n",events_on1.size());
			events_on1.clear();
			events_on2.clear();
		}
	}

	CCDPipeline::DumpFinalCoicReport( pPipeline1, pPipeline2,
												 events_on1, events_on2,
												 szLog1, szLog2, FALSE );

	if( strstr(basename,AVER_VERIF_EVENTS) ){
		if( gCCDParams.m_bSaveSumEventsToDB ){
			printf("CCDPipeline::DumpEventsOnSumFrame : Saving aververif-events to database ...\n");fflush(0);
			CCDPipeline::InsertEventsToDB( events_on1, events_on2, pPipeline1, pPipeline2, TRUE  );
			printf("CCDPipeline::DumpEventsOnSumFrame : Sum final events ( %d/%d ) saved to DB\n",events_on1.size(),events_on2.size());
		}
	}


	if( events_on1.size()>0 && gCCDParams.m_bSaveFramesWithEventsOnSum && bSaveParts ){
		Table2D<BIG_ELEM_TYPE>* pCurrSumFrame1 = pPipeline1->GetCurrSum();
	   Table2D<BIG_ELEM_TYPE>* pPrevSumFrame1 = pPipeline1->GetPrevSum();
		Table2D<BIG_ELEM_TYPE>* pCurrSumFrame2 = pPipeline2->GetCurrSum();
	   Table2D<BIG_ELEM_TYPE>* pPrevSumFrame2 = pPipeline2->GetPrevSum();

		int nPrevFrame1 = (pPipeline1->GetDayFrameCounter()-gCCDParams.m_bKeepSumOfPrevNFrames );
		int nPrevFrame2 = (pPipeline2->GetDayFrameCounter()-gCCDParams.m_bKeepSumOfPrevNFrames );

		pPipeline1->SaveParts( pPrevSumFrame1, events_on1, 
					  (pPipeline1->GetCurrent())[0],
					  0, "prev_average", nPrevFrame1, szEvtSubDir );
		pPipeline2->SaveParts( pPrevSumFrame2, events_on2, 
					  (pPipeline2->GetCurrent())[0], 1, "prev_average", 
						nPrevFrame2, szEvtSubDir  );
		
		pPipeline1->SaveParts( pCurrSumFrame1, events_on1, 
					  (pPipeline1->GetCurrent())[0], 0, "curr_average", 
						pPipeline1->GetDayFrameCounter(), szEvtSubDir );
		pPipeline2->SaveParts( pCurrSumFrame2, events_on2,
					  (pPipeline2->GetCurrent())[0], 1, "curr_average", 
						pPipeline2->GetDayFrameCounter(), szEvtSubDir );
	}
}

void CCDPipeline::SaveSNSumFrameEventsOnNext()
{
	if( m_BrightenOnSumVerifList.size()>0 ){
		printf("CCDPipeline::SaveSNSumFrameEventsOnNext saving next SumFrame SN-events ...\n");
		CCDEventList prev_frame_events;
		m_BrightenOnSumVerifList.GetEvents( prev_frame_events, m_DayFrameCounter-gCCDParams.m_bKeepSumOfPrevNFrames );

		if( prev_frame_events.size()>0 ){
			printf("CCDPipeline::SaveSNSumFrameEventsOnNext %d events to be saved\n",prev_frame_events.size());
			SaveSNSumFrameEvents( "next_average", prev_frame_events );
		}else{
			printf("CCDPipeline::SaveSNSumFrameEventsOnNext nothing to save\n");
		}

		// currently saving only on next average frame :
		// m_BrightenOnSumVerifList.clear();
		
		// clean only events not usfull for track fit :
		
	}
}


void CCDPipeline::SaveSumFrameEventsOnNext()
{
	if( m_OldEventsOnSumedFrame.size()>0 ){
		printf("CCDPipeline::SaveSumFrameEventsOnNext saving next SumFrame events ...\n");
		CCDEventList prev_frame_events;
		m_OldEventsOnSumedFrame.GetEvents( prev_frame_events, m_DayFrameCounter-gCCDParams.m_bKeepSumOfPrevNFrames );

		if( prev_frame_events.size()>0 ){
			printf("CCDPipeline::SaveSumFrameEventsOnNext %d events to be saved\n",prev_frame_events.size());
			SaveSumFrameEvents( "next_average", prev_frame_events );
		}else{
			printf("CCDPipeline::SaveSumFrameEventsOnNext nothing to save\n");
		}

		// currently saving only on next average frame :
		// m_OldEventsOnSumedFrame.clear();
		
		// clean only events not usfull for track fit :
		
	}
}

void CCDPipeline::SaveSumFrameEvents( const char* szFileName, CCDEventList& events )
{
	Table2D<BIG_ELEM_TYPE>* pCurrSumFrame1 = GetCurrSum();

	SaveParts( pCurrSumFrame1, events, (GetCurrent())[0], m_PipelineIndex, szFileName, m_DayFrameCounter );
}


void CCDPipeline::SaveSNSumFrameEvents( const char* szFileName, CCDEventList& events )
{
	Table2D<BIG_ELEM_TYPE>* pCurrSumFrame1 = GetCurrSum();

	SaveParts( pCurrSumFrame1, events, (GetCurrent())[0], m_PipelineIndex, szFileName, 
				  m_DayFrameCounter, "SN_SumEvents" );
}


// pozadnie klucze , i zeby nie obcianal 1 kolumny :
void CCDPipeline::SaveParts( Table2D<BIG_ELEM_TYPE>* pFrame, CCDEventList& evt_list,
									CCDMatrix& currFrame,
									int CamIdx, const char* szFileName,
									int saveFrameIndex, const char* szEvtSubDir/*=SumEvents*/ )
{
		char tmp[128];
		CCDMatrix part( 0 , 0 );
		part.GetKeyTab() = currFrame.GetKeyTab();

		eFITSCompressionType compr=gCCDParams.m_eCompressFITS;
		mystring szBaseName,szError;
		szBaseName << gCCDParams.GetOutputDir() << "/" << szEvtSubDir  << "/";


		// save parts :
		for(int i=0;i<evt_list.size();i++){
			CccdReport& evt = evt_list[i];
			char szSubDir[128];
		   sprintf(szSubDir,"Frame%.5d",evt.m_DayFrameIndex);

			mystring szFName,szDir,szBaseFileName;
			// szBaseFileName << szFileName << i << ".fit";			
			mystring szFormat;
			szFormat << EVENT_PART_FILE_NAME << ".fit";
			sprintf(tmp,szFormat.c_str(),evt.m_DayFrameIndex,0,evt.EvtIdx,saveFrameIndex);
			szBaseFileName << tmp;
	
			if( gCCDParams.m_eCompressFITS != eFITSComprNone)
		      szBaseFileName << "c";
			szDir << szBaseName << "/" << szSubDir << "/Cam" << CamIdx;
			szFName << szDir << "/" << szBaseFileName;

			int low_x=0,low_y=0,up_x=0,up_y=0;                                                                                                                                                                                                                                              
		   if(gCCDParams.m_bSaveEventSize>0){
		      low_x = (int)(MAX((evt.m_MaxPoint.x-gCCDParams.m_bSaveEventSize),0));
      		low_y = (int)(MAX((evt.m_MaxPoint.y-gCCDParams.m_bSaveEventSize),0));
				up_x = (int)(MIN((evt.m_MaxPoint.x+gCCDParams.m_bSaveEventSize),(pFrame->GetXSize()-1)));
				up_y = (int)(MIN((evt.m_MaxPoint.y+gCCDParams.m_bSaveEventSize),(pFrame->GetYSize()-1)));
				if( (up_y-low_y+1)%2 ){
					up_y--;					
				}
				if( (up_x-low_x+1)%2 ){
					up_x--;					
				}

				CCDUtil::GetPartFromBigElemToElem( *pFrame, part, low_x, low_y, up_x, up_y );

				CSafeKeyTab& keytab = part.GetKeyTab();

				/*char savearea[50];
			   sprintf(savearea,"%d %d %d %d",low_x,low_y,up_x,up_y);
   			part.GetKeyTab().Set( SAVEAREA, savearea );*/
		
				CCDMatrix::AddPartKeys( keytab, low_x,low_y,up_x,up_y,
												(int)evt.m_MaxPoint.x,
												(int)evt.m_MaxPoint.y,
												pFrame->GetXSize(), pFrame->GetYSize(),
												evt.m_DayFrameIndex, saveFrameIndex );
				
				mystring szList;
				szList << szDir << "/" << "list" << evt.EvtIdx;

				szBaseFileName.DumpToFile( szList.c_str() );

				part.WriteToFITSFile( szFName.c_str() );
			}
		}
	
}
									
									
/*void CCDPipeline::SetSynchroMode( BOOL_T bSynchroMode ){ 
	gCCDParams.m_bTakeInSynchroMode = bSynchroMode;

	gCCDParams.m_bAsynchroMode = !bSynchroMode;
}*/

void CCDPipeline::GetSynchroMode( BOOL_T& bSynchroMode, BOOL_T& bAsynchroMode )
{
	bSynchroMode = gCCDParams.m_bTakeInSynchroMode;
	bAsynchroMode = gCCDParams.m_bAsynchroMode;
}


// this is function for updating event list - for analysis of parts :
void CCDPipeline::UpdateEvents( CCD_Analyser* pAnal ) 
{
	CCDEventList& evtlist = GetNewEvents();
	CCDMatrix& frame = (GetCurrent())[0];

	CCDEventList::iterator i;
	for(i=evtlist.begin();i!=evtlist.end();i++){
		(i->m_Point).x += frame.m_X_On_Big;
		(i->m_Point).y += frame.m_Y_On_Big;				

		(i->m_MaxPoint).x += frame.m_X_On_Big;
      (i->m_MaxPoint).y += frame.m_Y_On_Big;
	}
	pAnal->CalcAdditionalInfoForEvents( *this, evtlist );

	CCDEventList& evtlist2 = m_allFoundEvents.back()[0];

	for(i=evtlist2.begin();i!=evtlist2.end();i++){
		(i->m_Point).x += frame.m_X_On_Big;
		(i->m_Point).y += frame.m_Y_On_Big;				

		(i->m_MaxPoint).x += frame.m_X_On_Big;
      (i->m_MaxPoint).y += frame.m_Y_On_Big;
	}
	pAnal->CalcAdditionalInfoForEvents( *this, evtlist2 );	
}

int CCDPipeline::SaveParts( const char* list_file, CCDEventList& evt_list, 
									 int CamIdx )
{
	int ret=0;
	if( MyFile::DoesFileExist( list_file ) ){
		if( evt_list.size()<1 )
			return 0;
		int start_frame = evt_list[0].m_DayFrameIndex-2*gCCDParams.m_bKeepSumOfPrevNFrames;
		int end_frame   = evt_list[0].m_DayFrameIndex+2*gCCDParams.m_bKeepSumOfPrevNFrames;

		printf("Saving single frames for sum events on frame : %d ( %d-%d )\n",
					GetDayFrameCounter(),start_frame,end_frame);

		CCDMatrix& dark = (*m_pDarkFrame)[0];
		CListFile list( list_file );
		CMyStrTable& tab = list.GetListTable();

		for(int i=0;i<tab.size();i++){
			mystring file = list[i];
		
			char szTmp[128];
			strcpy( szTmp, file.c_str() );
			int len = strlen( szTmp );
			for(int i=0;i<len;i++){
				if( szTmp[i] == '.' ){
					szTmp[i] = '\0';
					break;
				}
			}		

			char cam;
			int dt,frame;
			sscanf( szTmp, FITS_FILE_BASE_FORMAT, &cam, &dt, &frame );


			if( frame>=start_frame && frame<=end_frame ){
				CCDMatrix in_out(0,0);
				mystring szIn;
				szIn << gCCDParams.m_szSampleFramesDir << "/" << file.c_str();
				if( in_out.ReadFITSFile( szIn.c_str() ) ){
					in_out.Subtract( dark, in_out, TRUE );
		
					for(int j=0;j<evt_list.size();j++){
						CccdReport& evt = evt_list[j];
	

						if( frame>=start_frame && frame<=end_frame ){
							int low_x = (int)MAX( (evt.m_MaxPoint.x-gCCDParams.m_bSaveEventSize) , 0 );
							int low_y = (int)MAX( (evt.m_MaxPoint.y-gCCDParams.m_bSaveEventSize) , 0 );
							int up_x = (int)MIN( (evt.m_MaxPoint.x+gCCDParams.m_bSaveEventSize) , (in_out.GetXSize()-1) );
							int up_y = (int)MIN( (evt.m_MaxPoint.y+gCCDParams.m_bSaveEventSize) , (in_out.GetYSize()-1) );

							printf("getting part (%d,%d,%d,%d) of file : %s\n",low_x,low_y,up_x,up_y,file.c_str() );
						
							CCDMatrix part( 0,0 );
	
							part.GetKeyTab() = in_out.GetKeyTab();
							eFITSCompressionType compr=gCCDParams.m_eCompressFITS;
							mystring szBaseName,szError;
							szBaseName << gCCDParams.GetOutputDir() << "/SumEvents/";


							char szSubDir[128],szFrameNo[64];
							sprintf(szSubDir,"Frame%.5d/Cam%d/Event%.5d/",evt.m_DayFrameIndex,CamIdx,evt.EvtIdx);
							sprintf(szFrameNo,"%.5d",frame);

							mystring szFName,szDir,szBaseFileName;
							szDir << szBaseName << szSubDir;

							mystring szFormat;
							szFormat << EVENT_PART_FILE_NAME << "_" << szFrameNo << ".fit";
							char tmp[128];
		 					sprintf(tmp,szFormat.c_str(),evt.m_DayFrameIndex,0,evt.EvtIdx,frame);
							szBaseFileName << tmp;
	
							if( gCCDParams.m_eCompressFITS != eFITSComprNone)
				      		szBaseFileName << "c";
							szFName << szDir << "/" << szBaseFileName;

							int len_x = (up_x-low_x+1);
							int len_y = (up_y-low_y+1);
							in_out.GetImage( low_x, low_y, len_x, len_y, part );

	
							CSafeKeyTab& keytab = part.GetKeyTab();		
							CCDMatrix::AddPartKeys( keytab, low_x,low_y,up_x,up_y,
														(int)evt.m_MaxPoint.x,
														(int)evt.m_MaxPoint.y,
														in_out.GetXSize(), in_out.GetYSize(),
														evt.m_DayFrameIndex, frame );
				
							mystring szList;
							szList << szDir << "/" << "list" << evt.EvtIdx;
	
							szBaseFileName.DumpToFile( szList.c_str() );
			
							part.WriteToFITSFile( szFName.c_str() );
							ret++;
						}
					}
				}else{
					printf("ERROR : could not read file : %s\n",file.c_str());
				}
			}
		}
	}		
	return ret;
}

int CCDPipeline::GetMinMaxAver( int& min_frame, int& max_frame ){
	min_frame=0;
	max_frame=0;
	int ret=0;
	if( GetCurrent().GetCount()>0 ){
		CCDMatrix& curr_image = GetCurrent()[0];
		min_frame = GetDayFrameCounter();
		max_frame = min_frame;
		const char* szAVERF0 = curr_image.GetKeyTab().getKeyVal( "AVERF0" );

		if( szAVERF0 && szAVERF0[0] ){
			mystring szKey="AVERF1";
			const char* ptr;
			int i=1;
			mystring szMaxAver;
			while( ptr = curr_image.GetKeyTab().getKeyVal( szKey.c_str() ) ){
				szMaxAver = ptr;
				szKey="AVERF";
				i++;
				szKey << i;				
			}
	
			char c;
			int night;
			if( sscanf( szAVERF0, "k2%c_%d_%d.fit",&c,&night,&min_frame)==3 ){
				if( sscanf( szMaxAver.c_str(), "k2%c_%d_%d.fit",&c,&night,&max_frame)==3 ){
					ret=1;
				}
			}
		}
	}

	return ret;
}

int CCDPipeline::GetCameraID()
{
	mystring szName,szID;
   int idx=0,camid=2;

	if( strlen( m_szCameraID.c_str() ) ){
		return atol( m_szCameraID.c_str() );
	}
	
//	if( m_pCCDDeviceInterface ){
//		m_pCCDDeviceInterface->GetDeviceName( szName, m_szCameraID , idx );
//		return atol( m_szCameraID.c_str() );			
//	}else{
		if( GetCurrent().GetCount()>0 ){
			CCDMatrix& curr_image = GetCurrent()[0];
			const char* szCAMID = curr_image.GetKeyTab().getKeyVal( CAMID );
			if( szCAMID && szCAMID[0] ){
				m_szCameraID = szCAMID;
				return atol( m_szCameraID.c_str() );	
			}
		}
//	}

	return 2;
}

int CCDPipeline::GetNightGlobal()
{
	if( m_PipelineList.size()>0 ){
		return (m_PipelineList[0])->GetNight();
	}

	return 0;
}

int CCDPipeline::GetNight()
{
	mystring szNIGHT;
	time_t curr_time = GetCurrentFrameTime();
	get_night_date_gmt( curr_time, szNIGHT );

	return atol( szNIGHT.c_str() );
}


int CCDPipeline::GetDayFrameCounter()
{
   // return m_FrameCounter+m_PipelineCfg.m_DayFramesCounter;
   return m_DayFrameCounter;
}
                                                                                
