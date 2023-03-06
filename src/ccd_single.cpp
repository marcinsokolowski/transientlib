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
#include "ccd_image_creator.h"
#include "ccd_single.h"
#include "ccd_pipeline.h"
#include "ccd_analyse.h"
#include "ccd_globals.h"
#include "ccd_trace.h"
#include "ccd_controller.h"
#include <mymacros.h>
#include <random.h>

int CCDSingle::m_RunCount=0;

CCDSingle::CCDSingle()
: CCDSystemBase(), m_pPipeline(NULL)
{
	if(!m_pPipeline){
		long x_size = atol(gCCDParams.GetParam("CCD_SIZE_X"));
	   long y_size = atol(gCCDParams.GetParam("CCD_SIZE_Y"));

  	   m_pPipeline = new CCDPipeline(0,1,x_size,y_size);
	}
}

CCDSingle::~CCDSingle()
{
	_TRACE_PRINTF_3("~CCDSingle() ???\n");

	Clean();
}

void CCDSingle::Clean()
{
	if(m_pPipeline){
		delete m_pPipeline;
		m_pPipeline = NULL;
	}
}


void CCDSingle::Run( BOOL_T bPutImage, const char* szMagnitude, BOOL_T bSimulSuperNew )
{
	int nFrames = -1;
	m_RunCount++;

	struct timeb start_time,curr_time;
	if(bPutImage){
		CRandom::Initialize();
		m_pPipeline->InitGenObj();
	}
	
	// CImageCreator::InitParams();
	mystring szFrameFile;
	
	
	printf("WELCOME TO ANALYSIS PROGRAM - DOUBLE PIPELINE \n");
	my_printf_now("Initilizing pipelines ...");
	m_pPipeline->InitPipeline();
	my_printf_now("OK\n");

//	CCDController* pController = m_pPipeline->GetCCDInterface();
	CCD_Analyser* pAnal1 = m_pPipeline->GetAnalPtr();
	int nFramesForAstrometry=0;

	// here wait for signal from PI-MAN , exec ASAS astrometry and continue :
	if( gCCDParams.m_bDoASASAstroOnStartup>0 && !m_pPipeline->IsAstrometryOK()  ){
		SetFirstAstrometryMode();

		// wait for signal from pi-man and execute astrometry :
		// HandleRequests();
		
		// now do astrometry 
		// ExecAstrometry();
	}

	// init mode must be after so that after taking 
//	my_printf_now("Initializing daq working mode ...");
//	CCDPipeline::InitMode();
//	my_printf_now("OK\n");



/*	if(!gCCDParams.GetMC()){
		// initialize devices :
		my_printf_now("Initializing devices ...\n");
		if(pController->InitDevice( FALSE )<=0){ // do not wait for Temp
			my_printf_now("Could not initialize device, exiting ...\n");
			exit(-1);
		}
		printf_now("Waiting for temperatures ...\n");
		pController->WaitForTemp();

		my_printf_now("Devices initialized OK\n");

	}*/

	ftime(&start_time);

	mystring szDT = get_date_time_string();
	printf_now2("DAQ started at : %s\n",szDT.c_str());

	my_printf_now("Begining processing ...\n");

	long i = 0;	

	BOOL_T bContinue=TRUE;
	clock_t t1 = clock();	
	BOOL_T bStartMode=TRUE;
	int prevTn=1000000;
	int prevDayFrameIndex=-100;

	while(bContinue && m_bContinue && ( nFrames<0 || i<nFrames ) ){
		time_t anal_start=0,anal_end=0,start_time_t=get_dttm();

		// checking parameters consistancy :
		CCDPipeline::CheckConsistancy();

		time_t get_start=get_dttm();


		CCDPipeline::m_RestartDumpCounter = 0;

		if( m_pPipeline->m_bRestartOnNext ){
			// in case one pipeline needs restart - both are forced to restart here :
			my_printf_now("Forced both pipelines to be restarted in CCDDouble\n");

			// first dump events on both pipelines to have m_ConfirmedEvents list filled 
			// for DumpFinalEvents called later by RestartPipeline

			// 20070309 : when position is changed, forcing dump of all events -
			//  also not checked by tracks
			m_pPipeline->DumpNewEvents( &(m_pPipeline->m_TriggerEventsList), TRUE );
			m_pPipeline->DumpFinalEvents( &(m_pPipeline->m_FinalEventsLog),m_pPipeline->m_ConfirmedEvents, TRUE );

			if( gCCDParams.m_bCheckForSUPERNEW ){
				CCDPipeline::DumpFinalSN( m_pPipeline, NULL );
			}
			
			m_pPipeline->RestartPipeline( FALSE );
		}

		BOOL_T bNoMoreFrames=FALSE;

		if(!m_pPipeline->GetNextFrame( &m_szFName1, gCCDParams.m_bKeepAllInMemory, FALSE )){
			printf("\n\n##############################\n");
			printf("No more frames available from device : %d ...\n",m_pPipeline->GetPipelineIndex());
			printf("##############################\n\n\n");
			fflush(0);
//			CCDPipeline::StopAll();
			bNoMoreFrames = TRUE;
		}else{
			_TRACE_PRINTF_("Frame %d build from file : %s\n",(LONG_T)m_pPipeline->GetFrameIndex(),m_szFName1.c_str());fflush(0);
		}		
		
		if( bNoMoreFrames ){
			_STDOUT_TRACE_1("CCDDouble - no more frames !\n");
			m_bContinue = FALSE;
			break;						
		}

		RestartPipelinesIfNeeded();


		/* in case asdtrometry to slow and trace stars are lost ... :
		if( bStartMode ){
			if( (!m_pPipeline->m_bAstrometryExecuted || m_pPipeline->IsAstrometryOK() ) &&
				 (!m_pPipeline2->m_bAstrometryExecuted || m_pPipeline2->IsAstrometryOK() ) ){
				printf("Astrometry on first frame done on both reseting pipeline...\n");
				m_pPipeline->RestartPipeline();
				m_pPipeline2->RestartPipeline();
				bStartMode = FALSE;
			}
		}*/

		time_t get_end=get_dttm();


		if( gCCDParams.m_bDoNotAnalyse ){
			printf("No analysis mode, only collecting frames ...\n");
			continue;
		}
		

		// for new frame - in case pipeline is already filled 
		// analysis is performed here - loop over steps number :
		if(m_pPipeline->GetCount()==m_pPipeline->GetPipelineSize() ){
			if(bPutImage){
//				if(!PutSample( szMagnitude )){
//					printf("SAMPLE ERROR : could not put sample !\n"); 
//				}
				if(!PutSamples( szMagnitude, gCCDParams.m_nSamplesToPutOnFrame )){
					printf("SAMPLE ERROR : could not put sa ple !\n");
				}
	
				if( gCCDParams.m_bSaveImageWithSample ){
					mystring szOutFile;
					szOutFile << "image" << m_pPipeline->GetDayFrameCounter() << ".fitc";
					printf("saving image with sample to file : %s\n",szOutFile.c_str());
					if( !((m_pPipeline->GetCurrent())[0]).WriteToFITSFile( szOutFile.c_str() ) ){
						printf("Error in writing to fits file : %s\n",szOutFile.c_str());
					}
				}
			}

			anal_start=get_dttm();
			BOOL_T bAnalRes1 = m_pPipeline->AnalyzeNewFrame( FALSE, TRUE  );

// NEW - check number of stars and reject events if #stars < 1000 :
//			if( gCCDParams.m_bDoASASPhotAstr ){
				// only when performing astrometry - so that number of stars is known !!!
				if( m_pPipeline->m_nPhotometryStarsCount < 1000 && m_pPipeline->GetNewEvents().size()>0 ){
					printf("Number of stars on image is %d < 1000 -> all events rejected (dome probably closed)\n", m_pPipeline->m_nPhotometryStarsCount);
					m_pPipeline->GetNewEvents().clear();
				}
//			}
// NEW - check number of events after Tn ( star# indicator ) :
			if( pAnal1->backgrStat.size()>=1 ){
				int size = pAnal1->backgrStat.size();
				int currTn = (pAnal1->backgrStat).back().nTnewCut;
				int currDayFrameIndex = m_pPipeline->GetDayFrameCounter();
				printf("Number of events after Tn, prev=%d, curr=%d\n",prevTn,currTn);
				if( prevTn<5000 && currTn>5000 && (currTn-prevTn)>5000 ){
					printf("Big change of number of stars detected - probably dome opened, pipeline restart should be done and events skiped\n");
					if( (currDayFrameIndex - prevDayFrameIndex) <= 2 ){
						m_pPipeline->GetNewEvents().clear();
						m_pPipeline->m_bRestartOnNext = TRUE;
					}else{
						printf("Previous check on image=%d, current image=%d, ignored\n",prevDayFrameIndex,currDayFrameIndex);
					}
				}
				prevTn = currTn;
				prevDayFrameIndex = currDayFrameIndex;
			}
			

			CCDEventList verified_events;
			int verif_count = pAnal1->VerifySingleCameraEvents( *m_pPipeline, m_pPipeline->GetNewEvents(), verified_events );

			if( (m_pPipeline->m_allFoundEvents).size()>0 ){
				CCDEventList& verifEvents = (m_pPipeline->m_allFoundEvents).back()[0];
				verifEvents = verified_events;
			}

			// dump before coicydence :
			//m_pPipeline->DumpCurrentEvents();
			m_pPipeline->DumpNewEvents( &(m_pPipeline->m_TriggerEventsList) );
			m_pPipeline->DumpFinalEvents( &(m_pPipeline->m_FinalEventsLog),m_pPipeline->m_ConfirmedEvents, TRUE );
			int nEvents = (m_pPipeline->m_ConfirmedEvents).size();

			// m_pPipeline->ClearNewEvents();
			// m_pPipeline2->ClearNewEvents();

			// currently I used events stored in CCDPipeline.m_allFoundEvents
			// but maybe some data is incomplete - check VerifyTrack
			// if it updates correctly data in m_allFoundEvents
			// but without this information about m_LastDumpedFrame 
			// (set in SaveEvents was lost so the best way is to use
			// m_allFoundEvents when verifying events without removeing 
			// m_allFoundEvents.clear , but simply by updating then



			// if uncommented add if( !gCCDParams.m_bCCDDouble ) in ccd_analse.cpp
			if( gCCDParams.m_bCheckTracks ){
         	pAnal1->VerifyTracks( *m_pPipeline );
         }
			if( gCCDParams.m_bCheckPlaneTracks ){
				pAnal1->VerifyPlaneTracks( *m_pPipeline );
			}

			
			// after coicydence 
			// int nCount1 = m_pPipeline->DumpFoundEventsLog();

			// dumping SN events !
         if( gCCDParams.m_bCheckForSUPERNEW ){
            CCDPipeline::DumpFinalSN( m_pPipeline, NULL );
         }

			// CCDPipeline::m_nAlertsCount += nEvents;
			if( gCCDParams.m_bLogFrameStat ){
				pAnal1->LogEventRates( m_pPipeline );
			}

			// events on averaged frame :
			if( gCCDParams.m_bAnalyzeSumOfPrevNFrames ){
				// dump events :

				if( m_pPipeline->m_EventsOnSumedFrame.size()>0 ){
					if( gCCDParams.m_bRejectTracksOnSumedFrames ){
						// verify tracks on sumed frames now :
						pAnal1->VerifySumFrameTracks( m_pPipeline, (m_pPipeline->GetCurrent())[0] );
					}
					if( gCCDParams.m_bCheckNormalTracksOnSumEvents ){
						pAnal1->CheckSumEventsToNormalTracks( m_pPipeline );
					}
					// CCD_Analyser::CombineCamAnalResults( m_pPipeline->m_EventsOnSumedFrame , m_pPipeline2->m_EventsOnSumedFrame );
				}
				CCDPipeline::DumpEventsOnSumFrame( m_pPipeline, NULL );
				pAnal1->LogSumEventsRates( m_pPipeline );
			}

			if( gCCDParams.m_nCompareToOldFreqInSec>0 ){
				if( m_pPipeline->m_EventsFromCompareToOld.size()>0 ){
					CCDPipeline::DumpEventsOnOldFrame( m_pPipeline, NULL );
				}
			}
			
		
			// call this always - no not only when coicEvents>0
			// we must also dump 10 events after , also in case there
			// was no event on that frame :
			// m_pPipeline->SaveEventsToFITS();
			
			anal_end=get_dttm();						
		}else{
			printf("Frame number %d, not analysing yet ...\n",(LONG_T)m_pPipeline->GetFrameIndex());
		}

		// re-calculating maximus - using new frame
		// but not for potential hit
		time_t post_start=get_dttm();
		m_pPipeline->PostNewFrame();


		// sending info that frame is analysied :
      m_pPipeline->NewFrameAnalyseEnd();
		time_t post_end=get_dttm();

		time_t end_time_t=get_dttm();
		_TRACE_PRINTF_0("\n\nFULL ANALYSIS OF FRAME : %d TOOK %d sec \n",m_pPipeline->GetFrameIndex(),(end_time_t-start_time_t));
		_TRACE_PRINTF_0("ALGORITHM on both took : %d\n",(anal_end-anal_start));
		_TRACE_PRINTF_0("POSTING took           : %d\n",(post_end-post_start));
		_TRACE_PRINTF_0("GET_NEXT took          : %d\n\n",(get_end-get_start));
		_TRACE_PRINTF_0("Pure GET_NEXT ccd1     : %d\n",(m_pPipeline->GetPureGetNextTime()));
	}

	// 20070309 : when position is changed, forcing dump of all events -
   //  also not checked by tracks
	m_pPipeline->DumpNewEvents( &(m_pPipeline->m_TriggerEventsList), TRUE );
	m_pPipeline->DumpFinalEvents( &(m_pPipeline->m_FinalEventsLog),m_pPipeline->m_ConfirmedEvents, TRUE );

	m_pPipeline->End();

	// dumping events list :

	// commented 20040616 - already done in GetNextFrame :
	// m_pPipeline->DumpFoundEventsLog( TRUE );
   // m_pPipeline2->DumpFoundEventsLog( TRUE );

	// dumping final events log 

	// NEW - commented out on 20050407 :
	// int nEvents = CCDPipeline::DumpFinalCoicReport( m_pPipeline, NULL );

	// CCDPipeline::DumpSumEvents( m_pPipeline, m_pPipeline2 );


	if( gCCDParams.m_RunUpTo>0 && get_dttm()>gCCDParams.m_RunUpTo && !gCCDParams.GetMC() ){
		printf("Exiting due to end of analysis time scheduled up to : %d\n",gCCDParams.m_RunUpTo );
		CCDPipeline::ExecExitRequest();
	}

	// CCDPipeline::StopAll();
}



void CCDSingle::StopAnal()
{
	m_bContinue=FALSE;
	m_pPipeline->StopAnal();
}

void CCDSingle::SetFirstAstrometryMode()
{
      printf_now("Forcing astrometry on first frame \n");
      m_pPipeline->m_bForceAstroNow = TRUE;
      m_pPipeline->m_bForceSynchroMode = TRUE;
      m_pPipeline->m_nAstrometryRetry = gCCDParams.m_nASASAstrometryReTry;
}

void CCDSingle::RestartPipelinesIfNeeded()
{
	if( m_pPipeline->m_bRestartOnNext ){
		printf("CCDSingle::RestartPipelinesIfNeeded : Restarting pipeline now (%d)\n",m_pPipeline->m_bRestartOnNext);

		// OLD VERSION - commented on 20070309 :
		// first dump events on both pipelines to have m_ConfirmedEvents list filled
      // for DumpFinalEvents called later by RestartPipeline
      // m_pPipeline->DumpFoundEventsLog( TRUE );
		// CCDPipeline::DumpFinalCoicReport();

		// 20070309 :
		// when position is changed, forcing dump of all events -
      // also not checked by tracks
		m_pPipeline->DumpNewEvents( &(m_pPipeline->m_TriggerEventsList), TRUE );
		m_pPipeline->DumpFinalEvents( &(m_pPipeline->m_FinalEventsLog),m_pPipeline->m_ConfirmedEvents, TRUE );

		if( gCCDParams.m_bCheckForSUPERNEW ){
			CCDPipeline::DumpFinalSN( m_pPipeline, NULL );
		}
		

		m_pPipeline->RestartPipelineKeepCurrent( TRUE );
	}
}

BOOL_T CCDSingle::PutSample( const char* szMagnitude, BOOL_T bSimulSuperNew,
									  BOOL_T bRePut )
{
	if( bRePut ){
		// flag for putting many samples, in loop - re-put old samples
		// only for first put sample (i==0) :
		m_pPipeline->GetGenObj()->RePutPrevImages( *m_pPipeline );
	}

	m_pPipeline->GetGenObj()->AddObject( m_pPipeline->GetCurrent(), *m_pPipeline, NULL, szMagnitude, TRUE );

	if( gCCDParams.m_bKeepSumOfPrevNFrames>0 ){
		m_pPipeline->UpdateSumOfFrames( &m_szFName1 );
	}

	if(bSimulSuperNew){
		int star_x, star_y;
		double init_mag,res_mag;
		m_pPipeline->GetGenObj()->SimulateSuperNew( m_pPipeline->GetCurrent()[0],
													  m_pPipeline->GetFrameIndex(),
													  gCCDParams.m_szMinMagSUPERNEW.c_str(), 
													  gCCDParams.m_szMaxMagSUPERNEW.c_str(), 
													  gCCDParams.m_dRatioSUPERNEW, 
													  star_x, star_y, init_mag,res_mag );
	}
}

BOOL_T CCDSingle::PutSamples( const char* szMagnitude, int nSamples )
{
	int good=0;
	for(int i=0;i<nSamples;i++){
		if( PutSample( szMagnitude, FALSE, (i==0) ) ){
			good++;
		}
	}	
	return (good==nSamples);
}
