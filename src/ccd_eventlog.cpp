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
#include "ccd_eventlog.h"
#include "ccd_globals.h"
#include "ccd_pipeline.h"
#include "ccd_controller.h"
#include <cexcp.h>
#include <mymacros.h>


CCDEventLog::CCDEventLog( CCDPipeline* pPipeline, const char* szFileName,
								  eLogFileFormat _fmt )
: m_szFileName(szFileName),m_bFirstDump(TRUE),m_pPipeline(pPipeline),
  m_LogFormat(_fmt)
{

}
   


void CCDEventLog::SetFileName( const char* szFileName )
{
	m_szFileName = szFileName;
}


void CCDEventLog::DumpFinalEvents( CCDEventList& newEvents )
{
	CCDEventList::DumpFinalEvents( m_szFileName.c_str(), newEvents );
}

void CCDEventLog::SaveVerifiedEventsToDB( CCDEventList& newEvents, eLogFileFormat evtType )
{
	if( m_pPipeline ){
		int camid = m_pPipeline->GetCameraID();

//		for(int i=0;i<newEvents.size();i++){
//			(m_pPipeline->m_PIDB).AddEventSafe( newEvents[i], camid, evtType );
//		}
	}
}

void CCDEventLog::DumpNewEvents( CCDEventList& newEvents, const char* szIndex,
                  	 	         const char* msg, BOOL_T bToFile )
{
	int nDumped =0;
	if(msg && msg[0]){
		// passing desc :
		nDumped = CCDEventList::DumpEventReport( m_szFileName.c_str(), newEvents, m_pPipeline->GetPipelineIndex(),
															 gCCDParams.m_bStdoutOn, 
															 gCCDParams.m_bDumpNewEvents && bToFile,	
															 m_bFirstDump, FALSE, msg );
	}else{
		// using default value of description :
		nDumped = CCDEventList::DumpEventReport( m_szFileName.c_str(), newEvents, m_pPipeline->GetPipelineIndex(),
															 gCCDParams.m_bStdoutOn, 
															 gCCDParams.m_bDumpNewEvents && bToFile,	
															 m_bFirstDump, FALSE );
	}

	printf("CCDEventLog::DumpNewEvents m_LogFormat=%d\n",m_LogFormat);

	if( gCCDParams.m_bUseDB && gCCDParams.m_bSaveVerifEventsToDB && m_LogFormat==eVerif ){
		// saving verified events to database :
		printf("Saving verified events to database\n");
		SaveVerifiedEventsToDB( newEvents, eVerif );
	}

	if( gCCDParams.m_bUseDB && gCCDParams.m_bSaveAllEventsToDB && m_LogFormat==eAllEventsLog ){
		// saving verified events to database :
		printf("Saving all events to database\n");
		SaveVerifiedEventsToDB( newEvents, eAllEventsLog );
	}

	/*if(gCCDParams.m_bSaveFramesWithEvents){
		Assert(m_pPipeline!=NULL,"No pipeline object set for this CCDEventLog(%s) object",m_szFileName.c_str());
		m_pPipeline->SaveEventsToFITS( newEvents );
		int last = m_pPipeline->GetFrameIndex()-gCCDParams.m_nSaveFramesBeforeAndAfter;
		m_pPipeline->SaveEventsToFITS( last );
	}*/

	if(nDumped>0)
		m_bFirstDump = FALSE;
//	if(msg && strlen(msg))
//		_CCDLIB_PRINTF_("%s\n",msg);		

}


void CCDEventLog::Init( BOOL_T bHeader, eLogFileFormat fmt )
{
	m_bFirstDump = TRUE;

	Init( m_szFileName.c_str(), bHeader, fmt );
}

void CCDEventLog::Init( const char* fname, BOOL_T bHeader, eLogFileFormat fmt )
{
	// cleaning file if existed :
	mystring szFName;
	szFName << gCCDParams.GetOutputDir() << "/" << fname;
	MyOFile out( szFName.c_str() , "w" );
	if( bHeader ){
		if( fmt==eVerif ){
			out.Printf("%s",CccdReport::m_EventDescHeader.c_str());
		}
		if( fmt==eFinal ){
			out.Printf("%s",CccdReport::m_FinalEventHeader.c_str());
		}		
	}

	out.Open();
	out.Close();

	_TRACE_PRINTF_1("::Init Cleaned log file : %s\n",szFName.c_str());
}


void CCDEventLog::Init2()
{
	m_bFirstDump = TRUE;

	// cleaning file if existed :
	MyOFile out( m_szFileName.c_str() , "w" );
	out.Open();
	out.Close();

	_TRACE_PRINTF_1("::Init2 Cleaned log file : %s\n",m_szFileName.c_str());
}


