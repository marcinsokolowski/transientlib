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
#include "ccd_runstat.h"
#include "ccd_globals.h"
#include "ccd_pipeline.h"
#include <mydate.h>
#include <cfg.h>
#include <cexcp.h>
#include "mcinput.h"

LONG_T CRunStat::m_RunIndex=0;
BACKGROUND_COUNT_TYPE CRunStat::m_BackgroundCalcType=COUNT_NORMAL;

CRunStat::CRunStat( CCDPipeline* pPipeline, int sizeX, int sizeY )
: CCD_Analyser( sizeX, sizeY ),
  m_nAllAnalysed(0),m_nBackground(0),m_TotalGenerated(0),m_TotalIdentified(0),
  m_pInputFile(NULL),m_ParamGenerated(0),m_ParamIdentified(0),m_ParamBackground(0),
  m_ParamAnalysed(0),m_nDifferentBackground(0), m_pPipeline(pPipeline),
  m_TotalBackground(0)
{
}

CRunStat::~CRunStat()
{
}

LONG_T CRunStat::CalcBackgroundEvents( LONG_T bacgr_count, LONG_T total_count )
{	
	LONGLONG_T ret = bacgr_count;
	switch( m_BackgroundCalcType ){
		case COUNT_DIVIDE_BY_NUMBER_OF_FRAMES :
			if(total_count>0)
				ret = (bacgr_count/total_count);
			else
				ret = 0;
			break;
	}
	return ret;
}

void CRunStat::StartNewRun( CInputMC* pInputFile, const char* runname/*=NULL*/ )
{
	m_RunIndex++;
	m_szRunName = "";
	if(runname && runname[0])
		m_szRunName << runname;
	else
		m_szRunName << get_date_time_string() << "_" << m_RunIndex;

	m_szOutPutName << gCCDParams.GetOutputDir() << "/" << m_szRunName << ".txt";
	printf("Output file is : %s\n",m_szOutPutName.c_str());	
	m_pInputFile = pInputFile;	
	
	PrintRunOutputHeader( m_pInputFile );
}

// calls AnalyseNewFrame from CCD_Analyser class, but also ensures 
// that statisctics is upgraded !
BOOL_T CRunStat::AnalyseNewFrameOpt(CCDPipeline& ccd_pipeline,BOOL_T bReport,LONG_T idx)
{
	BOOL_T bRet = CCD_Analyser::AnalyseNewFrameOpt( ccd_pipeline,bReport,idx );
	m_nAllAnalysed++;		
	m_ParamAnalysed++;	

	/*if(gCCDParams.GetMC()){
		// updating statistical information 
		CCDEventList& newEvents = ccd_pipeline.GetNewEvents();
		CCDEventList::iterator i;	
		for(i=newEvents.begin();i!=newEvents.end();i++){		
			if(i->m_bGenerated){
				m_TotalGenerated++;						
				m_ParamGenerated++;
				if(i->m_bIdentified){
					m_TotalIdentified++;
					m_ParamIdentified++;
				}			
				m_IdentStatistics.AddOrUpdateMagInfo( i->m_Magnitude , i->m_bIdentified );
			}else{
				if(i->m_bIdentified){
					m_nBackground++;
					m_ParamBackground++;
				}			
			}
		}
	}*/
	return bRet;
}


void CRunStat::GetMagsReport( mystring& szMagsReport )
{
	szMagsReport = "";
	vector<CMagStat>::iterator i;
	
	szMagsReport << " -------------    Details    -------------\n"; 
	szMagsReport << " Magnitude     Generated#       Identified#\n";
	for( i = (m_IdentStatistics.GetTab())->begin() ; i!=(m_IdentStatistics.GetTab())->end(); i++){
		szMagsReport << " " << i->m_Magnitude << "\t\t" << i->m_MagGenerated
                   << "\t\t" << i->m_MagIdentified << "\n";
	}	
}

void CRunStat::DumpRunReport()
{
	mystring szBuf,szParams,szMagsReport;
	gCCDParams.GetParams(szParams);

	// m_pPipeline->CompileEventsReport( m_TotalGenerated, m_TotalBackground, m_TotalIdentified );
	
	szBuf << " ------------- Run statistics -------------\n";
	szBuf << "# of all anaylsed frames : " << m_nAllAnalysed << "\n"; 

	
	if(gCCDParams.GetMC()){
		// m_TotalBackground = CalcBackgroundEvents(m_nBackground,m_nAllAnalysed);
		szBuf << "# of generated events : " << m_TotalGenerated << "\n";
		szBuf << "# of identified events : " << m_TotalIdentified << "\n";
		szBuf << "# of background events : " << m_TotalBackground << "\n";
	}


	GetMagsReport( szMagsReport );
	szBuf << szMagsReport;	
	
	szBuf << "\n\n\n\n ------------- Run parameters -------------\n";
	szBuf << szParams;
	PrintToFile( szBuf.c_str() );
}

void CRunStat::PrintMagSeparator()
{
	mystring szMagSep;
	szMagSep << MAG_SEPARATOR;
	AppendToFile( szMagSep.c_str() );
}

void CRunStat::PrintParamResults( CInputMC* pInputFile )
{
	CMyStrTable VaryParams;
	mystring szVaryParams;

	if(pInputFile)
	   pInputFile->GetVaryParamsValues( VaryParams, szVaryParams );
	
	mystring szRecOut;
	szRecOut << szVaryParams;
	szRecOut << PARAMS_SEPARATOR;


	CCDEventList genList;
	m_pPipeline->CompileEventsReportFromFile( m_ParamGenerated , m_ParamBackground, m_ParamIdentified, genList );

	if( gCCDParams.m_bKeepSumOfPrevNFrames>0 ){
		char szLog1[512],szLog2[512];
                                                                                
	   sprintf(szLog1,AVER_FRAME_EVENTS,0);
   	sprintf(szLog2,AVER_FRAME_EVENTS,1);
		
		CompileAndSaveFile( szLog1, szLog2, genList );
	}
	m_TotalGenerated += m_ParamGenerated;	
	m_TotalIdentified += m_ParamIdentified;
	m_TotalBackground += m_ParamBackground;

	
	/*CCDEventList::iterator genEvt;
	for(genEvt=genList.begin();genEvt!=genList.end();genEvt++){
		m_IdentStatistics.AddOrUpdateMagInfo( *genEvt );
	}

	

	vector<CMagStat>::iterator i;
	vector<CMagStat>& MagReport = (*(m_IdentStatistics.GetTab()));
	for(i=MagReport.begin();i!=MagReport.end();i++){
		szRecOut << "\n";
		szRecOut << MAGNITUDE << "=" << i->m_Magnitude << "\n";
		szRecOut << GENERATED << "=" << i->m_MagGenerated << "\n";
		szRecOut << IDENTIFIED << "=" << i->m_MagIdentified << "\n"; 

		if(MagReport.size()==1){
			// only one magintude checked - total should be equal to mag result :
			Assert( i->m_MagGenerated==m_ParamGenerated,"Number of generated wired %d!=%d",i->m_MagGenerated,(int)m_ParamGenerated);	
			Assert( i->m_MagIdentified==m_ParamIdentified,"Number of identified wired %d!=%d",i->m_MagIdentified,(int)m_ParamIdentified);
		}
	}*/

	if(genList.size()){
		szRecOut << "\n";
		szRecOut << MAGNITUDE << "=" << genList[0].m_Magnitude << "\n";
		szRecOut << GENERATED << "=" << m_ParamGenerated << "\n";
		szRecOut << IDENTIFIED << "=" << m_ParamIdentified << "\n";
	}
		
	szRecOut << TOTAL_SEPARATOR;
	szRecOut << GENERATED << "=" << m_ParamGenerated << "\n";
	szRecOut << IDENTIFIED << "=" << m_ParamIdentified << "\n"; 
	szRecOut << BACKGROUND << "=" << m_ParamBackground << "\n";
	szRecOut << RECORD_SEPARATOR;

	AppendToFile( szRecOut.c_str() );
	Reset();
}

void CRunStat::PrintTotalReport( CInputMC* pInputFile )
{
	mystring szRecOut;

	szRecOut << END_SEPARATOR;
	szRecOut << GENERATED << "=" << m_TotalGenerated << "\n";
	szRecOut << IDENTIFIED << "=" << m_TotalIdentified << "\n"; 	
	szRecOut << BACKGROUND << "=" << m_TotalBackground << "\n";  
	/*if(m_BackgroundCalcType==COUNT_NORMAL) 
		szRecOut << CalcBackgroundEvents(m_nBackground,m_nAllAnalysed) << "\n";
	else
		szRecOut << m_nDifferentBackground << "\n";*/

	AppendToFile( szRecOut.c_str() );
	printf("\n\nDumping run %s report to file : %s\n",m_szRunName.c_str(),m_szOutPutName.c_str());
}

// clears current parameter run statisctics
void CRunStat::ClearCurrentRun()
{
	m_IdentStatistics.Reset();
	m_ParamGenerated = 0;
	m_ParamIdentified = 0;
	m_ParamBackground = 0;
	m_ParamAnalysed = 0;
	GetNewEvents().clear();
}


void CRunStat::Reset()
{	
	// calculates only different background events 
	// in case all samaples are put on same frame (changing parameters only)
	// all bacground events are same on each step - so here calculate only
   // once 
	m_nDifferentBackground += CalcBackgroundEvents(m_ParamBackground,m_ParamAnalysed);

	m_IdentStatistics.Reset();
	m_ParamGenerated = 0;
	m_ParamIdentified = 0;
	m_ParamBackground = 0;
	m_ParamAnalysed = 0;
}

void CRunStat::PrintRunOutputHeader( CInputMC* pInputFile )
{
	CMyStrTable VaryParams;

	if(pInputFile)
		pInputFile->GetVaryParams( VaryParams );

	mystring szHeader;
	szHeader << "#\t\tCCD Varying Params Analysis Report\n";
	szHeader << "#\t\tStarted at : " << CMyDate::getdate() << "\n";
	szHeader << "#\t\tThe following parameters were varied\n";
	szHeader << "#" << SEP << "NAME" << SEP << "START" << SEP 
            << "END" << SEP << "STEP" << SEP << "TYPE\n";
	szHeader << "PARAM_VARIATION_DESCRIPTION\n";
	

	if(pInputFile){
		vector<CVaryParamDesc>::iterator i;
		vector<CVaryParamDesc>& VaryParamDescTab = pInputFile->GetVaryParamsDesc();
		for(i=VaryParamDescTab.begin();i!=VaryParamDescTab.end();i++){
			mystring szParam;
			szParam << i->m_ParamName << SEP << i->m_LowEnd << SEP 
         	     << i->m_UpEnd << SEP << i->m_Step << SEP 
            	  << CVaryParamDesc::GetParamTypeStr(i->m_Type) << "\n";
			szHeader << szParam;
		}
	}
	szHeader << RECORD_SEPARATOR;
	PrintToFile( szHeader.c_str(), "w" );	
}

void CRunStat::AppendToFile( const char* buff )
{
	PrintToFile( buff , "a" );
}

void CRunStat::PrintToFile( const char* buff, const char* mode )
{
	MyOFile out(m_szOutPutName.c_str(),"a");
	printf("\n\nDumping run %s report to file : %s\n",m_szRunName.c_str(),m_szOutPutName.c_str());
	printf("%s\n\n",buff);
	out.Printf("%s",buff);
	out.Close();
}


void CRunStat::CompileAndSaveFile( const char* szLog1, const char* szLog2, CCDEventList& genList )
{
	mystring szLogFile1,szLogFile2;
	szLogFile1 << gCCDParams.GetOutputDir() << "/" << szLog1;
	szLogFile2 << gCCDParams.GetOutputDir() << "/" << szLog2;

	if ( MyFile::DoesFileExist( szLogFile1.c_str() ) && MyFile::DoesFileExist( szLogFile2.c_str() ) ){
		LONG_T nGen, nBackgr, nIdent;
		m_pPipeline->CompileEventsReportFromFileNew( szLogFile1, nGen , nBackgr, nIdent,
																	genList );
		if( gCCDParams.m_bSumAllMethodsInNparTest ){
			m_ParamBackground += nBackgr;
         m_ParamIdentified += nIdent;
		}else{
			m_ParamGenerated = nGen;
			m_ParamBackground = nBackgr;
			m_ParamIdentified = nIdent;
		}
	}			
}
