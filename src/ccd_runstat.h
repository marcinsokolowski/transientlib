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
#ifndef _CCD_RUNSTAT_H__
#define _CCD_RUNSTAT_H__

#include<mytypes.h>
#include<basestructs.h>
#include "ccd_analyse.h"

enum BACKGROUND_COUNT_TYPE { COUNT_NORMAL=0, COUNT_DIVIDE_BY_NUMBER_OF_FRAMES };

class CCDPipeline;
class CInputMC;


class CRunStat : public CCD_Analyser
{
protected:
	// 
	static BACKGROUND_COUNT_TYPE m_BackgroundCalcType;
	
	// 
	mystring m_szRunName;
	mystring m_szOutPutName;
	CInputMC* m_pInputFile;

	void AppendToFile( const char* buff );	
	void PrintToFile( const char* buff, const char* mode="w" );	
	void PrintRunOutputHeader( CInputMC* pInputFile );
	LONG_T CalcBackgroundEvents( LONG_T bacgr_count, LONG_T total_count );

	void CompileAndSaveFile( const char* szLog1, const char* szLog2, CCDEventList& genList );
public :
	CCDPipeline* m_pPipeline;
	CRunStat( CCDPipeline* pPipeline, int sizeX, int sizeY  );
	~CRunStat();
	
	static void SetCalcBacground( BACKGROUND_COUNT_TYPE type ) { m_BackgroundCalcType=type; }
	
	void PrintParamResults( CInputMC* pInputFile );
	void PrintTotalReport( CInputMC* pInputFile );
	void PrintMagSeparator();
	
	//
	void StartNewRun( CInputMC* pInputFile, const char* runname=NULL );

	// calls AnalyseNewFrame from CCD_Analyser class, but also ensures 
	// that statisctics is upgraded !
	BOOL_T AnalyseNewFrameOpt(CCDPipeline& ccd_pipeline,BOOL_T bReport=FALSE,LONG_T idx=0);

	void GetMagsReport( mystring& szMagsReport );
		
	void DumpRunReport();

	void Reset();
	
	// clears current parameter run statisctics
	void ClearCurrentRun();
	
	// identification statistics of single run - by magnitude 
	CIdentStat m_IdentStatistics;

	// parameters values set :
	vector<CEnvVar> m_CCDParams;
	
	// total bacground events at given values of parameters :
	LONG_T m_nBackground;
	LONG_T m_nDifferentBackground;
	
	// total frames anaysed 
	LONG_T m_nAllAnalysed;		
	
	// total generated events 
	LONG_T m_TotalGenerated;
	
	// total identified events 
	LONG_T m_TotalIdentified;		
	
	// total bacground :
	LONG_T m_TotalBackground;
	
	
	// generated for given params set
	LONG_T m_ParamGenerated;	
	LONG_T m_ParamIdentified;	
	LONG_T m_ParamBackground; // bacground events in event set
	LONG_T m_ParamAnalysed; // all analysed frames count

	static LONG_T m_RunIndex;	
};


#endif

