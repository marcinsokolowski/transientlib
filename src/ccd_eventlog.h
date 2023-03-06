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
#ifndef _CCD_EVENT_LOG_H__
#define _CCD_EVENT_LOG_H__

#include <mystring.h>
#include "ccd_report.h"

class CCDPipeline;

class CCDEventLog
{
public:
	CCDEventLog( CCDPipeline* pPipeline, const char* szFileName="", 
					 eLogFileFormat _fmt=eUnknownLogFormat );
	void SetFileName( const char* szFileName );

	// dumping and saving events :
	void DumpNewEvents( CCDEventList& newEvents, const char* szIndex, 
			 			     const char* msg, BOOL_T bToFile=TRUE );		
			 			     
	void DumpFinalEvents( CCDEventList& newEvents );
	
	// saving events to database :
	void SaveVerifiedEventsToDB( CCDEventList& newEvents, eLogFileFormat evtType=eVerif );
	


	// should be removed - after making OK with GetOutputDir !!!
	void Init( BOOL_T bHeader=FALSE, eLogFileFormat fmt=eUnknownLogFormat );
	static void Init( const char* fname, BOOL_T bHeader=FALSE, eLogFileFormat fmt=eUnknownLogFormat );


	// 
	void Init2();

	mystring m_szFileName;
	BOOL_T m_bFirstDump;
	CCDPipeline* m_pPipeline;
	eLogFileFormat m_LogFormat;	
	
};


#endif

