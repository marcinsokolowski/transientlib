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
#include "ccd_log.h"
#include <mydate.h>
#include <myfile.h>
#include "ccd_pipeline.h"

CCDLog gCCDErrLog( "%s\t%s", "Code\tDesc", "Errors", "errorlog" );

CCDLog::CCDLog( const char* szFmt, const char* szHeader,
					 const char* szSubDir, const char* szBase )
: m_szFmt(szFmt), m_szHeader(szHeader), m_szSubDir( szSubDir ), m_szBase( szBase )
{
}


CCDLog::~CCDLog()
{
	
}

void CCDLog::DumpLine( int Dummy, ... )
{
	va_list plist;
   va_start(plist,Dummy);
	
	DumpToFileBaseNew( NULL, m_szBase.c_str(), plist );

	va_end(plist);
}

void CCDLog::DumpToFile1(CCDPipeline* pPipeline , ...)
{
	va_list plist;
   va_start(plist,pPipeline);

	DumpToFileBase( pPipeline, m_szSubDir.c_str(), m_szBase.c_str(), plist );

	va_end(plist);
	
}

void CCDLog::DumpToFile2(CCDPipeline* pPipeline , int frame_index, ...)
{
	va_list plist;
   va_start(plist,frame_index);

	DumpToFileBase( pPipeline, frame_index, m_szSubDir.c_str(), 
						 m_szBase.c_str(), plist );

	va_end(plist);	
}


// same as DumpToFile2(CCDPipeline* pPipeline , int frame_index, ...)
// but different name 
void CCDLog::DumpToFile3(CCDPipeline* pPipeline , int frame_index, ...)
{
	va_list plist;
   va_start(plist,frame_index);

	DumpToFileBase( pPipeline, frame_index, m_szSubDir.c_str(), 
						 m_szBase.c_str(), plist );

	va_end(plist);	
}


void PrintError( CCDPipeline* pPipeline , int frame_index,
              	  const char* errcode, const char* errdesc )
{
	gCCDErrLog.DumpToFile3( pPipeline , frame_index, errcode, errdesc );
}

void CCDLog::DumpToFileBase( CCDPipeline* pPipeline , const char* szSubDir,
             		           const char* szBase, va_list plist )
{
	int nFrame=0;
	if( pPipeline )
		nFrame = pPipeline->GetDayFrameCounter();
	DumpToFileBase( pPipeline, nFrame, szSubDir, szBase, plist );
}

void CCDLog::DumpToFileBase( CCDPipeline* pPipeline , int frame_index,
									  const char* szSubDir,
             		           const char* szBase, va_list plist )
{
	mystring szFileName;
	CCDPipeline::GetOutFileName( szFileName, szSubDir, szBase, pPipeline , 0);

	BOOL_T bHeader=FALSE;
   if(!MyFile::DoesFileExist( szFileName.c_str() )){
      bHeader = TRUE;
   }
	
	MyOFile out( szFileName.c_str(), "a+" );
	if(bHeader){
		mystring szHeader;
		szHeader << "# DTM\t\tCCD#\tFRAME#\t" << m_szHeader << "\n";
     	out.Printf("%s",szHeader.c_str());
   }

	mystring szFmt;
	szFmt << get_gmtime_string() << "\t";
	if(pPipeline){
		szFmt << pPipeline->GetPipelineIndex() << "\t";		
	}else{
		szFmt << "0\t";
	 }
	 szFmt << frame_index << "\t" << m_szFmt;

	int len = strlen(szFmt.c_str());
	if(szFmt[len-1]!='\n')
		szFmt << "\n";

	out.VPrintf( szFmt.c_str(), plist );
	
}

void CCDLog::DumpToFileBaseNew( CCDPipeline* pPipeline ,
             		           const char* szFileName, va_list plist )
{
	BOOL_T bHeader=FALSE;
   if(!MyFile::DoesFileExist( szFileName )){
      bHeader = TRUE;
   }
	
	MyOFile out( szFileName, "a+" );
	if(bHeader){
		mystring szHeader;
		szHeader << "# DTM\t\tCCD#\tFRAME#\t" << m_szHeader << "\n";
     	out.Printf("%s",szHeader.c_str());
   }

	mystring szFmt;
	szFmt << get_gmtime_string() << "\t";
	if(pPipeline){
		szFmt << pPipeline->GetPipelineIndex() << "\t";		
	}else{
		szFmt << "0\t";
	 }
	 szFmt << 0 << "\t" << m_szFmt;

	int len = strlen(szFmt.c_str());
	if(szFmt[len-1]!='\n')
		szFmt << "\n";

	out.VPrintf( szFmt.c_str(), plist );
	
}



void CCDLog::GetOutPutFileName( mystring& szFileName, CCDPipeline* pPipeline,
										  const char* szSubDir, const char* szBase )
{
	CCDPipeline::GetOutFileName( szFileName, szSubDir, szBase, pPipeline , 0);	
}

void CCDLog::DumpToFile( CCDPipeline* pPipeline , const char* szSubDir, 
								 const char* szBase, ...)
{
	va_list plist;
   va_start(plist,szBase);

	DumpToFileBase( pPipeline ,  szSubDir, szBase, plist );

	va_end(plist);
}

void CCDLog::DumpToFile2( CCDPipeline* pPipeline , int frame_index,
								  const char* szSubDir, 
								  const char* szBase, ...)
{
	va_list plist;
   va_start(plist,szBase);

	DumpToFileBase( pPipeline , frame_index, szSubDir, szBase, plist );

	va_end(plist);
}

