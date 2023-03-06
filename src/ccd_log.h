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
#ifndef _CCD_LOG_H__
#define _CCD_LOG_H__


#include <mystring.h>
#include <stdarg.h>

class CCDPipeline;

class CCDLog
{
public:
	static void GetOutPutFileName( mystring& szFileName, CCDPipeline* pPipeline,
											 const char* szSubDir, const char* szBase );

	CCDLog( const char* szFmt, const char* szHeader, 
			  const char* szSubDir="", const char* szBase="" );
	~CCDLog();

	void DumpToFileBase( CCDPipeline* pPipeline , const char* szSubDir, 
								const char* szBase, va_list plist );
	void DumpToFileBase( CCDPipeline* pPipeline , int frame_index,
								const char* szSubDir, 
								const char* szBase, va_list plist );
	void DumpToFile(CCDPipeline* pPipeline , const char* szSubDir, const char* szBase, ...);
	void DumpToFile2( CCDPipeline* pPipeline , int frame_index,
						   const char* szSubDir, const char* szBase, ...);
						   
	void DumpToFile1(CCDPipeline* pPipeline , ... ); 	
	void DumpToFile2( CCDPipeline* pPipeline , int frame_index, ... );

	void DumpToFile3( CCDPipeline* pPipeline , int frame_index, ... );


	void DumpLine( int Dummy, ... );
	void DumpToFileBaseNew( CCDPipeline* pPipeline ,
	                        const char* szFileName, va_list plist );

	mystring m_szFmt;
	mystring m_szHeader;	
	mystring m_szSubDir;
	mystring m_szBase;
};

void PrintError( CCDPipeline* pPipeline , int frame_index,
                 const char* errcode, const char* errdesc );
                    

extern CCDLog gCCDErrLog;

#endif

