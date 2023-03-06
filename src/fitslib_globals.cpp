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
#include "fitslib_globals.h"
// #include <fitsio.h>
#include "mysafekeytab.h"
#include "mymacros.h"
#include "mylock.h"

BOOL_T gUseMyFITSIO=TRUE;
BOOL_T gDoAutoAllocData=FALSE;

BOOL_T m_bStaticInitialized=FALSE;
std::vector<string> m_HeaderKeys;


int getHeader( CSafeKeyTab& keyTab, void* ptr )
{
	keyTab.Clear();

/*	fitsfile* fitsPointer = (fitsfile*)ptr;
	char keyval[FLEN_CARD],comment[FLEN_COMMENT],keyword[FLEN_KEYWORD];
   int status=0;
   keyTab.Clear();
   for(int i=0;i<m_HeaderKeys.size();i++){
		strcpy(keyword,m_HeaderKeys[i].c_str());
		status=0;
      ffgkey(fitsPointer,keyword,keyval,comment,&status);
      if(status==0){
			// printf("%s = %s\n",keyword,keyval);
         keyTab.Add( keyword , keyval );
		}else{
			_TRACE_PRINTF_6("error = %d, reading header info key : %s not found in FITS header\n",status,keyword);
			//ffgkey(fitsPointer,"UT-START ",keyval,comment,&status);
			//printf("status = %d\n",status);
		}
   }*/

	return keyTab.GetCount();
}


// it seems that CCfits/cfitsio are not multi-threaded libraries
// so we need to have lock to protect from paralel Read/Write of different
// pictures :
CMyMutex gFITSLock;


int GetBitPixSign( int argtype )
{
	return 1;
}

int GetBitPixSign( LONG_T argtype )
{
	return 1;
}

int GetBitPixSign( float argtype )
{
	return -1;
}

int GetBitPixSign( double argtype )
{
	return -1;
}

