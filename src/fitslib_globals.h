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
#ifndef _FITSLIB_GLOBALS_H__
#define _FITSLIB_GLOBALS_H__

#include "mytypes.h"
#include <vector>
#include <string>
#include "mylock.h"

// to use CCFitsio - put FALSE :
extern BOOL_T gUseMyFITSIO;
extern BOOL_T gDoAutoAllocData;

int GetBitPixSign( int argtype );
int GetBitPixSign( LONG_T argtype );
int GetBitPixSign( float argtype );
int GetBitPixSign( double argtype );


using namespace std;

extern BOOL_T m_bStaticInitialized;
extern std::vector<string> m_HeaderKeys;

class CSafeKeyTab;
// typedef struct fitsfile;

int getHeader( CSafeKeyTab& keyTab, void* ptr );


// it seems that CCfits/cfitsio are not multi-threaded libraries
// so we need to have lock to protect from paralel Read/Write of different 
// pictures :
extern CMyMutex gFITSLock;

// defines :
#define STD_COMMENT_LINE_1 "FITS (Flexible Image Transport System) format is defined in 'Astronomy"
#define STD_COMMENT_LINE_2 "and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H "

#endif
