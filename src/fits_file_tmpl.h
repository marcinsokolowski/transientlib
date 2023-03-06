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
#ifndef _FITS_FILE_TEMPLATE_NEW_H__
#define _FITS_FILE_TEMPLATE_NEW_H__


//----------------------------------------------------------------
// ATTENTION - PROBLEMS :
// odd sizes - not allowed - coredump
// standard header keys must not be added by addKey !!!! - coredump also !!!!
// ...
// ...
//
//----------------------------------------------------------------


#include <mytypes.h>
#include <tab2D.h>
#include <mysafekeytab.h>
#include "fits_file_new.h"


template<class ARG_TYPE>
class CFITSFileTemplate : public CFITSFileNew
{
protected :

public :
	CFITSFileTemplate();
	~CFITSFileTemplate(){};
	
	
	BOOL_T WriteToFITSFile( void* data, 
			        long low_x, long low_y, long up_x, long up_y,
			        long FrameSizeX, long FrameSizeY, mystring& szError, 
                                const char* fname );


	BOOL_T WriteToFITSFile( void* data, long SizeX, long SizeY,
	                        mystring& szError,
	                        const char* fname ); 
	BOOL_T ReadFITSFile( Table2D<ARG_TYPE>& matrix, mystring& szError,
	                      const char* fname );


	BOOL_T ReadFITSHeader( CSafeKeyTab& keys, mystring& szError, const char* fname );
};


#endif
