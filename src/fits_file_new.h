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
#ifndef _FITS_FILE_NEW_H__
#define _FITS_FILE_NEW_H__

#include <mytypes.h>
#include <tab2D.h>
#include "fits_file_base.h"

// FITS standard defines :
// #define FLEN_CARD 81 /* length of a FITS header card */

// struct fitsfile;

class CFITSFileNew : public CFITSFileBase
{
public :
	CFITSFileNew();

	/* Pointers to input and output FITS files */
   fitsfile* m_FilePtr;
   
   /*Variables used to describe FITS files*/
   /*
   bitpix - how many bits are used to describe one pixel (def=16)
   datatype - data type code used in FITS files
   naxis - how many dimensions has the image
           (we use only 2 dimensions)
           naxes[i] - how many pixels has the image in 'i' dimension
           total_pix - total number of pixels in the image
           bscale - how we scale each pixel value when reading file (def=1.0)
           bzero - what we take as a zero level value (def=32768.0)
           */
	int status;
   int hdutype, bytepix, naxis, nkeys, datatype, any_null;
	// int bitpix;
   // long naxes[9];
   // long total_pix;
   double bscale, bzero, null_val;
   char card[FLEN_CARD];   /* Standard string lengths defined in fitsio.h - in fact 81 */
                                  

	// functions :
	BOOL_T open(const char* filename, const int iomode);
	BOOL_T load( void* data_ptr, int xSize, int ySize, int bitpix );
	BOOL_T createFITS (const char* filename, int xSize, int ySize, int bitpix );
	BOOL_T saveFITS( void* data_ptr, int xSize, int ySize  );
	void close();
	void close_file (fitsfile *file_ptr);
	int errors();
	int listHeaderKeywords(int single=1);
	int showHeaderKeyword(char *keyword);
	int modifyHeaderKeyword( char *keyword, char *new_value, char* new_comment);
	int insertHeaderKeyword(int keynum, char *keyword, char *new_value, char* new_comment);
	int renameHeaderKeyword( char *keyword, char *new_keyword);
	int setBZero();	
	int headerPI2ASAS();
	
                                  
};


#endif

