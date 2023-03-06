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
#ifndef _CCD_UTIL_H__
#define _CCD_UTIL_H__

#include <mystring.h>
#include <tab2D.h>
#include "ccd_defines.h"

class CCDMatrix;
class CMyStrTable;

class CCDUtil 
{
public :
	static BOOL_T bVerbose;

	CCDUtil();
	~CCDUtil();

	static mystring GetKeyValue( const char* fits_file, const char* key_name );

	static BOOL_T SumFrames( Table2D<BIG_ELEM_TYPE>& sumFrame, CCDMatrix& newFrame );
	
	static BOOL_T SumFrames( CMyStrTable& list, CCDMatrix& out,
									 Table2D<BIG_ELEM_TYPE>& tmp, CCDMatrix* pDark=NULL );
	static BOOL_T SaveSumFrame( CMyStrTable& list, const char* szOutName,
	                            CCDMatrix* pDark=NULL );
	
	static BOOL_T BuildMedian( CCDMatrix* pFrames , int nCount, CCDMatrix& DarkFrame );

	static BOOL_T BuildMedian( mystring szListFile, CCDMatrix& DarkFrame, 
							  mystring& szError, const char* szDark=NULL );
							  
	static BOOL_T BuildMedian( CMyStrTable& list, CCDMatrix& DarkFrame,
	                    mystring& szError, const char* szDark=NULL );		  

	static BOOL_T BuildMedian( Table2D<double>* pCCDTab, int nCount, Table2D<double>& DarkFrame );

	static BOOL_T Divide( CCDMatrix& frame,  Table2D<double>& flat, int border=0 );
	static BOOL_T Divide( CCDMatrix& frame,  Table2D<float>& flat, int border=0 );
	                    
	static BOOL_T BuildMedian( CMyStrTable& list, Table2D<double>& DarkFrame,
							  int xSize, int ySize,
           			     mystring& szError );
	                    

	static void CopyFromBigElemToElem(	Table2D<BIG_ELEM_TYPE>& bigElem, Table2D<ELEM_TYPE>& Elem );
	static void CopyFromElemToBigElem(  const Table2D<ELEM_TYPE>& Elem, Table2D<BIG_ELEM_TYPE>& bigElem );
	static void GetPartFromBigElemToElem( Table2D<BIG_ELEM_TYPE>& bigElem,
													  Table2D<ELEM_TYPE>& Elem,
													  int low_x, int low_y, int up_x, int up_y );



	static BOOL_T ReadFITSFile( Table2D<short>& sample, const char* fname);
	static BOOL_T ReadFITSFile( Table2D<BIG_ELEM_TYPE>& sample, const char* fname);
	static BOOL_T ReadFITSFile( Table2D<float>& sample, const char* fname);
	
	static BOOL_T WriteToFITSFile( Table2D<short>& sample, const char* fname, CKeyTab* pHduTab=NULL);

	static BOOL_T WriteToFITSFile( Table2D<BIG_ELEM_TYPE>& sample, const char* fname, CKeyTab* pHduTab=NULL);

	static BOOL_T WriteToFITSFile( Table2D<float>& float_image, const char* fname, CKeyTab* pHduTab=NULL);


	//static BOOL_T WriteToFITSFileInt( Table2D<int>& sample, const char* fname, CKeyTab* pHduTab=NULL);

	static BOOL_T Normalize( CCDMatrix& image, Table2D<float>& normalized );	

	static BOOL_T CalibrateFrame( CCDMatrix& image, Table2D<float>& calibrated );

	// Adds image to average image kept in out_frame :
	// RETURN VALUES :
	//    TRUE - if OK
	//    FALSE - means error, for example : different sizes of images 
	static BOOL_T AddFrame( Table2D<BIG_ELEM_TYPE>& out_frame, CCDMatrix& in_frame );

	static void NormalizeFrame( Table2D<BIG_ELEM_TYPE>& in_frame, CCDMatrix& out_frame, int cnt );

	static void CalcAndSaveAver( CMyStrTable& file_tab, CMyStrTable& aver_file_list,
							 int first, int i,
							 Table2D<BIG_ELEM_TYPE>* aver_frame,
							 CCDMatrix& sav, float* flat_data, ELEM_TYPE* dark_data, 
							 time_t minTime, time_t maxTime,
							 mystring& szProcessedList, mystring& szOutDir,
							 int& gOutNo,BOOL_T bDARK,BOOL_T bFLAT,
					       BOOL_T bReduced, int x_size, int y_size,
		                BOOL_T gCutBorder, int* gBorderValues,
			             BOOL_T bFH, BOOL_T bFV,mystring& szOutList,
			             const char* szFramesDir="", BOOL_T bSaveFlip=FALSE,
			             CSafeKeyTab* pAddKeys=NULL );


	static int CalcImageShift( const char* szImage1, const char* szImage2, CCDMatrix& dark,
	                    double& dx,double& dy,double& sigma_dx,double& sigma_dy, 
	                    BOOL_T bVerb=FALSE );
};



#endif
