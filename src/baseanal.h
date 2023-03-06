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
#ifndef _BASE_ANAL_H__
#define _BASE_ANAL_H__

#include "basedefines.h"
#include "mytypes.h"
#include "mypoints.h"


#define G1_NAME "g1"
#define G2_NAME "g2"
#define G4_NAME "g4"
#define G5_NAME "g5"
#define G54_NAME "g54"
#define G985_NAME "g985"
#define G987_NAME "g987"
#define G84_NAME "g84"
#define G810_NAME "g810"
#define G58_NAME "g58"
#define S_NAME "s"
#define G412_NAME "g4_12"
#define G412FAR_NAME "g4_12Far"

class CLongPoint;

class CBaseAnal 
{
public :
	CBaseAnal( int sizeX, int sizeY );
	
	int m_SizeX;
	int m_SizeY;

/*	static double CalcG54( LONG_T x, LONG_T y, LONG_T xSize, ELEM_TYPE** p_data );


	static LONG_T CalcLaplaceSum( LONG_T x, LONG_T y, LONG_T xSize, ELEM_TYPE** p_data,
										eLaplaceType_T laplaceType );										

	static LONG_T CalcLaplaceSumBig( LONG_T x, LONG_T y, LONG_T xSize, BIG_ELEM_TYPE** p_data,
										eLaplaceType_T laplaceType );										

	static LONG_T CalcLaplaceSum( LONG_T x, LONG_T y, LONG_T xSize, ELEM_TYPE** p_data,
										eLaplaceType_T laplaceType, LONG_T& plus_sum, LONG_T& minus_sum );

	static LONG_T CalcLaplaceSumBig( LONG_T x, LONG_T y, LONG_T xSize, BIG_ELEM_TYPE** p_data,
										eLaplaceType_T laplaceType, LONG_T& plus_sum, LONG_T& minus_sum );

	static LONG_T CalcLaplaceSumEstimateOnEdge( LONG_T x, LONG_T y, 
															  LONG_T xSize, LONG_T ySize,
   								               		ELEM_TYPE** p_data,
   								               		CLongPoint* plus_list, LONG_T plus_cnt,
   								               		CLongPoint* minus_list, LONG_T minus_cnt );

	static void GetPlusMinusList( eLaplaceType_T laplaceType,
											CLongPoint* plus_list, LONG_T& plus_cnt,
					  					   CLongPoint* minus_list, LONG_T& minus_cnt );*/
					  					   
					  					   
	static const char* GetLaplaceName( eLaplaceType_T laplaceType );					  					   
};


#endif
