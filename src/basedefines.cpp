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
#include "basedefines.h"
#include <string.h>

ImageStat::ImageStat()
:Average(0),Min(0),Max(0),MinNonZero(0),Sum(0),m_pValDistrib(NULL),
 m_ValCount(0),MaxX(0),MaxY(0),RMS(0),MeanFit(0),SigmaFit(0),
 bFitOK(FALSE)
{
}


ImageStat::~ImageStat()
{
	if(m_pValDistrib)
		delete [] m_pValDistrib;
}



ImageStat::ImageStat( const ImageStat& right )
{
	if(m_pValDistrib)
      delete [] m_pValDistrib;
	if(right.m_pValDistrib && right.m_ValCount){
		m_ValCount = right.m_ValCount;
		m_pValDistrib = new long[m_ValCount];
		memcpy(m_pValDistrib,right.m_pValDistrib,sizeof(long)*m_ValCount);
	}
}

void ImageStat::AllocDistribTab( long maxValue )
{
	if(m_pValDistrib)
		delete [] m_pValDistrib;
	m_ValCount = maxValue;
	m_pValDistrib = new long[m_ValCount];
	memset( m_pValDistrib, 0, sizeof(long)*m_ValCount );
}

LONG_T ImageStat::GetMaxOccurValue()
{
	LONG_T MaxCounter=0;
	LONG_T MaxOccurVal=0;
	for(int i=0;i<m_ValCount;i++){
		if(m_pValDistrib[i]>MaxCounter){
			MaxCounter = m_pValDistrib[i];
			MaxOccurVal = i;		
		}				
	}
	return MaxOccurVal;
}


const char* GetFlipDesc( eDriverReverseImage_T fliptype )
{
	if( fliptype==eReverseImageFull ){
		return "FULL";
	}
	if( fliptype==eReverseImageHor ){
		return "FH";
	}
	if( fliptype==eReverseImageVert ){
		return "FV";
	}
	if( fliptype==eReverseImageHorVert ){
		return "FVH";
	}
	return "NONE";
}