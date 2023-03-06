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
#ifndef _CCD_H__
#define _CCD_H__

#include "ccd_matrix.h"
#include <mystring.h>


class CCDPipeline;

class cCCD {
protected:
	CCDMatrix* m_pCCD;
	long m_SizeX;
	long m_SizeY;
	long m_Count;
	double m_Average;
	LONGLONG_T m_FrameIndex;
	mystring m_szFileNames;
	CCDPipeline* m_pPipeline;

	// processing :
	BOOL_T m_bIntersting;
public :
	cCCD();
	~cCCD();
	cCCD( long xSize, long ySize, long count,
			LONGLONG_T idx=NOT_DEFINED, BOOL_T bLaplaceOfCurrent=FALSE,
			CCDPipeline* pPipeline=NULL, BOOL_T bAllocHere=TRUE,
			BOOL_T bDoAllocEvents=TRUE );
	cCCD(const cCCD& right);

	void cCCD_InitConstructor( long xSize, long ySize, long count, 
			LONGLONG_T idx=NOT_DEFINED, BOOL_T bLaplaceOfCurrent=FALSE,
			CCDPipeline* pPipeline=NULL, BOOL_T bAllocHere=TRUE,
			BOOL_T bDoAllocEvents=TRUE );

	inline CCDMatrix* GetFrames(){ return m_pCCD; } 
	BOOL_T Init( long xSize, long ySize, long nCCD, 
					 BOOL_T bLaplaceOfCurrent=FALSE, BOOL_T bAllocHere=TRUE, BOOL_T bDoAllocEvents=TRUE );
	BOOL_T ClearState();
	void InitFrame();

	void SetFileNames( const char* filenames );
	const char* GetFileNames() { return m_szFileNames.c_str(); }
	
	// testing functions - can be commented out for realese version
	void InitRandom(long min_value=0,long max_value=255);
	void InitRandomDarkFrame();

	inline void SetFrameIndex( LONGLONG_T idx ){ m_FrameIndex = idx; }
	inline LONGLONG_T GetFrameIndex() { return m_FrameIndex; }
	inline long GetCount() const { return m_Count; }
	void SetData(ELEM_TYPE val);
	void Normalize();
	CCDMatrix& operator[](const int pos);
	CCDMatrix* GetMatrix(const int pos);
	cCCD& operator=(const cCCD& right);		


	void Subtract(cCCD& right,cCCD& result,BOOL_T bZero=FALSE);
	void Divide(cCCD& right,cCCD& result);
	
	
	void AdjustFrames();
	

	// flaging intersting frames :
	inline void SetInteresting(BOOL_T bFlag=TRUE){ m_bIntersting=bFlag; }	
	inline BOOL_T GetInteresting(){ return m_bIntersting; }
	
	// reading frame :
	BOOL_T ReadFrame( const char* dark );
	
	// writing frame to files :
	void WriteFrameToFile(const char* basename,BOOL_T bDumpEvents=FALSE);
};


class cCCDCache 
{
public :
	cCCDCache();
	~cCCDCache();
	void Init( long xSize, long ySize, long count, LONG_T size );	
	cCCD* GetCurrent();
	cCCDCache& operator++(int);
	void Reset();
	
	cCCD* m_pCache;
	LONG_T m_CurIdx;		
	LONG_T m_Count;
};

#endif
