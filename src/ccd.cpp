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
#include "ccd.h"
#include "ccd_trace.h"
#include "ccd_globals.h"
#include "ccd_pipeline.h"
#include <cfg.h>
#include <cexcp.h>
#include<mydate.h>
#include<errno.h>

#include <limits.h>

cCCD::cCCD():m_pCCD(NULL),m_FrameIndex(NOT_DEFINED),m_pPipeline(NULL)
{
	// gcc4.0 : requires commenting lines below, so that no 
	// allocation is done in this constructor !
	// this is due to fact that m_bAllocHere will be set to default TRUE value
	// while it should be set in InitConstructor to FALSE !!!
	/*m_SizeX = gCCDParams.m_SizeX;
	m_SizeY = gCCDParams.m_SizeY;
	m_Count = gCCDParams.m_nCamNo;
	Init( m_SizeX, m_SizeY, m_Count );*/
}

void cCCD::cCCD_InitConstructor( long xSize, long ySize, long count, 
			   LONGLONG_T idx,
				BOOL_T bLaplaceOfCurrent/*=FALSE*/,CCDPipeline* pPipeline/*=NULL*/,
				BOOL_T bAllocHere/*=TRUE*/, BOOL_T bDoAllocEvents/*=TRUE*/ )
{
	m_FrameIndex = idx;
 	m_pPipeline = pPipeline;

	Init( xSize, ySize, count, bLaplaceOfCurrent, bAllocHere, bDoAllocEvents );
}


cCCD::cCCD(long xSize, long ySize, long count, 
			   LONGLONG_T idx,
				BOOL_T bLaplaceOfCurrent/*=FALSE*/,CCDPipeline* pPipeline/*=NULL*/,
				BOOL_T bAllocHere/*=TRUE*/, BOOL_T bDoAllocEvents/*=TRUE*/ )
:m_SizeX(xSize),m_SizeY(ySize),m_Count(count),m_pCCD(NULL),m_FrameIndex(idx),
 m_pPipeline(pPipeline)
{
	Init( m_SizeX, m_SizeY, m_Count, bLaplaceOfCurrent, bAllocHere, bDoAllocEvents );
}

cCCD::cCCD(const cCCD& right):m_pCCD(NULL)
{
	(*this) = right;
}

cCCD::~cCCD()
{
	delete [] m_pCCD;
}

void cCCD::SetFileNames( const char* filenames )
{
	m_szFileNames = filenames;
}

void cCCD::InitFrame()
{
	for(int i=0;i<m_Count;i++){
		m_pCCD[i].SetData( 0 );
	}
}

BOOL_T cCCD::Init( long xSize, long ySize,long nCCD, 
						 BOOL_T bLaplaceOfCurrent/*=FALSE*/,
						 BOOL_T bAllocHere/*=TRUE*/,
						 BOOL_T bDoAllocEvents/*=TRUE*/ )
{
	m_SizeX = xSize;
	m_SizeY = ySize;
	m_Count = nCCD;
	if(m_pCCD)
		delete [] m_pCCD;
	if(nCCD){ 
		int nInitEventBufferSize = gCCDParams.GetEventsBufferSize();
		if(!bDoAllocEvents){
			nInitEventBufferSize = 0;
		}
		
//		m_pCCD = new CCDMatrix[nCCD](m_SizeX,m_SizeY,TRUE,NOT_DEFINED,this,
//											  nInitEventBufferSize,
//											  bLaplaceOfCurrent, NULL, bAllocHere );
// version for gcc4.0 - NEW :
		m_pCCD = new CCDMatrix[nCCD];

		for(int i=0;i<nCCD;i++){
//			BOOL_T bOK = m_pCCD[i].Alloc(m_SizeX,m_SizeY);
// version for gcc4.0 - NEW :
			BOOL_T bOK = m_pCCD[i].CCDMatrix_InitConstructor(
												m_SizeX,m_SizeY,TRUE,NOT_DEFINED,this,			
												nInitEventBufferSize,
                                    bLaplaceOfCurrent, NULL, bAllocHere );

			m_pCCD[i].SetIndex(i); 				
			
			CCcdCfg* pMatrixParams = NULL;
	      if(m_pPipeline){
   	      pMatrixParams = m_pPipeline->GetCamParamSet(i); 
      	}
			m_pCCD[i].SetParamSet( pMatrixParams );		
		
			if (!bOK){
				printf("Memory alloction problem occured, when allocating memory for CCD number %d\n",i);
				printf("Errno = %d\n",errno);
				printf("Size  = %dMB\n",(m_SizeX*m_SizeY)*sizeof(ELEM_TYPE)/1000000);
				printf("Desc  = %s\n",strerror(errno));
				MYTRACE1(gCCDTrace,"Memory allocation problem occured !!!");
				exit(-1);
			}else{
				MYTRACE3(gCCDTrace,"Allocated matrix number " << i);
			}
		}
	}	
	return TRUE;
}

cCCD& cCCD::operator=(const cCCD& right)
{
	if(right.m_SizeX!=m_SizeX || right.m_SizeY!=m_SizeY || 
      right.m_Count!=m_Count)
		Init( right.m_SizeX, right.m_SizeY, right.m_Count );	
	for(int i=0;i<m_Count;i++){
		m_pCCD[i] = right.m_pCCD[i];
	}
	return (*this);
}

void cCCD::Subtract(cCCD& right,cCCD& result,BOOL_T bZero)
{
	Assert(m_Count==right.m_Count && right.m_Count==result.m_Count,"Number of CCD matrixes must be equal in subtract");

	clock_t t1 = clock();
	for(int i=0;i<m_Count;i++){
		m_pCCD[i].Subtract(right.m_pCCD[i],result.m_pCCD[i],bZero);
	}
	clock_t t2 = clock();
	mystring msg = get_clock_in_sec_string( t2-t1 );
   MYTRACE2(gCCDTrace,"Subtraction of images took : " << msg);
}

void cCCD::Divide(cCCD& right,cCCD& result)
{
	Assert(m_Count==right.m_Count && right.m_Count==result.m_Count,"Number of CCD matrixes must be equal in divide");

	clock_t t1 = clock();
	for(int i=0;i<m_Count;i++){
		m_pCCD[i].Divide(right.m_pCCD[i],result.m_pCCD[i],USHRT_MAX);
	}
	clock_t t2 = clock();
	mystring msg = get_clock_in_sec_string( t2-t1 );
   MYTRACE2(gCCDTrace,"Dividing of images took : " << msg);
}

void cCCD::Normalize()
{
	clock_t t1 = clock();
	for(int i=0;i<m_Count;i++){
		m_pCCD[i].Normalize();
	}		
	clock_t t2 = clock();
	MYTRACE1(gCCDTrace,"Normalization took :  " << (t2-t1) << " ticks (=" << (t2-t1)/CLOCKS_PER_SEC << " sec)");
}

void cCCD::InitRandom(long min_value,long max_value)
{
	time_t t1 = get_dttm();
	for(int i=0;i<m_Count;i++){
		m_pCCD[i].InitRandom(min_value,max_value);
	}	
	time_t t2 = get_dttm();
	MYTRACE3(gCCDTrace,"Random initialization of ccd frame took " << (t2-t1) << "sec");
}


void cCCD::InitRandomDarkFrame()
{
	time_t t1 = get_dttm();
	for(int i=0;i<m_Count;i++){
		m_pCCD[i].InitRandom(50,70);
	}
	Normalize();
	time_t t2 = get_dttm();	
	MYTRACE3(gCCDTrace,"Random initialization dark frame took  " << (t2-t1) << "sec");
}

CCDMatrix* cCCD::GetMatrix(const int pos)
{
	Assert(pos>=0 && pos<m_Count,"CCD aparatus has only %d cameras, you want to access camers number %d",m_Count,pos);
	return &(m_pCCD[pos]);
}

CCDMatrix& cCCD::operator[](const int pos)
{
	Assert(pos>=0 && pos<m_Count,"CCD aparatus has only %d cameras, you want to access camers number %d",m_Count,pos);
	return m_pCCD[pos];
}

void cCCD::SetData(ELEM_TYPE val)
{
	time_t t1 = get_dttm();
	for(int i=0;i<m_Count;i++){
		m_pCCD[i].SetData(val);
	}
	time_t t2 = get_dttm();	
	MYTRACE3(gCCDTrace,"Random initialization dark frame took  " << (t2-t1) << "sec");

}

void cCCD::WriteFrameToFile(const char* basename,BOOL_T bDumpEvents)
{
	mystring szFullName;

	for(int i=0;i<m_Count;i++){
		szFullName = basename;
		szFullName << "_Camera_" << i;
		m_pCCD[i].WriteToFITSFile( szFullName.c_str(), bDumpEvents );
	}
}

BOOL_T cCCD::ClearState()
{
	m_bIntersting=FALSE;
	for(int i=0;i<m_Count;i++){
		m_pCCD[i].ClearState();
	}		
	return TRUE;
}


BOOL_T cCCD::ReadFrame( const char* dark )
{
	mystring szBase = dark;
	szBase.env2str();
	BOOL_T bRet=FALSE;
	if(m_Count==1){
		bRet = m_pCCD[0].ReadFITSFile( szBase.c_str() );
	}else{
		for(int i=0;i<m_Count;i++){
			mystring szTrace;
			szTrace << "In case frame consists of more then one image use syntax with %% sign";
	
			Assert(strstr(dark,"%")!=NULL,szTrace.c_str());
			char path[1000];
			sprintf(path,szBase.c_str(),i);
			bRet = bRet && m_pCCD[0].ReadFITSFile( path );
		}		
	}
	return bRet;
}


void cCCD::AdjustFrames()
{
	for(int i=0;i<m_Count;i++){
		m_pCCD[0].AdjustFrame();
	}
}


cCCDCache::cCCDCache()
: m_pCache(NULL),m_CurIdx(0),m_Count(0)
{
}

cCCDCache::~cCCDCache()
{
	if(m_pCache)
		delete [] m_pCache;
}

void cCCDCache::Init( long xSize, long ySize, long count, LONG_T size )
{
	if(m_pCache && size!=m_Count){
		delete [] m_pCache;
		m_pCache = NULL;
	}
	m_Count = size;
	if(!m_pCache){
		m_pCache = new cCCD[ m_Count ];
		for(int i=0;i<m_Count;i++){
			m_pCache[i].cCCD_InitConstructor( xSize, ySize, count, NOT_DEFINED, 
						FALSE,
						NULL, 
						TRUE,
						FALSE );
		}
	}
}

cCCD* cCCDCache::GetCurrent()
{
	if(m_pCache && m_CurIdx<m_Count){
		return &(m_pCache[m_CurIdx]);
	}
	return NULL;
}

cCCDCache& cCCDCache::operator++(int)
{
	m_CurIdx++;	
	return (*this);
}

void cCCDCache::Reset()
{
	m_CurIdx = 0;	
}
