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
#ifndef _CCD_SAMPLES_H__
#define _CCD_SAMPLES_H__

#include <stdlib.h>
#include <myfile.h>
#include <tab2D.h>
#include <vector>
#include <basestructs.h>

using namespace std;

class CMyStrTable;
class CDescFile;

#define MAX_MAGNITUDE 30


class CMagSamples
{
protected:
	CEnvVar m_SamplesList;
	CMyStrTable* m_pSampleFilesTab;
	Table2D<ELEM_SAMPLE_TYPE>* m_pSamples;		
	LONG_T m_Count;
	BOOL_T m_bRead;

public:
	mystring m_szMagnitude;

	CMagSamples();
	~CMagSamples();
	Table2D<ELEM_SAMPLE_TYPE>* GetSamples(){ return m_pSamples; }
	LONG_T GetCount(){ return m_Count; }
	void SetMagnitude( const char* mag ){ m_szMagnitude=mag; }

	LONG_T ReadSamples( const CEnvVar& sampleDesc,const char* list_file );
	
	Table2D<ELEM_SAMPLE_TYPE>& GetRandomSample( int* pIdx=NULL );		
	BOOL_T GetRandomSample_CheckXY( int x, int y, Table2D<ELEM_SAMPLE_TYPE>& sample_out,
											  int nFrame=-1 );

	BOOL_T GetRandomSample_ByName( int x, int y, 
											 Table2D<ELEM_SAMPLE_TYPE>& sample_in,
											 Table2D<ELEM_SAMPLE_TYPE>& sample_out,
											 int nFrame=-1 );

	static BOOL_T CheckXY( int x, int y, int sample_x, int sample_y, int nFrame=-1 );
};


class CccdSamples 
{
protected:
	CMagSamples* m_pMagSamples;
	
	BOOL_T m_bRead;
public:
	CDescFile m_DescFile;		
	double* m_pMagTab;
		
	CccdSamples();
	~CccdSamples();
	int GetCount(){ return m_DescFile.GetDescTab().GetCount(); }
	BOOL_T ReadSamples( const char* samples_dir );
	CMagSamples* GetMagSamples( const char* szMag, int* pPos=NULL );		
	CMagSamples* GetMagSamples( int pos );		
	Table2D<ELEM_SAMPLE_TYPE>& GetRandomSample( const char* szMag=NULL, mystring* szSampleMag=NULL );
	// vector<CMagSamples>& GetMagSamples(){ return m_MagSamples; }	
};


#endif
