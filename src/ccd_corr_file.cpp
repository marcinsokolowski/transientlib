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
#include "ccd_corr_file.h"
#include <myfile.h>
#include <myparser.h>
#include <mystrtable.h>

#define STEP_COUNT "STEP_COUNT"


CCDFastPhotoCorrFile::CCDFastPhotoCorrFile()
: m_CorrTab(NULL), StepCount(-1), StepSize(-1)
{
/*	if( !Read( m_szFileName.c_str() ) ){
		printf("could not read correction file : %s , exiting\n",m_szFileName.c_str());
		exit(-1);		
	}*/
}

CCDFastPhotoCorrFile::~CCDFastPhotoCorrFile()
{
	if( m_CorrTab ){
		delete [] m_CorrTab;
	}
}
                                                                                
BOOL_T CCDFastPhotoCorrFile::Read( const char* szCorrFile )
{
	if( !szCorrFile || strlen( szCorrFile )==0 || !MyFile::DoesFileExist( szCorrFile ) ){
		printf("file %s not found\n",szCorrFile);
		return FALSE;
	}

	printf("Reading fast photometry correction file %s ...\n",szCorrFile);fflush(0);

	m_szFileName = szCorrFile;

	MyIFile in( szCorrFile );
   const char* pLine=NULL;
	int count=0;
                                                                                
   while( pLine = in.GetLine( TRUE ) ){
      if( mystring::get_first_non_white( pLine )=='#' ){
			if(strstr( pLine, STEP_COUNT ) ){
				const char* ptr=pLine;
				while( (*ptr) && ( (*ptr)=='#' || (*ptr)==' ' ) )ptr++;
				MyParser pars = ptr;

				mystring szName,szValue;
				pars.GetVarAndValue(	szName, szValue );
				if( strcmp( szName.c_str() , STEP_COUNT )==0 ){
					StepCount = atol( szValue.c_str() );
					StepSize = 1.00/(double(StepCount));					
					m_CorrTab = new float[StepCount*StepCount];
				}
			}
			continue;
		}

		float x,y,corr;

		CMyStrTable items;
		MyParser pars=pLine;
		pars.GetItems( items );
	
		if( items.size()>=3 ){
			x = atof( items[0].c_str() )+StepSize/5.00;
			y = atof( items[1].c_str() )+StepSize/5.00;
			corr = atof( items[2].c_str() );

			int pos_x = ( x / StepSize );
			int pos_y = ( y / StepSize );
			
			int pos = pos_y*StepCount+pos_x;
			m_CorrTab[pos] = corr;
		}

		count++;
	}
	
	BOOL_T bOK=FALSE;
	if( count==(StepCount*StepCount) && StepCount>0 ){
		bOK = TRUE;
	}else{
		if( count!=(StepCount*StepCount) ){
			printf("Number of found numbers %d different then expected %d\n",count,(StepCount*StepCount));
		}
		if( StepCount<=0 ){
			printf("Item %s not found in comment of file !\n",STEP_COUNT);
		}
	}


	int size = 0;
	for(int i=0;i<size;i++){
		if( m_CorrTab[i] <= 0 ){
			printf("Error 0 value in corr table !\n");
			exit(0);
		}
	}

	printf("read %d numbers , expected %d , OK=%d\n",count,(StepCount*StepCount),bOK);
	return bOK;
}

float CCDFastPhotoCorrFile::Correct( float x, float y, float adu )
{
	int x_pixel = (int)x;
	int y_pixel = (int)y;
	float x_partial = (x-x_pixel);
	float y_partial = (y-y_pixel);

	int pos_x = ( x_partial / StepSize );
	int pos_y = ( y_partial / StepSize );
			
	int pos = pos_y*StepCount+pos_x;
	
	if( pos>=(StepCount*StepCount) ){
		printf("(%.2f,%.2f) - exceeded correction array size !\n",x,y);
		exit(0);
	}

	float corrected = ( m_CorrTab[pos]*adu );

	if( corrected == 0 ){
		printf("ERROR ! corrected = 0 !!\n");
		exit(0);
	}	

	return corrected;
}

