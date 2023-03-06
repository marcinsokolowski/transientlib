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
#include "ccd_pipeline.h"
#include "ccd_controller.h"
#include <Astroangle.h>
#include <mycmnglobals.h>
#include <myutil.h>
#include "ccd_photometry.h"

// flags events which occured in same place as some recent GRB :
void CCDPipeline::CheckExternalTriggers( CCDEventList& finallist1, CCDEventList& finallist2 )
{
}

int CCDPipeline::InsertEventsToDB( CCDEventList& finallist1, CCDEventList& finallist2,
											  CCDPipeline* pPipeline1, CCDPipeline* pPipeline2,
											  BOOL_T bForce )
{
	return 0;	
}

int CCDPipeline::InsertEventsToDB( CCDEventList& finallist1, CCDPipeline* pPipeline1,
											  BOOL_T bForce )
{
	return 0;	
}



int CCDPipeline::RunFastPhotometry( CCDMatrix& matrix, double tresh, 
												const char* szMagFile, BOOL_T bSaveMag /*=TRUE*/ )
{
	int bRet = 0;
	if( bSaveMag ){
		bRet = CPiPhotometry::ExecPhotometryAndSave( matrix, szMagFile, 
						tresh, tresh, m_PipelineCfg.m_nAsasBorderSize, 12 );
	}else{
		printf("Fast photometry without dumping mag file\n");
		bRet = CPiPhotometry::ExecPhotometryNoSave( matrix, szMagFile, 
					tresh, tresh, m_PipelineCfg.m_nAsasBorderSize, 12 );
	}
	return bRet;
}

int CCDPipeline::ReadSingleFrameEvents()
{
	return 0;
}

int CCDPipeline::ReadSingleFrameTracks()
{
	return 0;
}

