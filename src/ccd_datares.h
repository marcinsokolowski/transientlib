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
#ifndef _CCD_DATA_RESULTS_H__
#define _CCD_DATA_RESULTS_H__

#include "ccd_defines.h"
#include <laplace_info.h>

// TABLES ENUMERATED BY :
// eSinglePoint=0, eTwoPoints, eFourPoints, eFivePoints, eFivePlusFourMin,
// eNineEightFive, eNineEightSevenVeryBig,
// eEightFour, eEightTen, eFiveEight 


class CCDMatrix;
class CCDProcState;
class CCDParams;
class CCDConfig;

class CCDDataResults
{
public :
	CCDDataResults();
	~CCDDataResults();
	static void Initalize();

	static double GetTreshold( eLaplaceType_T laplace_type, double nSigmaAbove, InfoTable2D* pMatrixInfoTable  );
	static double GetTresholdAver( eLaplaceType_T laplace_type, double nSigmaAbove, InfoTable2D* pMatrixInfoTable );
	static double GetTresholdAver( eLaplaceType_T laplace_type, double nSigmaAbove, 
											 InfoTable2D* pMatrixInfoTable, double n );
	
	static double GetSigmaBackground( eLaplaceType_T laplace_type, CCDMatrix* pMatrix=NULL, 
												 InfoTable2D* pMatrixInfoTable=NULL );
												 
	static double GetSigmaBackgroundHomeo( eLaplaceType_T laplace_type, double alpha, CCDMatrix* pMatrix=NULL,
													   InfoTable2D* pMatrixInfoTable=NULL );
	
	static double GetSigmaBackgroundAverageOfPrevN( eLaplaceType_T laplace_type, long prevN, CCDMatrix* pMatrix=NULL,
																	InfoTable2D* pMatrixInfoTable=NULL );

	static double GetMagFromS( double s );
	static double GetSFromMag( double mag );


	static double GetG54FromMag( double mag );	

	static BOOL_T CalcStarPosition( double start_x, double start_y, int nSteps,
	                                double& curr_x, double& curr_y,
	                                double rotCenterX, double rotCenterY, double dAlfa, 
	                                double dAlfaPerSec, double prevTime, double currTime );	

	static BOOL_T CalcStarPosition( double start_x, double start_y,
	                                double prevTime, double currTime, int frame_no,
	                                double& curr_x, double& curr_y,
	                                double rotCenterX, double rotCenterY,
	                                double dAlfaPerFrame, double dAlfaPerSec );
	                                                                                  
	static BOOL_T CalcStarPositionAuto( double start_x, double start_y,
	                                    double prevTime, double currTime, int frame_no,
	                                    double& curr_x, double& curr_y,
	                                    double rotCenterX, double rotCenterY,
	                                    double frameDX, double frameDY, 
	                                    BOOL_T bUseRot, double dAlfaPerFrame, CCDProcState* pCCDInfo,
	                                    double dAlfaPerSec, CCDConfig* pCCDParams=NULL );
	                                                                                  
	static BOOL_T CalcStarPositionAuto( double start_x, double start_y,
	                                    double prevTime, double currTime, int frame_no,
	                                    double& curr_x, double& curr_y,
	                                    CCDProcState* pCCDInfo, CCDConfig* pCCDParams );
	                                                                                  
	static BOOL_T CalcStarPositionsFromFormula( double start_x, double start_y,
													int* prevX, int* prevY,
													int backStart, int backTo,
													const double* PrevFramesTime, int time_sign,
													CCDProcState* pCCDInfo, CCDConfig* pCCDParams );
	                                                                                  
											  
};


#endif
