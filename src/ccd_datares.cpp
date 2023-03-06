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
#include "ccd_datares.h"
#include "ccd_dataconst.h"
#include "ccd_globals.h"
#include <tab2Ddesc.h>
#include <math.h>
#include <ccd_astro_const.h>
#include <mathdefs.h>
#include <mymacros.h>
#include <myutil.h>
#include <mathfunc.h>
#include <cfg.h>

#include "ccd_astro.h"
#include "ccd_procstate.h"

// TABLES ENUMERATED BY :
// eSinglePoint=0, eTwoPoints, eFourPoints, eFivePoints, eFivePlusFourMin,
// eNineEightFive, eNineEightSevenVeryBig,
// eEightFour, eEightTen, eFiveEight, eRawS



CCDDataResults::CCDDataResults()
{}


CCDDataResults::~CCDDataResults()
{}


void CCDDataResults::Initalize()
{
}

double CCDDataResults::GetTreshold( eLaplaceType_T laplace_type, double nSigmaAbove, InfoTable2D* pMatrixInfoTable  )
{
	double mean,sigma;
	double ret = 100000.00;
	if(pMatrixInfoTable->GetMeanAndSigmaBackgr( 0, 0 , mean, sigma, laplace_type )){
		ret = mean + nSigmaAbove*sigma;
	}
	return ret;
}

double CCDDataResults::GetTresholdAver( eLaplaceType_T laplace_type, double nSigmaAbove, InfoTable2D* pMatrixInfoTable )
{
	double n = MIN( gCCDParams.m_nPrevFramesToCheckEdge, gCCDParams.m_nPipelineSize );
	

	return GetTresholdAver( laplace_type, nSigmaAbove, pMatrixInfoTable, n );	
}


double CCDDataResults::GetTresholdAver( eLaplaceType_T laplace_type, double nSigmaAbove, InfoTable2D* pMatrixInfoTable, double n )
{
	double mean,sigma;
	double ret = 100000.00;

	if(pMatrixInfoTable->GetMeanAndSigmaBackgr( 0, 0 , mean, sigma, laplace_type )){
		ret = mean + nSigmaAbove*(sigma/sqrt(n));
	}
	return ret;
	
}

double CCDDataResults::GetSigmaBackground( eLaplaceType_T laplace_type, CCDMatrix* pMatrix, InfoTable2D* pMatrixInfoTable )
{
	/*if(laplace_type==eSinglePoint){
		return 28.26;
	}
	if(laplace_type==eTwoPoints){
		return 44.44;
	}
	if(laplace_type==eFourPoints){
		return 72.67;
	}
	if(laplace_type==eFivePoints){
		return 85.26;
	}
	if(laplace_type==eFivePlusFourMin){
		return 87.00;
	}
	if(laplace_type==eNineEightFive){
		return 110.00;
	}
	if(laplace_type==eNineEightSevenVeryBig){
		return 112.00;
	}
	if(laplace_type==eEightFour){
		return 126.00;
	}
	if(laplace_type==eEightTen){
		return 97.36;
	}*/
	return gSigmaBacgroundG[(long)laplace_type];	
}

double CCDDataResults::GetSigmaBackgroundHomeo( eLaplaceType_T laplace_type, double alpha, 
																CCDMatrix* pMatrix, InfoTable2D* pMatrixInfoTable )
{
	double sigmaG = GetSigmaBackground( laplace_type, pMatrix, pMatrixInfoTable );
	double sigma_homeo = sqrt(alpha/(2-alpha))*sigmaG;
	return sigma_homeo;
}


double CCDDataResults::GetSigmaBackgroundAverageOfPrevN( eLaplaceType_T laplace_type, long prevN, 
																			CCDMatrix* pMatrix, InfoTable2D* pMatrixInfoTable )
{
	double sigmaG = GetSigmaBackground( laplace_type, pMatrix, pMatrixInfoTable );
	double sigma_average_of_prev_n = sigmaG/sqrt((double)prevN);
	return sigma_average_of_prev_n;
}

double CCDDataResults::GetMagFromS( double s )
{
	if(s>0){
		double m = (-1.00/CCDDataConst::Alpha)*log(s/CCDDataConst::S0); 
		return m;
	}
	return 0.00;
}

double CCDDataResults::GetSFromMag( double mag )
{
	if(mag > 0){
		double s = CCDDataConst::S0*exp(-CCDDataConst::Alpha*mag);
		return s;
	}
	return 0.00;
}


double CCDDataResults::GetG54FromMag( double mag )
{
	if(mag > 0){
      double s = CCDDataConst::S0_G54*exp(-CCDDataConst::Alpha_G54*mag);
      return s;
   }
   return 0.00;
}

BOOL_T CCDDataResults::CalcStarPosition( double start_x, double start_y, int nSteps,
                       		              double& curr_x, double& curr_y,
													  double rotCenterX, double rotCenterY, double dAlfa,
													  double dAlfaPerSec, double prevTime, double currTime )
{
	double x0_orginal_prim = start_x-rotCenterX;
   double y0_orginal_prim = start_y-rotCenterY;
	double r1 = sqrt(x0_orginal_prim*x0_orginal_prim+y0_orginal_prim*y0_orginal_prim);

	double sin_alfa_0 = y0_orginal_prim/r1;
	double alfa_0 = CMyMathFunc::GetAngleFromSin( fabs(sin_alfa_0), x0_orginal_prim, y0_orginal_prim );

	// calculation of current position :
	double x0_prim = curr_x-rotCenterX;
   double y0_prim = curr_y-rotCenterY;
// assuming small moves r does not change to much :
//   double r = sqrt(x0_prim*x0_prim + y0_prim*y0_prim);

	double dAlfaTotal=0;
	if(gCCDParams.m_bUseRotPerSec){
		dAlfaTotal = (currTime-prevTime)*dAlfaPerSec;
	}else{
		dAlfaTotal = nSteps*dAlfa;
	}
	double alfa_new = dAlfaTotal+alfa_0;

//	double x0_rot = rotCenterX + r*cos( alfa_new );
// double y0_rot = rotCenterY + r*sin( alfa_new );
	double x0_rot = rotCenterX + r1*cos( alfa_new );
	double y0_rot = rotCenterY + r1*sin( alfa_new );
	if(x0_rot!=curr_x || y0_rot!=curr_y){
   	curr_x = x0_rot;
      curr_y = y0_rot;
		return TRUE;
   }

	return FALSE;	
}


BOOL_T CCDDataResults::CalcStarPosition( double start_x, double start_y, 
													  double prevTime, double currTime, int frame_no,
                       		              double& curr_x, double& curr_y,
													  double rotCenterX, double rotCenterY,
													  double dAlfaPerFrame, double dAlfaPerSec  )
{
	double x0_orginal_prim = start_x-rotCenterX;
   double y0_orginal_prim = start_y-rotCenterY;

	double r = sqrt(x0_orginal_prim*x0_orginal_prim+y0_orginal_prim*y0_orginal_prim);
	double sin_alfa_0 = y0_orginal_prim/(r);
   double alfa_0 = CMyMathFunc::GetAngleFromSin( fabs(sin_alfa_0), x0_orginal_prim, y0_orginal_prim );	

	// calculation of current position :
//	  double x0_prim = start_x-rotCenterX;
//   double y0_prim = start_y-rotCenterY;
//	double r = sqrt(x0_orginal_prim*x0_orginal_prim + y0_prim*y0_prim);


//	double dAlfa = (currTime-prevTime)*sinGeoLatitude*EARTH_ROTATION_OMEGA;
	double dAlfa = 0;
	double dt = ( currTime - prevTime );

	if(gCCDParams.m_bUseRotPerSec){
		dAlfa = dAlfaPerSec*dt;
	}else{
		dAlfa = dAlfaPerFrame*frame_no;
	}


//	_TRACE_PRINTF_6("currTime=%f, prevTime=%f, omega=%f, sin_val=%f ==> dAlfa =%f\n",currTime,prevTime,EARTH_ROTATION_OMEGA,sinGeoLatitude,dAlfa);

	double alfa_new = dAlfa+alfa_0;
	double x0_rot = rotCenterX + r*cos(alfa_new);
   double y0_rot = rotCenterY + r*sin(alfa_new);
	if(x0_rot!=start_x || y0_rot!=start_y){
   	curr_x = x0_rot;
      curr_y = y0_rot;
		return TRUE;
   }

	return FALSE;	
}



BOOL_T CCDDataResults::CalcStarPositionAuto( double start_x, double start_y,
	                                    double prevTime, double currTime, int frame_no,
	                                    double& curr_x, double& curr_y,
	                                    CCDProcState* pCCDInfo, CCDConfig* pCCDParams )
{
	if(gCCDParams.m_bShiftUsesAstroFormulas){
		// using formulas :
		/*pCCDInfo->xyAfterTime( start_x, start_y, (currTime-prevTime), curr_x, curr_y,
												 gCCDParams.m_DAObs, gCCDParams.m_RAObs,
												 gCCDParams.m_TransformCCDOrientation,
												(gCCDParams.m_SizeX/2),(gCCDParams.m_SizeY/2),
												 gCCDParams.m_TransformCCDFocus, gCCDParams.m_TransformCCDPixelSize );*/

		pCCDInfo->xyAfterTimeNEW( start_x, start_y, (currTime-prevTime), curr_x, curr_y );

		//double curr_x_2,curr_y_2;
		//CalcStarPosition( start_x, start_y, prevTime, currTime, frame_no, curr_x_2, curr_y_2, 
		//						pCCDParams->m_RotCenterX, pCCDParams->m_RotCenterY, sinGeoLatitude, 
		//						pCCDParams->m_RotValueDAlfa, pCCDParams->m_RotValueDAlfaPerSec );
		// printf("Astro forulas give (%2.f,%2.f) -> (%2.f,%2.f) - while from data it is (%2.f,%2.f) dT=%d\n",start_x, start_y, curr_x, curr_y, curr_x_2, curr_y_2, (int)(currTime-prevTime) );

	}else{
		// double curr_x_tmp,curr_y_tmp;
		if(pCCDParams && pCCDParams->m_bTransformMatrixOK && gCCDParams.m_bUseTransformMatrix){
			(pCCDParams->m_FrameTransformMatrix).TransformPoint( (currTime-prevTime), start_x, start_y, curr_x, curr_y );
		}else{
			if(pCCDParams->m_bUseRotInAverageOfPrevN && !pCCDParams->m_bDoNotUseRotation){
				CalcStarPosition( start_x, start_y, prevTime, currTime, frame_no, curr_x, curr_y, 
										pCCDParams->m_RotCenterX, pCCDParams->m_RotCenterY, 
										pCCDParams->m_RotValueDAlfa, pCCDParams->m_RotValueDAlfaPerSec );
			}else{
				// curr_x = start_x+(frameDX*(currTime-prevTime));
				// curr_y = start_y+(frameDY*(currTime-prevTime));
				if(gCCDParams.m_bUseRotPerSec){
					curr_x = start_x+(pCCDParams->m_FrameDXPerSec*(currTime-prevTime));
					curr_y = start_y+(pCCDParams->m_FrameDYPerSec*(currTime-prevTime));	
				}else{
					curr_x = my_round(start_x+pCCDParams->m_FrameDX*frame_no);
					curr_y = my_round(start_y+pCCDParams->m_FrameDY*frame_no);
				}
			}
		}
	}
	return TRUE;

}


BOOL_T CCDDataResults::CalcStarPositionAuto( double start_x, double start_y, 
													  double prevTime, double currTime, int frame_no,
                       		              double& curr_x, double& curr_y,
													  double rotCenterX, double rotCenterY,
													  double frameDX, double frameDY,
		                                   BOOL_T bUseRot, double dAlfaPerFrame, CCDProcState* pCCDInfo,
													  double dAlfaPerSec, CCDConfig* pCCDParams  )
{
	if(gCCDParams.m_bShiftUsesAstroFormulas){
		// using formulas :
		/*pCCDInfo->xyAfterTime( start_x, start_y, (currTime-prevTime), curr_x, curr_y,
												 gCCDParams.m_DAObs, gCCDParams.m_RAObs,
												 gCCDParams.m_TransformCCDOrientation,
												(gCCDParams.m_SizeX/2),(gCCDParams.m_SizeY/2),
												 gCCDParams.m_TransformCCDFocus, gCCDParams.m_TransformCCDPixelSize );*/

		pCCDInfo->xyAfterTimeNEW( start_x, start_y, (currTime-prevTime), curr_x, curr_y );

		/*double curr_x_2,curr_y_2;
		CalcStarPosition( start_x, start_y, prevTime, currTime, frame_no, curr_x_2, curr_y_2, 
								rotCenterX, rotCenterY, dAlfaPerFrame, dAlfaPerSec );*/
		// printf("Astro forulas give (%2.f,%2.f) -> (%2.f,%2.f) - while from data it is (%2.f,%2.f) dT=%d\n",start_x, start_y, curr_x, curr_y, curr_x_2, curr_y_2, (int)(currTime-prevTime) );

	}else{
		// double curr_x_tmp,curr_y_tmp;
		if(pCCDParams && pCCDParams->m_bTransformMatrixOK && gCCDParams.m_bUseTransformMatrix){
			(pCCDParams->m_FrameTransformMatrix).TransformPoint( (currTime-prevTime), start_x, start_y, curr_x, curr_y );
		}else{
			if(bUseRot){
				CalcStarPosition( start_x, start_y, prevTime, currTime, frame_no, curr_x, curr_y, 
										rotCenterX, rotCenterY, dAlfaPerFrame, dAlfaPerSec );
			}else{
				// curr_x = start_x+(frameDX*(currTime-prevTime));
				// curr_y = start_y+(frameDY*(currTime-prevTime));
				curr_x = my_round(start_x+frameDX*frame_no);
				curr_y = my_round(start_y+frameDY*frame_no);
			}
		}
	}
	return TRUE;
}

BOOL_T CCDDataResults::CalcStarPositionsFromFormula( double start_x, double start_y,
													int* prevX, int* prevY,
													int backStart, int backTo,
													const double* PrevFramesTime,int time_sign,
													CCDProcState* pCCDInfo, CCDConfig* pCCDParams )
{
	// using formulas :
	pCCDInfo->xyAfterTimeNEW( start_x, start_y, prevX, prevY, backStart, backTo, PrevFramesTime, time_sign );
	return TRUE;
}
