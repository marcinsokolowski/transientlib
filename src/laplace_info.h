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
#ifndef _LAPLACE_INFO_H__
#define _LAPLACE_INFO_H__

#include "mytypes.h"
#include "mypoints.h"
#include "basedefines.h"


int InitGlobalLaplaceInfoTable();

class CLaplaceInfo
{
public:
	CLaplaceInfo();

	CLongPoint m_PlusList[MAX_PIXELS_IN_LAPLACE];
	CLongPoint m_MinusList[MAX_PIXELS_IN_LAPLACE];		
	LONG_T m_PlusCount;
	LONG_T m_MinusCount;
	double m_MinusFactor;
	double m_TypicalSigma;
		
};

extern CLaplaceInfo gLaplaceInfo[MAX_LAPLACE_DEFINED];
extern double gSigmaGCoef[MAX_LAPLACE_DEFINED];
extern double gSigmaBacgroundG[MAX_LAPLACE_DEFINED];
extern double gAverageBacgroundG[MAX_LAPLACE_DEFINED];

// sigma dependence on number of averaged frames :
struct sigmaDep
{
	double mean;
	double sigma;
	double mean_aver;
	double a;
	double b;
	double c;
};

extern sigmaDep gSigmaAverOfPrevNOMAX[MAX_LAPLACE_DEFINED];
extern sigmaDep gSigmaAverOfPrevMAX[MAX_LAPLACE_DEFINED];

double get_sigma( eLaplaceType_T lap, int nPrev, int bDoMax3x3, double sigmaNew );

double get_sigma_and_mean( eLaplaceType_T lap, int nPrev, int bDoMax3x3, 
						double sigmaNew, double& mean );

double get_gt_gauss_prob( double mean, double sigma, double val );
double get_lt_gauss_prob( double mean, double sigma, double val );
double get_significance( double lapNew, double lapPrev, double sigmaNew, double meanNew,
                         eLaplaceType_T lap, int nPrev, int bDoMax3x3, int bCoic );         

#endif
