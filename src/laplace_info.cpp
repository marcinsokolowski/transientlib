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
#include "laplace_info.h"
#include "tab2D.h"
#include  "mymacros.h"
#include <math.h>
#include "mathdefs.h"

CLaplaceInfo gLaplaceInfo[MAX_LAPLACE_DEFINED];

double gSigmaBacgroundG[MAX_LAPLACE_DEFINED] = { 
										 28.26, 44.44, 72.67, 85.26, 87.00, 
										 110.00, 112.00, 165.00, 136.36, 100.00, 
										 80.00, 80.00, 80.00, 80.00, 30, 30, 30, 30, 30, 30 };

double gAverageBacgroundG[MAX_LAPLACE_DEFINED] = 
										{ 0, 0, 0, 0, 0, 
										  0, 0, 0, 0, 0, 
									 1680 , 0, 0, 0, 0,
										  0, 0, 0, 0, 0 };
										  

double gSimgaGCoef[MAX_LAPLACE_DEFINED] = 
									{ 5.00/4.00, 3.00, 8.00, 45.00/4.00, 45.00/4.00,
                             153.00/8.00, 153.00/8.00, 72.00, 328.00/25.00, 
									  65.00/8.00, 1.00, 0.3333333, 0.33333,
                             0, 0, 0, 0, 0, 0, 0 };



sigmaDep gSigmaAverOfPrevNOMAX[MAX_LAPLACE_DEFINED]=
{ 
	{   -2.5100,  32.6100,  -1.9693,  19.5415,  0.0000,   8.8520 },
	{   -5.3300,  53.2100,  -4.8460,  31.9486,  0.0000,   14.4320 },
	{   -10.5700, 92.2200,  -10.5980, 57.2245,  0.0000,   26.1203 },
	{   -5.0800,  98.4000,  -4.9127,  56.4237,  0.0000,   24.9643 },
	{  -13.8000, 114.6200, -14.3620, 82.2722, -0.0241,  33.3958 },
	{ -23.2100, 143.8500, -21.1807, 122.5045, -0.2176,  32.5463 },
	{ -31.0200, 156.1700, -30.1700, 99.0706 , -0.0774,  60.2787 },
	{  -20.5800, 168.2000, -21.2853, 99.0836,  0.0880,   65.0320 },
	{ -30.9400, 142.6200, -30.2427, 87.8115 , -0.1321,  59.9863 },
	{  -19.4100, 102.8800, -19.5313, 68.0852,  0.0181,   34.3616 },
	{ 1772.9700,   35.5500,  1764.0633,   13.0518,  0.3713,   19.2951 },
	{   -10.7900, 73.9200,  -9.6233,  44.9498,  0.0000,  20.3901 },
	{   -15.1700, 80.3200,  -14.4787, 59.1177,  -0.0176,  21.7035 }
};


sigmaDep gSigmaAverOfPrevMAX[MAX_LAPLACE_DEFINED]=
{
	{ -2.5100 , 32.6100 , 46.0453 , 6.0438 , 0.5824 , 9.7271 },
	{ -5.3300 , 53.2100 , 67.9853 , 21.7396 , 0.0950 , 8.6023 },
	{ -10.5700 , 92.2200 , 110.3487 , 42.8720 , -0.0567 , 13.8393 },
	{ -5.0800 , 98.4000 , 136.3607 , 31.4776 , 0.0000 , 14.2954 },
	{ -13.8000 , 114.6200 , 127.6500 , 46.9597 , -0.0795 , 25.3439 },
	{ -23.2100 , 143.8500 , 167.1253 , 53.5523 , -0.0175 , 29.0915 },
	{ -31.0200 , 156.1700 , 153.0300 , 49.0228 , 0.1812 , 44.0812 },
	{ -20.5800 , 168.2000 , 178.2220 , 54.2184 , 0.1398 , 40.2966 },
	{ -30.9400 , 142.6200 , 128.7927 , 46.8879 , 0.0607 , 46.0225 },
	{ -19.4100 , 102.8800 , 102.5760 , 40.1110 , 0.1052 , 27.6774 },
	{ 1772.9700 , 35.5500 , 1806.3220 , 5.3800 , 0.4404 , 22.4234 },
	{ -10.7900 , 73.9200 , 89.4007 , 34.2933 , -0.0202 , 11.7484 },
	{ -15.1700 , 80.3200 , 82.9167 , 37.0713 , -0.0259 , 17.8688 }
};


double sigma_func_aver_n( double a, double b, double c, double n )
{
	double bottom = sqrt( n - b );
	double ret = (a/bottom)+c;
	return ret;
}


double sigma_func_aver_n( double a, double b, double c, double n, 
								  double sigma, double sigmaNew )
{
	a = (a/sigma)*sigmaNew;
	b = (b/sigma)*sigmaNew;
	c = (c/sigma)*sigmaNew;

	double bottom = sqrt( n - b );
	double ret = (a/bottom)+c;
	return ret;
}


double get_sigma_and_mean( eLaplaceType_T lap, int nPrev, int bDoMax3x3, 
									double sigmaNew, double& mean )
{
	if(bDoMax3x3){
		if(nPrev==1){
			mean = gSigmaAverOfPrevMAX[lap].mean;
			return gSigmaAverOfPrevMAX[lap].sigma;
		}else{
			mean = gSigmaAverOfPrevMAX[lap].mean_aver;
			double ret = sigma_func_aver_n( gSigmaAverOfPrevMAX[lap].a, gSigmaAverOfPrevMAX[lap].b, gSigmaAverOfPrevMAX[lap].c, (double)nPrev,
													  gSigmaAverOfPrevMAX[lap].sigma, sigmaNew );
			return ret;
		}
	}else{
		mean = gSigmaAverOfPrevNOMAX[lap].mean;
		if(nPrev==1){
			return gSigmaAverOfPrevNOMAX[lap].sigma;
		}else{
			double ret = sigma_func_aver_n( gSigmaAverOfPrevNOMAX[lap].a, gSigmaAverOfPrevNOMAX[lap].b, gSigmaAverOfPrevNOMAX[lap].c, (double)nPrev,
													  gSigmaAverOfPrevMAX[lap].sigma, sigmaNew );
			return ret;
		}
	}
}




double get_sigma( eLaplaceType_T lap, int nPrev, int bDoMax3x3, double sigmaNew )
{
	if(bDoMax3x3){
		if(nPrev==1){
			return gSigmaAverOfPrevMAX[lap].sigma;
		}else{
			double ret = sigma_func_aver_n( gSigmaAverOfPrevMAX[lap].a, gSigmaAverOfPrevMAX[lap].b, gSigmaAverOfPrevMAX[lap].c, (double)nPrev,
													  gSigmaAverOfPrevMAX[lap].sigma, sigmaNew );
			return ret;
		}
	}else{
		if(nPrev==1){
			return gSigmaAverOfPrevNOMAX[lap].sigma;
		}else{
			double ret = sigma_func_aver_n( gSigmaAverOfPrevNOMAX[lap].a, gSigmaAverOfPrevNOMAX[lap].b, gSigmaAverOfPrevNOMAX[lap].c, (double)nPrev,
													  gSigmaAverOfPrevMAX[lap].sigma, sigmaNew );
			return ret;
		}
	}
}

double get_significance( double lapNew, double lapPrev, double sigmaNew, double meanNew,
								 eLaplaceType_T lap, int nPrev, int bDoMax3x3, int bCoic )
{
	double sigmaPrev = get_sigma( lap, nPrev, bDoMax3x3, sigmaNew );
	double meanPrev  = gSigmaAverOfPrevMAX[lap].mean_aver;

	double probLapNew = get_lt_gauss_prob( meanNew, sigmaNew, lapNew );
	double probLapPrev = get_gt_gauss_prob( meanPrev, sigmaPrev, lapPrev);


	double ret = probLapNew*probLapPrev;
	//if(bCoic)
	//	ret = 2.00*ret;
	
	//ret = ( 1.00 - ret );
	return ret;
}

double get_gt_gauss_prob( double mean, double sigma, double val )
{
	double beta0 = ( (val-mean) / (SQRT_2*sigma) );
	double ret = 0.5*erfc(beta0);
	return ret;
}

double get_lt_gauss_prob( double mean, double sigma, double val )
{
	double ret = (1.00 - get_gt_gauss_prob( mean,sigma,val ));
	return ret;		
}



CLaplaceInfo::CLaplaceInfo()
: m_PlusCount(0), m_MinusCount(0), m_TypicalSigma(0), m_MinusFactor(0)
{

}

int InitGlobalLaplaceInfoTable()
{
	int lastVal = (int)eLastEnumElem;
	for(int i=0;i<lastVal;i++){
		eLaplaceType_T lap = (eLaplaceType_T)i;
		gLaplaceInfo[i].m_TypicalSigma = gSigmaBacgroundG[i];
		Table2D<ELEM_TYPE>::GetPlusMinusList( lap, 
														  gLaplaceInfo[i].m_PlusList, gLaplaceInfo[i].m_PlusCount,
														  gLaplaceInfo[i].m_MinusList, gLaplaceInfo[i].m_MinusCount );		
		if(gLaplaceInfo[i].m_MinusCount>0)
			gLaplaceInfo[i].m_MinusFactor = double(gLaplaceInfo[i].m_PlusCount)/double(gLaplaceInfo[i].m_MinusCount);
		else
			gLaplaceInfo[i].m_MinusFactor = 0;
	}
	_TRACE_PRINTF_6("Laplace info initialized - OK\n");
	for(int i=0;i<lastVal;i++){
		_TRACE_PRINTF_6("%s\t%.3f\n",Table2D<ELEM_TYPE>::GetLaplaceName((eLaplaceType_T)i),gLaplaceInfo[i].m_MinusFactor);
	}
	return lastVal;
}


