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
#ifndef _MY_HISTO_H__
#define _MY_HISTO_H__

#include <math.h>
#include <stdlib.h>
#include "mytypes.h"
#include "mystring.h"
#include "cexcp.h"

class CMyHisto
{
public :
	static BOOL_T m_DoRetryFit;

	double m_SkipBigger;	
	double m_MinValue;
	double m_MaxValue;
	int m_BinNo;
	double m_BinWidth;
	mystring m_szHistoName;

	double m_TotalSum;
	int m_AllCounts;

	double m_MininumValue;
	double m_MaximumValue;
	double m_MeanValue;
	double m_MaxOfDistr;	
	double m_RMSValue;
	double m_NormValue;
	BOOL_T m_bCharCalc;
	double m_Sum;
	double m_Sum2;
	
	double m_SigmaFit;
	double m_MeanFit;
	double m_NormFit;
	BOOL_T m_bPrevFitOK;		
	int m_ErrorCode;
	
	int* m_pCountTab;
	
	void GetStatValues( double& mean, double& rms, double& max_val );
	
	BOOL_T WasPrevFitOK(){ return m_bPrevFitOK; }
	double GetSigmaFit(){ return m_SigmaFit; }
	
	BOOL_T IsCharCalculated(){ return m_bCharCalc; }
	void SetCharCalculated(){ m_bCharCalc=TRUE; }
	
	void SetCharValues( double mean, double rms, double max_val );
	void SetRMS( double rms ){ m_RMSValue = rms; }

	void CalcValues( double typicalSigma );
	void CalcAllValues( double typicalSigma );	
	double GetMeanFromTotalSum();
	double GetMaxValue();

	inline BOOL_T IsInHisto(double value )
		{ 
			return (value>=m_MinValue && value<=m_MaxValue); 
		}				

	CMyHisto( const char* szName );	
	CMyHisto( const char* szName, double min_value, double max_value, int bin_no );

	void Fill(double value);
	inline void Fill2(double value){ Fill(value); }


	// return ok = Yes/No
	BOOL_T Fit( eFitType_T fittype, double* par, int par_count );

	BOOL_T CheckIfReasonableFit( double& mean_fit, double& sigma_fit, double& mean, double& rms );

	
	BOOL_T FitGauss( double& mean, double& rms, double& norm, int nSteps=1000, 
						  BOOL_T bCheckIfResonable=TRUE, BOOL_T bUseMeanOnFitFailed=FALSE );
	
	BOOL_T RetryFitGauss( double MeanStart, double SigmaStart, double NormStart,
	                      int nStepsSav, int bin_no, int _FirstBin=0, BOOL_T bWithLoop=TRUE );
	
	double GetBinCenter( int bin_no );		
	
	void Init();
	void Init( double min_value, double max_value, int bin_no );
	void Clear();
	~CMyHisto();

	int get_values( double* x_values, double* y_values, int max_count );
	void DumpToFile( const char* szFileName );
	mystring DumpToFile( int frame_index, const char* szInName="" );
	BOOL_T Read( const char* szFileName );
		
	const char* GetName() { return m_szHistoName.c_str(); }	
	
	inline int GetBinValue( int pos ){
		Assert( pos>=0 && pos<m_BinNo, "Outside of bin range : %d\n",pos);
		return m_pCountTab[pos];
	}
	int GetBinNo( double value );
	int GetMaxBin();
	double GetRMS( double& mean_val, int first_bin=0, int bin_no=0 );
};


#endif
