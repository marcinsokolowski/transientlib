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
	static double m_ReasonableMeanFit;
	static double m_ReasonableSigmaFit;

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
						  BOOL_T bCheckIfResonable=TRUE, BOOL_T bUseMeanOnFitFailed=FALSE,
						  BOOL_T bDumpFailedFit=TRUE );
	
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
