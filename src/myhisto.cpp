#include "myhisto.h"
#include <string.h>
#include "cexcp.h"
#include "myfile.h"
#include <myfits.h>
#include <mathfunc.h>
#include "myparser.h"
#include "basestructs.h"
#include "mymacros.h"

BOOL_T CMyHisto::m_DoRetryFit=FALSE;
double CMyHisto::m_ReasonableMeanFit=10000.000;
double CMyHisto::m_ReasonableSigmaFit=10000.000;

void CMyHisto::Init()
{
   if(m_BinNo){
      if(!m_pCountTab){
         m_pCountTab = new int[m_BinNo];
      }
      Clear();
   }
}


void CMyHisto::Init( double min_value, double max_value, int bin_no )
{
	if(m_BinNo!=bin_no){
		delete [] m_pCountTab;
		m_pCountTab = NULL;		
	}	
	m_MinValue = min_value;
	m_MaxValue = max_value;
	m_BinNo = bin_no;
	m_SkipBigger = -1;

	m_MininumValue = 10000;
	m_MaximumValue = -10000;

	Assert(m_BinNo>0,"Number of bins must be greater then 0");
	if( m_MaxValue<=m_MinValue ){
		printf("error\n");
	}
//	Assert(m_MaxValue>m_MinValue,"Max value must be greater then minimal value max_val=%.2f, min_val=%.2f",
//				m_MaxValue,m_MinValue);
	
	if( m_MaxValue <= m_MinValue ){
		printf("ERROR ( ANALYSIS ) : Max value must be greater then minimal value max_val=%.2f, min_val=%.2f", 
		m_MaxValue,m_MinValue);
	}

	m_BinWidth = (m_MaxValue-m_MinValue)/m_BinNo;
	m_Sum = 0;
	m_Sum2 = 0;

	Init();
}


CMyHisto::CMyHisto( const char* szName, double min_value, double max_value, int bin_no )
: m_MinValue(min_value),m_MaxValue(max_value),m_pCountTab(NULL),m_BinNo(bin_no),
  m_MeanValue(0), m_MaxOfDistr(-10000), m_RMSValue(0), m_bCharCalc(FALSE),m_TotalSum(0),
  m_AllCounts(0), m_NormValue(0), m_bPrevFitOK(FALSE), m_SigmaFit(0), m_MeanFit(0), m_NormFit(0),
  m_ErrorCode(0),m_MininumValue(10000),m_MaximumValue(-10000),m_SkipBigger(-1)
{
	m_szHistoName = szName;

	Init( min_value,  max_value, bin_no );
}

CMyHisto::CMyHisto( const char* szName )
: m_MinValue(0),m_MaxValue(0),m_pCountTab(NULL),m_BinNo(0),
  m_MeanValue(0), m_MaxOfDistr(-10000), m_RMSValue(0), m_bCharCalc(FALSE),m_TotalSum(0),
  m_AllCounts(0), m_NormValue(0), m_bPrevFitOK(FALSE), m_SigmaFit(0), m_MeanFit(0), m_NormFit(0),
  m_ErrorCode(0), m_MininumValue(10000),m_MaximumValue(-10000)
{
}

double CMyHisto::GetMaxValue()
{ 
	double max_val = -10000.000;
	int max_pos = -1;
	for(register int i=0;i<m_BinNo;i++){
		if(m_pCountTab[i]>max_val){
			max_val = m_pCountTab[i];
			max_pos = i;
		}
	}
	m_MaxOfDistr = max_val;
	return m_MaxOfDistr;
}

/*BOOL_T CMyHisto::IsInHisto(double value )
{
	return (value>=m_MinValue && value<=m_MaxValue);
}*/

void CMyHisto::Clear()
{
	if(m_pCountTab){
		memset(m_pCountTab,0,sizeof(int)*m_BinNo);
	}
	m_MeanValue = 0;
	m_MaxOfDistr = -10000;
	m_RMSValue = 0;
	m_bCharCalc	= FALSE;
	m_TotalSum = 0;
	m_AllCounts = 0;
	m_NormValue = 0;
	m_Sum=0;
	m_Sum2=0;
}


CMyHisto::~CMyHisto()
{
	if(m_pCountTab)
		delete [] m_pCountTab;
}

double CMyHisto::GetMeanFromTotalSum()
{
	m_MeanValue = 0.00;
	if(m_AllCounts>0){
		m_MeanValue = m_TotalSum/m_AllCounts;
	}
	return m_MeanValue;
}

int CMyHisto::GetBinNo( double value )
{
	double from_min = (value-m_MinValue);
   int bin_no = (int)(from_min/m_BinWidth);
	return bin_no;
}

/*void CMyHisto::Fill2(double value)
{
	if(value>=m_MinValue && value<m_MaxValue){
		Fill( value );
		
		m_Sum = m_Sum + value;
		m_Sum2 = m_Sum2 + (value*value);
	}
}*/

void CMyHisto::Fill(double value)
{
	if( m_SkipBigger>0 ){
		if( value >= m_SkipBigger ){
			return;
		}
	}

	if(value>=m_MinValue && value<m_MaxValue){
		double from_min = (value-m_MinValue);
		int bin_no = (int)(from_min/m_BinWidth);
		// Assert(bin_no>=0 && bin_no<m_BinNo,"Bin for value %f , %d outside of range %d-%d, histogram range %f-%f !!!",value,bin_no,0,m_BinNo,m_MinValue,m_MaxValue);
		if(bin_no<0 || bin_no>=m_BinNo){
			Assert(FALSE,"Bin for value %f , %d outside of range %d-%d, histogram range %f-%f !!!",value,bin_no,0,m_BinNo,m_MinValue,m_MaxValue);
		}
		m_pCountTab[bin_no]++;

		m_TotalSum += value;
		m_AllCounts++;
	
		m_Sum = m_Sum + value;
      m_Sum2 = m_Sum2 + (value*value);
	}
	
	if(value>m_MaximumValue){
		m_MaximumValue = value;
	}
	if(value<m_MininumValue){
		m_MininumValue = value;
	}
}

double CMyHisto::GetBinCenter( int bin_no )
{
	double ret = m_MinValue + bin_no*m_BinWidth + 0.5*m_BinWidth;
	return ret;
}

mystring CMyHisto::DumpToFile( int frame_index,  const char* szInName )
{
	char szIndex[16];
	sprintf(szIndex,"%.5d",frame_index);

	mystring szName;
   szName << "HISTO/";
	if(szInName && szInName[0])
		szName << szInName << "_";
	szName << "histo_" << m_szHistoName << "_" << szIndex << ".txt";	
	DumpToFile( szName.c_str() );

	return szName;
}

void CMyHisto::DumpToFile( const char* szFileName )
{
	MyOFile out( szFileName );
	mystring szBuf;
	szBuf << "HISTO_NAME=" << m_szHistoName << "\n";
	szBuf << "MIN_VALUE=" << m_MinValue << "\n";
	szBuf << "MAX_VALUE=" << m_MaxValue << "\n";
	szBuf << "BIN_NO=" << m_BinNo << "\n";
	szBuf << "RMS=" << m_RMSValue << "\n";
	szBuf << "MEAN=" << m_MeanValue << "\n";
	szBuf << "MAX_OF_DISTR=" << m_NormValue << "\n";
	szBuf << "FIT_OK=" << m_bPrevFitOK << "\n";
	szBuf << "ALL_HITS=" << m_AllCounts << "\n";
	if(m_bPrevFitOK){
		szBuf << "MEAN_FIT=" << m_MeanFit << "\n";
		szBuf << "SIGMA_FIT=" << m_SigmaFit << "\n";
		szBuf << "NORM_FIT=" << m_NormFit << "\n";
	}
	
	char szBufTmp[128];
	for(register int i=0;i<m_BinNo;i++){
		sprintf(szBufTmp,"%.4f\t%d\n",GetBinCenter( i ),m_pCountTab[i]);
		szBuf << szBufTmp;
		// szBuf << GetBinCenter( i ) << "\t" << m_pCountTab[i] << "\n"; 	 	
	}
	out.Printf("%s",szBuf.c_str());
	out.Close();
	
	printf("DEBUG : CMyHisto::DumpToFile dumped histogram to file %s\n",szFileName);
}

int CMyHisto::get_values( double* x_values, double* y_values, int max_count )
{
	if( m_BinNo>max_count ){
		return 0;
	}
	for(register int i=0;i<m_BinNo;i++){
		x_values[i] = GetBinCenter( i );
		y_values[i] = m_pCountTab[i];
	}
	return m_BinNo;
}

BOOL_T CMyHisto::Read( const char* szFileName )
{
	if(!MyFile::DoesFileExist( szFileName ))
		return FALSE;
	

	MyIFile in( szFileName );
	
	const char* pLine;
	BOOL_T bOK=FALSE;
	while( pLine = in.GetLine(TRUE) ){
		MyParser pars = pLine;
		
		CEnvVar tmp;
      if(strstr(pars.c_str(),"=") && pars.GetVarAndValueFromLine(tmp.szName,tmp.szValue))
		{
			if(tmp.szName.Strcmp("HISTO_NAME")==0){
				m_szHistoName = tmp.szValue;
			}
			if(tmp.szName.Strcmp("MIN_VALUE")==0){
				m_MinValue = atof( tmp.szValue.c_str() );
			}
			if(tmp.szName.Strcmp("MAX_VALUE")==0){
				m_MaxValue = atof( tmp.szValue.c_str() );
			}
			if(tmp.szName.Strcmp("BIN_NO")==0){
				m_BinNo = atol( tmp.szValue.c_str() );
				m_BinWidth = (m_MaxValue-m_MinValue)/m_BinNo;
				Init();
			}
			if(tmp.szName.Strcmp("RMS")==0){
				m_RMSValue = atof( tmp.szValue.c_str() );
			}
			if(tmp.szName.Strcmp("MEAN")==0){
				m_MeanValue = atof( tmp.szValue.c_str() );
			}
			if(tmp.szName.Strcmp("MAX_OF_DISTR")==0){
				m_NormValue = atof( tmp.szValue.c_str() );
			}
			if(tmp.szName.Strcmp("FIT_OK")==0){
				m_bPrevFitOK = (atol(tmp.szValue.c_str())>0);
			}
			if(tmp.szName.Strcmp("ALL_HITS")==0){
				m_AllCounts = atol( tmp.szValue.c_str() );
			}
		}else{
			pars.Reset();
			mystring szBin,szBinVal;
			szBin = pars.GetNextItem();
			szBinVal = pars.GetNextItem();
			int bin_no = GetBinNo( atof(szBin.c_str()) );
			m_pCountTab[ bin_no ] = atol( szBinVal.c_str() );
			bOK = TRUE;
		}
	}
	return bOK;
}

BOOL_T CMyHisto::Fit( eFitType_T fittype, double* par, int par_count )
{
	return FALSE;
}

BOOL_T CMyHisto::CheckIfReasonableFit( double& mean_fit, double& sigma_fit,
													double& mean, double& rms )
{
	if( fabs(mean_fit) > CMyHisto::m_ReasonableMeanFit || fabs(sigma_fit) > CMyHisto::m_ReasonableSigmaFit )
		return FALSE;

	return TRUE;
}

BOOL_T CMyHisto::RetryFitGauss( double MeanStart, double SigmaStart, double NormStart, 
										  int nStepsSav, int bin_no, int _FirstBin, BOOL_T bWithLoop/*=TRUE*/ )
{
	m_SigmaFit = SigmaStart;
	m_MeanFit  = MeanStart;
	m_NormFit = NormStart;
	int nSteps = nStepsSav;
	printf("FIT correction first attempt : %.6f, %.6f, %.6f\n",m_SigmaFit,m_MeanFit,m_NormFit);
	BOOL_T bFit = CMyFit::FitGauss( m_MinValue, m_SigmaFit, m_MeanFit, m_NormFit,
									 bin_no, m_BinWidth, m_pCountTab, nSteps, _FirstBin );

	if(bFit){
		if(!CheckIfReasonableFit( m_MeanFit, m_SigmaFit, m_MeanValue, m_RMSValue ))
		{
			bFit = FALSE;
			printf("RetryFitGauss : fit no sens (%.2f, %.2f, %.2f)\n",m_SigmaFit, m_MeanFit,m_NormFit);
		}
	}
		
	if(!bFit && bWithLoop){
		double sigma_start = (int)(SigmaStart*0.5);
		double sigma_end = (int)(SigmaStart*1.5);
		double step = 5.00;
		for(double sigma=sigma_start;sigma<sigma_end;sigma+=step){
			double sigma_fit = sigma;
			m_MeanFit  = MeanStart;
			m_NormFit = NormStart;				
			nSteps = nStepsSav;
			_TRACE_PRINTF_3("starting from : %.6f, %.6f, %.6f\n",sigma,m_MeanFit,m_NormFit);

			
			bFit = CMyFit::FitGauss( m_MinValue, sigma_fit, m_MeanFit, m_NormFit,
									 bin_no, m_BinWidth, m_pCountTab, nSteps, _FirstBin );
			if(bFit){
				printf("OK - after correction : %.6f, %.6f, %.6f\n",sigma,m_MeanFit,m_NormFit);
				m_SigmaFit = sigma_fit;
				break;
			}else{
				printf("error = %d\n",nSteps);
			}
		}
		if(!CheckIfReasonableFit( m_MeanFit, m_SigmaFit, m_MeanValue, m_RMSValue ))
		{
			bFit = FALSE;
			printf("RetryFitGauss - 2 : fit no sens (%.2f, %.2f, %.2f)\n",m_SigmaFit, m_MeanFit,m_NormFit);
		}
	}
	return bFit;
}
										  

BOOL_T CMyHisto::FitGauss( double& mean, double& rms, double& norm, int nSteps, 
									BOOL_T bCheckIfResonable, BOOL_T bUseMeanOnFitFailed /* =FALSE */,
									BOOL_T bDumpFailedFit /*=TRUE*/ )
{
	int nStepsSav = nSteps;
	double SigmaStart = CMyMathFunc::round( m_RMSValue , 0 );
	double MeanStart = CMyMathFunc::round( m_MeanValue , 0 );
	double NormStart = CMyMathFunc::round( m_NormValue , 0 );

	// saving mean,rms and norm;
	double mean_sav = mean;
	double rms_sav = rms;
	double norm_sav = norm;

	m_bPrevFitOK = FALSE;

	m_SigmaFit = SigmaStart;
	m_MeanFit  = MeanStart;
	m_NormFit = NormStart;


	_TRACE_PRINTF_3("starting from : %.6f, %.6f, %.6f\n",m_SigmaFit,m_MeanFit,m_NormFit);


	BOOL_T bFit = CMyFit::FitGauss( m_MinValue, m_SigmaFit, m_MeanFit, m_NormFit,
				   			 		    m_BinNo, m_BinWidth, m_pCountTab, nSteps );


	if(!bFit && CMyHisto::m_DoRetryFit){
		int max_bin = GetMaxBin();
		int first_bin = MAX(0, max_bin-10);
		SigmaStart = GetRMS( MeanStart, first_bin, 20 );
		

		bFit = RetryFitGauss( MeanStart, SigmaStart, NormStart, nStepsSav, m_BinNo );

		if(!bFit){
			printf("trying re-bining ...\n");
			

			int re_bin = (m_BinNo/2);

			if(re_bin<5)
				return FALSE;

			int re_bin_2 = re_bin*2;
			CMyHisto tmp("test",m_MinValue,m_MaxValue,m_BinNo/2);
			for(int j=0;j<m_BinNo;j++){
				int pos = (j/2);
				if(pos<re_bin){
					tmp.m_pCountTab[pos] += m_pCountTab[j];
				}
			}
			tmp.CalcAllValues( SigmaStart );
			nSteps = nStepsSav;
			bFit = tmp.FitGauss( tmp.m_MeanFit, tmp.m_SigmaFit, tmp.m_NormFit, nSteps, tmp.m_BinNo );
			tmp.DumpToFile("tmp_rebin.txt");
			printf("re-bining result : %d\n",bFit);
			if(bFit){
				m_SigmaFit = tmp.m_SigmaFit;
				m_MeanFit = tmp.m_MeanFit;
				m_NormFit = tmp.m_NormFit/2;
			}
		}

		if(!bFit){
			first_bin = MAX(0, max_bin-15);
	
			printf("Step 0 : Last stage of re-try - changing fit range ...\n");
			bFit = RetryFitGauss( MeanStart, SigmaStart, NormStart, nStepsSav, 30, first_bin, FALSE );
			/*if(bFit){
				bFit = RetryFitGauss( m_MeanFit, m_SigmaFit, m_NormFit, nStepsSav, m_BinNo, 0, FALSE );
			}
		
			if(!bFit){
				printf("Step 1 : Last stage of re-try - changing fit range ...\n");
				first_bin = MAX(0, max_bin-5);
				bFit = RetryFitGauss( MeanStart, SigmaStart, NormStart, nStepsSav, 10, first_bin, FALSE );
				if(bFit){
					first_bin = MAX(0, max_bin-10);
					bFit = RetryFitGauss( m_MeanFit, m_SigmaFit, m_NormFit, nStepsSav, 20, first_bin, FALSE );
				}				
			}*/


			if(!bFit){
				printf("Step 2 : Last stage of re-try - changing fit range ...\n");
				double sigma_start = SigmaStart;
				double step_norm=(NormStart/10.00);
				double max_norm=(NormStart*2.00);


				for(double norm=(NormStart-step_norm);norm>=step_norm;norm-=step_norm){
					bFit = RetryFitGauss( MeanStart, sigma_start, norm, nStepsSav, 20, first_bin,  FALSE );
					if(bFit){					
						break;
					}	
				}

				if( !bFit ){
					for(double norm=(NormStart+step_norm);norm<=max_norm;norm+=step_norm){
						bFit = RetryFitGauss( MeanStart, sigma_start, norm, nStepsSav, 20, first_bin, FALSE );
						if(bFit){					
							break;
						}	
					}
				}

				if(!bFit){
					printf("Step 3 : Last stage of re-try - changing fit range ...\n");
					for(double norm=(NormStart-step_norm);norm>=step_norm;norm-=step_norm){
						bFit = RetryFitGauss( MeanStart, sigma_start, norm, nStepsSav, 20, first_bin );
						if(bFit){					
							break;
						}	
					}

					for(double norm=(NormStart+step_norm);norm<=max_norm;norm+=step_norm){
						bFit = RetryFitGauss( MeanStart, sigma_start, norm, nStepsSav, 20, first_bin );
						if(bFit){					
							break;
						}	
					}
				}
			}
		}

	}
	

	if(bFit){
		mean = m_MeanFit;
		rms = fabs(m_SigmaFit);
		norm = m_NormFit;
		m_SigmaFit = fabs(m_SigmaFit);

		if( m_SigmaFit<0 && rms<0 ){
			DumpToFile( 0, "histo_fit_error.txt" );
			printf("Error when fitting gauss !\n");fflush(0);
			return FALSE;
		}
		m_ErrorCode = nSteps;
	
		_TRACE_PRINTF_3("Fitted : mean=%.2f sigma=%.2f norm=%.2f\n",m_MeanFit,rms,norm);
		if(bCheckIfResonable){
			if(!CheckIfReasonableFit( mean, rms, m_MeanValue, m_RMSValue ))
			{
				if( bDumpFailedFit ){
					DumpToFile( 0, "failed_fit" );
				}
				_TRACE_PRINTF_1("FIT result does not make any sens - fit FAILED\n");
				return FALSE;			
			}
		}

		m_bPrevFitOK = TRUE;
		return TRUE;
	}else{
		_TRACE_PRINTF_1("\n\n/**********************************/\n");
		_TRACE_PRINTF_1("Fit of Gauss FAILED !!!!!!!!!!!!!!\n");
		_TRACE_PRINTF_1("/**********************************/\n\n");
		if( bDumpFailedFit ){
			DumpToFile( 0, "failed_fit2" );
		}

		if( bUseMeanOnFitFailed ){
			mean = mean_sav;
	      rms = fabs(rms_sav);			
			norm = norm_sav;
			printf("Fit failed, returing mean=%.2f, rms=%.2f, norm=%.2f\n",mean,rms,norm);
		}
	}	
	m_ErrorCode = nSteps;
	return FALSE;						 	
}

void CMyHisto::SetCharValues( double mean, double rms, double max_val )
{
	m_MeanValue = mean;
	m_RMSValue = rms;
	m_NormValue = m_MaxOfDistr = max_val;	
	m_bCharCalc = TRUE;	
}

void CMyHisto::CalcValues( double typicalSigma )
{
	double mean_sum=0;
	double all_count=0;
	for(register int i=0;i<m_BinNo;i++){
		mean_sum += GetBinCenter( i )*m_pCountTab[i];
		all_count += m_pCountTab[i];
	}		
	m_MeanValue = (mean_sum/all_count);

	double disp=0;
	for(register int i=0;i<m_BinNo;i++){
		disp += m_pCountTab[i]*CMyMathFunc::mysqr(GetBinCenter( i ) - m_MeanValue);
	}
	m_RMSValue = sqrt(disp/all_count);
}

void CMyHisto::CalcAllValues( double typicalSigma )
{
	CalcValues( typicalSigma );
	GetMaxValue();
	SetCharValues( m_MeanValue, m_RMSValue, m_MaxOfDistr );
}


double CMyHisto::GetRMS( double& mean_val, int first_bin, int bin_no )
{
	if(!bin_no)
		bin_no = m_BinNo;

	double mean_sum=0;
	double all_count=0;

	int last_bin=MIN( (first_bin+bin_no) , m_BinNo );
	for(register int i=first_bin;i<last_bin;i++){
		mean_sum += GetBinCenter( i )*m_pCountTab[i];
		all_count += m_pCountTab[i];
	}		
	mean_val = (mean_sum/all_count);


	double disp=0;
	for(register int i=first_bin;i<last_bin;i++){
		disp += m_pCountTab[i]*CMyMathFunc::mysqr(GetBinCenter( i ) - mean_val);		
	}
	double ret = sqrt(disp/all_count);
	return ret;
}



int CMyHisto::GetMaxBin()
{
	int ret=0;
	int max_count=m_pCountTab[0];
	for(register int i=1;i<m_BinNo;i++){
		if( m_pCountTab[i] > max_count )
		{
			max_count = m_pCountTab[i];
			ret = i;
		}
	}
	return ret;
}

void CMyHisto::GetStatValues( double& mean, double& rms, double& max_val )
{
	mean = 0;
	rms = 0;
	max_val = 0;
	if( m_AllCounts>0 ){
		mean = (m_Sum/m_AllCounts);
		rms = sqrt( m_Sum2/m_AllCounts - mean*mean );
		max_val = GetMaxValue();
		SetCharValues( mean, rms, max_val );
	}
}