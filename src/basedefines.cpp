#include "basedefines.h"
#include <string.h>

ImageStat::ImageStat()
:Average(0),Min(0),Max(0),MinNonZero(0),Sum(0),m_pValDistrib(NULL),
 m_ValCount(0),MaxX(0),MaxY(0),RMS(0),MeanFit(0),SigmaFit(0),
 bFitOK(FALSE)
{
}


ImageStat::~ImageStat()
{
	if(m_pValDistrib)
		delete [] m_pValDistrib;
}



ImageStat::ImageStat( const ImageStat& right )
{
	if(m_pValDistrib)
      delete [] m_pValDistrib;
	if(right.m_pValDistrib && right.m_ValCount){
		m_ValCount = right.m_ValCount;
		m_pValDistrib = new long[m_ValCount];
		memcpy(m_pValDistrib,right.m_pValDistrib,sizeof(long)*m_ValCount);
	}
}

void ImageStat::AllocDistribTab( long maxValue )
{
	if(m_pValDistrib)
		delete [] m_pValDistrib;
	m_ValCount = maxValue;
	m_pValDistrib = new long[m_ValCount];
	memset( m_pValDistrib, 0, sizeof(long)*m_ValCount );
}

LONG_T ImageStat::GetMaxOccurValue()
{
	LONG_T MaxCounter=0;
	LONG_T MaxOccurVal=0;
	for(int i=0;i<m_ValCount;i++){
		if(m_pValDistrib[i]>MaxCounter){
			MaxCounter = m_pValDistrib[i];
			MaxOccurVal = i;		
		}				
	}
	return MaxOccurVal;
}


const char* GetFlipDesc( eDriverReverseImage_T fliptype )
{
	if( fliptype==eReverseImageFull ){
		return "FULL";
	}
	if( fliptype==eReverseImageHor ){
		return "FH";
	}
	if( fliptype==eReverseImageVert ){
		return "FV";
	}
	if( fliptype==eReverseImageHorVert ){
		return "FVH";
	}
	return "NONE";
}