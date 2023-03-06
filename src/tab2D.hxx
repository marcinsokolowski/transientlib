#ifndef _TAB2D_HXX_
#define _TAB2D_HXX_


#include <math.h>
#include "tab2D.h"
#include "cexcp.h"
#include <errno.h>
#include "random.h"
#include "mydate.h"
#include "tab2Ddesc.h"
#include "myutil.h" 
#include "mypoints.h"
#include "mymacros.h"
#include "mycmnglobals.h"
#include "baseanal.h"
#include "myhisto.h"
#include "mytrace.h"
#include "basestructs.h"
#include "laplace_info.h"
#include "ccdrowcollist.h"
#include "mathconst.h"
#include "myprogress.h"
#include "mymacros.h"

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::SetDataPtr( ARG_TYPE* pData )
{
	Assert( m_bAllocHere==FALSE, "Only outside allocation tables can use this function");
	m_pData = pData;
	ReCreateIndex();
}

template<class ARG_TYPE>
Table2D<ARG_TYPE>::Table2D()
: m_pData(NULL),m_Index(0),m_pFastData(NULL), m_X_On_Big(0),
  m_Y_On_Big(0),m_PrevShiftX(0),m_PrevShiftY(0),m_StartOnOrigX(0),m_StartOnOrigY(0),
  m_SizeX(0),m_SizeY(0),m_Size(0),m_MaxVal(0),m_MaxX(0),m_MaxY(0),
  m_totalShiftX(0),m_totalShiftY(0),m_bAllocHere(FALSE),
  m_pDistribTab(NULL),m_MaxValue(0),m_pValuesTable(NULL),m_pWrkTable(NULL),
  m_MedianValue(0),m_MeanValue(0)
{


}

template<class ARG_TYPE>
BOOL_T Table2D<ARG_TYPE>::InitConstructor(long x_size,long y_size,long idx,BOOL_T bAllocHere)
{
	BOOL_T bRet=TRUE;
	
	m_Index = idx;
   m_bAllocHere = bAllocHere;

	if(m_bAllocHere)
		bRet = Alloc(x_size,y_size);

	// lines below are MUST !!!
	// in case not bAllocHere , size of image must be properly set !
	m_SizeX = x_size;
	m_SizeY = y_size;
	m_Size = (m_SizeX*m_SizeY);

	return bRet;
}

template<class ARG_TYPE>
Table2D<ARG_TYPE>::Table2D(long x_size,long y_size,long idx,BOOL_T bAllocHere)
: m_pData(NULL),m_Index(idx),m_pFastData(NULL), m_X_On_Big(0),
  m_Y_On_Big(0),m_PrevShiftX(0),m_PrevShiftY(0),m_StartOnOrigX(0),m_StartOnOrigY(0),
  m_SizeX(0),m_SizeY(0),m_Size(0),m_MaxVal(0),m_MaxX(0),m_MaxY(0),
  m_totalShiftX(0),m_totalShiftY(0),m_bAllocHere(bAllocHere),
  m_pDistribTab(NULL),m_MaxValue(0),m_pValuesTable(NULL),m_pWrkTable(NULL),
  m_MedianValue(0),m_MeanValue(0)
{
	if(m_bAllocHere)
		Alloc(x_size,y_size);
	m_SizeX = x_size;
	m_SizeY = y_size;
	m_Size = (m_SizeX*m_SizeY);
}

template<class ARG_TYPE>
Table2D<ARG_TYPE>::Table2D(const Table2D& right)
: m_pData(NULL),m_Index(0),m_pFastData(NULL), m_X_On_Big(0),
  m_Y_On_Big(0),m_PrevShiftX(0),m_PrevShiftY(0),m_StartOnOrigX(0),m_StartOnOrigY(0),
  m_SizeX(0),m_SizeY(0),m_Size(0),m_MaxVal(0),m_MaxX(0),m_MaxY(0),
  m_totalShiftX(0),m_totalShiftY(0),m_pDistribTab(NULL),m_MaxValue(0),
  m_pValuesTable(NULL),m_pWrkTable(NULL)
{
	(*this) = right;
}

template<class ARG_TYPE>
Table2D<ARG_TYPE>::~Table2D()
{
	// printf("Desctructor !!!! ????\n");
	if(m_bAllocHere){
		if(m_pData)
			delete [] m_pData;
	}
	if(m_pFastData)
		delete [] m_pFastData;

	if(m_pDistribTab)
		delete [] m_pDistribTab;

	if(m_pValuesTable)
		delete [] m_pValuesTable;

	if(m_pWrkTable)
		delete [] m_pWrkTable;
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::ClearState()
{
	m_PrevShiftX = 0;
	m_PrevShiftY = 0;
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::ReCalcIndex()
{
	long pos=0;
	for(register int y=0;y<m_SizeY;y++){
		m_pFastData[y] = (m_pData+pos);
		pos += m_SizeX;
	}
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::ReCreateIndex()
{
	if(m_pFastData){
		delete [] m_pFastData;
      m_pFastData = NULL;
   }	
	m_pFastData = new ARG_TYPE*[m_SizeY];
printf("DEBUG : ::ReCreateIndex m_SizeY=%d (ptr = %x)\n",m_SizeY,m_pFastData);
	ReCalcIndex();		
}

template<class ARG_TYPE>
BOOL_T Table2D<ARG_TYPE>::Alloc(long x_size,long y_size)
{
	// printf(" In Alloc ??? this=%x\n",this);
	if(m_SizeX==x_size && m_SizeY==y_size)
		return TRUE;

	if(m_bAllocHere){
		// can only delete if allocated HERE !!!
		if (m_pData)
			delete [] m_pData;
	}

	if(m_pFastData){
		delete [] m_pFastData;
		m_pFastData = NULL;
	}
	m_SizeX = x_size;
	m_SizeY = y_size;
	m_Size = (m_SizeX*m_SizeY);
	if (m_SizeX>0 && m_SizeY>0){
		 try{
			 if(m_bAllocHere){
				 long size_of_elem = sizeof(ARG_TYPE);
			    size_t nBytes = m_Size*size_of_elem;
			    // m_pData = (ARG_TYPE*)malloc(nBytes);
				 m_pData = (ARG_TYPE*)new char[nBytes];
				 ReCreateIndex();
				 if(m_pData)
	            Init();
			 }
		 }catch(...){
			 return FALSE;			
		 }
		 if (!m_pData)
			return FALSE;
		return TRUE;
	}
	return FALSE;
}

template<class ARG_TYPE>
BOOL_T Table2D<ARG_TYPE>::ReSize(long x_size,long y_size)
{
	if(m_SizeX==x_size && m_SizeY==y_size)
		return TRUE;

	if( x_size*y_size > m_SizeX*m_SizeY ){
		// new image is larger -> buffer size must be increased :
		return Alloc( x_size , y_size );
	}
	

	if(m_pFastData){
		delete [] m_pFastData;
		m_pFastData = NULL;
	}
	m_SizeX = x_size;
	m_SizeY = y_size;
	m_Size = (m_SizeX*m_SizeY);
	if (m_SizeX>0 && m_SizeY>0){
		 try{
			 if(m_bAllocHere){
				 ReCreateIndex();
				 if(m_pData)
	            Init();
			 }
		 }catch(...){
			 return FALSE;			
		 }
		 if (!m_pData)
			return FALSE;
		return TRUE;
	}
	return FALSE;
	
}

template<class ARG_TYPE>
ARG_TYPE Table2D<ARG_TYPE>::getval(long x,long y)
{
	ARG_TYPE ret;
	ret = get(x,y);
	return ret;
}

template<class ARG_TYPE>
ARG_TYPE& Table2D<ARG_TYPE>::get(long x,long y)
{
    AssertNULL(m_pData);
    Assert((x<m_SizeX && x>=0),"X size range (%d) exceeded, requested x=%d",m_SizeX,x);	 	
    Assert((y<m_SizeY && y>=0),"Y size range (%d) exceeded, requested y=%d",m_SizeY,y);	
	 return m_pData[y*m_SizeX+x]; 	
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::Init(ARG_TYPE val){
	if(m_bAllocHere){
		if(m_pData){
			memset(m_pData,'\0',sizeof(ARG_TYPE)*m_Size); 
		}
	}
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::SetData(ARG_TYPE val)
{
	if(m_pData){
		// memset(m_pData,(int)val,sizeof(ARG_TYPE)*m_Size); 		
		for(register int i=0;i<m_Size;i++){
			m_pData[i] = val;
		}
	}
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::InitRandom(long min_value,long max_value)
{
	if(m_pData){
		CRandom::Initialize(DEFAULT_SEED/*+get_dttm()*/);
		CRandom rnd;
		for(register int i=0;i<m_Size;i++){
			LONG_T r = rnd.GetRandomInteger(min_value,max_value);
			m_pData[i] = r;
		}		
	}	
}

// buffer must already be allocated !!!!
template<class ARG_TYPE>
void Table2D<ARG_TYPE>::copy( ARG_TYPE* src )
{
	memcpy(m_pData,src,m_SizeX*m_SizeY*sizeof(ARG_TYPE));
}

template<class ARG_TYPE>
Table2D<ARG_TYPE>& Table2D<ARG_TYPE>::Assign(const Table2D<ARG_TYPE>& right)
{
//printf("here::Assign 1\n %d,%d := %d,%d",m_SizeX,m_SizeY,right.m_SizeX,right.m_SizeY);fflush(0); 
	if(m_SizeX!=right.m_SizeX || m_SizeY!=right.m_SizeY)
		Alloc(right.m_SizeX,right.m_SizeY);
//printf("here::Assign 2\n");fflush(0);

	ARG_TYPE* tmp_ptr = m_pData; // needed for tests  

	if( m_pData != right.m_pData ){
		if( m_bAllocHere ){
			printf("DEBUG : Table2D<ARG_TYPE>::Assign : memcpy( 0x%x , 0x%x )\n",m_pData,right.m_pData);
			memcpy(m_pData,right.m_pData,m_SizeX*m_SizeY*sizeof(ARG_TYPE));
		}else{
			// MS : 2009-09-28 due to bug found that after GetData was done
			// the image 00875 was overwriten somehow by 00874-DARK 
			// it turned out that it is because of  CCDPipeline::RestartPipelineKeepCurrent
			// and (GetNext()) = m_pCCD[curr_pos];
			// which made a copy of previous image to memory already containing 
			// new image from GetData !!!
			// This is unacceptable !!!
			printf("DEBUG : Table2D<ARG_TYPE>::Assign m_pData(= 0x%x ) := SetDataPtr( 0x%x )\n",m_pData,right.m_pData);
			SetDataPtr( right.m_pData );
		}

		// TEMPORARY - to be commented out - corresponds to ccd_controller.cpp / GetData
		printf("OLD ( 0x%x ) : ",tmp_ptr);
		for(int ii=(2062*100+1000);(ii<=(2062*100+1005) && ii<m_Size);ii++){
			printf("%d ",tmp_ptr[ii]);
		}
		printf("\nNEW ( 0x%x ) : ",m_pData);
		for(int ii=(2062*100+1000);(ii<=(2062*100+1005) && ii<m_Size);ii++){
			printf("%d ",m_pData[ii]);
		}
		printf("\n");
	}
//printf("here::Assign 3\n");fflush(0);
	return (*this);
}

template<class ARG_TYPE>
Table2D<ARG_TYPE>& Table2D<ARG_TYPE>::OperatorEq(const Table2D<ARG_TYPE>& right)
{
	(*this) = right;
	return (*this);
}

template<class ARG_TYPE>
Table2D<ARG_TYPE>& Table2D<ARG_TYPE>::operator=(const Table2D<ARG_TYPE>& right)
{
//printf("here 1?\n");fflush(0);
	Assign( right );
//printf("here 2?\n");fflush(0);
	m_KeyTab = right.m_KeyTab;
//printf("here 3?\n");fflush(0);
	m_MaxX = right.m_MaxX;
	m_MaxY = right.m_MaxY;
	m_MaxVal = right.m_MaxVal;
	m_MedianValue = right.m_MedianValue;
	m_MeanValue = right.m_MeanValue;
	m_StartOnOrigX = right.m_StartOnOrigX;
   m_StartOnOrigY = right.m_StartOnOrigY;

//printf("here 4?\n");fflush(0);
	return (*this);
}

template<class ARG_TYPE>
ARG_TYPE* Table2D<ARG_TYPE>::GetPos(long x,long y)
{
	// Assert(x>=0 && x<m_SizeX,"GetPos : x range exceeded x=%d, m_Size=%d",x,m_SizeX);
	if(x<0 || x>=m_SizeX){
		Assert(FALSE,"GetPos : x range exceeded x=%d, m_Size=%d",x,m_SizeX);		
	}


	Assert(y>=0 && y<m_SizeY,"GetPos : y range exceeded");
 
	ARG_TYPE* pos = m_pData + y*m_SizeX + x;
	return pos;
}

template<class ARG_TYPE>
ARG_TYPE& Table2D<ARG_TYPE>::GetPixel(long x,long y)
{
	return (*GetPos(x,y));
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::GetXY_FromPos( LONG_T pos, LONG_T& x ,LONG_T& y )
{
	x = (pos % m_SizeX);
	y = (pos / m_SizeX);	
}

template<class ARG_TYPE>
LONG_T Table2D<ARG_TYPE>::GetPos_FromXY( LONG_T& x ,LONG_T& y )
{
	LONG_T pos = (y*m_SizeX + x);
	return pos;
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::CalcMaxAndPos()
{
	m_MaxVal = GetMaxValueAndPos( 0, 0, m_SizeX, m_SizeY, m_MaxX, m_MaxY );
}

template<class ARG_TYPE>
int Table2D<ARG_TYPE>::GetMaxValue( int& max_x, int& max_y )
{
	CalcMaxAndPos();
	max_x = m_MaxX;
	max_y = m_MaxY;
	return m_MaxVal;
}


template<class ARG_TYPE>
LONG_T Table2D<ARG_TYPE>::GetMaxValueAndPos( long start_x, long start_y,
                        					      long end_x, long end_y,
					                              long& max_x, long& max_y )
{
	long x0 = MAX(start_x,0);
	long y0 = MAX(start_y,0);
	long x1 = MIN(end_x,m_SizeX);
	long y1 = MIN(end_y,m_SizeY);

	long max=-10000;
	for(register long y=y0;y<y1;y++){
		for(register long x=x0;x<x1;x++){
			if(m_pFastData[y][x]>max){
				max_x = x;
				max_y = y;
				max = (long)(m_pFastData[y][x]);
			}
		}
	}
	return max;
}

template<class ARG_TYPE>
int Table2D<ARG_TYPE>::GetPixelsAbove( double treshold, CPointList& pointList, int start_x, int end_x, int start_y, int end_y )
{
	int minX = MAX(start_x,0);
	int minY = MAX(start_y,0);
	int maxX = MIN(end_x,m_SizeX);
	int maxY = MIN(end_y,m_SizeY);

	for(register int y=minY;y<maxY;y++){
		for(register int x=minX;x<maxX;x++){
			if(m_pFastData[y][x]>treshold){
				CPoint tmp;
				tmp.x = x;
				tmp.y = y;
				pointList.push_back( tmp );
			}
		}
	}
	return pointList.size();
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::GetSafeStatistics( LONG_T start_x, LONG_T start_y,
                               		       LONG_T end_x, LONG_T end_y,
														 double& mean, double& RMS )
{
	if( start_x<0 ){
		start_x = 0;
	}
	if( start_y<0 ){
		start_y = 0;
	}
	if( end_x>=m_SizeX ){
		end_x = m_SizeX - 1;
	}
	if( end_y>=m_SizeY ){
		end_y = m_SizeY - 1;
	}
	GetStatistics( start_x, start_y, end_x, end_y, mean, RMS );
}


template<class ARG_TYPE>
int Table2D<ARG_TYPE>::GetBackgroundValue( LONG_T start_x, LONG_T start_y, 
														 LONG_T end_x, LONG_T end_y,
														 LONG_T* buffer, int nSize )
{
	if( start_x<0 ){
		start_x = 0;
	}
	if( start_y<0 ){
		start_y = 0;
	}
	if( end_x>=m_SizeX ){
		end_x = m_SizeX - 1;
	}
	if( end_y>=m_SizeY ){
		end_y = m_SizeY - 1;
	}

	register int i=0;
	for(register int y=start_y;y<=end_y;y++){
		for(register int x=start_x;x<=end_x;x++){
			if(i<nSize){
				buffer[i] = m_pFastData[y][x];
				i++;
			}else{
				printf("BUFFER TO SMALL IN Table2D<ARG_TYPE>::GetBackgroundValue\n");
				break;
			}
		}
	}
	my_qsort( buffer, i );
	int quarter_pos = i/4;
	int end_pos = (quarter_pos + (i*0.1))*1.2;
	if(end_pos<(i/2)){
		end_pos = i/2;
	}
	int pos = (end_pos/2);

	return buffer[pos];
}


template<class ARG_TYPE>
void Table2D<ARG_TYPE>::GetStatistics( LONG_T start_x, LONG_T start_y,
                                       LONG_T end_x, LONG_T end_y,
													double& mean, double& RMS )
{
	ImageStat statInfo;
	GetStatistics( statInfo, FALSE, start_x,  start_y, end_x, end_y );
	mean = statInfo.Average;
	RMS = statInfo.RMS;
//	mean_fit = statInfo.MeanFit;
//	sigma_fit = statInfo.SigmaFit
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::GetSafeStatistics( ImageStat& statInfo ,
											  			 BOOL_T bDistrib,
                                       	 LONG_T start_x, LONG_T start_y,
                                       	 LONG_T end_x, LONG_T end_y,
														 BOOL_T bDoFit /*=FALSE*/ )
{
	if( start_x<0 ){
		start_x = 0;
	}
	if( start_y<0 ){
		start_y = 0;
	}
	if( end_x>=m_SizeX ){
		end_x = m_SizeX - 1;
	}
	if( end_y>=m_SizeY ){
		end_y = m_SizeY - 1;
	}
	GetStatistics( statInfo, bDistrib, start_x, start_y, end_x, end_y, bDoFit );
}


template<class ARG_TYPE>
void Table2D<ARG_TYPE>::GetStatistics( ImageStat& statInfo ,
													BOOL_T bDistrib,
                                       LONG_T start_x, LONG_T start_y,
                                       LONG_T end_x, LONG_T end_y,
													BOOL_T bDoFit /*=FALSE*/,
													double min_value, double max_value )
{
	LONG_T size = (m_SizeX*m_SizeY);

	// do not change : this are defaults - in case 
	// nothing passed whole frame is used !!!
	if(end_x<0)
		end_x = m_SizeX;
	if(end_y<0)
		end_y =m_SizeY;
	// do not change 


	statInfo.Sum=0;
	statInfo.Max=-10000;
	statInfo.Min=66000;
	statInfo.MinNonZero=20000;	
	statInfo.Average=0;	
	statInfo.RMS=0;
	double sum2=0;
	int count=0;

	for( register int y=start_y;y<=end_y;y++){
		for( register int x=start_x;x<=end_x;x++){
			if( x>=0 && y>=0 && x<m_SizeX && y<m_SizeY ){
				int val = m_pFastData[y][x];

				if( val >= min_value && val <= max_value ){
					statInfo.Sum += val;
					if(val<statInfo.Min){
						statInfo.Min = val;
					}
					if(val<statInfo.MinNonZero && val>0){
						statInfo.MinNonZero = val;
					}		
					if(val>statInfo.Max){
						statInfo.Max = val;		
						statInfo.MaxX = x;
						statInfo.MaxY = y;
					}
					sum2 += (double(val)*double(val));
					count++;
				}
			}
		}
	}
	statInfo.Average = (statInfo.Sum/count);	
	double mean2 = (sum2/count);
	double avg2 = double(statInfo.Average)*double(statInfo.Average);
	if( mean2 >= avg2 ){
		statInfo.RMS = sqrt( mean2 - avg2 );
	}else{
		printf("ERROR in Table2D<ARG_TYPE>::GetStatistics mean2<avg2 (%.2f < %.2f) \n",mean2,avg2);
		statInfo.RMS = 0.00;		
	}
	if(bDistrib)
		GetDistrib( statInfo, (LONG_T)statInfo.Max );
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::GetPartStat( ImageStat& statInfo,
					                      BOOL_T bDistrib,
					                      LONG_T start_x, LONG_T start_y,
                					       LONG_T end_x, LONG_T end_y )
{
	
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::GetPartDistrib( LONG_T* pDistribTab, LONG_T maxValue,
					                         LONG_T start_x, LONG_T start_y,
               					          LONG_T end_x, LONG_T end_y )
{
	if(end_x<0 || end_x>m_SizeX)
		end_x = m_SizeX;
	if(end_y<0 || end_y>m_SizeY)
		end_y = m_SizeY;
	if(start_x<0)
		start_x = 0;
	if(start_y<0)
		start_y = 0;
	
	// printf("looking in (%d,%d)-(%d,%d)\n",start_x,start_y,end_x,end_y);
	InitDistribTab( maxValue );
	for(register int y=start_y;y<end_y;y++){
      for(register int x=start_x;x<end_x;x++){
			// printf("%d ",m_pFastData[y][x]);
			if(m_pFastData[y][x]<maxValue && m_pFastData[y][x]>=0){
         	pDistribTab[ (long) m_pFastData[y][x] ]++;
	      }
		}
	}
	// printf("\n");			
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::ClearDistribTab()
{	
	if(m_pDistribTab)
		memset((void*)m_pDistribTab,0,sizeof(LONG_T)*m_MaxValue);
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::InitDistribTab( LONG_T MaxValue )
{
	if(m_MaxValue!=MaxValue || !m_pDistribTab){
		if(m_pDistribTab)
			delete [] m_pDistribTab;
		m_MaxValue = MaxValue;
		m_pDistribTab = new LONG_T[m_MaxValue];
	}
	ClearDistribTab();
}

template<class ARG_TYPE>
LONG_T Table2D<ARG_TYPE>::GetMostPopular( LONG_T MaxValue,
														LONG_T start_x, LONG_T start_y,
		                                    LONG_T end_x, LONG_T end_y )
{
	LONG_T MostPopular=0;
   LONG_T MaxCount=0;
	InitDistribTab( MaxValue );
	GetPartDistrib( m_pDistribTab, MaxValue, start_x, start_y, end_x, end_y );
	for(register int i=0;i<MaxValue;i++){
		// printf("%d ",m_pDistribTab[i]);
		if(m_pDistribTab[i]>MaxCount){
			MaxCount = m_pDistribTab[i];
			MostPopular = i;
		}
	}
	// printf("\n");
	return MostPopular;
}

template<class ARG_TYPE>
BOOL_T Table2D<ARG_TYPE>::CalcStat( ARG_TYPE* data, int sizeX, int sizeY,
                          				CWindowList& rowcollist, double& rms,
			                           double& mean, BOOL_T bFitGauss,
			                           int frame_index )
{
	ARG_TYPE** pFastData = new ARG_TYPE*[sizeY];
	long pos=0;
   for(register int y=0;y<sizeY;y++){
      pFastData[y] = (data+pos);
      pos += sizeX;
   }
	
	BOOL_T bRet = CalcStat( pFastData, sizeX, sizeY, rowcollist, rms, 
									mean, bFitGauss, frame_index );
	delete [] pFastData;

	return bRet;
}


template<class ARG_TYPE>
BOOL_T Table2D<ARG_TYPE>::CalcStat( ARG_TYPE** data, int sizeX, int sizeY,
                          				CWindowList& rowcollist, double& rms,
			                           double& mean, BOOL_T bFitGauss,
			                           int frame_index )
{
	double sum=0,sum2=0;

	int count=0;
	double min_value=100000;
	double max_value=-100000;
	for(register int i=0;i!=rowcollist.size();i++){
		CCDWindow& desc = rowcollist[i];
		
		for(register int y=desc.m_LowY;y<=desc.m_UpY;y++){
			for(register int x=desc.m_LowX;x<=desc.m_UpX;x++){
				if(x>=0 && y>=0 && x<sizeX && y<sizeY){
					double prev_sum=sum,prev_sum2=sum2;
					double newval = ((double)(data[y][x]));
					double newval2 = newval*newval;

					sum += newval;
					sum2 += newval2;

					if( sum < prev_sum ){
						printf("ERROR : in Table2D<ARG_TYPE>::CalcStat : sum=%.8f < prev_sum=%.8f after adding %d\n",sum,prev_sum,(int)data[y][x]);
					}
					if( sum2 < prev_sum2 ){
						printf("ERROR : in Table2D<ARG_TYPE>::CalcStat : sum2=%.8f < prev_sum2=%.8f after adding %d^2 = %.8f\n",sum2,prev_sum2,(int)data[y][x],newval2);
					}

					count++;
					if(data[y][x]<min_value){
						min_value = data[y][x];
					}
					if(data[y][x]>max_value){
						max_value = data[y][x];
					}				
				}
			}
		}
	}		
	mean = 0.00;
	rms = 0.00;

	if(count>0){
		mean = (sum/count);
		double mean2 = (sum2/count);
		if( mean2 >= (mean*mean) ){
			rms = sqrt( (mean2) - (mean*mean) );
		}else{
			rms = 0.00;
			printf("WARNING : rms calculation failed ! sum=%.8f sum2=%.8f count=%d mean2=%.8f < mean^2=%.8f^2=%.8f !!!\n",
						sum,sum2,count,mean2,mean,(mean*mean));fflush(stdout);
		}
	}

	if( rowcollist.size() ){
		CCDWindow& desc = rowcollist[0];		

		printf("WINDOW = (%d,%d)-(%d,%d) min=%.2f max=%.2f mean=%.2f rms=%.2f\n",
				desc.m_LowX,desc.m_LowY,
				desc.m_UpX,desc.m_UpY,min_value,max_value,mean,rms);
	}


	BOOL_T bRet=FALSE;
	if(bFitGauss && count>0){
		// in case gauss fit required , fill histo and fits gauss :
		CMyHisto Histo( "calc_stat", min_value-(5+rms), max_value+(5+rms),  50 );
		for(register int i=0;i!=rowcollist.size();i++){
			CCDWindow& desc = rowcollist[i];
		
			for(register int y=desc.m_LowY;y<=desc.m_UpY;y++){
				for(register int x=desc.m_LowX;x<=desc.m_UpX;x++){
					if(x>=0 && y>=0 && x<sizeX && y<sizeY){
						Histo.Fill( data[y][x] );
					}
				}
			}
		}
		int nSteps=500;
		double mean_fit=mean,sigma_fit=rms,norm_fit=Histo.GetMaxValue();
		Histo.SetCharValues( mean_fit, sigma_fit, norm_fit );
		if(Histo.FitGauss( mean_fit,sigma_fit,norm_fit, nSteps, FALSE )){
			// successfull fit :
			rms = sigma_fit;
			mean = mean_fit;			
			bRet = TRUE;
		}else{
			printf("WARNING : fit in CalcStat FAILED mean=%.2f rms=%.2f\n",mean,rms);
		}
		if(gDoDumpAllHisto || !bRet){
         Histo.DumpToFile( frame_index );
      }
	}

	return bRet;		

}

template<class ARG_TYPE>
BOOL_T Table2D<ARG_TYPE>::CalcStat( CWindowList& rowcollist, double& rms,
                    double& mean, BOOL_T bFitGauss, int frame_index )
{
	return CalcStat( m_pFastData, m_SizeX, m_SizeY, rowcollist, rms,
						  mean, bFitGauss, frame_index );
}

template<class ARG_TYPE>
BOOL_T Table2D<ARG_TYPE>::CalcStat( CRowColList& rowcollist, double& rms, 
												double& mean, BOOL_T bFitGauss, int frame_index )
{
	double sum=0,sum2=0;

	int count=0;
	double min_value=100000;
	double max_value=-100000;
	for(register int i=0;i!=rowcollist.size();i++){
		CRowColDesc& desc = rowcollist[i];
		if(desc.type==eCCDRow){
			// row :
			if(desc.num>=0 && desc.num<m_SizeY){
				for(register int x=0;x<m_SizeX;x++){
					double newval = ((double)(m_pFastData[desc.num][x]));
					double newval2 = newval*newval;

					sum += newval;
					sum2 += newval2;
					count++;
					if(m_pFastData[desc.num][x]<min_value){
						min_value = m_pFastData[desc.num][x];
					}
					if(m_pFastData[desc.num][x]>max_value){
						max_value = m_pFastData[desc.num][x];
					}
				}
			}
		}		
		if(desc.type==eCCDCol){
			// column :
			if(desc.num>=0 && desc.num<m_SizeX){
				for(register int y=0;y<m_SizeY;y++){
					double newval = ((double)(m_pFastData[y][desc.num]));
					double newval2 = newval*newval;

					sum += newval;
					sum2 += newval2;
					count++;
					if(m_pFastData[y][desc.num]<min_value){
						min_value = m_pFastData[y][desc.num];
					}
					if(m_pFastData[y][desc.num]>max_value){
						max_value = m_pFastData[y][desc.num];
					}
				}
			}
		}
	}		

	mean = 0.00;
	rms = 0.00;

	if(count>0){
		mean = (sum/count);
		double mean2 = (sum2/count);
		rms = sqrt( (mean2) - (mean*mean) );
	}

	if(bFitGauss && count>0){
		// in case gauss fit required , fill histo and fits gauss :
		CMyHisto Histo( "calc_stat", min_value, max_value,  100 );
		for(register int i=0;i!=rowcollist.size();i++){
			CRowColDesc& desc = rowcollist[i];
			if(desc.type==eCCDRow){
				// row :
				if(desc.num>=0 && desc.num<m_SizeY){
					for(register int x=0;x<m_SizeX;x++){
						Histo.Fill( m_pFastData[desc.num][x] );
					}
				}		
			}		
			if(desc.type==eCCDCol){
				// column :
				if(desc.num>=0 && desc.num<m_SizeX){
					for(register int y=0;y<m_SizeY;y++){
						Histo.Fill( m_pFastData[y][desc.num] );
					}
				}
			}
		}
		int nSteps=500;
		double mean_fit=mean,sigma_fit=rms,norm_fit=Histo.GetMaxValue();
		Histo.SetCharValues( mean_fit, sigma_fit, norm_fit );
		if(Histo.FitGauss( mean_fit,sigma_fit,norm_fit, nSteps )){
			// successfull fit :
			rms = sigma_fit;
			mean = mean_fit;			
		}
	}

	return TRUE;		
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::GetMeanAndRMS( double& mean, double& rms )
{
	register int size=m_SizeX*m_SizeY;
	double sum=0;
	double sum2=0;
	for( register int i=0;i<size;i++){
		double newval = ((double)(m_pData[i]));
		double newval2 = newval*newval;

		sum += newval;
		sum2 += newval2;
	}
	mean = (sum/size);
	rms = sqrt( sum2/size - mean*mean );
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::GetMean( double& mean, int border /*=30*/ )
{
	register int size=m_SizeX*m_SizeY;
	double sum=0;
	double sum2=0;
	for( register int y=0;y<m_SizeY;y++){
		for( register int x=0;x<m_SizeX;x++){		
			if( x>=border && y>=border && x<=(m_SizeX-border) && y<=(m_SizeY-border) ){
				sum = sum + m_pFastData[y][x];
			}
		}
	}
	mean = (sum/size);
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::GetMean( int border /*=30*/ )
{
	GetMean( m_MeanValue , border );
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::GetVariableMeanAndSigma( eLaplaceType_T laplaceType,
                        					          double& mean, double& sigma, 
																 LONG_T start_x, LONG_T start_y,
																 LONG_T end_x, LONG_T end_y )
{
	Area2DInfo info;
	GetVariableMeanAndSigma( laplaceType, mean, sigma, info,
									 start_x, start_y, end_x, end_y );
}


template<class ARG_TYPE>
void Table2D<ARG_TYPE>::GetVariableSigma( eLaplaceType_T laplaceType,double* typicalSigmaTable,
														double mean, double& sigma, Area2DInfo& info,
														LONG_T start_x, LONG_T start_y,
                                          LONG_T end_x, LONG_T end_y )
{
	double typicalSigma = typicalSigmaTable[(long)laplaceType];
	const char* laplace_name = GetLaplaceName( laplaceType );	
	double typical3Sigma = 3*typicalSigma;

	if(end_x<0 || end_x>m_SizeX)
		end_x = m_SizeX;
   if(end_y<0 || end_y>m_SizeY)
		end_y = m_SizeY;
   if(start_x<0)
		start_x = 0;
   if(start_y<0)
		start_y = 0;

	LONG_T xSize = (end_x-start_x);
	LONG_T ySize = (end_y-start_y);
	LONG_T total_size = (xSize*ySize);
	Assert(total_size<=m_SizeX*m_SizeY,"To small buffer initialized for working table ...");

	if(laplaceType!=eRawS){
		register int low_y = 4;
		register int up_y = ySize-4;		
		register int low_x = 4;
		register int up_x = xSize-4;		
		register double count=0.00;

		sigma = 0.00;
		count = 0;
		for(register int y=low_y;y<=up_y;y++){
			for(register int x=low_x;x<=up_x;x++){
				double laplace = CalcLaplaceSum( x, y, m_SizeX, m_pFastData, laplaceType );
				if(laplace<=typical3Sigma && laplace>=-typical3Sigma){
					sigma += (laplace-mean)*(laplace-mean);
					count = count + 1;
				}
			}
		}
		sigma = sqrt(sigma/count);
		// printf("mean = %f, sigma=%f, sigma2=%f, sigma_n=%f count=%f\n",mean,sigma2,sigma,(sigma/count),count);
	}else{
		double typical3Sigma = 3*typicalSigma;
		sigma = 0.00;
		register int count = 0;
		for(register int y=0;y<ySize;y++){
			for(register int x=0;x<xSize;x++){
				if(m_pFastData[y][x]<=mean+typical3Sigma && m_pFastData[y][x]>=mean-typical3Sigma){
					sigma += (m_pFastData[y][x]-mean)*(m_pFastData[y][x]-mean);
					count++;
				}
			}
		}
		sigma = sqrt(sigma/count);
	}		

	//static int i=0;
	//i++;
	//mystring szName;
	//szName << "histo_" << i << ".txt";
	//if(gCmnTrace.GetTrLevel()>=3)
	//	pHisto->DumpToFile( szName.c_str() );

}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::FillHisto( eLaplaceType_T laplaceType, int low_x, int up_x,
                          int low_y, int up_y,
								  double& mean, double& rms, double& max_val,
								  CMyHisto* pHisto )
{
	mean = 0.00;
	register int count=0;


	if( laplaceType!=eRawS ){
		for(register int y=low_y;y<=up_y;y++){
			for(register int x=low_x;x<=up_x;x++){
				double laplace = CalcLaplaceSum( x, y, m_SizeX, m_pFastData, laplaceType  );

				pHisto->Fill(laplace);
			}
		}	
	}else{
		
		for(register int y=low_y;y<=up_y;y++){
			for(register int x=low_x;x<=up_x;x++){				
				pHisto->Fill( m_pFastData[y][x] );
//				printf("%d\n",m_pFastData[y][x] );
			}
		}	
//		exit(0);
		if( pHisto->m_AllCounts<10 ){
			// change range and re-fill :
			pHisto->Init( pHisto->m_MininumValue-500, pHisto->m_MaximumValue+500, pHisto->m_BinNo );
			for(register int y=low_y;y<=up_y;y++){
				for(register int x=low_x;x<=up_x;x++){
					pHisto->Fill( m_pFastData[y][x] );
				}
			}				
		}
	}
	pHisto->GetStatValues( mean, rms, max_val );

/*	mean = pHisto->GetMeanFromTotalSum();

	rms = 0.00;
	count = 0;
	if( laplaceType!=eRawS ){
		for(register int y=low_y;y<=up_y;y++){
			for(register int x=low_x;x<=up_x;x++){
				double laplace = CalcLaplaceSum( x, y, m_SizeX, m_pFastData, laplaceType  );
				if(pHisto->IsInHisto( laplace )){
					rms += (laplace-mean)*(laplace-mean);
					count++;
				}
			}
		}
	}else{
		for(register int y=low_y;y<=up_y;y++){
			for(register int x=low_x;x<=up_x;x++){
				double laplace = m_pFastData[y][x];		
				if(pHisto->IsInHisto( laplace )){
					rms += (laplace-mean)*(laplace-mean);
					count++;
				}
			}
		}

	}
	if(count>0){
		rms = sqrt(rms/count);
	}
	max_val = pHisto->GetMaxValue();*/
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::FillHistoFromData( eLaplaceType_T laplaceType, 
								  int low_x, int up_x,
                          int low_y, int up_y,
								  double& mean, double& rms, double& max_val,
								  CMyHisto* pHisto, BIG_ELEM_TYPE** p_data )
{
	mean = 0.00;
	register int count=0;


	if( laplaceType!=eRawS ){
		for(register int y=low_y;y<=up_y;y++){
			for(register int x=low_x;x<=up_x;x++){
				pHisto->Fill( p_data[y][x] );
			}
		}	
	}else{
		for(register int y=low_y;y<=up_y;y++){
			for(register int x=low_x;x<=up_x;x++){
				pHisto->Fill( m_pFastData[y][x] );
			}
		}
	}	
	pHisto->GetStatValues( mean, rms, max_val );


/*	mean = pHisto->GetMeanFromTotalSum();

	rms = 0.00;
	count = 0;


	if( laplaceType!=eRawS ){
		for(register int y=low_y;y<=up_y;y++){
			for(register int x=low_x;x<=up_x;x++){
				register double laplace = p_data[y][x];
				if(pHisto->IsInHisto( laplace )){
					rms += (laplace-mean)*(laplace-mean);
					count++;
				}
			}
		}
	}else{
		for(register int y=low_y;y<=up_y;y++){
			for(register int x=low_x;x<=up_x;x++){
				register double laplace = p_data[y][x];
				if(pHisto->IsInHisto( laplace )){
					rms += (laplace-mean)*(laplace-mean);
					count++;
				}
			}
		}

	}
	rms = sqrt(rms/count);
	max_val = pHisto->GetMaxValue();*/
}


template<class ARG_TYPE>
BOOL_T Table2D<ARG_TYPE>::GetVariableMeanAndSigma( eLaplaceType_T laplaceType,
																 double& mean, double& sigma, 
																 Area2DInfo& info,
                        					          LONG_T start_x, LONG_T start_y,
					                                  LONG_T end_x, LONG_T end_y, 
																 int frame_index,
																 BIG_ELEM_TYPE** p_lap_data,
																 BOOL_T bUseMeanOnFitFailed /*=FALSE*/ )
{
	
	const char* laplace_name = GetLaplaceName( laplaceType );
	CMyHisto* pHisto = info.m_HistoList.Find( laplace_name );

	double typicalSigma = gSigmaBacgroundG[(long)laplaceType];
	if(pHisto){
		if(pHisto->WasPrevFitOK()){
			typicalSigma = pHisto->GetSigmaFit();
		}else{
			typicalSigma = pHisto->m_RMSValue;
		}
	}
	
	int binNo=60;
	double minVal = -12*typicalSigma;
   double maxVal = +12*typicalSigma;

	// raw value caulculation :
	double mean_tmp=0;
	double sum2_tmp=0;
	double sigma_tmp=0;

	// laplace :	
	double mean_lap=0;
	double sum2_lap=0;
	double sigma_lap=0;

	if( laplaceType == eRawS ){
		// calculating mean/sigma to declare histogram :
		minVal=100000;
		maxVal=-100000;
		int size=m_SizeX*m_SizeY;
		int count_tmp=0;

		for( int y=start_y;y<end_y;y++ ){
         for( int x=start_x;x<end_x;x++ ){
				double newval = ((double)(m_pFastData[y][x]));

				mean_tmp += newval;
				sum2_tmp += newval*newval;
				count_tmp++;
		
				if( m_pFastData[y][x] < minVal ){
					minVal = m_pFastData[y][x];
				}
				if( m_pFastData[y][x] > maxVal ){
					maxVal = m_pFastData[y][x];
				}
				// printf("%d\n",m_pData[i]);


//				if( frame_index==0 ){
//					printf("TEMPORARY_DUMP = %d %d %d\n",x,y,m_pFastData[y][x]);
//				}
			}
		}

		mean_tmp = (mean_tmp/count_tmp);
		double mean2_tmp = (sum2_tmp/count_tmp);
		sigma_tmp = sqrt( mean2_tmp - mean_tmp*mean_tmp );
		// exit(0);


		if( mean_tmp <= 0 ){
			printf("ERROR :Empty frame ?\n");			
			mean = mean_tmp;
			sigma = 0;
			return FALSE;
		}
	}else{
		if( !p_lap_data ){
			printf("ERROR IN CONFIG : TURNED OFF CALCULATION OF BACKGROUND FOR LAPLACE %d\n",laplaceType);
			return FALSE;
		}

		int count_lap=0;
		
		for( int y=start_y;y<end_y;y++ ){
			for( int x=start_x;x<end_x;x++ ){
				if( p_lap_data[y][x]>=minVal && p_lap_data[y][x]<=maxVal ){
					double newval=p_lap_data[y][x];

					mean_lap += newval;
					sum2_lap += newval*newval;
		
					count_lap++;

//					if( frame_index==0 ){
//						printf("TEMPORARY_DUMP = %d %d %d\n",x,y,p_lap_data[y][x]);
//					}
				}
			}
		}
		mean_lap = (mean_lap/count_lap);
	   double mean2_lap = (sum2_lap/count_lap);
   	sigma_lap = sqrt( mean2_lap - mean_lap*mean_lap );

		minVal = mean_lap - 5*sigma_lap;
		maxVal = mean_lap + 5*sigma_lap;
	}

	BOOL_T bDoAdd=TRUE;
	if( pHisto ){
		bDoAdd = FALSE;
	}

	if( laplaceType == eRawS ){	
		if( !pHisto ){
			double max_histo_val = MAX((2*mean_tmp),(mean_tmp+5*sigma_tmp));
			pHisto = new CMyHisto( laplace_name, 0, max_histo_val, binNo );
		}else{
			pHisto->Init( 0, 2*mean_tmp, binNo );
		}
	}else{
		// temporary - NEW change 
		if( minVal >= maxVal ){
			printf("Bad ranges ! %.2f , %.2f\n",minVal, maxVal);

			if( fabs(mean_lap)<0.1 && fabs(sigma_lap)<0.1 ){
				printf("WARNING : all values in window (%d,%d)-(%d,%d) are 0 ???!!!\n",start_x,start_y,end_x,end_y);

				if( bUseMeanOnFitFailed ){
					printf("WARNING : setting large value sigma=1000.00 ADU\n");

					mean = 0.00;
					sigma = 1000.00;
					
					info.m_DataInfo[ laplaceType ].m_Average = mean;
			  	   info.m_DataInfo[ laplaceType ].m_Sigma = sigma;
			  	   info.m_DataInfo[ laplaceType ].m_bFitOK = FALSE;

					return TRUE;
				}
			}

			return FALSE;
		}
		
		if( !pHisto ){
		   pHisto = new CMyHisto( laplace_name, minVal, maxVal,  binNo );	
		}else{
			pHisto->Init( minVal, maxVal,  binNo );
		}
	}
	if( bDoAdd ){
	   info.m_HistoList.Add( pHisto );
	}

	pHisto->Clear();



	if(end_x<0 || end_x>m_SizeX)
		end_x = m_SizeX;
   if(end_y<0 || end_y>m_SizeY)
		end_y = m_SizeY;
   if(start_x<0)
		start_x = 0;
   if(start_y<0)
		start_y = 0;

	LONG_T xSize = (end_x-start_x);
	LONG_T ySize = (end_y-start_y);
	LONG_T total_size = (xSize*ySize);
	Assert(total_size<=m_SizeX*m_SizeY,"To small buffer initialized for working table ...");


	mean = 0.00;
	register int low_y = start_y+4;
	register int up_y = end_y-4;
	register int low_x = start_x + 4;
	register int up_x = end_x - 4;
	register double count=0.00;

	double norm;
	if( p_lap_data ){
		FillHistoFromData( laplaceType, low_x, up_x, low_y, up_y, mean, sigma,
								 norm, pHisto, p_lap_data );
	}else{
		// skiping saturated pixels :
		pHisto->m_SkipBigger = MAX_VALUE_ALLOWED;

		FillHisto( laplaceType, low_x, up_x, low_y, up_y, mean, 
					  sigma, norm, pHisto );
	}
		
	if(!pHisto->IsCharCalculated()){
		// not first time :
		pHisto->SetCharValues( mean, sigma, norm );
	}
	
	// now fit gauss using mean,rms,max :
	double mean_fit = mean;
	double sigma_fit = sigma;
	double norm_fit = norm;

	static int i=0;
	i++;

	PROFILER_START
	int nSteps=500;
	BOOL_T bFitGaussResult = pHisto->FitGauss(  mean_fit, sigma_fit, norm_fit, 
															  nSteps, TRUE, bUseMeanOnFitFailed, 
															  FALSE );

	// INFO : newer root fits much better, maybe it is not needed at all :
	// try to re-fit :
	if( !bFitGaussResult ){
		if( pHisto ){
			printf("WARNING : gauss fit failed in window (%d,%d)-(%d,%d), retrying in different binning\n",low_x, up_x, low_y, up_y);
			printf("WARNING : retrying in different binning ...\n");fflush(stdout);
			if( laplaceType == eRawS ){	
				// typicaly fit failes, because there is only 1-3 points > 0 and very short spike, this 
				// means that histogram range (initial RMS/sigma) is too high, so i try to make 
				// histogram range smaller and re-try fit, that is why there is n_sigma=1 here :
				double n_sigma=1.00; 

				double min_histo_val = mean_tmp-n_sigma*sigma_tmp;
				double max_histo_val = mean_tmp+n_sigma*sigma_tmp;
				pHisto->Init( min_histo_val, max_histo_val, binNo );
				printf("Histogram re-initialized with range %.2f - %.2f ( from %.2f +/- %.2f * %.2f ) , binning = %d\n",
							min_histo_val,max_histo_val,mean_tmp,n_sigma,sigma_tmp,binNo);
			}else{
				double n_sigma=3.00;
				double min_histo_val = mean_lap-n_sigma*sigma_lap;
         	double max_histo_val = mean_lap+n_sigma*sigma_lap;
				pHisto->Init( min_histo_val, max_histo_val, binNo );
				printf("Histogram re-initialized with range %.2f - %.2f ( from %.2f +/- %.2f * %.2f ) , binning = %d\n",
							min_histo_val,max_histo_val,mean_tmp,n_sigma,sigma_tmp,binNo);
			}

			double norm;
			if( p_lap_data ){
				FillHistoFromData( laplaceType, low_x, up_x, low_y, up_y, mean, sigma,
									 norm, pHisto, p_lap_data );
			}else{
				// skiping saturated pixels :
				pHisto->m_SkipBigger = MAX_VALUE_ALLOWED;

				FillHisto( laplaceType, low_x, up_x, low_y, up_y, mean, 
						  sigma, norm, pHisto );
			}	
		
			if(!pHisto->IsCharCalculated()){
				// not first time :
				pHisto->SetCharValues( mean, sigma, norm );
			}
	
			int nSteps=500;
			bFitGaussResult = pHisto->FitGauss(  mean_fit, sigma_fit, norm_fit, nSteps, TRUE, bUseMeanOnFitFailed );
			if ( bFitGaussResult ){
				printf("SUCCESS : re-fit of gauss ok !\n");

				if( TRUE ){
					char szPrefix[128];
			      sprintf(szPrefix,"refitted_%.5d_%.5d",start_x,start_y);   
					mystring szName;
			      szName = pHisto->DumpToFile( frame_index, szPrefix );
					printf("DEBUG : refitted histogram dumped to file %s\n",szName.c_str());
				}
			}else{
				printf("WARNING : re-fit of gauss failed !\n");
			}
		}else{
			printf("ERROR : cannot re-fit pHisto=NULL error in code ???\n");
		}
	} // end of re-fit 


	mystring szName;
	if((!bFitGaussResult && gDoDumpBadFit) || gDoDumpAllHisto){
		char szPrefix[128];
      sprintf(szPrefix,"err_%.5d_%.5d",start_x,start_y);
      szName = pHisto->DumpToFile( frame_index, szPrefix );
	}

	if( !bFitGaussResult ){
		printf("DEBUG : gauss fit result = %d ( printf_level = %d )\n",bFitGaussResult,gPrintfLevel);fflush(stdout);
	
		printf("---------------------------------\n");
		printf("ERROR_GAUSS : could not fit GAUSS !!! mean=%f, sigma=%f, norm=%f\n",mean_fit, sigma_fit, norm_fit);
		printf("ERROR_GAUSS_CODE : %d\n",pHisto->m_ErrorCode);
		printf("USING RMS/MEAN : %f / %f\n",sigma,mean);
		if(strlen(szName.c_str())){
			printf("File name : %s\n",szName.c_str() );
		}
		printf("---------------------------------\n");fflush(stdout);

		if( bUseMeanOnFitFailed ){
			if( sigma_fit < 50 ){
				printf("WARNING !!! : RMS < 50 changed to 500 !!!\n");fflush(stdout);
				sigma_fit = 500.00;
			}
			printf("GAUSS FIT FAILED - ACCEPTING MEAN=%.2f and RMS=%.2f\n",mean_fit, sigma_fit);
			bFitGaussResult = TRUE;
		}
		info.m_DataInfo[ laplaceType ].m_eFitResType = eFitResAVG;
	}else{
		info.m_DataInfo[ laplaceType ].m_eFitResType = eFitResGauss;
	}


	if( bFitGaussResult ){
		mean = mean_fit;
		sigma = sigma_fit;
		norm = norm_fit;
		_TRACE_PRINTF_3("%s FITTED GAUSS TO BACKGROUND IN WINDOW (%d,%d)-(%d,%d): sigma=%.5f, mean=%.5f, norm=%.5f\n",
					laplace_name,low_x,low_y,up_x,up_y,sigma,mean,norm);		
		if(gDoDumpAllHisto){
			pHisto->DumpToFile( frame_index );
		}
	}
	if( gPrintfLevel>=3 ){
		PROFILER_END_CMN("Fit of gauss took :");
	}
	info.m_DataInfo[ laplaceType ].m_Average = mean;
	info.m_DataInfo[ laplaceType ].m_Sigma = sigma;
	info.m_DataInfo[ laplaceType ].m_bFitOK = bFitGaussResult;

	return bFitGaussResult;
}

template<class ARG_TYPE>
LONG_T Table2D<ARG_TYPE>::GetMedian( LONG_T MaxValue,
					                      LONG_T start_x, LONG_T start_y,
                					       LONG_T end_x, LONG_T end_y )
{
	if(end_x<0 || end_x>m_SizeX)
		end_x = m_SizeX;
   if(end_y<0 || end_y>m_SizeY)
		end_y = m_SizeY;
   if(start_x<0)
		start_x = 0;
   if(start_y<0)
		start_y = 0;

	LONG_T xSize = (end_x-start_x);
	LONG_T ySize = (end_y-start_y);
	LONG_T total_size = (xSize*ySize);
	LONG_T median_pos=(total_size/2);

	if(!m_pValuesTable)
		m_pValuesTable = new LONG_T[total_size];
	register int i=0;
	for(register int y=start_y;y<end_y;y++){
      for(register int x=start_x;x<end_x;x++){
			m_pValuesTable[i] = (LONG_T)m_pFastData[y][x];
			i++;			
		}
	}			
	my_qsort(m_pValuesTable,total_size);
	LONG_T median = m_pValuesTable[median_pos];
	m_MedianValue = median;
	return m_MedianValue;
}


template<class ARG_TYPE>
LONG_T Table2D<ARG_TYPE>::GetFastMedian( LONG_T MaxValue,
					                  	     LONG_T start_x, LONG_T start_y,
                					      	  LONG_T end_x, LONG_T end_y )
{
	if(end_x<0 || end_x>m_SizeX)
		end_x = m_SizeX;
   if(end_y<0 || end_y>m_SizeY)
		end_y = m_SizeY;
   if(start_x<0)
		start_x = 0;
   if(start_y<0)
		start_y = 0;

	LONG_T xSize = (end_x-start_x);
	LONG_T ySize = (end_y-start_y);
	LONG_T total_size = (xSize*ySize);
	LONG_T median_pos=(total_size/2);

	CMyHisto histo( "median", 0, 65555, 65556 );

	register int i=0;
	for(register int y=start_y;y<end_y;y++){
      for(register int x=start_x;x<end_x;x++){
			histo.Fill( m_pFastData[y][x] );
		}
	}	
	
	int median=0;
	int sum=0;
	int all_counts=histo.m_AllCounts;		
	int all_counts_div2=(all_counts/2);
	for(register int j=0;j<histo.m_BinNo;j++){
		sum += histo.m_pCountTab[j];
		// printf("j=%d, sum =%d, bin=%d\n",j,sum,histo.m_pCountTab[j]);
		if(sum>all_counts_div2){
			double val1 = histo.GetBinCenter( j-1 );			
			double val2 = histo.GetBinCenter( j );
		 	m_MedianValue = (LONG_T)((val1 + val2)/2);
			break;
		}
	}

	return m_MedianValue;
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::GetDistrib( LONG_T* pDistribTab, LONG_T maxValue )
{
	InitDistribTab( maxValue );
	for(register int i=0;i<m_Size;i++){
		if(m_pData[i]<=maxValue && m_pData[i]>=0){
			pDistribTab[ (long) m_pData[i] ]++;
		}
	}	
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::GetDistrib( ImageStat& statInfo, LONG_T maxValue )
{
	statInfo.AllocDistribTab( maxValue );
	GetDistrib( statInfo.m_pValDistrib, maxValue );
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::GetBackgrDescStr( InfoTable2D& info, mystring& szDesc )
{
	//szDesc = "G54 : ";
	//for(register int x=0;x<info.m_X_count;x++){
   //  for(register int y=0;y<info.m_Y_count;y++){
	//		szDesc << " mean= " << info.GetElem(x,y).m_AverageG54 << " RMS=" << info.GetElem(x,y).m_SigmaG54 << "\n";
	//	}
	//}
	szDesc = "";
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::CalcTableMapInfo( InfoTable2D& info, 
								  LONG_T maxValue, 
								  BOOL_T* pBackgrFlagTable,
								  eLaplaceType_T currLaplace,
								  BIG_ELEM_TYPE** p_lap_data,
								  BOOL_T bCalcMostPopular /*=FALSE*/, 
	                       BOOL_T bCalcMedian /*=FALSE*/,
	                       int frame_index /*=0*/,
								  int border/*=0*/ )
{
	if(m_SizeX!=info.m_SizeX || m_SizeY!=info.m_SizeY){
		// NEW change due to SLT - 20041022 
		// info.Init( m_SizeX, m_SizeY, m_SizeX/10, m_SizeY/10 );
		info.Init( m_SizeX, m_SizeY, m_SizeX, m_SizeY );
	}
	for(register int x=0;x<info.m_X_count;x++){
		for(register int y=0;y<info.m_Y_count;y++){
			Area2DInfo& elem = info.GetElem(x,y);

			/*if(bCalcMostPopular)
				info.GetElem(x,y).m_MostPopularValue = GetMostPopular( maxValue,
																x*info.m_dX, y*info.m_dY,
																(x+1)*info.m_dX, (y+1)*info.m_dY );*/
			if(bCalcMedian)
				info.GetElem(x,y).m_MedianValue = GetMedian( maxValue, x*info.m_dX, y*info.m_dY,
																			(x+1)*info.m_dX, (y+1)*info.m_dY );
			for(int i=0;i<MAX_LAPLACE_DEFINED;i++){
				eLaplaceType_T lap_type = (eLaplaceType_T)i;
				int laptype_int = (int)lap_type;
				if( pBackgrFlagTable[i] ){
					BIG_ELEM_TYPE** p_laplace_data = NULL;
					if( lap_type==currLaplace && p_lap_data ){
						p_laplace_data =  p_lap_data;
					}
					GetVariableMeanAndSigma( lap_type, 
													 elem.m_DataInfo[laptype_int].m_Average,
													 elem.m_DataInfo[laptype_int].m_Sigma,
													 elem,
													 (int)elem.m_LowLeft.x, (int)elem.m_LowLeft.y,
													 (int)elem.m_UpRight.x, (int)elem.m_UpRight.y, 
													 frame_index, p_laplace_data );
					_TRACE_PRINTF_2("Background of laplace : %s, mean=%.2f, sigma=%.2f\n",
							  GetLaplaceName( lap_type ), 
							  elem.m_DataInfo[laptype_int].m_Average,
							  elem.m_DataInfo[laptype_int].m_Sigma );
				}
			}
		}
	}	
}

template<class ARG_TYPE>
LONGLONG_T Table2D<ARG_TYPE>::CalcSum( LONG_T* pixels, LONG_T cnt )
{
   LONGLONG_T sum=0;
   for(register int i=0;i<cnt && pixels[i]!=NEIGHB_LIST_END;i++){
		sum += (LONGLONG_T)(m_pData[ pixels[i] ]);
   }
   return sum;
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::GetXYRotCorrected( LONG_T x, LONG_T y, LONG_T dx, LONG_T dy,
                        					    LONG_T& x_corr, LONG_T& y_corr )
{
	x_corr = x + dx;
	y_corr = y + dy;	
}

template<class ARG_TYPE>
LONGLONG_T Table2D<ARG_TYPE>::CalcSumRotCorrected( LONG_T* pixels, LONG_T cnt,
                      		         	            double dx, double dy )
{
	LONGLONG_T sum=0;
   for(register int i=0;i<cnt && pixels[i]!=NEIGHB_LIST_END;i++){
		LONG_T pos = pixels[i];
		LONG_T x = (pos % m_SizeX);
		LONG_T y = (pos / m_SizeX);
		LONG_T real_x = my_round(x+dx);
		LONG_T real_y = my_round(y+dy);

		if(real_x>=0 && real_x<m_SizeX && real_y>=0 || real_y<m_SizeY)
	      sum += (LONGLONG_T)(m_pFastData[real_y][real_x]);
   }
   return sum;
}


template<class ARG_TYPE>
void Table2D<ARG_TYPE>::Dump( long _x0, long _y0, long _x1, long _y1 )
{
	long x0 = MAX(0,_x0);
	long y0 = MAX(0,_y0);
	long x1 = MIN((m_SizeX-1),_x1);
	long y1 = MIN((m_SizeY-1),_y1);

	for(register int y=y1;y>=y0;y--){
		for(register int x=x0;x<=x1;x++){
			_TRACE_PRINTF_6("(%d,%d)=%d ",y,x,m_pFastData[y][x]);
		}
		_TRACE_PRINTF_6("\n");
	}
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::Dump()
{
	for(register int y=0;y<m_SizeY;y++){
		for(register int x=0;x<m_SizeX;x++){
			printf("%d ",getval(x,y));
		}
		printf("\n");
	}
}


template<class ARG_TYPE>
const char* Table2D<ARG_TYPE>::getKey(const char* keyname)
{
	CEnvVar* pKey = m_KeyTab.Find( keyname );
	if(pKey)
		return (pKey->szValue).c_str();
	return NULL;
}

template<class ARG_TYPE>
BOOL_T Table2D<ARG_TYPE>::Compare( const Table2D<ARG_TYPE>& right )
{
	Assert(m_SizeX==right.m_SizeX,"X - sizes must be equal");
	Assert(m_SizeY==right.m_SizeY,"Y - sizes must be equal");

	BOOL_T bRes = TRUE;
	long size = (m_SizeX*m_SizeY);
	for(int i=0;i<size;i++){
		if(m_pData[i] != right.m_pData[i]){
			printf("Difference at (%d,%d) %d!=%d\n",(i % m_SizeX),(i / m_SizeX),m_pData[i],right.m_pData[i]);fflush(0);
			bRes=FALSE;
		}
	}
	return bRes;
}

template<class ARG_TYPE>
Table2D<ARG_TYPE>& Table2D<ARG_TYPE>::operator-=(const Table2D<ARG_TYPE>& right)
{
	for(register long pos=0;pos<m_Size;pos++){
		m_pData[pos] -= right.m_pData[pos];
	}	
	return (*this);
}

template<class ARG_TYPE>
void  Table2D<ARG_TYPE>::SubtractBacgroundInFrame( long nFrameSize ) 
{
	long border_sum=0;
	long count=0;
	for(register long pos=0;pos<m_Size;pos++){
		long x,y;
		GetXY_FromPos( pos , x , y );
		if( (x<nFrameSize) || (y<nFrameSize) || (x>(m_SizeX-1-nFrameSize)) || (y>(m_SizeY-1-nFrameSize)) ){
            border_sum += (long)m_pData[pos];
            count++;
		}
	}
	long backgr_level = (border_sum/count);
	// printf("Calculated border sum=%d, background_level=%d\n",border_sum,backgr_level);
   SubtractConst( backgr_level );
}

template<class ARG_TYPE>
void  Table2D<ARG_TYPE>::SubtractConst( long value, BOOL_T bZero )
{
	long size = (m_SizeX*m_SizeY);
   for(register int i=0;i<size;i++){
		int val = ( m_pData[i] - value );
		if( val<0 && bZero ){
			val = 0;
		}
		m_pData[i] = val;
   }
	CalcMaxAndPos();
}

template<class ARG_TYPE>
void  Table2D<ARG_TYPE>::MultiplyByConst( double mult )
{
	long size = (m_SizeX*m_SizeY);
   for(register int i=0;i<size;i++){
		m_pData[i] = (ARG_TYPE)(m_pData[i]*mult); 
   }
	CalcMaxAndPos();
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::Subtract( Table2D<ARG_TYPE>& right,Table2D<ARG_TYPE>& result, BOOL_T bZero)
{
	Assert(m_SizeX==right.m_SizeX,"X - sizes must be equal");
	Assert(m_SizeY==right.m_SizeY,"Y - sizes must be equal");
	result.Alloc(m_SizeX,m_SizeY);

	register int size = (m_SizeX*m_SizeY);
	for(register int i=0;i<size;i++){
		register ARG_TYPE original = m_pData[i];
		result.m_pData[i] = original - right.m_pData[i];
		if(bZero){
			if(original < right.m_pData[i])
				result.m_pData[i] = 0;
				// Assert(FALSE,"Error in analysis value in pixel <0");
		}
	}
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::Multiply(Table2D<ARG_TYPE>& right,Table2D<ARG_TYPE>& result,
											int nNorm, int max_val)
{
	Assert(m_SizeX==right.m_SizeX,"X - sizes must be equal");
	Assert(m_SizeY==right.m_SizeY,"Y - sizes must be equal");
	result.Alloc(m_SizeX,m_SizeY);

	register int size = (m_SizeX*m_SizeY);
	for(register int i=0;i<size;i++){
		int val=0;
	
		double ratio=(double)right.m_pData[i];
		if( nNorm>0 )
			ratio = ratio / ( (double)nNorm );
		if(right.m_pData[i]!=0){
			val = (ARG_TYPE)(m_pData[i]*ratio);
		}

		// in case limit of data type exceeded, setmax value :
		if(val>max_val)
			val = max_val;
		result.m_pData[i] = val;
	}
}


template<class ARG_TYPE>
BOOL_T Table2D<ARG_TYPE>::GetImageWithZeros( long x0, long y0,long len_x, long len_y,
                    Table2D<ARG_TYPE>& Image )
{
	long x1 = x0 + len_x - 1;
	long y1 = y0 + len_y - 1;
	BOOL_T bRet = TRUE;

	Image.Alloc( len_x, len_y );

	ARG_TYPE** p_out = Image.get_data_buffer_fast();
	for(register long y=y0;y<=y1;y++){
		for(register long x=x0;x<=x1;x++){
			if(y>=0 && y<m_SizeY && x>=0 && x<m_SizeX)
				p_out[y-y0][x-x0] = m_pFastData[y][x];
			else{
				p_out[y-y0][x-x0] = 0;
				bRet = FALSE;
			}
		}
	}
	return bRet;		
	
}

template<class ARG_TYPE>
BOOL_T Table2D<ARG_TYPE>::MarkStarWithCircle( int x0, int y0, int radius, 
															 int width /* =2 */,
															 int value /* =50000 */ )
{
	int start_x = MAX((x0-radius),0);
	int start_y = MAX((y0-radius),0);	
	int end_x = MIN((x0+radius),(m_SizeX-1));
	int end_y = MIN((y0+radius),(m_SizeY-1));

	BOOL_T bRet=FALSE;
	for(int x=start_x;x<=end_x;x++){
		for(int y=start_y;y<=end_y;y++){
			double dist = sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0) );
			if( ((dist - radius))<=width && dist>radius ){
				m_pFastData[y][x] = value;
				bRet = TRUE;
			}
		}
	}

	return bRet;
}

template<class ARG_TYPE>
BOOL_T Table2D<ARG_TYPE>::GetImage( long x0, long y0,long len_x, long len_y,
                            Table2D<ARG_TYPE>& Image )
{
	long inital_len_x = len_x;
	long inital_len_y = len_y;
	long x0_in = x0;
	long y0_in = y0; 


	if(x0<0)
		x0_in = 0;
	if(y0<0)
		y0_in = 0; 

	Image.m_X_On_Big = x0_in;
	Image.m_Y_On_Big = y0_in;	
	Image.m_StartOnOrigX = x0_in;
	Image.m_StartOnOrigY = y0_in;
	long x1 = x0_in + len_x - 1;
	long y1 = y0_in + len_y - 1;

	if(x1<0 || y1<0)
		return FALSE;
	if(x1>=m_SizeX)
		x1 = m_SizeX-1;
	if(y1>=m_SizeY)
		y1 = m_SizeY-1;
	len_x = (x1-x0_in+1);
	len_y = (y1-y0_in+1);
	
	Image.Alloc( len_x, len_y );

	ARG_TYPE* left = GetPos(x0_in,y0_in);
	ARG_TYPE* right = Image.GetPos(0,0);

	for(int i=0;i<len_y;i++){
		memcpy(right,left,sizeof(ARG_TYPE)*len_x);
		left = left + m_SizeX;
		right = right + len_x;
	}		
	return (inital_len_x==len_x && inital_len_y==len_y);		
}

template<class ARG_TYPE>
BOOL_T Table2D<ARG_TYPE>::GetImageELEMTYPE( long x0, long y0,long len_x, long len_y,
                            Table2D<ELEM_TYPE>& Image )
{
	long inital_len_x = len_x;
	long inital_len_y = len_y;
	long x0_in = x0;
	long y0_in = y0; 


	if(x0<0)
		x0_in = 0;
	if(y0<0)
		y0_in = 0; 

	Image.m_X_On_Big = x0_in;
	Image.m_Y_On_Big = y0_in;	
	Image.m_StartOnOrigX = x0_in;
   Image.m_StartOnOrigY = y0_in;
	long x1 = x0_in + len_x - 1;
	long y1 = y0_in + len_y - 1;

	if(x1<0 || y1<0)
		return FALSE;
	if(x1>=m_SizeX)
		x1 = m_SizeX-1;
	if(y1>=m_SizeY)
		y1 = m_SizeY-1;
	len_x = (x1-x0_in+1);
	len_y = (y1-y0_in+1);
	
	Image.Alloc( len_x, len_y );

	ARG_TYPE* left = GetPos(x0_in,y0_in);
	ELEM_TYPE* right = Image.GetPos(0,0);

	for(register int i=0;i<len_y;i++){
		for(register int j=0;j<len_x;j++){
			right[j] = (ELEM_TYPE)left[j];
		}
		left = left + m_SizeX;
		right = right + len_x;
	}		
	return (inital_len_x==len_x && inital_len_y==len_y);		
}

template<class ARG_TYPE>
BOOL_T Table2D<ARG_TYPE>::GetImageBIGELEMTYPE( long x0, long y0,long len_x, long len_y,
                            Table2D<BIG_ELEM_TYPE>& Image )
{
	long inital_len_x = len_x;
	long inital_len_y = len_y;
	long x0_in = x0;
	long y0_in = y0; 


	if(x0<0)
		x0_in = 0;
	if(y0<0)
		y0_in = 0; 

	Image.m_X_On_Big = x0_in;
	Image.m_Y_On_Big = y0_in;	
	Image.m_StartOnOrigX = x0_in;
   Image.m_StartOnOrigY = y0_in;
	long x1 = x0_in + len_x - 1;
	long y1 = y0_in + len_y - 1;

	if(x1<0 || y1<0)
		return FALSE;
	if(x1>=m_SizeX)
		x1 = m_SizeX-1;
	if(y1>=m_SizeY)
		y1 = m_SizeY-1;
	len_x = (x1-x0_in+1);
	len_y = (y1-y0_in+1);
	
	Image.Alloc( len_x, len_y );

	ARG_TYPE* left = GetPos(x0_in,y0_in);
	BIG_ELEM_TYPE* right = Image.GetPos(0,0);

	for(register int i=0;i<len_y;i++){
		for(register int j=0;j<len_x;j++){
			right[j] = (BIG_ELEM_TYPE)left[j];
		}
		left = left + m_SizeX;
		right = right + len_x;
	}		
	return (inital_len_x==len_x && inital_len_y==len_y);		
}


template<class ARG_TYPE>
void Table2D<ARG_TYPE>::ShiftFrameWeighted( LONG_T FrameCounter, double frame_dx, double frame_dy )
{
	LONG_T nSteps = ( FrameCounter - 1 );
	if(nSteps<0)
      nSteps = 0;
	m_totalShiftX += frame_dx;
	m_totalShiftY += frame_dy;
	
	double realStepX = m_totalShiftX-(double)m_PrevShiftX;
	double realStepY = m_totalShiftY-(double)m_PrevShiftY;
	LONG_T newStepX = my_round( realStepX );
	LONG_T newStepY = my_round( realStepY );
	Table2D<ARG_TYPE> WrkTable( m_SizeX,m_SizeY );

	_TRACE_PRINTF_5("Shift values : steps=%d,totalX=%f,totalY=%f,prevX=%d,prevY=%d,StepX=%d,StepY=%d\n",
			nSteps,m_totalShiftX,m_totalShiftY,m_PrevShiftX,m_PrevShiftY,newStepX,newStepY);	
	
	register long new_x=0;
	register long new_y=0;
	double real_new_x = 0; 
	double real_new_y = 0; 

		
	if(newStepX<=0 && newStepY<=0){
		// only for moving left/down !

		for( register int y=0;y<m_SizeY;y++){
			new_y = y + newStepY;
			for( register int x=0;x<m_SizeX;x++){
				new_x = x + newStepX;

										
				if(new_x>=0 && new_y>=0)	
					m_pFastData[new_y][new_x] = m_pFastData[y][x];
			}
		}
	}else{
		Assert(FALSE,"Translation not implemented for this values of m_FrameDX=%f and m_FrameDY=%f, newStepX=%d, newStepY=%d",frame_dx,frame_dy,newStepX,newStepY);
	}
	
	m_PrevShiftX += newStepX;
	m_PrevShiftY += newStepY;


}


template<class ARG_TYPE>
void Table2D<ARG_TYPE>::RotateFrame( LONG_T FrameCounter, double dAlpha, double rotCenterX, double rotCenterY ) 
{

	/*for(register long y0=0;y0<m_SizeY;y0++){
		double y0_prim = y0-rotCenterY;
		double y0_prim_sq = (y0_prim*y0_prim);
		for(register long x0=0;x0<m_SizeX;x0++){
			double x0_prim = x0-rotCenterX;
		   double r = sqrt(x0_prim*x0_prim + y0_prim_sq);

			double x0_rot = rotCenterX + r*cos(nSteps*dAlfa+alfa_0);
		   double y0_rot = rotCenterY + r*sin(nSteps*dAlfa+alfa_0);

	}*/

	
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::GetTotalShift( LONG_T FrameCounter, double frame_dx,
                                       double frame_dy, 
													double& totalShiftX, double& totalShiftY, LONG_T& nSteps )
{
	nSteps = (FrameCounter - 1);
   if(nSteps<0)
   	nSteps = 0;		
	totalShiftX = nSteps*frame_dx;
	totalShiftY = nSteps*frame_dy;
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::CalcTotalShift( LONG_T FrameCounter, double frame_dx, double frame_dy )
{
	LONG_T nSteps;
	GetTotalShift( FrameCounter, frame_dx, frame_dy, m_totalShiftX, m_totalShiftY, nSteps );
}


template<class ARG_TYPE>
void Table2D<ARG_TYPE>::ShiftFrame( LONG_T FrameCounter, double frame_dx, 
												double frame_dy, BOOL_T bUseTotal/*=TRUE*/ )
{
	register LONG_T newStepX = 0;
	register LONG_T newStepY = 0;	

	if(bUseTotal){
		LONG_T nSteps = (FrameCounter - 1);
		GetTotalShift( FrameCounter, frame_dx, frame_dy, m_totalShiftX, m_totalShiftY, nSteps );
	
		double realStepX = m_totalShiftX-(double)m_PrevShiftX;
		double realStepY = m_totalShiftY-(double)m_PrevShiftY;
		newStepX = my_round( realStepX );
		newStepY = my_round( realStepY );
		_TRACE_PRINTF_5("Shift values : steps=%d,totalX=%f,totalY=%f,prevX=%d,prevY=%d,StepX=%d,StepY=%d\n",
				nSteps,m_totalShiftX,m_totalShiftY,m_PrevShiftX,m_PrevShiftY,newStepX,newStepY);	

	}else{
		newStepX = my_round( frame_dx );
		newStepY = my_round( frame_dy );
		_TRACE_PRINTF_3("Shift values : StepX=%d,StepY=%d\n",newStepX,newStepY);
	}

	
	register long new_x=0;
	register long new_y=0;
	double real_new_x = 0; 
	double real_new_y = 0; 

		
	if(newStepX<=0 && newStepY<=0){
		// only for moving left/down !

		for( register int y=0;y<m_SizeY;y++){
			new_y = y + newStepY;
			if(new_y>=0){
				for( register int x=0;x<m_SizeX;x++){
					new_x = x + newStepX;

					/*real_new_x = x + realStepX;
					real_new_y = y + realStepY;*/
										
					if(new_x>=0)	
						m_pFastData[new_y][new_x] = m_pFastData[y][x];
				}
			}
		}
	}else{
		if(newStepX<=0 && newStepY>0){
			for( register int y=m_SizeY-1;y>=0;y--){
				new_y = y + newStepY;
				if(new_y<m_SizeY){
					for( register int x=0;x<m_SizeX;x++){
						new_x = x + newStepX;

						if(new_x>=0)
							m_pFastData[new_y][new_x] = m_pFastData[y][x];
					}
				}
			}
		}else{
			Assert(FALSE,"Translation not implemented for this values of m_FrameDX=%f and m_FrameDY=%f, newStepX=%d, newStepY=%d",frame_dx,frame_dy,newStepX,newStepY);
		}
	}
	
	m_PrevShiftX += newStepX;
	m_PrevShiftY += newStepY;
}



template<class ARG_TYPE>
void Table2D<ARG_TYPE>::CalcOverlapParts2( double dx, double dy, ShiftInfo& pShiftInfo )
{
	CalcOverlapParts( dx, dy, pShiftInfo.m_LocalS0, pShiftInfo.m_LocalS1,
			  pShiftInfo.m_LocalS2, pShiftInfo.m_LocalS3,
			  pShiftInfo.overlapledPoints );

}


template<class ARG_TYPE>
void Table2D<ARG_TYPE>::CalcOverlapParts( double dx, double dy, double& s0,
                                          double& s1, double& s2, double& s3,
														CLongPoint* overlapedPixels )
{
	overlapedPixels[0].x = (long)floor(dx);
	overlapedPixels[0].y = (long)floor(dy);

	overlapedPixels[1].x = (long)floor(dx);
	overlapedPixels[1].y = (long)floor(dy+1);

	overlapedPixels[2].x = (long)floor(dx+1);
	overlapedPixels[2].y = (long)floor(dy+1);

	overlapedPixels[3].x = (long)floor(dx+1);
	overlapedPixels[3].y = (long)floor(dy);


	s0 = (floor(dx)-dx+1)*(floor(dy)-dy+1);
	s1 = (floor(dx)-dx+1)*(dy-floor(dy));
	s2 = (dx-floor(dx))*(dy-floor(dy));
	s3 = (dx-floor(dx))*(floor(dy)-dy+1);

	// printf("s0=%f s1=%f s2=%f s3=%f\n",s0,s1,s2,s3);
}


// LAPLACE CALCULATIONS :
/*template<class ARG_TYPE>
double Table2D<ARG_TYPE>::CalcG54( LONG_T x, LONG_T y, LONG_T xSize, ARG_TYPE** p_data )
{
//	register int plusSum = (int)(p_data[y][x]+p_data[y-1][x]+p_data[y+1][x]+p_data[y][x-1]+p_data[y][x+1]);
//	register int minusSum = (int)(p_data[y-2][x-2]+p_data[y-2][x+2]+p_data[y+2][x-2]+p_data[y+2][x+2]);
//
//	double ret = plusSum - (1.25*minusSum);
	return ( ( p_data[y][x]+p_data[y-1][x]+p_data[y+1][x]+p_data[y][x-1]+p_data[y][x+1])
				-1.25*( p_data[y-2][x-2]+p_data[y-2][x+2]+p_data[y+2][x-2]+p_data[y+2][x+2] ) );
//	return ( ( p_data[y][x]+p_data[y-1][x]+p_data[y+1][x]+p_data[y][x-1]+p_data[y][x+1])
//				- (5.00*(p_data[y-2][x-2]+p_data[y-2][x+2]+p_data[y+2][x-2]+p_data[y+2][x+2]))/4.00  );
}*/


template<class ARG_TYPE>
LONG_T Table2D<ARG_TYPE>::CalcLaplaceSum( LONG_T x, LONG_T y, LONG_T xSize,
                                 ARG_TYPE** p_data, eLaplaceType_T laplaceType,
                                 LONG_T& plus_sum, LONG_T& minus_sum,
											BOOL_T bMedianBkg /*=FALSE*/ )
{
	plus_sum = 0;
	for(register int i=0;i<gLaplaceInfo[laplaceType].m_PlusCount;i++){
		plus_sum += (int)p_data[y+gLaplaceInfo[laplaceType].m_PlusList[i].y][x+gLaplaceInfo[laplaceType].m_PlusList[i].x];
	}

	LONG_T bkg_values[100];
	minus_sum=0;
	for(register int i=0;i<gLaplaceInfo[laplaceType].m_MinusCount;i++){
		bkg_values[i] = (int)p_data[y+gLaplaceInfo[laplaceType].m_MinusList[i].y][x+gLaplaceInfo[laplaceType].m_MinusList[i].x];
		minus_sum += bkg_values[i];
	}
	
	double bkg_value = gLaplaceInfo[laplaceType].m_MinusFactor*minus_sum;
	if( bMedianBkg ){
		my_qsort( bkg_values, gLaplaceInfo[laplaceType].m_MinusCount );
		bkg_value = bkg_values[(int)(gLaplaceInfo[laplaceType].m_MinusCount/2)];	
		bkg_value = bkg_value*gLaplaceInfo[laplaceType].m_PlusCount;
	}

	return (int)(plus_sum - bkg_value);
}

template<class ARG_TYPE>
LONG_T Table2D<ARG_TYPE>::CalcLaplaceSumBig( LONG_T x, LONG_T y, LONG_T xSize,
                                 BIG_ELEM_TYPE** p_data, eLaplaceType_T laplaceType,
                                 LONG_T& plus_sum, LONG_T& minus_sum )
{
	plus_sum = 0;
	for(register int i=0;i<gLaplaceInfo[laplaceType].m_PlusCount;i++){
		plus_sum += p_data[y+gLaplaceInfo[laplaceType].m_PlusList[i].y][x+gLaplaceInfo[laplaceType].m_PlusList[i].x];
	}
	minus_sum=0;
	for(register int i=0;i<gLaplaceInfo[laplaceType].m_MinusCount;i++){
		minus_sum += p_data[y+gLaplaceInfo[laplaceType].m_MinusList[i].y][x+gLaplaceInfo[laplaceType].m_MinusList[i].x];
	}

	return (int)(plus_sum - gLaplaceInfo[laplaceType].m_MinusFactor*minus_sum);
}


/*template<class ARG_TYPE>
LONG_T Table2D<ARG_TYPE>::CalcLaplaceSumOLD( LONG_T x, LONG_T y, LONG_T xSize, 
                                 ARG_TYPE** p_data, eLaplaceType_T laplaceType,
											LONG_T& plus_sum, LONG_T& minus_sum )
{
	register long x_1 = x-1;
	register long y_1 = y-1;
	register long xp1 = x+1;
	register long yp1 = y+1;
	register long x_2 = x-2;
	register long y_2 = y-2;
	register long xp2 = x+2;
	register long yp2 = y+2;

	int laplace = 0;
	plus_sum = 0;
	minus_sum = 0;

	if(laplaceType==eRawS){
		minus_sum = 0;
		plus_sum = p_data[y][x];		
		return plus_sum;
	}

	if(laplaceType==eSinglePoint || laplaceType==eFivePoints){
		minus_sum = (p_data[yp1][x_1] + p_data[yp1][xp1] + 
				       p_data[y_1][x_1] + p_data[y_1][xp1]);

		if(laplaceType==eSinglePoint){
			plus_sum = p_data[y][x];
			minus_sum = (minus_sum >> 2); // div by 4
		}else{
			plus_sum = p_data[y][x]+p_data[y][x_1]+p_data[y][xp1]+
			 			  p_data[y_1][x]+p_data[yp1][x];
			minus_sum = (int)((5*minus_sum)/4.00);
		}
	}else{
		if(laplaceType==eTwoPoints){
			plus_sum = p_data[y][x]+p_data[y][xp1];
			minus_sum = (p_data[yp1][x_1]+p_data[y_1][x_1]+
                      p_data[yp1][xp1+1]+p_data[y_1][xp1+1]);
			minus_sum = (minus_sum >> 1); // div by  2				
		}
		if(laplaceType==eFourPoints || laplaceType==eFourTwelve || laplaceType==eFourTwelveFar){
			plus_sum = p_data[y][x]   + p_data[y][xp1]+
                    p_data[y_1][x] + p_data[y_1][xp1];
			minus_sum = (p_data[yp1][x_1]+p_data[yp1][xp1+1]+
							 p_data[y_1-1][x_1]+p_data[y_1-1][xp1+1]);
			if(laplaceType==eFourTwelve){
				minus_sum += ( p_data[yp1][x]+p_data[yp1][xp1]+p_data[y_2][x]+p_data[y_2][xp1]+
 									p_data[y][x_1]+p_data[y_1][x_1]+p_data[y][xp2]+p_data[y_1][xp2] );
				minus_sum = (minus_sum)*0.75;
			}
			if(laplaceType==eFourTwelveFar){
				minus_sum += ( p_data[yp2][x]+p_data[yp2][xp1]+p_data[y_2-1][x]+p_data[y_2-1][xp1]+
									p_data[y][x_2]+p_data[y_1][x_2]+p_data[y][xp2+1]+p_data[y_1][xp2+1] );
				minus_sum = (minus_sum)*0.75;
			}
		}
		if(laplaceType==eFivePlusFourMin){
			plus_sum = p_data[y][x]+p_data[y][x_1]+p_data[y][xp1]+
			 			  p_data[y_1][x]+p_data[yp1][x];
			minus_sum = p_data[yp1+1][x_1-1]+p_data[yp1+1][xp1+1]+
							p_data[y_1-1][x_1-1]+p_data[y_1-1][xp1+1];
			minus_sum = (int)((5*minus_sum)/4.00);
			
		}
		if(laplaceType==eNineEightFive){
			plus_sum = p_data[y][x]+p_data[y][x_1]+p_data[y][xp1]+
						  p_data[y_1][x]+p_data[y_1][x_1]+p_data[y_1][xp1]+
						  p_data[yp1][x]+p_data[yp1][x_1]+p_data[yp1][xp1];
			minus_sum = p_data[y_1][x_1-1]+p_data[y_1][xp1+1]+
							p_data[y_1-1][x_1]+p_data[y_1-1][xp1]+
							p_data[yp1][x_1-1]+p_data[yp1][xp1+1]+
							p_data[yp1+1][x_1]+p_data[yp1+1][xp1];
			minus_sum = (long)((9.00*minus_sum)/8.00);
		}
		if(laplaceType==eNineEightSevenVeryBig){
			plus_sum = p_data[y][x]+p_data[y][x_1]+p_data[y][xp1]+
                    p_data[y_1][x]+p_data[y_1][x_1]+p_data[y_1][xp1]+
                    p_data[yp1][x]+p_data[yp1][x_1]+p_data[yp1][xp1];
			minus_sum = p_data[y][x_2-1]+p_data[y][xp2+1]+
							p_data[y_2][x_2]+p_data[y_2][xp2]+p_data[y_2-1][x]+
							p_data[yp2][x_2]+p_data[yp2][xp2]+p_data[yp2+1][x];
			minus_sum = (long)((9.00*minus_sum)/8.00);
		}
		if(laplaceType==eEightFour || laplaceType==eEightTen){
			plus_sum = p_data[y][x_1]+p_data[y][x]+p_data[y][xp1]+p_data[y][xp2]+
						  p_data[y_1][x]+p_data[y_1][xp1]+
						 p_data[yp1][x]+p_data[yp1][xp1];
			if(laplaceType==eEightFour){
				minus_sum = p_data[y_2][x_2]+p_data[y_2][xp2+1]+
							   p_data[yp2][x_2]+p_data[yp2][xp2+1];
				minus_sum = (minus_sum << 1); // times 2
			}
			if(laplaceType==eEightTen){
				minus_sum = p_data[y_2][x_2]+p_data[y_2][xp2+1]+
								p_data[yp2][x_2]+p_data[yp2][xp2+1]+
								p_data[y][x_2-1]+p_data[y][xp2+2]+
								p_data[yp2+1][x]+p_data[yp2+1][xp1]+
								p_data[y_2-1][x]+p_data[y_2-1][xp1];
				minus_sum = (long)(minus_sum*0.8);
			}
		}
		if(laplaceType==eFiveEight){
			plus_sum = p_data[y][x]+p_data[y][x_1]+p_data[y][xp1]+
			 			  p_data[y_1][x]+p_data[yp1][x];
			minus_sum = p_data[yp1+1][x_1-1]+p_data[yp1+1][xp1+1]+
							p_data[y_1-1][x_1-1]+p_data[y_1-1][xp1+1]+
							p_data[yp2+1][x]+p_data[y_2-1][x]+
							p_data[y][xp2+1]+p_data[y][x_2-1];
			minus_sum = (long)(minus_sum*FIVE_EIGHTS);
		}
	}


	laplace = (plus_sum - minus_sum);		
	return laplace;

}

template<class ARG_TYPE>
LONG_T Table2D<ARG_TYPE>::CalcLaplaceSumBigOLD( LONG_T x, LONG_T y, LONG_T xSize, 
                                 BIG_ELEM_TYPE** p_data, eLaplaceType_T laplaceType,
											LONG_T& plus_sum, LONG_T& minus_sum )
{
	register long x_1 = x-1;
	register long y_1 = y-1;
	register long xp1 = x+1;
	register long yp1 = y+1;
	register long x_2 = x-2;
	register long y_2 = y-2;
	register long xp2 = x+2;
	register long yp2 = y+2;

	int laplace = 0;
	plus_sum = 0;
	minus_sum = 0;

	if(laplaceType==eRawS){
		minus_sum = 0;
		plus_sum = p_data[y][x];		
		return plus_sum;
	}

	if(laplaceType==eSinglePoint || laplaceType==eFivePoints){
		minus_sum = (p_data[yp1][x_1] + p_data[yp1][xp1] + 
				       p_data[y_1][x_1] + p_data[y_1][xp1]);

		if(laplaceType==eSinglePoint){
			plus_sum = p_data[y][x];
			minus_sum = (minus_sum >> 2); // div by 4
		}else{
			plus_sum = p_data[y][x]+p_data[y][x_1]+p_data[y][xp1]+
			 			  p_data[y_1][x]+p_data[yp1][x];
			minus_sum = (int)((5*minus_sum)/4.00);
		}
	}else{
		if(laplaceType==eTwoPoints){
			plus_sum = p_data[y][x]+p_data[y][xp1];
			minus_sum = (p_data[yp1][x_1]+p_data[y_1][x_1]+
                      p_data[yp1][xp1+1]+p_data[y_1][xp1+1]);
			minus_sum = (minus_sum >> 1); // div by  2				
		}
		if(laplaceType==eFourPoints || laplaceType==eFourTwelve || laplaceType==eFourTwelveFar){
			plus_sum = p_data[y][x]   + p_data[y][xp1]+
                    p_data[y_1][x] + p_data[y_1][xp1];
			minus_sum = (p_data[yp1][x_1]+p_data[yp1][xp1+1]+
							 p_data[y_1-1][x_1]+p_data[y_1-1][xp1+1]);
			if(laplaceType==eFourTwelve){
				minus_sum += ( p_data[yp1][x]+p_data[yp1][xp1]+p_data[y_2][x]+p_data[y_2][xp1]+
 									p_data[y][x_1]+p_data[y_1][x_1]+p_data[y][xp2]+p_data[y_1][xp2] );
				minus_sum = (minus_sum)*0.75;
			}
			if(laplaceType==eFourTwelveFar){
				minus_sum += ( p_data[yp2][x]+p_data[yp2][xp1]+p_data[y_2-1][x]+p_data[y_2-1][xp1]+
									p_data[y][x_2]+p_data[y_1][x_2]+p_data[y][xp2+1]+p_data[y_1][xp2+1] );
				minus_sum = (minus_sum)*0.75;
			}
		}	
		if(laplaceType==eFivePlusFourMin){
			plus_sum = p_data[y][x]+p_data[y][x_1]+p_data[y][xp1]+
			 			  p_data[y_1][x]+p_data[yp1][x];
			minus_sum = p_data[yp1+1][x_1-1]+p_data[yp1+1][xp1+1]+
							p_data[y_1-1][x_1-1]+p_data[y_1-1][xp1+1];
			minus_sum = (int)((5*minus_sum)/4.00);
			
		}
		if(laplaceType==eNineEightFive){
			plus_sum = p_data[y][x]+p_data[y][x_1]+p_data[y][xp1]+
						  p_data[y_1][x]+p_data[y_1][x_1]+p_data[y_1][xp1]+
						  p_data[yp1][x]+p_data[yp1][x_1]+p_data[yp1][xp1];
			minus_sum = p_data[y_1][x_1-1]+p_data[y_1][xp1+1]+
							p_data[y_1-1][x_1]+p_data[y_1-1][xp1]+
							p_data[yp1][x_1-1]+p_data[yp1][xp1+1]+
							p_data[yp1+1][x_1]+p_data[yp1+1][xp1];
			minus_sum = (long)((9.00*minus_sum)/8.00);
		}
		if(laplaceType==eNineEightSevenVeryBig){
			plus_sum = p_data[y][x]+p_data[y][x_1]+p_data[y][xp1]+
                    p_data[y_1][x]+p_data[y_1][x_1]+p_data[y_1][xp1]+
                    p_data[yp1][x]+p_data[yp1][x_1]+p_data[yp1][xp1];
			minus_sum = p_data[y][x_2-1]+p_data[y][xp2+1]+
							p_data[y_2][x_2]+p_data[y_2][xp2]+p_data[y_2-1][x]+
							p_data[yp2][x_2]+p_data[yp2][xp2]+p_data[yp2+1][x];
			minus_sum = (long)((9.00*minus_sum)/8.00);
		}
		if(laplaceType==eEightFour || laplaceType==eEightTen){
			plus_sum = p_data[y][x_1]+p_data[y][x]+p_data[y][xp1]+p_data[y][xp2]+
						  p_data[y_1][x]+p_data[y_1][xp1]+
						 p_data[yp1][x]+p_data[yp1][xp1];
			if(laplaceType==eEightFour){
				minus_sum = p_data[y_2][x_2]+p_data[y_2][xp2+1]+
							   p_data[yp2][x_2]+p_data[yp2][xp2+1];
				minus_sum = (minus_sum << 1); // times 2
			}
			if(laplaceType==eEightTen){
				minus_sum = p_data[y_2][x_2]+p_data[y_2][xp2+1]+
								p_data[yp2][x_2]+p_data[yp2][xp2+1]+
								p_data[y][x_2-1]+p_data[y][xp2+2]+
								p_data[yp2+1][x]+p_data[yp2+1][xp1]+
								p_data[y_2-1][x]+p_data[y_2-1][xp1];
				minus_sum = (long)(minus_sum*0.8);
			}
		}
		if(laplaceType==eFiveEight){
			plus_sum = p_data[y][x]+p_data[y][x_1]+p_data[y][xp1]+
			 			  p_data[y_1][x]+p_data[yp1][x];
			minus_sum = p_data[yp1+1][x_1-1]+p_data[yp1+1][xp1+1]+
							p_data[y_1-1][x_1-1]+p_data[y_1-1][xp1+1]+
							p_data[yp2+1][x]+p_data[y_2-1][x]+
							p_data[y][xp2+1]+p_data[y][x_2-1];
			minus_sum = (long)(minus_sum*FIVE_EIGHTS);
		}		
	}


	laplace = (plus_sum - minus_sum);		

	return laplace;
}*/

//   ATTENTION :
//
// does not check if points around are still in frame - must 
// be ensured be calling function !!!!!!!!!!!!!!!!
//
template<class ARG_TYPE>
LONG_T Table2D<ARG_TYPE>::CalcLaplaceSum( LONG_T x, LONG_T y,
										 LONG_T xSize, ARG_TYPE** p_data,
										 eLaplaceType_T laplaceType,
										 BOOL_T bMedianBkg /*=FALSE*/ )
{
/*	LONG_T laplace = 0;
	LONG_T plus_sum = 0;
	LONG_T minus_sum = 0;
	laplace = CalcLaplaceSum( x , y, xSize, p_data, laplaceType, 
									  plus_sum, minus_sum );
	return laplace;*/

	register int plus_sum = 0;
	for(register int i=0;i<gLaplaceInfo[laplaceType].m_PlusCount;i++){
		plus_sum += (int)p_data[y+gLaplaceInfo[laplaceType].m_PlusList[i].y][x+gLaplaceInfo[laplaceType].m_PlusList[i].x];
	}

	LONG_T bkg_values[100];
	register int minus_sum=0;
	for(register int i=0;i<gLaplaceInfo[laplaceType].m_MinusCount;i++){
		bkg_values[i] = (int)p_data[y+gLaplaceInfo[laplaceType].m_MinusList[i].y][x+gLaplaceInfo[laplaceType].m_MinusList[i].x];
		minus_sum += bkg_values[i];
	}

	double bkg_value = gLaplaceInfo[laplaceType].m_MinusFactor*minus_sum;
	if( bMedianBkg ){
		my_qsort( bkg_values, gLaplaceInfo[laplaceType].m_MinusCount );
		bkg_value = bkg_values[(int)(gLaplaceInfo[laplaceType].m_MinusCount/2)];	
		bkg_value = bkg_value*gLaplaceInfo[laplaceType].m_PlusCount;
	}

	return (int)(plus_sum - bkg_value);
}

template<class ARG_TYPE>
LONG_T Table2D<ARG_TYPE>::CalcLaplaceSumBig( LONG_T x, LONG_T y,
										 LONG_T xSize, BIG_ELEM_TYPE** p_data,
										 eLaplaceType_T laplaceType )
{
/*	LONG_T laplace = 0;
	LONG_T plus_sum = 0;
	LONG_T minus_sum = 0;
	laplace = CalcLaplaceSumBig( x , y, xSize, p_data, laplaceType,
										  plus_sum, minus_sum );
	return laplace;*/

	register int plus_sum = 0;
	for(register int i=0;i<gLaplaceInfo[laplaceType].m_PlusCount;i++){
		plus_sum += (int)p_data[y+gLaplaceInfo[laplaceType].m_PlusList[i].y][x+gLaplaceInfo[laplaceType].m_PlusList[i].x];
	}
	register int minus_sum=0;
	for(register int i=0;i<gLaplaceInfo[laplaceType].m_MinusCount;i++){
		minus_sum += (int)p_data[y+gLaplaceInfo[laplaceType].m_MinusList[i].y][x+gLaplaceInfo[laplaceType].m_MinusList[i].x];
	}

	return (int)(plus_sum - gLaplaceInfo[laplaceType].m_MinusFactor*minus_sum);
}


template<class ARG_TYPE>
LONG_T Table2D<ARG_TYPE>::CalcLaplaceSumEstimateOnEdge( LONG_T x, LONG_T y, 
															   LONG_T xSize, LONG_T ySize,
                 											ARG_TYPE** p_data,
																CLongPoint* plus_list, LONG_T plus_count,
				                                    CLongPoint* minus_list, LONG_T minus_count )
{
	register int plus_cnt=0;
	register int minus_cnt = 0;
	register int plus_sum=0;
	register int minus_sum=0;

	for(register int p=0;p<plus_count;p++){
		register int xx = x+plus_list[p].x;
		register int yy = y+plus_list[p].y;
		
		if(xx>=0 && xx<xSize && yy>=0 && yy<ySize){
			plus_cnt++;
			plus_sum += (int)p_data[yy][xx];
		}
	}
	for(register int p=0;p<minus_count;p++){
		register int xx = x+minus_list[p].x;
		register int yy = y+minus_list[p].y;
		
		if(xx>=0 && xx<xSize && yy>=0 && yy<ySize){
			minus_cnt++;
			minus_sum += (int)p_data[yy][xx];
		}
	}
	register int laplace = (int)(plus_sum - (double(plus_cnt)/double(minus_cnt))*double(minus_sum));
	return laplace;	
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::GetPlusMinusList( eLaplaceType_T laplaceType,
											 CLongPoint* plus_list, LONG_T& plus_cnt,
                                  CLongPoint* minus_list, LONG_T& minus_cnt )
{	
	plus_cnt = 0;
	minus_cnt = 0;
	if(laplaceType==eRawS){
		plus_cnt = 1;
		minus_cnt = 0;
		plus_list[0].x = 0;
		plus_list[0].y = 0;
		return;
	}

	if(laplaceType==eSinglePoint){
		plus_cnt = 1;
		minus_cnt = 4;
		plus_list[0].x = 0;
		plus_list[0].y = 0;

		minus_list[0].x = -1;
		minus_list[0].y = -1;

		minus_list[1].x = 1;
		minus_list[1].y = -1;

		minus_list[2].x = -1;
		minus_list[2].y = 1;

		minus_list[3].x = 1;
		minus_list[3].y = 1;		
		return;
	}		
	if(laplaceType==eTwoPoints){
		plus_cnt = 2;
		minus_cnt = 4;

		plus_list[0].x = 0;
		plus_list[0].y = 0;

		plus_list[1].x = 1;
		plus_list[1].y = 0;

		minus_list[0].x = -1;
		minus_list[0].y = -1;

		minus_list[1].x = 2;
		minus_list[1].y = -1;

		minus_list[2].x = -1;
		minus_list[2].y = 1;

		minus_list[3].x = 2;
		minus_list[3].y = 1;		
		return;
	}		
	if(laplaceType==eFourPoints || laplaceType==eFourTwelve || laplaceType==eFourTwelveFar){
		plus_cnt = 4;
		minus_cnt = 4;

		plus_list[0].x = 0;
		plus_list[0].y = 0;

		plus_list[1].x = 1;
		plus_list[1].y = 0;

		plus_list[2].x = 0;
		plus_list[2].y = -1;

		plus_list[3].x = 1;
		plus_list[3].y = -1;

		minus_list[0].x = -1;
		minus_list[0].y = -2;

		minus_list[1].x = 2;
		minus_list[1].y = -2;

		minus_list[2].x = -1;
		minus_list[2].y = 1;

		minus_list[3].x = 2;
		minus_list[3].y = 1;		
		
		if(laplaceType==eFourTwelve){
			minus_list[4].x = 0;
	      minus_list[4].y = 1;

			minus_list[5].x = 1;
	      minus_list[5].y = 1;

			minus_list[6].x = 0;
	      minus_list[6].y = -2;

			minus_list[7].x = 1;
	      minus_list[7].y = -2;

			minus_list[8].x = -1;
	      minus_list[8].y = 0;

			minus_list[9].x = -1;
	      minus_list[9].y = -1;

			minus_list[10].x = +2;
	      minus_list[10].y = 0;

			minus_list[11].x = +2;
	      minus_list[11].y = -1;
		
			minus_cnt = 12;
		}
		if(laplaceType==eFourTwelveFar){
			minus_list[4].x = 0;
	      minus_list[4].y = 2;

			minus_list[5].x = 1;
	      minus_list[5].y = 2;

			minus_list[6].x = 0;
	      minus_list[6].y = -3;

			minus_list[7].x = 1;
	      minus_list[7].y = -3;

			minus_list[8].x = -2;
	      minus_list[8].y = 0;

			minus_list[9].x = -2;
	      minus_list[9].y = -1;

			minus_list[10].x = +3;
	      minus_list[10].y = 0;

			minus_list[11].x = +3;
	      minus_list[11].y = -1;

			minus_cnt = 12;
		}
		return;
	}		
	if(laplaceType==eFivePoints){
		plus_cnt = 5;
		minus_cnt = 4;

		plus_list[0].x = 0;
		plus_list[0].y = 0;

		plus_list[1].x = 1;
		plus_list[1].y = 0;

		plus_list[2].x = -1;
		plus_list[2].y = 0;

		plus_list[3].x = 0;
		plus_list[3].y = -1;

		plus_list[4].x = 0;
		plus_list[4].y = 1;

		minus_list[0].x = -1;
		minus_list[0].y = -1;

		minus_list[1].x = 1;
		minus_list[1].y = -1;

		minus_list[2].x = -1;
		minus_list[2].y = 1;

		minus_list[3].x = 1;
		minus_list[3].y = 1;		
		return;
	}		
	if(laplaceType==eFivePlusFourMin || laplaceType==eFiveEight){
		plus_cnt = 5;
		minus_cnt = 4;

		plus_list[0].x = 0;
		plus_list[0].y = 0;

		plus_list[1].x = 1;
		plus_list[1].y = 0;

		plus_list[2].x = -1;
		plus_list[2].y = 0;

		plus_list[3].x = 0;
		plus_list[3].y = -1;

		plus_list[4].x = 0;
		plus_list[4].y = 1;

		minus_list[0].x = -2;
		minus_list[0].y = -2;

		minus_list[1].x = 2;
		minus_list[1].y = -2;

		minus_list[2].x = -2;
		minus_list[2].y = 2;

		minus_list[3].x = 2;
		minus_list[3].y = 2;		

		if(laplaceType==eFiveEight){
			minus_cnt = 8;

			minus_list[4].x = 0;
	      minus_list[4].y = -3;

			minus_list[5].x = 0;
	      minus_list[5].y = 3;

			minus_list[6].x = -3;
	      minus_list[6].y = 0;

			minus_list[7].x = 3;
	      minus_list[7].y = 0;
		}
	}		
	if(laplaceType==eNineEightFive || laplaceType==eNineEightSevenVeryBig){
		plus_cnt = 9;
		minus_cnt = 8;
		int pos=0;
		for(register int y=-1;y<=1;y++){
			for(register int x=-1;x<=1;x++){
				plus_list[pos].y = y;
				plus_list[pos].x = x;
				pos++;
			}
		}

		if(laplaceType==eNineEightFive){
			minus_list[0].x = -1;
			minus_list[0].y = -2;

			minus_list[1].x = 1;
			minus_list[1].y = -2;
	
			minus_list[2].x = -2;
			minus_list[2].y = -1;

			minus_list[3].x = 2;
			minus_list[3].y = -1;

			minus_list[4].x = -1;
			minus_list[4].y = 2;

			minus_list[5].x = 1;
			minus_list[5].y = 2;

			minus_list[6].x = -2;
			minus_list[6].y = 1;

			minus_list[7].x = 2;
			minus_list[7].y = 1;
		}else{
			minus_list[0].x = -2;
			minus_list[0].y = -2;

			minus_list[1].x = 2;
			minus_list[1].y = -2;
	
			minus_list[2].x = -2;
			minus_list[2].y = 2;

			minus_list[3].x = 2;
			minus_list[3].y = 2;

			minus_list[4].x = 0;
			minus_list[4].y = -3;

			minus_list[5].x = 0;
			minus_list[5].y = 3;

			minus_list[6].x = -3;
			minus_list[6].y = 0;

			minus_list[7].x = 3;
			minus_list[7].y = 0;			
		}
	}
	if(laplaceType == eEightFour || laplaceType == eEightTen){
			plus_cnt = 8;
			minus_cnt = 4;

			plus_list[0].x = -1;
			plus_list[0].y = 0;

			plus_list[1].x = 0;
			plus_list[1].y = 0;

			plus_list[2].x = 1;
			plus_list[2].y = 0;

			plus_list[3].x = 2;
			plus_list[3].y = 0;

			plus_list[4].x = 0;
			plus_list[4].y = 1;

			plus_list[5].x = 1;
			plus_list[5].y = 1;

			plus_list[6].x = 0;
			plus_list[6].y = -1;

			plus_list[7].x = 1;
			plus_list[7].y = -1;
		
			minus_list[0].x = -2;
         minus_list[0].y = -2;		

			minus_list[1].x = 3;
         minus_list[1].y = -2;		

			minus_list[2].x = -2;
         minus_list[2].y = 2;		

			minus_list[3].x = 3;
         minus_list[3].y = 2;		
		
			if( laplaceType == eEightTen ){
				minus_list[4].x = 0;
	         minus_list[4].y = -3;

				minus_list[5].x = 1;
	         minus_list[5].y = -3;

				minus_list[6].x = 0;
	         minus_list[6].y = 3;

				minus_list[7].x = 1;
	         minus_list[7].y = 3;

				minus_list[8].x = -3;
	         minus_list[8].y = 0;

				minus_list[9].x = 4;
	         minus_list[9].y = 0;
				minus_cnt = 10;
			}
	}  
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::GetLaplacePlusMinusValues( ARG_TYPE** p_data,
														int* plusValues, int* minusValues,
														int& plusValCount,int& minusValCount,
														LONG_T x, LONG_T y,
														LONG_T xSize, LONG_T ySize,
														eLaplaceType_T laplaceType )
{
	plusValCount =0;
	for(register int i=0;i<gLaplaceInfo[laplaceType].m_PlusCount;i++){
		register int xx = x+gLaplaceInfo[laplaceType].m_PlusList[i].x;
		register int yy = y+gLaplaceInfo[laplaceType].m_PlusList[i].y;


		if(CheckRange(x,y,xSize,ySize)){
			plusValues[plusValCount] = (int)(p_data[yy][xx]);
			plusValCount++;
		}
	}

	minusValCount = 0;
	for(register int i=0;i<gLaplaceInfo[laplaceType].m_MinusCount;i++){
		register int yy = y+gLaplaceInfo[laplaceType].m_MinusList[i].y;
		register int xx = x+gLaplaceInfo[laplaceType].m_MinusList[i].x;


		if(CheckRange(x,y,xSize,ySize)){
			minusValues[minusValCount] = (int)p_data[yy][xx];
			minusValCount++;
		}
	}
	
}


template<class ARG_TYPE>
LONG_T Table2D<ARG_TYPE>::GetLaplacePlusMinusList( LONG_T* plusTab, LONG_T& plusCnt,
                                              LONG_T* minusTab, LONG_T& minusCnt,
															 LONG_T x, LONG_T y, LONG_T pos,
															 LONG_T xSize, eLaplaceType_T laplaceType )
{
	CLongPoint minus_list[15],plus_list[15];
	LONG_T plus_cnt,minus_cnt;
	GetPlusMinusList( laplaceType,  plus_list, plus_cnt, minus_list, minus_cnt );
	for(register int i=0;i<plus_cnt;i++){
		plusTab[i] = pos + plus_list[i].y*xSize + plus_list[i].x;
		minusTab[i] = pos + minus_list[i].y*xSize + minus_list[i].x;
	}
	return 0;
}


template<class ARG_TYPE>
const char* Table2D<ARG_TYPE>::GetLaplaceName( eLaplaceType_T laplaceType )
{
	if(laplaceType==eSinglePoint)
		return G1_NAME;
	if(laplaceType==eTwoPoints)
		return G2_NAME;
	if(laplaceType==eFourPoints)
		return G4_NAME;
	if(laplaceType==eFivePoints)
		return G5_NAME;
	if(laplaceType==eFivePlusFourMin)
		return G54_NAME;
	if(laplaceType==eNineEightFive)
		return G985_NAME;
	if(laplaceType==eNineEightSevenVeryBig)
		return G987_NAME;
	if(laplaceType==eEightFour)
		return G84_NAME;
	if(laplaceType==eEightTen)
		return G810_NAME;
	if(laplaceType==eFiveEight)
		return G58_NAME;
	if(laplaceType==eRawS)
		return S_NAME;
	if(laplaceType==eFourTwelve)
		return G412_NAME;
	if(laplaceType==eFourTwelveFar)
      return G412FAR_NAME;



	return "UNKNOWN";
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::Laplace( BIG_ELEM_TYPE** p_laplace_fast, eLaplaceType_T laplaceType,
											long x0, long y0, long x1, long y1 )
{
	for(register long y=y0;y<=y1;y++){
		for(register long x=x0;x<=x1;x++){
			p_laplace_fast[y][x] = CalcLaplaceSum( x, y, m_SizeX, m_pFastData, 
																laplaceType );
		}
	}		
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::LaplaceSafe( BIG_ELEM_TYPE** p_laplace_fast,
												 CLongPoint* LaplacePlusList, int LaplacePlusCount,
												 CLongPoint* LaplaceMinusList, int LaplaceMinusCount,
												 long x0, long y0, long x1, long y1 )
{
	for(register long y=y0;y<=y1;y++){
		for(register long x=x0;x<=x1;x++){
			// p_laplace_fast[y][x] = CBaseAnal::CalcLaplaceSumSafe( x, y, m_SizeX, m_pFastData );
			p_laplace_fast[y][x] = CalcLaplaceSumEstimateOnEdge( x, y,
                                            m_SizeX, m_SizeY,
                                            m_pFastData,
                                            LaplacePlusList, LaplacePlusCount,
                                            LaplaceMinusList, LaplaceMinusCount );

		}
	}
}


template<class ARG_TYPE>
void Table2D<ARG_TYPE>::Laplace( Table2D<ARG_TYPE>& out, eLaplaceType_T laplaceType, 
											int ignore_edge, 
											long x0, long y0, long x1, long y1 )
{
	out.Alloc( m_SizeX, m_SizeY );

	ARG_TYPE* p_data = m_pData;
	ARG_TYPE* p_laplace = out.m_pData;
	ARG_TYPE** p_laplace_fast = out.m_pFastData;



	long upY = m_SizeY-ignore_edge;
	long upX = m_SizeX-ignore_edge;
	long bottomX = ignore_edge;
	long bottomY = ignore_edge;

	if(x0>=0 && x0>=bottomX && x0<=upX){
		bottomX = x0;
	}
	if(x1>=0 && x1>=bottomX && x1<=upX){
		upX = x1;
	}
	if(y0>=0 && y0>=bottomY && y0<=upY){
		bottomY = y0;
	}
	if(y1>=0 && y1>=bottomY && y1<=upY){
		upY = y1;
	}
		

	for(register long y=bottomY;y<upY;y++){
		for(register long x=bottomX;x<upX;x++){
			long laplace = CalcLaplaceSum(	x, y, m_SizeX, m_pFastData, laplaceType );
			p_laplace_fast[y][x] = laplace;
		}
	}
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::FH()
{
	FlipImage( eReverseImageHor );
}
                            
template<class ARG_TYPE>                                                    
void Table2D<ARG_TYPE>::FV()
{
	FlipImage( eReverseImageVert );
}


template<class ARG_TYPE>
void Table2D<ARG_TYPE>::FlipImage( eDriverReverseImage_T flipType )
{
	if(flipType!=eReverseImageNone){
		if(flipType==eReverseImageHor){
			ARG_TYPE* buff = new ARG_TYPE[m_SizeX];
			for(register int y=0;y<m_SizeY;y++){
				for(register int x=0;x<m_SizeX;x++){
					buff[x] = m_pFastData[y][m_SizeX-x-1];
				}
				memcpy( (m_pData+(y*m_SizeX)), buff, sizeof(ARG_TYPE)*m_SizeX );			
			}			
			delete [] buff;
		}
		if(flipType==eReverseImageVert){
			ARG_TYPE* buff = new ARG_TYPE[m_SizeY];
			for(register int x=0;x<m_SizeX;x++){
				for(register int y=0;y<m_SizeY;y++){
					buff[y] = m_pFastData[m_SizeY-y-1][x];
				}
				for(register int i=0;i<m_SizeY;i++){
					m_pFastData[i][x] = buff[i];
				}
			}
			delete [] buff;
		}
	}
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::FlipImage( ARG_TYPE* data, int sizeX, int sizeY, 
											  eDriverReverseImage_T flipType )
{
	if(flipType!=eReverseImageNone){
		if(flipType==eReverseImageHor){
			ARG_TYPE* buff = new ARG_TYPE[sizeX];
			for(register int y=0;y<sizeY;y++){
				register int y_start = y*sizeX;
				register int y_end = (y+1)*sizeX-1;
				for(register int x=0;x<sizeX;x++){
					buff[x] = data[y_end-x];
				}
				memcpy( (data+y_start), buff, sizeof(ARG_TYPE)*sizeX );			
			}			
			delete [] buff;
		}
		if(flipType==eReverseImageVert){
			ARG_TYPE* buff = new ARG_TYPE[sizeY];
			register int size=sizeX*sizeY;
			register int sizeY_1 = (sizeY-1);
			for(register int x=0;x<sizeX;x++){
				int pos=x;
				int y=sizeY_1;
				while( pos<size && y>=0 ){
					buff[y] = data[pos];
					pos += sizeX;
					y--;
				} 
				pos=x;
				y=0;
				while( pos<size && y<sizeY ){
					data[pos] = buff[y];
					pos += sizeX;
					y++;
				}
			}
			delete [] buff;
		}
	}
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::ClearBorder( int border )
{
	int upX=(m_SizeX-border);
	int upY=(m_SizeY-border);
	for(register int y=0;y<m_SizeY;y++){
		for(register int x=0;x<m_SizeX;x++){
			if( x<border || y<border || x>=upX || y>=upY ){
				m_pFastData[y][x] = 0;
			}
		}
	}
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::Laplace( Table2D<BIG_ELEM_TYPE>* pFrameLaplace, 
											eLaplaceType_T laplaceType,
											BOOL_T bEstimateOnEdge/*=TRUE*/ )
{
	if(pFrameLaplace){
		pFrameLaplace->Alloc( m_SizeX, m_SizeY );

		ARG_TYPE* p_data = m_pData;
		BIG_ELEM_TYPE* p_laplace = pFrameLaplace->get_data_buffer();
		BIG_ELEM_TYPE** p_laplace_fast = pFrameLaplace->get_data_buffer_fast();

		register int ignore_edge = 5;
		register int upY = m_SizeY-ignore_edge;
		register int upX = m_SizeX-ignore_edge;

		BOOL_T bDone=FALSE;
		if( laplaceType==eFivePlusFourMin){
			for(register int y=ignore_edge;y<upY;y++){
            for(register int x=ignore_edge;x<upX;x++){
					// p_laplace_fast[y][x] = (m_pFastData[y][x]+m_pFastData[y-1][x]+m_pFastData[y+1][x]+m_pFastData[y][x-1]+m_pFastData[y][x+1])-1.25*(m_pFastData[y-2][x-2]+m_pFastData[y-2][x+2]+m_pFastData[y+2][x-2]+m_pFastData[y+2][x+2] );
					p_laplace_fast[y][x] = (int)CalcG54( x, y, m_SizeX, m_pFastData );
				}
			}
			bDone=TRUE;
		}

		if(!bDone){
			// printf("NOT OPTIMIZED ?\n");
			LONG_T plus_sum,minus_sum;
			for(register int y=ignore_edge;y<upY;y++){
				for(register int x=ignore_edge;x<upX;x++){
					// p_laplace_fast[y][x] = CCD_Analyser::CalcLaplaceSum( x, y, m_SizeX, m_pFastData );;
					p_laplace_fast[y][x] = CalcLaplaceSum( x, y, m_SizeX, m_pFastData, 
																	  laplaceType, plus_sum,minus_sum );
																				  
				}
			}
		}


		//if(bEstimateOnEdge){
		//	LaplaceEdgeOnly( ignore_edge );
		//}
	}
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::GenEmptyImage( double sigma, double mean )
{
	printf("Generating image using mean=%.2f, sigma=%.2f, (image size=%d)\n",mean,sigma,m_Size);fflush(stdout);
	CMyProgressBar bar(0,m_Size);
	for(int i=0;i<m_Size;i++){
		// without root use :
		// m_pData[i] = CRandom::GetGauss( sigma, mean );

		m_pData[i] = CRandom::GetFastGauss( sigma, mean );
		bar.SetValue(i);
		bar.Update();
	}
	printf("Image generated OK\n");fflush(stdout);
}

template<class ARG_TYPE>
void Table2D<ARG_TYPE>::CutBorderBase( int left, int up, int right, int bottom, Table2D<ARG_TYPE>& out_image )
{
	int size_x = m_SizeX - (left+right);
   int size_y = m_SizeY - (bottom+up); 

   out_image.InitConstructor( size_x, size_y );
   ARG_TYPE** out_data = out_image.get_data_buffer_fast();

   for(int y=bottom;y<(m_SizeY-up);y++){
      for(int x=left;x<(m_SizeX-right);x++){
         out_data[y-bottom][x-left] = m_pFastData[y][x];
      }
   }   
       
       
   char savearea[50];
   sprintf(savearea,"%d %d %d %d",left,bottom,(m_SizeX-right),(m_SizeY-up));
//   out_image.GetKeyTab().Set( SAVEAREA , savearea );
//   out_image.GetKeyTab().Set( ASTTHIS , 0 );
}

#endif



