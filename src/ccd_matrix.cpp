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
#include "ccd_matrix.h"
#include <cfg.h>
#include <cexcp.h>
#include <mydate.h>
#include "ccd_trace.h"
#include <fits_file.h>
#include <myfile.h>
#include <mykeytab.h>
#include <math.h>
#include <tab2Ddesc.h>
#include <myutil.h>
#include <mymacros.h>
#include <mycmnglobals.h>
#include "ccd_globals.h"
#include "ccd_analyse.h"
#include "ccd_procstate.h"
#include "ccd_pipeline.h"
#include <ccd_fits_header_defs.h>
#include <ccd_piman_interf_defs.h>
#include <fitslib_globals.h>
#include <myparser.h>
#include "ccd_util.h"
#include "ccd_log.h"
#include <Astroangle.h>
#include "ccd_asastransform.h"
#include <myfits.h>
#include <mathfunc.h>
#include <myprogress.h>
#include <AstroCCD.h>
#include <algorithm> // sort 

#include <errno.h>

// #include <vector>
// using namespace std;

BOOL_T CCDMatrix::CCDMatrix_InitConstructor(long x_size,long y_size,BOOL_T bAutoInit,LONG_T idx,
							cCCD* pFrame,int nAllocEvents,BOOL_T bKeepLaplaceOfCurrent,
							CCcdCfg* pParamsSet/*=NULL*/,BOOL_T bAllocHere/*=TRUE*/)
{
	BOOL_T bRet = Table2D<ELEM_TYPE>::InitConstructor(x_size,y_size,idx,bAllocHere);
	m_Average = 0;
	m_bIntersting = FALSE;
	m_pFrame = pFrame;
	m_FoundEventList.reserve(nAllocEvents);
	m_GeneratedEventList.reserve(nAllocEvents);
	m_pFrameLaplace = NULL;
	m_pMatrixParamsSet = pParamsSet;
	m_pHeaderPtr = NULL;
	m_dx = 0;
	m_dy = 0;

	if(bAutoInit){
		if (x_size==0 && y_size==0){
			 long size_x = gCCDParams.m_SizeX;
			 long size_y = gCCDParams.m_SizeY;
			 bRet = Alloc( size_x, size_y );
		}		
	}
	if(bKeepLaplaceOfCurrent){		
		m_pFrameLaplace = new Table2D<BIG_ELEM_TYPE>( m_SizeX, m_SizeY );
	}

	return bRet;
}

CCDMatrix::CCDMatrix()
: Table2D<ELEM_TYPE>(), m_Average(0), m_bIntersting(FALSE),
	m_pFrame(NULL),
	m_pFrameLaplace(NULL), m_pMatrixParamsSet( NULL ), m_pHeaderPtr(NULL),
	m_dx(0),m_dy(0)
{

}

CCDMatrix::CCDMatrix(long x_size,long y_size,BOOL_T bAutoInit,LONG_T idx,
                     cCCD* pFrame,int nAllocEvents,BOOL_T bKeepLaplaceOfCurrent,
							CCcdCfg* pParamsSet/*=NULL*/,BOOL_T bAllocHere/*=TRUE*/)
: Table2D<ELEM_TYPE>(x_size,y_size,idx,bAllocHere), m_Average(0), m_bIntersting(FALSE),
	m_pFrame(pFrame),m_FoundEventList(nAllocEvents),m_GeneratedEventList(nAllocEvents),
	m_pFrameLaplace(NULL), m_pMatrixParamsSet( pParamsSet ), m_pHeaderPtr(NULL),
	m_dx(0),m_dy(0)
{
	if(bAutoInit){
		if (x_size==0 && y_size==0){
			 long size_x = gCCDParams.m_SizeX;
			 long size_y = gCCDParams.m_SizeY;
			 Alloc( size_x, size_y );
		}		
	}
	if(bKeepLaplaceOfCurrent){		
		m_pFrameLaplace = new Table2D<BIG_ELEM_TYPE>( m_SizeX, m_SizeY );
	}
}

CCDMatrix::~CCDMatrix()
{
	if(m_pFrameLaplace)
		delete m_pFrameLaplace;
}


BIG_ELEM_TYPE* CCDMatrix::get_frame_laplace()
{
	if(m_pFrameLaplace)
		return m_pFrameLaplace->get_data_buffer();
	return NULL;
}

BIG_ELEM_TYPE** CCDMatrix::get_frame_laplace_fast()
{
	if(m_pFrameLaplace)
		return m_pFrameLaplace->get_data_buffer_fast();
	return NULL;
}

void CCDMatrix::UpdateGenEventsByArea( int low_x, int up_x, int low_y, int up_y )
{
	CCDEventList::iterator i;
	for(i=m_GeneratedEventList.begin();i!=m_GeneratedEventList.end();i++){
		if( i->m_MaxPoint.x<low_x || i->m_MaxPoint.y<low_y ||
			 i->m_MaxPoint.x>=up_x || i->m_MaxPoint.y>=up_y ){
			(i->m_PixelAnalResults).m_bOutSideAcceptanceRegion = TRUE;
		}
	}
}

BOOL_T CCDMatrix::IsGenerated( int x, int y )
{
	CCDEventList::iterator i;
	for(i=m_GeneratedEventList.begin();i!=m_GeneratedEventList.end();i++){
		if(fabs(i->m_MaxPoint.x-x)<=3 && fabs(i->m_MaxPoint.y-y)<=3){
			return TRUE;
		}
	}
	return FALSE;
}

void CCDMatrix::UpdateGenEvents( CPixelAnalyseIn& in, CPixelAnalyseOut& out)
{
	CCDEventList::iterator i;
	for(i=m_GeneratedEventList.begin();i!=m_GeneratedEventList.end();i++){
		register int new_x = in.x;
		register int new_y = in.y;
		if(out.cluster_cnt){
			new_x = out.m_PixelOut.x0;
			new_y = out.m_PixelOut.y0;
		}	
		
		long gen_start_x = (long)i->m_MaxPoint.x-1;
		long gen_start_y = (long)i->m_MaxPoint.y-1;
		long gen_end_x = (long)i->m_MaxPoint.x+1;
		long gen_end_y = (long)i->m_MaxPoint.y+1;
		if(new_x >= gen_start_x && new_x<=gen_end_x &&
     	   new_y >= gen_start_y && new_y<=gen_end_y){
			// printf("newSum=%d\n",out.m_PixelOut.newSum);
			
			if((out.m_PixelOut.newSum>(i->m_PixelAnalResults).newSum) || (out.cluster_cnt>i->m_ClusterCount)){
				// filling only first time !
				// i->m_MaxPoint.x = new_x;
				// i->m_MaxPoint.y = new_y;
				MYTRACE4(gCCDTrace,"gen event rejection at (" << new_x << "," 
                        << new_y << ") : newSum=" << out.m_PixelOut.newSum 
								<< " AND otherSum=" << out.m_PixelOut.otherSum << " newClusterSum="
							   << out.newClusterSum << " and prev_max_cluster_sum=" 
                        << my_find_max_value_long(out.prevClusterSum, out.prevFramesCount) );
				i->SetIdentificationData( in, out );
			}
		}
	}
}


CccdReport* CCDMatrix::AddGeneratedEvent( long start_x, long start_y, 
                                   long end_x, long end_y, const char* mag,
                                   int otherMaxX, int otherMaxY, int frame_index,
											  int& OutMaxX, int& OutMaxY,
											  eEventType event_type /*=eFlash*/,
											  int max_lap )
{
	CPoint LowLeft(start_x,start_y);
	CPoint TopRight(end_x,end_y);
	LONG_T x_s = (start_x+end_x)/2;
	LONG_T y_s = (start_y+end_y)/2;

	// generating events - not suposed to be multi-threaded :
	static CPixelAnalyseIn in;
	static CPixelAnalyseOut out;

	in.x = x_s;
	in.y = y_s;
	in.frame_index = frame_index;

	out.m_PixelOut.x0 = x_s;
	out.m_PixelOut.y0 = y_s;

	CccdReport tmp(m_Index,	in, out, TRUE, FALSE, &LowLeft,&TopRight,mag);

	long maxX,maxY;
	// tmp.m_MaxValue = GetMaxValueAndPos( start_x,start_y, end_x, end_y,
	//				  							   maxX, maxY );
	// _TRACE_PRINTF_5("Other Image : (%d,%d)\n",OtherImage.m_MaxX,OtherImage.m_MaxY);

	tmp.m_MaxPoint.x = start_x + otherMaxX;
   tmp.m_MaxPoint.y = start_y + otherMaxY;
	tmp.m_MaxValue = m_pFastData[(long)tmp.m_MaxPoint.y][(long)tmp.m_MaxPoint.x];
	
	tmp.m_PixelAnalResults.eventType = event_type;
	tmp.m_PixelAnalResults.laplaceSum = max_lap;

	tmp.m_Time = getObsTime( (BOOL_T)FALSE );
	if(tmp.m_Time<=0)
		tmp.m_Time = get_dttm();

	_TRACE_PRINTF_4("\n\n PUT IMAGE at (%d,%d) MAX:(max_x,max_y) = (%f,%f)\n\n",start_x,start_y,tmp.m_MaxPoint.x,tmp.m_MaxPoint.y);

	OutMaxX = (int)tmp.m_MaxPoint.x;
	OutMaxY = (int)tmp.m_MaxPoint.y;
	
	m_GeneratedEventList.Add( tmp );
	return &(m_GeneratedEventList[ m_GeneratedEventList.size()-1 ]);
}


void CCDMatrix::AddGeneratedEvent( long x_s, long y_s,
                                   const char* mag, int frame_index,
											  eEventType event_type /*=eFlash*/)
{
	CPoint LowLeft( x_s-5 , y_s-5);
   CPoint TopRight( x_s+5, y_s+5 );

	// generating events - not suposed to be multi-threaded :
	static CPixelAnalyseIn in;
	static CPixelAnalyseOut out;

	in.x = x_s;
	in.y = y_s;
	in.frame_index = frame_index;

	out.m_PixelOut.x0 = x_s;
	out.m_PixelOut.y0 = y_s;

	CccdReport tmp(m_Index,	in, out, TRUE, FALSE, &LowLeft,&TopRight,mag);

	long maxX,maxY;

	tmp.m_MaxPoint.x = x_s;
   tmp.m_MaxPoint.y = y_s;
	tmp.m_MaxValue = m_pFastData[(long)tmp.m_MaxPoint.y][(long)tmp.m_MaxPoint.x];
	
	tmp.m_PixelAnalResults.eventType = event_type;

	_TRACE_PRINTF_4("\n\n PUT IMAGE at MAX:(max_x,max_y) = (%f,%f)\n\n",tmp.m_MaxPoint.x,tmp.m_MaxPoint.y);
	
	m_GeneratedEventList.Add( tmp );
}



void CCDMatrix::CalcAdditionalInfo( ELEM_TYPE* p_data, long x, long y, long pos,
											   CPixelAnalyseOut& out,
												CccdReport& event )
{
	long pixel_cnt=0;
	double val = CCD_Analyser::CalcClusterRatio( p_data, pos, out.cluster, out.cluster_cnt, pixel_cnt,
																FALSE, 0, TRUE );		
	event.m_AdditionalInfo[CccdReport::eClusterChar] = val;

	val = CCD_Analyser::CalcClusterRatio( p_data, pos, out.neighb_list,  out.ncnt, pixel_cnt );		
	event.m_AdditionalInfo[CccdReport::eCrossChar] = val;

	val = CCD_Analyser::CalcClusterRatio( p_data, pos, out.neighb_list,  out.ncnt, 
													  pixel_cnt, TRUE, gCCDParams.m_ConfTresholdPerPixel);
	event.m_AdditionalInfo[CccdReport::eCrossAboveTreshChar] = val;

	event.m_nAddInfo=3;

	long max_val = -10000,max_pos=pos;	
	/*for(register long i=0;i<out.cluster_cnt;i++){
		if(m_pData[out.cluster[i]]>max_val){
			max_val = m_pData[out.cluster[i]];
			max_pos = out.cluster[i];
		}	
	}*/
	/*if(max_pos>=0){
		event.m_MaxPoint.Set( (max_pos%m_SizeX), (max_pos/m_SizeX) );
		event.m_MaxValue = m_pData[max_pos];
	}*/
	if( out.m_PixelOut.m_MaxClusterPixel.x>0 && out.m_PixelOut.m_MaxClusterPixel.y>0 ){
		max_pos=(int)(out.m_PixelOut.m_MaxClusterPixel.y*m_SizeX+out.m_PixelOut.m_MaxClusterPixel.x);		
		if( !gCCDParams.m_bUseOriginalXY ){ // NEW 20041102
			event.m_MaxPoint = out.m_PixelOut.m_MaxClusterPixel;
		}
	}
	event.m_MaxValue = m_pData[ max_pos ];
}

void CCDMatrix::AddFoundEvent( CccdReport& tmp )
{
	m_FoundEventList.push_back( tmp );
}

void CCDMatrix::SetEventTime()
{
	time_t ut_time = getObsTime( TRUE );
	for(int i=0;i<m_FoundEventList.size();i++){
		m_FoundEventList[i].m_Time = ut_time;
	}
}

CccdReport& CCDMatrix::AddFoundEvent( CCDEventList& event_list,
												  CPixelAnalyseIn& in, CPixelAnalyseOut& out, 
												  BOOL_T bCheckUnique )
{
	CccdReport tmp( m_Index, in, out, FALSE, TRUE, NULL, NULL, 0, TRUE );
	tmp.m_PipelineIndex = (in.pPipeline)->GetPipelineIndex();

	CalcAdditionalInfo( in.p_data, in.x, in.y, in.pos, out, tmp );	

	
	// temporary :
	memcpy(out.m_PixelOut.PrevFramesX,in.PrevFramesX,sizeof(int)*(gCCDParams.m_nMaxOfAverageOfPrevN+1));
	memcpy(out.m_PixelOut.PrevFramesY,in.PrevFramesY,sizeof(int)*(gCCDParams.m_nMaxOfAverageOfPrevN+1));

	tmp.SetIdentificationData( in, out );
	
	if( bCheckUnique ){
		CccdReport* pEvt = event_list.FindEventByMaxPoint( (tmp.m_MaxPoint).x, 
																			(tmp.m_MaxPoint).y, 1 );
		if ( pEvt ){
			return (*pEvt);
		}
	}

	event_list.push_back( tmp );
	// delete tmp;
	return event_list.back();	
}

CccdReport& CCDMatrix::AddFoundEvent( CPixelAnalyseIn& in, CPixelAnalyseOut& out )
{
	return AddFoundEvent( m_FoundEventList, in, out );	
}


void CCDMatrix::AddFoundEvent( long x, long y, CLongList& cluster)
{
	CLongList::iterator i;
	LONG_T tmp_cluster[MAX_CLUSTER_SIZE];
   long cnt=0;
	for(i=cluster.begin();i!=cluster.end();i++){
   	tmp_cluster[cnt] = (*i);
		cnt++;
   }

	CPixelAnalyseIn in;
   CPixelAnalyseOut out;

   in.x = x;
   in.y = y;

   out.m_PixelOut.x0 = x;
   out.m_PixelOut.y0 = y;

	CccdReport tmp(m_Index, in, out, FALSE,TRUE,NULL,NULL,0 );
	m_FoundEventList.push_back( tmp );
}


int CCDMatrix::CountFoundEvents( int x, int y, int radius )
{
	int ret=0;
	CCDEventList::iterator i;
	for(i=m_FoundEventList.begin();i!=m_FoundEventList.end();i++){
		if( (fabs(i->m_MaxPoint.x-x)<radius && fabs(i->m_MaxPoint.y-y)<radius) ){
			ret++;
		}
	}	
	return ret;
}

int CCDMatrix::CountOverlaps( long x, long y,
                               LONG_T* cluster, LONG_T cluster_cnt,
                               LONG_T& prev_evt_x, LONG_T& prev_evt_y )
{
		int ret=0;
		// assumes analysis goes from lowest y to higer (0 --> SizeY)
		register int size = m_FoundEventList.size();
		for(register int i=size-1;i>=0;i--){
			if(abs((long)m_FoundEventList[i].m_Point.x-x)<=gCCDParams.m_OverlapRedial &&
  				abs((long)m_FoundEventList[i].m_Point.y-y)<=gCCDParams.m_OverlapRedial ){
				prev_evt_x = (long)m_FoundEventList[i].m_Point.x;
				prev_evt_y = (long)m_FoundEventList[i].m_Point.y;
				ret++;
			}else{
				if(y-m_FoundEventList[i].m_Point.y>=gCCDParams.m_OverlapRedial){
					// if already at y below overlap redial - no need to continue !
					// because not overlaping FOR SURE
					return ret;
				}
			}
		}
		return ret;
}

BOOL_T CCDMatrix::CheckIfOverlapsFast( long x, long y, 
                               LONG_T* cluster, LONG_T cluster_cnt,
                               LONG_T& prev_evt_x, LONG_T& prev_evt_y )
{
		BOOL_T bRet = CheckIfOverlapsFast( m_FoundEventList, x, y, cluster,
													  cluster_cnt, prev_evt_x, prev_evt_y );
		return bRet;
}

BOOL_T CCDMatrix::CheckIfOverlapsFast( CCDEventList& found_events, 
										 long x, long y, 
                               LONG_T* cluster, LONG_T cluster_cnt,
                               LONG_T& prev_evt_x, LONG_T& prev_evt_y )
{
		// assumes analysis goes from lowest y to higer (0 --> SizeY)
		register int size = found_events.size();
		register int size_1 = (size-1);
		for(register int i=size_1;i>=0;i--){
			if(abs((long)found_events[i].m_Point.x-x)<=gCCDParams.m_OverlapRedial &&
  				abs((long)found_events[i].m_Point.y-y)<=gCCDParams.m_OverlapRedial ){
				prev_evt_x = (long)found_events[i].m_Point.x;
				prev_evt_y = (long)found_events[i].m_Point.y;
				return TRUE;			
			}else{
				if(y-found_events[i].m_Point.y>=gCCDParams.m_OverlapRedial){
					// if already at y below overlap redial - no need to continue !
					// because not overlaping FOR SURE
					return FALSE;
				}
			}
		}

		return FALSE;
}

BOOL_T CCDMatrix::CheckIfOverlapsFastMaxPoint( CCDEventList& found_events, 
										 long x, long y, 
                               LONG_T* cluster, LONG_T cluster_cnt,
                               LONG_T& prev_evt_x, LONG_T& prev_evt_y )
{
		// assumes analysis goes from lowest y to higer (0 --> SizeY)
		register int size = found_events.size();
		register int size_1 = (size-1);
		for(register int i=size_1;i>=0;i--){
			if(abs((long)found_events[i].m_MaxPoint.x-x)<=gCCDParams.m_OverlapRedial &&
  				abs((long)found_events[i].m_MaxPoint.y-y)<=gCCDParams.m_OverlapRedial ){
				prev_evt_x = (long)found_events[i].m_MaxPoint.x;
				prev_evt_y = (long)found_events[i].m_MaxPoint.y;
				return TRUE;			
			}else{
				if(y-found_events[i].m_MaxPoint.y>=gCCDParams.m_OverlapRedial){
					// if already at y below overlap redial - no need to continue !
					// because not overlaping FOR SURE
					return FALSE;
				}
			}
		}

		return FALSE;
}


BOOL_T CCDMatrix::CheckIfOverlaps( LONG_T* cluster, LONG_T cluster_cnt,
                  			        LONG_T& prev_evt_x, LONG_T& prev_evt_y )
{
		CCDEventList::iterator pFoundEvt;				
		for(pFoundEvt=m_FoundEventList.begin();pFoundEvt!=m_FoundEventList.end();pFoundEvt++){
			//if(my_table_overlap_check( pFoundEvt->m_Cluster, pFoundEvt->m_ClusterCount,
         //                           cluster, cluster_cnt ))
			//	return TRUE;

			// cluster is already sorted :
			if(my_sorted_table_overlap_check( pFoundEvt->m_Cluster, pFoundEvt->m_ClusterCount,
                      		                cluster, cluster_cnt )){
				prev_evt_x = (LONG_T)pFoundEvt->m_Point.x;
				prev_evt_y = (LONG_T)pFoundEvt->m_Point.y;
				return TRUE;
			}
			
		}
		return FALSE;
}


void CCDMatrix::Divide(CCDMatrix& right,CCDMatrix& result,int max_val)
{
	Assert(m_SizeX==right.m_SizeX,"X - sizes must be equal");
	Assert(m_SizeY==right.m_SizeY,"Y - sizes must be equal");
	result.Alloc(m_SizeX,m_SizeY);

	clock_t t1 = clock();
	register int size = (m_SizeX*m_SizeY);
	for(register int i=0;i<size;i++){
		int val=0;
		if(right.m_pData[i]!=0)
			val = m_pData[i] / right.m_pData[i];

		// in case limit of data type exceeded, setmax value :
		if(val>max_val)
			val = max_val;

		result.m_pData[i] = val;
	}
	
	clock_t t2 = clock();
	mystring msg = get_clock_in_sec_string( t2-t1 );

	MYTRACE1(gCCDTrace,"Subtraction of images took : " << msg );
}


int CCDMatrix::Subtract(CCDMatrix& right,BOOL_T bZero)
{
	return Subtract( right, *this, bZero );
}


int CCDMatrix::Subtract(CCDMatrix& right,CCDMatrix& result,BOOL_T bZero)
{
	Assert(m_SizeX==right.m_SizeX,"X - sizes must be equal");
	Assert(m_SizeY==right.m_SizeY,"Y - sizes must be equal");
	result.Alloc(m_SizeX,m_SizeY);

	clock_t t1 = clock();
	register int size = (m_SizeX*m_SizeY);
	int ret=0;
	for(register int i=0;i<size;i++){
		register int original = m_pData[i];
		result.m_pData[i] = original - right.m_pData[i];
		if(bZero){
			if(original < right.m_pData[i])
				result.m_pData[i] = 0;
				// Assert(FALSE,"Error in analysis value in pixel <0");
		}
		if( result.m_pData[i] > 0 )		
			ret++;
	}
	
	clock_t t2 = clock();
	mystring msg = get_clock_in_sec_string( t2-t1 );

	MYTRACE1(gCCDTrace,"Subtraction of images took : " << msg );
	
	return ret;
}

void CCDMatrix::Normalize()
{
	long sum=0;
	long size = (m_SizeX*m_SizeY);
	for(int i=0;i<size;i++){
		sum += m_pData[i];
	}	
	m_Average = (sum/size);
	if(m_Average!=0){
		for(int i=0;i<size;i++){
			m_pData[i] = (ELEM_TYPE)(((double)m_pData[i])/m_Average);
		}	
	}
}

void CCDMatrix::AddImage( long x0, long y0,
                          Table2D<ELEM_TYPE>& OtherImage,
                          long start_x,long start_y,long end_x,long end_y,
                          BOOL_T bGen,const char* mag,int frame_index,
								  int* maxX/*=NULL*/, int* maxY/*=NULL*/ )
{
	PutImage(x0,y0,OtherImage,start_x,start_y,end_x,end_y,TRUE,bGen,mag,frame_index,maxX,maxY);
}

void CCDMatrix::PutImage( long x0, long y0,
                			  Table2D<ELEM_TYPE>& OtherImage,
                          long start_x,long start_y,long end_x,long end_y,
                          BOOL_T bAdd,BOOL_T bGen,const char* mag,int frame_index,
								  int* maxX/*=NULL*/, int* maxY/*=NULL*/)
{
	if(end_x==-1)
		end_x = OtherImage.GetXSize() - 1;
	if(end_y==-1)
		end_y = OtherImage.GetYSize() - 1;

	long x_length = (end_x-start_x)+1; // x - size of image to put
	long y_length = (end_y-start_y)+1; // y - size of image to put
	long x0_in = x0;
	long y0_in = y0; 

	if(x0<-x_length || x0>=m_SizeX)
		return;
	if(y0<-y_length || y0>=m_SizeY)
		return;

	if(x0<0){
		x0_in = 0;
		start_x = start_x - x0;
	}
	if(y0<0){
		y0_in = 0; 
		start_y = start_y - y0;
	}
	x_length = (end_x-start_x)+1;
	y_length = (end_y-start_y)+1;
	
	if(x0_in+x_length-1>=m_SizeX){
		x_length = m_SizeX - x0_in;
	}
	if(y0_in+x_length-1>=m_SizeY){
		y_length = m_SizeY - y0_in;
	}

	// Assert(X_of_left_down_corner+x_length-1<m_SizeX,"Image to be pasted too big in X direction"); 	
	// Assert(Y_of_left_down_corner+y_length-1<m_SizeY,"Image to be pasted too big in Y direction"); 	



	long x1 = x0_in + x_length - 1;
	long y1 = y0_in + y_length - 1;
		
	
	ELEM_TYPE* left = GetPos(x0_in,y0_in);
	ELEM_TYPE* right = OtherImage.GetPos(start_x,start_y);


	for(int i=0;i<y_length;i++){
		if(!bAdd){
			for(register int kk=0;kk<x_length;kk++){
				register int newVal = right[kk];
            if(newVal<0)
               newVal = 0;

            left[kk] = newVal;
			}
		}else{
			for(register int kk=0;kk<x_length;kk++){
				register int newVal = left[kk] + right[kk];
            if(newVal<0)
               newVal = 0;

            left[kk] = newVal;
			}
		}
		left = left + m_SizeX;
		right = right + OtherImage.GetXSize();
	}		
	
	int ret_max_x,ret_max_y;
	if(bGen)
		AddGeneratedEvent( x0_in , y0_in , x1 , y1, mag, OtherImage.m_MaxX, 
								 OtherImage.m_MaxY, frame_index, ret_max_x, ret_max_y );

	int maxVal = m_pFastData[ret_max_y][ret_max_x];	
	for(register int yy=y0_in;yy<y_length;yy++){
		for(register int xx=x0_in;xx<x_length;xx++){
			if(m_pFastData[yy][xx]>maxVal){
				maxVal = m_pFastData[yy][xx];
				ret_max_x = xx;
				ret_max_y = yy;
			}
		}
	}

	if(maxX)
		(*maxX) = ret_max_x;
	if(maxY)
		(*maxY) = ret_max_y;

	
	if(m_pFrameLaplace){
		register int edgeSize = gCCDParams.GetEdgeSizeInLaplace();

		int x_rec_start = MAX(edgeSize,x0_in);
		int x_rec_end = MIN((m_SizeX-edgeSize-1),x1);
		int y_rec_start = MAX(edgeSize,y0_in);
      int y_rec_end = MIN((m_SizeY-edgeSize-1),y1);

		Laplace( x_rec_start , y_rec_start , x_rec_end , y_rec_end );
	}
	MYTRACE3(gCCDTrace,"Image Inserted at position : x=" << x0_in << ",y=" << y0_in);
}

CccdReport* CCDMatrix::AddSample( long x0, long y0,
                          Table2D<ELEM_SAMPLE_TYPE>& OtherImage,
                          long start_x,long start_y,long end_x,long end_y,
                          BOOL_T bGen,const char* mag,int frame_index,
								  int* maxX/*=NULL*/, int* maxY/*=NULL*/ )
{
	return PutSample(x0,y0,OtherImage,start_x,start_y,end_x,end_y,TRUE,bGen,mag,frame_index,maxX,maxY);
}

CccdReport* CCDMatrix::PutSample( long x0, long y0,
                			  Table2D<ELEM_SAMPLE_TYPE>& OtherImage,
                          long start_x,long start_y,long end_x,long end_y,
                          BOOL_T bAdd,BOOL_T bGen,const char* mag,int frame_index,
								  int* maxX/*=NULL*/, int* maxY/*=NULL*/)
{
	// move so that max point of sample is in :
	// OtherImage.m_MaxX , OtherImage.m_MaxY :
	x0 = (x0 - OtherImage.m_MaxX);
	y0 = (y0 - OtherImage.m_MaxY);
	if( x0<0 ) x0 = 0;
	if( y0<0 ) y0 = 0;

	if(end_x==-1)
		end_x = OtherImage.GetXSize() - 1;
	if(end_y==-1)
		end_y = OtherImage.GetYSize() - 1;

	long x_length = (end_x-start_x)+1; // x - size of image to put
	long y_length = (end_y-start_y)+1; // y - size of image to put
	long x0_in = x0;
	long y0_in = y0; 

	if(x0<-x_length || x0>=m_SizeX)
		return NULL;
	if(y0<-y_length || y0>=m_SizeY)
		return NULL;

	if(x0<0){
		x0_in = 0;
		int tmp_start_x = start_x - x0;
		if(tmp_start_x<0)
			tmp_start_x = 0;
		if(tmp_start_x>=OtherImage.GetXSize())
			tmp_start_x = start_x;
		start_x = tmp_start_x;
	}
	if(y0<0){
		y0_in = 0; 
		int tmp_start_y = start_y - y0;
		if(tmp_start_y<0)
			tmp_start_y = 0;
		if(tmp_start_y>=OtherImage.GetYSize())
			tmp_start_y = start_y;
		start_y = tmp_start_y;
	}
	
	

	x_length = (end_x-start_x)+1;
	y_length = (end_y-start_y)+1;
	
	if(x0_in+x_length-1>=m_SizeX){
		x_length = m_SizeX - x0_in;
	}
	if(y0_in+x_length-1>=m_SizeY){
		y_length = m_SizeY - y0_in;
	}

	// Assert(X_of_left_down_corner+x_length-1<m_SizeX,"Image to be pasted too big in X direction"); 	
	// Assert(Y_of_left_down_corner+y_length-1<m_SizeY,"Image to be pasted too big in Y direction"); 	



	long x1 = x0_in + x_length - 1;
	long y1 = y0_in + y_length - 1;
		
	
	ELEM_TYPE* left = GetPos(x0_in,y0_in);
	ELEM_SAMPLE_TYPE* right = OtherImage.GetPos(start_x,start_y);


	for(int i=0;i<y_length;i++){
		if(!bAdd){
			for(register int kk=0;kk<x_length;kk++){
				register int newVal = right[kk];
				if(newVal<0)
					newVal = 0;

				left[kk] = newVal;
			}
		}else{
			for(register int kk=0;kk<x_length;kk++){
				register int newVal = left[kk] + right[kk];
            if(newVal<0)
               newVal = 0;				

				left[kk] = newVal;
			}
		}
		left = left + m_SizeX;
		right = right + OtherImage.GetXSize();
	}		

	int max_lap=0;
	if(m_pFrameLaplace){
		register int edgeSize = gCCDParams.GetEdgeSizeInLaplace();

		int x_rec_start = MAX(edgeSize,x0_in);
		int x_rec_end = MIN((m_SizeX-edgeSize-1),x1);
		int y_rec_start = MAX(edgeSize,y0_in);
      int y_rec_end = MIN((m_SizeY-edgeSize-1),y1);

		Laplace( x_rec_start , y_rec_start , x_rec_end , y_rec_end );

		// long max_x,max_y;
		// max_lap = m_pFrameLaplace->GetMaxValueAndPos( x_rec_start, y_rec_start, x_rec_end, y_rec_end, max_x, max_y );
	}

	
	int ret_max_x,ret_max_y;
	CccdReport* ret=NULL;
	if(bGen){
		ret = AddGeneratedEvent( x0_in , y0_in , x1 , y1,  mag,OtherImage.m_MaxX, 
										 OtherImage.m_MaxY, frame_index, ret_max_x, ret_max_y,
										 eFlash, max_lap );
	}

	if(maxX)
		(*maxX) = ret_max_x;
	if(maxY)
		(*maxY) = ret_max_y;

	
	MYTRACE3(gCCDTrace,"Image Inserted at position : x=" << x0_in << ",y=" << y0_in);

	return ret;
}



/*BOOL_T CCDMatrix::GetImage( long x0, long y0,long len_x, long len_y,
                            Table2D<ELEM_TYPE>& Image )
{
	long x0_in = x0;
	long y0_in = y0; 

	Image.m_X_On_Big = x0_in;
	Image.m_Y_On_Big = y0_in;

	if(x0<0)
		x0_in = 0;
	if(y0<0)
		y0_in = 0; 
	
	long x1 = x0_in + len_x - 1;
	long y1 = y0_in + len_y - 1;

	if(x1<0 || y1<0)
		return FALSE;
	if(x1>=m_SizeX)
		x1 = m_SizeX-1;
	if(y1>=m_SizeY)
		y1 = m_SizeY-1;
	len_x = (x1-x0+1);
	len_y = (y1-y0+1);
	
	Image.Alloc( len_x, len_y );

	ELEM_TYPE* left = GetPos(x0_in,y0_in);
	ELEM_TYPE* right = Image.GetPos(0,0);

	for(int i=0;i<len_y;i++){
		memcpy(right,left,sizeof(ELEM_TYPE)*len_x);
		left = left + m_SizeX;
		right = right + len_x;
	}		
	return TRUE;		
}*/

CCDMatrix& CCDMatrix::operator=(const CCDMatrix& right)
{
	// Table2D<ELEM_TYPE>::Assign(right);
	OperatorEq( right );
	if(right.m_pFrameLaplace && m_pFrameLaplace){
		(*m_pFrameLaplace) = (*right.m_pFrameLaplace);
	}
	return (*this);
}

CCDMatrix& CCDMatrix::operator+=(const CCDMatrix& right)
{
	clock_t t1 = clock();
	
	Assert(m_SizeX==right.m_SizeX,"X - sizes must be equal");
	Assert(m_SizeY==right.m_SizeY,"Y - sizes must be equal");

	long size = (m_SizeX*m_SizeY);
	for(int i=0;i<size;i++){
		m_pData[i] += right.m_pData[i];
	}
	
	clock_t t2 = clock();
	mystring msg = get_clock_in_sec_string( t2-t1 );

	MYTRACE1(gCCDTrace,"Addition of images took : " << msg );	
	return (*this);
}

void CCDMatrix::GetEventDirAndList( CccdReport& event, mystring& szEvtDir, mystring& szEvtList )
{
	szEvtDir = "";
	szEvtList = "";

	char szSubDir[128];
	sprintf(szSubDir,"Frame%.5d",event.m_DayFrameIndex);
	
	szEvtDir << gCCDParams.GetOutputDir() << "/Events/" << szSubDir << "/";
	if( CCDPipeline::GetPipelineCount()>1 )
		szEvtDir << "Cam" << event.m_PipelineIndex << "/";

	szEvtList << "list" << event.EvtIdx;	
}

void CCDMatrix::WriteEventToFITSFile( CccdReport& event,
												  int CurrentFrameIndex,
												  int EventNo, CCDPipeline* pPipeline, 
												  mystring& szEventDir,
												  const char* szMainSubDir )
{
	if(gCCDParams.m_bSaveEventDescOnly){
		// not writing to FITS - only description :
		return;
	}

	szEventDir = "";
	mystring szEventFile,szListFile;
	mystring szBaseFileName;


	char szSubDir[128];
	sprintf(szSubDir,"Frame%.5d",event.m_DayFrameIndex);
	
	szEventDir << gCCDParams.GetOutputDir() << "/" << szMainSubDir << "/" 
				  << szSubDir << "/";
	if( pPipeline->GetPipelineCount()>1 )
		szEventDir << "Cam" << pPipeline->GetPipelineIndex() << "/";
	szBaseFileName << CCDEventList::GetEventFileName( event, (int)CurrentFrameIndex, EventNo ) << ".fit";
	szEventFile << szEventDir << szBaseFileName;
		

	szListFile << szEventDir << "list" << EventNo;
	MyOFile list_file( szListFile.c_str(), "a" );

	mystring szAbsolutePath,szCurDir;
	//if(szEventFile[0]!='/')
	//	szAbsolutePath << MyFile::GetCWD( szCurDir) << "/" << szEventFile;
	//else
	//	szAbsolutePath << szEventFile;
	szAbsolutePath << szBaseFileName;
	if( gCCDParams.m_eCompressFITS != eFITSComprNone)
		szAbsolutePath << "c";
	list_file.Printf("%s\n",szAbsolutePath.c_str());
	list_file.Close();


	int low_x=0,low_y=0,up_x=0,up_y=0;
	

	
	if(gCCDParams.m_bSaveEventSize>0){
		low_x = (int)(MAX((event.m_MaxPoint.x-gCCDParams.m_bSaveEventSize),0));
		low_y = (int)(MAX((event.m_MaxPoint.y-gCCDParams.m_bSaveEventSize),0));
		up_x = (int)(MIN((event.m_MaxPoint.x+gCCDParams.m_bSaveEventSize),(m_SizeX-1)));
		up_y = (int)(MIN((event.m_MaxPoint.y+gCCDParams.m_bSaveEventSize),(m_SizeY-1)));

		if(((up_x-low_x+1)%2)!=0){
			up_x--;
			Assert(up_x>low_x,"Cannot resize window to even size in X");
		}
		if(((up_y-low_y+1)%2)!=0){
			up_y--;
			Assert(up_y>low_y,"Cannot resize window to even size in Y");
		}	
	}

	CKeyTab eventHDUInfo;

	for(int i=0;i<m_KeyTab.GetCount();i++){
		CEnvVar& key = m_KeyTab[i];
		if( gUseMyFITSIO || !IsStandardKey( key.szName.c_str() ) ){
			mystring szValue = key.szValue.TrimApostro();
			eventHDUInfo.Add( key.szName.c_str(), szValue.c_str() );
		}
	}

	// add save area here ( skiped in loop ):
	char savearea[50];
   sprintf(savearea,"%d %d %d %d",low_x,low_y,up_x,up_y);
   eventHDUInfo.Set( SAVEAREA, savearea );

	mystring szEventFrame;
	szEventFrame << event.m_FrameIndex;
	eventHDUInfo.Set(EVENT_FRAME,szEventFrame.c_str());

	mystring szCurrFrame;
	szCurrFrame << CurrentFrameIndex;
	eventHDUInfo.Set(EVENT_CURRENT_FRAME,szCurrFrame.c_str());

	mystring szOrginalX;
	szOrginalX << (long)event.m_MaxPoint.x;
	eventHDUInfo.Set(EVENT_ORGINAL_X,szOrginalX.c_str());

	mystring szOrginalY;
	szOrginalY << (long)event.m_MaxPoint.y;
	eventHDUInfo.Set(EVENT_ORGINAL_Y,szOrginalY.c_str());

	mystring szXOnSmall,szYOnSmall;
	szXOnSmall << ((long)event.m_MaxPoint.x-low_x);
	szYOnSmall << ((long)event.m_MaxPoint.y-low_y);

	eventHDUInfo.Set(EVENT_X_ON_SMALL,szXOnSmall.c_str());
	eventHDUInfo.Set(EVENT_Y_ON_SMALL,szYOnSmall.c_str());


	_TRACE_PRINTF_3( "saving event file : %s ..." , szEventFile.c_str() );fflush(0);
	if(gCCDParams.m_bSaveEventSize>0){
		// writing part of Image :
		eventHDUInfo.Set( TAKEN_AT_X, low_x );
		eventHDUInfo.Set( TAKEN_AT_Y, low_y );
		eventHDUInfo.Set( ORGSIZEX , (int)GetXSize() );
		eventHDUInfo.Set( ORGSIZEY , (int)GetYSize() );
		WriteToFITSFile( szEventFile.c_str() , low_x, low_y, up_x, up_y, 
							  &eventHDUInfo );
	}else{
		// Writing full frame with events :
		eventHDUInfo.Set( TAKEN_AT_X, (int)0 );
      eventHDUInfo.Set( TAKEN_AT_Y, (int)0 );
		WriteToFITSFile( szEventFile.c_str() , FALSE, &eventHDUInfo );
	}
		
	// WriteToFITSFile( szEventFile.c_str() , FALSE, &eventHDUInfo );		
	event.m_LastDumpedFrame = CurrentFrameIndex;

	_TRACE_PRINTF_3("OK\n");
}

void CCDMatrix::AddPartKeys( CSafeKeyTab& eventHDUInfo, 
									  int low_x, int low_y, int up_x, int up_y,
									  int evt_x, int evt_y, 
									  int big_size_x, int big_size_y,
									  int evt_frame_index,
									  int curr_frame_index )
{
	char savearea[50];
   sprintf(savearea,"%d %d %d %d",low_x,low_y,up_x,up_y);
   eventHDUInfo.Set( SAVEAREA, savearea );

	mystring szEventFrame;
	szEventFrame << evt_frame_index;
	eventHDUInfo.Set(EVENT_FRAME,szEventFrame.c_str());

	mystring szCurrFrame;
	szCurrFrame << curr_frame_index;
	eventHDUInfo.Set(EVENT_CURRENT_FRAME,szCurrFrame.c_str());

	mystring szOrginalX;
	szOrginalX << evt_x;
	eventHDUInfo.Set(EVENT_ORGINAL_X,szOrginalX.c_str());

	mystring szOrginalY;
	szOrginalY << evt_y;
	eventHDUInfo.Set(EVENT_ORGINAL_Y,szOrginalY.c_str());

	/*mystring szXOnSmall,szYOnSmall;
	szXOnSmall << ((long)event.m_MaxPoint.x-low_x);
	szYOnSmall << ((long)event.m_MaxPoint.y-low_y);

	eventHDUInfo.Set(EVENT_X_ON_SMALL,szXOnSmall.c_str());
	eventHDUInfo.Set(EVENT_Y_ON_SMALL,szYOnSmall.c_str());*/

	
	eventHDUInfo.Set( EVTX0, low_x );
	eventHDUInfo.Set( EVTY0, low_y );
	eventHDUInfo.Set( ORGSIZEX , big_size_x );
	eventHDUInfo.Set( ORGSIZEY , big_size_y );	
}

void CCDMatrix::SaveEventFrameSum( CccdReport& event, const char* dir_path,
                                   int startFrame, int CurrentFrameIndex,
                                   int EventNo, CCDPipeline* pPipeline )
{
	// in fact this can be merged into single version 
	// the one with gCCDParams.m_ConfirmEventsOnNextNFrames>0 
	// because in case it is =0 both versions are equal
	// this split just for SAFTY reasons not to make it not working

	if( gCCDParams.m_ConfirmEventsOnNextNFrames>0 ){
		_TRACE_PRINTF_3("saving averages - CONFIRMATION ON NEXT MODE\n");fflush(0);
		if( CurrentFrameIndex==(event.m_FrameIndex+gCCDParams.m_ConfirmEventsOnNextNFrames) ){
			// now save average of frames before event :
			CMyStrTable listBefore;
			mystring szOutFile;
			mystring szList,szFileName;

			// save 1 
			GetFrameList( dir_path, event, EventNo, listBefore, szOutFile,
						  startFrame,
						  (CurrentFrameIndex-1-gCCDParams.m_ConfirmEventsOnNextNFrames) );	
			if(!CCDUtil::SaveSumFrame( listBefore, szOutFile.c_str() )){
				gCCDErrLog.DumpToFile1( NULL, "ERR_SUM_PARTS", "Could not save sum of event pictures" );
			}	
			mystring szExt;
			szFileName = getfname( szOutFile.c_str(), szExt );
			szList << dir_path << "/" << "avlist" << EventNo;
			szFileName << "." << szExt;
			if( gCCDParams.m_eCompressFITS != eFITSComprNone)
	      	szFileName << "c";
			szFileName.DumpToFile( szList.c_str() );

			int first = (CurrentFrameIndex-5);
			if(first<=0)
				first = 1;
			if(first<startFrame)
				first=startFrame;

			// save 2 
			GetFrameList( dir_path, event, EventNo, listBefore, szOutFile,
							  first, (CurrentFrameIndex-1-gCCDParams.m_ConfirmEventsOnNextNFrames) );	
			if( !CCDUtil::SaveSumFrame( listBefore, szOutFile.c_str() ) ){
				gCCDErrLog.DumpToFile1( NULL, "ERR_SUM_PARTS", "Could not save sum of event pictures" );
			}	
		
			szFileName = getfname( szOutFile.c_str(), szExt );
			szList = "";
			szList << dir_path << "/" << "avlist" << EventNo;
			szFileName << "." << szExt;
			if( gCCDParams.m_eCompressFITS != eFITSComprNone)
	      	szFileName << "c";
			szFileName.DumpToFile( szList.c_str() );
		}else{
			int after = ( CurrentFrameIndex - ( event.m_FrameIndex + gCCDParams.m_ConfirmEventsOnNextNFrames ) );
			if( after>0 && (after % 5)==0 ){
				// save average of after frames :
				CMyStrTable listAfter;
		      mystring szOutFile;

				GetFrameList( dir_path, event, EventNo, listAfter, szOutFile,
						  (event.m_FrameIndex+1+gCCDParams.m_ConfirmEventsOnNextNFrames), 
						  CurrentFrameIndex );	
				if( !CCDUtil::SaveSumFrame( listAfter, szOutFile.c_str() ) ){
					gCCDErrLog.DumpToFile1( NULL, "ERR_SUM_PARTS", "Could not save sum of event pictures" );
				}
				mystring szList,szFileName,szExt;
				szFileName = getfname( szOutFile.c_str(), szExt );
				szFileName << "." << szExt;
				if( gCCDParams.m_eCompressFITS != eFITSComprNone)
			      szFileName << "c";

				szList << dir_path << "/" << "avlist" << EventNo;
				szFileName.DumpToFile( szList.c_str() );
			}
		}
	}else{
		_TRACE_PRINTF_3("saving averages - NORMAL\n");fflush(0);
		if( CurrentFrameIndex==event.m_FrameIndex ){
			// now save average of frames before event :
			CMyStrTable listBefore;
			mystring szOutFile;
			mystring szList,szFileName;

			// save 1 
			GetFrameList( dir_path, event, EventNo, listBefore, szOutFile,
						  startFrame,
						  (CurrentFrameIndex-1) );	
			if(!CCDUtil::SaveSumFrame( listBefore, szOutFile.c_str() )){
				gCCDErrLog.DumpToFile1( NULL, "ERR_SUM_PARTS", "Could not save sum of event pictures" );
			}	
			mystring szExt;
			szFileName = getfname( szOutFile.c_str(), szExt );
			szList << dir_path << "/" << "avlist" << EventNo;
			szFileName << "." << szExt;
			if( gCCDParams.m_eCompressFITS != eFITSComprNone)
	      	szFileName << "c";
			szFileName.DumpToFile( szList.c_str() );

			int first = (CurrentFrameIndex-5);
			if(first<=0)
				first = 1;
			if(first<startFrame)
				first=startFrame;

			// save 2 
			GetFrameList( dir_path, event, EventNo, listBefore, szOutFile,
							  first, (CurrentFrameIndex-1) );	
			if( !CCDUtil::SaveSumFrame( listBefore, szOutFile.c_str() ) ){
				gCCDErrLog.DumpToFile1( NULL, "ERR_SUM_PARTS", "Could not save sum of event pictures" );
			}	
		
			szFileName = getfname( szOutFile.c_str(), szExt );
			szList = "";
			szList << dir_path << "/" << "avlist" << EventNo;
			szFileName << "." << szExt;
			if( gCCDParams.m_eCompressFITS != eFITSComprNone)
	      	szFileName << "c";
			szFileName.DumpToFile( szList.c_str() );
		}else{
			int after = ( CurrentFrameIndex - event.m_FrameIndex );		
			if( after>0 && (after % 5)==0 ){
				// save average of after frames :
				CMyStrTable listAfter;
		      mystring szOutFile;

				GetFrameList( dir_path, event, EventNo, listAfter, szOutFile,
						  (event.m_FrameIndex+1), 
						  CurrentFrameIndex );	
				if( !CCDUtil::SaveSumFrame( listAfter, szOutFile.c_str() ) ){
					gCCDErrLog.DumpToFile1( NULL, "ERR_SUM_PARTS", "Could not save sum of event pictures" );
				}
				mystring szList,szFileName,szExt;
				szFileName = getfname( szOutFile.c_str(), szExt );
				szFileName << "." << szExt;
				if( gCCDParams.m_eCompressFITS != eFITSComprNone)
			      szFileName << "c";

				szList << dir_path << "/" << "avlist" << EventNo;
				szFileName.DumpToFile( szList.c_str() );
			}
		}
	}
}

void CCDMatrix::GetFrameList( const char* dir_path, CccdReport& event,
										int EventNo,
										CMyStrTable& out_list, mystring& szOutFile,
										int first, int last )
{
	out_list.clear();
	if( first<=0 )
		first = 1;

	// make sure event frame is not added :
	if( first<event.m_FrameIndex ){
		if( last >= event.m_FrameIndex )
	      last = (event.m_FrameIndex-1);
	}
	if( last>event.m_FrameIndex ){
		if( first <= event.m_FrameIndex )
			first = (event.m_FrameIndex+1);
	}

	for(int i=first;i<=last;i++){
		mystring szFileName = dir_path;
		szFileName << "/" << CCDEventList::GetEventFileName( event, i, EventNo ) << ".fit";
		if( gCCDParams.m_eCompressFITS != eFITSComprNone )
			szFileName << "c";
		out_list.push_back( szFileName );					
	}
	szOutFile = dir_path;
	szOutFile << "/" << "avg_evt" << EventNo << "_fr" << first 
				 << "to" << last << ".fit";

}


BOOL_T CCDMatrix::WriteToFITSFile( const char* fname,
                			         int low_x, int low_y, int up_x, int up_y,
                        		   CKeyTab* pHduTab, BOOL_T bNoCompress  )
{
	mystring szError,szFName = fname;
	szFName.env2str();	

	CFITSFile<ELEM_TYPE> out;

	eFITSCompressionType compr=gCCDParams.m_eCompressFITS;
   if( bNoCompress ){
      compr = eFITSComprNone;
   }

	if( gCCDParams.m_bCheckFrame && gCCDParams.m_bRepairFrame && !gCCDParams.GetMC() ){
		printf("CCDMatrix::WriteToFITSFile - repairing of bad frames before save is required !\n");fflush(0);
		int bad_line;
		if( !CheckFrame( bad_line ) ){
			printf("Error detected in frame %s, bad line at y=%d, reparing frame ...\n",fname,bad_line);fflush(0);
			RepairFrame( bad_line );
		}else{
			printf("Frame %s OK - no pixel shift detected\n",fname);fflush(0);
		}
	}
	
	if(!strstr(szFName.c_str(),".fit"))
		szFName << ".fit";

	for(int i=0;i<m_KeyTab.GetCount();i++){
		//if( strcmp( m_KeyTab[i].szName.c_str(), fitsHeaderKeyDefTab[eUT_START].szKeyName)==0 ||
		//	 strcmp( m_KeyTab[i].szName.c_str(), fitsHeaderKeyDefTab[eDATE_OBS].szKeyName)==0 )
		if( gUseMyFITSIO || !IsStandardKey( m_KeyTab[i].szName.c_str() ) ){
			out.AddHdu( m_KeyTab[i].szName.c_str(), m_KeyTab[i].szValue.c_str() );
		}
	}

	if(pHduTab){
		CKeyTab::iterator i;
		for(i=pHduTab->begin();i!=pHduTab->end();i++){
	      if( gUseMyFITSIO || !IsStandardKey( i->szName.c_str() ) ){
				out.SetHdu(i->szName.c_str(),i->szValue.c_str());
			}
		}
	}

	if( low_x>0 && low_y>0 && up_x>0 && up_y>0 ){
		char savearea[50];
   	sprintf(savearea,"%d %d %d %d",low_x,low_y,up_x,up_y);
	   out.SetHdu( SAVEAREA, savearea );
	}

	if (!out.WriteToFITSFile( m_pData,  low_x, low_y, up_x, up_y ,
									  m_SizeX, m_SizeY, 
                             szError, szFName.c_str() , 
									  compr ) ){
		printf("Error while writing to FITS file %s : %s\n",szFName.c_str(),szError.c_str());
		return FALSE;
	}else{
		if(!szError.empty()){
			printf("Warning while writing to FITS file %s : %s\n",szFName.c_str(),szError.c_str());
		}
	}
	m_szFileName = fname;


	return TRUE;	
}

BOOL_T  CCDMatrix::WriteToFITSFile( const char* fname, BOOL_T bDumpEvents,
                                  CKeyTab* pHduTab, BOOL_T bNeverCompress )
{
	mystring szError,szFName = fname;
	szFName.env2str();	

	CFITSFile<ELEM_TYPE> out;
	
	if(!strstr(szFName.c_str(),".fit"))
		szFName << ".fit";

	if( gCCDParams.m_bCheckFrame && gCCDParams.m_bRepairFrame && !gCCDParams.GetMC() ){
		printf("CCDMatrix::WriteToFITSFile - repairing of bad frames before save is required !\n");fflush(0);
		int bad_line;
		if( !CheckFrame( bad_line ) ){
			printf("Error detected in frame %s, bad line at y=%d, reparing frame ...\n",fname,bad_line);fflush(0);
			RepairFrame( bad_line );
		}else{
			printf("Frame %s OK - no pixel shift detected\n",fname);fflush(0);
		}
	}

	for(int i=0;i<m_KeyTab.GetCount();i++){
		//if( strcmp( m_KeyTab[i].szName.c_str(), fitsHeaderKeyDefTab[eUT_START].szKeyName)==0 ||
		//	 strcmp( m_KeyTab[i].szName.c_str(), fitsHeaderKeyDefTab[eDATE_OBS].szKeyName)==0 )
		// printf("Header = %s %s\n",m_KeyTab[i].szName.c_str(), m_KeyTab[i].szValue.c_str() );fflush(0);
		if( gUseMyFITSIO || !IsStandardKey( m_KeyTab[i].szName.c_str() ) ){
			out.AddHdu( m_KeyTab[i].szName.c_str(), m_KeyTab[i].szValue.c_str() );
		}
	}
	if(pHduTab){
		CKeyTab::iterator i;
		for(i=pHduTab->begin();i!=pHduTab->end();i++){
			if( gUseMyFITSIO || !IsStandardKey( i->szName.c_str() ) ){
				out.SetHdu(i->szName.c_str(),i->szValue.c_str());
			}
		}
	}

	mystring szBaseName=szFName.c_str();
	getbasename_new( szFName, szBaseName );
	out.SetHdu( fitsHeaderKeyDefTab[eFILENAME].szKeyName, szBaseName.c_str() );

	eFITSCompressionType compr=gCCDParams.m_eCompressFITS;
	if( bNeverCompress ){
		compr = eFITSComprNone;
	}
	if (!out.WriteToFITSFile( m_pData,  m_SizeX, m_SizeY, 
                             szError, szFName.c_str(), 
									  compr ) ){
		printf("Error while writing to FITS file %s , error is : %s\n",szFName.c_str(),szError.c_str());
		return FALSE;
	}else{
		if(bDumpEvents){
			DumpEventReport( fname );
		}
		if(!szError.empty()){
			printf("Warning while writing to FITS file %s , error is %s",szFName.c_str(),szError.c_str());
		}
	}
	m_szFileName = fname;

	return TRUE;
}

void CCDMatrix::UpdateGenEvent( CccdReport& genEvent, CccdReport& foundEvent )
{
	genEvent.m_GenPoint = genEvent.m_Point;
	genEvent.m_Point = foundEvent.m_Point;
	genEvent.m_bIdentified = foundEvent.m_bIdentified;
}

void CCDMatrix::UpdateFoundEvent( CccdReport& foundEvent, CccdReport& genEvent )
{
	foundEvent.m_GenPoint = genEvent.m_Point;	
	foundEvent.m_LowLeft = genEvent.m_LowLeft;
	foundEvent.m_TopRight = genEvent.m_TopRight;
	foundEvent.m_Magnitude = genEvent.m_Magnitude;
	foundEvent.m_bGenerated = genEvent.m_bGenerated;
}


void CCDMatrix::MoveGeneratedToBegin( CCDEventList& ccdEvents )
{
	if(ccdEvents.size()>1){
		CCDEventList tmpList;
		for(int i=0;i<ccdEvents.size();i++){
			if(ccdEvents[i].m_bGenerated){
				tmpList.Add( ccdEvents[i] );
			}
		}
		for(int i=0;i<ccdEvents.size();i++){
			if(!ccdEvents[i].m_bGenerated){
				tmpList.Add( ccdEvents[i] );
			}
		}
		ccdEvents = tmpList;
	}
}

BOOL_T CCDMatrix::IsGenEventIdentified( const CccdReport& genEvent,
                                        const CccdReport& foundEvent )
{
	const CPoint& maxPointFound = foundEvent.m_MaxPoint;
	const CPoint& maxPointGen = genEvent.m_MaxPoint;

	_TRACE_PRINTF_6("GEN : (%d,%d), FOUND : (%d,%d)\n",(int)maxPointGen.x,(int)maxPointGen.y,(int)maxPointFound.x,(int)maxPointFound.y);
	if( fabs(maxPointFound.x-maxPointGen.x)<=gCCDParams.m_bGenEventRedial &&
  	    fabs(maxPointFound.y-maxPointGen.y)<=gCCDParams.m_bGenEventRedial ){
		return TRUE;
	}
	return FALSE;
}

BOOL_T CCDMatrix::IsGenEventIdentified( int gen_x, int gen_y, int found_x, int found_y )
{
	if( fabs(found_x-gen_x)<=gCCDParams.m_bGenEventRedial &&
  	    fabs(found_y-gen_y)<=gCCDParams.m_bGenEventRedial ){
		return TRUE;
	}
	return FALSE;
}



int CCDMatrix::RejectIfMoreThen( CCDPipeline* pPipeline ){
	return RejectIfMoreThen( m_FoundEventList, pPipeline  );
}

int CCDMatrix::RejectIfMoreThen( CCDEventList& foundEventList, 
											CCDPipeline* pPipeline,
											BOOL_T bCheckConfirmedOnly/*=FALSE*/,
											BOOL_T bAfterCoicCheck/*=FALSE*/ )
{	
	int allCount=foundEventList.size();

	if(gCCDParams.m_MaxNumberOfEventsOnFrame>0){
		if(foundEventList.size() > gCCDParams.m_MaxNumberOfEventsOnFrame){
			printf("\n\nNumber of events on single frame : %d exceeds limit %d, all rejected !!!\n",foundEventList.size(),gCCDParams.m_MaxNumberOfEventsOnFrame);
			int ret = foundEventList.size();
			foundEventList.clear();
			return ret;
		}
	}

	if(gCCDParams.m_bSkipIfMoreThen<=0)
		return 0;	

	for(register int i=0;i<foundEventList.size();i++){
		foundEventList[i].m_bRemoveFromList = FALSE;
	}

	int RejecetedCount=0;
	for(register int i=0;i<foundEventList.size();i++){
		CccdReport& evt_i = foundEventList[i];

		// now checking if event is flash - SUPER NEW are independent and are accepted anyway :
		if( evt_i.m_PixelAnalResults.eventType == eBrighten){
			// skiping SUPER-NEW - means no MORE_THEN rejection
			continue;
		}

		if(bCheckConfirmedOnly){
			if(evt_i.m_PixelAnalResults.m_bRejectedByNextFrames || 
				evt_i.m_PixelAnalResults.m_bRejectedDueToTrack ||
				evt_i.m_PixelAnalResults.m_bRejectedByCoic ){
				continue;
			}
		}

		// if(!foundEventList[i].m_bRemoveFromList){

			int nNearEventsCount=1;
			for(register int j=i+1;j<foundEventList.size();j++){		
				CccdReport& evt_j = foundEventList[j];

				if( evt_j.m_PixelAnalResults.eventType==eBrighten )
					continue;

				if(bCheckConfirmedOnly){
					if(evt_j.m_PixelAnalResults.m_bRejectedByNextFrames || 
						evt_j.m_PixelAnalResults.m_bRejectedDueToTrack ||
						evt_j.m_PixelAnalResults.m_bRejectedByCoic ){
						continue;
					}
				}
				

				register int x_dist = (int)fabs(evt_j.m_Point.x-evt_i.m_Point.x);
				register int y_dist = (int)fabs(evt_j.m_Point.y-evt_i.m_Point.y);

 				if( x_dist<=gCCDParams.m_bSkipIfMoreThenRedial && y_dist<=gCCDParams.m_bSkipIfMoreThenRedial && 
					 (x_dist > gCCDParams.m_bSkipIfMoreThenMinDist || y_dist > gCCDParams.m_bSkipIfMoreThenMinDist ) ){
					nNearEventsCount++;
				}

				if((evt_j.m_Point.y-evt_i.m_Point.y)>gCCDParams.m_bSkipIfMoreThenRedial){
					break;
				}
			}

			if( CCD_Analyser::m_pAnalFoundEventFunc ){
				if( bAfterCoicCheck ){
					(*CCD_Analyser::m_pAnalFoundEventFunc)( &nNearEventsCount, eHistoIfMoreAfterCoic, pPipeline->GetPipelineIndex(), NULL );					
				}else{
					(*CCD_Analyser::m_pAnalFoundEventFunc)( &nNearEventsCount, eHistoIfMore, pPipeline->GetPipelineIndex(), NULL );
				}
			}
	
			if(nNearEventsCount>gCCDParams.m_bSkipIfMoreThen){
				// flag events to be skiped :


				_TRACE_PRINTF_5("REJECTED_ at (%d,%d): !!!!!!!!!!!!!!!!!!!!!\n",(int)evt_i.m_MaxPoint.x,(int)evt_i.m_MaxPoint.y);
				//exit(0);

				if( gCCDParams.m_bRejectIfMoreVerb ){
					printf("IF_MORE, events near by to (%d,%d) rejected\n",
							(int)evt_i.m_Point.x,(int)evt_i.m_Point.y);fflush(0);
				}

				for(register int j=i;j<foundEventList.size();j++){		
					CccdReport& evt_j = foundEventList[j];

					if(bCheckConfirmedOnly){
						if(evt_j.m_PixelAnalResults.m_bRejectedByNextFrames || 
							evt_j.m_PixelAnalResults.m_bRejectedDueToTrack ||
							evt_j.m_PixelAnalResults.m_bRejectedByCoic ){
							continue;
						}
					}

					if( evt_j.m_PixelAnalResults.eventType==eBrighten )
						continue;

					register int x_dist = (int)fabs(evt_j.m_Point.x-evt_i.m_Point.x);
					register int y_dist = (int)fabs(evt_j.m_Point.y-evt_i.m_Point.y);

 					if( x_dist<=gCCDParams.m_bSkipIfMoreThenRedial && y_dist<=gCCDParams.m_bSkipIfMoreThenRedial ){
						evt_j.m_bRemoveFromList = TRUE;
						RejecetedCount++;
					}
				}				
			}else{
				// now check if skip some overlaping events :
				int overlapCount=0;
				// int overlapRedial = MAX_func(gCCDParams.m_OverlapRedial,gCCDParams.m_bSkipIfMoreThenRedial);
				int overlapRedial = gCCDParams.m_OverlapRedial;
				for(register int j=i+1;j<foundEventList.size();j++){		
					CccdReport& evt_j = foundEventList[j];										

					if( evt_j.m_PixelAnalResults.eventType != eBrighten ){
	 					if( (fabs(evt_j.m_Point.x-evt_i.m_Point.x)<=overlapRedial) &&
 							 (fabs(evt_j.m_Point.y-evt_i.m_Point.y)<=overlapRedial) ){
							overlapCount++;
							evt_j.m_bRemoveFromList = TRUE;	
							RejecetedCount++;
						}

						if((evt_j.m_Point.y-evt_i.m_Point.y)>overlapRedial){
							break;
						}
					}
				}
			}
		// }
	}

	vector<CccdReport*> rejectedList;

	int left=0;
	for(register int i=0;i<foundEventList.size();i++){
		if(!foundEventList[i].m_bRemoveFromList){
			left++;
		}else{
			rejectedList.push_back( (CccdReport*)(&(foundEventList[i])) );
		}
	}

	/*if( pPipeline && rejectedList.size()>2 && 
		 gCCDParams.m_bCheckRejectIfMoreTracks && rejectedList.size()<50 ){
		if( pPipeline->GetAnalPtr() ){
			if( pPipeline->GetAnalPtr()->FitLine( rejectedList, pPipeline ) ){
				// line was fitted to events - possibly plane :
				// now add to list of "plane tracks"
			}
		}
	}*/

	
	int Rejected=0;	
	if( left > 0 ){
		// now go throw table and pick only survivied events :
		CCDEventList survList( left );
		for(register int i=0;i<foundEventList.size();i++){
			if( !foundEventList[i].m_bRemoveFromList ){
				survList.Add( foundEventList[i] );
			}
		}			
		Rejected = (allCount-left);
		foundEventList.clear();
		foundEventList = survList;
	}else{
		foundEventList.clear();
		Rejected=allCount;
	}

	return Rejected;
}

LONG_T CCDMatrix::CompileEventReport( CCDEventList& ccdEvents, LONG_T idx, 
												  BOOL_T bAdd, BOOL_T bEraseFoundFromGen,
												  BOOL_T bProfile )
{
	return CompileEventReport( ccdEvents, m_FoundEventList, m_GeneratedEventList,
										this, idx, bAdd, bEraseFoundFromGen, bProfile );
}

LONG_T CCDMatrix::CompileEventReportPtr( CCDEventList* ccdEvents, LONG_T idx, 
												  BOOL_T bAdd, BOOL_T bEraseFoundFromGen,
												  BOOL_T bProfile )
{
	return CompileEventReportPtr( ccdEvents, m_FoundEventList, m_GeneratedEventList,
										this, idx, bAdd, bEraseFoundFromGen, bProfile );
}


LONG_T CCDMatrix::CompileEventReport( CCDEventList& ccdEvents,
												  CCDEventList& foundEvents,
												  CCDEventList& genEvents,
												  CCDMatrix* pMatrix,
												  LONG_T idx, BOOL_T bAdd, 
												  BOOL_T bEraseFoundFromGen,
												  BOOL_T bProfile )
{

	return CompileEventReportPtr( &ccdEvents, foundEvents, genEvents,
											pMatrix, idx, bAdd, bEraseFoundFromGen, bProfile );


/*	CCDEventList::iterator pGenEvt;
	PROFILER_START

	if(!bAdd)
		ccdEvents.clear();
	
	pGenEvt=genEvents.begin();
	while(pGenEvt!=genEvents.end() ){
		CCDEventList::iterator pFoundEvt;				
		BOOL_T bFound=FALSE;
		for(pFoundEvt=foundEvents.begin();pFoundEvt!=foundEvents.end();pFoundEvt++){
			CPoint& pEvt = pFoundEvt->m_Point;
			if(pFoundEvt->IsIdentified()){
				if(IsGenEventIdentified( *pGenEvt, *pFoundEvt )){
					_TRACE_PRINTF_6("\n\n !!! FOUND EVENT at (%d,%d)\n",(long)pEvt.x,(long)pEvt.y);
					if(pMatrix){
						long max_x,max_y;			
						long maxValue = pMatrix->GetMaxValueAndPos( (long)pGenEvt->m_LowLeft.x,(long)pGenEvt->m_LowLeft.y,
																	  (long)pGenEvt->m_TopRight.x, (long)pGenEvt->m_TopRight.y,
																	  max_x,max_y);
						_TRACE_PRINTF_6("max value=%d, at (%d,%d)\n",maxValue,max_x,max_y);fflush(0);
						pMatrix->Dump( (long)pGenEvt->m_LowLeft.x, (long)pGenEvt->m_LowLeft.y, 
   		               (long)pGenEvt->m_TopRight.x ,(long)pGenEvt->m_TopRight.y );
					}

					UpdateFoundEvent( *pFoundEvt, *pGenEvt );
					bFound = TRUE;
					break;
				}
			}
		}
		if(bFound){
			if(bEraseFoundFromGen){
				genEvents.erase( pGenEvt );
			}else{
				UpdateGenEvent( *pGenEvt, *pFoundEvt );
				pGenEvt++;
			}
		}else{
			pGenEvt++;
		}			
	}

	CCDEventList::iterator pEvt;

	// add only not identified here (identified are included in found events list :
   for(pEvt=genEvents.begin();pEvt!=genEvents.end();pEvt++){
		pEvt->m_FrameIndex = idx;
		if(!pEvt->IsIdentified())
			ccdEvents.push_back( *pEvt );		
	}

	
	for(pEvt=foundEvents.begin();pEvt!=foundEvents.end();pEvt++){		
		pEvt->m_FrameIndex = idx;
		if(pEvt->IsIdentified())
			ccdEvents.push_back( *pEvt );
	}
	MoveGeneratedToBegin( foundEvents );
	if(bProfile){
		PROFILER_END("CompileEventReport took:")
	}

	return ccdEvents.size();*/
}


LONG_T CCDMatrix::CompileEventReportPtr( CCDEventList* ccdEvents,
												  CCDEventList& foundEvents,
												  CCDEventList& genEvents,
												  CCDMatrix* pMatrix,
												  LONG_T idx, BOOL_T bAdd, 
												  BOOL_T bEraseFoundFromGen,
												  BOOL_T bProfile )
{
	CCDEventList::iterator pGenEvt;
	PROFILER_START

	if(ccdEvents){
		if(!bAdd)
			ccdEvents->clear();
	}
	
	pGenEvt=genEvents.begin();
	while(pGenEvt!=genEvents.end() ){
		CCDEventList::iterator pFoundEvt;				
		BOOL_T bFound=FALSE;
		for(pFoundEvt=foundEvents.begin();pFoundEvt!=foundEvents.end();pFoundEvt++){
			CPoint& pEvt = pFoundEvt->m_Point;
			if(IsGenEventIdentified( *pGenEvt, *pFoundEvt )){
				_TRACE_PRINTF_6("\n\n !!! FOUND EVENT at (%d,%d)\n",(long)pEvt.x,(long)pEvt.y);
				if(pMatrix){
					long max_x,max_y;			
					long maxValue = pMatrix->GetMaxValueAndPos( (long)pGenEvt->m_LowLeft.x,(long)pGenEvt->m_LowLeft.y,
																  (long)pGenEvt->m_TopRight.x, (long)pGenEvt->m_TopRight.y,
																  max_x,max_y);
					_TRACE_PRINTF_6("max value=%d, at (%d,%d)\n",maxValue,max_x,max_y);fflush(0);
					pMatrix->Dump( (long)pGenEvt->m_LowLeft.x, (long)pGenEvt->m_LowLeft.y, 
   	               (long)pGenEvt->m_TopRight.x ,(long)pGenEvt->m_TopRight.y );
				}

				UpdateFoundEvent( *pFoundEvt, *pGenEvt );
				bFound = TRUE;
				break;
			}
		}
		if(bFound){
			if(bEraseFoundFromGen){
				genEvents.erase( pGenEvt );
			}else{
				UpdateGenEvent( *pGenEvt, *pFoundEvt );
				pGenEvt++;
			}
		}else{
			pGenEvt++;
		}			
	}

	int ret =0;

	CCDEventList::iterator pEvt;
	// add only not identified here (identified are included in found events list :
   for(pEvt=genEvents.begin();pEvt!=genEvents.end();pEvt++){
		pEvt->m_FrameIndex = idx;

		if(!pEvt->m_bIdentified){
			if(ccdEvents){
				ccdEvents->push_back( *pEvt );		
			}
			ret++;
		}
	}

	
	for(pEvt=foundEvents.begin();pEvt!=foundEvents.end();pEvt++){
		pEvt->m_FrameIndex = idx;

		if(ccdEvents)
			ccdEvents->push_back( *pEvt );
		ret++;
	}

	MoveGeneratedToBegin( foundEvents );
	if(bProfile){
		PROFILER_END("CompileEventReport took:")
	}

	if(ccdEvents)
		return ccdEvents->size();
	else
		return ret;
}


void CCDMatrix::DumpEventReport( const char* fname ){
	CCDEventList ccdReport;
	DumpEventReport( fname , ccdReport );
}

void CCDMatrix::DumpEventReport( const char* fname , CCDEventList& ccdReport )
{
	// first compile the report - what was identified vs what was generated
   // then write such report to a file 
	ccdReport.clear();
	CompileEventReport( ccdReport );
	CCDEventList::DumpEventReport( fname , ccdReport );
}

BOOL_T CCDMatrix::ParseSaveArea( const char* szVal, 
											int& x0, int& y0, int& x1, int& y1 )
{
	x0 = 0;
	y0 = 0;
	x1 = 0;
	y1 = 0;
	if( szVal && szVal[0] ){
   	CMyStrTable items;
      MyParser szPars = szVal;
      szPars.GetItems(items);
      if( items.size()>=2 ){
      	if( atol( items[0].c_str() )>0 )
         	x0 = atol( items[0].c_str() );
	      if( atol( items[1].c_str() )>0 )
   	      y0 = atol( items[1].c_str() );
      }
      if( items.size()>=4 ){
      	if( atol( items[2].c_str() )>0 )
         	x1 = atol( items[2].c_str() );
	      if( atol( items[3].c_str() )>0 )
   	      y1 = atol( items[3].c_str() );
      }
		
		return TRUE;
	}
	return FALSE;
}

BOOL_T  CCDMatrix::ReadFITSHeader(  const char* fname )
{
	mystring szError,szFName = fname;
	szFName.env2str();

	CFITSFile<ELEM_TYPE> input;
/*	if(!input.ReadFITSHeader( *this, szError, szFName.c_str() )){
		printf("could not read header file from file  : %s\n",szFName.c_str());
		return FALSE;
	}
	m_szFileName = getfname( szFName );		*/

	return TRUE;
	
}

BOOL_T  CCDMatrix::ReadFITSFile(  const char* fname, BOOL_T bMustBeAllocated )
{
	if(bMustBeAllocated){
		Assert(m_pData!=NULL,"Memory for image not allocted - cannot ReadFITSFile");
	}

	mystring szError,szFName = fname;
	szFName.env2str();

	int sizeX = m_SizeX;
	int sizeY = m_SizeY;
	
	clock_t t1=clock();
	CFITSFile<ELEM_TYPE> input;
	if(!input.ReadFITSFile( *this, szError, szFName.c_str() )){
		printf("could not read frame : %s\n",szFName.c_str());
		return FALSE;
	}else{
		if(!szError.empty()){
			MYTRACE1(gCCDTrace,"Warning while reading FITS file : " << szFName
									  << ", error is : " << szError);
		}
		if( sizeX>0 && sizeY>0 && gCCDParams.m_bCheckSizeWhenRead ){
			Assert(sizeX==m_SizeX && sizeY==m_SizeY,"Wrong initial values for m_SizeX=%d and m_SizeY=%d, fit file %s has : (%d,%d)",sizeX,sizeY,szFName.c_str(),m_SizeX,m_SizeY);
		}

		const char* keyX = fitsHeaderKeyDefTab[eTAKEN_AT_X].szKeyName;
		const char* keyY = fitsHeaderKeyDefTab[eTAKEN_AT_Y].szKeyName;

		CEnvVar* pKeyX = GetKeyTab().Find( keyX );
		CEnvVar* pKeyY = GetKeyTab().Find( keyY );
		if( pKeyX && pKeyY ){
			m_X_On_Big = atol( mystring::get_number( pKeyX->szValue ).c_str() );
			m_Y_On_Big = atol( mystring::get_number( pKeyY->szValue ).c_str() );		
		}else{
			CEnvVar* pSaveArea = GetKeyTab().Find( SAVEAREA );
         if( pSaveArea ){
				int  x0,y0,x1,y1;
				ParseSaveArea( pSaveArea->szValue.c_str(), x0, y0, x1, y1 );
				m_X_On_Big = x0;
				m_Y_On_Big = y0;
			}
		}

	}
	m_szFileName = getfname( szFName );		
	clock_t t2=clock();
	mystring msg=get_clock_in_sec_string( t2-t1 );
	_TRACE_PRINTF_3("Reading of file : %s, took %s\n",szFName.c_str(),msg.c_str());

	return TRUE;
}

BOOL_T CCDMatrix::GetPictureParams( const char* fname, int& sizeX, int& sizeY, int& x, int& y )
{
	mystring szError,szFName = fname;
	szFName.env2str();
	
	CFITSFile<ELEM_TYPE> input;
	if(!input.ReadFITSFile( *this, szError, szFName.c_str() )){
		MYTRACE1(gCCDTrace,"Error while reading FITS file : " << szFName
                         << ", error is : " << szError);
		return FALSE;
	}else{
		if(!szError.empty()){
			MYTRACE1(gCCDTrace,"Warning while reading FITS file : " << szFName
									  << ", error is : " << szError);
		}
	}
	m_szFileName = getfname( szFName );		


	const char* keyX = fitsHeaderKeyDefTab[eNAXIS1].szKeyName;
	const char* keyY = fitsHeaderKeyDefTab[eNAXIS2].szKeyName;

	CEnvVar* pKeyX = GetKeyTab().Find( keyX );
	CEnvVar* pKeyY = GetKeyTab().Find( keyY );

	if(pKeyX && pKeyY){
		sizeX = atol( pKeyX->szValue.c_str() );
		sizeY = atol( pKeyY->szValue.c_str() );
	}

	keyX = fitsHeaderKeyDefTab[eTAKEN_AT_X].szKeyName;
	keyY = fitsHeaderKeyDefTab[eTAKEN_AT_Y].szKeyName;
	pKeyX = GetKeyTab().Find( keyX );
	pKeyY = GetKeyTab().Find( keyY );

	if(pKeyX && pKeyY){
		x = atol( pKeyX->szValue.c_str() );
		y = atol( pKeyY->szValue.c_str() );
	}

	return TRUE;	
}


void CCDMatrix::ClearState()
{
	 m_bIntersting = FALSE;
	 m_FoundEventList.clear();	 	
 	 m_GeneratedEventList.clear();	 
 	 // m_IdentifiedEvents.clear();
}

void CCDMatrix::AdjustFrame()
{
	
}


LONGLONG_T CCDMatrix::CalcSumAround( LONG_T x, LONG_T y , LONG_T r0 )
{
	LONGLONG_T sum=0;
	LONG_T start_x = MAX(0,x-r0);
	LONG_T start_y = MAX(0,y-r0);
	LONG_T end_x = MIN(m_SizeX-1,x+r0);
	LONG_T end_y = MIN(m_SizeY-1,y+r0);

	for(LONG_T x0=start_x;x0<end_x;x0++){
		for(LONG_T y0=start_y;y0<end_y;y0++){
			LONG_T pos = x0 + y0*m_SizeX;
			double r = sqrt((x0-x)*(x0-x)+(y0-y)*(y0-y));
			if(r <= r0){
				sum += m_pData[pos];
			}
		}
	}
	return sum;
}


Table2D<BIG_ELEM_TYPE>* CCDMatrix::InitLaplaceFrame()
{
	if(!m_pFrameLaplace){
		m_pFrameLaplace = new Table2D<BIG_ELEM_TYPE>( m_SizeX, m_SizeY );
	}
	return m_pFrameLaplace;
}


void CCDMatrix::Laplace( long x0, long y0, long x1, long y1 )
{
	if(!m_pFrameLaplace){
		InitLaplaceFrame();
	}
	
	if(m_pFrameLaplace){
		BIG_ELEM_TYPE** p_laplace_fast = m_pFrameLaplace->get_data_buffer_fast();
		Table2D<ELEM_TYPE>::Laplace( p_laplace_fast, gCCDParams.m_eLaplaceType, x0, y0, x1, y1 );
	}		
}

void CCDMatrix::LaplaceSafe( long x0, long y0, long x1, long y1 )
{
	if(m_pFrameLaplace){
		BIG_ELEM_TYPE** p_laplace_fast = m_pFrameLaplace->get_data_buffer_fast();
		Table2D<ELEM_TYPE>::LaplaceSafe( p_laplace_fast,
						 gCCDParams.m_LaplacePlusList, gCCDParams.m_LaplacePlusCount,						 
						 gCCDParams.m_LaplaceMinusList, gCCDParams.m_LaplaceMinusCount,
						 x0, y0, x1, y1 );
	}		
}


void CCDMatrix::Laplace( CCDMatrix& out, long x0, long y0, long x1, long y1 )
{
	Table2D<ELEM_TYPE>::Laplace( out, gCCDParams.m_eLaplaceType, 
				gCCDParams.GetEdgeSizeInLaplace(),
				x0, y0, x1, y1 );
}

void CCDMatrix::LaplaceEdgeOnly( long ignore_edge )
{
	if(m_pFrameLaplace){
		ELEM_TYPE* p_data = m_pData;
		BIG_ELEM_TYPE* p_laplace = m_pFrameLaplace->get_data_buffer();
		BIG_ELEM_TYPE** p_laplace_fast = m_pFrameLaplace->get_data_buffer_fast();
		long upY = m_SizeY-ignore_edge;
      long upX = m_SizeX-ignore_edge;

		for(register long y=0;y<ignore_edge;y++){
			for(register long x=0;x<m_SizeX;x++){
				register int laplace = CalcLaplaceSumEstimateOnEdge( x, y, 
														m_SizeX, m_SizeY,
														m_pFastData,
														gCCDParams.m_LaplacePlusList, gCCDParams.m_LaplacePlusCount,
														gCCDParams.m_LaplaceMinusList, gCCDParams.m_LaplaceMinusCount );
				p_laplace_fast[y][x] = laplace;
			}
		}
		for(register long y=upY;y<m_SizeY;y++){
			for(register long x=0;x<m_SizeX;x++){
				register int laplace = CalcLaplaceSumEstimateOnEdge( x, y, 
														m_SizeX, m_SizeY,
														m_pFastData,
														gCCDParams.m_LaplacePlusList, gCCDParams.m_LaplacePlusCount,
														gCCDParams.m_LaplaceMinusList, gCCDParams.m_LaplaceMinusCount );
				p_laplace_fast[y][x] = laplace;
			}
		}
		for(register long y=0;y<m_SizeY;y++){
			for(register long x=0;x<ignore_edge;x++){
				register int laplace = CalcLaplaceSumEstimateOnEdge( x, y, 
														m_SizeX, m_SizeY,
														m_pFastData,
														gCCDParams.m_LaplacePlusList, gCCDParams.m_LaplacePlusCount,
														gCCDParams.m_LaplaceMinusList, gCCDParams.m_LaplaceMinusCount );
				p_laplace_fast[y][x] = laplace;
			}
		}
		for(register long y=0;y<m_SizeY;y++){
			for(register long x=upX;x<m_SizeX;x++){
				register int laplace = CalcLaplaceSumEstimateOnEdge( x, y, 
														m_SizeX, m_SizeY,
														m_pFastData,
														gCCDParams.m_LaplacePlusList, gCCDParams.m_LaplacePlusCount,
														gCCDParams.m_LaplaceMinusList, gCCDParams.m_LaplaceMinusCount );
				p_laplace_fast[y][x] = laplace;
			}
		}
	}

}


void CCDMatrix::LaplaceEdgeOnly( long ignore_edge, 
											CLongPoint* plus_list, LONG_T plus_count,
											CLongPoint* minus_list, LONG_T minus_count )
{
	if(m_pFrameLaplace){
		ELEM_TYPE* p_data = m_pData;
		BIG_ELEM_TYPE* p_laplace = m_pFrameLaplace->get_data_buffer();
		BIG_ELEM_TYPE** p_laplace_fast = m_pFrameLaplace->get_data_buffer_fast();
		long upY = m_SizeY-ignore_edge;
      long upX = m_SizeX-ignore_edge;

		for(register long y=0;y<ignore_edge;y++){
			for(register long x=0;x<m_SizeX;x++){
				register int laplace = CalcLaplaceSumEstimateOnEdge( x, y, 
														m_SizeX, m_SizeY,
														m_pFastData,
														plus_list, plus_count,
														minus_list, minus_count );
				p_laplace_fast[y][x] = laplace;
			}
		}
		for(register long y=upY;y<m_SizeY;y++){
			for(register long x=0;x<m_SizeX;x++){
				register int laplace = CalcLaplaceSumEstimateOnEdge( x, y, 
														m_SizeX, m_SizeY,
														m_pFastData,
														plus_list, plus_count,
                                          minus_list, minus_count );
				p_laplace_fast[y][x] = laplace;
			}
		}
		for(register long y=0;y<m_SizeY;y++){
			for(register long x=0;x<ignore_edge;x++){
				register int laplace = CalcLaplaceSumEstimateOnEdge( x, y, 
														m_SizeX, m_SizeY,
														m_pFastData,
														plus_list, plus_count,
                                          minus_list, minus_count );
				p_laplace_fast[y][x] = laplace;
			}
		}
		for(register long y=0;y<m_SizeY;y++){
			for(register long x=upX;x<m_SizeX;x++){
				register int laplace = CalcLaplaceSumEstimateOnEdge( x, y, 
														m_SizeX, m_SizeY,
														m_pFastData,
														plus_list, plus_count,
                                          minus_list, minus_count );
				p_laplace_fast[y][x] = laplace;
			}
		}
	}

}


void CCDMatrix::Laplace( BOOL_T bEstimateOnEdge/*=TRUE*/ )
{
	if(!m_pFrameLaplace){
      InitLaplaceFrame();
   }

	if(m_pFrameLaplace){
		m_pFrameLaplace->Alloc( m_SizeX, m_SizeY );

		ELEM_TYPE* p_data = m_pData;
		BIG_ELEM_TYPE* p_laplace = m_pFrameLaplace->get_data_buffer();
		BIG_ELEM_TYPE** p_laplace_fast = m_pFrameLaplace->get_data_buffer_fast();

		register int ignore_edge = gCCDParams.GetEdgeSizeInLaplace();
		register int upY = m_SizeY-ignore_edge;
		register int upX = m_SizeX-ignore_edge;

		BOOL_T bDone=FALSE;
		if( gCCDParams.m_eLaplaceType==eFivePlusFourMin){
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
																	  gCCDParams.m_eLaplaceType, plus_sum,minus_sum );
																				  
				}
			}
		}


		if(bEstimateOnEdge){
			LaplaceEdgeOnly( ignore_edge );
		}
	}
}

void CCDMatrix::Laplace( eLaplaceType_T laplace_type, BOOL_T bEstimateOnEdge )
{
	// calculation of plus and minus lists :
	CLongPoint plus_list[64],minus_list[64];
	LONG_T plus_count, minus_count;
	GetPlusMinusList( laplace_type, plus_list, plus_count,
              			minus_list, minus_count );

	if(!m_pFrameLaplace){
      InitLaplaceFrame();
   }

	if(m_pFrameLaplace){
		m_pFrameLaplace->Alloc( m_SizeX, m_SizeY );

		ELEM_TYPE* p_data = m_pData;
		BIG_ELEM_TYPE* p_laplace = m_pFrameLaplace->get_data_buffer();
		BIG_ELEM_TYPE** p_laplace_fast = m_pFrameLaplace->get_data_buffer_fast();

		register int ignore_edge = gCCDParams.GetEdgeSizeInLaplace( laplace_type );
		register int upY = m_SizeY-ignore_edge;
		register int upX = m_SizeX-ignore_edge;

		BOOL_T bDone=FALSE;
		if( laplace_type == eFivePlusFourMin){
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
																	   laplace_type, plus_sum,minus_sum );
																				  
				}
			}
		}


		if(bEstimateOnEdge){
			LaplaceEdgeOnly( ignore_edge, 
								  plus_list, plus_count,
		                    minus_list, minus_count  );
		}
	}
}



long CCDMatrix::CalcGaussFilter( Table2D<double>& gaussMatrix, CCDMatrix& out, 
											long threshold,
							     			long x0, long y0, long x1, long y1 )
{
	out.Alloc( m_SizeX, m_SizeY );

	// Assert( (gaussMatrix.GetXSize() % 2 ) == 1," 
	Assert( gaussMatrix.GetXSize()==gaussMatrix.GetYSize(),"Gauss matrix must be a square");
	long radius = (gaussMatrix.GetXSize()-1)/2;
	long center = radius;
	long gmsize = gaussMatrix.GetXSize();


	long upY = m_SizeY-radius;
	long upX = m_SizeX-radius;
	long bottomX = radius;
	long bottomY = radius;

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
	
  long pixcount=0; //stores 'how many pixel were calculated' information

     /* For each pixel that is above 'threshold' multiply surrounding
        pixels by corresponding 'gauss_matrix' values and store
	  	  the sum in 'workcopy' table (float) */

	ELEM_TYPE** p_data = get_data_buffer_fast();
	ELEM_TYPE** p_data_out = out.get_data_buffer_fast();	
	double** p_gauss_matrix = gaussMatrix.get_data_buffer_fast();

	for (register long y=bottomY; y < upY ; y++){
		for (register long x=bottomX; x < upX ; x++){
			p_data_out[y][x] = 0;
			if ( p_data[y][x] > threshold){
			   for (register int jy=0; jy<gmsize; jy++){
				   for (register int ix=0; ix<gmsize; ix++){
						p_data_out[y][x] += (int)((p_data[y-center+jy][x-center+ix])*
												  p_gauss_matrix[jy][ix]);
					}
		      }
			   pixcount++;
			}else{
				_TRACE_PRINTF_3("!!!!!!!! (%d,%d) - below treshold => 0\n",x,y);fflush(0);
			}
		}
	}

	return pixcount;	
}

long CCDMatrix::CalcGaussFilter( Table2D<double>& gaussMatrix, Table2D<double>& out, 
											long threshold,
							     			long x0, long y0, long x1, long y1 )
{
	out.Alloc( m_SizeX, m_SizeY );

	// Assert( (gaussMatrix.GetXSize() % 2 ) == 1," 
	Assert( gaussMatrix.GetXSize()==gaussMatrix.GetYSize(),"Gauss matrix must be a square");
	long radius = (gaussMatrix.GetXSize()-1)/2;
	long center = radius;
	long gmsize = gaussMatrix.GetXSize();


	long upY = m_SizeY-radius;
	long upX = m_SizeX-radius;
	long bottomX = radius;
	long bottomY = radius;

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
	
  long pixcount=0; //stores 'how many pixel were calculated' information

     /* For each pixel that is above 'threshold' multiply surrounding
        pixels by corresponding 'gauss_matrix' values and store
	  	  the sum in 'workcopy' table (float) */

	ELEM_TYPE** p_data = get_data_buffer_fast();
	double** p_data_out = out.get_data_buffer_fast();	
	double** p_gauss_matrix = gaussMatrix.get_data_buffer_fast();

	for (register long y=bottomY; y < upY ; y++){
		for (register long x=bottomX; x < upX ; x++){
			p_data_out[y][x] = 0;
			if ( p_data[y][x] > threshold){
			   for (register int jy=0; jy<gmsize; jy++){
				   for (register int ix=0; ix<gmsize; ix++){
						p_data_out[y][x] += (p_data[y-center+jy][x-center+ix])*
												  p_gauss_matrix[jy][ix];
					}
		      }
			   pixcount++;
			}else{
				_TRACE_PRINTF_3("!!!!!!!! (%d,%d) - below treshold => 0\n",x,y);fflush(0);
			}
		}
	}

	return pixcount;	
}



LONGLONG_T CCDMatrix::CalcSum( LONG_T* pixels, LONG_T cnt )
{
	return Table2D<ELEM_TYPE>::CalcSum(pixels, cnt);
}

LONGLONG_T CCDMatrix::CalcSum( CLongList& pixels )
{
	CLongList::iterator i;
	LONGLONG_T sum=0;
	for(i=pixels.begin();i!=pixels.end();i++){
		sum += m_pData[*i];
	}
	return sum;
}



void CCDMatrix::Transform( eConfShape_T shapeType, double Redial, CCDMatrix& out )
{
	CLongPoint shapePoints[100];
	int nPoints = CCD_Analyser::CalcShapePoints( shapePoints, shapeType, Redial );
	register long i;

	_TRACE_PRINTF_3("\n\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	_TRACE_PRINTF_3("SHAPE POINTS :\n");
	for(i=0;i<nPoints;i++){
		_TRACE_PRINTF_3("(%d,%d) ",(long)shapePoints[i].x,(long)shapePoints[i].y);
	}
	_TRACE_PRINTF_3("\n\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");

	Assert(nPoints<100,"Maximum number of points for transformation exceeded");
	ELEM_TYPE** p_out = out.get_data_buffer_fast();
	ELEM_TYPE** p_in  = get_data_buffer_fast();

	for(register long x=0;x<m_SizeX;x++){
		for(register long y=0;y<m_SizeY;y++){


			double out_sum=0;
			for(i=0;i<nPoints;i++){
					long x_pos = x+shapePoints[i].x;
					long y_pos = y+shapePoints[i].y;
					if(x_pos>=0 && x_pos<m_SizeX && y_pos>=0 && y_pos<m_SizeY){
						out_sum += p_in[y_pos][x_pos];
					}
			}
			p_out[y][x] = (long)out_sum;
		}
	}	
	
}



const char* CCDMatrix::getKeyValue( const char* keyname )
{
	for(register int i=0;i<m_KeyTab.GetCount();i++){
		if(strcmp(m_KeyTab[i].szName.c_str(),keyname)==0)
			return m_KeyTab[i].szValue.c_str();
	}
	return "";
}

void CCDMatrix::setKeyValue( const char* keyname, int value )
{
	for(register int i=0;i<m_KeyTab.GetCount();i++){
		if(strcmp(m_KeyTab[i].szName.c_str(),keyname)==0){
			char szValue[256];
			sprintf(szValue,"%d",value);
			m_KeyTab[i].szValue = szValue;
			return;
		}
	}
	m_KeyTab.Add( keyname, value );
}


void CCDMatrix::setKeyValue( const char* keyname, double value )
{
	for(register int i=0;i<m_KeyTab.GetCount();i++){
		if(strcmp(m_KeyTab[i].szName.c_str(),keyname)==0){
			m_KeyTab[i].szValue = value;
			return;
		}
	}
	m_KeyTab.Add( keyname, value );
}

void CCDMatrix::setKeyValue( const char* keyname, const char* value )
{
	for(register int i=0;i<m_KeyTab.GetCount();i++){
		if(strcmp(m_KeyTab[i].szName.c_str(),keyname)==0){
			m_KeyTab[i].szValue = value;
			return;
		}
	}
	m_KeyTab.Add( keyname, value );
}


const char* CCDMatrix::getKeyValue( const char* keyname, int start_pos )
{
	if( start_pos<0 || start_pos>=m_KeyTab.GetCount() ){
		start_pos = 0;
	}
	for(;start_pos<m_KeyTab.GetCount();start_pos++){
		if(strcmp(m_KeyTab[start_pos].szName.c_str(),keyname)==0)
			return m_KeyTab[start_pos].szValue.c_str();
	}
	return getKeyValue( keyname );
	
}


BOOL_T CCDMatrix::IsDark()
{
	const char* szObject = getKeyValue( OBJECT );
	if( szObject && szObject[0] ){
		if( strstr( szObject, "DARK" ) ){
			return TRUE;
		}
	}		
	const char* szExpTime = getKeyValue( EXPTIME );
	if( szExpTime && szExpTime[0] ){
		return ( atol(szExpTime)<0 );
	}
	return FALSE;
}

time_t CCDMatrix::getMinTime()
{
	time_t ret=0;
	const char* szMINTIME = getKeyValue( "MINTIME" );
	if( szMINTIME && szMINTIME[0] ){
		ret = atol( szMINTIME );
		return ret;
	}

	ret = (time_t)getObsTime( TRUE );
	return ret;
}

time_t CCDMatrix::getMaxTime()
{
	time_t ret=0;
	const char* szMAXTIME = getKeyValue( "MAXTIME" );
	if( szMAXTIME && szMAXTIME[0] ){
		ret = atol( szMAXTIME );
		return ret;
	}

	ret = (time_t)getObsTime( TRUE );
	return ret;

}

BOOL_T CCDMatrix::getObsTime( struct tm* _tm )
{
	memset( _tm,'\0',sizeof(struct tm));
	
	const char* szDATE_TMP     = getKeyValue( (gCCDParams.m_DateObs).c_str() );
	const char* szTIME_TMP     = NULL;

	if( strlen( (gCCDParams.m_TimeObs).c_str() ) ){
		szTIME_TMP = getKeyValue( (gCCDParams.m_TimeObs).c_str() );
	}

	mystring szDateString, szTimeString;

	if( szDATE_TMP && szDATE_TMP[0] ){
		szDateString = szDATE_TMP;
		szDateString.SkipApostrophs();
	}	
	if( szTIME_TMP && szTIME_TMP[0] ){
		szTimeString = szTIME_TMP;
		szTimeString.SkipApostrophs();
	}

	const char* szDATE = szDateString.c_str();
	const char* szTIME = szTimeString.c_str();

	const char* szDT_FMT = gCCDParams.m_DateTimeObsFormat.c_str();
	

	mystring szDateTime,szFMT;

	if(szTIME && strlen(szTIME)){
		szDateTime << szDATE << "T" << szTIME;
		szFMT << DATE_FMT  << "T" << TIME_FMT;
	}else{
		szDateTime << szDATE;
		szFMT << DATE_OBS_FMT;
	}

	if(strlen(szDateTime.c_str())<=1){
		return FALSE;
	}


	_TRACE_PRINTF_6("getObsTime : %s\n",szDateTime.c_str());

	int msec = 0;
	int match = sscanf(szDateTime.c_str(),szDT_FMT,&_tm->tm_year,&_tm->tm_mon,&_tm->tm_mday,&_tm->tm_hour,&_tm->tm_min,&_tm->tm_sec,&msec);

	if(match<6){
		return FALSE;
	}
	if(match==6)
		msec = 0;

	return TRUE;
}

double CCDMatrix::getObsTime( BOOL_T bIgnoreError )
{
	const char* szUtTime = getKeyValue( TIME_UT );
	if( szUtTime && strlen(szUtTime) && atol( szUtTime )>0){
		double ret = (double)atol( szUtTime );
		return ret;
	}

// when testing dirver - uncomment this line :
// there is no FITS header passing throw shared memory yet (cannot read
// DATE/TIME from header )
//return 0;


	const char* szDATE_TMP     = getKeyValue( (gCCDParams.m_DateObs).c_str() );
	const char* szTIME_TMP     = NULL;

	if( strlen( (gCCDParams.m_TimeObs).c_str() ) ){
		szTIME_TMP = getKeyValue( (gCCDParams.m_TimeObs).c_str() );
	}

	mystring szDateString, szTimeString;

	if( szDATE_TMP && szDATE_TMP[0] ){
		szDateString = szDATE_TMP;
		szDateString.SkipApostrophs();
	}	
	if( szTIME_TMP && szTIME_TMP[0] ){
		szTimeString = szTIME_TMP;
		szTimeString.SkipApostrophs();
	}

	const char* szDATE = szDateString.c_str();
	const char* szTIME = szTimeString.c_str();


	mystring szDT_FMT = gCCDParams.m_DateTimeObsFormat.c_str();
	

	mystring szDateTime,szFMT;

	if(szTIME && strlen(szTIME)){
		szDateTime << szDATE << "T" << szTIME;
		szFMT << DATE_FMT  << "T" << TIME_FMT;
	}else{
		szDateTime << szDATE;
		szFMT << DATE_OBS_FMT;
	}

	if(strlen(szDateTime.c_str())<=1){
		if(!bIgnoreError)
			Assert(FALSE,"Date and time not defined");
		else
			return -1.00;
	}


	_TRACE_PRINTF_6("getObsTime : %s\n",szDateTime.c_str());

	// time_t ret = CMyDate::getTime( "2003-02-02T18:04:33.024", DATE_OBS_FMT );

	// double ret = (double)CMyDate::getTime( szDT, szDT_FMT );

	struct tm _tm;	
	memset( &_tm,'\0',sizeof(_tm));

	// strptime( szDateTime.c_str() , szFMT.c_str() , &_tm );
	int msec = 0;
	int match = sscanf(szDateTime.c_str(),szDT_FMT.c_str(),&_tm.tm_year,&_tm.tm_mon,&_tm.tm_mday,&_tm.tm_hour,&_tm.tm_min,&_tm.tm_sec,&msec);

	if(match<6){
		// retry with apostrohps :
		mystring szDT_FMT_tmp;
		szDT_FMT_tmp << "'" << szDT_FMT << "'";
		match = sscanf(szDateTime.c_str(),szDT_FMT_tmp.c_str(),&_tm.tm_year,&_tm.tm_mon,&_tm.tm_mday,&_tm.tm_hour,&_tm.tm_min,&_tm.tm_sec,&msec);

		if( match<6 ){
			if(!bIgnoreError){
				Assert(FALSE,"Could not get time from %s using format : %s",szDateTime.c_str(),szDT_FMT.c_str());
			}else{
				return -1.00;
			}
		}
	}
	if(match==6)
		msec = 0;


	// in header we have data presented in the normal way :
	// YYYY-MM-DD HH24:MI:SS 
	// and months are counted from 1 ( not as required in tm structer from 0 ) :
	convert_date_from_fits_header_format( &_tm );

	//_tm.tm_year = 2003;
	//_tm.tm_mday = 10;
	//_tm.tm_mon  = 10;
	//_tm.tm_hour = 10;
	//_tm.tm_min  = 10;
	//_tm.tm_sec  = 10;
	// strcpy(_tm.tm_zone,"GMT");
	// double ret = (double)mktime( &_tm );		
	double ret = (double)timegm( &_tm );
	

	if( ret<0 ){
		if(!bIgnoreError){
			printf("ERROR : %s\n",strerror(errno));
			Assert(FALSE,"Could not get DATE-TIME from FITS header");
		}else{
			return -1.00;
		}
	}

	ret += (double(msec)/1000.00);

	return ret;		
}


double CCDMatrix::getObsTime( CSafeKeyTab& keys, BOOL_T bIgnoreError )
{
	const char* szUtTime = keys.getKeyVal( TIME_UT );
	if( szUtTime && strlen(szUtTime) && atol( szUtTime )>0){
		double ret = (double)atol( szUtTime );
		return ret;
	}

// when testing dirver - uncomment this line :
// there is no FITS header passing throw shared memory yet (cannot read
// DATE/TIME from header )
//return 0;


	const char* szDATE_TMP     = keys.getKeyVal( (gCCDParams.m_DateObs).c_str() );
	const char* szTIME_TMP     = NULL;

	if( strlen( (gCCDParams.m_TimeObs).c_str() ) ){
		szTIME_TMP = keys.getKeyVal( (gCCDParams.m_TimeObs).c_str() );
	}

	mystring szDateString, szTimeString;

	if( szDATE_TMP && szDATE_TMP[0] ){
		szDateString = szDATE_TMP;
		szDateString.SkipApostrophs();
	}	
	if( szTIME_TMP && szTIME_TMP[0] ){
		szTimeString = szTIME_TMP;
		szTimeString.SkipApostrophs();
	}

	const char* szDATE = szDateString.c_str();
	const char* szTIME = szTimeString.c_str();


	mystring szDT_FMT = gCCDParams.m_DateTimeObsFormat.c_str();
	

	mystring szDateTime,szFMT;

	if(szTIME && strlen(szTIME)){
		szDateTime << szDATE << "T" << szTIME;
		szFMT << DATE_FMT  << "T" << TIME_FMT;
	}else{
		szDateTime << szDATE;
		szFMT << DATE_OBS_FMT;
	}

	if(strlen(szDateTime.c_str())<=1){
		if(!bIgnoreError)
			Assert(FALSE,"Date and time not defined");
		else
			return -1.00;
	}


	_TRACE_PRINTF_6("getObsTime : %s\n",szDateTime.c_str());

	// time_t ret = CMyDate::getTime( "2003-02-02T18:04:33.024", DATE_OBS_FMT );

	// double ret = (double)CMyDate::getTime( szDT, szDT_FMT );

	struct tm _tm;	
	memset( &_tm,'\0',sizeof(_tm));

	// strptime( szDateTime.c_str() , szFMT.c_str() , &_tm );
	int msec = 0;
	int match = sscanf(szDateTime.c_str(),szDT_FMT.c_str(),&_tm.tm_year,&_tm.tm_mon,&_tm.tm_mday,&_tm.tm_hour,&_tm.tm_min,&_tm.tm_sec,&msec);

	if(match<6){
		// retry with apostrohps :
		mystring szDT_FMT_tmp;
		szDT_FMT_tmp << "'" << szDT_FMT << "'";
		match = sscanf(szDateTime.c_str(),szDT_FMT_tmp.c_str(),&_tm.tm_year,&_tm.tm_mon,&_tm.tm_mday,&_tm.tm_hour,&_tm.tm_min,&_tm.tm_sec,&msec);

		if( match<6 ){
			if(!bIgnoreError){
				Assert(FALSE,"Could not get time from %s using format : %s",szDateTime.c_str(),szDT_FMT.c_str());
			}else{
				return -1.00;
			}
		}
	}
	if(match==6)
		msec = 0;


	// in header we have data presented in the normal way :
	// YYYY-MM-DD HH24:MI:SS 
	// and months are counted from 1 ( not as required in tm structer from 0 ) :
	convert_date_from_fits_header_format( &_tm );

	//_tm.tm_year = 2003;
	//_tm.tm_mday = 10;
	//_tm.tm_mon  = 10;
	//_tm.tm_hour = 10;
	//_tm.tm_min  = 10;
	//_tm.tm_sec  = 10;
	// strcpy(_tm.tm_zone,"GMT");
	// double ret = (double)mktime( &_tm );		
	double ret = (double)timegm( &_tm );
	

	if( ret<0 ){
		if(!bIgnoreError){
			printf("ERROR : %s\n",strerror(errno));
			Assert(FALSE,"Could not get DATE-TIME from FITS header");
		}else{
			return -1.00;
		}
	}

	ret += (double(msec)/1000.00);

	return ret;		
}




BOOL_T CCDMatrix::CalcAzimutalCoord( int x, int y, double& azim, double& altit,
												 CCDProcState& ccdInfo )
{
	time_t ut_time = (int)getObsTime();
	if(1){
		if(ccdInfo.m_pAstroForms){
			BOOL_T bRet = ccdInfo.CalcAzimutalCoord( x, y , ut_time, azim, altit );
			return bRet;
		}else{
			return FALSE;
		}
	}else{
		printf("Could not determine frame time\n");
		return FALSE;
	}
	return TRUE;				
}

BOOL_T CCDMatrix::CalcAstroCoordinates( time_t ut_time, int x, int y, double& azim, double& altit, double& dec,
                                        double& ra, double& ha, CCDProcState& ccdInfo )
{
	if(ccdInfo.m_pAstroForms){
		BOOL_T bRet = ccdInfo.CalcAstroCoord( x, y , ut_time, azim, altit, dec, ra, ha );
		return bRet;
	}
	return FALSE;
}

BOOL_T CCDMatrix::CalcAstroCoordinates( int x, int y, double& azim, double& altit, double& dec,
                  		                double& ra, double& ha, CCDProcState& ccdInfo )
{
	time_t ut_time = (int)getObsTime();
	if (1){
		BOOL_T bRet = CalcAstroCoordinates( ut_time, x, y , azim, altit, dec, ra, ha, ccdInfo );
		return bRet;
	}else{
		printf("Could not determine frame time\n");
		return FALSE;
	}
	return TRUE;	
}

void CCDMatrix::CutBorder( int left, int up, int right, int bottom,
                			   CCDMatrix& out_image )
{
	int size_x = m_SizeX - (left+right);
	int size_y = m_SizeY - (bottom+up);

	out_image.CCDMatrix_InitConstructor( size_x, size_y );
	ELEM_TYPE** out_data = out_image.get_data_buffer_fast();

	for(int y=bottom;y<(m_SizeY-up);y++){
      for(int x=left;x<(m_SizeX-right);x++){
         out_data[y-bottom][x-left] = m_pFastData[y][x];
      }
   }
	
	
	char savearea[50];
   sprintf(savearea,"%d %d %d %d",left,bottom,(m_SizeX-right),(m_SizeY-up));
	out_image.GetKeyTab().Set( SAVEAREA	, savearea );
}

void CCDMatrix::CatBorder( int border, BOOL_T bPutMean ){
	if(border<=0)
		return;

	double mean,rms;
	GetMeanAndRMS( mean, rms );

	ELEM_TYPE** lap_data = get_data_buffer_fast();
	ImageStat stat;

	if( bPutMean ){
		int total=50;
		int half = (total/2);

		for(register int y=0;y<gCCDParams.m_SizeY;y++){
         for(register int x=0;x<gCCDParams.m_SizeX;x++){
				double mean,rms;
				if(x<border || y<border || x>=(gCCDParams.m_SizeX-border) || y>=(gCCDParams.m_SizeY-border)){		
					if(x<border){
						GetStatistics( stat, FALSE, 
													 border,y-half,border+total,y+half );					
					}else{
						if( y<border ){
							GetStatistics( stat, FALSE,
											x-half,border,x+half,border+total );					
						}else{
							if( x>=(gCCDParams.m_SizeX-border) ){
								GetStatistics( stat, FALSE,
											(gCCDParams.m_SizeX-border)-total,y-half,(gCCDParams.m_SizeX-border),y+half );
							}else{
								if( y>=(gCCDParams.m_SizeY-border) ){
									GetStatistics( stat, FALSE,
										x-half,(gCCDParams.m_SizeY-border)-total,x+half,(gCCDParams.m_SizeY-border) );
								}
							}
						}
					}
					lap_data[y][x] = stat.Average;
				}				
			}
		}
	}else{
		for(register int y=0;y<gCCDParams.m_SizeY;y++){
			for(register int x=0;x<gCCDParams.m_SizeX;x++){
				if(x<border || y<border || x>=(gCCDParams.m_SizeX-border) || y>=(gCCDParams.m_SizeY-border)){
						lap_data[y][x] = 0;
				}
			}
		}		
	}
}

void CCDMatrix::CatBorderFast( int border ){
	if(border<=0)
		return;

	double mean,rms;
	GetMeanAndRMS( mean, rms );

	ELEM_TYPE** lap_data = get_data_buffer_fast();
	ImageStat stat;

	int total=50;
	int half = (total/2);
	int size=10000;
	LONG_T tab[10000];
	int step=10;

	for( register int y=border;y<(m_SizeY-border);y=y+step){
		// left border
		int backgr_value = GetBackgroundValue( border, y-half, border+total, y+half, 
															tab,size);
		for(register int x=0;x<border;x++){
			for(register int i=0;i<step;i++){
				lap_data[y+i][x] = backgr_value;
			}
		}

		// right border
		backgr_value = GetBackgroundValue( m_SizeX-(border+total), y-half, m_SizeX-border,  y+half, tab,size );
		for(register int x=(m_SizeX-border);x<m_SizeX;x++){
			for(register int i=0;i<step;i++){
				lap_data[y+i][x] = backgr_value;
			}
		}		
	}

	for( register int x=border;x<(m_SizeX-border);x=x+step){
		// down border
		int backgr_value = GetBackgroundValue( x-half, border, x+half, border+total ,tab,size );
		for(register int y=0;y<border;y++){
			for(register int i=0;i<step;i++){
				lap_data[y][x+i] = backgr_value;
			}
		}

		// up border
		backgr_value = GetBackgroundValue( x-half, m_SizeY-(border+total), x+half, m_SizeY-border, tab,size  );
		for(register int y=(m_SizeY-border);y<m_SizeY;y++){
			for(register int i=0;i<step;i++){
				lap_data[y][x+step] = backgr_value;
			}
		}		
	}

	// now corners :	
	int ll_val = (lap_data[border+1][0]+lap_data[0][border+1])/2;
	int lu_val = (lap_data[m_SizeY-border-1][0]+lap_data[m_SizeY-1][border+1])/2;
	for( register int x=0;x<border;x++){
		for(register int y=0;y<border;y++){
			lap_data[y][x] = ll_val;
		}
		for(register int y=(m_SizeY-border);y<m_SizeY;y++){
			lap_data[y][x] = lu_val;
		}	
	}

	int rl_val = (lap_data[0][m_SizeX-border-1]+lap_data[border+1][m_SizeX-1])/2;
	int ru_val = (lap_data[m_SizeY-border-1][m_SizeX-1]+lap_data[m_SizeY-1][m_SizeX-border-1])/2;
	for( register int x=(m_SizeX-border);x<m_SizeX;x++){
		for(register int y=0;y<border;y++){
			lap_data[y][x] = rl_val;
		}
		for(register int y=(m_SizeY-border);y<m_SizeY;y++){
			lap_data[y][x] = ru_val;
		}	
	}
}

const int CCD_PIX_HI = 65535;

const int BG_MAP_X_RES = 30; //30
const int BG_MAP_Y_RES = 30; //30

int CCDMatrix::createBackgroundLuminosityMap( int hibad_range,
															 int x_start, int y_start,
															 int x_end, int y_end )
{
	int average_sky=0,sigma_sky=0;

	 vector< vector<int> > background_map;
	 vector<short> background_map_x;        // x-coord. of bg map points
    vector<short> background_map_y;

	background_map.resize(BG_MAP_X_RES);
   background_map_x.resize(BG_MAP_X_RES);
   background_map_y.resize(BG_MAP_Y_RES);

   for (int j=0; j < BG_MAP_X_RES ; j++)   background_map_x.push_back(0);
   for (int j=0; j < BG_MAP_Y_RES ; j++)   background_map_y.push_back(0);

   for (int i=0; i < BG_MAP_X_RES ; i++)
   {
      background_map[i].reserve(BG_MAP_Y_RES);
      for (int j=0; j < BG_MAP_Y_RES ; j++)
			background_map[i].push_back(5);
   }


    const int points_per_area=50; // tu chyba powinno byc 30x30 = 90 ???
    const int ppa_low = 15; //this is used to quickly find typical (lowest ?) value for area
	 const int square_size=(BG_MAP_X_RES*BG_MAP_Y_RES);
    
    int istep, istep_low, x0, y0, x_step, y_step;
    short start_x, start_y;
    int real_ppa, middle, x, y;
    bool clouds = false;
    
    int n, m=0, hibad;
    double median = 0.0, sigma = 0.0, diff;

    vector<int> sky_vals;
    sky_vals.reserve(points_per_area);

	start_x = x_start;	// x_start defines image begining, start_x defines current area begining
	start_y = y_start;	// the same for 'y'
	int end_x = x_end;
	int end_y = y_end;

	short x_range = ( end_x - start_x );
	short y_range = ( end_y - start_y );

	// ustalenie krokow: dla wspolrzednych x i y na bg_mapie i dla pixeli wewnatrz jednego obszaru

	short bg_map_x_step = x_range / BG_MAP_X_RES;
	short bg_map_y_step = y_range / BG_MAP_Y_RES;
	
	x_step = (int) bg_map_x_step;
	y_step = (int) bg_map_y_step;
	
    	int temp_step = int( bg_map_x_step*bg_map_y_step / points_per_area );
	if ( temp_step > 1 ) istep = temp_step;
	  else		     istep = 1;

	istep_low = (points_per_area / ppa_low) * istep;	// 5 * istep

	for (int area_col=0; area_col < BG_MAP_X_RES; area_col++)
	   background_map_x[area_col] = int( bg_map_x_step * ( area_col+0.5 ) + 0.5 );
	   
	for (int area_row=0; area_row < BG_MAP_Y_RES; area_row++)
	   background_map_y[area_row] = int( bg_map_y_step * ( area_row+0.5 ) + 0.5 );

	// ponizsze mozna usunac w czasie optymalizacji

	for (int area_col=0; area_col < BG_MAP_X_RES; area_col++)
	   background_map_x[area_col] += x_start;
	   
	for (int area_row=0; area_row < BG_MAP_Y_RES; area_row++)
	   background_map_y[area_row] += y_start;
	   
	   	   
	//lecimy po wszystkich obszarach
	
	for (int area_col=0; area_col < BG_MAP_X_RES; area_col++, start_x += bg_map_x_step)
	{
	   x0 = (int) start_x;
	   start_y = y_start;
	   
	   for (int area_row=0; area_row < BG_MAP_Y_RES; area_row++, start_y += bg_map_y_step)
	    {
		n = 0;
		y0 = (int) start_y;
		hibad = CCD_PIX_HI;

		// wczytanie do talbicy sky_vals wartosci 15 pixli 
		// z badanego kwadratu 30x30 		
		for (long i=0; n < ppa_low ; i += istep_low, n++)
		{
			x = i/y_step + x0;
			y = i%y_step + y0;
			sky_vals.push_back( m_pFastData[y][x] );
		}
		
		sort(sky_vals.begin(), sky_vals.end());
		
		// ustalenie gornego poziomu, powyzej ktorego pixele sa ignorowane
		// to jest 3 wartosc z 15 na 30x30 square + hibad_range( default = 125 )
		hibad = sky_vals[2] + hibad_range;
		
		sky_vals.clear();		

		// szukanie sredniej wartosci tla nieba dla obszaru			
		n = 0;
//		for (long i=0; n < points_per_area ; i += istep, n++)
		for (long i=0; n < square_size ; i += istep, n++)						
	  	{
	    		x = i/y_step + x0;
			y = i%y_step + y0;
			if( y<end_y && x<end_x ){
				if ( m_pFastData[y][x] < hibad )
					sky_vals.push_back( m_pFastData[y][x] );
				if( sky_vals.size()>= points_per_area ){
					break;
				}
			}
	   	}

		/* 'p.p.a.' does not always equal 'sky_vals.size()', so get real value  */
		
		real_ppa = sky_vals.size();
		
	
		if ( real_ppa < points_per_area / 2 )
		{
			printf("Warning! Low pixel count per area %d / %d\n",real_ppa,points_per_area);
			clouds = true;
		}

		/* sort and take middle value */
		
		sort(sky_vals.begin(), sky_vals.end());

		middle = real_ppa / 2;
		
		median += sky_vals[middle];
		
		// szukanie odchylenia standardowego

		for (int i=0; i < real_ppa ; i += 2)
		{
			diff =  sky_vals[i] - sky_vals[middle];
			sigma += diff * diff;
			m++;
		}

		// tworzenie mapy tla (mediana):
		
		background_map[area_col][area_row] = sky_vals[middle] ;
		
		// czyszczenie zbioru pixeli dla aktualnego obszaru

		sky_vals.clear();
	   }
	}

	// srednia jasnosc srednia tla nieba i standardowe odchylenie:
	
	int areas = BG_MAP_X_RES * BG_MAP_Y_RES;
	sigma = sqrt( sigma / m );
	median = median / areas;
	
	average_sky = (int) median;
	sigma_sky = (int) sigma;

	printf("BGSKY (avrg,sig) : ( %.2f , %.2f )\n",median,sigma);
	
	/*if (clouds) return 0;
	else
	return average_sky;*/
	return average_sky;
}

void CCDMatrix::CalcColumnStat( int sigma_col,
										  double& median, double& sigma34_14, 
										  double& mean, double& sigma )
{
	int size_y = GetYSize();
	LONG_T* col_data = new LONG_T[size_y];	
	ELEM_TYPE** data = get_data_buffer_fast();
	double sum=0;
   double sum2=0;
   int count=0;
	
	for(int i=30;i<(size_y-30);i++){
  		sum = sum + data[i][sigma_col];
      sum2 = sum2 + (data[i][sigma_col])*(data[i][sigma_col]);
      col_data[count] = data[i][sigma_col];
      count++;
   }
   my_qsort( col_data , count );
	median = col_data[count/2];
   sigma34_14 = col_data[(int)(0.75*count)] - col_data[(int)(0.25*count)];
   mean=sum/count;
   double mean2=(sum2/count);
   sigma=sqrt(mean2 - mean*mean);

	delete [] col_data;
}

int CCDMatrix::CalcColumnStat( int col, int y_start, int y_end,
								  double& mean, double& sigma )
{
	int size_y = GetYSize();
   ELEM_TYPE** data = get_data_buffer_fast();
   double sum=0;
   double sum2=0;
	int count=0;

	for( int y=y_start;(y<y_end && y<size_y);y++){
		sum = sum + data[y][col];
      sum2 = sum2 + (data[y][col])*(data[y][col]);
      count++;
	}
	
	mean=sum/count;
   double mean2=(sum2/count);
   sigma=sqrt(mean2 - mean*mean);

	return count;
}

BOOL_T CCDMatrix::FindShift( int col , int step_size_adu, int star_tresh,
									  int& bad_line )
{
	bad_line = -1;
	int size_y = GetYSize();
	ELEM_TYPE** data = get_data_buffer_fast();
	double sum=0;
   double sum2=0;
	BOOL_T bShift=FALSE;

	int max_diff=-1;
	int max_diff_y=-1;

	int start_row=3;
	int prev_value = data[start_row-1][col];
	int count=0;

//	for(int y=10;y<(size_y-10);y++){
//	for(int y=8;y<(size_y-10);y++){
	for(int y=start_row;y<(size_y-10);y++){
		sum = sum + data[y][col];
      sum2 = sum2 + (data[y][col])*(data[y][col]);
		count++;

		int new_val = data[y][col];			
		// find step - but not due to star 
		int max_val = MAX(new_val,prev_value);
		if( abs( new_val - prev_value ) > step_size_adu && bad_line<0 ){
			if( max_val<star_tresh){
				printf("Pixel shift detected at rows %d -> %d\n",(y-1),y);
				printf("value[%d][%d] = %d\n",(y-1),col,prev_value);
				printf("value[%d][%d] = %d\n",y,col,new_val);
				bad_line = y;
				bShift = TRUE;
			}
		}
		if( abs( new_val - prev_value ) > max_diff ){
			max_diff = abs( new_val - prev_value );
			max_diff_y = y;
		}
		prev_value = new_val;
	}

	double mean=sum/count;
   double mean2=(sum2/count);
   double sigma=sqrt(mean2 - mean*mean);
	printf("sigma( column=%d ) = %.2f ADU, max_diff(y=%d) = %d, tresh = %d\n",
				col,sigma,max_diff_y,max_diff,step_size_adu);


	return bShift;
}

BOOL_T CCDMatrix::RepairFrame( int start_bad_line )
{
	printf("Repairing pixel shift starting at y=%d\n",start_bad_line);fflush(0);
	// shift by single pixel 
	
		
	int size=m_SizeX*m_SizeY;
	int last_pixel=(size-1);
	int start_bad_index=((start_bad_line+1)*m_SizeX);	

	printf("start_bad_line = %d\n",start_bad_line);
	int i=last_pixel;
	while( i >= start_bad_index ){
		m_pData[i] = m_pData[i-1];
		i--;
	}
	GetKeyTab().Set( "QUALITY" , "REPAIRED_PIXEL_SHIFT" );
	
	return TRUE;
}

BOOL_T CCDMatrix::CheckFrame( int& bad_line, BOOL_T bForce /*=FALSE*/ )
{
	bad_line = -1;

	// in case calculated treshold is > then this value it will be set 
	// to value max_step_tresh
	const double max_step_tresh=2000.00;

	if( !bForce ){
		const char* szQuality = getKeyValue( QUALITY );
		if( szQuality && szQuality[0] && strlen( szQuality )>0  ){
			printf("Frame already checked/fixed, skiped\n");
			return TRUE;
		}
	}

	int size_y = GetYSize();
	printf("Checking frame for pixel shift , checking column=%d for value change greater then tresh=%d\n",gCCDParams.m_ShiftCol,gCCDParams.m_ShiftVal);fflush(0);

	double mean,sigma,median,sigma34_14;
	double mean_sigma=0,mean_median=0,normal_sigma=0;

	CalcColumnStat( 30, median,sigma34_14 , mean,sigma );
	mean_sigma = mean_sigma + sigma34_14;
	mean_median = mean_median + median;
	normal_sigma = normal_sigma + sigma;
	CalcColumnStat( 40, median,sigma34_14 , mean,sigma );
   mean_sigma = mean_sigma + sigma34_14;
	mean_median = mean_median + median;
	normal_sigma = normal_sigma + sigma;
	CalcColumnStat( 50, median,sigma34_14 , mean,sigma );
   mean_sigma = mean_sigma + sigma34_14;
	mean_median = mean_median + median;
	normal_sigma = normal_sigma + sigma;
	CalcColumnStat( 60, median,sigma34_14 , mean,sigma );
   mean_sigma = mean_sigma + sigma34_14;
	mean_median = mean_median + median;
	normal_sigma = normal_sigma + sigma;

	mean_sigma = ( mean_sigma / 4.0 );	
	mean_median = ( mean_median / 4.0 );
	normal_sigma = ( normal_sigma / 4.0 );
			
	
	int adu_step_size = mean_sigma*gCCDParams.m_ShiftVal;
	if( gCCDParams.m_ShiftVal >= 200 )
		adu_step_size = gCCDParams.m_ShiftVal;
	if( adu_step_size > max_step_tresh ){
		printf("adu_step_size = %d > %d , set %d\n",(int)adu_step_size,
				(int)max_step_tresh,(int)max_step_tresh);
		adu_step_size = max_step_tresh;
	}
	int star_tresh = (int)(mean_median + 5*mean_sigma);

	printf("ShiftVal=%d Avg_median(col=30,40,50,60) = %.2f Avg_Sigma(col=30,40,50,60) = %.2f, Step_Tresh[ADU] = %d ,Star_Tresh=%d[ADU], Normal_sigma=%.2f\n",
				gCCDParams.m_ShiftVal,
				mean_median,mean_sigma,(int)adu_step_size,star_tresh,normal_sigma);


	int bad_line15,bad_line16,bad_line23,bad_line2055,bad_line0;

//	BOOL_T bShift2055 = FindShift( 2055, adu_step_size, star_tresh, bad_line2055 );
	BOOL_T bShift2055 = FindShift( 2060, adu_step_size, star_tresh, bad_line2055 );	
	BOOL_T bZeroAt2055=FALSE;
	if( bShift2055 ){
		if( m_pFastData[bad_line2055][2060]==0 && bad_line2055>0 && m_pFastData[bad_line2055-1][2060]>500 ){
			printf("ZERO at column 2060, change %d  -> %d\n",
				m_pFastData[bad_line2055-1][2060],m_pFastData[bad_line2055][2060]);
			bZeroAt2055 = TRUE;			

			// in case 0 detected threshold for step is smaller :
			if( gCCDParams.m_ShiftVal < 200 && gCCDParams.m_ShiftVal>3 ){		
				int step_size_sav = adu_step_size;
				adu_step_size = mean_sigma*3;	
				printf("ZERO at column 2060, step_size changed : %d => %d\n",step_size_sav,(int)adu_step_size);
			}
		}
	}


	BOOL_T bShift0 = FindShift( 0, adu_step_size, 100000, bad_line0 );
	BOOL_T bShift15 = FindShift( gCCDParams.m_ShiftCol, adu_step_size, 100000, bad_line15 );
	BOOL_T bShift16 = FindShift( gCCDParams.m_ShiftCol+1, adu_step_size, 100000, bad_line16 );
	BOOL_T bShift23 = FindShift( 23, adu_step_size, star_tresh, bad_line23 );

	int bad_lines[10];
	int bad_count=0;
	if( bShift0 ){
		bad_lines[bad_count] = bad_line0;
      bad_count++;
   }
	if( bShift15 ){
		bad_lines[bad_count] = bad_line15;
		bad_count++;
	}
	if( bShift16 ){
		bad_lines[bad_count] = bad_line16;
		bad_count++;
	}
	if( bShift23 ){
		bad_lines[bad_count] = bad_line23;
		bad_count++;
	}
	if( bShift2055 ){
		bad_lines[bad_count] = bad_line2055;
		bad_count++;
	}

	BOOL_T bOK=TRUE;


	if( bad_count>=2 ){
		int coic_count=0;
		for(int i=0;i<bad_count;i++){
			coic_count=1;
			for(int j=(i+1);j<bad_count;j++){
				if( abs(bad_lines[i]-bad_lines[j]) <= 1 ){
					coic_count++;
				}
			}
			if( coic_count >= 2 ){
       		bOK = FALSE;
				bad_line = bad_lines[i];
				break;
	      }
		}
	}

	if(!bOK ){
		double mean_low,mean_up,sigma_low,sigma_up;
		int low_count = CalcColumnStat( 15, 10, bad_line-2, mean_low, sigma_low );
		int up_count = CalcColumnStat( 15, bad_line+2, (size_y-10), mean_up, sigma_up );
		printf("below %d (mean,sigma)=(%.2f,%.2f), above (mean,sigma)=(%.2f,%.2f)\n",
					bad_line,mean_low, sigma_low,mean_up, sigma_up );
		
		double sig = MAX(sigma_low, sigma_up );
//		if( fabs( mean_low - mean_up ) < normal_sigma*gCCDParams.m_ShiftVal ){
//		if( fabs( mean_low - mean_up ) < normal_sigma*3 ){
//			bOK = TRUE;
//		}		
	}

	if( !bOK ){
		GetKeyTab().Set( "QUALITY" , "PIXEL_SHIFT" );
		GetKeyTab().Set( "ERRLINE" , bad_line );			
	}else{
		GetKeyTab().Set( "QUALITY" , "GOOD" );
	}

	return bOK;	
}


CMyMutex gGenImageLock;
int CCDMatrix::GenSkyImage( double sigma, double mean,
									 double ra, double dec, double fi, double pixscale,
									 BOOL_T bCheckAlt /*=TRUE*/  )
{
	printf("IMAGE GENERATOR (ra,dec)=(%.2f [h],%.2f [deg]) fi=%.2f[deg] pixscale=%.2f (mean,sigma)=(%.2f,%.2f)\n",ra,dec,fi,pixscale,mean,sigma);
	// params :
	double alt_limit=10.00;
	
	gGenImageLock.Lock();

	double i0=150000000.0;
	time_t ut = get_dttm();
	double azim,alt;
	int radius_default=20;
	GenEmptyImage( sigma, mean );
	BOOL_T bAddStars=TRUE;
	int added=0;

	int SizeX = m_SizeX;
	int SizeY = m_SizeY;
	if( SizeX<=0 ){
		SizeX = gCCDParams.m_SizeX;
	}
	if( SizeY<=0 ){
		SizeY = gCCDParams.m_SizeY;
	}

	double alt_in_deg = 0.00;
	if( bCheckAlt ){
		double longit = gCCDParams.m_GeoLongitude;
   	double lat  = gCCDParams.m_GeoLatitude;
		AstroCCD::calculateHorizontalCoordinatesFromEq( AstroAngle::deg2rad( dec ),
																		AstroAngle::hours2rad( ra ),
																		ut, longit, lat, alt, azim );
		alt_in_deg = AstroAngle::rad2deg( alt );
		printf("Image (azim,alt)=(%.2f,%.2f) [deg]\n",AstroAngle::rad2deg(azim),alt_in_deg);
		if( alt_in_deg < alt_limit ){
			bAddStars = FALSE;
		}
	}

/*	if( bAddStars ){															
		vector<CCatalogStar> starList;
		double ra_rad = AstroAngle::hours2rad( ra );
		double dec_rad = AstroAngle::deg2rad( dec );
		double radius_rad = AstroAngle::deg2rad( gCCDParams.m_fFrameSize );

		gStarCatDefault.getStarList( ra_rad, dec_rad, radius_rad, starList );
		printf("Read %d stars from catalog around position (%.2f,%.2f)\n",starList.size(),ra,dec);fflush(stdout);
		printf("adding stars to image of size = (%d,%d) fi=%.2f pixscale=%.2f...\n",SizeX,SizeY,fi,pixscale);fflush(stdout);

		CCDAsasTransform asas_transform( &gCCDParams, NULL );
		asas_transform.InitSimple( ra, dec, fi, pixscale );	

		CMyProgressBar bar(0,starList.size());
		int j=0;
		for( vector<CCatalogStar>::iterator i=starList.begin();i!=starList.end();i++){
			double mag = i->mag;
//			double intensity = pow( 10, 0.4*(mag-21.00) );
			double intensity = i0*pow( 10, -0.4*mag );
			int x,y;
			asas_transform.ad2xy_raw( i->ra, i->dec, x, y );


			if( x>=8 && x<(SizeX-8) && y>=8 && y<(SizeY-8) ){
//				printf("%.2f %.2f %.2f\n",i->ra,i->dec,i->mag);
//				CMyFit::CreateGaussFunc( x+0.5, y+0.5, 20.00 );
				CMyMathFunc::m_X = x+0.5;
				CMyMathFunc::m_Y = y+0.5;

				int radius = radius_default;

				// first check effective radius : 
				// new version 0m 16.760s 
				// old version 2m 6sec !!!! - zysk jest ewidentny !!!
				// po porownaniu klatek - sa identyczne !
				for (int xx=(x-radius);xx<=x;xx++){
					double g = 0.00;
//					double g = CMyFit::GaussIntegral( xx, y, xx+1, y+1 );
					double pix_value = intensity*g;

					if( pix_value < 1.00 ){
						radius = abs( x - xx );
					}
				}
	
				for (int yy=(y-radius);yy<=(y+radius);yy++){
					for (int xx=(x-radius);xx<=(x+radius);xx++){
						if( xx>=0 && xx<SizeX && yy>=0 && yy<SizeY ){
							double g = 0.00;
//							double g = CMyFit::GaussIntegral( xx, yy, xx+1, yy+1 );
							double pix_value = intensity*g;

//							printf("%.4f/%.2f ",g,pix_value);

							double val = m_pFastData[yy][xx] + pix_value;
							if( val >= USHRT_MAX ){
								val = USHRT_MAX;
							}
							m_pFastData[yy][xx] = val;		
	
						}
					}
//					printf("\n");
				}
				added++;
//				exit(0);
			}


			bar.SetValue(j);
			bar.Update();
		
			j++;
		}
		printf("%d stars added correctly to image\n",added);fflush(stdout);
	}else{
		printf("h=%.2f deg , lower then dome edge, no stars added\n",alt_in_deg);fflush(stdout);
	}*/


	gGenImageLock.UnLock();
	return 1;
}


void CCDMatrix::Data_Correction(double Tczyt, double Teksp)
{
	ELEM_TYPE** obrazek = m_pFastData;

	// correction (for each column), initial value "0":
   std::vector<double> poprawka(m_SizeX,0);

   for(int y(0);y<m_SizeY;++y){
		for(int x(0);x<m_SizeX;++x){
			//correcting data:
			obrazek[y][x]=obrazek[y][x]-poprawka[x];
			//computing next correction value, precision of an ELEM_TYPE:
			poprawka[x]+=static_cast<ELEM_TYPE>(( (static_cast<double>(obrazek[y][x])) /Teksp)*Tczyt);
		}
	}
}
