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
#include <mypixellist.h>
#include "ccd_analyse.h"
#include "ccd_pipeline.h"
#include "ccd_trace.h"
#include "ccd_globals.h"
#include "ccd_datares.h"
#include "ccd_excp.h"
#include "ccd_log.h"
#include <math.h>
#include <mymacros.h>
#include <tab2Ddesc.h>
#include <myutil.h>
#include <myvalcounter.h>
#include <mathconst.h>
#include <mymacros.h>
#include <myfits.h>
#include <ccd_astro_const.h>
#include <mathdefs.h>
#include <mygraphcut.h>
#include <AstroCCD.h>
#include <sat_info.h>
#include <algorithm> // lower_bound on Ubuntu 

#include "ccd_image_creator.h"
#include "ccd_util.h"
#include "ccd_fastcalc.h"
#include "ccd_starcat.h"
#include "ccd_asastransform.h"


// global - to optimize but later it will be controled by precompiler
// directives 
LONG_T gNeighbList[MAX_CLUSTER_SIZE];

BOOL_T gFirstDump=TRUE;
BOOL_T CCD_Analyser::m_bInitialized=FALSE;

// processing 
CCDEventList CCD_Analyser::m_AllRunEvents;
CLongPoint CCD_Analyser::m_Aureola[MAX_CLUSTER_SIZE];
LONG_T CCD_Analyser::m_nPixelsInAureola=-1;

// shape to confirm type
CLongPoint CCD_Analyser::m_ShapeToConfirmMap[MAX_CLUSTER_SIZE];
LONG_T CCD_Analyser::m_ConfShapeCount=0;

// points to sum - calculated once at the begining !
CLongPoint CCD_Analyser::m_NeighbPoints[MAX_CLUSTER_SIZE];
LONG_T CCD_Analyser::m_SelfPos=0;
LONG_T CCD_Analyser::m_FirstDy=0;

CLongPoint CCD_Analyser::m_singlePoint[1];


// function pointers :
BOOL_T (*CCD_Analyser::m_pAnalFoundEventFunc)( void* event_info, int type, int ccd_idx, void* event_info2 ) = NULL;

// calculations :


void CCD_Analyser::SetCustomEventAnalFunc( BOOL_T (*pAnalFoundEventFunc)( void* event_info, int type, int ccd_idx, void* event_info2  ) )
{ 
	m_pAnalFoundEventFunc = pAnalFoundEventFunc; 
	CMyFit::m_pFillFunc = m_pAnalFoundEventFunc;
}



// 
// function calculates points in the shape 
// depending on the shapeType it will be :
//   1 - square (Redial is half of BOK ??? )
//   2 - circle 
//   3 - star ( example - cross when 5 pixels )
//  

int compareLongPoint( const void* left, const void* right)
{
	CLongPoint* pLeft = (CLongPoint*)left;
	CLongPoint* pRight = (CLongPoint*)right;

	if(pLeft->y < pRight->y){
		return -1;
	}else{
		if(pLeft->y > pRight->y){
			return 1;
		}else{
			if(pLeft->x < pRight->x){
				return -1;
			}else{
				if(pLeft->x > pRight->x){
					return 1;
				}else{
					return 0;
				}
			}
		}
	}
	
		
}

double CCD_Analyser::CalcClusterRatio( ELEM_TYPE* p_data, LONG_T pos,
                                       LONG_T* pixel_list, LONG_T list_count, LONG_T& pixel_cnt,
                                       BOOL_T bAboveTresholdOnly/*=FALSE*/, double Treshold/*=0*/,
													BOOL_T bUseMax/*=FALSE*/)
{
	double sum = 0, ret=0;
	pixel_cnt=0;

	if( p_data ){
		for(int p=0;p<list_count;p++){
			if(bAboveTresholdOnly){				
				if( p_data[pixel_list[p]] > Treshold ){
					sum += p_data[pixel_list[p]];
					pixel_cnt++;
				}			
			}else{
				sum += p_data[pixel_list[p]];
				pixel_cnt++;				
			}
		}				
		if(pixel_cnt){
			double value = p_data[pos];
			if(bUseMax){
				for( register int i=0;i<list_count;i++){
					if(p_data[pixel_list[i]] > value)
						value = p_data[pixel_list[i]];
				}
			}
			ret = (value-(sum/pixel_cnt))/(sum/pixel_cnt);
		}
	}
	return ret;
}

LONG_T CCD_Analyser::CalcShapePoints( CLongPoint* shapePoints, 
		  										  eConfShape_T shapeType,
												  double Redial )
{
	LONG_T Redial_Long = (LONG_T)Redial;
	int x,y;
	
	if(!shapePoints)
		shapePoints = m_NeighbPoints;

	int cnt = 0;
	for(x=-Redial_Long;x<=Redial_Long;x++){
		for(y=-Redial_Long;y<=Redial_Long;y++){		
			if( shapeType == shapeSquare || shapeType == shapeStar){
				shapePoints[cnt].x = x;
				shapePoints[cnt].y = y;	
				cnt++;
			}					
			if( shapeType == shapeCircle ){
				double r = sqrt(x*x+y*y);
				if(r<=Redial){
					shapePoints[cnt].x = x;
					shapePoints[cnt].y = y;	
					cnt++;
				}
			}
		}
	}


	// for shape type add also edge points - and sort 
	if( shapeType == shapeStar ){
		LONG_T step=1;		
		if(Redial_Long==0){
			shapePoints[1].x = -1;
			shapePoints[1].y = 0;

			shapePoints[2].x = 1;
         shapePoints[2].y = 0;

			shapePoints[3].x = 0;
			shapePoints[3].y = -1;

			shapePoints[4].x = 0;
         shapePoints[4].y = 1;
			cnt = 5;
		}else{
			LONG_T to_add = (2*Redial_Long+1)-2;

			// while there are still points to add
			while(to_add>0){		
				for(x=-Redial_Long+step;x<=Redial_Long-step;x++){				
					// add upper points
					shapePoints[cnt].x = x;
					shapePoints[cnt].y = Redial_Long + step;
					cnt++;

					// add lower points
					shapePoints[cnt].x = x;
					shapePoints[cnt].y = -Redial_Long - step;
					cnt++;
				}

				for(y=-Redial_Long+step;y<=Redial_Long-step;y++){				
					// add left points
					shapePoints[cnt].x = -Redial_Long - step;
					shapePoints[cnt].y = y;
					cnt++;

					// add lower points
					shapePoints[cnt].x = Redial_Long + step;
					shapePoints[cnt].y = y;
					cnt++;
				}	

				step++;
				to_add = to_add - 2;
			}
		}
	}				
	qsort( shapePoints, cnt, sizeof(CLongPoint), compareLongPoint );
	shapePoints[cnt].x = NEIGHB_LIST_END;
   shapePoints[cnt].y = NEIGHB_LIST_END;

	_TRACE_PRINTF_6("Dumping points in calculated shape :\n");
	for(int i=0;i<cnt;i++){
		if(shapePoints[i].x==0 && shapePoints[i].y==0)
			m_SelfPos = i;
		_TRACE_PRINTF_6("(%d,%d)\n",shapePoints[i].x,shapePoints[i].y);
	}

	// exit(0);	
	Assert(cnt<MAX_CLUSTER_SIZE,"Max size of cluster table exceeded");

	if(shapePoints==m_NeighbPoints){
		gCCDParams.m_nNeighbToSumCount = cnt;
		gCCDParams.RecalculateParams();
	}
	return cnt;
}

// prepare list of points around - with point itself :
void CCD_Analyser::CalcShapeToConfirm()
{
	m_ConfShapeCount = CalcShapePoints( m_ShapeToConfirmMap, gCCDParams.m_ConfShape, gCCDParams.m_ConfRedial );
}

// list of points around without point itself
void CCD_Analyser::CalcAurola()
{
	if(m_nPixelsInAureola<0){
		LONG_T PointsAroundToConfirm = gCCDParams.m_nPixelsAroundToConfirm;

		m_nPixelsInAureola = 0;
		for(int x=-PointsAroundToConfirm;x<=PointsAroundToConfirm;x++){
			for(int y=-PointsAroundToConfirm;y<=PointsAroundToConfirm;y++){
				double r = sqrt(x*x+y*y);
				if(r<=PointsAroundToConfirm && (x!=0 || y!=0)){
					// do not include center point : (x!=0 || y!=0)
					m_Aureola[m_nPixelsInAureola].x = x;
					m_Aureola[m_nPixelsInAureola].y = y;
					m_nPixelsInAureola++;
				}
			}
		}
	}
}

CCD_Analyser::CCD_Analyser( int sizeX, int sizeY ) 
: CBaseAnal( sizeX, sizeY ), m_pNeighMap(NULL),m_AllEvents(gCCDParams.GetEventsBufferSize()),
  m_pMaxLaplacePrev(NULL),m_pWrkTableInX(NULL), m_pCurrSatList(NULL), m_pAllSatList(NULL),
	m_eCurrFitTrackType( eNormalTrack )
{
	Initialize();
	if(gCCDParams.m_bAnalyseMaxFromPrevInPipeline && gCCDParams.m_bCheckLaplaceCondition){
		m_pMaxLaplacePrev = new Table2D<BIG_ELEM_TYPE>(m_SizeX,m_SizeY);
		m_pWrkTableInX = new BIG_ELEM_TYPE[ m_SizeX ];
	}
	m_pCurrSatList = new CSatList();	
	m_pAllSatList = new CSatList();
}

CCD_Analyser::~CCD_Analyser()
{
	if(m_pMaxLaplacePrev)
		delete m_pMaxLaplacePrev;
	if(m_pWrkTableInX)
		delete [] m_pWrkTableInX;
	if(m_pNeighMap)
		delete m_pNeighMap;
	if(m_pCurrSatList)
		delete m_pCurrSatList;
	if(m_pAllSatList)
		delete m_pAllSatList;
}

void CCD_Analyser::RefreshParams()
{
	m_bInitialized = FALSE;
	Initialize();
}

void CCD_Analyser::Initialize(){
	if(!m_bInitialized){
		gCCDParams.InitParams();
		m_bInitialized = TRUE;
	}
}

void CCD_Analyser::Init()
{
	Initialize();
	memset(gNeighbList,MAX_CLUSTER_SIZE*sizeof(LONG_T),'\0');

	if(gCCDParams.m_nPixelsAroundToConfirm){
		CalcAurola();
	}
	if(gCCDParams.m_bConfirmReq && gCCDParams.m_ConfShape!=shapeCluster){
		CalcShapeToConfirm();
	}

	// now shape is well steered from parameters level
	// points are calculated here :
	gCCDParams.m_nNeighbToSumCount = CalcShapePoints( m_NeighbPoints, gCCDParams.m_eNeigbShape, gCCDParams.m_dNeighbRedial );
	

	if(gCCDParams.m_bKeepNeighbMap){
	   MYTRACE3(gCCDTrace,"Initializing neighbours map table");
      m_pNeighMap = new CDescTab2D( m_SizeX, m_SizeY, gCCDParams.m_nNeighbToSumCount );
      m_pNeighMap->InitNeigb2D();
      MYTRACE3(gCCDTrace,"Neighbours table initialized");
   }	
	m_FirstDy = ((m_NeighbPoints[0].y)*m_SizeX);

	m_singlePoint[0].x = 0;
	m_singlePoint[0].y = 0;
}

void CCD_Analyser::ClearState()
{
	m_AllEvents.clear();
}

LONG_T* CCD_Analyser::GetAnalNeighbFromMap( long pos, long& ncnt )
{
	return (m_pNeighMap->GetDesc( pos, ncnt ));	
}


void CCD_Analyser::VerifyPrevValue(LONG_T& val)
{
	if(val>=gCCDParams.m_MaxAllowedVal){
		// do not accept values exceeding limit
		val = gCCDParams.m_MaxAllowedVal;
	}
}

void CCD_Analyser::VerifyValue(LONG_T& val,BOOL_T bMax/*=FALSE*/)
{
	if(val>=gCCDParams.m_MaxAllowedVal){
		// do not accept values exceeding limit
		if(bMax)
			val = gCCDParams.m_MaxAllowedVal;
		else
			val = 0;
	}
}


#define SET_NEXT_ELEM( tab, pos, val ) { tab[pos]=val;pos++; }

long CCD_Analyser::GetAnalNeighbOpt( long x,long y,long xSize,long ySize,
                   		             LONG_T* neighb_list, LONG_T& ncnt,
                                     LONG_T* out_list, LONG_T& ocnt )
{
	Initialize();
	LONG_T pos = x+y*xSize;
	ncnt = 0;
	ocnt = 0;
	switch(gCCDParams.m_nNeighbToSumCount){
		case 5:
		case 4 :
				/*    x
                 xxx
                  x    */
				{
					// add neighbours sorted !
					if(y>0){	
						neighb_list[ ncnt ] = pos-xSize;
						out_list[ ocnt ] = pos-xSize;						
						ncnt++;
						ocnt++;
					}
					if(x>0){
						neighb_list[ ncnt ] = pos-1;
						out_list[ ocnt ] = pos-1;						
						ncnt++;
						ocnt++;						
					}
					
					// curent point - only to neighb - not outer
					neighb_list[ ncnt ] = pos;
					ncnt++;

					if(x<(xSize-1)){
						neighb_list[ ncnt ] = pos+1;
						out_list[ ocnt ] = pos+1;						
						ncnt++;
						ocnt++;						
					}
					if(y<(ySize-1)){
						neighb_list[ ncnt ] = pos+xSize;
						out_list[ ocnt ] = pos+xSize;
						ncnt++;
						ocnt++;						
					}
				}
				break;
		case 8 :
		case 9 :
				/*   xxx
                 xxx
                 xxx    */
				{
					// add neighbours sorted !
					if(y>0){
						if(x>0){	
							SET_NEXT_ELEM( neighb_list, ncnt, pos-xSize-1)
						}
						SET_NEXT_ELEM( neighb_list, ncnt, pos-xSize )
						if(x<(xSize-1)){
							SET_NEXT_ELEM( neighb_list, ncnt, pos-xSize+1 )	
						}							
					}
					if(x>0){
						SET_NEXT_ELEM( neighb_list, ncnt, pos-1)
					}
					LONG_T PointPos = ncnt;
					SET_NEXT_ELEM( neighb_list, ncnt, pos)
					if(x<(xSize-1)){
						SET_NEXT_ELEM( neighb_list, ncnt, pos+1)
					}
					if(y<(ySize-1)){
						if(x>0){
							SET_NEXT_ELEM( neighb_list, ncnt, pos+xSize-1)
						}
						SET_NEXT_ELEM( neighb_list, ncnt, pos+xSize)
						if(x<(xSize-1)){
							SET_NEXT_ELEM( neighb_list, ncnt, pos+xSize+1)
						}
					}
					if(PointPos)
						memcpy(out_list,neighb_list,PointPos*sizeof(LONG_T));
					memcpy(out_list+PointPos,neighb_list+PointPos+1,(ncnt-PointPos-1)*sizeof(LONG_T));
					ocnt = ncnt-1;
				}
				break;
			
		default :
				{
					// only point itself
					SET_NEXT_ELEM( neighb_list, ncnt, pos )
					SET_NEXT_ELEM( out_list, ocnt, pos)
				}
				break;
	}
	return ncnt;
}


void CCD_Analyser::AutoCalculateIgnoreEdges()
{
	double dx = gCCDParams.m_FrameDX;
	double dy = gCCDParams.m_FrameDY;

	int nPipelineSize = gCCDParams.m_nPipelineSize;

	double left = gCCDParams.m_nIgnoreEdgeLeft;
	double right = gCCDParams.m_nIgnoreEdgeRight;
	double up = gCCDParams.m_nIgnoreEdgeUp;
	double bottom = gCCDParams.m_nIgnoreEdgeBottom;

/*	if(fabs(dx)>0.5){
		double edge_val = fabs(dx*nPipelineSize);
		if(dx>0){
			left = MAX( left, edge_val );
		}else{
			right = MAX( right, edge_val );
		}
			
	}
	if(fabs(dy)>0.5){
		double edge_val = fabs(dy*nPipelineSize);
		if(dy>0){
			bottom = MAX( bottom, edge_val );
		}else{
			up = MAX( up, edge_val );
		}			
	}
*/

	gCCDParams.SetIgnoreEdges( left, right , bottom, up );
}

void CCD_Analyser::AutoCalculateIgnoreEdges( CCDPipeline* pPipeline )
{
	Assert(pPipeline!=NULL,"Pipeline pointer must not be NULL");

	double dx = pPipeline->GetPipelineCfg().m_FrameDX;
	double dy = pPipeline->GetPipelineCfg().m_FrameDY;

	int nPipelineSize = gCCDParams.m_nPipelineSize;
	if(pPipeline && pPipeline->GetPipelineSize()>0)
		nPipelineSize = pPipeline->GetPipelineSize();

	double left = gCCDParams.m_nIgnoreEdgeLeft;
	double right = gCCDParams.m_nIgnoreEdgeRight;
	double up = gCCDParams.m_nIgnoreEdgeUp;
	double bottom = gCCDParams.m_nIgnoreEdgeBottom;

	left = pPipeline->GetPipelineCfg().m_nIgnoreEdgeLeft;
	right = pPipeline->GetPipelineCfg().m_nIgnoreEdgeRight;
	bottom = pPipeline->GetPipelineCfg().m_nIgnoreEdgeBottom;
	up = pPipeline->GetPipelineCfg().m_nIgnoreEdgeUp;		
	
/*	if(fabs(dx)>0.5){
		double edge_val = fabs(dx*nPipelineSize);
		if(dx>0){
			left = MAX( left, edge_val );
		}else{
			right = MAX( right, edge_val );
		}
			
	}
	if(fabs(dy)>0.5){
		double edge_val = fabs(dy*nPipelineSize);
		if(dy>0){
			bottom = MAX( bottom, edge_val );
		}else{
			up = MAX( up, edge_val );
		}			
	}
*/

	pPipeline->GetPipelineCfg().SetIgnoreEdges( left, right , bottom, up );
}

long CCD_Analyser::GetAnalNeighbOpt( long x,long y, long pos,long xSize,long ySize,
                   		             LONG_T* neighb_list, LONG_T& ncnt,
                                     BOOL_T bAddSelf )
{
	ncnt=0;
	if( gCCDParams.m_dNeighbRedial < gCCDParams.m_nIgnoreEdge ){
		// optimized version - no checking if x>0, x<xSize etc ... 
		// cause IgnoreEdge ensures this !		
		for(register int i=0;i<gCCDParams.m_nNeighbToSumCount;i++){
			if(bAddSelf || i!=m_SelfPos){
				neighb_list[ncnt] = 	pos + ((m_NeighbPoints[i].y)*xSize) + m_NeighbPoints[i].x;
				ncnt++;
			}
		}
	}else{		
		for(register int i=0;i<gCCDParams.m_nNeighbToSumCount;i++){
			if(bAddSelf || i!=m_SelfPos){
				LONG_T x0 = x+m_NeighbPoints[i].x;
				LONG_T y0 = y+m_NeighbPoints[i].y;
				if( x0>=0 && x0<xSize && y0>=0 && y0<ySize ){
					neighb_list[ncnt] = 	pos + ((m_NeighbPoints[i].y)*xSize) + m_NeighbPoints[i].x;
					ncnt++;
				}
			}
		}
	}
	return ncnt;
}


long CCD_Analyser::GetAnalNeighbWithSelf( long x,long y,long pos, long xSize,long ySize,
                    		          LONG_T* neighb_list, LONG_T& ncnt )
{
	ncnt=0;

	// optimized version - no checking if x>0, x<xSize etc ... 
	// cause IgnoreEdge ensures this !
	LONG_T d_pos_y = m_FirstDy;		
	for(register int i=0;i<gCCDParams.m_nNeighbToSumCount;i++){			
		LONG_T x0 = x+m_NeighbPoints[i].x;
      LONG_T y0 = y+m_NeighbPoints[i].y;
		
		if(i){
			if(m_NeighbPoints[i-1].y!=m_NeighbPoints[i].y)
				d_pos_y += xSize;
		}
		if(x0>=0 && x0<xSize && y0>=0 && y0<ySize ){
			neighb_list[ncnt] = 	pos + d_pos_y + m_NeighbPoints[i].x;
			ncnt++;
		}
	}
	return ncnt;
}

long CCD_Analyser::GetAllNeighbNoSelf( long x0,long y0,long pos, long xSize,long ySize,
                                        LONG_T* neighb_list, LONG_T& ncnt )
{
	ncnt=0;
	for(register int y=y0-1;y<=y0+1;y++){
		for(register int x=x0-1;x<=x0+1;x++){
			if(x>=0 && y>=0 && x<xSize && y<ySize && ( x!=x0 || y!=y0 ) ){
				neighb_list[ncnt] = y*xSize+x;
				ncnt++;
			}
		}
	}
	return ncnt;
}

long CCD_Analyser::GetAnalNeighbNoSelf( long x,long y,long pos, long xSize,long ySize,
                          	             LONG_T* neighb_list, LONG_T& ncnt )
{
	ncnt=0;

	// optimized version - no checking if x>0, x<xSize etc ... 
	// cause IgnoreEdge ensures this !
	LONG_T d_pos_y = m_FirstDy;		
	for(int i=0;i<gCCDParams.m_nNeighbToSumCount;i++){			
		if(i){
			if(m_NeighbPoints[i-1].y!=m_NeighbPoints[i].y)
				d_pos_y += xSize;
		}

		if(i!=m_SelfPos){
			LONG_T x0 = x+m_NeighbPoints[i].x;
	      LONG_T y0 = y+m_NeighbPoints[i].y;

			if(x0>=0 && x0<xSize && y0>=0 && y0<ySize ){		
				neighb_list[ncnt] = 	pos + d_pos_y + m_NeighbPoints[i].x;
				ncnt++;
			}
		}
	}
	return ncnt;
}





long CCD_Analyser::GetAnalNeighb( long x,long y,long xSize,long ySize,
                   		          CLongList& Neighbours, CLongList& OuterNeighbours )
{
	Initialize();
	LONG_T pos = x+y*xSize;
	Neighbours.Clear();
	switch(gCCDParams.m_nNeighbToSumCount){
		case 5:
		case 4 :
				/*    x
                 xxx
                  x    */
				{
					// add neighbours sorted !
					if(y>0)	
					   Neighbours.Add( pos-xSize );
					if(x>0)	
						Neighbours.Add( pos-1 );
					OuterNeighbours = Neighbours;
					Neighbours.Add( pos );
					if(x<(xSize-1)){
						Neighbours.Add( pos+1 );
						OuterNeighbours.Add( pos+1 );
					}
					if(y<(ySize-1)){
						Neighbours.Add( pos+xSize );
						OuterNeighbours.Add( pos+xSize );
					}
				}
				break;
		case 8 :
		case 9 :
				/*   xxx
                 xxx
                 xxx    */
				{
					// add neighbours sorted !
					if(y>0){
						if(x>0)
							Neighbours.Add( pos-xSize-1 );
						Neighbours.Add( pos-xSize );
						if(x<(xSize-1))
							Neighbours.Add( pos-xSize+1 );
					}
					if(x>0)	
						Neighbours.Add( pos-1 );
					OuterNeighbours = Neighbours;					
					Neighbours.Add( pos );					
					if(x<(xSize-1)){
						Neighbours.Add( pos+1 );
						OuterNeighbours.Add( pos+1 );
					}
					if(y<(ySize-1)){
						if(x>0){
							Neighbours.Add( pos+xSize-1 );
							OuterNeighbours.Add( pos+xSize-1 );
						}
						Neighbours.Add( pos+xSize );
						OuterNeighbours.Add( pos+xSize );
						if(x<(xSize-1)){
							Neighbours.Add( pos+xSize+1 );
							OuterNeighbours.Add( pos+xSize+1 );
						}
					}
				}
				break;
			
		default :
				{
					// only point itself
					Neighbours.Add( pos );		
					OuterNeighbours = Neighbours;			
				}
				break;
	}
	return Neighbours.size();
}


long CCD_Analyser::GetNeighbours( long pos, long xSize,long ySize, 
                                  CLongList& Neighbours, BOOL_T bAddSelf )
{
	Initialize();
	long x = ( pos % xSize );
   long y =  ( pos / xSize );	
	long ret = GetNeighbPositions(x ,y , xSize, ySize, Neighbours, bAddSelf );
	return ret;
}

long CCD_Analyser::GetNeighbPositions( long x,long y,long xSize,long ySize,
                   		       CLongList& Neighbours, BOOL_T bAddSelf )
{
	Initialize();
	LONG_T pos = x+y*xSize;
	Neighbours.RemoveAll();
	if(y>0)
      Neighbours.Add(pos-xSize);
	if(x>0)
		Neighbours.Add(pos-1);
	if(bAddSelf)
		Neighbours.Add(pos);
	if(x<(xSize-1))
		Neighbours.Add(pos+1);
	if(y<(ySize-1))
		Neighbours.Add(pos+xSize);
	return Neighbours.size();
}

long CCD_Analyser::GetMaxLaplaceOnPrevRotCorrect( const CPixelAnalyseIn& in )
{
	long MaxLaplace=-100000;


	// loop is starting from OLDEST frame !
	for(register int i=0;i<in.PrevMatrixPtrCnt;i++){		
		// double frame_dx = -(in.PrevMatrixPtrCnt-i)*gCCDParams.m_FrameDX;
		// double frame_dy = -(in.PrevMatrixPtrCnt-i)*gCCDParams.m_FrameDY;
		double frame_dx = in.PrevFramesShiftTab[i].frameDX;
		double frame_dy = in.PrevFramesShiftTab[i].frameDY;

		ELEM_TYPE** p_data_fast = in.PrevMatrixPtr[i]->get_data_buffer_fast();
		BIG_ELEM_TYPE** p_laplace = ((CCDMatrix*)in.PrevMatrixPtr[i])->get_frame_laplace_fast();

		if(p_laplace){
			for(register int v=0;v<gCCDParams.m_nVetoPointsCount;v++){
		      register long x = (long)(in.x+(gCCDParams.m_VetoArea)[v].x + frame_dx);
   	      register long y = (long)(in.y+(gCCDParams.m_VetoArea)[v].y + frame_dy);

				// printf("frame : %d at (%d,%d)=%d\n",i,x,y,p_laplace[y][x]);
	
   	      if(MaxLaplace<p_laplace[y][x]){
      	         MaxLaplace = p_laplace[y][x];
         	}
			}
		}else{
			for(register int v=0;v<gCCDParams.m_nVetoPointsCount;v++){
		      register long x = (long)(in.x+(gCCDParams.m_VetoArea)[v].x + frame_dx);
   	      register long y = (long)(in.y+(gCCDParams.m_VetoArea)[v].y + frame_dy);
	
				long laplace = Table2D<ELEM_TYPE>::CalcLaplaceSum( x, y, in.xSize, p_data_fast, gCCDParams.m_eLaplaceType );
	
   	      if(MaxLaplace<laplace){
      	         MaxLaplace = laplace;
         	}
			}
		}
	}
	/*if(MaxLaplace!=in.p_max_prev_lap[in.y][in.x]){
		printf("(x,y)=(%d,%d) %d!=%d\n",in.x,in.y,(long)MaxLaplace,(long)in.p_max_prev_lap[in.y][in.x]);
		for(register int i=0;i<in.PrevMatrixPtrCnt;i++){
			printf("\n\nframe : %d\n",i);
			double frame_dx = in.PrevFramesShiftTab[i].frameDX;
	      double frame_dy = in.PrevFramesShiftTab[i].frameDY;
			ELEM_TYPE** p_laplace = ((CCDMatrix*)in.PrevMatrixPtr[i])->get_frame_laplace_fast();

			// for(register int v=0;v<gCCDParams.m_nVetoPointsCount;v++){


        //  register long x = in.x+(gCCDParams.m_VetoArea)[v].x + frame_dx;
        //  register long y = in.y+(gCCDParams.m_VetoArea)[v].y + frame_dy;

			long x_pos = in.x+frame_dx;
			long y_pos = in.y+frame_dy;
			for(register long y=y_pos-5;y<y_pos+5;y++){
				for(register long x=x_pos-5;x<x_pos+5;x++){
					if(x>=0 && y>=0){
						// printf("%d ",p_laplace[y][x]);
						printf("(%d,%d)=%d ",x,y,p_laplace[y][x]);
						if(p_laplace[y][x]==MaxLaplace)
							printf("   !!!!!   ");
					}
				}
				printf("\n");
			}
		}

	 	Assert(FALSE,"In correctly calculated max laplace prev matrix %d!=%d at (%d,%d)",(long)MaxLaplace,(long)in.p_max_prev_lap[in.y][in.x],in.x,in.y);	
	}*/

	return MaxLaplace;	
}

LONG_T CCD_Analyser::GetMaxAverageInVetoArea( CPixelAnalyseIn& in, BOOL_T bLaplace )
{
	register long maxAverage = -10000;
	for(register int v=0;v<gCCDParams.m_nVetoPointsCount;v++){
       register long x = (long)(in.x+(gCCDParams.m_VetoArea)[v].x);
       register long y = (long)(in.y+(gCCDParams.m_VetoArea)[v].y);
	
		 register long prevAverage = CalcAverageOfPrevN( x, y, in,  
																		 gCCDParams.m_nMaxOfAverageOfPrevN,
																		 bLaplace );
		 if( prevAverage>maxAverage ){
		     maxAverage = prevAverage;
		 }
	}		
	return maxAverage;
}

LONG_T CCD_Analyser::CalcAverageOfPrevN( CPixelAnalyseIn& in, long prevCount, 
                                         BOOL_T bLaplace/*=FALSE*/,
													  CLongPoint* shapePoints/*=NULL*/, long pointCnt/*=0*/ )
{
	return CalcAverageOfPrevN(in.x, in.y, in, prevCount, bLaplace );
}


LONG_T CCD_Analyser::CalcWeightedAverageOfPrevN( const CPixelAnalyseIn& in, long prevCount,
																InfoTable2D* pTotalShiftInfo, BOOL_T bLaplace )
{
	long prev_x = in.x;
	long prev_y = in.y;

	register long sum = 0;
	register long f=1;
	register long cnt=0;
	

	static CLongPoint overlapedPixels[4];


	for(f=1;f<=prevCount;f++){
		Area2DInfo& pShiftInfo = pTotalShiftInfo[f].GetAreaDesc( in.x , in.y );

		double dx_to_1 = pShiftInfo.m_LocalShiftInfo.m_LocalDX;
		double dy_to_1 = pShiftInfo.m_LocalShiftInfo.m_LocalDY;

		long prev_x_0 = (long)(in.x + dx_to_1);
		long prev_y_0 = (long)(in.y + dy_to_1);
	

		long prev_x_0_1 = prev_x_0+1;
		long prev_y_0_1 = prev_y_0+1;

		double pixel_sum = 0;

		if(prev_x_0>=0 && prev_x_0<in.xSize && prev_y_0>=0 && prev_y_0<in.ySize )
			pixel_sum += ((in.PrevMatrixPtr[f])->m_pFastData[prev_y_0][prev_x_0])*pShiftInfo.m_LocalShiftInfo.m_LocalS0;
		if(prev_x_0>=0 && prev_x_0<in.xSize && prev_y_0_1>=0 && prev_y_0_1<in.ySize )
			pixel_sum += ((in.PrevMatrixPtr[f])->m_pFastData[prev_y_0_1][prev_x_0])*pShiftInfo.m_LocalShiftInfo.m_LocalS1;
		if(prev_x_0_1>=0 && prev_x_0_1<in.xSize && prev_y_0_1>=0 && prev_y_0_1<in.ySize )
			pixel_sum += ((in.PrevMatrixPtr[f])->m_pFastData[prev_y_0_1][prev_x_0_1])*pShiftInfo.m_LocalShiftInfo.m_LocalS2;
		if(prev_x_0_1>=0 && prev_x_0_1<in.xSize && prev_y_0>=0 && prev_y_0<in.ySize )
			pixel_sum += ((in.PrevMatrixPtr[f])->m_pFastData[prev_y_0][prev_x_0_1])*pShiftInfo.m_LocalShiftInfo.m_LocalS3;


		sum += (long)pixel_sum;
		cnt++;
	}
	if(cnt)
		sum = (sum / cnt);
	return sum;
	
}

void CCD_Analyser::GetPrevPixelPos( int& prev_x, int& prev_y, int x, int y,
												int stepsBack, const CPixelAnalyseIn& in )
{
	if( gCCDParams.m_bUseFrameShift ){
		double tmp_x = x;
		double tmp_y = y;
		for(int i=0;i<stepsBack;i++){
			CCDMatrix* pImage = in.PrevMatrixPtr[i];
			if( pImage ){
				tmp_x = tmp_x - pImage->m_dx;
				tmp_y = tmp_y - pImage->m_dy;
			}
		}
		prev_x = my_round(tmp_x);
		prev_y = my_round(tmp_y);
		return;
	}

	double curr_x=x, curr_y=y;
	CCDDataResults::CalcStarPositionAuto( (double)x, (double)y, 0,
                                     -in.PrevFramesTime[stepsBack], -stepsBack,
                                     curr_x, curr_y,
                                     in.pCCDInfo, &((in.pPipeline)->m_PipelineCfg) );
	
	prev_x = my_round(curr_x);
	prev_y = my_round(curr_y);
}	

void CCD_Analyser::GetPrevPixelPos( int x, int y, int* PrevFramesX, int* PrevFramesY,
												int backStart, int backTo, const CPixelAnalyseIn& in )
{
	if( gCCDParams.m_bShiftUsesAstroFormulas ){
		CCDDataResults::CalcStarPositionsFromFormula( x, y, PrevFramesX, PrevFramesY, backStart, backTo,
															in.PrevFramesTime, -1, in.pCCDInfo, 
															&((in.pPipeline)->m_PipelineCfg) );
	}else{
		// normal way :
		for(int f=backStart;f<=backTo;f++){
			GetPrevPixelPos( PrevFramesX[f], PrevFramesY[f], x, y, f, in );	
		}
	}
}


/*
void CCD_Analyser::GetPrevPixelPosFromFormula( int& prev_x, int& prev_y, long x, long y,
												int stepsBack, const CPixelAnalyseIn& in, int sec )
{
	double prev_x_d,prev_y_d;
	(in.pCCDInfo)->xyAfterTimeNEW( (double)x, (double)y, sec, prev_x_d,prev_y_d );
	prev_x = my_round(prev_x_d);
	prev_y = my_round(prev_y_d);
}*/

void CCD_Analyser::CalcPrevPositionWithRot( long& prev_x, long& prev_y, long x0, long y0, 
														  int stepsBack, const CPixelAnalyseIn& in )
{
   double x0_prim = x0-gCCDParams.m_RotCenterX;
   double y0_prim = y0-gCCDParams.m_RotCenterY;	
	double cos_alfa_0 = x0_prim/(sqrt(x0_prim*x0_prim+y0_prim*y0_prim));
   double alfa_0 = acos(cos_alfa_0 );

	if(y0_prim<0){
		alfa_0 = (TWO_PI - alfa_0);
	}

	double r = sqrt(x0_prim*x0_prim + y0_prim*y0_prim);

	// double dAlfa1 = gCCDParams.m_RotValueDAlfa*stepsBack;
	double dAlfa = EARTH_ROTATION_OMEGA*in.PrevFramesTime[stepsBack]*gCCDParams.m_SinOfGeoLatitude;

	prev_x = my_round( gCCDParams.m_RotCenterX + r*cos(alfa_0-dAlfa) );
	prev_y = my_round( gCCDParams.m_RotCenterY + r*sin(alfa_0-dAlfa) );		


	// printf("(%d,%d) Frame time = %f dAlfaOld=%f dAlfaNew=%f\n",x0,y0,frameTime,dAlfa1,dAlfa);
	
}


// very slow :
int CCD_Analyser::CalcAverageOfPrevN( int x, int y, 
												  int xSize, int ySize,
												  int prevCount, 
												  BOOL_T bLaplace, CCDPipeline* pPipeline )
{
	CPixelAnalyseIn in;
	in.xSize = xSize;
	in.ySize = ySize;
	in.treshold_for_max = -1000;
	in.pPipeline = pPipeline;
	in.pipeline_size_minus_1 = (pPipeline->GetPipelineSize()-1);
	in.Matrix = &((pPipeline->GetCurrent())[0]);
	in.ccd_index = 0;	
	in.pCCDInfo = &(pPipeline->GetCCDInfoTab()[in.ccd_index]);
	in.p_data = (in.Matrix)->get_data_buffer();
	in.p_data_fast = (in.Matrix)->get_data_buffer_fast();
	in.p_curr_data_laplace = (in.Matrix)->get_frame_laplace_fast();
	in.p_curr_data_laplace_normal = (in.Matrix)->get_frame_laplace();
	in.pCamCfg = (CCcdCfg*)((in.pPipeline)->GetCamCfgTab()[in.ccd_index]);
	in.x = x;
	in.y = y;
	in.pos = y*xSize+x;
	in.PrevMatrixPtrCnt = (in.pPipeline)->GetAllMatrixPtrsChronologicalInt( 0, in, TRUE, prevCount+1 );
	
	int ret = CalcAverageOfPrevN( x, y, in, prevCount, bLaplace );
	return ret;
}


LONG_T CCD_Analyser::CalcAverageOfPrevNInPixel( const long x, const long y,
                        CPixelAnalyseIn& in, long prevCount,
                        BOOL_T bLaplace/*=FALSE*/ )
{
   register int sum = 0;
   register int f=1;

   for(f=1;f<=prevCount;f++){
      if(bLaplace){
         sum += (in.PrevLaplacePtr[f])[y][x];
      }else{
         sum += (in.PrevMatrixPtr[f])->m_pFastData[y][x];
      }
   }
   if( prevCount>=1 ){
      sum = (sum / prevCount);
   }

   return sum;
}


// BE CARFULL HERE :
// 
// in case you want to use table PrevValues for average calculation
// do not sum in loop from 1 to prevCount because some points 
// can be outside the frame and are not counted !!!!!!!!
// thus you should make loop up to cnt !!!!
// assuming on all frames older then cnt it is also outside the image 
// or prepare class - CFastLongTable
// which will store PrevValues and counter so that it is really ease
// and sure 
//
LONG_T CCD_Analyser::CalcAverageOfPrevN( const long x, const long y,
														CPixelAnalyseIn& in, long prevCount,
														BOOL_T bLaplace/*=FALSE*/ )
{
	register long prev_x = x;
	register long prev_y = y;

	register int sum = 0;
	register int f=1;
	register int cnt=0;
	double curr_x=x, curr_y=y;
	

	register int deciding_value = (in.PrevMatrixPtr[0])->m_pFastData[in.y][in.x];
	if( bLaplace ){
		deciding_value = (in.PrevLaplacePtr[0])[y][x];
		//deciding_value = CBaseAnal::CalcLaplaceSum( in.x, in.y, in.xSize,
		//								(in.PrevMatrixPtr[0])->m_pFastData,
		//								gCCDParams.m_eLaplaceType );
	}

	// in.PrevFramesMaxVal[0] = 0;

	in.PrevFramesX[0] = x;
	in.PrevFramesY[0] = y;
	in.PrevValues.counter = 0;


	GetPrevPixelPos( x, y, in.PrevFramesX, in.PrevFramesY, 1, prevCount, in );
	for(f=1;f<=prevCount;f++){
			// GetPrevPixelPos( prev_x, prev_y, x, y, f, in );
			/*CCDDataResults::CalcStarPositionAuto( (double)x, (double)y, 0, 
												 -in.PrevFramesTime[f], -f,
												 curr_x, curr_y, (in.pCamCfg->m_CCDParams).m_RotCenterX,
												 (in.pCamCfg->m_CCDParams).m_RotCenterY, gCCDParams.m_SinOfGeoLatitude, 
												 (in.pCamCfg->m_CCDParams).m_FrameDX, (in.pCamCfg->m_CCDParams).m_FrameDY,
												 (in.pCamCfg->m_CCDParams).m_bUseRotInAverageOfPrevN, 
												 (in.pCamCfg->m_CCDParams).m_RotValueDAlfa,
												 in.pCCDInfo, (in.pCamCfg->m_CCDParams).m_RotValueDAlfaPerSec, 
												 &((in.pPipeline)->m_PipelineCfg.m_CCDParams) );*/
			/*CCDDataResults::CalcStarPositionAuto( (double)x, (double)y, 0,
                                     -in.PrevFramesTime[f], -f,
                                     curr_x, curr_y,
												 in.pCCDInfo, &((in.pPipeline)->m_PipelineCfg.m_CCDParams) );	
			prev_x = my_round(curr_x);
			prev_y = my_round(curr_y);*/

			// in.PrevFramesX[f] = prev_x;
			// in.PrevFramesY[f] = prev_y;
			prev_x = in.PrevFramesX[f];
			prev_y = in.PrevFramesY[f];
	

			if(prev_x>=0 && prev_y>=0 && prev_x<in.xSize && prev_y<in.ySize){
				if(gCCDParams.m_bCalcMaxNeighbHomeo){
					register long max_val = (in.PrevLaplacePtr[f])[prev_y][prev_x];

					register long low_x = prev_x;
					register long up_x = prev_x;
					register long low_y = prev_y;
					register long up_y = prev_y;
					if( deciding_value > in.treshold_for_max){
						low_x -= (in.pPipeline)->m_PipelineCfg.m_nCalcMaxNieghbRedial;
						low_y -= (in.pPipeline)->m_PipelineCfg.m_nCalcMaxNieghbRedial;
						up_x  += (in.pPipeline)->m_PipelineCfg.m_nCalcMaxNieghbRedial;
						up_y  += (in.pPipeline)->m_PipelineCfg.m_nCalcMaxNieghbRedial;
					}


					if(bLaplace){
						for(register int xx=low_x;xx<=up_x;xx++){
							for(register int yy=low_y;yy<=up_y;yy++){	
								if(xx>=0 && yy>=0 && xx<in.xSize && yy<in.ySize){
									if((in.PrevLaplacePtr[f])[yy][xx]>max_val)
										max_val = (in.PrevLaplacePtr[f])[yy][xx];
								}
							}
						}
					}else{
						for(register int xx=low_x;xx<=up_x;xx++){
							for(register int yy=low_y;yy<=up_y;yy++){	
								if(xx>=0 && yy>=0 && xx<in.xSize && yy<in.ySize){
									if((in.PrevMatrixPtr[f])->m_pFastData[yy][xx] > max_val ){
										max_val = (in.PrevMatrixPtr[f])->m_pFastData[yy][xx];								
									}
								}
							}
						}
					}
					in.PrevValues.Add( max_val );
					sum += max_val;	
					cnt++;
				}else{
					if(bLaplace){				
						sum += (in.PrevLaplacePtr[f])[prev_y][prev_x];
						in.PrevValues.Add( (in.PrevLaplacePtr[f])[prev_y][prev_x] );
					}else{
						sum += (in.PrevMatrixPtr[f])->m_pFastData[prev_y][prev_x];
						in.PrevValues.Add( (in.PrevMatrixPtr[f])->m_pFastData[prev_y][prev_x] );
					}
					cnt++;
				}
			}
	}
	if(cnt)
		sum = (sum / cnt);


/*	TESTED - this can be uncommented and test_val can be used , in 
   option to rejecte MAX_AND_MIN of prev added (odrzucenie skrajnych)
	int test_val = 0;
	for(register int k=0;k<in.PrevValues.counter;k++){
		test_val += in.PrevValues.values[k];
	}
	test_val = (test_val/in.PrevValues.counter);
	if(test_val!=sum){
		printf("ERROR !!!! old=%d, new=%d\n",sum,test_val);
		printf("NEW table (%d): ",in.PrevValues.counter);
		for(register int k=0;k<in.PrevValues.counter;k++){
			printf("%d,",in.PrevValues.values[k]);
		}
		printf("\n\n");
		exit(0);
	}	*/

	return sum;
}



LONG_T CCD_Analyser::CalcAverageOfPrevNRecalc( const long x, const long y,
														CPixelAnalyseIn& in, long prevCount,
														BOOL_T bLaplace/*=FALSE*/,
													   CLongPoint* shapePoints/*=NULL*/, long pointCnt/*=0*/,
														BOOL_T bPrintMax  )
{
	long prev_x = x;
	long prev_y = y;

	register long sum = 0;
	register long f=1;
	register long cnt=0;
	
	// if(!shapePoints){
	shapePoints = m_singlePoint;
	pointCnt = 1;
	//}

	if(x<3 || y<3 || y>(in.ySize-3) || x>(in.xSize-3))
		return 0;

	register long deciding_value = (in.PrevMatrixPtr[0])->m_pFastData[in.y][in.x];
	if( bLaplace ){
		//deciding_value = (in.PrevLaplacePtr[0])[y][x];
		deciding_value = Table2D<ELEM_TYPE>::CalcLaplaceSum( in.x, in.y, in.xSize,
										(in.PrevMatrixPtr[0])->m_pFastData,
										gCCDParams.m_eLaplaceType );
	}

	for(f=1;f<=prevCount;f++){
		for(register long p=0;p<pointCnt;p++){

			prev_x = my_round(x + in.PrevFramesShiftTab[f].frameDX + shapePoints[p].x);
			prev_y = my_round(y + in.PrevFramesShiftTab[f].frameDY + shapePoints[p].y);
	

			if(prev_x>=0 && prev_y>=0 && prev_x<in.xSize && prev_y<in.ySize){
				if(gCCDParams.m_bCalcMaxNeighbHomeo){
					register long max_val = (in.PrevLaplacePtr[f])[prev_y][prev_x];
					register int max_x=0;
					register int max_y=0;

					register long low_x = prev_x;
					register long up_x = prev_x;
					register long low_y = prev_y;
					register long up_y = prev_y;
					if( deciding_value > in.treshold_for_max){
						low_x -= (in.pPipeline)->m_PipelineCfg.m_nCalcMaxNieghbRedial;
						low_y -= (in.pPipeline)->m_PipelineCfg.m_nCalcMaxNieghbRedial;
						up_x  += (in.pPipeline)->m_PipelineCfg.m_nCalcMaxNieghbRedial;
						up_y  += (in.pPipeline)->m_PipelineCfg.m_nCalcMaxNieghbRedial;
					}
					for(register int xx=low_x;xx<=up_x;xx++){
						for(register int yy=low_y;yy<=up_y;yy++){
	
							if(xx>=3 && yy>=3 && xx<(in.xSize-3) && yy<(in.ySize-3)){
								if(bLaplace){				
									/*Assert( in.PrevLaplacePtr[f]!=NULL,"laplace frame of frame %f not calculated !",f);
									if( (in.PrevLaplacePtr[f])[yy][xx] > max_val ){
										max_val = (in.PrevLaplacePtr[f])[yy][xx];
									}*/

									//register int laplace = (in.PrevLaplacePtr[f])[yy][xx];	

									register int laplace = Table2D<ELEM_TYPE>::CalcLaplaceSum( xx, yy, in.xSize,
																			 (in.PrevMatrixPtr[f])->m_pFastData,
																			 gCCDParams.m_eLaplaceType );
									if(laplace>max_val){
										max_val = laplace;
										max_x = (xx-prev_x);
										max_y = (yy-prev_y);
									}
								}else{
									if((in.PrevMatrixPtr[f])->m_pFastData[yy][xx] > max_val ){
										max_val = (in.PrevMatrixPtr[f])->m_pFastData[yy][xx];								
										max_x = (xx-prev_x);
                              max_y = (yy-prev_y);
									}
								}
							}
						}
					}
					sum += max_val;
					cnt++;

					if(bPrintMax && deciding_value > in.treshold_for_max)
						_TRACE_PRINTF_2("MAX_POS : %d %d\n",max_x,max_y);
				}else{
					if(bLaplace){				
						/*Assert( in.PrevLaplacePtr[f]!=NULL,"laplace frame of frame %f not calculated !",f);*/
						// sum += (in.PrevLaplacePtr[f])[prev_y][prev_x];
						register int laplace = Table2D<ELEM_TYPE>::CalcLaplaceSum( prev_x, prev_y, in.xSize,
																(in.PrevMatrixPtr[f])->m_pFastData,
																gCCDParams.m_eLaplaceType );
						sum += laplace;
					}else{
						sum += (in.PrevMatrixPtr[f])->m_pFastData[prev_y][prev_x];	
					}
					cnt++;
				}
			} // if prev_x,prev_y in frames
		}
	}
	if(cnt)
		sum = (sum / cnt);
	return sum;
}

BOOL_T CCD_Analyser::VerifyIfNotMoreThenNExceedsTreshold( const CPixelAnalyseIn& in,
                                                   CPixelAnalyseOut& out,
																	LONG_T prevCount,
																	LONG_T nAllowedToExceed,
                                                   BOOL_T bLaplace,
																	CLongPoint* shapePoints/*=NULL*/, long pointCnt/*=0*/ )
{
	register int prev_x = in.x;
	register int prev_y = in.y;
	register long f=1;
   register long cntAbove=0;
	
	if(!shapePoints){
   	shapePoints = m_singlePoint;
	   pointCnt = 1;
   }

	if(bLaplace){
		for(f=1;f<=prevCount;f++){
			prev_x = my_round(in.x + in.PrevFramesShiftTab[f].frameDX);
			prev_y = my_round(in.y + in.PrevFramesShiftTab[f].frameDY);

			for(register long p=0;p<pointCnt;p++){
				register int xx = prev_x + shapePoints[p].x;
				register int yy = prev_y + shapePoints[p].y;

				long laplace = Table2D<ELEM_TYPE>::CalcLaplaceSum( xx, yy, in.xSize,
																		(in.PrevMatrixPtr[f])->m_pFastData,
																		gCCDParams.m_eLaplaceType );
				if(laplace > in.treshold_for_not_more_then_n_exceeds){
					cntAbove++;
					break;
				}
			}
		}		
	}else{
		for(f=1;f<=prevCount;f++){
			prev_x = my_round(in.x + in.PrevFramesShiftTab[f].frameDX);
			prev_y = my_round(in.y + in.PrevFramesShiftTab[f].frameDY);

			for(register long p=0;p<pointCnt;p++){
				register int xx = prev_x + shapePoints[p].x;
				register int yy = prev_y + shapePoints[p].y;

				if((in.PrevMatrixPtr[f])->m_pFastData[yy][xx] > in.treshold_for_not_more_then_n_exceeds){
					cntAbove++;
					break;
				}
			}
		}		
	}
	out.m_PixelOut.m_nCountAbove = cntAbove;
	if(cntAbove>gCCDParams.m_nCheckIfNoneOfNExceedsTresh){
		out.m_PixelOut.m_bNotMoreThenNAboveRejected = TRUE;
		return FALSE;
	}else{
		out.m_PixelOut.m_bNotMoreThenNAboveRejected = FALSE;
		return TRUE;
	}
}

LONGLONG_T CCD_Analyser::CalcMaxSumOnPrevRotCorrect( LONG_T* neighb_list, LONG_T ncnt,
              const CPixelAnalyseIn& in )
{
	LONGLONG_T MaxSum=-100000,newSum=-1;
	LONGLONG_T MaxIdx=-1;

	// loop is starting from OLDEST frame !
	for(register int i=0;i<in.PrevMatrixPtrCnt;i++){		
		//double frame_dx = -(frames_back-i)*gCCDParams.m_FrameDX;
		//double frame_dy = -(frames_back-i)*gCCDParams.m_FrameDY;

		LONGLONG_T Sum = (in.PrevMatrixPtr[i])->CalcSumRotCorrected( neighb_list, 
										ncnt, in.PrevFramesShiftTab[i].frameDX, 
                             in.PrevFramesShiftTab[i].frameDY );
		if(Sum>MaxSum){
			MaxSum = Sum;
		}
	}
	return MaxSum;
}

LONGLONG_T CCD_Analyser::CalcMaxSumOnPrev( LONG_T* neighb_list, LONG_T ncnt, 
													    Table2D<ELEM_TYPE>** PrevMatrixPtr, LONG_T frames_back )
{
	LONGLONG_T MaxSum=-1,newSum=-1;
	LONGLONG_T MaxIdx=-1;

	for(int i=0;i<frames_back;i++){	
		LONGLONG_T Sum = (PrevMatrixPtr[i])->CalcSum( neighb_list, ncnt);
		if(Sum>MaxSum){
			MaxSum = Sum;
		}
	}
	return MaxSum;
}

double CCD_Analyser::CalcGaussMatrix( Table2D<double>& gaussMatrix, double fwhm, long radius, int mode )
{
	int gmsize = 2*radius + 1;
	int gmsize2 = gmsize * gmsize;
	double gmsize2_1 = 1.0 / gmsize2;


	int center = radius;
	double sig = 0.4246*fwhm;
	double dsig2 = 2*sig*sig;
	double dpisig2 = 3.1415926*dsig2; //=2*pi*sig^2
	double sum1 = 0.0;


	  /*if (radius>7) {
	 
	  }*/

	gaussMatrix.Alloc( gmsize, gmsize );
	double* p_data = gaussMatrix.get_data_buffer();
	double** p_data_fast = gaussMatrix.get_data_buffer_fast();


	/* Creation of matrix with Gauss/normal distribution values */
	for (register long y=0; y<gmsize; y++) {
		for (register long x=0; x<gmsize; x++){
			p_data_fast[y][x] = exp( -((x-center)*(x-center) + (y-center)*(y-center)) / dsig2)/dpisig2;
			sum1 += p_data_fast[y][x];
		}
	}



	/* Additional operations to make matrix work better for us */
	double sum0 = 0.0;
	double dif = (1.0 - sum1) / gmsize2 ;

	for (register long y=0; y<gmsize; y++) {
		for (register long x=0; x<gmsize; x++){
			if (mode == 0)
				p_data_fast[y][x] = p_data_fast[y][x]  - gmsize2_1;

			if (mode == 1)
				p_data_fast[y][x] = p_data_fast[y][x]/sum1 - gmsize2_1;

			if (mode == 2)
				p_data_fast[y][x] = p_data_fast[y][x] + dif - gmsize2_1;

			/*
				if (mode == 3)   gauss_matrix[x][y] = gauss_matrix[x][y];
				if (mode == 4)   gauss_matrix[x][y] = gauss_matrix[x][y] / sum1;
				if (mode == 5)   gauss_matrix[x][y] = gauss_matrix[x][y] + dif;
			*/ 
 
			sum0 += p_data_fast[y][x];

		}
	}

	return sum0;		
}



LONGLONG_T CCD_Analyser::CalcSumOpt( LONG_T* neighb_list, LONG_T ncnt, 
				     const ELEM_TYPE* pData )
{
	LONGLONG_T sum=0;
	for(register int i=0;i<ncnt && neighb_list[i]!=NEIGHB_LIST_END;i++){
		sum += pData[neighb_list[i]];
	}	
/*#ifdef _DEBUG
	if(sum<0){
		printf("sum of neighbour pixels <0\n");
		for(i=0;i<ncnt;i++){
			printf("data[%d]=%d\n", neighb_list[i],pData[ neighb_list[i] ]);
		}
		Assert(FALSE,"Negative value of neighbours sum");
	}
#endif*/
	return sum;
}

LONGLONG_T CCD_Analyser::CalcSum( CLongList& List, const ELEM_TYPE* pData,
                                  CPixelList* pUsedList/*=NULL*/ )
{
	LONGLONG_T sum=0;
	CLongList::iterator pPix;	
	for(pPix=List.begin();pPix!=List.end();pPix++){
		sum += pData[ *pPix ];
	}
/*#ifdef _DEBUG
	if(sum<0){
		printf("sum of neighbour pixels <0\n");
		for(pPix=List.begin();pPix!=List.end();pPix++){			
			printf("data[%d]=%d\n", *pPix,pData[ *pPix ]);
		}
		Assert(FALSE,"Negative value of neighbours sum");
	}
#endif*/
	return sum;
}

LONGLONG_T CCD_Analyser::CalcSum( ELEM_TYPE* p_data, long x, long y, 
                                  long xSize, long ySize)
{
	LONGLONG_T sum = 0;
	long pos = x+y*xSize;
	switch(gCCDParams.m_nNeighbToSumCount){
		case 5:
		case 4 :
				/*    x
                 xxx
                  x    */
				{
					LONG_T left = (x>0 ? p_data[pos-1] : 0 );
					LONG_T right = (x<(xSize-1) ? p_data[pos+1] : 0);
					LONG_T up = (y>0 ? p_data[pos-xSize] : 0);	
					LONG_T down = (y<(ySize-1) ? p_data[pos+xSize] : 0);
					LONG_T point = p_data[pos];
					VerifyValue(left);
					VerifyValue(right);
					VerifyValue(up);
					VerifyValue(down);
					VerifyValue(point);
					sum = sum + point + left + right + up + down;
				}
		default :
				{
					LONG_T val = p_data[x+y*xSize];
					VerifyValue(val);
					sum += val;
				}
	}
	return sum;
}

BOOL_T CCD_Analyser::AnalyseSumOfNeighbours(CCDPipeline& ccd_pipeline)
{
	if(ccd_pipeline.GetCount()!=ccd_pipeline.GetPipelineSize())
		return FALSE;
	BOOL_T bRet = FALSE;
	if(gCCDParams.m_bAnalyseSumAround){
		clock_t t1 = clock();	
		Assert(gCCDParams.m_FramesBack<ccd_pipeline.GetCount(),"Analysis can not use %d frames, because pipeline size is %d",gCCDParams.m_FramesBack,ccd_pipeline.GetCount());
		cCCD& newFrame = ccd_pipeline.GetCurrent();
		long size = newFrame.GetCount();
		for(int i=0;i<size;i++){
			CCDMatrix& Matrix = newFrame[i];			
			ELEM_TYPE* p_data = Matrix.get_data_buffer();

			long xSize = Matrix.GetXSize();
			long ySize = Matrix.GetYSize();
			long yUpperLimit = (ySize-gCCDParams.m_nIgnoreEdge);
	      long xUpperLimit = (xSize-gCCDParams.m_nIgnoreEdge);

			CPixelList pixel_list(xSize*ySize);
			CLongList Pixels;
			for(int y=gCCDParams.m_nIgnoreEdge;y<yUpperLimit;y++){
            for(int x=gCCDParams.m_nIgnoreEdge;x<xUpperLimit;x++){
					long pos = x+y*xSize;
					if(!pixel_list.CheckPixel(pos)){
						// calculate sums only for pixels not included in
                  // any previous sum
						GetNeighbPositions(x,y,xSize,ySize,Pixels);
						LONGLONG_T newSum = CalcSum( Pixels, p_data );
						cCCD* pFrame;
						long num=0;
						BOOL_T bEvent=TRUE;
						ccd_pipeline.GetCurr();
						for(pFrame = ccd_pipeline.GetPrev();pFrame!=NULL && num<gCCDParams.m_FramesBack  
                      ;pFrame= ccd_pipeline.GetPrev()){
							CCDMatrix& PrevMatrix = (*pFrame)[i];
							LONGLONG_T PrevSum = CalcSum( Pixels, PrevMatrix.get_data_buffer() );
							if(newSum-PrevSum<=gCCDParams.m_SumTresholdForNewFrame){
								bEvent = FALSE;
								break;
							}
						}
						if(bEvent){
							// potential event found
							CLongList cluster;
							FindSumCluster( Matrix, ccd_pipeline, x, y, cluster, pixel_list );
							long x0,y0;
							CalcCenterOfHit( p_data, cluster, xSize, x0, y0 );
							Assert(x0>=0 && x0<xSize,"X=%d coordinate of hit center out of range [%d-%d]",x0,0,xSize);
							Assert(y0>=0 && y0<ySize,"X=%d coordinate of hit center out of range [%d-%d]",y0,0,ySize);
							Matrix.AddFoundEvent( x0, y0, cluster );
							Matrix.SetInteresting();
							newFrame.SetInteresting();
							MYTRACE2(gCCDTrace,"Sum Algorithm : Potential event found x=" << x0 << ", y="<<y0);
						}
					}
				}
			}
		}
	}
	return bRet;
}

Table2D<ELEM_TYPE>** CCD_Analyser::GetPrevMatrixPtrs( LONG_T index, CCDPipeline& ccd_pipeline )
{
	Table2D<ELEM_TYPE>** pMatrixPtrTab = new Table2D<ELEM_TYPE>*[ccd_pipeline.GetPipelineSize()];
	CCDPipelineIterator i(&ccd_pipeline);
	
	int cnt=0;
	for(;!i.curr();i++){
		pMatrixPtrTab[cnt] = i->GetMatrix( index );
		cnt++;
	}
	return pMatrixPtrTab;
}


//
// only analysis of previous pipeline frames 
//
BOOL_T CCD_Analyser::FullNewFrameAnalysePrev(CCDPipeline& ccd_pipeline)
{
	Initialize();
	LONG_T pipeline_size = ccd_pipeline.GetPipelineSize();
	LONG_T pipeline_size_minus_1 = pipeline_size-1;
	if(ccd_pipeline.GetCount()!=pipeline_size)
		return FALSE;

	BOOL_T bRet = FALSE;
	LONG_T ncnt = gCCDParams.m_nNeighbToSumCount;
	LONG_T NeighbList[MAX_CLUSTER_SIZE];
	LONG_T* neighb_list=NeighbList;

	if(gCCDParams.m_bAnalyseMaxFromPrevInPipeline){
		PROFILER_START
		cCCD& newFrame = ccd_pipeline.GetCurrent();
		long size = newFrame.GetCount();
		if(size<1){
			return FALSE;
		}
		long xSize = newFrame[0].GetXSize();
      long ySize = newFrame[0].GetYSize();
		long yUpperLimit = (ySize-gCCDParams.m_nIgnoreEdge);
		long xUpperLimit = (xSize-gCCDParams.m_nIgnoreEdge);		

		for(int i=0;i<size;i++){
			CCDMatrix& Matrix = newFrame[i];			
			ELEM_TYPE* p_data = Matrix.get_data_buffer();
			Table2D<ELEM_TYPE>** PrevMatrixPtr = GetPrevMatrixPtrs( i, ccd_pipeline );;
			CPixelList pixel_list(xSize*ySize);
			
			long y_pos = (gCCDParams.m_nIgnoreEdge-1)*xSize+gCCDParams.m_nIgnoreEdge;
			for(int y=gCCDParams.m_nIgnoreEdge;y<yUpperLimit;y++){
				y_pos = y_pos + xSize;
				long pos = y_pos;
				for(int x=gCCDParams.m_nIgnoreEdge;x<xUpperLimit;x++,pos++){
					if(!pixel_list.CheckPixel(pos)){						
						// calculate sums only for pixels not included in
                  // any previously found cluster

						// to optimized - on preprocesor level - OPTPREPROC
						if(m_pNeighMap){
							neighb_list = m_pNeighMap->GetDesc( pos, ncnt );
						}else{
							GetAnalNeighbWithSelf(x,y,pos,xSize,ySize,neighb_list,ncnt);
						}

						LONGLONG_T newSum = CalcSumOpt( neighb_list, ncnt, p_data );
						LONGLONG_T maxSum = CalcMaxSumOnPrev( neighb_list, ncnt, 
                                                        PrevMatrixPtr, pipeline_size_minus_1 );

						long x_y_val = p_data[pos];
//#ifdef _DO_NOT_ALLOW_NEGATIVE_
//              Assert(x_y_val>=0,"Negative value found in matrix analysis = %d at position (%d,%d)",x_y_val,x,y);			
//#endif

						if(x_y_val>=gCCDParams.m_MaxAllowedVal){
							// skiping pixel exceeding maximum allowed values in
							continue;
						}
						if( CheckCondition(newSum,maxSum) ){
							CLongList cluster,star_cluster;

							LONG_T x0,y0;
							double r_max;
							FindClusterAboveTreshold( Matrix, x, y, pos, x0, 
                                               y0, r_max, pixel_list, cluster );
	
							if(cluster.size()<gCCDParams.m_MinClusterSize){
							 	MYTRACE4(gCCDTrace,"Potential event at x=" << x << ",y=" << y << ", but cluster size=" << cluster.size() << "<" << gCCDParams.m_MinClusterSize);
							 	continue;	
							}
							LONG_T pos0 = x0 + y0*xSize;
							CLongList ClusterWithMore;
							GetClusterWithPointsAround(cluster,ClusterWithMore,
																xSize,ySize,gCCDParams.m_nPixelsAroundToConfirm);
							if(!ConfirmEvent_InCluster( i , ccd_pipeline, Matrix, 
                                                   pos0, xSize, ySize, ClusterWithMore )){
								MYTRACE4(gCCDTrace,"Potential event at x=" << x0 << ",y=" << y0 << ", but rejected by prev sum rejection procedure !");
								continue;									
							}
							bRet = TRUE;

							Matrix.AddFoundEvent( x0, y0, cluster );
							Matrix.SetInteresting();
							newFrame.SetInteresting();
							MYTRACE2(gCCDTrace,"ACCEPTED : Potential event found x=" << x0 << ", y="<<y0<<" ccd_index="<<i<<", sum="<<newSum<<", max_sum="<<maxSum);
						}	
					}
				}
			}

			if(gCCDParams.m_bSkipIfMoreThen>0 && gCCDParams.m_eCheckIfMorePoint!=eAfterCoic){
	         // requires rejection of all events in case in the redial of gCCDParams.m
   	      // there is more then m_bSkipIfMoreThen events :
      	   Matrix.RejectIfMoreThen( &ccd_pipeline );
	      }
			


			delete [] PrevMatrixPtr;
		}
		PROFILER_END("Full analysis of new frame took : ")
	}
	return bRet;
}


BOOL_T CCD_Analyser::CheckSumEventsToNormalTracks( CCDPipeline* pPipeline )
{
	my_printf_now("###### AVG_ANAL : checking against normal track #######\n");
	printf_now2("Cam%d :\n",pPipeline->GetPipelineIndex());
	// NEW - 20041021 , no velocity check is performed for sumed frames events
	// due to fackt that it often fails
	BOOL_T bDoCheckVelocity = TRUE;
	BOOL_T bRet = CheckNormalTracks( pPipeline, pPipeline->m_EventsOnSumedFrame, 
												bDoCheckVelocity );
	if( gCCDParams.m_bCheckForSUPERNEW  ){
		BOOL_T bRet2 = CheckNormalTracks( pPipeline, pPipeline->m_BrightenOnSumList, 
													bDoCheckVelocity );

	}

	printf_now2("######### AVG_ANAL : track check (cam%d) end #######\n",pPipeline->GetPipelineIndex());
	return bRet;
}

BOOL_T CCD_Analyser::CheckNormalTracks( CCDPipeline* pPipeline, CCDEventList& events,
													 BOOL_T bDoCheckVelocity ){
	CCDEventList::iterator i;
	int frameno = pPipeline->GetDayFrameCounter();

	for(i=events.begin();i!=events.end();i++){
		if( i->IsIdentified() ){
			CTrackDesc track;
			// BOOL_T bDoCheckVelocity = gCCDParams.m_bCheckVelocity;
			// now no check is performed 
			// BOOL_T bDoCheckVelocity = FALSE; // NEW - 20041021
			BOOL_T bRet = CheckIfEventBelongsToTrack( (*i), frameno, pPipeline,
																	bDoCheckVelocity,
						                                 gCCDParams.m_fVelocityError, track );
			if( gDebugTracks ){
				if( bRet ){
					printf("SumFramesEvent (%d,%d) rejected by normal track %d\n",
							(int)(i->m_MaxPoint).x,(int)(i->m_MaxPoint).y, track.track_id);
				}
			}
		}
	}
}

BOOL_T CCD_Analyser::CheckIfEventBelongsToTrack( CccdReport& event, int confirmedFrameIndex,
																 CCDPipeline* pPipeline,
									                      BOOL_T bCheckVelocity, double fVelocityError,
																 CTrackDesc& track )
{
	CTrackList& oldTracks = pPipeline->GetTrackList();
	CTrackList::iterator it;
	/*for(it=oldTracks.begin();it!=oldTracks.end();it++){
		BOOL_T bOK=TRUE;
		int index = (it->m_EventsOnTrack).back().m_FrameIndex;
		if( (confirmedFrameIndex-index)>gCCDParams.m_nCheckFramesIfNotOlderThenNFrames )
		{
			bOK=FALSE;
		}
		printf("%d(%d,%d,%d),",it->track_id,(int)bOK,confirmedFrameIndex,index);
	}*/
	printf("\n");

	CTrackList::iterator t;
	for(t=oldTracks.begin();t!=oldTracks.end();t++){
		// if( (confirmedFrameIndex-(t->m_EventsOnTrack)[0].m_FrameIndex)>gCCDParams.m_nCheckFramesIfNotOlderThenNFrames )
		// [NEW] 20040927 : now check if last event 
		// was not added longer then :
		if( (confirmedFrameIndex-(t->m_EventsOnTrack).back().m_FrameIndex)>gCCDParams.m_nCheckFramesIfNotOlderThenNFrames )
		{
			// do not re-fit to tracks older then 200 frames :		
			continue;
		}

		BOOL_T bReFit=FALSE;	
		double new_a,new_b;

		if(CheckIfEventBelongsToTrack( event, (*t) , bReFit, 
												new_a, new_b, bCheckVelocity,
												fVelocityError )){
			event.m_PixelAnalResults.m_bRejectedDueToTrack = TRUE;
			event.m_PixelAnalResults.m_TrackType = eNormalTrack;
			track = (*t);
			return TRUE;
		}						
	}
	return FALSE;
}


BOOL_T CCD_Analyser::CheckIfEventBelongsToTrack( CccdReport& event, CTrackDesc& track, BOOL_T& bReFitted,
																 double& new_a, double& new_b,
																 BOOL_T bCheckVelocity, double velError,
																 BOOL_T bDoFit /* = TRUE */  )
{
	int cam_idx = m_pPipeline->GetPipelineIndex();

	double minDist = track.FindMinDist( event.m_Point.x, event.m_Point.y );
	/* [NEW] - 20040927 
		when adding new points to exisitng track do not check distance 
		if( minDist<=gCCDParams.m_nMinDistOfEventsInTrack ){
		if( gDebugTracks ){
         track.LogEventCheck( event, cam_idx, 0, 0, 0, 0 ,0 , FALSE, minDist, "min_dist" );
		}
		return FALSE;
	}*/
	event.m_PixelAnalResults.minDist = minDist;

	bReFitted=FALSE;
	double chi2 = CMyFit::CalcDist2FromLine2Par( track.a, track.b, event.m_Point.x, event.m_Point.y );

	if( m_pAnalFoundEventFunc ){
		if( m_eCurrFitTrackType == eNormalTrack ){
			(*m_pAnalFoundEventFunc)( &chi2, eChi2ToOld, cam_idx, NULL );
		}else{
			if( m_eCurrFitTrackType == eTrackOnSumedFrame ){
				(*m_pAnalFoundEventFunc)( &chi2, eChi2ToOldSum, cam_idx, NULL );
			}
		}
	}

	if( chi2 < event.m_PixelAnalResults.minChi2 ){
		event.m_PixelAnalResults.minChi2 = chi2;
		event.m_PixelAnalResults.bestTrackID = track.track_id;
	}

	if(chi2<gCCDParams.m_MaxChi2ForPointToMatchLine){
		double* x_values = new double[track.m_EventsOnTrack.size()+1];
	   double* y_values = new double[track.m_EventsOnTrack.size()+1];

   	int cnt = track.get_xy_values( x_values, y_values );
	   x_values[cnt] = event.m_Point.x;
   	y_values[cnt] = event.m_Point.y;
	   cnt++;

		double a=track.a,b=track.b;
	   double new_line_chi2 = 0;
		if( bDoFit ){
			new_line_chi2 = CMyFit::FitLineChi2( x_values, y_values, cnt, a, b);
		}

		if( gDebugTracks ){
			_TRACE_PRINTF_2( "Event %d-(%d,%d) chi2=%.2f OK ( track=%d ) , ",event.m_DayFrameIndex,
						(int)event.m_MaxPoint.x,(int)event.m_MaxPoint.y,
						chi2, track.track_id );
		}else{
			_TRACE_PRINTF_2("Track debug is off\n");
		}

		BOOL_T bOK=TRUE;
		double vx=0,vy=0;
		if( bCheckVelocity ){
			time_t dt = ( event.m_Time - track.m_EventsOnTrack.back().m_Time );
			vx = ( event.m_MaxPoint.x - track.m_EventsOnTrack.back().m_MaxPoint.x )/dt;
			vy = ( event.m_MaxPoint.y - track.m_EventsOnTrack.back().m_MaxPoint.y )/dt;


			if( dt==0 ){
				// skiping if any event from same frame already added :
				bOK = FALSE;
				printf("Event from same frame already added %d-(%d,%d) rejected from track\n",event.m_DayFrameIndex,(int)event.m_MaxPoint.x,(int)event.m_MaxPoint.y);
			}else{
				if( !CMyFit::CheckVelocityCondition( track.vx, track.vy,
																  vx, vy, velError,
																  event.m_PixelAnalResults.rx, 
																  event.m_PixelAnalResults.ry ,
																  DEFAULT_CCD_IDX_FOR_CHECK_COND,
																  m_eCurrFitTrackType ) ){
					event.m_PixelAnalResults.veloCheckOK=0;
					bOK = FALSE;

					_TRACE_PRINTF_2( "Event %d-(%d,%d) velocity check failed for track=%d",event.m_DayFrameIndex,
						(int)event.m_MaxPoint.x,(int)event.m_MaxPoint.y,
						track.track_id );					
				}else{
					event.m_PixelAnalResults.veloCheckOK=1;
					if(gDebugTracks){
						printf("Event %d-(%d,%d) added to track %d, (vx,vy)=(%.2f,%.2f) (rx,ry)=(%.2f,%.2f) err=%.2f\n",
							event.m_DayFrameIndex,
							(int)event.m_MaxPoint.x,(int)event.m_MaxPoint.y,track.track_id,
							vx,vy,event.m_PixelAnalResults.rx,event.m_PixelAnalResults.ry,velError);
					}
				}
				if( m_pAnalFoundEventFunc ){
					FillVelocityHisto( track.vx, track.vy, vx, vy, m_pPipeline->GetPipelineIndex() );
				}
			}
		}else{
			if( gDebugTracks ){
				printf("velocity check skiped\n");
			}
		}
		/*if( bOK && new_line_chi2>gCCDParams.m_MaxChi2InTrack){
			bOK = FALSE;
		}*/

		if( bOK ){
			new_a = a;
   	   new_b = b;
      	bReFitted = TRUE;
		}
	
		if( gDebugTracks ){
			mystring szReason="OK";
			if(!bOK){
				szReason = "velocity";
			}
			track.LogEventCheck( event, cam_idx, chi2, vx,vy,
										event.m_PixelAnalResults.rx,
										event.m_PixelAnalResults.ry, bOK,
										minDist, szReason.c_str() );
			if( bCheckVelocity ){
				if( gDebugTracks>1 ){
					printf("velocity check %s : (rx,ry)=(%.2f,%.2f)\n",
								GetOK( bOK ),
								event.m_PixelAnalResults.rx,
								event.m_PixelAnalResults.ry);
				}
			}
		}

      delete [] x_values;
      delete [] y_values;
		fflush(0); // printf
		return bOK;
	}else{
		if( gDebugTracks ){
			track.LogEventCheck( event, cam_idx, chi2, 0, 0, 0,0,FALSE, minDist, "chi2");
		}		
	}

	fflush(0); // printf
	return FALSE;
}

void CCD_Analyser::FillVelocityHisto( double vx, double vy,
                                      double vx_new, double vy_new, 
												  int ccd_idx )
{
	if( m_pAnalFoundEventFunc ){
		double vx_min=MIN(fabs(vx),fabs(vx_new));
	   double vy_min=MIN(fabs(vy),fabs(vy_new));
   	double vx_max=MAX(fabs(vx),fabs(vx_new));
	   double vy_max=MAX(fabs(vy),fabs(vy_new));

		double r = (vx_min/vx_max);
		(*m_pAnalFoundEventFunc)( &r, eHistoVXRatioToOld, ccd_idx, NULL );

		r = (vy_min/vy_max);
		(*m_pAnalFoundEventFunc)( &r, eHistoVYRatioToOld, ccd_idx, NULL );
	}
}

void CCD_Analyser::UpdateBestFit( const CPixelAnalyseIn& in, int startFrame, int EndFrame, 
											 CTrackDesc& old_track, CTrackDesc& new_track, 
											 const CccdReport& newEvent, CTrackList& tracks,
											 BOOL_T bAddNewEventOnly )
{
	for(int back1=startFrame;back1<EndFrame;back1++){
		CTrackList::iterator t1;
 		for(t1=tracks.begin();t1!=tracks.end();t1++){
			if( (*t1) == old_track ){
				(t1->m_EventsOnTrack).push_back( (CEventBaseInfo&)newEvent );
				if(!bAddNewEventOnly){
					t1->a = new_track.a;
					t1->b = new_track.b;
					t1->m_cam_idx = new_track.m_cam_idx;
				}
				t1->calcAverageVelocity();
			}
		}			
	}					
		
}

BOOL_T CCD_Analyser::CheckIfInternalTrigger( CccdReport& newEvent )
{
	// TODO :
	// no criteria was desinged yet ...
	if(!IsVerifiedOK( newEvent ))
		return FALSE;

	return FALSE;
}

BOOL_T CCD_Analyser::CheckIfInternalTrigger( CccdReport& coicEvent1, CccdReport& coicEvent2 )
{
	// TODO :
	// create some criteria for coicyding events to claim triggers :
	if(!IsVerifiedOK( coicEvent1 ) || !IsVerifiedOK( coicEvent2 ))
		return FALSE;

	return FALSE;
}

BOOL_T CCD_Analyser::VerifyIfEventsNotOnTrack( const CPixelAnalyseIn& in, CCDMatrix& Matrix )
{
	int cam_idx = (in.pPipeline)->GetPipelineIndex();
	m_eCurrFitTrackType =  eNormalTrack;

	printf("########### CAM%d TRACKS #############\n",cam_idx);
	BOOL_T bRet = VerifyIfEventsNotOnTrackBase( in, Matrix, (in.pPipeline)->GetTrackList(),
														  gCCDParams.m_MaxEventRate,
														  gCCDParams.m_nNumBackFramesForTracks,
														  gCCDParams.m_bCheckVelocity,
														  gCCDParams.m_fVelocityError,
														  gCCDParams.m_nMaxEventsToFitTrack,
														  gCCDParams.m_nNumBackFramesHasEventsForTracks,
														  gCCDParams.m_MaxChi2InTrack,
														  gCCDParams.m_nMinEventNoInTrack,
														  gCCDParams.m_MaxChi2ForPointToMatchLine,
														  gCCDParams.m_nCheckFramesIfNotOlderThenNFrames );
	printf("####### END OF CAM%d TRACKS #########\n",cam_idx);

	return bRet;	
}

void CCD_Analyser::AddNewNormalTracks( CCDPipeline* pPipeline )
{
	CTrackList& normalTracks = pPipeline->GetTrackList();
	CTrackList& planeTracks  = pPipeline->m_PlaneTracks;
	
	CTrackList::iterator i;
	int frameno = pPipeline->GetFrameIndex();

	my_printf_now("Updating plane tracks by new normal tracks :");
	for(i = normalTracks.begin();i!=normalTracks.end();i++){
		if( i->frame_index == frameno ){
			planeTracks.push_back( *i );
			printf("%d,",i->track_id);
		}
	}
	printf("\n");
}

BOOL_T CCD_Analyser::VerifyPlaneTracks( const CPixelAnalyseIn& in, CCDMatrix& Matrix,
													 BOOL_T bDoCheckForNewTracks,
													 BOOL_T bDoCheckVelocity,
			                               double fVelError )
{
	// tu by trzeba dodac najpierw dodanie nowych trackow zwyklych z
	// ostatniej klatki do listy m_PlaneTracks
	// tak ze proboujemy tez sie dopasowac do normal trackow 
	// o to moze byc 1 punkt wyjscia - tylko dopasowac sie 
	// do istniejacych normal trackow - nie szukac jeszcze nowych !

	BOOL_T bRet=FALSE;
	if( gCCDParams.m_bCheckPlaneTracks ){
		// AddNewNormalTracks( (in.pPipeline) );		

		// 
		// int numBack = gCCDParams.m_nCheckFramesIfNotOlderThenNFrames;
		int numBack = gCCDParams.m_nNumBackFramesForPlaneTrack; // only 5 frames back when checking plane tracks

		/*BOOL_T bDoCheckVelocity=FALSE;
		double fVelError=0.00;*/

		int nNumBackFramesHasEventsForTracks = 0;
		bRet = VerifyIfEventsNotOnTrackBase( in, Matrix, (in.pPipeline)->m_PlaneTracks, 		
														  1000, gCCDParams.m_nNumBackFramesForTracks,
														  bDoCheckVelocity, fVelError, 
															1000, nNumBackFramesHasEventsForTracks,
														  gCCDParams.m_MaxChi2InPlaneTrack,
														  gCCDParams.m_nMinPointNoOnPlaneTrack,
														  gCCDParams.m_MaxChi2InPlaneTrack,
														  numBack,
														  ePlaneTrack, 
														  bDoCheckForNewTracks,
														  gCCDParams.m_MinFramesOnPlaneTrack );
	}
	return bRet;
}


BOOL_T CCD_Analyser::VerifyIfEventsNotOnTrackBase( const CPixelAnalyseIn& in, CCDMatrix& Matrix, 
																CTrackList& oldTracks,
																int MaxEventRate, int nNumBackFramesForTracks,
																BOOL_T bCheckVelocity, double fVelocityError, 
																int nMaxEventsToFitTrack, 
                                                int nNumBackFramesHasEventsForTracks,
																double MaxChi2InTrack, int nMinEventNoInTrack,
																double MaxChi2ForPointToMatchLine,
																int nCheckFramesIfNotOlderThenNFrames,
																eTrackCheckType_T track_type,
																BOOL_T bDoCheckForNewTracks /*=TRUE*/,
																int nMinFramesWithEvent	)
{
	if( (in.pPipeline->m_allFoundEvents).size()==0 ){
		return FALSE;			
	}

	m_eCurrFitTrackType =  track_type;

	if( gDebugTracks )printf("1 in CCD_Analyser::VerifyIfEventsNotOnTrackBase\n");

	// in.pPipeline->VerifyIfEventsNotOnTrack( in, Matrix );	
	deque<CFrameEvents>& allPipelineEvents = in.pPipeline->m_allFoundEvents;
	int nConfirmOnNext=gCCDParams.m_ConfirmEventsOnNextNFrames;
	int nNewFrameIndex=(allPipelineEvents.size()-1);
	int confirmedFrameIndex = ((in.pPipeline)->GetFrameIndex()-gCCDParams.m_ConfirmEventsOnNextNFrames);	
	BOOL_T bRet = FALSE;
	BOOL_T bHasTracks=TRUE;
	CCDEventList* newEvents = &(allPipelineEvents.back()[in.ccd_index]);
	CImageCreator* pGenObj = (in.pPipeline)->GetGenObj();
	CFrameIdentStat* pBackgrStat=NULL;
	CFrameIdentStat* pSampleStat=NULL;
	if( gCCDParams.m_bLogFrameStat && backgrStat.size() ){
		pBackgrStat = &(backgrStat.back());
	}
	if( gCCDParams.m_bLogFrameStat && samplesStat.size() ){
		pSampleStat = &(samplesStat.back());
	}

	// checking if more frames then required :
	int checkBack = nNumBackFramesForTracks-1;

	if( nConfirmOnNext>0 ){
		nNewFrameIndex = (allPipelineEvents.size()-1)-nConfirmOnNext;
		if( nNewFrameIndex<0 )
			return FALSE;
		newEvents = &(allPipelineEvents[nNewFrameIndex][in.ccd_index]);
	}
	if(newEvents->size()>=MaxEventRate){
		printf("Number of events %d to big to fit track on cam%d\n",newEvents->size(),(in.pPipeline)->GetPipelineIndex());fflush(0);
		return FALSE;
	}
	
//	if(allPipelineEvents.size()>=nNumBackFramesForTracks){
		// in case already analaysied at least checkBack back frames 

		int allEventsCount=0;
		int allFramesWithEvents=0;


		int firstFrameToCheck=allPipelineEvents.size()-nNumBackFramesForTracks;
		if( firstFrameToCheck<0 ){
			firstFrameToCheck=0;
		}
		int startFrom=(allPipelineEvents.size()-1-nConfirmOnNext);
		if( startFrom<0 ){
			return FALSE;
		}

		if( gDebugTracks )printf("2 in CCD_Analyser::VerifyIfEventsNotOnTrackBase\n");


		if(oldTracks.size()){
			// there are already found tracks - verify if new events belong to one of them :
			if( gDebugTracks>1 ){
				printf("Tracks to be checked : ");
				CTrackList::iterator it;
				for(it=oldTracks.begin();it!=oldTracks.end();it++){
					BOOL_T bOK=TRUE;
					int index = (it->m_EventsOnTrack).back().m_FrameIndex;
					if( (confirmedFrameIndex-index)>nCheckFramesIfNotOlderThenNFrames )
					{
						bOK=FALSE;
					}
					printf("%d(%d,%d,%d),",it->track_id,(int)bOK,confirmedFrameIndex,index);
				}
				printf("\n");
			}

			for(int k=0;k<newEvents->size();k++){
				CccdReport& evt = (*newEvents)[k];
				if( evt.m_PixelAnalResults.m_bRejectedByNextFrames || 
					 evt.m_PixelAnalResults.m_bRejectedDueToTrack )
					continue;
				if( track_type == ePlaneTrack && !evt.IsIdentified() ){
					// if plane tracks check , skip already rejected by other 
					// criterias :
					continue;
				}
				
				CTrackList::iterator t;
				for(t=oldTracks.begin();t!=oldTracks.end();t++){
					// if( (confirmedFrameIndex-(t->m_EventsOnTrack)[0].m_FrameIndex)>nCheckFramesIfNotOlderThenNFrames )
					// [NEW] 20040927 : now check if last event 
					// was not added longer then :
					if( (confirmedFrameIndex-(t->m_EventsOnTrack).back().m_FrameIndex)>nCheckFramesIfNotOlderThenNFrames )
					{
						// do not re-fit to tracks older then 200 frames :		
						continue;
					}

					BOOL_T bReFit=FALSE;	
					double new_a,new_b;

					if(CheckIfEventBelongsToTrack( (*newEvents)[k], (*t) , bReFit, 
										new_a, new_b, bCheckVelocity,
										fVelocityError )){
						CTrackDesc new_track( (*t) );
						new_track.a = new_a;
						new_track.b = new_b;

						(*newEvents)[k].m_PixelAnalResults.m_bRejectedDueToTrack = TRUE;
						(*newEvents)[k].m_PixelAnalResults.m_TrackType = track_type;
						
						// count to statistics :
						if( !gCCDParams.GetMC() || !pGenObj || !pGenObj->IsGenEvent( (*newEvents)[k].m_MaxPoint ) ){
							if( pBackgrStat ){
								pBackgrStat->DecTrackAndCoicEvents();
							}
						}else{
							if( gCCDParams.GetMC() && pGenObj && pGenObj->IsGenEvent( (*newEvents)[k].m_MaxPoint ) ){
								if(pSampleStat){
									pSampleStat->DecTrackAndCoicEvents();
								}
							}
						}

						BOOL_T bFound=TRUE;
						if(nConfirmOnNext==0){
							bFound = Matrix.GetFoundEvents().RejectEventDueToTrack( (int)((*newEvents)[k].m_Point.x) , (int)((*newEvents)[k].m_Point.y) );
						}else{
							allPipelineEvents[startFrom][in.ccd_index].RejectEventDueToTrack( (int)((*newEvents)[k].m_Point.x) , (int)((*newEvents)[k].m_Point.y) );
						}
						Assert(bFound,"Event (%d,%d) on frame %d not found on Matrix.GetFoundEvents list",(int)((*newEvents)[k].m_Point.x) , (int)((*newEvents)[k].m_Point.y), (in.pPipeline)->GetFrameIndex());

						(new_track.m_EventsOnTrack).push_back( ((*newEvents)[k]) );
						new_track.calcAverageVelocity();
						if( !oldTracks.Find( (*t) ) ){
							// adding only first track or it is not update of already present track :
							oldTracks.Add( new_track );
						}
									
						if(bReFit){
							// upadting after refiting 
							(t->m_EventsOnTrack).push_back( ((*newEvents)[k]) );
							t->a = new_track.a;
							t->b = new_track.b;
							t->m_cam_idx = new_track.m_cam_idx;
						}else{
							(t->m_EventsOnTrack).push_back( ((*newEvents)[k]) );
						}
						t->calcAverageVelocity();
						LogNewTrack( new_track, FALSE, track_type );
						LogNewTrack_RADEC( new_track, FALSE, track_type );

						// flaging events as saved to DB 
						t->FlagEventsSaved();

						break;						
					}
				}
			}
		}

		if( !bDoCheckForNewTracks ){
			return bRet;
		}

		
		for(int back=startFrom;back>=firstFrameToCheck;back--){
			// int notMatched = allPipelineEvents[back].GetNotMatchedEventsCount( in.ccd_index );

			// skip generated events ( if we put events they are to many for the procedure
			int notMatched = allPipelineEvents[back].GetNotMatchedEventsCountSkipGen( in.ccd_index );

			allEventsCount += notMatched;
			allFramesWithEvents += (notMatched>0);
		}

  	   if( gDebugTracks )printf("3 in CCD_Analyser::VerifyIfEventsNotOnTrackBase\n");

		if( allEventsCount>0 && allFramesWithEvents>nMinFramesWithEvent && bHasTracks ){
			int nLoops=0;
			while(allEventsCount>0 && allFramesWithEvents>nMinFramesWithEvent && bHasTracks && nLoops<10){
				// trying to find all tracks :
				bHasTracks = FALSE;
				if(allEventsCount<nMaxEventsToFitTrack){
					// not even trying if to much events - maybe something very bad ...

					if(allFramesWithEvents>=nNumBackFramesHasEventsForTracks){
						// no verify if on track - fit line etc ...
						vector<CccdReport*> eventsPtrList;
				
						for(int back=firstFrameToCheck;(back<=startFrom && back<allPipelineEvents.size());back++){
							CCDEventList::iterator pEvt;
							for(pEvt=allPipelineEvents[back][in.ccd_index].begin();pEvt!=allPipelineEvents[back][in.ccd_index].end();pEvt++){
							// if( pEvt->IsIdentifiedForTrack() && 
								if( ( pEvt->IsIdentified() || 
										( gCCDParams.m_bCheckPlaneTracks && pEvt->IsIdentifiedForTrack() ) ) && 
									!pEvt->m_bGenerated && 
									( CccdReport::FindMinDist( eventsPtrList, *pEvt )>gCCDParams.m_nMinDistOfEventsInTrack ) ){
									eventsPtrList.push_back( (CccdReport*)(&(*pEvt)) );
								}
							}
						}
	
						if( eventsPtrList.size()>0 ){
							CTrackDesc track;
							bHasTracks = FlagEventsOnTrack( eventsPtrList, track, in.pPipeline,
																	  bCheckVelocity, fVelocityError,
								                             MaxChi2InTrack, nMinEventNoInTrack,
								                             MaxChi2ForPointToMatchLine, track_type );
							if(bHasTracks){
								BOOL_T bPlaneTrack=FALSE;
								if( gCCDParams.m_bCheckPlaneTracks ){
									// now check if this is not plane track !:
									bPlaneTrack = CheckIfPlaneTrack( eventsPtrList, track, in.pPipeline, 3 );
								}

								bRet = TRUE;
								track.m_cam_idx = in.ccd_index;
								if( bPlaneTrack ){
									if( gCCDParams.m_bLogTracks ){
										LogNewTrack( track, TRUE, ePlaneTrack );
										LogNewTrack_RADEC( track, TRUE, ePlaneTrack  );
									}
									(in.pPipeline)->m_PlaneTracks.push_back( track );	
								}else{
									if( gCCDParams.m_bLogTracks ){
         							LogNewTrack( track, TRUE, track_type );
										LogNewTrack_RADEC( track, TRUE, track_type );
         						}
									oldTracks.push_back( track );
								}
							}
						}
					}else{
						if( gDebugTracks ){
							printf("Number of frames with events %d smaller then minumum of %d required, no track fit\n",allFramesWithEvents,nNumBackFramesHasEventsForTracks);
						}
					}
				}else{
					if( gDebugTracks ){
						printf("Number of events %d > then max rate allowed (%d), no track fit\n",allEventsCount,nMaxEventsToFitTrack);
					}
				}
				nLoops++;

				if(bHasTracks){
					// calculating not matched tracks :	
					allEventsCount = 0;
					allFramesWithEvents = 0;
					for(int back=allPipelineEvents.size()-1;back>=firstFrameToCheck;back--){
						// skip generated - to much mess ...
						int notMatched = allPipelineEvents[back].GetNotMatchedEventsCountSkipGen( in.ccd_index );
						allEventsCount += notMatched;
						allFramesWithEvents += (notMatched>0);
					}
				}
			}		

			if(nLoops>=10){
				printf("Too many loops in track searching ...\n");
			}
		}else{
			if( gDebugTracks ){
				printf("To small number of events to check for tracks (events# = %d, frames with events = %d)\n",allEventsCount,allFramesWithEvents);
			}
		}
//	}	

	return bRet;	
}


BOOL_T CCD_Analyser::VerifySingleCamTracks( const CPixelAnalyseIn& in, CCDMatrix& Matrix,
                                              CCDEventList& newEvents, CCDEventList& oldEvents )
{
	return VerifySingleCamTracks( in, Matrix, (in.pPipeline)->m_SingleCamTracks, 
											newEvents, oldEvents, 
											gCCDParams.m_MaxEventRate, 1, TRUE, 0.5, 
                                 200, 1, 2, 5, 2 );											
}

BOOL_T CCD_Analyser::VerifySumFrameTracks( CCDPipeline* pPipeline, CCDMatrix& Matrix )
{
	m_eCurrFitTrackType = eTrackOnSumedFrame;

	CPixelAnalyseIn in;
	in.pPipeline = pPipeline;
	in.ccd_index = 0;
	CCDEventList oldEvents;
	char szLog[512];
   sprintf(szLog,AVER_FRAME_EVENTS,pPipeline->GetPipelineIndex());
	mystring szPath;
	szPath << gCCDParams.GetOutputDir() << "/" << szLog;

	//(pPipeline->m_OldEventsOnSumedFrame).GetEventsSinceFrame( oldEvents, ( pPipeline->GetDayFrameCounter()-gCCDParams.m_nNumBackFramesForTracksOnSum ) );
	//printf("Found %d events newer then frame %d\n",oldEvents.size(),(pPipeline->GetDayFrameCounter()-gCCDParams.m_nNumBackFramesForTracksOnSum));


	// printf("CCD_Analyser::VerifySumFrameTrack reading old events from file : %s ...",szPath.c_str());
	// oldEvents.ReadEvents( szPath.c_str(), FALSE, pPipeline->GetDayFrameCounter()-gCCDParams.m_nNumBackFramesForTracksOnSum );
	// printf("read %d old events ( later then frame : %d )\n",oldEvents.size(),(pPipeline->GetDayFrameCounter()-gCCDParams.m_nNumBackFramesForTracksOnSum));
	

	BOOL_T bRet = VerifySingleCamTracks( in, Matrix, pPipeline->m_TracksOnSumedFrames,
													 pPipeline->m_EventsOnSumedFrame,
													 pPipeline->m_OldEventsOnSumedFrame, 
													 gCCDParams.m_MaxEventRate, 0,
													 FALSE, 0.5, 
													 gCCDParams.m_MaxEventRate,
													 1, gCCDParams.m_MaxChi2InPlaneTrack,
													 3, gCCDParams.m_MaxChi2InPlaneTrackToOld,
													 eTrackOnSumedFrame );
	return bRet;													 
}

BOOL_T CCD_Analyser::VerifySingleCamTracks( const CPixelAnalyseIn& in, CCDMatrix& Matrix, 
																CTrackList& oldTracks,
															 CCDEventList& newEvents, CCDEventList& oldEvents,
																int MaxEventRate, int nNumBackFramesForTracks,
																BOOL_T bCheckVelocity, double fVelocityError, 
																int nMaxEventsToFitTrack, 
                                                int nNumBackFramesHasEventsForTracks,
																double MaxChi2InTrack, int nMinEventNoInTrack,
																double MaxChi2ForPointToMatchLine,
																eTrackCheckType_T track_type )
{
	// in.pPipeline->VerifyIfEventsNotOnTrack( in, Matrix );	
	// eTrackCheckType_T track_type=eSingleCamTrack;

	BOOL_T bRet = FALSE;
	BOOL_T bHasTracks=TRUE;
	int currFrame = (in.pPipeline)->GetFrameIndex();

	if(newEvents.size()>=MaxEventRate){
		printf("Number of events %d to big to fit track on cam%d\n",newEvents.size(),(in.pPipeline)->GetPipelineIndex());fflush(0);
		return FALSE;
	}
	


	if(oldTracks.size()){
		// there are already found tracks - verify if new events belong to one of them :


		for(int k=0;k<newEvents.size();k++){
			CTrackList::iterator t;
			for(t=oldTracks.begin();t!=oldTracks.end();t++){
				if( (currFrame-(t->m_EventsOnTrack)[0].m_FrameIndex)>gCCDParams.m_nCheckFramesIfNotOlderThenNFrames )
				{
					// do not re-fit to tracks older then 200 frames :
					continue;
				}

				BOOL_T bReFit=FALSE;	
				double new_a,new_b;

				if(CheckIfEventBelongsToTrack( (newEvents)[k], (*t) , bReFit, 
									new_a, new_b, bCheckVelocity,
									fVelocityError )){
					CTrackDesc new_track( (*t) );
					new_track.a = new_a;
					new_track.b = new_b;

					(newEvents)[k].m_PixelAnalResults.m_bRejectedDueToTrack = TRUE;
					(newEvents)[k].m_PixelAnalResults.m_bRejByTrackOnSingleCam = TRUE;
					(newEvents)[k].m_PixelAnalResults.m_TrackType = track_type;

					(new_track.m_EventsOnTrack).push_back( ((newEvents)[k]) );
					new_track.calcAverageVelocity();
					if( !oldTracks.Find( (*t) ) ){
						// adding only first track or it is not update of already present track :
						oldTracks.Add( new_track );
					}
									
					if(bReFit){
						// upadting after refiting 
						(t->m_EventsOnTrack).push_back( ((newEvents)[k]) );
						t->a = new_track.a;
						t->b = new_track.b;
						t->m_cam_idx = new_track.m_cam_idx;
					}else{
						(t->m_EventsOnTrack).push_back( ((newEvents)[k]) );
					}
					t->calcAverageVelocity();
					LogNewTrack( new_track, FALSE, track_type );
					LogNewTrack_RADEC( new_track, FALSE, track_type );
					t->FlagEventsSaved();

					break;						
				}
			}
		}
	}


	int allEventsCount=newEvents.size();
	int allFramesWithEvents = ( newEvents.size()>0 );		

	int notMatched = oldEvents.GetNotMatchedEventsCountSkipGen();
	allEventsCount += notMatched;
	allFramesWithEvents += (notMatched>0);


	int nLoops=0;
	while(allEventsCount>0 && allFramesWithEvents>nNumBackFramesForTracks && bHasTracks && nLoops<10){
		// trying to find all tracks :
		bHasTracks = FALSE;
		if(allEventsCount<nMaxEventsToFitTrack){
			// not even trying if to much events - maybe something very bad ...

			if(allFramesWithEvents>=nNumBackFramesHasEventsForTracks){
				// no verify if on track - fit line etc ...
				vector<CccdReport*> eventsPtrList;
			
				for(int i=0;i<newEvents.size();i++){
					if( newEvents[i].IsIdentifiedForTrack() &&
						!(newEvents[i].m_PixelAnalResults.m_bRejByTrackOnSingleCam) &&
						!newEvents[i].m_bGenerated && 
						( CccdReport::FindMinDist( eventsPtrList, newEvents[i] )>gCCDParams.m_nMinDistOfEventsInTrack ) ){
						eventsPtrList.push_back( (CccdReport*)(&(newEvents[i])) );
					}
				}	
				for(int i=0;i<oldEvents.size();i++){
					if( oldEvents[i].IsIdentifiedForTrack() && 
						!(oldEvents[i].m_PixelAnalResults.m_bRejByTrackOnSingleCam) &&
						!oldEvents[i].m_bGenerated && 
						( CccdReport::FindMinDist( eventsPtrList, oldEvents[i] )>gCCDParams.m_nMinDistOfEventsInTrack ) ){
						eventsPtrList.push_back( (CccdReport*)(&(oldEvents[i])) );
					}
				}	

				if( eventsPtrList.size()>0 ){
					CTrackDesc track;
					bHasTracks = FlagEventsOnTrack( eventsPtrList, track, in.pPipeline,
															  bCheckVelocity, fVelocityError,
						                             MaxChi2InTrack, nMinEventNoInTrack,
						                             MaxChi2ForPointToMatchLine, track_type );
					if(bHasTracks){
						bRet = TRUE;
						track.m_cam_idx = in.ccd_index;
						if( gCCDParams.m_bLogTracks ){
                 		LogNewTrack( track, TRUE, track_type );
							LogNewTrack_RADEC( track, TRUE, track_type );
                  }
						oldTracks.push_back( track );
					}
				}
			}
		}
		nLoops++;

		if(bHasTracks){
			// calculating not matched tracks :	
			allEventsCount=newEvents.GetNotMatchedEventsCountSkipGen();
			int notMatched = oldEvents.GetNotMatchedEventsCountSkipGen();

			allEventsCount += notMatched;
		}
	}		

	if(nLoops>=10){
		printf("Too many loops in track searching ...\n");
	}

	return bRet;	
}




BOOL_T CCD_Analyser::CheckLineForEventsOnCurrFrame( CCDEventList& frameEvents )
{
	if( frameEvents.size()<3 )
		return FALSE;

	double* x_values = new double[frameEvents.size()];
	double* y_values = new double[frameEvents.size()];
	
	for( register int i=0;i<frameEvents.size();i++){
      CccdReport* pEvt = &(frameEvents[i]);
		x_values[i] = (pEvt->m_Point).x;
		y_values[i] = (pEvt->m_Point).y;
	}

	double a,b;
	double chi2 = CMyFit::FitLineChi2( x_values, y_values, frameEvents.size(), a, b );

	double chi2_point = (chi2/frameEvents.size());

	// clean 
	delete [] x_values;
	delete [] y_values;

	if( m_pAnalFoundEventFunc ){
		(*m_pAnalFoundEventFunc)( &chi2_point , eChi2OnCurrFrame, m_pPipeline->GetPipelineIndex(), NULL );
	}
	if( chi2_point<gCCDParams.m_fChi2ForLineOnSingleFrame ){
		CTrackDesc new_track( a, b, 0 );
		new_track.frame_index = m_pPipeline->GetFrameIndex();
		for(int i=0;i<frameEvents.size();i++){
			new_track.m_EventsOnTrack.push_back( (frameEvents[i]) );
		}
		// (pPipeline->m_RejectedByIfMoreTracks).push_back( new_track );

		// TODO : add to track list 
		
		// LogNewTrack( new_track, TRUE, "rejectifmore_newtracks" , "rejectifmore_tracks" );
		return TRUE;
	}
	return FALSE;
}

BOOL_T CCD_Analyser::FitLineToSingleFrameEvents( CCDEventList& recentEvents )
{
	int ccd_index = m_pPipeline->GetPipelineIndex();
	double a,b;
	BOOL_T bRet = FitLineToSingleFrameEvents( recentEvents, ccd_index, 
															gCCDParams.m_fChi2ForLineOnSingleFrame,
															a, b  );
	return bRet;
}

BOOL_T CCD_Analyser::FitLineToSingleFrameEvents( CCDEventList& recentEvents, int ccd_index,
																double fChi2ForLineOnSingleFrame,
																 double& a_line, double& b_line,
																BOOL_T bCheckVelo/*=FALSE*/, double fVeloError/*=0.00*/  )
{
	if( recentEvents.size()<3 )
		return FALSE;

	sEventDesc* events = new sEventDesc[recentEvents.size()];
	sEventDesc* on_line = new sEventDesc[recentEvents.size()];
	sEventDesc* rejected = new sEventDesc[recentEvents.size()];
	int rejected_cnt=0;
	double minChi2_On3PointsPerPoint;

	int cnt_on_line=0;
	double a,b,maxchi2_out;
	BOOL_T bRet = FALSE;

	int count = recentEvents.size();
	for( register int i=0;i<recentEvents.size();i++){
		CccdReport* pEvt = &(recentEvents[i]);
		events[i].x = pEvt->m_Point.x;
		events[i].y = pEvt->m_Point.y;
		events[i].frame = i; // this is to workaround , that only frames from 	
									// different frames are used 
		events[i].timeUT = pEvt->m_Time;
	}
	if(CMyFit::FindPointsOnLine3New( events, count, fChi2ForLineOnSingleFrame,
										  on_line, cnt_on_line, a, b, maxchi2_out,
										  rejected, rejected_cnt,
										  bCheckVelo, fVeloError ,
										  minChi2_On3PointsPerPoint, ccd_index, TRUE )){
		bRet = TRUE;
		a_line = a;
		b_line = b;
		// now try adding rejected points to line , criteria is less strict
		// for each point gCCDParams.m_MaxChi2ForPointToMatchLine is required
		// and for whole track not to exceed same value per point
		if(rejected_cnt>0){
			double maxchi2_out_prim=0;
			int first_on_line_cnt = cnt_on_line;
			if(CMyFit::TryToAddNewPointsNew( rejected, rejected_cnt,
												gCCDParams.m_fChi2ForLineOnSingleFrame,
												on_line, cnt_on_line, a, b, 				
												maxchi2_out_prim, count, 
											 	NULL, ccd_index )){
				a_line = a;
				b_line = b;
				printf("additional points added");
			}

			// now log found track :
			/*if( gCCDParams.m_bLogTracks ){
				LogNewTrack( track, TRUE );
			}*/
		}
	}
	if( m_pAnalFoundEventFunc ){
		(*m_pAnalFoundEventFunc)( &minChi2_On3PointsPerPoint , eChi2OnCurrFrame3Points, ccd_index, NULL );
	}

	if( m_pAnalFoundEventFunc ){
		// TODO 
	}

	delete [] events;
	delete [] on_line;
	delete [] rejected;

	// maybe try to fit to out of tracks points - maybe there are multiple
	// tracks (many airplanes ...) - so this TO BE DONE here :

	return bRet;
}


BOOL_T CCD_Analyser::FitLine( vector<CccdReport*>& rejectedList, CCDPipeline* pPipeline )
{
	double* x_values = new double[rejectedList.size()];
	double* y_values = new double[rejectedList.size()];
	
	for( register int i=0;i<rejectedList.size();i++){
      CccdReport* pEvt = rejectedList[i];
		x_values[i] = (pEvt->m_Point).x;
		y_values[i] = (pEvt->m_Point).y;
	}

	double a,b;
	double chi2 = CMyFit::FitLineChi2( x_values, y_values, rejectedList.size(), a, b );

	double chi2_point = (chi2/rejectedList.size());

	// clean 
	delete [] x_values;
	delete [] y_values;

	if( chi2_point<gCCDParams.m_MaxChi2ForRejectIfMoreTrack){
		CTrackDesc new_track( a, b, 0 );
		new_track.frame_index = pPipeline->GetFrameIndex();
		for(int i=0;i<rejectedList.size();i++){
			new_track.m_EventsOnTrack.push_back( *(rejectedList[i]) );
		}
		
		LogNewTrack( new_track, TRUE, eRejIfMoreTrack  );
		LogNewTrack_RADEC( new_track, TRUE, eRejIfMoreTrack  );
		(pPipeline->m_CurrFrameTracks).push_back( new_track );
		return TRUE;
	}
	return FALSE;
}


BOOL_T CCD_Analyser::IsSingleFrameTrack( eTrackCheckType_T track_type )
{
	return ( track_type==ePlaneTrack || track_type==eTrackOnSumedFrame || track_type==eSingleCamTrack );
}

BOOL_T CCD_Analyser::DoAddFromSameFrame( eTrackCheckType_T track_type )
{
	return ( track_type==ePlaneTrack || track_type==eTrackOnSumedFrame );
}

BOOL_T CCD_Analyser::CheckIfPlaneTrack(  vector<CccdReport*>& recentEvents,
													  CTrackDesc& track,
	                                      CCDPipeline* pPipeline,
													  int minAddedForPlane )
{
	CTrackDesc tmp_track = track;
	int frame_index = pPipeline->GetFrameIndex();
	int nAdded=0;

	for( int i=0;i<recentEvents.size();i++){
		CccdReport* pEvt = recentEvents[i];

		if(!pEvt->m_PixelAnalResults.m_bRejectedDueToTrack && 
			 pEvt->m_FrameIndex >= (frame_index-gCCDParams.m_nNumBackFramesForPlaneTrack)  ){
			// check if event belongs to track 
			BOOL_T bReFitted=FALSE;
			double new_a,new_b;
			if( CheckIfEventBelongsToTrack( (*pEvt) , tmp_track, bReFitted, 
													  new_a, new_b, FALSE, 0.00 ) ){
				// event belongs to track :
				pEvt->m_PixelAnalResults.m_bRejectedDueToTrack = TRUE;
				tmp_track.UpdateTrack( new_a, new_b, (*pEvt) );
				nAdded++;
			}
		}
	}		
	if( nAdded>=minAddedForPlane ){
		track = tmp_track;
	}	

	return ( nAdded>=minAddedForPlane );	
}

BOOL_T CCD_Analyser::FlagEventsOnTrack( vector<CccdReport*>& recentEvents, 
													CTrackDesc& track,
													CCDPipeline* pPipeline,
													BOOL_T bCheckVelocity, double fVelocityError,
													double MaxChi2InTrack, int nMinEventNoInTrack,
													double MaxChi2ForPointToMatchLine,
													eTrackCheckType_T track_type )
{
	if( recentEvents.size()<=2 )
		return FALSE;

	
	int ccd_index = pPipeline->GetPipelineIndex();
	CImageCreator* pGenObj = pPipeline->GetGenObj();
	sEventDesc* events = new sEventDesc[recentEvents.size()];
	sEventDesc* on_line = new sEventDesc[recentEvents.size()];
	sEventDesc* rejected = new sEventDesc[recentEvents.size()];
	int rejected_cnt=0;
	double minChi2_On3PointsPerPoint;

	int cnt_on_line=0;
	double a,b,maxchi2_out;
	BOOL_T bRet = FALSE;

	int count = recentEvents.size();
	for( register int i=0;i<recentEvents.size();i++){
		CccdReport* pEvt = recentEvents[i];
		events[i].x = pEvt->m_Point.x;
		events[i].y = pEvt->m_Point.y;
		events[i].frame = pEvt->m_FrameIndex;
		events[i].timeUT = pEvt->m_Time;
	}

	// TODO - sort events 
	sort_event_list( events , count );

	if(gDebugTracks)printf("CCD_Analyser::FlagEventsOnTrack Checking for tracks in %d of events ...",count);

	int original_cnt_on_line = 0;
	double original_chi2=0;
	BOOL_T bContinue=TRUE;
	
	while( bContinue ){
		if(CMyFit::FindPointsOnLine3New( events, count, MaxChi2InTrack,
										  on_line, cnt_on_line, a, b, maxchi2_out,
										  rejected, rejected_cnt,
										  bCheckVelocity, fVelocityError,
										  minChi2_On3PointsPerPoint, ccd_index,
										  DoAddFromSameFrame(track_type), track_type )){

			original_cnt_on_line = cnt_on_line;
			original_chi2 = minChi2_On3PointsPerPoint;

			if( gDebugTracks ){ 
				mystring szNewTrack;
				szNewTrack << "3 points found, " << cnt_on_line << "points in track=" << track.track_id << " :";
				show_points(szNewTrack.c_str(), on_line, cnt_on_line );
			}

			if( cnt_on_line>=nMinEventNoInTrack ){
				for( int l=0;l<cnt_on_line;l++){
					for( int i=0;i<recentEvents.size();i++){
						CccdReport* pEvt = recentEvents[i];
						if( ((int)pEvt->m_Point.x==(int)on_line[l].x) &&
							 ((int)pEvt->m_Point.y==(int)on_line[l].y) && 
							 ((int)pEvt->m_FrameIndex==(int)on_line[l].frame) ){
							pEvt->m_PixelAnalResults.m_bRejectedDueToTrack = TRUE;
							if( IsSingleFrameTrack( track_type ) )
								pEvt->m_PixelAnalResults.m_bRejByTrackOnSingleCam = TRUE;
							pEvt->m_PixelAnalResults.m_TrackType = track_type;
							if( !gCCDParams.GetMC() || !pEvt->m_bGenerated ){
								backgrStat.DecTrackAndCoicEvents( pEvt->m_DayFrameIndex );
							}else{
								if( gCCDParams.GetMC() && pEvt->m_bGenerated ){
									samplesStat.DecTrackAndCoicEvents( pEvt->m_DayFrameIndex );
								}
							}
							track.m_EventsOnTrack.push_back( *pEvt );
							bRet = TRUE;
							break;
						}
					}
				}
				track.calcAverageVelocity();			

				// now try adding rejected points to line , criteria is less strict
				// for each point gCCDParams.m_MaxChi2ForPointToMatchLine is required
				// and for whole track not to exceed same value per point
				if(rejected_cnt>0){
					double maxchi2_out_prim=0;
					int first_on_line_cnt = cnt_on_line;
					if(CMyFit::TryToAddNewPointsNew( rejected, rejected_cnt,
														MaxChi2ForPointToMatchLine,
														on_line, cnt_on_line, a, b, 				
														maxchi2_out_prim, count, 
													 	m_pAnalFoundEventFunc, pPipeline->GetPipelineIndex(),
														DoAddFromSameFrame(track_type), track_type )){
						for( int l=first_on_line_cnt;l<cnt_on_line;l++){
							for( int i=0;i<recentEvents.size();i++){
								CccdReport* pEvt = recentEvents[i];

								if(!pEvt->m_PixelAnalResults.m_bRejectedDueToTrack){
									if( ((int)pEvt->m_Point.x==(int)on_line[l].x) && 
										 ((int)pEvt->m_Point.y==(int)on_line[l].y) && 
										 ((int)pEvt->m_FrameIndex!=(int)on_line[l].frame || DoAddFromSameFrame(track_type) ) && 
										(!bCheckVelocity || track.CheckVelocity( pEvt, fVelocityError, ccd_index )) ){
										pEvt->m_PixelAnalResults.m_bRejectedDueToTrack = TRUE;
										if( IsSingleFrameTrack( track_type ) )
											pEvt->m_PixelAnalResults.m_bRejByTrackOnSingleCam = TRUE;
										pEvt->m_PixelAnalResults.m_TrackType = track_type;
										if( !gCCDParams.GetMC() || !pGenObj || !(pGenObj->IsGenEvent( pEvt->m_MaxPoint )) ){
											backgrStat.DecTrackAndCoicEvents( pPipeline->GetDayFrameCounter() );
										}else{
											if( gCCDParams.GetMC() || pGenObj || (pGenObj->IsGenEvent( pEvt->m_MaxPoint )) ){
												samplesStat.DecTrackAndCoicEvents( pPipeline->GetDayFrameCounter() );
											}
										}
										track.m_EventsOnTrack.push_back( *pEvt );
										break;
									}
								}
							}
						}
					}
				}	

				track.calcAverageVelocity();
				track.SetTrackDesc( a, b, 0 );


				// now log found track :
				//if( gCCDParams.m_bLogTracks ){
				//	LogNewTrack( track, TRUE, track_type );
				//}
			}else{
				// to small number of points in case limit >3 we have to 
				// check if there is no track with worse chi2 but 
				// better number of points :
				if( nMinEventNoInTrack>3 ){
					for( int l=0;l<cnt_on_line;l++){
						for( register int i=0;i<count;i++){
							if( ((int)events[i].x==(int)on_line[l].x) &&
		                   ((int)events[i].y==(int)on_line[l].y) &&
      		             ((int)events[i].frame==(int)on_line[l].frame) ){
								// remove event from events list :
								if( i != (count-1) ){
									memcpy( &(events[i]), &(events[count-1]),sizeof(sEventDesc));
								}
								count--;
								break;
							}
						}
					}
					continue; // continue big loop to skip bContinue=FALSE 
								 // and retry finding track without track that
								 // were found now but had to small number of 
								 // points in track
				}
			}
		}else{
			printf("no satisfying 3 points found\n");		
			original_cnt_on_line = cnt_on_line;
      	original_chi2 = minChi2_On3PointsPerPoint;			
		}

		// track found or not found but no point to continue looking for track
		bContinue=FALSE;
	} // end of while(bContinue)


	if( gShowChi2_NPoints &&  original_chi2<1000.00 ){
		printf("FindPointsOnLine3New : %d %.8f\n",original_cnt_on_line,original_chi2);
	}

	if( m_pAnalFoundEventFunc ){
		(*m_pAnalFoundEventFunc)( &minChi2_On3PointsPerPoint, eMinChi2_On3Points, m_pPipeline->GetPipelineIndex(), NULL );
	}

	delete [] events;
	delete [] on_line;
	delete [] rejected;

	// maybe try to fit to out of tracks points - maybe there are multiple
	// tracks (many airplanes ...) - so this TO BE DONE here :

	return bRet;
}

void CCD_Analyser::GetTrackLogName( eTrackCheckType_T track_type,
												mystring& szNewTrackLog, mystring& szTrackLog,
												mystring& szType, mystring& szSubDir )
{
	szSubDir = "Tracks";
	switch (track_type ){
		case eNormalTrack :
			szNewTrackLog = "newtrackslog";
			szTrackLog = "tracks";
			szType = "NORMAL";
			break;
		case ePlaneTrack :
			szNewTrackLog = "planetrackslog";
			szTrackLog = "planetracks";
			szType = "PLANE";
			break;
		case eSingleCamTrack:
			szNewTrackLog = "singlecamtrackslog";
			szTrackLog = "singlecamtracks";
			szType = "SINGLECAM";
			break;
		case eRejIfMoreTrack :
			szNewTrackLog = "rejifmoretrackslog";
			szTrackLog = "rejifmoretracks";
			szType = "IFMORE";
			break;
		case eTrackOnSumedFrame :
			szNewTrackLog = "newtrackslog";
         szTrackLog = "tracks";
         szType = "SUMFRAMES";
			szSubDir = "TracksOnSum";
			break;
		default :
			szNewTrackLog = "newtrackslog";
			szTrackLog = "tracks";
			szType = "NORMAL";			
         break;	
	}
}

void CCD_Analyser::LogNewTrack( CTrackDesc& track, BOOL_T bNew, eTrackCheckType_T track_type )
{
	if( gCCDParams.m_bLogTracks ){
		if(track.m_EventsOnTrack.size()>0){
			mystring szTrackDesc( 2048, 2048 );
			mystring szNewTrackLog = "newtrackslog";
			mystring szTrackLog    = "tracks";
			mystring szType,szSubDir="Tracks";

			GetTrackLogName( track_type, szNewTrackLog, szTrackLog, szType, szSubDir );
			
			vector<CEventBaseInfo>::iterator evt;
			szTrackDesc << track.a << " " << track.b << " " << track.track_id 
							<< " " << track.vx << " " << track.vy;
			if(bNew){
				szTrackDesc << " N\t";
			}else{
				szTrackDesc << " U\t";
			}
			/*for(evt=track.m_EventsOnTrack.begin();evt!=track.m_EventsOnTrack.end();evt++){
				szTrackDesc << (int)evt->m_FrameIndex << "-(" << (int)evt->m_MaxPoint.x << ","
								<< (int)evt->m_MaxPoint.y << "),";
			}*/
			mystring szEvtList;
			// track.GetSortedEventsList( szEvtList );
			track.SortEventListByTime();
			track.GetEventsList( szEvtList );
			szTrackDesc << szEvtList;
			CCDLog trackLog( "%s\n", "a b ID vx vy Oper Events_on_Track",szSubDir.c_str(), szNewTrackLog );
			trackLog.DumpToFile1( m_pPipeline , szTrackDesc.c_str() );

			if( bNew ){
				// tracking each track only once here :
				CCDLog trackLog( "%.4f %.4f %d %.4f %.4f %s\n", "a b ID vx vy Type",szSubDir.c_str(),szTrackLog);
				trackLog.DumpToFile1( m_pPipeline , track.a, track.b, 
											 track.track_id, track.vx, track.vy,
											 szType.c_str() );
			}
		}
	}

//	if( gCCDParams.m_bSaveTracksToDB && m_pPipeline ){
//		int night = m_pPipeline->GetNight();
//		SaveTrack( track, night, bNew, track_type );
//	}
}

void CCD_Analyser::LogNewTrack_RADEC( CTrackDesc& track, BOOL_T bNew, eTrackCheckType_T track_type, BOOL_T bDoConvert )
{
	 LogNewTrack_RADEC( track, track_type, "", bNew, m_pPipeline, FALSE,
   	                  this, TRUE, bDoConvert );  
}	


void CCD_Analyser::LogNewTrack_RADEC( CTrackDesc& track, eTrackCheckType_T track_type,
												  const char* szFileName, BOOL_T bNew,
												  CCDPipeline* pPipeline, BOOL_T bLocal,
												  CCD_Analyser* pAnal, BOOL_T bDoFit,
												  BOOL_T bDoConvert )
{
	if( gCCDParams.m_bLogTracks ){
		if(track.m_EventsOnTrack.size()>0){
			if( bDoFit ){
				double* ra_tab = new double[track.m_EventsOnTrack.size()];
				double* dec_tab = new double[track.m_EventsOnTrack.size()];
				int cnt = track.m_EventsOnTrack.size();

				track.get_radec_values( ra_tab, dec_tab );

				// converting RAD -> DEG , to have line fit in DEGREES :
				for( register int i=0;i<cnt;i++){
					ra_tab[i]  = ( ra_tab[i]*RAD_TO_DEG );
					dec_tab[i] = ( dec_tab[i]*RAD_TO_DEG );
				}

				CMyFit::FitLineChi2( ra_tab, dec_tab, track.m_EventsOnTrack.size(), track.radec_a, track.radec_b );

				double ra_end = AstroAngle::rad2deg( track.m_EventsOnTrack.back().m_AstroCoord.RA );
				double ra_start = AstroAngle::rad2deg( track.m_EventsOnTrack.front().m_AstroCoord.RA );
				double dec_end = track.m_EventsOnTrack.back().m_AstroCoord.Dec;
				double dec_start = track.m_EventsOnTrack.front().m_AstroCoord.Dec;

				track.v_ra = ( ra_end - ra_start ) / (track.m_EventsOnTrack.back().m_Time-track.m_EventsOnTrack.front().m_Time);
				track.v_dec = ( dec_end - dec_start ) / (track.m_EventsOnTrack.back().m_Time-track.m_EventsOnTrack.front().m_Time);
				delete [] ra_tab;	
				delete [] dec_tab;
			}

			mystring szTrackDesc( 2048, 2048 );
			mystring szNewTrackLog;			
			mystring szEvtList;
			mystring szTrackLog="tracks_radec_ccd0.txt";
			mystring szType="N";
			if( szFileName && szFileName[0] )
				szNewTrackLog = szFileName;

			// GetTrackLogName( track_type, szNewTrackLog, szTrackLog, szType, szSubDir );
			
			// track.GetSortedEventsList( szEvtList );
			char szInfo[1024];
			if( sprintf(szInfo,"%.4f %.4f %d %.4f %.4f %.4f %.4f %.4f %.4f",
							track.a,track.b,track.track_id,track.vx,track.vy,
 							track.radec_a,track.radec_b,track.v_ra,track.v_dec)>=1024 ){
				printf("buffer size execeeded in LogNewTrack_RADEC\n");
				exit(0);
			}
			szTrackDesc << szInfo;
			if(bNew){
				szTrackDesc << " N\t";
			}else{
				szTrackDesc << " U\t";
			}
			track.GetEventsList_RADEC( szEvtList, bDoConvert );
			szTrackDesc << szEvtList;

			if( bLocal || !pAnal ){
				CCDLog trackLog( "%s\n", "a b ID vx vy radec_a radec_b v_ra v_dec Oper Events_on_Track", "RADEC", szNewTrackLog.c_str() );
				trackLog.DumpLine( 0, szTrackDesc.c_str() );

				CCDLog trackLog2( "%.4f %.4f %d %.4f %.4f %.4f %.4f %.4f %.4f %s\n", 
									  "a b ID vx vy radec_a radec_b v_ra v_dec Type ",
									  NULL, szTrackLog.c_str() );
				trackLog2.DumpLine( 1 , 
								 track.a, track.b,
								 track.track_id, track.vx, track.vy, 
								 track.radec_a, track.radec_b,
								 track.v_ra, track.v_dec,
								 szType.c_str() );
			}else{
				mystring szTrackLog    = "tracks_radec";
				mystring szType="NORMAL",szSubDir="Tracks";
				mystring szNewTrackLog2 = "newtrackslog_radec";

				// if( pAnal ){
				//	pAnal->GetTrackLogName( track_type, szNewTrackLog2, szTrackLog, szType, szSubDir );

					CCDLog trackLog( "%s\n", "a b ID vx vy radec_a radec_b v_ra v_dec Oper Events_on_Track", szSubDir.c_str(), szNewTrackLog2.c_str() );
					trackLog.DumpToFile1( pPipeline , szTrackDesc.c_str() );

					if( bNew ){
						CCDLog trackLog2( "%.4f %.4f %d %.4f %.4f %.4f %.4f %.4f %.4f %s\n", 
										  "a b ID vx vy radec_a radec_b v_ra v_dec Type ",
										  szSubDir.c_str(), szTrackLog.c_str() );
						trackLog2.DumpToFile1(  pPipeline, 
									 track.a, track.b,
									 track.track_id, track.vx, track.vy, 
									 track.radec_a, track.radec_b,
									 track.v_ra, track.v_dec,
									 szType.c_str() );
					}
				// }
			}
		}
	}

}


BOOL_T CCD_Analyser::CheckBiggerArea( const CPixelAnalyseIn& in,
                                      LONG_T* cluster, LONG_T& cluster_cnt,
												  LONG_T* ClusterWithMore, LONG_T& ClusterWithMoreCnt,
                                      CPixelList& pixel_list,BOOL_T bTrace)
{
	LONG_T clusterSum=0;
	LONG_T x0,y0,pos0;
	LONG_T prevClusterSum[MAX_PREV_SUMS];

	if( gCCDParams.m_ConfShape==shapeCluster ){
		clusterSum = FindClusterAboveTresholdOpt( &in, in.p_data, in.x, in.y, in.pos, in.xSize, in.ySize,
   	      	                      				x0, y0, pixel_list, cluster, cluster_cnt );
		if(cluster_cnt<gCCDParams.m_MinClusterSize){
			if(bTrace){
				MYTRACE4(gCCDTrace,"Potential event at x=" << in.x << ",y=" 
   	         	         	  << in.y << ", but cluster size=" 
      	         	           << cluster_cnt << "<" << gCCDParams.m_MinClusterSize);
			}
			return FALSE;
		}
		pos0 = x0 + y0*in.xSize;
		ClusterWithMoreCnt=0;
				
		// cluster is also sorted inside :
		GetClusterWithPointsAroundOpt(cluster,cluster_cnt,
												ClusterWithMore,ClusterWithMoreCnt,
												in.xSize,in.ySize,
												gCCDParams.m_nPixelsAroundToConfirm);
			

		// DumpCluster( ClusterWithMore , ClusterWithMoreCnt , p_data, xSize, ySize, "" );
	}else{
		GetClusterFromShapeMap( in.x, in.y, in.xSize, in.ySize, ClusterWithMore,ClusterWithMoreCnt );
		cluster_cnt = ClusterWithMoreCnt;
		memcpy(cluster,ClusterWithMore,sizeof(LONG_T)*cluster_cnt);
	}

	if(!ConfirmEvent_InClusterOpt_RotCorrected( in, ClusterWithMore, ClusterWithMoreCnt,
											              clusterSum, prevClusterSum )){
		if(bTrace){
			MYTRACE4(gCCDTrace,"Potential event at x=" << x0 << ",y=" << y0 << ", but rejected by prev sum rejection procedure !");
		}
		return FALSE;
	}
	return TRUE;
}


BOOL_T CCD_Analyser::AnalysePixel( LONG_T x, LONG_T y, CCDMatrix& Matrix, CCDPipeline* pPipeline, LONG_T ccd_index,
                     			     mystring& szOutput )
{
	LONG_T ncnt = gCCDParams.m_nNeighbToSumCount;
	CPixelAnalyseIn in;
   CPixelAnalyseOut out;
	
	in.xSize = Matrix.GetXSize();
   in.ySize = Matrix.GetYSize();
	CPixelList pixel_list(in.xSize*in.ySize);
	

	in.x = x;
	in.y = y;
	in.pos = (in.y*in.xSize + in.x);
	in.Matrix = &Matrix;
	in.pPipeline = pPipeline;
	in.ccd_index = ccd_index;
	in.pCCDInfo = &(pPipeline->GetCCDInfoTab()[ccd_index]);
	pPipeline->GetAllMatrixPtrsChronologicalInt( ccd_index, in, FALSE );
	in.pipeline_size_minus_1 = pPipeline->GetPipelineSize()-1;	
	in.p_data = Matrix.get_data_buffer();
	in.p_data_fast = Matrix.get_data_buffer_fast();
   if(gCCDParams.m_bKeepLaplaceFrame){
      CManyTab2D<BIG_ELEM_TYPE>* pLaplaceFrame = pPipeline->GetLaplaceFrame();
		in.p_laplace_data = (*pLaplaceFrame)[0].get_data_buffer();
		in.p_laplace_data_fast = (*pLaplaceFrame)[0].get_data_buffer_fast();
	}		

		
	out.ncnt = gCCDParams.m_nNeighbToSumCount;
	out.pixel_list = &pixel_list;
	out.m_PixelOut.eventType = eNone;

	if(gCCDParams.m_bKeepHomeopaticSum && pPipeline->GetHomeopaticFrame()){
		in.p_homeo_data = (*(pPipeline->GetHomeopaticFrame()))[ccd_index].get_data_buffer();
   }

	BOOL_T bRet = AnalysePixel( in, out, TRUE );

	szOutput = "";
	szOutput << "Point (x,y)=(" << in.x << "," << in.y << ") was ";	
	if(bRet)
		szOutput << " ACCEPTED\n";
	else
		szOutput << " REJECTED\n";

	mystring szMethod;
	if(gCCDParams.m_bAnalyseMaxFromPrevInPipeline){
		szMethod = " PREV_SUM ";
   }
   if(gCCDParams.m_bCheckHomeoRawCond){
		szMethod = " HOMEOPATIC ";
   }
	szOutput << "At first level trigger " << szMethod << " method was used, using " << out.ncnt << " pixels \n";
 	szOutput <<  "For values : newSum=" << out.m_PixelOut.newSum << ", otherSum=" << out.m_PixelOut.otherSum << "\n";	
	szOutput << "condition was : ";
	szOutput << "newSum>" << (out.ncnt*gCCDParams.m_TresholdPerPixel)
            << " AND otherSum<" << (out.ncnt*gCCDParams.m_MaxPrevPerPixel) << "\n";
	if(gCCDParams.m_bConfirmReq && (out.cluster_cnt || out.ClusterWithMoreCnt)){
		szOutput << "Second level trigger used cluster of " << MAX(out.ClusterWithMoreCnt,out.cluster_cnt) << " pixels";
		szOutput << " new_cluster_sum=" << out.newClusterSum << "\n";
		szOutput << "Previous sums were : ";
		for(int i=0;i<gCCDParams.m_FramesBack && i<MAX_PREV_SUMS;i++){
			szOutput << (out.prevClusterSum)[i] << ",";
		}
		szOutput << "\n";
		szOutput << "Condition : PREV<" << (out.cluster_cnt*(gCCDParams.m_ConfMaxPrevPerPixel)) << " AND NEW>"
					<< (out.cluster_cnt*(gCCDParams.m_ConfTresholdPerPixel)) << "\n";				
	}			
	return bRet;
}



BOOL_T CCD_Analyser::AnalysePixel( CPixelAnalyseIn& in,
                                   CPixelAnalyseOut& out,
                                   BOOL_T bTrace )
{
	GetAnalNeighbWithSelf(in.x,in.y,in.pos,in.xSize,in.ySize,
								 out.neighb_list,out.ncnt);

	out.m_PixelOut.Init();
	out.m_PixelOut.x0 = in.x;
   out.m_PixelOut.y0 = in.y;
   out.m_PixelOut.pos0 = in.pos;
	out.m_PixelOut.treshold_for_not_more_then_n_exceeds = in.treshold_for_not_more_then_n_exceeds;

	if(gCCDParams.m_bCheckLaplaceCondition){
		if(in.p_curr_data_laplace){
			out.m_PixelOut.laplaceSum = in.p_curr_data_laplace[in.y][in.x];
			//Assert( in.p_curr_data_laplace[in.y][in.x]==out.laplaceSum,"Error in laplacjan calculation procedure at (%d,%d), laplaceSum=%f != %f !!!"
			//			,in.x, in.y,(double)out.laplaceSum,(double)in.p_curr_data_laplace[in.y][in.x]);
		}else{
			 out.m_PixelOut.laplaceSum = Table2D<ELEM_TYPE>::CalcLaplaceSum( in.x, in.y , in.xSize, in.p_data_fast, gCCDParams.m_eLaplaceType );			
		}

		out.m_PixelOut.newSum = out.m_PixelOut.laplaceSum;

		// try to reject as ealrly as possible :
		if(out.m_PixelOut.laplaceSum<=((in.pPipeline)->GetPipelineCfg()).m_nNewLaplace)
			return FALSE;

		if(in.p_laplace_data){
			out.m_PixelOut.prevLaplaceSum = in.p_laplace_data[in.pos];
		}
		if(gCCDParams.m_bAnalyseMaxFromPrevInPipeline){
			out.m_PixelOut.prevLaplaceSum = GetMaxLaplaceOnPrevRotCorrect( in );						
		}
		if(gCCDParams.m_nMaxOfAverageOfPrevN>0){
			out.m_PixelOut.maxAverageOfPrev = GetMaxAverageInVetoArea( in, TRUE );
		}

		out.m_PixelOut.otherSum = out.m_PixelOut.prevLaplaceSum;
		out.m_PixelOut.maxSum = out.m_PixelOut.prevLaplaceSum;
	}else{
		if(gCCDParams.m_bAnalyseMaxFromPrevInPipeline || gCCDParams.m_bCheckHomeoRawCond){
			out.m_PixelOut.newSum = CalcSumOpt( out.neighb_list, out.ncnt, in.p_data );
			out.m_PixelOut.newRawSum = out.m_PixelOut.newSum;
		}
		if(gCCDParams.m_bAnalyseMaxFromPrevInPipeline){			
			out.m_PixelOut.maxSum = CalcMaxSumOnPrevRotCorrect( out.neighb_list, out.ncnt, in );
			out.m_PixelOut.otherSum = out.m_PixelOut.maxSum; 
		}
		if(gCCDParams.m_bKeepHomeopaticSum && in.p_homeo_data){
			out.m_PixelOut.homeoSum = (in.p_homeo_data)[in.pos];
			out.m_PixelOut.otherSum = out.m_PixelOut.homeoSum;
		}				
		if(gCCDParams.m_nMaxOfAverageOfPrevN>0){
			out.m_PixelOut.maxAverageOfPrev = GetMaxAverageInVetoArea( in );
		}
	}		

	long x_y_val = (in.p_data)[in.pos];

	if(x_y_val>=gCCDParams.m_MaxAllowedVal){
		// skiping pixel exceeding maximum allowed values in
		return FALSE;
	}

	BOOL_T bChecked = FALSE;
	BOOL_T bFirstCheck = TRUE;

	if(gCCDParams.m_bCheckLaplaceCondition){
		bFirstCheck = (bFirstCheck && CheckLaplaceCondition( in, out ));
		if(!bFirstCheck)
			return FALSE;
		bChecked=TRUE;
		if( bFirstCheck && gCCDParams.m_bConfirmLaplaceMedianMinus){
			bFirstCheck = CheckLaplaceConditionMedianMinus( in, out );
			if(!bFirstCheck)
				return FALSE;
		}		
		if( bFirstCheck && gCCDParams.m_nCheckIfNoneOfNExceedsTresh>=0){
			bFirstCheck = VerifyIfNotMoreThenNExceedsTreshold( in , out, gCCDParams.m_nNotMoreThenNExceedsBackFramesCount,
																				gCCDParams.m_nCheckIfNoneOfNExceedsTresh,
																				TRUE, gCCDParams.mNotMoreThenNExceedsArea, 
																			   gCCDParams.m_nNotMoreThenNExceedsPointsCount );
			if(!bFirstCheck)
				return FALSE;
		}
	}else{
		if(gCCDParams.m_bAnalyseMaxFromPrevInPipeline){
			bFirstCheck = (bFirstCheck && CheckCondition(out.m_PixelOut.newSum,out.m_PixelOut.maxSum) );
			bChecked=TRUE;
		}
		if(gCCDParams.m_nMaxOfAverageOfPrevN>0){
			bFirstCheck = (bFirstCheck && CheckCondition(out.m_PixelOut.newSum,out.m_PixelOut.maxAverageOfPrev) );
			bChecked=TRUE;
		}		
	}

	if(gCCDParams.m_bCheckHomeoRawCond){
		// printf("Checking homeopatic raw condition !!!\n");
		
		bFirstCheck = (bFirstCheck && CheckCondition(out.m_PixelOut.newSum,out.m_PixelOut.homeoSum));
		bChecked=TRUE;
	}


	if(gCCDParams.m_bLocalMaxReq){
		LONG_T max_pos=-1;
		if(bFirstCheck){
			BOOL_T bLocalMaxCheck = CheckLocalMaxCondition(in, out, max_pos); 		
			bFirstCheck = (bFirstCheck && bLocalMaxCheck);
			if(!bLocalMaxCheck){
				if(max_pos>=0){
					MYTRACE3( gCCDTrace,"(x,y)=" << in.x << "," << in.y
      	                        << " - rejected by local MAX requirement max_pos=(" 
											<< (max_pos % in.xSize) << "," << (max_pos / in.xSize) << ")");
				}
			}
		}
	}

	/*if( bFirstCheck && bChecked ){
		if( m_pAnalFoundEventFunc ){
			bFirstCheck = (*m_pAnalFoundEventFunc)( *(in.Matrix), (long)in.x, (long)in.y, NULL );
			if(!bFirstCheck){
				out.m_PixelOut.m_bRejectedByCustomEventAnal = TRUE;
			}
		}
	}*/

	if( bFirstCheck && bChecked ){
		out.m_PixelOut.eventType = eFlash;
	}

	if(gCCDParams.m_bCheckForSUPERNEW){
		if( !bFirstCheck || !bChecked ){
			// check if it was brightening !
	   	if( out.m_PixelOut.otherSum >= gCCDParams.m_MinPrevTotal ){
	      	if( CheckBrighteningCond( out.m_PixelOut.newSum, out.m_PixelOut.otherSum ) ){
					out.m_PixelOut.eventType = eBrighten;
				}
			}		
		}
	}

	if( out.m_PixelOut.eventType!=eNone ){
		if(bTrace)
			MYTRACE4( gCCDTrace,"First Level trigger accepted event at x="
                   << in.x << " ,y=" << in.y 
   	             << " condition :" << GetCondDesc(out.m_PixelOut.newSum,out.m_PixelOut.maxSum,out.m_PixelOut.homeoSum));

		// initializing cluster to neighbouring points :
		out.ClusterWithMoreCnt = out.ncnt;
		memcpy(out.ClusterWithMore,out.neighb_list,out.ncnt*sizeof(LONG_T));
		out.cluster_cnt = out.ncnt;
		memcpy(out.cluster,out.neighb_list,out.ncnt*sizeof(LONG_T));

		if(( gCCDParams.m_bConfirmReq && gCCDParams.m_ConfShape==shapeCluster) 
			  || gCCDParams.m_bCalcClusterReq){
				out.newClusterSum = FindClusterAboveTresholdOpt( &in, in.p_data, 
                                      in.x, in.y, in.pos, in.xSize, in.ySize,
  	                      				  out.m_PixelOut.x0, out.m_PixelOut.y0, *(out.pixel_list), 
                                      out.cluster, out.cluster_cnt );
				out.ClusterWithMoreCnt = out.cluster_cnt;
				memcpy(out.ClusterWithMore,out.cluster,out.cluster_cnt*sizeof(LONG_T));
				out.m_PixelOut.pos0 = out.m_PixelOut.x0 + out.m_PixelOut.y0*in.xSize;

				if(out.cluster_cnt<gCCDParams.m_MinClusterSize){
					if(bTrace)
					 	MYTRACE4(gCCDTrace,"Potential event at x=" << in.x << ",y=" 
  		       	   	               << in.y << ", but cluster size=" 
     	         	   	            << out.cluster_cnt << "<" << gCCDParams.m_MinClusterSize);
					return FALSE;
				}
		}

		if(gCCDParams.m_bConfirmReq){
			BOOL_T bRet = FALSE;

			if(out.m_PixelOut.eventType==eFlash){
				// required confirmation - 1.5 level of trigger :-)
				LONGLONG_T totalSum=0;

				if( gCCDParams.m_ConfShape==shapeCluster ){
					if(gCCDParams.m_bUseClusterWithMore){
						out.ClusterWithMoreCnt=0;				
						// cluster is also sorted inside :
						GetClusterWithPointsAroundOpt(out.cluster,out.cluster_cnt,
																out.ClusterWithMore,out.ClusterWithMoreCnt,
																in.xSize,in.ySize,
																gCCDParams.m_nPixelsAroundToConfirm);			
					}
				}else{
					GetClusterFromShapeMap( in.x, in.y, in.xSize, in.ySize,
                                       out.ClusterWithMore,out.ClusterWithMoreCnt );
					out.cluster_cnt = out.ClusterWithMoreCnt;
					memcpy(out.cluster,out.ClusterWithMore,sizeof(LONG_T)*out.cluster_cnt);
				}

				if(!ConfirmEvent_InClusterOpt_RotCorrected( in, out.ClusterWithMore, out.ClusterWithMoreCnt,
																	  out.newClusterSum, out.prevClusterSum )){
					if(bTrace)
						MYTRACE4(gCCDTrace,"Potential event at x=" << out.m_PixelOut.x0 
								<< ",y=" << out.m_PixelOut.y0 << ", but rejected by prev sum rejection procedure !");
					return FALSE;
				}
				bRet = TRUE;
			}

			if(out.m_PixelOut.eventType==eBrighten){
				ConfirmEvent_InClusterOpt_RotCorrected( in, out.cluster,out.cluster_cnt,
																    out.newClusterSum, out.prevClusterSum );
				LONGLONG_T maxPrevSum = my_find_max_value_long( out.prevClusterSum, in.pipeline_size_minus_1);
				BOOL_T bCheckOK = CheckBrighteningCond( out.newClusterSum, maxPrevSum );
				if(bCheckOK){
					out.ClusterWithMoreCnt = out.cluster_cnt;
	            memcpy(out.ClusterWithMore,out.cluster,sizeof(LONG_T)*out.cluster_cnt);
					bRet = TRUE;
				}
			}


				// mark all pixels in cluster :
			if( gCCDParams.m_ConfShape==shapeCluster ){	
				// only in custer above NOISE_TRESHOLD
				for(int ii=0;ii<out.cluster_cnt;ii++){
					(out.pixel_list)->HitPixel( (out.cluster)[ii] );
				}
			}else{
				// in this case ClueterWithMore contains full cluster 
				// which is - square, circle, etc
				for(int ii=0;ii<out.ClusterWithMoreCnt;ii++){
					(out.pixel_list)->HitPixel( (out.ClusterWithMore)[ii] );
				}
			}
			return bRet;
		}

		// no second level is required :
		return TRUE;
	}
	return FALSE;

}

BOOL_T CCD_Analyser::AnalysePixelOpt3( const CPixelAnalyseIn& in,
	                                    CPixelAnalyseOut& out,
                                       BOOL_T bTrace )
{
	/*out.ncnt = 5;
	out.neighb_list[0] = in.pos-in.xSize;
	out.neighb_list[1] = in.pos-1;
	out.neighb_list[2] = in.pos;
	out.neighb_list[3] = in.pos+1;
	out.neighb_list[4] = in.pos+in.xSize;*/

	GetAnalNeighbWithSelf(in.x,in.y,in.pos,in.xSize,in.ySize,
		                   out.neighb_list,out.ncnt);

	out.m_PixelOut.Init();
	out.m_PixelOut.x0 = in.x;
	out.m_PixelOut.y0 = in.y;
	out.m_PixelOut.pos0 = in.pos;

	if(gCCDParams.m_bCheckLaplaceCondition){
		// out.laplaceSum = CalcLaplaceSum( in.x, in.y , in.xSize, in.p_data_fast );
		
		if(in.p_curr_data_laplace){
			out.m_PixelOut.laplaceSum = in.p_curr_data_laplace[in.y][in.x];
			//Assert( in.p_curr_data_laplace[in.y][in.x]==out.laplaceSum,"Error in laplacjan calculation procedure at (%d,%d), laplaceSum=%f != %f !!!"
			//			,in.x, in.y,(double)out.laplaceSum,(double)in.p_curr_data_laplace[in.y][in.x]);
		}

		out.m_PixelOut.newSum = out.m_PixelOut.laplaceSum;
		if(!CheckLaplaceCondition( in, out ))
			return FALSE;
		if(in.p_laplace_data){
			out.m_PixelOut.prevLaplaceSum = in.p_laplace_data[in.pos];
		}else{
			if(gCCDParams.m_bAnalyseMaxFromPrevInPipeline){
				// out.prevLaplaceSum = GetMaxLaplaceOnPrevRotCorrect( in );				
				out.m_PixelOut.prevLaplaceSum = in.p_max_prev_lap[in.y][in.x];
			}
		}
		out.m_PixelOut.otherSum = out.m_PixelOut.prevLaplaceSum;
		out.m_PixelOut.maxSum = out.m_PixelOut.prevLaplaceSum;
	}else{
		/*if((gCCDParams.m_bAnalyseMaxFromPrevInPipeline && gCCDParams.m_bCheckHomeoRawCond) ||
			gCCDParams.m_bAnalyseMaxFromPrevInPipeline){
			GetAnalNeighbWithSelf(in.x,in.y,in.pos,in.xSize,in.ySize,
		                         out.neighb_list,out.ncnt);
		}*/

		if(gCCDParams.m_bAnalyseMaxFromPrevInPipeline && gCCDParams.m_bCheckHomeoRawCond){
			out.m_PixelOut.newSum = CalcSumOpt( out.neighb_list, out.ncnt, in.p_data );
		}
		if(gCCDParams.m_bAnalyseMaxFromPrevInPipeline){			
			out.m_PixelOut.maxSum = CalcMaxSumOnPrevRotCorrect( out.neighb_list, out.ncnt, in );
			out.m_PixelOut.otherSum = out.m_PixelOut.maxSum; 
		}
		if(gCCDParams.m_bKeepHomeopaticSum && in.p_homeo_data){
			out.m_PixelOut.homeoSum = (in.p_homeo_data)[in.pos];
			out.m_PixelOut.otherSum = out.m_PixelOut.homeoSum;
		}				
	}		

	long x_y_val = (in.p_data)[in.pos];

	if(x_y_val>=gCCDParams.m_MaxAllowedVal){
		// skiping pixel exceeding maximum allowed values in
		return FALSE;
	}

	BOOL_T bChecked = FALSE;
	BOOL_T bFirstCheck = TRUE;

	if(gCCDParams.m_bCheckLaplaceCondition){
		bFirstCheck = (bFirstCheck && CheckLaplaceCondition( in, out ));
		bChecked=TRUE;
	}else{
		if(gCCDParams.m_bAnalyseMaxFromPrevInPipeline){
			bFirstCheck = (bFirstCheck && CheckCondition(out.m_PixelOut.newSum,out.m_PixelOut.maxSum) );
			bChecked=TRUE;
		}
	}

	if(gCCDParams.m_bCheckHomeoRawCond){
		// printf("Checking homeopatic condition !!!\n");
		
		bFirstCheck = (bFirstCheck && CheckCondition(out.m_PixelOut.newSum,out.m_PixelOut.homeoSum));
		bChecked=TRUE;
	}


	if(gCCDParams.m_bLocalMaxReq){
		LONG_T max_pos=-1;
		if(bFirstCheck){
			BOOL_T bLocalMaxCheck = CheckLocalMaxCondition(in, out, max_pos); 		
			bFirstCheck = (bFirstCheck && bLocalMaxCheck);
			if(!bLocalMaxCheck){
				if(max_pos>=0){
					MYTRACE3( gCCDTrace,"(x,y)=" << in.x << "," << in.y
      	                        << " - rejected by local MAX requirement max_pos=(" 
											<< (max_pos % in.xSize) << "," << (max_pos / in.xSize) << ")");
				}
			}
		}
	}

	if( bFirstCheck && bChecked ){
		out.m_PixelOut.eventType = eFlash;
	}

	if(gCCDParams.m_bCheckForSUPERNEW){
		if( !bFirstCheck || !bChecked ){
			// check if it was brightening !
	   	if( out.m_PixelOut.otherSum >= gCCDParams.m_MinPrevTotal ){
	      	if( CheckBrighteningCond( out.m_PixelOut.newSum, out.m_PixelOut.otherSum ) ){
					out.m_PixelOut.eventType = eBrighten;
				}
			}		
		}
	}
	/*if( out.eventType!=eNone ){
		// initializing cluster to neighbouring points :
		out.ClusterWithMoreCnt = out.ncnt;
		memcpy(out.ClusterWithMore,out.neighb_list,out.ncnt*sizeof(LONG_T));
		out.cluster_cnt = out.ncnt;
		memcpy(out.cluster,out.neighb_list,out.ncnt*sizeof(LONG_T));
	}*/

	return (out.m_PixelOut.eventType!=eNone);
}


BOOL_T CCD_Analyser::FullAnalysePixel( CPixelAnalyseIn& in,
                                       CPixelAnalyseOut& out,
                                       BOOL_T bTrace/*=FALSE*/												
											 )
{
	out.m_PixelOut.eventType = eNone;
	cCCD& newFrame = (in.pPipeline)->GetCurrent();
	BOOL_T bCheckPixel = AnalysePixel( in, out, bTrace );

	if(bCheckPixel){
		// check overlaping with ealier detected event
		LONG_T prev_x,prev_y;
		if( (in.Matrix)->CheckIfOverlaps( out.cluster, out.cluster_cnt, prev_x, prev_y ) ){
			MYTRACE4(gCCDTrace,"Cluster overlaps with already detected event at (" << prev_x << "," << prev_y << ")" );
			bCheckPixel = FALSE;
		}
	}else{
		// in case not considered as flash - check if exisiting
		// object flash - variable objects 
		if( gCCDParams.m_bCheckForSUPERNEW ){
		 	if( out.m_PixelOut.otherSum >= gCCDParams.m_MinPrevTotal ){
					BOOL_T bCheckOK = FALSE;
					LONG_T increase = (out.m_PixelOut.newSum - out.m_PixelOut.otherSum);
					if(increase>0){
						if( gCCDParams.m_eBrigtheningCheckType==ePercentIncrease )
						{
							if( (double(increase)/double(out.m_PixelOut.otherSum)) > gCCDParams.m_IncreaseTreshPercent ){
								bCheckOK=TRUE;
							}
						}else{
							if( increase>gCCDParams.m_IncreaseTreshADU )
								bCheckOK=TRUE;
						}
					}
					if(bCheckOK){
						bCheckPixel=TRUE;
						out.m_PixelOut.eventType = eBrighten;
					}
			}
		}
	}						

	if(bCheckPixel){
			MYTRACE2(gCCDTrace,"ACCEPTED : Potential event found x=" << out.m_PixelOut.x0 << ", y="<< out.m_PixelOut.y0
								<<" ccd_index="<<in.ccd_index<<", sum="<<out.m_PixelOut.newSum<<", max_sum="<<out.m_PixelOut.maxSum);

			(in.Matrix)->AddFoundEvent( in, out );
			MYTRACE2(gCCDTrace,"ACCEPTED : Potential " 
							<< GetEventTypeDesc(out.m_PixelOut.eventType) <<  " event found x=" 
							<< in.x << ", y="<<in.y<<" ccd_index="<<in.ccd_index<<", sum="<<out.m_PixelOut.newSum
							<<", max_sum="<<out.m_PixelOut.maxSum);
			(in.Matrix)->SetInteresting();
			newFrame.SetInteresting();
	}
#ifdef _MONTE_CARLO_
	else{			
			// posibility to turn off - in case no event detailed report
			// is nessesary for us - for big parameters anlysis can make it slow
			if(gCCDParams.m_bGenNotFoundReport && gCCDParams.GetMC()){
				// in simulation version - update generated event info
				// - why it was refused !
				(in.Matrix)->UpdateGenEvents( in, out );
			}
	}
#endif

	return bCheckPixel;
}


BOOL_T CCD_Analyser::CheckBrighteningCond( LONGLONG_T newSum, LONGLONG_T otherSum )
{
	LONG_T increase = (newSum - otherSum);
	BOOL_T bCheckOK = FALSE;

	if(increase>0){
		if( gCCDParams.m_eBrigtheningCheckType==ePercentIncrease )
		{
			if( (double(increase)/double(otherSum)) > gCCDParams.m_IncreaseTreshPercent ){
				bCheckOK=TRUE;
			}
		}else{
			if( increase>gCCDParams.m_IncreaseTreshADU )
				bCheckOK=TRUE;
		}
	}
	return bCheckOK;
}

BOOL_T CCD_Analyser::CheckBrightening( const CPixelAnalyseIn& in,
                                   CPixelList& pixel_list,CCDMatrix& Matrix,
											  LONG_T newSum, LONG_T maxSum,
											  LONG_T homeoSum,LONG_T otherSum,
											  LONG_T& newClusterSum,LONG_T* prevClusterSum,
											  LONG_T x0, LONG_T y0,LONG_T pos0,
											  LONG_T* cluster, LONG_T& cluster_cnt,
											  LONG_T* ClusterWithMore, LONG_T& ClusterWithMoreCnt )
{
	BOOL_T bCheckOK = FALSE;
	if( otherSum >= gCCDParams.m_MinPrevTotal ){
		bCheckOK = CheckBrighteningCond( newSum, otherSum );


		if(bCheckOK){
			FindClusterAboveTresholdOpt( &in, in.p_data, in.x, in.y, in.pos, in.xSize, in.ySize,
   	                           x0, y0, pixel_list, cluster, cluster_cnt );
			if( gCCDParams.m_bConfirmReq ){
				ConfirmEvent_InClusterOpt_RotCorrected( in, cluster,cluster_cnt, newClusterSum, prevClusterSum );
				LONGLONG_T maxPrevSum = my_find_max_value_long( prevClusterSum, in.pipeline_size_minus_1);
				bCheckOK = CheckBrighteningCond( newClusterSum, maxPrevSum );
				if(bCheckOK){
					ClusterWithMoreCnt = cluster_cnt;
	            memcpy(ClusterWithMore,cluster,sizeof(LONG_T)*cluster_cnt);
				}
			}
		}
	}
	return bCheckOK;
}

// 
//  all algorithm applied here :
//
BOOL_T CCD_Analyser::FullNewFrameAnalyse(CCDPipeline& ccd_pipeline)
{
	Initialize();
	LONG_T pipeline_size = ccd_pipeline.GetPipelineSize();
	LONG_T ClusterWithMoreCnt=0;
	LONG_T n_found_events = 0;
	LONG_T n_accepted_on_level_1=0;
	LONG_T n_accepted_on_level_2=0;
   LONG_T star_cluster[MAX_CLUSTER_SIZE];
	LONG_T nFramesBack = gCCDParams.m_FramesBack;

	if(gCCDParams.m_nMaxOfAverageOfPrevN>nFramesBack)
		nFramesBack = gCCDParams.m_nMaxOfAverageOfPrevN;
	if(gCCDParams.m_nNotMoreThenNExceedsBackFramesCount>nFramesBack)
		nFramesBack = gCCDParams.m_nNotMoreThenNExceedsBackFramesCount;

	// to have also current frame
	nFramesBack++;			


	if(ccd_pipeline.GetCount()!=pipeline_size)
		return FALSE;

	BOOL_T bRet = FALSE;

	PROFILER_START
	cCCD& newFrame = ccd_pipeline.GetCurrent();
	long size = newFrame.GetCount();
	if(size<1){
		return FALSE;
	}
	CPixelAnalyseIn in;
   CPixelAnalyseOut out;


	in.PrevMatrixPtrCnt = nFramesBack;
	in.frame_index = ccd_pipeline.GetFrameIndex();
	in.xSize = newFrame[0].GetXSize();
   in.ySize = newFrame[0].GetYSize();	
	long yUpperLimit = (in.ySize-gCCDParams.m_nIgnoreEdgeUp);
	long xUpperLimit = (in.xSize-gCCDParams.m_nIgnoreEdgeRight);		
	m_FirstDy = ((m_NeighbPoints[0].y)*in.xSize);
	in.treshold_for_max = (long)(gCCDParams.m_nCalcMaxForAboveNSigma*CCDDataResults::GetSigmaBackground( gCCDParams.m_eLaplaceType ));
	in.treshold_for_not_more_then_n_exceeds = (long)(gCCDParams.m_TreshForCheckIfNoneOfNExceedsInSigma*CCDDataResults::GetSigmaBackground( gCCDParams.m_eLaplaceType ));

	out.m_PixelOut.treshold_for_not_more_then_n_exceeds = in.treshold_for_not_more_then_n_exceeds;

	// printf("TRESHOLD = %d\n",(long)in.treshold_for_not_more_then_n_exceeds);

	in.pPipeline = &ccd_pipeline;
	in.pipeline_size_minus_1 = pipeline_size-1;

	// treshold for raw cluster :
	in.treshold_for_raw_cluster = (int)CCDDataResults::GetTreshold( eRawS, gCCDParams.m_nSigmaAboveMeanInRawCluster, (in.pPipeline)->GetBackgroundMap() );
	in.treshold_for_raw_cluster_aver = (int)CCDDataResults::GetTresholdAver( eRawS, gCCDParams.m_nSigmaAboveMeanInRawCluster, (in.pPipeline)->GetBackgroundMap() );


	out.ncnt = gCCDParams.m_nNeighbToSumCount;

	CManyTab2D<BIG_ELEM_TYPE>* pHomeoFrame = NULL;
	if(gCCDParams.m_bKeepHomeopaticSum)
		pHomeoFrame = ccd_pipeline.GetHomeopaticFrame();
	CManyTab2D<BIG_ELEM_TYPE>* pLaplaceFrame = NULL;
	if(gCCDParams.m_bKeepLaplaceFrame)
		pLaplaceFrame = ccd_pipeline.GetLaplaceFrame();

	clock_t total_in_check=0;
	for(register int i=0;i<size;i++){
		in.ccd_index = i;
		in.pCCDInfo = &(ccd_pipeline.GetCCDInfoTab()[in.ccd_index]);
		in.Matrix = &(newFrame[i]);
		in.p_data = (in.Matrix)->get_data_buffer();
		in.p_data_fast = (in.Matrix)->get_data_buffer_fast();
		in.p_curr_data_laplace = (in.Matrix)->get_frame_laplace_fast();
		in.p_curr_data_laplace_normal = (in.Matrix)->get_frame_laplace();

		if(pHomeoFrame){
			in.p_homeo_data = (*pHomeoFrame)[i].get_data_buffer();
			in.p_fast_homeo_data = (*pHomeoFrame)[i].get_data_buffer_fast();
		}
		if(pLaplaceFrame){
			in.p_laplace_data = (*pLaplaceFrame)[i].get_data_buffer();
			in.p_laplace_data_fast = (*pLaplaceFrame)[i].get_data_buffer_fast();
		}

		// in.PrevMatrixPtr = ((Table2D<ELEM_TYPE>**)ccd_pipeline.GetAllMatrixPtrsChronological( i, FALSE ));
		in.PrevMatrixPtrCnt = ccd_pipeline.GetAllMatrixPtrsChronologicalInt( i, in,
									  TRUE, nFramesBack );
																					


		CPixelList pixel_list(in.xSize*in.ySize);
		CPixelList bright_list(in.xSize*in.ySize);
		out.pixel_list = &pixel_list;
		in.pPixelList = &pixel_list;
		
		long y_pos = (gCCDParams.m_nIgnoreEdgeBottom-1)*in.xSize+gCCDParams.m_nIgnoreEdgeLeft;
		for(in.y=gCCDParams.m_nIgnoreEdgeBottom;in.y<yUpperLimit;in.y++){
			y_pos += in.xSize;
			in.pos = y_pos;
			for(in.x=gCCDParams.m_nIgnoreEdgeLeft;in.x<xUpperLimit;in.x++,in.pos++){
				out.m_PixelOut.eventType	= eNone;



				if(!pixel_list.CheckPixel(in.pos)){						
				  BOOL_T bCheckPixel = AnalysePixel( in, out, TRUE);
				  
					if(gCCDParams.m_bSkipOverlaps){						
						clock_t t1_=clock();
						if(bCheckPixel && gCCDParams.m_bCheckForFlashes){
							// check overlaping with ealier detected event
							LONG_T prev_x,prev_y;
							if( (in.Matrix)->CheckIfOverlapsFast( in.x, in.y, out.cluster, out.cluster_cnt, prev_x, prev_y ) ){
								// printf("Event at (%d,%d) overlaps with (%d,%d)\n",in.x,in.y,prev_x,prev_y);
								MYTRACE4(gCCDTrace,"Event at (" << in.x << "," 
                           <<  in.y << ") overlaps with already detected event at (" << prev_x << "," << prev_y << ")" );
								bCheckPixel = FALSE;
							}
						}
						clock_t t2_=clock();
						total_in_check += (t2_ - t1_);
					}

					 // always checking sphericity :
   	         // LONG_T x0,y0;
	            // FindClusterAboveTresholdOpt2( in, x0,y0, out.cluster, out.cluster_cnt );
      	      // double max_redial = CCDFastCalc::FindMaxRedial( out.cluster, out.cluster_cnt, x0, y0, in.xSize );
         	   // out.m_PixelOut.m_Sphericity = CCDFastCalc::CalcSphericity( max_redial, out.cluster_cnt );

					if(gCCDParams.m_CheckEventShape>=0){
						// lastly checking shape of the event :
						(in.pPixelList)->Clear();
						if(!CheckEventShape( in, out ))
							bCheckPixel = FALSE; 					
					}
					

					if(bCheckPixel){
							MYTRACE2(gCCDTrace,"ACCEPTED : Potential event found x=" << out.m_PixelOut.x0 << ", y="<<out.m_PixelOut.y0
										<<" ccd_index="<<i<<", sum="<<out.m_PixelOut.newSum<<", max_sum="<<out.m_PixelOut.maxSum);

							n_found_events++;
							n_accepted_on_level_1++;



							(in.Matrix)->AddFoundEvent( in, out );
							MYTRACE2(gCCDTrace,"ACCEPTED : Potential " 
									<< GetEventTypeDesc(out.m_PixelOut.eventType) <<  " event found x=" 
									<< in.x << ", y="<<in.y<<" ccd_index="<<i<<", sum="
									<<out.m_PixelOut.newSum
									<<", max_sum="<<out.m_PixelOut.maxSum);
							(in.Matrix)->SetInteresting();
							newFrame.SetInteresting();
					}
#ifdef _MONTE_CARLO_
					else{			
						// posibility to turn off - in case no event detailed report
						// is nessesary for us - for big parameters anlysis can make it slow
						if(gCCDParams.m_bGenNotFoundReport && gCCDParams.GetMC()){
							// in simulation version - update generated event info
							// - why it was refused !
							(in.Matrix)->UpdateGenEvents( in, out );
						}
					}
#endif
				}
			}
		}			

		if(gCCDParams.m_bSkipIfMoreThen>0 && gCCDParams.m_eCheckIfMorePoint!=eAfterCoic){
	      // requires rejection of all events in case in the redial of gCCDParams.m
   	   // there is more then m_bSkipIfMoreThen events :
         (in.Matrix)->RejectIfMoreThen( &ccd_pipeline );
	   }
			


	}
	(in.Matrix)->UpdateGenEventsByArea( gCCDParams.m_nIgnoreEdgeLeft,xUpperLimit,	
													gCCDParams.m_nIgnoreEdgeBottom,yUpperLimit );
	PROFILER_END("Full analysis of new frame took : ")
	PROFILER_PUT("CheckIfOverlaps took : ",total_in_check);

	if(n_found_events)
		bRet = TRUE;
	return bRet;
}

void CCD_Analyser::CalcMaxLaplaceTab( const CPixelAnalyseIn& in )
{
	PROFILER_START
	Assert(in.pMaxLaplacePrev!=NULL,"Max prev laplace frame not initailized");

	BIG_ELEM_TYPE** p_prev_lap = (in.pMaxLaplacePrev)->get_data_buffer_fast();
	ELEM_TYPE** p_data_fast = (in.PrevMatrixPtr[0])->get_data_buffer_fast();
	long bottom = gCCDParams.GetEdgeSizeInLaplace();
	long upX = in.xSize-gCCDParams.GetEdgeSizeInLaplace();
	long upY = in.ySize-gCCDParams.GetEdgeSizeInLaplace();
	LONG_T plus_sum,minus_sum;
	register long x,y,v,iy,ix;
	register long max_x_redial = -1,max_y_redial=-1;

	// checking max shape redial :
	/*for(v=0;v<gCCDParams.m_nVetoPointsCount;v++){
      register long x = (gCCDParams.m_VetoArea)[v].x;
      register long y = (gCCDParams.m_VetoArea)[v].y;
	
		if(x>max_x_redial)
			max_x_redial = abs(x);
		if(y>max_y_redial)
      	max_y_redial = abs(y);
	}

	for(register long prev=0;prev<in.PrevMatrixPtrCnt;prev++){
		ELEM_TYPE** p_laplace = ((CCDMatrix*)in.PrevMatrixPtr[prev])->get_frame_laplace_fast();
		// long dxLong = my_round(in.PrevFramesShiftTab[prev].frameDX);
		// long dyLong = my_round(in.PrevFramesShiftTab[prev].frameDY);
		long dxLong = in.PrevFramesShiftTab[prev].frameDX;
		long dyLong = in.PrevFramesShiftTab[prev].frameDY;

	for(ix=bottom;ix<upX;ix++){
		for(iy=bottom;iy<upY;iy++){
			
				register long y_pos = iy+dyLong;
				register long x_pos = ix+dxLong;

				if(prev==0){
					p_prev_lap[iy][ix] = -10000;
					if(x_pos>=0 && y_pos>=0 && x_pos<in.xSize && y_pos<in.ySize){
						p_prev_lap[y_pos][x_pos] = p_laplace[y_pos][x_pos];						
					}else{
						p_prev_lap[y_pos][x_pos] = 0;
					}
				}
				for(register int v=0;v<gCCDParams.m_nVetoPointsCount;v++){					
	            register long x = x_pos+(gCCDParams.m_VetoArea)[v].x;
	            register long y = y_pos+(gCCDParams.m_VetoArea)[v].y;


					if(x>=0 && y>=0 && x<in.xSize && y<in.ySize){
						//if(ix==7 && iy==5)
						//	printf("(%d,%d)=%d\n",x,y,p_laplace[y][x]);
							
						if(p_laplace[y][x]>p_prev_lap[iy][ix]){
								p_prev_lap[iy][ix] = p_laplace[y][x];
						}
					}
				}
	
				//if(ix==7 && iy==5)
				//	printf("2 MAX = (%d,%d)=%d\n",ix,iy,p_prev_lap[iy][ix]);
			}
		}	
	}*/

	
	for(register long prev=0;prev<in.PrevMatrixPtrCnt;prev++){
		BIG_ELEM_TYPE** p_laplace = ((CCDMatrix*)in.PrevMatrixPtr[prev])->get_frame_laplace_fast();
		// long dxLong = my_round(in.PrevFramesShiftTab[prev].frameDX);
		// long dyLong = my_round(in.PrevFramesShiftTab[prev].frameDY);
		long dxLong = (long)in.PrevFramesShiftTab[prev].frameDX;
		long dyLong = (long)in.PrevFramesShiftTab[prev].frameDY;

		for(ix=bottom;ix<upX;ix++){
			register long x_pos = ix+dxLong;
			for(iy=bottom;iy<upY;iy++){
				// if(ix==7 && iy==5)
            //   printf("1 MAX = (%d,%d)=%d\n",ix,iy,p_prev_lap[iy][ix]);

				register long y_pos = iy+dyLong;

				if(prev==0){
					p_prev_lap[iy][ix] = -10000;
					if(x_pos>=0 && y_pos>=0 && x_pos<in.xSize && y_pos<in.ySize){
						p_prev_lap[y_pos][x_pos] = p_laplace[y_pos][x_pos];						
					}else{
						p_prev_lap[y_pos][x_pos] = 0;
					}
				}
				for(register int v=0;v<gCCDParams.m_nVetoPointsCount;v++){					
	            register long x = x_pos+(gCCDParams.m_VetoArea)[v].x;
	            register long y = y_pos+(gCCDParams.m_VetoArea)[v].y;


					if(x>=0 && y>=0 && x<in.xSize && y<in.ySize){
						//if(ix==7 && iy==5)
						//	printf("(%d,%d)=%d\n",x,y,p_laplace[y][x]);
							
						if(p_laplace[y][x]>p_prev_lap[iy][ix]){
								p_prev_lap[iy][ix] = p_laplace[y][x];
						}
					}
				}
	
				//if(ix==7 && iy==5)
				//	printf("2 MAX = (%d,%d)=%d\n",ix,iy,p_prev_lap[iy][ix]);
			}
		}	
	}
	PROFILER_END("Calculation of laplace MAX frame from prev took:")

/*	return;	

	// 
	long MaxLaplace = -10000;
	for(v=0;v<gCCDParams.m_nVetoPointsCount;v++){
      register long x = bottom+(gCCDParams.m_VetoArea)[v].x;
      register long y = bottom+(gCCDParams.m_VetoArea)[v].y;
	
      if(MaxLaplace<p_laplace[y][x]){
   	   MaxLaplace = p_laplace[y][x];
     	}
	}
	p_prev_lap[bottom][bottom] = MaxLaplace;
	p_prev_lap[bottom-1][bottom] = MaxLaplace;

	for(iy=bottom;iy<upY;iy++){
		// change of y check new level :
		register long x,y;
		//for(x=bottom-redial;x<=bottom+redial;x++){
		//	for(y=iy-redial;y<=iy+redial;y++){
		//		if(p_laplace[y][x]>p_prev_lap[iy][bottom])
		//			p_prev_lap[iy][bottom] = p_laplace[y][x];
		//	}
		//}
			

		//for(x=bottom-redial;x<=bottom+redial;x++){
		//	if(p_laplace[iy+redial][x]>p_prev_lap[iy][bottom])
		//		p_prev_lap[iy][bottom] = p_laplace[iy+redial][x];	
		//}

		//for(ix=bottom+1;ix<upX;ix++){
		//	p_prev_lap[iy][ix] = MAX(p_prev_lap[iy][ix-1],p_prev_lap[iy][ix]);
		//	for(y=bottom-redial;y<=bottom+redial;y++){
		//		if( p_laplace[y][ix+redial] > p_prev_lap[iy][ix] )
		//			p_prev_lap[iy][ix] = p_laplace[y][ix+redial];
		//	}			
		//}
		for(ix=bottom;ix<upX;ix++){
			for(x=ix-redial;x<=ix+redial;x++){
       	  for(y=iy-redial;y<=iy+redial;y++){
            if(p_laplace[y][x]>p_prev_lap[iy][ix])
               p_prev_lap[iy][ix] = p_laplace[y][x];
	        }
   	   }
		}
	}
	PROFILER_END("Calculation of laplace MAX frame from prev took:")*/
}

BOOL_T CCD_Analyser::StartOptimizedAnalyse(CCDPipeline& ccd_pipeline)
{
	if(gCCDParams.m_bCheckGenOnly){
		return CheckGenerated( ccd_pipeline );
	}else{
		if(gCCDParams.m_nSkipNFramesAfterChange>0){
			gCCDParams.m_nSkipNFramesAfterChange--;
			return FALSE;	
		}

		if( gCCDParams.m_bReadFirstLevelInfoFromLog ){
			printf("Using FLT info from old log files\n");
			return AnalyseNewFrameWithAverageOfPrevAlg_FromLogFiles( ccd_pipeline );
		}else{
			if(gCCDParams.m_bAverageOfPrev){			
				// test :
				// if( gCCDParams.GetMC() ){
				// 	return AnalyseNewFrameWithAverageOfPrevAlgMC( ccd_pipeline );				
				// }else{
				return AnalyseNewFrameWithAverageOfPrevAlg( ccd_pipeline );
				// }
			}
		}
	}
	Assert(FALSE,"No algorithm specified");
	return FALSE;
}


BOOL_T CCD_Analyser::CheckGenerated( CCDPipeline& ccd_pipeline )
{
	Initialize();
	LONG_T pipeline_size = ccd_pipeline.GetPipelineSize();
	LONG_T n_found_events = 0;
	LONG_T n_accepted_on_level_1=0;
	LONG_T n_accepted_on_level_2=0;

	if(ccd_pipeline.GetCount()!=pipeline_size)
		return FALSE;

	BOOL_T bRet = FALSE;

	PROFILER_START
	cCCD& newFrame = ccd_pipeline.GetCurrent();
	long size = newFrame.GetCount();
	if(size<1){
		return FALSE;
	}
	CPixelAnalyseIn in;
   CPixelAnalyseOut out;

	in.PrevMatrixPtrCnt = 0;
	in.frame_index = ccd_pipeline.GetFrameIndex();
	in.xSize = newFrame[0].GetXSize();
   in.ySize = newFrame[0].GetYSize();	
	register int yUpperLimit = (in.ySize-gCCDParams.m_nIgnoreEdgeUp);
	register int xUpperLimit = (in.xSize-gCCDParams.m_nIgnoreEdgeRight);		
	m_FirstDy = ((m_NeighbPoints[0].y)*in.xSize);
	in.treshold_for_max = (long)(gCCDParams.m_nCalcMaxForAboveNSigma*CCDDataResults::GetSigmaBackground( gCCDParams.m_eLaplaceType ));
	in.treshold_for_not_more_then_n_exceeds = (long)(gCCDParams.m_TreshForCheckIfNoneOfNExceedsInSigma*CCDDataResults::GetSigmaBackground( gCCDParams.m_eLaplaceType ));
	out.m_PixelOut.treshold_for_not_more_then_n_exceeds = in.treshold_for_not_more_then_n_exceeds;

	// printf("TRESHOLD = %d\n",(long)in.treshold_for_not_more_then_n_exceeds);

	in.pPipeline = &ccd_pipeline;
	in.pipeline_size_minus_1 = pipeline_size-1;


	CManyTab2D<BIG_ELEM_TYPE>* pLaplaceFrame = ccd_pipeline.GetLaplaceFrame();
	Assert(pLaplaceFrame!=NULL,"No homeopatic frame of laplace !");

	clock_t total_in_check=0;

	register int nNewLaplaceTreshold = ((in.pPipeline)->GetPipelineCfg()).m_nNewLaplace;
	register int nMaxLaplaceOnOther = ((in.pPipeline)->GetPipelineCfg()).m_nMaxLaplaceOnOther;

	for(register int i=0;i<size;i++){
		in.ccd_index = i;
		in.pCCDInfo = &(ccd_pipeline.GetCCDInfoTab()[in.ccd_index]);
		in.Matrix = &(newFrame[i]);
		in.p_data = (in.Matrix)->get_data_buffer();
		in.p_data_fast = (in.Matrix)->get_data_buffer_fast();
		in.p_curr_data_laplace = (in.Matrix)->get_frame_laplace_fast();
		in.p_curr_data_laplace_normal = (in.Matrix)->get_frame_laplace();			

		in.p_laplace_data = (*pLaplaceFrame)[i].get_data_buffer();
		in.p_laplace_data_fast = (*pLaplaceFrame)[i].get_data_buffer_fast();



		register int bOK=0;
		register int y_pos = (gCCDParams.m_nIgnoreEdgeBottom-1)*in.xSize+gCCDParams.m_nIgnoreEdgeLeft;
		register int pos = 0;
		register int low_y = gCCDParams.m_nIgnoreEdgeBottom;
		register int low_x = gCCDParams.m_nIgnoreEdgeLeft;
		register LONG_T max_pos=-1;
		for(register int y=low_y;y<yUpperLimit;y++){
			y_pos += in.xSize;
			pos = y_pos;
			for(register int x=low_x;x<xUpperLimit;x++,pos++){
				if(!(in.Matrix)->IsGenerated(x,y)){
					continue;
				}


				if(in.p_curr_data_laplace[y][x]<=nNewLaplaceTreshold){
					continue;
				}

				in.x = x;
            in.y = y;
            in.pos = pos;

				/*if(in.p_data_fast[y][x]>=gCCDParams.m_MaxAllowedVal){
					continue;
				}

				bOK = 1;
				for(register int v=0;v<gCCDParams.m_nVetoPointsCount;v++){
					register long xx = x+(gCCDParams.m_VetoArea)[v].x;
					register long yy = y+(gCCDParams.m_VetoArea)[v].y;

					if(in.p_laplace_data_fast[yy][xx]>=nMaxLaplaceOnOther ||
						in.p_laplace_data_fast[yy][xx]<=gCCDParams.m_nMinLaplaceOnOther){
						bOK = 0;
					   break;
					}
				}
				if(!bOK){
					continue;
				}

         	if(!CheckLocalMaxCondition(in, out, max_pos)){
					continue;
				}

				if(gCCDParams.m_bSkipOverlaps){						
					register LONG_T prev_x,prev_y;
					if( (in.Matrix)->CheckIfOverlapsFast( x, y, out.cluster, 0 , prev_x, prev_y ) ){
						continue;
					}
				}*/

				/*if(!AnalysePixel(in,out)){
					printf("ERROR IN CODE !!!! at (%d,%d)\n",x,y);
					BOOL_T bCheckAgain = AnalysePixel(in,out);
					Assert(FALSE,"ERROR IN CODE !!!! at (%d,%d)\n",x,y);
				}*/

				//if(!bAnalPixel)
				//	Assert(FALSE,"error at (%d,%d) !",x,y);

				out.m_PixelOut.laplaceSum = in.p_curr_data_laplace[y][x];
				// out.m_PixelOut.prevLaplaceSum =  
				out.m_PixelOut.eventType = eFlash;
				(in.Matrix)->AddFoundEvent( in, out );
				n_found_events++;
			}
		}
	}
	PROFILER_END("Full analysis of new frame took : ")
	PROFILER_PUT("CheckIfOverlaps took : ",total_in_check);

	if(n_found_events)
		bRet = TRUE;
	return bRet;

}


void CCD_Analyser::CalcTresholdsByMap( CPixelAnalyseIn& in, int& nNewLaplaceTreshold,
                   				         int& nMaxLaplaceOnOther )
{
	if( (in.pPipeline)->GetBackgroundMap() ){
		Area2DInfo& map_elem = (in.pPipeline)->GetBackgroundMap()->GetAreaDesc( in.x, in.y );
		
		nNewLaplaceTreshold = (int)(gCCDParams.m_nNewLaplaceInSigma*(map_elem.m_DataInfo[ gCCDParams.m_eLaplaceType ].m_Sigma )
			+ map_elem.m_DataInfo[ gCCDParams.m_eLaplaceType ].m_Average);		
		nMaxLaplaceOnOther = (int)(gCCDParams.m_nMaxLaplaceOnOtherInSigma*(map_elem.m_DataInfo[ gCCDParams.m_eLaplaceType ].m_Sigma )
			+ map_elem.m_DataInfo[ gCCDParams.m_eLaplaceType ].m_Average);
		in.treshold_for_cluster = (int)(gCCDParams.m_ClusterIfNSigmaAboveBackgr*(map_elem.m_DataInfo[ gCCDParams.m_eLaplaceType ].m_Sigma )
			+ map_elem.m_DataInfo[ gCCDParams.m_eLaplaceType ].m_Average);
	}
}

void CCD_Analyser::CalculateTresholds( CPixelAnalyseIn& in )
{
	// treshold for raw cluster :
	in.treshold_for_raw_cluster = (int)CCDDataResults::GetTreshold( eRawS, gCCDParams.m_nSigmaAboveMeanInRawCluster, (in.pPipeline)->GetBackgroundMap() );
	in.treshold_for_raw_cluster_aver = (int)CCDDataResults::GetTresholdAver( eRawS, gCCDParams.m_nSigmaAboveMeanInRawCluster, (in.pPipeline)->GetBackgroundMap() );
	in.treshold_for_cluster = GetMaxNoiseLevelLaplace( &in );

	in.nAver =  gCCDParams.m_nMaxOfAverageOfPrevN;
	double sigmaB,meanB,meanAver;
	CCDPipeline::GetBackgroundStat( in.pPipeline, in.ccd_index, meanB, sigmaB, gCCDParams.m_eLaplaceType );
	in.SigmaLap = sigmaB;
	in.MeanLap = meanB;
	double sigAver = sigmaB/sqrt((double)in.nAver);
	in.treshold_for_prev = meanAver + sigAver*gCCDParams.m_ClusterIfNSigmaAboveBackgr;
	in.treshold_for_hot = (int)(gCCDParams.m_nRejectHotPixelsTresholdInSigma*sigmaB);

	for(register int i=0;i<(int)eLastEnumElem;i++){
		CCDPipeline::GetBackgroundStat( in.pPipeline, in.ccd_index,
												  in.m_pDistrValues[i].mean, in.m_pDistrValues[i].sigma,
                                      (eLaplaceType_T)i );
	}

	// treshold for black pixel - pixels is black if is below this limit :
	in.treshold_for_black_pixel = in.m_pDistrValues[(int)eRawS].mean - gCCDParams.m_bBlackPixelsIfNSigmaBelow*in.m_pDistrValues[(int)eRawS].sigma;

	in.treshold_for_max = (long)(gCCDParams.m_nCalcMaxForAboveNSigma*sigmaB);
	in.treshold_for_not_more_then_n_exceeds = (long)(gCCDParams.m_TreshForCheckIfNoneOfNExceedsInSigma*sigmaB);
}

BOOL_T CCD_Analyser::CheckIfHotPixel( CPixelAnalyseIn& in, CPointList& hotList )
{
	CPointList::iterator i;
	for(i=hotList.begin();i!=hotList.end();i++){
		if( abs((int)i->x-in.x)<=2 && abs((int)i->y-in.y)<=2 )
			return TRUE;
	}
	return FALSE;
}

BOOL_T CCD_Analyser::CheckBlackPixels( int* plusValues, int* minusValues,
												  int plusValCount,int minusValCount, 
												  const CPixelAnalyseIn& in,
												  CPixelAnalyseOut& out )
{	
	/*for(register int i=0;i<minusValCount;i++){
		if(minusValues[i]<=0 || minusValues[i]<=in.treshold_for_black_pixel)
			return FALSE;
	}
	
	double two_sig=2*in.m_pDistrValues[(int)eRawS].sigma;
	int min_pos=-1;
	int min_of_minus = find_min_value( minusValues, minusValCount, min_pos );

	for(register int i=0;i<minusValCount;i++){
		if( i!=min_pos ){
			if( min_of_minus >= (minusValues[i]-two_sig) ){
				return TRUE;
			}
		}
	}*/

	int min_pos=-1;
	int min_of_minus = find_min_value( minusValues, minusValCount, min_pos );
	double sum=0.00;
	if( min_pos>=0 ){
		for(register int i=0;i<minusValCount;i++){	
			if( i!=min_pos ){
				sum += minusValues[i];
			}
		}
		sum = (sum/(minusValCount-1));
		out.m_PixelOut.m_fBlackRatio = double(min_of_minus)/sum;
		if(out.m_PixelOut.m_fBlackRatio<gCCDParams.m_fBlackPixelsRatio)
			return FALSE;
	}


	// in case it is OK - first TRUE is returned - otherwise 
	// minimum is smaller then others-2Sigma
	return TRUE;
}

double CCD_Analyser::GetBlackRatio( int* plusValues, int* minusValues,
 											   int plusValCount,int minusValCount )
{	
	double ret=1.00;
	int min_pos=-1;
	int min_of_minus = find_min_value( minusValues, minusValCount, min_pos );
	double sum=0.00;
	if( min_pos>=0 ){
		for(register int i=0;i<minusValCount;i++){	
			if( i!=min_pos ){
				sum += minusValues[i];
			}
		}
		sum = (sum/(minusValCount-1));
		ret = double(min_of_minus)/sum;
	}
	return ret;
}


BOOL_T CCD_Analyser::ApplyCuts( int x, int y, int pos,
										  CPixelAnalyseIn& in, CPixelAnalyseOut& out,
										  CFrameIdentStat& stat, CFrameIdentStat& sstat,
										  BOOL_T bIsGenObject, BOOL_T& bDoBreak,
										  CPixelList& allclusters )
{  
#ifndef _OPTIMIZED_VERSION_
	if( !bIsGenObject ){
		stat.nTprevCut++;
	}else{
		sstat.nTprevCut=1;
	}
#endif

	if( gCCDParams.m_bMinLaplaceOnOtherEnabled ){
		if( out.m_PixelOut.maxAverageOfPrev < gCCDParams.m_nMinLaplaceOnOther ){
			return FALSE;
		}
						
	}
#ifndef _OPTIMIZED_VERSION_
	if( !bIsGenObject ){
		if( out.m_PixelOut.eventType!=eBrighten ){
			stat.nMinPrevCut++;
		}
	}else{
		sstat.nMinPrevCut=1;
	}
#endif

	// checking rates after tv also :
	if( gCCDParams.m_bSkipIfMoreThen>0 && gCCDParams.m_eCheckIfMorePoint!=eAfterCoic ){
			if( gCCDParams.m_MaxNumberOfEventsOnFrameAfterTv>0 ){
				if( stat.nMinPrevCut > gCCDParams.m_MaxNumberOfEventsOnFrameAfterTv )
				{
					printf("REJECTING ALL DUE TO AFTER_TV CUT %d > %d\n",stat.nMinPrevCut,
								gCCDParams.m_MaxNumberOfEventsOnFrameAfterTv);
					(in.Matrix)->GetFoundEvents().clear();
                stat.bAllRejected=TRUE;
                sstat.bAllRejected=TRUE;
					 bDoBreak=TRUE;
					 return FALSE;
				}
			}
	}


	
	if( out.m_PixelOut.eventType!=eBrighten && gCCDParams.m_bLocalMaxReq){
		LONG_T max_pos=-1;
		if(!CheckLocalMaxCondition(in, out, max_pos)){
			return FALSE;
		}
	}
#ifndef _OPTIMIZED_VERSION_
	if( !bIsGenObject ){
		stat.nLocalMaxCut++;
	}else{
		sstat.nLocalMaxCut=1;
	}
#endif

				
/*				if( m_pAnalFoundEventFunc ){
					// sum customization analysisis - now cosmic analysis will be performed 
					// in this way :
					BOOL_T bEventAnal = (*m_pAnalFoundEventFunc)( *(in.Matrix), (long)in.x, (long)in.y, NULL );
					if(!bEventAnal)
						continue;
				}*/


	BOOL_T bClusterFound=FALSE;
	if(gCCDParams.m_bSkipOverlaps){						
		register LONG_T prev_x,prev_y;
		if( (in.Matrix)->CheckIfOverlapsFast( x, y, out.cluster, 0 , prev_x, prev_y ) ){
			return FALSE;
		}
	}
#ifndef _OPTIMIZED_VERSION_
	if( !bIsGenObject ){
		stat.nOverlap++;
	}else{
		sstat.nOverlap=1;
	}
#endif

	// always checking sphericity : 
	if( gCCDParams.m_CheckEventShape>=0){
		if( allclusters.CheckPixel( in.pos )){
			// pixels in already rejected cluster - skip
			_TRACE_PRINTF_5("IN ALREADY REJECTED SHAPE spher=%f, point_in_cluster=%d\n",out.m_PixelOut.m_Sphericity,out.cluster_cnt);				
			return FALSE;
		}										
	
		double x0,y0;
		double max_noise_level = in.treshold_for_cluster;// WAS 0
		//(in.pPixelList)->Clear(); // very slow !!!!!!!!!!!!!!!!! OPT - required !
		FindClusterAboveTresholdOpt3( in, out, x0,y0, max_noise_level );
		bClusterFound=TRUE;
			
   	double max_redial = CCDFastCalc::FindMaxRedial( out.cluster, out.cluster_cnt, x0, y0, in.xSize );
   	out.m_PixelOut.m_Sphericity = CCDFastCalc::CalcSphericity( max_redial, out.cluster_cnt );

		_TRACE_PRINTF_5("CLUSTER FOUND (%d,%d): spher=%f, point_in_cluster=%d\n",in.x,in.y,out.m_PixelOut.m_Sphericity,out.cluster_cnt);
		/*if(out.m_PixelOut.m_Sphericity < gCCDParams.m_CheckEventShape || 
			out.cluster_cnt > gCCDParams.m_MaxPixelsInClusterAllowed)
		{
			_TRACE_PRINTF_4("at (%d,%d) SHAPE REJECTED spher=%f, point_in_cluster=%d\n",x,y,out.m_PixelOut.m_Sphericity,out.cluster_cnt);
			for(register int c=0;c<out.cluster_cnt;c++){
				// marking rejected pixels :
				allclusters.HitPixel( out.cluster[c] );
			}
			return FALSE;
		}*/
	}
#ifndef _OPTIMIZED_VERSION_
	if( !bIsGenObject ){
		stat.nShapeCut++;
	}else{
		sstat.nShapeCut=1;
	}
#endif

	// Black pixels rejection :
	Table2D<ELEM_TYPE>::GetLaplacePlusMinusValues( in.p_data_fast, 
		out.m_PixelOut.laplacePlusValues, out.m_PixelOut.laplaceMinusValues,
		out.m_PixelOut.laplacePlusCount, out.m_PixelOut.laplaceMinusCount,
		x,y,in.xSize,in.ySize,gCCDParams.m_eLaplaceType );
		
	if( gCCDParams.m_bRejectBlackPixels ){
		// checking if no of minuses was below minmum limit :
		if(!CheckBlackPixels( out.m_PixelOut.laplacePlusValues, out.m_PixelOut.laplaceMinusValues,
								 	 out.m_PixelOut.laplacePlusCount, out.m_PixelOut.laplaceMinusCount, 
									in, out )){
			// returns FALSE if any of pixels used is black :
			return FALSE;	
		}
	}
#ifndef _OPTIMIZED_VERSION_
	if( !bIsGenObject ){
		stat.nBlackCut++;
	}else{
		sstat.nBlackCut=1;
	}
#endif

	// hot pixels rejection based on list :
	if( (in.pPipeline->m_PipelineCfg).m_bRejectHotByList ){
		if(CheckIfHotPixel( in, (in.pPipeline->m_PipelineCfg).m_HotPixels ) ){
			// it is hot pixel reject it :
			return FALSE;
		}
	}

	if( gCCDParams.m_bRejectHotPixelsByAverage ){
		BOOL_T bOK=TRUE;
		out.m_PixelOut.m_PrevLapInPixel = -100000;
		for(register int yy=(y-1);yy<=(y+1) && bOK;yy++){
			for(register int xx=(x-1);xx<=(x+1) && bOK;xx++){
				int new_val = CalcAverageOfPrevNInPixel( xx, yy, in, 
										gCCDParams.m_nMaxOfAverageOfPrevN,
										TRUE );
				if(new_val>out.m_PixelOut.m_PrevLapInPixel){
					out.m_PixelOut.m_PrevLapInPixel = new_val;
				}
			}
		}	
		if( out.m_PixelOut.m_PrevLapInPixel > in.treshold_for_hot ){
			// average from previous frames on this pixel is 
			// high - rejected as hot pixel :
			bOK = FALSE;
		}		
		if( !bOK ){
			return FALSE;
		}
	}
#ifndef _OPTIMIZED_VERSION_
	if( !bIsGenObject ){
		stat.nHotCut++;
	}else{
		sstat.nHotCut=1;
	}
#endif


	out.m_PixelOut.m_MaxClusterValue = in.p_curr_data_laplace[y][x];
	out.m_PixelOut.m_MaxClusterPixel.x = x;
   out.m_PixelOut.m_MaxClusterPixel.y = y;
	out.m_PixelOut.MaxPixelRawValue = in.p_data_fast[y][x];
	int maxPos = pos;
	if( bClusterFound ){
		// choose bigest value in cluster and set :
		for(register int k=0;k<out.cluster_cnt;k++){
			if( in.p_curr_data_laplace_normal[ out.cluster[k] ]>out.m_PixelOut.m_MaxClusterValue ){
				out.m_PixelOut.m_MaxClusterValue = in.p_curr_data_laplace_normal[ out.cluster[k] ];
				maxPos = out.cluster[k];
			}
		}
		if( maxPos!=pos ){
			out.m_PixelOut.m_MaxClusterPixel.x = ( maxPos % in.xSize );
			out.m_PixelOut.m_MaxClusterPixel.y = ( maxPos / in.xSize );
			out.m_PixelOut.MaxPixelRawValue = in.p_data_fast[(int)out.m_PixelOut.m_MaxClusterPixel.y][(int)out.m_PixelOut.m_MaxClusterPixel.x];
		}
	}
	out.m_PixelOut.laplaceSum = in.p_curr_data_laplace[y][x]; 

	// now check MaxPixel - overlap :
	/*if(gCCDParams.m_bSkipOverlaps && maxPos!=pos ){		
		register LONG_T prev_x,prev_y;
		if( (in.Matrix)->CheckIfOverlapsFastMaxPoint( (in.Matrix)->GetFoundEvents() , (int)out.m_PixelOut.m_MaxClusterPixel.x, 
					(int)out.m_PixelOut.m_MaxClusterPixel.y, out.cluster, 0 , prev_x, prev_y ) ){
			return FALSE;
		}
	}*/

			
	if(gCCDParams.GetMC()){
		if( (in.pPipeline)->CheckIfReGenEvent( x, y ) ){
			return FALSE;
		}
	}

	return TRUE;	
}

BOOL_T CCD_Analyser::HistoVariables( int x, int y, int pos,
										  CPixelAnalyseIn& in, CPixelAnalyseOut& out,
										  CFrameIdentStat& stat, CFrameIdentStat& sstat,
										  CVariableInfo& eventVariables, double SigmaB )
{  

	eventVariables.l_new = (in.p_curr_data_laplace[y][x]/SigmaB);
	eventVariables.l_prev = (out.m_PixelOut.maxAverageOfPrev/SigmaB);

	register LONG_T prev_x,prev_y;
	eventVariables.nOverlaps = (in.Matrix)->CountOverlaps( x, y, out.cluster, 0 , prev_x, prev_y );
	
	// SPHERICITY :
	double x0,y0;
	double max_noise_level = in.treshold_for_cluster;// WAS 0
	FindClusterAboveTresholdOpt3( in, out, x0,y0, max_noise_level );
   double max_redial = CCDFastCalc::FindMaxRedial( out.cluster, out.cluster_cnt, x0, y0, in.xSize );
   eventVariables.fSpericity = CCDFastCalc::CalcSphericity( max_redial, out.cluster_cnt );


	// Black pixels ratio :
	Table2D<ELEM_TYPE>::GetLaplacePlusMinusValues( in.p_data_fast, 
		out.m_PixelOut.laplacePlusValues, out.m_PixelOut.laplaceMinusValues,
		out.m_PixelOut.laplacePlusCount, out.m_PixelOut.laplaceMinusCount,
		x,y,in.xSize,in.ySize,gCCDParams.m_eLaplaceType );		
	CheckBlackPixels( out.m_PixelOut.laplacePlusValues, out.m_PixelOut.laplaceMinusValues,
						 	 out.m_PixelOut.laplacePlusCount, out.m_PixelOut.laplaceMinusCount, 
							in, out );
	eventVariables.fBlackRatio = out.m_PixelOut.m_fBlackRatio;

	// Hot pixels :
	BOOL_T bOK=TRUE;
	out.m_PixelOut.m_PrevLapInPixel = -100000;
	for(register int yy=(y-1);yy<=(y+1) && bOK;yy++){
		for(register int xx=(x-1);xx<=(x+1) && bOK;xx++){
			int new_val = CalcAverageOfPrevNInPixel( xx, yy, in, 
									gCCDParams.m_nMaxOfAverageOfPrevN,
									TRUE );
			if(new_val>out.m_PixelOut.m_PrevLapInPixel){
				out.m_PixelOut.m_PrevLapInPixel = new_val;
			}
		}
	}	
	eventVariables.averInPixel = (out.m_PixelOut.m_PrevLapInPixel/SigmaB);
	
	if( m_pAnalFoundEventFunc ){
		(*m_pAnalFoundEventFunc)( &varInfo, eCVariableInfo, (in.pPipeline)->GetPipelineIndex(), NULL );
	}
}


BOOL_T CCD_Analyser::AnalyseNewFrameWithAverageOfPrevAlg_FromLogFiles( CCDPipeline& ccd_pipeline )
{
	Initialize();
	BOOL_T bMC = gCCDParams.GetMC();
	LONG_T pipeline_size = ccd_pipeline.GetPipelineSize();
	int pipeline_index = ccd_pipeline.GetPipelineIndex();
	int nDayFrame = ccd_pipeline.GetDayFrameCounter();
	register LONG_T max_pos=-1;


	PROFILER_START

	int bContinue=1;

	// preparing statistics objects :
	CFrameIdentStat stat( nDayFrame ), sstat( nDayFrame );
		
	if(ccd_pipeline.GetCount()!=pipeline_size)
      return FALSE;

   cCCD& newFrame = ccd_pipeline.GetCurrent();
	CCDMatrix* current_image = &(newFrame[0]);


	BOOL_T bDoBreak=FALSE;


	int nEvents=0;

	printf("selecting events for current frame from list ...");fflush(0);
	CCDEventList& allEvents = ccd_pipeline.GetEventFromLogFile();
	if( allEvents.size() ){
		for(int i=0;i<allEvents.size();i++){
			CccdReport& evt = allEvents[i];
			
			if( evt.m_DayFrameIndex == nDayFrame ){
				current_image->AddFoundEvent( evt );
				nEvents++;
			}else{
				if( evt.m_DayFrameIndex > nDayFrame ){
					break;
				}
			}
		}
	}
	printf("selection of events - DONE\n");
	

	// log event rates :
	if(gCCDParams.m_bLogFrameStat && !gCCDParams.m_bCCDDouble){
		LogEventRates( &ccd_pipeline );
	}
		

	PROFILER_END("Full analysis of new frame took : ")

	printf("For Frame %d - read %d events from log file\n",nDayFrame,nEvents);

	return ( nEvents>0 );
}


// filling samples info, the following fields are filled here :
// evt.m_PixelAnalResults.maxAverageOfPrev 
// evt.m_PixelAnalResults.PixelRawValue    - maximum of raw value
// evt.m_PixelAnalResults.laplaceSum       - maximum of laplace value
int CCD_Analyser::FillSamplesInfo( CPixelAnalyseIn& in, int nNewLaplaceTreshold )
{
	if( in.pPipeline ){
		if( (in.pPipeline)->m_allGenEvents.size() && (in.pPipeline)->m_allGenEvents.back().size()){
			CCDEventList& genListUpdt = (in.Matrix)->GetGenEvents();

			for(int i=0;i<(((in.pPipeline)->m_allGenEvents).back()[0]).size();i++){
				// for each generated event fill : max_value, max_laplace value etc ...
				CccdReport& evt = (((in.pPipeline)->m_allGenEvents).back()[0])[i];
				printf("GEN_EVENT (%d,%d)\n",(int)evt.m_MaxPoint.x,(int)evt.m_MaxPoint.y);


				int max_raw=-10000;
				int max_laplace=-10000;
				int evt_x = (int)evt.m_MaxPoint.x;
				int evt_y = (int)evt.m_MaxPoint.y;
				for( int y=evt_y-5;y<=evt_y+5;y++){
					for( int x=evt_x-5;x<=evt_x+5;x++){
						if( y>=0 && x>=0 && y<in.ySize && x<in.xSize ){
							if( in.p_data_fast[y][x] > max_raw ){
								max_raw = in.p_data_fast[y][x];
							}
							
							
							if( in.p_curr_data_laplace[y][x] > max_laplace ){
								max_laplace = in.p_curr_data_laplace[y][x];
							}
						}
					}
				}
				evt.m_PixelAnalResults.laplaceSum = max_laplace;
				evt.m_PixelAnalResults.PixelRawValue = max_raw;				

				int maxPrev=-10000;
				for( int y=evt_y-gCCDParams.m_bGenEventRedial;y<=evt_y+gCCDParams.m_bGenEventRedial;y++){
					for( int x=evt_x-gCCDParams.m_bGenEventRedial;x<=evt_x+gCCDParams.m_bGenEventRedial;x++){
						if( in.p_curr_data_laplace[y][x] > nNewLaplaceTreshold ){
	
							//  out.m_PixelOut.maxAverageOfPrev = GetMaxAverageInVetoArea( in, TRUE );
							// out.m_PixelOut.maxAverageOfPrev = -10000;
							for(register int v=0;v<gCCDParams.m_nVetoPointsCount;v++){
						       register int xx = (long)(x+(gCCDParams.m_VetoArea)[v].x);
						       register int yy = (long)(y+(gCCDParams.m_VetoArea)[v].y);
	
								 register int prevAverage = CalcAverageOfPrevN( xx, yy, in,  
																		 gCCDParams.m_nMaxOfAverageOfPrevN,
																		 TRUE );
								 if( prevAverage>maxPrev ){
						   		  maxPrev = prevAverage;
								 }
							}		
						}
					}
				}
				evt.m_PixelAnalResults.maxAverageOfPrev = maxPrev;


				CccdReport* pEvt = genListUpdt.FindEvent( (int)evt.m_Point.x, (int)evt.m_Point.y );
				if( pEvt ){
					(pEvt->m_PixelAnalResults).laplaceSum = max_laplace;
					(pEvt->m_PixelAnalResults).PixelRawValue = max_raw;
					(pEvt->m_PixelAnalResults).maxAverageOfPrev = maxPrev;
				}
			}			

			return (((in.pPipeline)->m_allGenEvents).back()[0]).size();
		}
	}
	return 0;	
}



BOOL_T CCD_Analyser::AnalyseNewFrameWithAverageOfPrevAlg(CCDPipeline& ccd_pipeline)
{	
	Initialize();
	BOOL_T bMC = gCCDParams.GetMC();
	LONG_T pipeline_size = ccd_pipeline.GetPipelineSize();
	int pipeline_index = ccd_pipeline.GetPipelineIndex();
	LONG_T ClusterWithMoreCnt=0;
	LONG_T n_found_events = 0;
	LONG_T n_accepted_on_level_1=0;
	LONG_T n_accepted_on_level_2=0;
   LONG_T star_cluster[MAX_CLUSTER_SIZE];
	LONG_T nFramesBack = gCCDParams.m_FramesBack;
	int nDayFrame = ccd_pipeline.GetDayFrameCounter();
	register LONG_T max_pos=-1;

	if(gCCDParams.m_nMaxOfAverageOfPrevN>nFramesBack)
		nFramesBack = gCCDParams.m_nMaxOfAverageOfPrevN;
	if(gCCDParams.m_nNotMoreThenNExceedsBackFramesCount>nFramesBack)
		nFramesBack = gCCDParams.m_nNotMoreThenNExceedsBackFramesCount;

	// to have also current frame
	nFramesBack++;			


	if(ccd_pipeline.GetCount()!=pipeline_size)
		return FALSE;

	BOOL_T bRet = FALSE;

	PROFILER_START
	cCCD& newFrame = ccd_pipeline.GetCurrent();
	long size = newFrame.GetCount();
	if(size<1){
		return FALSE;
	}
	CPixelAnalyseIn in;
   CPixelAnalyseOut out;

	CPixelAnalyseIn* inForPrev = new CPixelAnalyseIn();
	CPixelAnalyseOut* outForPrev = new CPixelAnalyseOut();

	in.PrevMatrixPtrCnt = nFramesBack;
	in.frame_index = ccd_pipeline.GetFrameIndex();
	in.xSize = newFrame[0].GetXSize();
   in.ySize = newFrame[0].GetYSize();	
	int yUpperLimit = (in.ySize-gCCDParams.m_nIgnoreEdgeUp);
	int xUpperLimit = (in.xSize-gCCDParams.m_nIgnoreEdgeRight);		
	m_FirstDy = ((m_NeighbPoints[0].y)*in.xSize);


	CPixelList pixel_list(in.xSize*in.ySize),allclusters(in.xSize*in.ySize);
	CPixelList star_clusters(in.xSize*in.ySize),super_new_clusters(in.xSize*in.ySize);
	in.pPixelList = &pixel_list;


	in.pPipeline = &ccd_pipeline;
	in.pipeline_size_minus_1 = pipeline_size-1;


	// calculating of tresholds - using (0,0) square of 
	// calculated bacground , in case more precise quatisation
	// is needed - tresholds are overwritten later :
	CalculateTresholds( in );
	out.m_PixelOut.treshold_for_not_more_then_n_exceeds = in.treshold_for_not_more_then_n_exceeds;

	out.ncnt = gCCDParams.m_nNeighbToSumCount;

	CManyTab2D<BIG_ELEM_TYPE>* pHomeoFrame = NULL;
	if(gCCDParams.m_bKeepHomeopaticSum)
		pHomeoFrame = ccd_pipeline.GetHomeopaticFrame();
	CManyTab2D<BIG_ELEM_TYPE>* pLaplaceFrame = NULL;
	if(gCCDParams.m_bKeepLaplaceFrame)
		pLaplaceFrame = ccd_pipeline.GetLaplaceFrame();

	// graphical cut for SUPER-NEW stars area :
	CGraphCut graph_cut;

	// loosy cut :
   // graph_cut.AddLineByPoints( 5000, 1250, 18750, 15000, eRelSmaller );
	// strict cut :
	graph_cut.AddLineByPoints( 15000, 7500, 25000, 17500, eRelSmaller );
   // graph_cut.AddHorizontalLine( 1250, eRelGreater );
	graph_cut.AddHorizontalLine( 500, eRelGreater ); // small stars too 
   graph_cut.AddHorizontalLine( 17500, eRelSmaller );


	inForPrev->pPixelList = &star_clusters;
	(*outForPrev) = out;

	int super_new=0;

	register int nNewLaplaceTreshold = ((in.pPipeline)->GetPipelineCfg()).m_nNewLaplace;
	register int nMaxLaplaceOnOther = ((in.pPipeline)->GetPipelineCfg()).m_nMaxLaplaceOnOther;
//	CWindowList& chip_bad_parts = ((in.pPipeline)->GetPipelineCfg()).m_BadPartsOfChip;

	// supernova tresholds :
	int TvForSN = ((in.pPipeline)->GetPipelineCfg()).m_fSuperNovaTvInSigma*in.SigmaLap;
	int TnForSN = ((in.pPipeline)->GetPipelineCfg()).m_fSuperNovaTnInSigma*in.SigmaLap;
	int MinPrevForSN = ((in.pPipeline)->GetPipelineCfg()).m_fSuperNovaMinPrevValue*in.SigmaLap;

	CImageCreator* pGenObj = ccd_pipeline.GetGenObj();

	clock_t total_in_check=0;
	for(register int i=0;i<size;i++){
		in.ccd_index = i;
		in.pCCDInfo = &(ccd_pipeline.GetCCDInfoTab()[in.ccd_index]);
		in.Matrix = &(newFrame[i]);
		in.p_data = (in.Matrix)->get_data_buffer();
		in.p_data_fast = (in.Matrix)->get_data_buffer_fast();
		in.p_curr_data_laplace = (in.Matrix)->get_frame_laplace_fast();
		in.p_curr_data_laplace_normal = (in.Matrix)->get_frame_laplace();			

		in.pCamCfg = (CCcdCfg*)((in.pPipeline)->GetCamCfgTab()[in.ccd_index]);		
		(in.pCCDInfo)->ReCalcAlfaMultTab( ((in.pCamCfg)->m_CCDParams).m_RotValueDAlfa, ccd_pipeline.GetPipelineSize() );
		

		if(pHomeoFrame){
			in.p_homeo_data = (*pHomeoFrame)[i].get_data_buffer();
			in.p_fast_homeo_data = (*pHomeoFrame)[i].get_data_buffer_fast();
		}
		if(pLaplaceFrame){
			in.p_laplace_data = (*pLaplaceFrame)[i].get_data_buffer();
			in.p_laplace_data_fast = (*pLaplaceFrame)[i].get_data_buffer_fast();
		}

		in.PrevMatrixPtrCnt = ccd_pipeline.GetAllMatrixPtrsChronologicalInt( i, in, TRUE, nFramesBack );
		// inForPrev->SetPrevMatrix( in );
		(*inForPrev) = in;
																					

		register int nIgnoreEdgeBottom = ccd_pipeline.GetPipelineCfg().m_nIgnoreEdgeBottom;
		register int nIgnoreEdgeUp = ccd_pipeline.GetPipelineCfg().m_nIgnoreEdgeUp;
		register int nIgnoreEdgeRight = ccd_pipeline.GetPipelineCfg().m_nIgnoreEdgeRight;
		register int nIgnoreEdgeLeft = ccd_pipeline.GetPipelineCfg().m_nIgnoreEdgeLeft;

		register int y_pos = (nIgnoreEdgeBottom-1)*in.xSize + nIgnoreEdgeLeft;
		register int pos = 0;
		register int low_y = nIgnoreEdgeBottom;
		register int low_x = nIgnoreEdgeLeft;

		register int genX = -1;
		register int genY = -1;
		if( pGenObj ){
			genX = pGenObj->m_LastPutObjectX;
			genY = pGenObj->m_LastPutObjectY;
		}

		int bContinue=1;

		// preparing statistics objects :
		CFrameIdentStat stat( nDayFrame ), sstat( nDayFrame );
		stat.nTotal = (yUpperLimit-nIgnoreEdgeBottom)*(xUpperLimit-nIgnoreEdgeLeft);
		if( gCCDParams.GetMC() && gCCDParams.m_bPutSample && pGenObj && genX>0 && genY>0 ){
			sstat.nTotal = 1;
			stat.nTotal = 0;				
		}
		

		BOOL_T bIsGenObject=FALSE;
		BOOL_T bDoBreak=FALSE;

		_TRACE_PRINTF_3("Cam%d Analysing area : (%d,%d)-(%d,%d)\n",pipeline_index,
					nIgnoreEdgeLeft,nIgnoreEdgeBottom,(xUpperLimit-1),(yUpperLimit-1));

		for(register int y=nIgnoreEdgeBottom;y<yUpperLimit && bContinue;y++){
			y_pos += in.xSize;
			pos = y_pos;

			for(register int x=nIgnoreEdgeLeft;x<xUpperLimit && bContinue;x++,pos++){

//				if( nDayFrame==41 && abs(x-487)<=1 && abs(y-234)<=1 ){
					//if(abs(x-75)<=0 && abs(y-8)<=0 )
//					printf("odo");
//				}

// only in MC :
#ifndef _OPTIMIZED_VERSION_
				// MONTE CARLO - comment out in real analysis :
				bIsGenObject=FALSE;
				if( bMC && gCCDParams.m_bPutSample && pGenObj && pGenObj->m_pGenEvent){
					if(x==genX && y==genY){	
						(pGenObj->m_pGenEvent)->m_PixelAnalResults.laplaceSum = in.p_curr_data_laplace[y][x];
					}
					if( abs(x-genX)<=gCCDParams.m_bGenEventRedial && abs(y-genY)<=gCCDParams.m_bGenEventRedial ){
						bIsGenObject=TRUE;
					}else{
						stat.nTotal++;
					}
				}else{
					stat.nTotal++;
				}
#endif

				if(in.p_data_fast[y][x]>=gCCDParams.m_MaxAllowedVal){
					continue;
				}


				if( gCCDParams.m_bCalcTresholdsByBackgrMap ){
					CalcTresholdsByMap( in, nNewLaplaceTreshold, nMaxLaplaceOnOther );
				}

#ifndef _OPTIMIZED_VERSION_ 
				if( !bIsGenObject ){
					stat.nMaxAllowedValue++;
				}else{
					sstat.nMaxAllowedValue=1;
				}
#endif

				if(in.p_curr_data_laplace[y][x]<=nNewLaplaceTreshold)
					continue;

#ifndef _OPTIMIZED_VERSION_
				// MONTE CARLO - comment out in real analysis :
				if( bIsGenObject && pGenObj && pGenObj->m_pGenEvent){
					if( in.p_curr_data_laplace[y][x] > (pGenObj->m_pGenEvent)->m_PixelAnalResults.laplaceSum )
					{
						(pGenObj->m_pGenEvent)->m_PixelAnalResults.laplaceSum = in.p_curr_data_laplace[y][x];
					}
				}
#endif

//				if( nDayFrame==205 && abs(x-1160)<=2 && abs(y-1248)<=2 ){
//					printf("odo");
//				}

				
				in.x = x;
            in.y = y;
            in.pos = pos;
#ifndef _OPTIMIZED_VERSION_
				if( !bIsGenObject ){
					stat.nTnewCut++;
				}else{
					sstat.nTnewCut=1;
				}
#endif
			
	
				//  out.m_PixelOut.maxAverageOfPrev = GetMaxAverageInVetoArea( in, TRUE );
				out.m_PixelOut.maxAverageOfPrev = -10000;
				for(register int v=0;v<gCCDParams.m_nVetoPointsCount;v++){
			       register long xx = (long)(x+(gCCDParams.m_VetoArea)[v].x);
			       register long yy = (long)(y+(gCCDParams.m_VetoArea)[v].y);
	
					 register long prevAverage = CalcAverageOfPrevN( xx, yy, in,  
																		 gCCDParams.m_nMaxOfAverageOfPrevN,
																		 TRUE );
					 if( prevAverage>out.m_PixelOut.maxAverageOfPrev ){
					     out.m_PixelOut.maxAverageOfPrev = prevAverage;
					 }
				}		
	
#ifndef _OPTIMIZED_VERSION_
				// MONTE CARLO - comment out in real analysis :
				// if(x==genX && y==genY && pGenObj && pGenObj->m_pGenEvent){
				if( bIsGenObject && pGenObj && pGenObj->m_pGenEvent){
					(pGenObj->m_pGenEvent)->m_PixelAnalResults.maxAverageOfPrev = out.m_PixelOut.maxAverageOfPrev;
				}
#endif
		
				out.m_PixelOut.eventType = eFlash;
				if(out.m_PixelOut.maxAverageOfPrev>=nMaxLaplaceOnOther){
					if(!gCCDParams.m_bCheckForSUPERNEW){
						continue;
					}else{
						// checking for SN 
						if( out.m_PixelOut.maxAverageOfPrev >= MinPrevForSN ){
							out.m_PixelOut.eventType = eBrighten;
						}else{
							continue;
						}
					}
				}

// Check bad areas of the chip here - so that it is not too time consuming :
				/*for(int w=0;w<chip_bad_parts.size();w++){
					CCDWindow& win = chip_bad_parts[w];
					if( x>=win.m_LowX && x<=win.m_UpX && y>=win.m_LowY && y<=win.m_UpY ){
						// pixels belongs to bad area of the chip and will be skiped :
						continue;
					}
				}*/


//				if( nDayFrame==205 && y>=1240 ){
//					printf("odo");
//				}

				if( m_pAnalFoundEventFunc ){
					HistoVariables( x, y, pos, in, out, stat, sstat, varInfo, in.SigmaLap );
				}


				if(!ApplyCuts( x, y, pos, in, out, stat, sstat, bIsGenObject, 
									bDoBreak, allclusters ) ){
					if( bDoBreak ){
						bContinue=0;
						break;
					}
					continue;
				}

				if( gCCDParams.m_bCheckForSUPERNEW && out.m_PixelOut.eventType==eBrighten ){
					if( out.cluster_cnt ){
						// cluster must have been calculated :
						if( !CheckForSN( in, out, TnForSN, TvForSN, MinPrevForSN ) ){
							continue;
						}						
					}else{
						continue;
					}
				}

				// NEW 20050403 - new change , added due to big background 
				// in ccdsingle analysis :
				if( gCCDParams.m_bCheckPrevOfMaxPixel ){
					int max_x = (int)out.m_PixelOut.m_MaxClusterPixel.x;
					int max_y = (int)out.m_PixelOut.m_MaxClusterPixel.y;

					if( max_x!=x || max_y!=y ){
						// in case MAX_PIXEL different then current one 
						// check if not edge of STAR :

						int maxPrev = -10000;

						for(register int v=0;v<gCCDParams.m_nVetoPointsCount;v++){
					       register long xx = (long)(max_x+(gCCDParams.m_VetoArea)[v].x);
					       register long yy = (long)(max_y+(gCCDParams.m_VetoArea)[v].y);
	
							 register long prevAverage = CalcAverageOfPrevN( xx, yy, in,  
																			 gCCDParams.m_nMaxOfAverageOfPrevN,
																			 TRUE );
							 if( prevAverage>maxPrev ){
							     maxPrev = prevAverage;
							 }
						}		

						if( maxPrev >= nMaxLaplaceOnOther){
							if(!gCCDParams.m_bCheckForSUPERNEW){
								continue;
							}else{
								// checking for SN 
								if( maxPrev >= MinPrevForSN ){
									out.m_PixelOut.eventType = eBrighten;
									out.m_PixelOut.maxAverageOfPrev = maxPrev;
								}else{
									continue;
								}
							}
						}
					}
				}

/*				if( gCCDParams.m_bRejectEventsNearShift && gCCDParams.m_bCheckFrame ){
					if( in.pPipeline && in.pPipeline->m_LastBadFrame>0 && in.pPipeline->m_LastBadLineY>0 ){
						if( abs( in.pPipeline->GetDayFrameCounter() - in.pPipeline->m_LastBadFrame )<=1 ){						
							if( abs( y - in.pPipeline->m_LastBadLineY ) <= 5 ){
								printf("Event (%d,%d) rejected due to bad frame (pixel shift at line %d)\n",x,y,in.pPipeline->m_LastBadLineY);
								continue;
							}
						}
					}
				}*/

				out.m_PixelOut.PixelRawValue = in.p_data_fast[y][x];
				if( out.m_PixelOut.eventType != eBrighten ){
					(in.Matrix)->AddFoundEvent( in, out );
					n_found_events++;
				}else{
					// super nova events added to other list :
					(in.Matrix)->AddFoundEvent( (in.pPipeline)->m_BrightenList, in, out, TRUE );
				}

				if( gCCDParams.m_bSkipIfMoreThen>0 && gCCDParams.m_eCheckIfMorePoint!=eAfterCoic ){
					// check limit for total number of events - reject all
					// if exceeded :
					if(gCCDParams.m_MaxNumberOfEventsOnFrame>0){
						if( (in.Matrix)->GetFoundEvents().size() > gCCDParams.m_MaxNumberOfEventsOnFrame){
							// in case limit exceeded , skip all and do not continue 
							printf("REJECTING ALL DUE TO %d > %d\n",(in.Matrix)->GetFoundEvents().size(),
									gCCDParams.m_MaxNumberOfEventsOnFrame);
							(in.Matrix)->GetFoundEvents().clear();
							bContinue=0;
							stat.bAllRejected=TRUE;
							sstat.bAllRejected=TRUE;
							break;
						}
					}
				}
			}
		}

		int nFrameEvents = (in.Matrix)->GetFoundEventCount();

		if( gCCDParams.m_bUseFoundPosition ){
			CCDEventList& events = (in.Matrix)->GetFoundEvents();
			for(int i=0;i<nFrameEvents;i++){
				events[i].m_MaxPoint = events[i].m_Point;
			}			
		}

		
		if( m_pAnalFoundEventFunc ){
			int afterTv = stat.nMinPrevCut;

			(*m_pAnalFoundEventFunc)( &afterTv, eHistoAfterTv, (in.pPipeline)->GetPipelineIndex(), NULL );
			(*m_pAnalFoundEventFunc)( &nFrameEvents, eHistoFrameEvents, (in.pPipeline)->GetPipelineIndex(), NULL );
		}

		int nGenEvents = 0;
		if( gCCDParams.GetMC() && pGenObj && genX>0 && genY>0 ){
			nGenEvents = (in.Matrix)->CountFoundEvents( genX, genY, 10 );
			if( nGenEvents>1 )
				nGenEvents = 1;
		}

		BOOL_T bToAll=FALSE;
		/*if( gCCDParams.m_bFitLineToSingleFrameToAll ){
			// fiting line to events on single frame :
			bToAll = CheckLineForEventsOnCurrFrame( (in.Matrix)->GetFoundEvents() );
			if( bToAll ){
				printf("There are events on line on this frame, pipeline : %d\n",(in.pPipeline)->GetPipelineIndex());
			}			
		}*/
		if( !bToAll ){
			if( gCCDParams.m_bFitLineToSingleFrame ){
				if( FitLineToSingleFrameEvents( (in.Matrix)->GetFoundEvents() ) ){
					printf("line fit SUCCEDED !!!\n");
				}
			}
		}

		// adding additional information to events :
		// before VerifySingleCamTracks - to have event time :
		AddAditionalEventInfo( in, &stat );

		if( gCCDParams.m_bCheckTracksOnSingleCam ){
			VerifySingleCamTracks( in, *(in.Matrix), (in.Matrix)->GetFoundEvents(), 
										  (in.pPipeline)->m_SingleCamEvents );
		}

		stat.nIfMoreThen = nFrameEvents;
		if(gCCDParams.m_bSkipIfMoreThen>0 && gCCDParams.m_eCheckIfMorePoint!=eAfterCoic){
	      // requires rejection of all events in case in the redial of gCCDParams.m
   	   // there is more then m_bSkipIfMoreThen events :
      	int nRejectedIfMore = (in.Matrix)->RejectIfMoreThen( in.pPipeline );

			nGenEvents = (in.Matrix)->CountFoundEvents( genX, genY, 10 );
			if( nGenEvents>1 )
            nGenEvents = 1;
						
			stat.nIfMoreThen = ( (in.Matrix)->GetFoundEventCount() - nGenEvents);
			if( stat.nIfMoreThen < 0 ){
				stat.nIfMoreThen = 0;
			}
			sstat.nIfMoreThen = nGenEvents;
	   }
		stat.nTracks = MAX( ( (in.Matrix)->GetFoundEventCount() - nGenEvents ) , 0 );
		sstat.nTracks = nGenEvents;
		if( gCCDParams.m_bLogFrameStat ){
			backgrStat.push_back( stat );
			if( gCCDParams.GetMC() && pGenObj && pGenObj->m_pGenEvent){
				// add to sample only if really added :
				samplesStat.push_back ( sstat );
			}
		}

		// log event rates :
		if(gCCDParams.m_bLogFrameStat && !gCCDParams.m_bCCDDouble){
			LogEventRates( in.pPipeline );
		}
		

		// CURRENTLY MOVED TO CALLING ROUTINE - now skips generated events 
		// so must be called after CompileEventsReport function !!!!!
		// no verify tracks :
      //if(gCCDParams.m_bCheckTracks){
      //   if((in.Matrix)->GetFoundEvents().size()>0){
            // in such case verify if not on track :
      //      VerifyIfEventsNotOnTrack( in, *(in.Matrix) );
      //   }
      //}

		// now add new events to pipeline events list :
 		// ccd_pipeline.m_allFoundEvents.back().push_back( (in.Matrix)->GetFoundEvents() );

	}
	PROFILER_END("Full analysis of new frame took : ")
	PROFILER_PUT("CheckIfOverlaps took : ",total_in_check);

	if( (in.Matrix)->GetFoundEvents().size()>0 ){
		bRet = TRUE;
		(in.Matrix)->SetEventTime();
		(in.Matrix)->GetFoundEvents().SetEvtIdx();
	}

	if( gCCDParams.m_bCheckForSUPERNEW && (in.pPipeline)->m_BrightenList.size()>0 ){
		(in.pPipeline)->m_BrightenList.SetEventTime( (in.pPipeline)->m_FrameUnixTime );
	}

	if( gCCDParams.m_nSamplesToPutOnFrame>1 ){
		FillSamplesInfo( in, nNewLaplaceTreshold );
	}

	delete inForPrev;
	delete outForPrev;

	_TRACE_PRINTF_3("checked for SUPER_NEW %d times\n",super_new);

	PrintFrameAnalyseInfo( in, (in.Matrix)->GetFoundEvents().size() );

	return bRet;
}


BOOL_T CCD_Analyser::AnalyseCompareToOldFrame( CCDPipeline& ccd_pipeline )
{
	if( gCCDParams.m_nCompareToOldFreqInSec > 0 ){
		if( strlen( ccd_pipeline.m_szPrevSumOfNFileName.c_str() ) && strlen( ccd_pipeline.m_szCurrSumOfNFileName.c_str() ) 
			 && (  ccd_pipeline.m_FrameUnixTime - ccd_pipeline.m_PrevSumOfFramesTime ) > gCCDParams.m_nCompareToOldFreqInSec && 
			 (ccd_pipeline.GetFrameIndex() % gCCDParams.m_bKeepSumOfPrevNFrames)==0 ){
			printf("OLD_ANAL : Performing comparison to old frame\n");


			// assuming prev sum frame is not in use later , so that we can read old
			// frame to it now :
			Table2D<BIG_ELEM_TYPE>* pCurrSumFrame = ccd_pipeline.GetCurrSum();
			Table2D<BIG_ELEM_TYPE>* pPrevSumFrame = ccd_pipeline.GetPrevSum( ccd_pipeline.m_szPrevSumOfNFileName.c_str() );
			Table2D<BIG_ELEM_TYPE>* pCurrSumFrameLap = ccd_pipeline.GetCurrSumLaplace();
			Table2D<BIG_ELEM_TYPE>* pPrevSumFrameLap = ccd_pipeline.GetPrevLaplace( ccd_pipeline.m_szPrevSumOfNFileName.c_str() ) ;

			// now read old frame into : pPrevSumFrame, pPrevSumFrameLap :


			// and perform analysis ... :
			BOOL_T ret =  AnalyseSumOfPrevNFrames(  ccd_pipeline, ccd_pipeline.m_EventsFromCompareToOld,
																 pCurrSumFrame, pPrevSumFrame, 
																 pCurrSumFrameLap, pPrevSumFrameLap );
			return ret;
						
		}
	}
	return FALSE;
}

BOOL_T CCD_Analyser::AnalyseSumOfPrevNFrames( CCDPipeline& ccd_pipeline )
{
	Table2D<BIG_ELEM_TYPE>* pCurrSumFrame = ccd_pipeline.GetCurrSum();
	Table2D<BIG_ELEM_TYPE>* pPrevSumFrame = ccd_pipeline.GetPrevSum();
	Table2D<BIG_ELEM_TYPE>* pCurrSumFrameLap = ccd_pipeline.GetCurrSumLaplace();
	Table2D<BIG_ELEM_TYPE>* pPrevSumFrameLap = ccd_pipeline.GetPrevSumLaplace();

	BOOL_T ret =  AnalyseSumOfPrevNFrames(  ccd_pipeline, ccd_pipeline.m_EventsOnSumedFrame,
														 pCurrSumFrame, pPrevSumFrame, 
														 pCurrSumFrameLap, pPrevSumFrameLap );
	return ret;
}

BOOL_T CCD_Analyser::AnalyseSumOfPrevNFrames( CCDPipeline& ccd_pipeline,
															 CCDEventList& event_list,
															 Table2D<BIG_ELEM_TYPE>* pCurrSumFrame,
															 Table2D<BIG_ELEM_TYPE>* pPrevSumFrame,
															 Table2D<BIG_ELEM_TYPE>* pCurrSumFrameLap,
															 Table2D<BIG_ELEM_TYPE>* pPrevSumFrameLap  )
{	
	Initialize();
	CCDMatrix& matrix = (ccd_pipeline.GetCurrent())[0];
	BOOL_T bMC = gCCDParams.GetMC();
	LONG_T pipeline_size = ccd_pipeline.GetPipelineSize();
	LONG_T ClusterWithMoreCnt=0;
	LONG_T n_found_events = 0;
	LONG_T n_accepted_on_level_1=0;
	LONG_T n_accepted_on_level_2=0;
   LONG_T star_cluster[MAX_CLUSTER_SIZE];
	LONG_T nFramesBack = gCCDParams.m_FramesBack;
	int nDayFrame = ccd_pipeline.GetDayFrameCounter();
	register LONG_T max_pos=-1;

	if( ccd_pipeline.GetFrameIndex()< (2*gCCDParams.m_bKeepSumOfPrevNFrames))
		return FALSE;


	BOOL_T bRet = FALSE;

	PROFILER_START

	BIG_ELEM_TYPE** p_curr_sum_fast = pCurrSumFrame->get_data_buffer_fast();
	BIG_ELEM_TYPE* p_curr_sum = pCurrSumFrame->get_data_buffer();
	BIG_ELEM_TYPE** p_prev_sum_fast = pPrevSumFrame->get_data_buffer_fast();
	BIG_ELEM_TYPE* p_prev_sum = pPrevSumFrame->get_data_buffer();
	
	BIG_ELEM_TYPE** p_curr_lap_fast = pCurrSumFrameLap->get_data_buffer_fast();
	BIG_ELEM_TYPE** p_prev_lap_fast = pPrevSumFrameLap->get_data_buffer_fast();

	CPixelAnalyseIn in;
   CPixelAnalyseOut out;

	in.frame_index = ccd_pipeline.GetFrameIndex();
	in.xSize = pCurrSumFrame->GetXSize();
   in.ySize = pCurrSumFrame->GetYSize();	
	m_FirstDy = ((m_NeighbPoints[0].y)*in.xSize);
	in.ccd_index = 0;
	in.pCCDInfo = &(ccd_pipeline.GetCCDInfoTab()[0]);
	in.p_curr_data_laplace = p_curr_lap_fast;
	in.p_curr_data_laplace_normal = pCurrSumFrameLap->get_data_buffer();			
	// in.pCamCfg = &(ccd_pipeline.m_PipelineCfg);
	CPixelList pixel_list(in.xSize*in.ySize);
	in.pPixelList = &pixel_list;
	in.pPipeline = &ccd_pipeline;
	in.pipeline_size_minus_1 = pipeline_size-1;
	in.p_prev_lap_fast = p_prev_lap_fast;


	// calculating of tresholds - using (0,0) square of 
	// calculated bacground , in case more precise quatisation
	// is needed - tresholds are overwritten later :
	CalculateTresholds( in );


	int yUpperLimit = (in.ySize-gCCDParams.m_nIgnoreEdgeUp);
	int xUpperLimit = (in.xSize-gCCDParams.m_nIgnoreEdgeRight);
	int nIgnoreEdgeBottom = ccd_pipeline.GetPipelineCfg().m_nIgnoreEdgeBottom;
	int nIgnoreEdgeUp = ccd_pipeline.GetPipelineCfg().m_nIgnoreEdgeUp;
	int nIgnoreEdgeRight = ccd_pipeline.GetPipelineCfg().m_nIgnoreEdgeRight;
	int nIgnoreEdgeLeft = ccd_pipeline.GetPipelineCfg().m_nIgnoreEdgeLeft;
	CPixelList allclusters(in.xSize*in.ySize);

	int map_size=1;
   InfoTable2D info( in.xSize, in.ySize, map_size, map_size, FALSE, gCCDParams.m_nIgnoreEdge );

   my_printf_now("calculating bacground map ...");
	double mean,sigma;
	for( int i=0;i<map_size;i++){
      for( int j=0;j<map_size;j++){
         Area2DInfo& elem = info.GetElem(i,j);
         pCurrSumFrame->GetVariableMeanAndSigma( gCCDParams.m_eLaplaceType, mean, sigma, elem,
                                  (int)elem.m_LowLeft.x, (int)elem.m_LowLeft.y,
                                  (int)elem.m_UpRight.x, (int)elem.m_UpRight.y,
                                  0, pCurrSumFrameLap->get_data_buffer_fast() );
      }
   }
   my_printf_now("OK\n");	

	CLongPoint vetopoints[100];
	int nVetoCount=CalcShapePoints( vetopoints, gCCDParams.m_eVetoShape, gCCDParams.m_fVetoRadiusOnSumFrame );


	// currently map_size=1
	Area2DInfo& elem = info.GetElem( 0 ,0 ); 
	int nSigmaBacground = elem.m_DataInfo[gCCDParams.m_eLaplaceType].m_Sigma;

   int nNewLaplaceTreshold = (int)(elem.m_DataInfo[gCCDParams.m_eLaplaceType].m_Average
									+gCCDParams.m_nNewLaplaceInSigmaOnSum*elem.m_DataInfo[gCCDParams.m_eLaplaceType].m_Sigma);
   int prevTresh = (int)(elem.m_DataInfo[gCCDParams.m_eLaplaceType].m_Average
							+gCCDParams.m_nMaxLaplaceOnOtherInSigmaOnSum*elem.m_DataInfo[gCCDParams.m_eLaplaceType].m_Sigma);
	int clusterTreshold = nNewLaplaceTreshold;

	// supernova tresholds :
	int TvForSN = ((in.pPipeline)->GetPipelineCfg()).m_fSuperNovaTvInSigma*nSigmaBacground;
	int TnForSN = ((in.pPipeline)->GetPipelineCfg()).m_fSuperNovaTnInSigma*nSigmaBacground;
	int MinPrevForSN = ((in.pPipeline)->GetPipelineCfg()).m_fSuperNovaMinPrevValue*nSigmaBacground;


	double x0,y0;
   double max_noise_level = in.treshold_for_cluster;// WAS 0

	event_list.clear();

	int nAboveTn=0;

	// statisctics object 
	CFrameIdentStat stat( nDayFrame );

	for(register int y=nIgnoreEdgeBottom;y<yUpperLimit;y++){
		register int y_pos = y*in.xSize;
		for(register int x=nIgnoreEdgeLeft;x<xUpperLimit;x++){

				stat.nTotal++;
				if( p_curr_sum_fast[y][x]>=gCCDParams.m_MaxAllowedVal){
					continue;
				}
				stat.nMaxAllowedValue++;

				register int pos = y_pos + x;				
				if( in.pPixelList->CheckPixel( pos ) )
					continue;
				stat.nClusterOverlap++;

				if( p_curr_lap_fast[y][x]<=nNewLaplaceTreshold)
					continue;
				stat.nTnewCut++;

				nAboveTn++;

				register int prevAverage = p_prev_lap_fast[y][x];

				for(register int v=0;v<nVetoCount;v++){
			       register long xx = (long)(x+(vetopoints)[v].x);
			       register long yy = (long)(y+(vetopoints)[v].y);

					/*if(gCCDParams.m_bCalcMaxNeighbHomeo){
						register long low_x = xx;
						register long up_x = xx;
						register long low_y = yy;
						register long up_y = yy;
						if(  p_prev_lap_fast[y][x] > nNewLaplaceTreshold ){
							low_x -= ccd_pipeline.m_PipelineCfg.m_nCalcMaxNieghbRedial;
							low_y -= ccd_pipeline.m_PipelineCfg.m_nCalcMaxNieghbRedial;
							up_x  += ccd_pipeline.m_PipelineCfg.m_nCalcMaxNieghbRedial;
							up_y  += ccd_pipeline.m_PipelineCfg.m_nCalcMaxNieghbRedial;
						}


					  for(register int i=low_x;i<=up_x;i++){
					  	 for(register int j=low_y;j<=up_y;j++){	
							if(i>=0 && j>=0 && i<in.xSize && j<in.ySize){
								if(p_prev_lap_fast[j][i]>prevAverage)
									prevAverage = p_prev_lap_fast[j][i];
							}
						 }
					  }
					}*/
					if(xx>=0 && yy>=0 && xx<in.xSize && yy<in.ySize){
						if(p_prev_lap_fast[yy][xx]>prevAverage){
							prevAverage = p_prev_lap_fast[yy][xx];
						}
					}
				}		

				out.m_PixelOut.eventType = eFlash;
				if( prevAverage>prevTresh ){
					if(!gCCDParams.m_bCheckForSUPERNEW){
						continue;
					}else{
						 // checking for SN
                  if( prevAverage >= MinPrevForSN ){
                     out.m_PixelOut.eventType = eBrighten;
                  }else{
                     continue;
                  }
					}
				} 
				stat.nTprevCut++;
				in.x = x;
				in.y = y;
				in.pos = pos;


				if( gCCDParams.m_bMinLaplaceOnOtherEnabled ){
					if( prevAverage < gCCDParams.m_nMinLaplaceOnOther ){
						continue;
					}						
				}
				stat.nMinPrevCut;


				out.m_PixelOut.maxAverageOfPrev = prevAverage;
				out.m_PixelOut.laplaceSum = p_curr_lap_fast[y][x];				

				// now tn, tv cut passed , determine cluster size :
				FindClusterAboveTresholdOpt3( in, out, out.cluster, out.cluster_cnt,
														x0,y0, max_noise_level, clusterTreshold , 
														FALSE, TRUE , FALSE );
				if(gCCDParams.m_bSkipOverlaps){						
					register LONG_T prev_x,prev_y;
					if( matrix.CheckIfOverlapsFast( event_list, 
															  x, y, out.cluster, 0 , prev_x, prev_y ) ){
						continue;
						// return FALSE;
					}
				}
				stat.nOverlap++;

				// calculate sphericity - but no cut on this level :
				if( gCCDParams.m_CheckEventShape>=0){
			   	double max_redial = CCDFastCalc::FindMaxRedial( out.cluster, out.cluster_cnt, x0, y0, in.xSize );
			   	out.m_PixelOut.m_Sphericity = CCDFastCalc::CalcSphericity( max_redial, out.cluster_cnt );
				}
				stat.nShapeCut++;
				
				if( gCCDParams.m_bCheckForSUPERNEW && out.m_PixelOut.eventType==eBrighten ){				
					if( out.cluster_cnt ){
      			  // cluster must have been calculated :
			        if( !CheckForSN_OnSum( in, out, TnForSN, TvForSN, MinPrevForSN ) ){
			           continue;
			        }
     				}else{
			        continue;
			     }
				}

				if( out.m_PixelOut.eventType == eFlash ){
					matrix.AddFoundEvent( event_list , in, out );
				}else{
					// super nova events added to other list :
				  matrix.AddFoundEvent( (in.pPipeline)->m_BrightenOnSumList, in, out, TRUE );		
				}
		}
	}

	int evt_count = event_list.size();
	if( evt_count>0 ){
		bRet = TRUE;
		time_t ut_time = (time_t)matrix.getObsTime( TRUE );
		event_list.SetEventTime( ut_time );
		event_list.SetEvtIdx();
		for(int i=0;i<evt_count;i++){
			event_list[i].m_nFrameStarCount = nAboveTn;
		}
	}

	if( gCCDParams.m_bCheckForSUPERNEW && (in.pPipeline)->m_BrightenOnSumList.size()>0 ){
		(in.pPipeline)->m_BrightenOnSumList.SetEventTime( (in.pPipeline)->m_FrameUnixTime );
		printf("Cam%d Number of SN events on SUM = %d\n",
			(in.pPipeline)->GetPipelineIndex(),
			(in.pPipeline)->m_BrightenOnSumList.size());
	}


	sumBackgrStat.clear();
	sumBackgrStat.push_back( stat );
	
	printf("Found %d events on sum frame\n",event_list.size());
	PROFILER_END("Full analysis of sum frame frame took : ")

	return bRet;
}

BOOL_T CCD_Analyser::CheckIfEdgeOfBigStar( CPixelAnalyseIn& in, CPixelAnalyseOut& out,
														 BOOL_T bClusterFound, int nMaxLaplaceOnOther )
{
	// optional check if not fluctuating edge of big star :
	if(gCCDParams.m_bCheckIfNotEdgeOfBigStar){
		int x = in.x;
		int y = in.y;
		// using ClusterWithMore - not used anymore in this algorithm :
			
		int xx=-1,yy=-1;			
		double cluster_treshold=in.treshold_for_cluster;
						
		BOOL_T bClusterHasPoints=FALSE;
		if(gCCDParams.m_bEdgeOfBigByRawData){
			// double tresh_for_cluster_raw = in.treshold_for_raw_cluster;
			// double tresh_for_cluster_raw_aver = in.treshold_for_raw_cluster_aver;
			double tresh_for_cluster_raw = (in.pPipeline)->GetTreshold( x, y, in.ccd_index, eRawS, gCCDParams.m_nSigmaAboveMeanInRawCluster );
			double tresh_for_cluster_raw_aver = (in.pPipeline)->GetTreshold( x, y, in.ccd_index, eRawS,
															gCCDParams.m_nSigmaAboveMeanInRawCluster, TRUE, in.nAver );
			int max_pos=-1;
			double x0,y0;
			FindClusterAboveTresholdRawDataOpt3( in, out, x0, y0, 
															 out.ClusterWithMore, out.ClusterWithMoreCnt,
															 tresh_for_cluster_raw, FALSE);						
			FindMaxInCluster( in.p_data, in.xSize, out.ClusterWithMore, out.ClusterWithMoreCnt, max_pos );


			if( out.ClusterWithMoreCnt>0 ){
				int x_max = (max_pos % in.xSize);
				int y_max = (max_pos / in.xSize );

				// now find clusters on previous frames starting from pixel (x_max,y_max)
				// function FindAverageClusterOnPrev also initializes out.m_pCluster with zeros
				FindAverageClusterOnPrevRawData( in, out, out.ClusterWithMore, out.ClusterWithMoreCnt,						
														 (int)tresh_for_cluster_raw, x_max, y_max );
				xx = (x-x_max)+((out.m_pCluster->GetXSize())/2);
				yy = (y-y_max)+((out.m_pCluster->GetYSize())/2);
				bClusterHasPoints = TRUE;
			}
		}else{
			int max_pos=-1;
			double max_noise_level;
         if(!bClusterFound){
				double x0,y0;
				FindClusterAboveTresholdOpt3( in, out, x0,y0, max_noise_level, cluster_treshold );								
			}						
			FindMaxInClusterBig( in.p_curr_data_laplace_normal, in.xSize, out.cluster, out.cluster_cnt, max_pos );	

			if( out.cluster_cnt > 0 ){
				LONG_T x_max = (max_pos % in.xSize);
				LONG_T y_max = (max_pos / in.xSize );
				FindAverageClusterOnPrevLaplace( in, out, out.ClusterWithMore, out.ClusterWithMoreCnt,
															(int)cluster_treshold, x_max, y_max, in.nAver );
																  
				xx = (x-x_max)+((out.m_pCluster->GetXSize())/2);
				yy = (y-y_max)+((out.m_pCluster->GetYSize())/2);
				bClusterHasPoints = TRUE;
			}
		}

		// (xx,yy) - position of event on m_pCluster - centered in x_max,y_max :
		if(bClusterHasPoints && xx>=0 && yy>=0 && xx<(out.m_pCluster)->GetXSize() && yy<(out.m_pCluster)->GetYSize() ){
			/*int value = ((out.m_pCluster)->get_data_buffer_fast())[yy][xx];
			if( (gCCDParams.m_bEdgeOfBigByRawData && value>0) ||
				 (!gCCDParams.m_bEdgeOfBigByRawData && value>=in.treshold_for_prev)  ){
				// pixel belongs to cluster of BIG STAR :
				// value - is >=in.treshold_for_raw_cluster_aver and not 0
				printf("Frame# : %d at (%d,%d) Rejected by EDGE OF BIG STAR !!!\n",(in.pPipeline)->GetFrameIndex(),x,y);
				continue;
			}*/
			//mystring szOutName;
			//szOutName << "TMP/cluster_" << in.frame_index << "_" << x << "-" << y  << ".fit";
			//CCDUtil::WriteToFITSFile( *(out.m_pCluster), szOutName.c_str() );

			int max_x,max_y;
			int maxVal = (out.m_pCluster)->GetMaxValue( max_x, max_y );
			int laplaceValue = maxVal;
			BOOL_T bDoCheck = (maxVal>0 && max_x>2 && max_y>2 &&
									max_x<((out.m_pCluster)->GetXSize()-2) &&
									max_y<((out.m_pCluster)->GetYSize()-2));
			if( bDoCheck ){	
				if( gCCDParams.m_bEdgeOfBigByRawData ){
					_TRACE_PRINTF_1("cehcking : (%d,%d) for (%d,%d)\n",max_x,max_y,x,y);
					laplaceValue = Table2D<BIG_ELEM_TYPE>::CalcLaplaceSum( max_x, max_y, in.xSize, 
												(out.m_pCluster)->get_data_buffer_fast(),gCCDParams.m_eLaplaceType);
				}
				if( laplaceValue >= nMaxLaplaceOnOther ){
					// point belongs to cluster which is star and 
					// exceeds value of veto treshold 
					_TRACE_PRINTF_1("Frame# : %d at (%d,%d) Rejected by EDGE OF BIG STAR !!!\n",(in.pPipeline)->GetFrameIndex(),x,y);
					return FALSE;
				}else{	
					_TRACE_PRINTF_1("Frame# : %d at (%d,%d) - Accepted by EDGE OF BIG STAR !!!\n",(in.pPipeline)->GetFrameIndex(),x,y);
				}
			}
		}
	}
	return TRUE;
}


void CCD_Analyser::LogSumEventsRates( CCDPipeline* pPipeline )
{
	if(gCCDParams.m_bLogFrameStat){
		int k=0;

		// in case more then back frames for tracks we can log to file :			
		CFrameEventStatTab::iterator i;
		for(i=sumBackgrStat.begin();i!=sumBackgrStat.end();i++,k++){
			i->LogToFile( pPipeline, "sum_backgr_event_rates", "SumEvents" );
		}
		sumBackgrStat.clear();
		// sumBackgrStat.erase( backgrStat.begin(), pDelEnd+1 );
	}
}

void CCD_Analyser::LogEventRates( CCDPipeline* pPipeline )
{
	if(gCCDParams.m_bLogFrameStat){
		int k=0;
		if( backgrStat.size()>=(gCCDParams.m_nNumBackFramesForTracks+5)){
			// in case more then back frames for tracks we can log to file :			
			CFrameEventStatTab::iterator pDelEnd;
			CFrameEventStatTab::iterator i;
			for(i=backgrStat.begin();i!=backgrStat.end();i++,k++){
				if( k<=5 ){
					i->LogToFile( pPipeline, "backgr_event_rates" );
				}
				if(k==5){
					pDelEnd = i;
				}				
			}
			backgrStat.erase( backgrStat.begin(), pDelEnd+1 );
		}

		
		if( gCCDParams.GetMC() && samplesStat.size()>=(gCCDParams.m_nNumBackFramesForTracks+5)){
			// in case more then back frames for tracks we can log to file :			
			CFrameEventStatTab::iterator pDelEnd;
			CFrameEventStatTab::iterator i;
			k=0;
			for(i=samplesStat.begin();i!=samplesStat.end();i++,k++){
				if( k<=5 ){
					i->LogToFile( pPipeline, "sample_event_rates" );
				}
				if(k==5){
					pDelEnd = i;
				}				
			}
			samplesStat.erase( samplesStat.begin(), pDelEnd+1 );
		}
	}
}

void CCD_Analyser::LogAllEventRates( CCDPipeline* pPipeline )
{
	if(gCCDParams.m_bLogFrameStat){
		CFrameEventStatTab::iterator i;
		for(i=backgrStat.begin();i!=backgrStat.end();i++){
			i->LogToFile( pPipeline, "backgr_event_rates" );
		}
		backgrStat.clear();

		
		if( gCCDParams.GetMC() ){
			// in case more then back frames for tracks we can log to file :			
			CFrameEventStatTab::iterator i;
			for(i=samplesStat.begin();i!=samplesStat.end();i++){
				i->LogToFile( pPipeline, "sample_event_rates" );
			}			
		}
		samplesStat.clear();
	}	
}

void CCD_Analyser::CalcCoordForEvents( CCDPipeline& ccd_pipeline, CCDProcState* pFrameInfo,
													CCDEventList& events, CCDMatrix& matrix,
													double sigmaB )
{
	for(CCDEventList::iterator pEvt=events.begin();pEvt!=events.end();pEvt++){
		// celestial coordinates :
		
		pEvt->CalcAstroCoordinates( matrix, pFrameInfo );

		//pEvt->m_PixelAnalResults.m_Likehood = get_significance( pEvt->m_PixelAnalResults.m_MaxClusterValue, 
		//												pEvt->m_PixelAnalResults.maxAverageOfPrev,
		//										      sigmaB, meanB, gCCDParams.m_eLaplaceType, 
		//												gCCDParams.m_nMaxOfAverageOfPrevN,
		//												gCCDParams.m_bCalcMaxNeighbHomeo,
		//												gCCDParams.m_bIsCOIC );
		pEvt->m_PixelAnalResults.m_Significance  = ( pEvt->m_PixelAnalResults.laplaceSum - pEvt->m_PixelAnalResults.maxAverageOfPrev )/sigmaB;
	}
}

void CCD_Analyser::CalcAdditionalInfoForEvents( CCDPipeline& ccd_pipeline )
{
	if( gCCDParams.m_bReadFirstLevelInfoFromLog ){
		printf("re-analysing log files, no recalculation is performed : CCD_Analyser::CalcAdditionalInfoForEvents\n");
		return;
	}	

	cCCD& currFrame = ccd_pipeline.GetCurrent();
	int nCount=0;
	InfoTable2D* pBackgrMap = ccd_pipeline.GetBackgroundMap();
	double sigmaB,meanB;
	pBackgrMap->GetMeanAndSigmaBackgr( 0, 0, meanB, sigmaB, gCCDParams.m_eLaplaceType );

   for(register int i=0;i<currFrame.GetCount();i++){
		CCDMatrix& matrix = currFrame[i];
		CCDEventList& events = matrix.GetFoundEvents();
		CCDProcState* pFrameInfo = &((ccd_pipeline.GetCCDInfoTab())[i]);

		CalcCoordForEvents( ccd_pipeline, pFrameInfo, events, matrix, sigmaB );
	}


	if( gCCDParams.m_bAnalyzeSumOfPrevNFrames && ccd_pipeline.m_EventsOnSumedFrame.size()>0 ){	
		CCDProcState* pFrameInfo = &((ccd_pipeline.GetCCDInfoTab())[0]);
		CalcCoordForEvents(  ccd_pipeline, pFrameInfo,
									ccd_pipeline.m_EventsOnSumedFrame, currFrame[0], 1.00 );
		if( gCCDParams.m_bCheckForSUPERNEW  ){
			CalcCoordForEvents(  ccd_pipeline, pFrameInfo,
								 		ccd_pipeline.m_BrightenOnSumList, currFrame[0], 1.00 );
		}
	}

	if( gCCDParams.m_bCheckForSUPERNEW  ){
		CCDProcState* pFrameInfo = &((ccd_pipeline.GetCCDInfoTab())[0]);
		CalcCoordForEvents(  ccd_pipeline, pFrameInfo,
							 		ccd_pipeline.m_BrightenList, currFrame[0], 1.00 );
	}
}

void CCD_Analyser::CalcAdditionalInfoForEvents( CCDPipeline& ccd_pipeline, CCDEventList& events )
{
	cCCD& currFrame = ccd_pipeline.GetCurrent();
	int nCount=0;
	InfoTable2D* pBackgrMap = ccd_pipeline.GetBackgroundMap();
	double sigmaB,meanB;
	pBackgrMap->GetMeanAndSigmaBackgr( 0, 0, meanB, sigmaB, gCCDParams.m_eLaplaceType );

	CCDMatrix& matrix = currFrame[0];
	CCDProcState* pFrameInfo = &((ccd_pipeline.GetCCDInfoTab())[0]);

	CalcCoordForEvents( ccd_pipeline, pFrameInfo, events, matrix, sigmaB );
}


void CCD_Analyser::CalcAstroCoorinatesOfNewEvents( CCDPipeline& ccd_pipeline )
{
	cCCD& currFrame = ccd_pipeline.GetCurrent();
	int nCount=0;
   for(register int i=0;i<currFrame.GetCount();i++){
		CCDMatrix& matrix = currFrame[i];
		CCDEventList& events = matrix.GetFoundEvents();
		CCDProcState* pFrameInfo = &((ccd_pipeline.GetCCDInfoTab())[i]);

		for(CCDEventList::iterator pEvt=events.begin();pEvt!=events.end();pEvt++){
			pEvt->CalcAstroCoordinates( matrix, pFrameInfo );
		}
	}
}

void CCD_Analyser::PrintFrameAnalyseInfo( CPixelAnalyseIn& in, int found_cnt )
{
	mystring szPath;
	//szPath << gCCDParams.GetOutputDir() << "/FrameAnalyseInfo/Cam" << (in.pPipeline)->GetPipelineIndex()
	//	    << "/frames.log";
	CCDPipeline::GetOutFileName( szPath, "FrameAnalyseInfo", "frames.log", in.pPipeline, -1, FALSE );

	BOOL_T bFirstLine=FALSE;
	if(!MyFile::DoesFileExist( szPath.c_str() ))
		bFirstLine = TRUE;

	MyOFile out( szPath.c_str(), "a" );
	if(bFirstLine){
		out.Printf( (gCCDParams.m_FrameAnalHeader).c_str() );
	}
	
	BOOL_T bFitG=FALSE,bFitS=FALSE;
	int frame = (in.pPipeline)->GetFrameIndex();
	int tn = (int)((in.pPipeline)->GetPipelineCfg()).m_nNewLaplace;
	int tv = (int)((in.pPipeline)->GetPipelineCfg()).m_nMaxLaplaceOnOther;	
	int sigmaG = (int)(in.pPipeline)->GetBackgroundMap()->GetSigmaBackground( 0 ,0 , gCCDParams.m_eLaplaceType, bFitG );
	int meanG = (int)(in.pPipeline)->GetBackgroundMap()->GetMeanBackground( 0 ,0 , gCCDParams.m_eLaplaceType );
	int sigmaS = (int)(in.pPipeline)->GetBackgroundMap()->GetSigmaBackground( 0 , 0, eRawS, bFitS );
	int meanS = (int)(in.pPipeline)->GetBackgroundMap()->GetMeanBackground( 0 , 0, eRawS );
	int t_hot = in.treshold_for_hot;
	int t_cluster = (int)in.treshold_for_cluster;

	out.Printf( (gCCDParams.m_FrameAnalLogFmt).c_str(), frame, 
					tn,tv, meanG, sigmaG, meanS, sigmaS, found_cnt, t_hot,
					t_cluster, bFitG, bFitS  );
	out.Close();		
}


void CCD_Analyser::AddAditionalEventInfo( CPixelAnalyseIn& in, CFrameIdentStat* stat )
{
	CCDEventList& newEvents = (in.Matrix)->GetFoundEvents();
	
	struct tm frameDataTime;
	BOOL_T bDateTime = (in.Matrix)->getObsTime( &frameDataTime );

	// date from header must be converted :
	convert_date_from_fits_header_format( &frameDataTime );

	time_t ut_time = timegm( &frameDataTime ); // or use my_timegm - if timegm not present on your platform

	for( int i=0;i<newEvents.size();i++){
		if(bDateTime){
			newEvents[i].m_Time = ut_time;						
		}
		if( (in.pPipeline)->GetWorkingMode() == eDAQSatTriggerMode ){
			newEvents[i].m_PixelAnalResults.m_bInTriggerMode = TRUE;
		}
		if( stat ){
			newEvents[i].m_nFrameStarCount = stat->nTnewCut;
		}
	}
}

int CCD_Analyser::VerifyTracks( CCDPipeline& ccd_pipeline, int nConfirmOnNext, 
										  BOOL_T bForceCheck /*=FALSE*/ )
{
	CPixelAnalyseIn in;

	in.frame_index = ccd_pipeline.GetFrameIndex();	
	in.pPipeline = &ccd_pipeline;
	cCCD& frame = ccd_pipeline.GetCurrent();
	long size = frame.GetCount();
	int nTracks=0;

	for(register int i=0;i<size;i++){
		 // no verify tracks :		
		in.ccd_index = i;
		in.pCCDInfo = &(ccd_pipeline.GetCCDInfoTab()[in.ccd_index]);
		in.Matrix = &(frame[i]);

      if(gCCDParams.m_bCheckTracks){
			int nEventCount = (in.Matrix)->GetFoundEvents().size();
			if( nConfirmOnNext>0 ){
				deque<CFrameEvents>& allPipelineEvents = in.pPipeline->m_allFoundEvents;
				int confirmedIndex = (allPipelineEvents.size()-1-nConfirmOnNext);
				if( confirmedIndex<0 )
					return 0;
				nEventCount = allPipelineEvents[confirmedIndex][i].size();
			}

			// if there are events and there is less then upper rate limit for track fits :
         if( nEventCount>0 || bForceCheck ){
				if( nEventCount<gCCDParams.m_MaxEventRate ){
	            // in such case verify if not on track :
   	         BOOL_T bRet = VerifyIfEventsNotOnTrack( in, *(in.Matrix) );
					if(bRet){
						nTracks++;
					}
				}else{
					printf("EventRate=%d, exceeds limit %d, no track fit\n",nEventCount,gCCDParams.m_MaxEventRate);
				}
         }
      }
	}
	return (nTracks>0);
}

int CCD_Analyser::VerifyPlaneTracks( CCDPipeline& ccd_pipeline, BOOL_T bDoCheckForNewTracks,
												 BOOL_T bDoCheckVelocity,
			                            double fVelError  )
{
	CPixelAnalyseIn in;

	in.frame_index = ccd_pipeline.GetFrameIndex();	
	in.pPipeline = &ccd_pipeline;
	cCCD& frame = ccd_pipeline.GetCurrent();
	long size = frame.GetCount();

	int ret=0;
	for(register int i=0;i<size;i++){
		 // no verify tracks :		
		in.ccd_index = i;
		in.pCCDInfo = &(ccd_pipeline.GetCCDInfoTab()[in.ccd_index]);
		in.Matrix = &(frame[i]);

      if(gCCDParams.m_bCheckPlaneTracks){
			int nEventCount = (in.Matrix)->GetFoundEvents().size();

			// if there are events and there is less then upper rate limit for track fits :
         if( nEventCount>0 ){
            // in such case verify if not on track :
            if( VerifyPlaneTracks( in, *(in.Matrix), bDoCheckForNewTracks,
											  bDoCheckVelocity, fVelError ) ){
					ret++;
				}
         }
      }
	}
	return ret;
}


LONG_T CCD_Analyser::VerifyEventsOnImage( CCDEventList& eventToVerify, CPixelAnalyseIn& in )
{
	// CCDEventList& eventToVerify = FrameToVerify.GetFoundEvents();
	// CCDEventList verifiedEvents;
	CCDEventList::iterator pEvt;
	LONG_T nextClusterSum[MAX_PREV_SUMS];
	LONG_T nVerified=0;

	_TRACE_PRINTF_5("----------------------------------------------\n");
	_TRACE_PRINTF_5("Verification on next frames :\n");

	CCDEventList finallist1;
	for(pEvt = eventToVerify.begin();pEvt!=eventToVerify.end();pEvt++){
		pEvt->m_bCheckOnNext = TRUE;

		// starting from 1 , at 0 there is frame to be confirmed :
		for(int next=1;next<in.PrevMatrixPtrCnt;next++){
			double next_x,next_y;
			CCDDataResults::CalcStarPositionAuto( pEvt->m_MaxPoint.x, pEvt->m_MaxPoint.y,
															  in.PrevFramesTime[0], in.PrevFramesTime[next], next ,
															  next_x, next_y,
															  (in.pCamCfg->m_CCDParams).m_RotCenterX,
															  (in.pCamCfg->m_CCDParams).m_RotCenterY,
//															  (in.pCamCfg->m_CCDParams).m_FrameDXPerSec,(in.pCamCfg->m_CCDParams).m_FrameDYPerSec,
															 (in.pCamCfg->m_CCDParams).m_FrameDX,(in.pCamCfg->m_CCDParams).m_FrameDY,
															  (in.pCamCfg->m_CCDParams).m_bUseRotInAverageOfPrevN,
															 (in.pCamCfg->m_CCDParams).m_RotValueDAlfa,
															  in.pCCDInfo, (in.pCamCfg->m_CCDParams).m_RotValueDAlfaPerSec,
															  &((in.pPipeline)->m_PipelineCfg) );

															  															  


			if(pEvt->m_bIdentified){	
				if((pEvt->m_PixelAnalResults).eventType==eFlash){
					register int nextFrameX = my_round(next_x);
					register int nextFrameY = my_round(next_y);
					
					BOOL_T bVerif = CheckNextFrameCondition( (CCDMatrix*)in.PrevMatrixPtr[next], nextFrameX, 
																		  nextFrameY, &(*pEvt), next, in );
					if(!bVerif){
						pEvt->RejectByNextFrames();
					}else{
						// events is verfieid on next frame now check with cataloge :
						// now exec check in catalogue :
						BOOL_T bOK=TRUE;
						if( gCCDParams.m_bCheckIfStar && 
							 (in.pPipeline)->m_pAsasTransform->m_bTransformOK ){
							// check if star :
							mystring szStarDesc;
							BOOL_T bIsKnownStar = IsStar( *pEvt, NULL, szStarDesc );
							if( bIsKnownStar ){
								pEvt->m_EventType = EVENT_TYPE_STAR;
								pEvt->m_szSatName = szStarDesc;
								if( gCCDParams.m_bRejectStars ){
									if( gCCDParams.GetMC() && samplesStat.size() ){
										samplesStat.back().nStarRej++;
									}
								}
								bOK=FALSE;
							}
						}
						if(bOK){
							if( pEvt->IsIdentified() ){
								finallist1.push_back( *pEvt );
							}
						}
						nVerified++;
					}
				}
			}
		}
	}

	//if( m_pPipeline ){
	//	m_pPipeline->DumpFinalEvents( &(m_pPipeline->m_FinalEventsLog), finallist1 );
	//}


	if(gCCDParams.m_eCheckIfMorePoint!=eAfterCurrentFrameOnly){
		// verification of number of events - rejection if too many required 
		// also at this stage :
		CCDMatrix::RejectIfMoreThen( eventToVerify, in.pPipeline, TRUE );
	}

	_TRACE_PRINTF_5("----------------------------------------------\n");

	// eventToVerify.clear();
	// eventToVerify = verifiedEvents;
	return nVerified;
}

void CCD_Analyser::SetNextFrameValues( CCDEventList& eventToVerify, CPixelAnalyseIn& in )
{
	CCDEventList::iterator pEvt;
	LONG_T nextClusterSum[MAX_PREV_SUMS];
	LONG_T nVerified=0;


	CCDEventList finallist1;
	for(pEvt = eventToVerify.begin();pEvt!=eventToVerify.end();pEvt++){
		// starting from 1 , at 0 there is frame to be confirmed :
		for(int next=1;next<in.PrevMatrixPtrCnt;next++){
			double next_x,next_y;

			CCDDataResults::CalcStarPositionAuto( pEvt->m_MaxPoint.x, pEvt->m_MaxPoint.y,
															  in.PrevFramesTime[0], in.PrevFramesTime[next], next ,
															  next_x, next_y,
															  (in.pCamCfg->m_CCDParams).m_RotCenterX,
															  (in.pCamCfg->m_CCDParams).m_RotCenterY,
//															  (in.pCamCfg->m_CCDParams).m_FrameDXPerSec,(in.pCamCfg->m_CCDParams).m_FrameDYPerSec,
															 (in.pCamCfg->m_CCDParams).m_FrameDX,(in.pCamCfg->m_CCDParams).m_FrameDY,
															  (in.pCamCfg->m_CCDParams).m_bUseRotInAverageOfPrevN,
															 (in.pCamCfg->m_CCDParams).m_RotValueDAlfa,
															  in.pCCDInfo, (in.pCamCfg->m_CCDParams).m_RotValueDAlfaPerSec,
															  &((in.pPipeline)->m_PipelineCfg) );



			if(pEvt->m_bIdentified){	


				if((pEvt->m_PixelAnalResults).eventType==eFlash){
					register int nextFrameX = my_round(next_x);
					register int nextFrameY = my_round(next_y);

					int val = GetNextFrameValue( (CCDMatrix*)in.PrevMatrixPtr[next], nextFrameX, 
														  nextFrameY, &(*pEvt), next, in );
					(pEvt->m_PixelAnalResults).laplaceOnNext = val;
				}
			}
		}
	}
}


BOOL_T CCD_Analyser::AnalysePrevMaxInBuffer(CCDPipeline& ccd_pipeline)
{
	Initialize();
	LONG_T pipeline_size = ccd_pipeline.GetPipelineSize();
	LONG_T pipeline_size_minus_1 = pipeline_size-1;
	if(ccd_pipeline.GetCount()!=pipeline_size)
		return FALSE;

	BOOL_T bRet = FALSE;
	LONG_T neighb_list[MAX_CLUSTER_SIZE];
	LONG_T ncnt;

	if(gCCDParams.m_bAnalyseMaxFromPrevInPipeline){
		clock_t t1 = clock();	
		cCCD& newFrame = ccd_pipeline.GetCurrent();
		long size = newFrame.GetCount();
		if(size<1){
			return FALSE;
		}
		long xSize = newFrame[0].GetXSize();
      long ySize = newFrame[0].GetYSize();
		long yUpperLimit = (ySize-gCCDParams.m_nIgnoreEdge);
		long xUpperLimit = (xSize-gCCDParams.m_nIgnoreEdge);		

		for(int i=0;i<size;i++){
			CCDMatrix& Matrix = newFrame[i];			
			ELEM_TYPE* p_data = Matrix.get_data_buffer();
			Table2D<ELEM_TYPE>** PrevMatrixPtr = GetPrevMatrixPtrs( i, ccd_pipeline );
			CPixelList pixel_list(xSize*ySize);
			
			long y_pos = (gCCDParams.m_nIgnoreEdge-1)*xSize+gCCDParams.m_nIgnoreEdge;
			for(int y=gCCDParams.m_nIgnoreEdge;y<yUpperLimit;y++){
				y_pos = y_pos + xSize;
				long pos = y_pos;
				for(int x=gCCDParams.m_nIgnoreEdge;x<xUpperLimit;x++,pos++){
					if(!pixel_list.CheckPixel(pos)){						
						// calculate sums only for pixels not included in
                  // any previously found cluster
						GetAnalNeighbWithSelf(x,y,pos,xSize,ySize,neighb_list,ncnt);
						LONGLONG_T newSum = CalcSumOpt( neighb_list, ncnt, p_data );
						LONGLONG_T maxSum = CalcMaxSumOnPrev( neighb_list, ncnt, 
                                                        PrevMatrixPtr, pipeline_size_minus_1 );

						// temporary :
						// printf("x=%d, y=%d, NewSum=%f, MaxSum=%f\n",x,y,(double)newSum,(double)maxSum);


						long x_y_val = p_data[pos];
//#ifdef _DO_NOT_ALLOW_NEGATIVE_
//               Assert(x_y_val>=0,"Negative value found in matrix analysis = %d at position (%d,%d)",x_y_val,x,y);			
//#endif

						if(x_y_val>=gCCDParams.m_MaxAllowedVal){
							// skiping pixel exceeding maximum allowed values in
							continue;
						}
						if( CheckCondition(newSum,maxSum) ){
							CLongList cluster,star_cluster;

							LONG_T x0,y0;
							double r_max;
							FindClusterAboveTreshold( Matrix, x, y, pos, x0, 
                                               y0, r_max, pixel_list, cluster );
	
							if(cluster.size()<gCCDParams.m_MinClusterSize){
							 	MYTRACE4(gCCDTrace,"Potential event at x=" << x << ",y=" << y << ", but cluster size=" << cluster.size() << "<" << gCCDParams.m_MinClusterSize);
							 	continue;	
							}
							LONG_T pos0 = x0 + y0*xSize;
							CLongList ClusterWithMore;
							GetClusterWithPointsAround(cluster,ClusterWithMore,
																xSize,ySize,gCCDParams.m_nPixelsAroundToConfirm);
							if(!ConfirmEvent_InCluster( i , ccd_pipeline, Matrix, 
                                                   pos0, xSize, ySize, ClusterWithMore )){
								MYTRACE4(gCCDTrace,"Potential event at x=" << x0 << ",y=" << y0 << ", but rejected by prev sum rejection procedure !");
								continue;									
							}
							bRet = TRUE;

							Matrix.AddFoundEvent( x0, y0, cluster );
							Matrix.SetInteresting();
							newFrame.SetInteresting();
							MYTRACE2(gCCDTrace,"ACCEPTED : Potential event found x=" << x0 << ", y="<<y0<<" ccd_index="<<i<<", sum="<<newSum<<", max_sum="<<maxSum);
						}	
					}
				}
			}
		}
		clock_t t2 = clock();
   	mystring msg = get_clock_in_sec_string( t2-t1 );

	   MYTRACE1(gCCDTrace,"MaxOnPrevInPipeline Algorithm analysis took : " << msg );
	}

	return bRet;
}

BOOL_T CCD_Analyser::AnalysePrevMaxNew(CCDPipeline& ccd_pipeline)
{
	if(ccd_pipeline.GetCount()!=ccd_pipeline.GetPipelineSize())
		return FALSE;

	BOOL_T bRet = FALSE;
	if(gCCDParams.m_bAnalyseMaxFromPrevFrames){
		clock_t t1 = clock();	
		Assert(gCCDParams.m_FramesBack<ccd_pipeline.GetCount(),"Analysis can not use %d frames, because pipeline size is %d",gCCDParams.m_FramesBack,ccd_pipeline.GetCount());
		cCCD& newFrame = ccd_pipeline.GetCurrent();
		long size = newFrame.GetCount();
		cCCD* pMaxFrame = ccd_pipeline.GetMaxFrame();
		for(int i=0;i<size;i++){
			CCDMatrix& Matrix = newFrame[i];			
			CCDMatrix& MaxMatrix = (*pMaxFrame)[i];
			ELEM_TYPE* p_max_data = MaxMatrix.get_data_buffer();
			ELEM_TYPE* p_data = Matrix.get_data_buffer();
			long xSize = Matrix.GetXSize();
			long ySize = Matrix.GetYSize();
			long yUpperLimit = (ySize-gCCDParams.m_nIgnoreEdge);
         long xUpperLimit = (xSize-gCCDParams.m_nIgnoreEdge);
			CPixelList pixel_list(xSize*ySize);
			CLongList neighb_list,outer_neighb;
			for(int y=gCCDParams.m_nIgnoreEdge;y<yUpperLimit;y++){
				for(int x=gCCDParams.m_nIgnoreEdge;x<xUpperLimit;x++){
					long pos = x + y*xSize;
					if(!pixel_list.CheckPixel(pos)){						
						// calculate sums only for pixels not included in
                  // any previous sum
						GetAnalNeighb(x,y,xSize,ySize,neighb_list,outer_neighb);
						LONGLONG_T newSum = CalcSum( neighb_list, p_data );
						LONGLONG_T maxSum = p_max_data[pos]; // already sum stored
						long x_y_val = p_data[pos];
						long max_x_y_val = p_max_data[pos];
//#ifdef _DO_NOT_ALLOW_NEGATIVE_
//              Assert(x_y_val>=0,"Negative value found in matrix analysis = %d at position (%d,%d)",x_y_val,x,y);			
//#endif

						if(x_y_val>=gCCDParams.m_MaxAllowedVal || max_x_y_val>=gCCDParams.m_MaxAllowedVal){
							// skiping pixel exceeding maximum allowed values in
							continue;
						}
						if( CheckCondition( newSum, maxSum ) ){
							CLongList cluster,star_cluster;
							// FindMaxClusterNew( Matrix, MaxMatrix, outer_neighb ,
                     //                   cluster, pos, pixel_list );

							LONG_T x0,y0;
							double r_max;
							FindClusterAboveTreshold( Matrix, x, y, pos, x0, 
                                               y0, r_max, pixel_list, cluster );
	
							if(cluster.size()<gCCDParams.m_MinClusterSize){
							 	MYTRACE4(gCCDTrace,"Potential event at x=" << x << ",y=" << y << ", but cluster size=" << cluster.size() << "<" << gCCDParams.m_MinClusterSize);
							 	continue;	
							}
							LONG_T pos0 = x0 + y0*xSize;
							CLongList ClusterWithMore;
							GetClusterWithPointsAround(cluster,ClusterWithMore,
																xSize,ySize,gCCDParams.m_nPixelsAroundToConfirm);
							if(!ConfirmEvent_InCluster( i , ccd_pipeline, Matrix, 
                                                   pos0, xSize, ySize, ClusterWithMore )){
								MYTRACE4(gCCDTrace,"Potential event at x=" << x0 << ",y=" << y0 << ", but rejected by prev sum rejection procedure !");
								continue;									
							}
							/*if(!ConfirmEvent_MaxPrevSums( i , ccd_pipeline, Matrix, 
                                                   pos0, xSize, ySize, r_max+3 )){
								MYTRACE4(gCCDTrace,"Potential event at x=" << x0 << ",y=" << y0 << ", but rejected by prev sum rejection procedure !");
								continue;									
							}*/
							/*if(!ConfirmEvent_MaxAlgorithm( Matrix, MaxMatrix, pos, xSize, ySize )){
								MYTRACE4(gCCDTrace,"Potential event at x=" << x << ",y=" << y << ", but rejected by confirmation procedure !");
								continue;									
							}*/
							bRet = TRUE;
							// long x0,y0;
							// CalcCenterOfHit( p_data, cluster, xSize, x0, y0 );
							// Assert(x0>=0 && x0<xSize,"X=%d coordinate of hit center out of range [%d-%d]",x0,0,xSize);
							// Assert(y0>=0 && y0<ySize,"X=%d coordinate of hit center out of range [%d-%d]",y0,0,ySize);
							Matrix.AddFoundEvent( x0, y0, cluster );
							Matrix.SetInteresting();
							newFrame.SetInteresting();
							MYTRACE2(gCCDTrace,"ACCEPTED : Potential event found x=" << x0 << ", y="<<y0<<" ccd_index="<<i<<", sum="<<newSum<<", max_sum="<<maxSum);
						}	
					}
				}
			}
		}
		clock_t t2 = clock();
   	mystring msg = get_clock_in_sec_string( t2-t1 );

	   MYTRACE1(gCCDTrace,"SumMax Algorithm analysis took : " << msg );
	}

	return bRet;
}



// still remains to be changed !!!!!!!!!!!!!!!!!!!!!!!!!!
// it should check already existing events (if any) and 
// require that only those accepted by latest criteria remain on the list 
void CCD_Analyser::UpdateEventReport( CCDPipeline& ccd_pipeline, LONG_T idx )
{
	PrepareEventReport( ccd_pipeline, idx );
}

void CCD_Analyser::PrepareEventReport( CCDPipeline& ccd_pipeline, LONG_T idx )
{
	cCCD& newFrame = ccd_pipeline.GetCurrent();
   long size = newFrame.GetCount();

	m_AllEvents.clear();
	long total_size = 0;
	int i;
	
	for(i=0;i<size;i++){
		CCDMatrix& Matrix = newFrame[i];
 	   total_size += Matrix.GetFoundEvents().size()+Matrix.GetGenEvents().size();
	}
	if(m_AllEvents.capacity()<total_size){
		printf("reserving more memory for %d events...\n",(total_size+10));
		m_AllEvents.reserve( total_size+10  );
	}
   for(i=0;i<size;i++){
		CCDMatrix& Matrix = newFrame[i];
		// CCDEventList mEvents( Matrix.GetFoundEvents().size()+20 );
		Matrix.CompileEventReport( m_AllEvents, idx, TRUE );
		// m_AllEvents += mEvents;
	}
	
	/*if(gCCDParams.m_bDumpAllEvents){
		m_AllRunEvents += m_AllEvents;
		if(m_AllRunEvents.size()>=gCCDParams.m_DumpEventsFreq){		
			mystring mode="a";
			if(gFirstDump)	
				mode="w";
			m_AllRunEvents.DumpEventReport( NULL, m_AllRunEvents, -1, FALSE, TRUE, mode.c_str() );
			m_AllRunEvents.clear();
			gFirstDump=FALSE;
		}
	}*/
}

int CCD_Analyser::InitSatList( CCDPipeline& ccd_pipeline, time_t ut_time )
{
	int ret=0;
	if( gCCDParams.m_bCheckIfSatelite && m_pCurrSatList ){
		CCDMatrix& newFrame = ccd_pipeline.GetCurrent()[0];

		
		if( ut_time == 0 ){
			if( gCCDParams.m_bReadFirstLevelInfoFromLog ){
				if( gCCDParams.m_bFromLogFile_WithFrames ){
					ut_time = (time_t)newFrame.getObsTime();
				}else{
					ut_time = ccd_pipeline.m_FrameUnixTime;
				}
			}else{
				ut_time = (time_t)newFrame.getObsTime();
			}
		}
	
		int all_count = CSatInfo::GetSatInfo( *m_pAllSatList, ut_time, FALSE );

		// ret = CSatInfo::GetSatInfo( *m_pCurrSatList, ut_time, TRUE );
		ret = CSatInfo::GetVisibleOnly( *m_pAllSatList, *m_pCurrSatList );
		printf("SAT-DB on frame %d all satellites# = %d, visible# = %d\n",
						ccd_pipeline.GetDayFrameCounter(),all_count,ret);
	}
	return ret;
}


BOOL_T CCD_Analyser::AnalyseNewFrameOpt(CCDPipeline& ccd_pipeline, BOOL_T bReport, LONG_T idx)
{
	if(ccd_pipeline.GetCount()!=ccd_pipeline.GetPipelineSize())
		return FALSE;
	
	ClearState();
	BOOL_T bRet=TRUE;

	InitSatList( ccd_pipeline );

	// test - calling optimized function !
	if(gCCDParams.m_bSuperOptimized)
		bRet = StartOptimizedAnalyse( ccd_pipeline );
	else
		bRet = FullNewFrameAnalyse( ccd_pipeline );

	if( gCCDParams.m_bAnalyzeSumOfPrevNFrames ){
   	if( ccd_pipeline.m_nTotalCollectedSumFrames >= 2 &&
      	 (ccd_pipeline.m_nFramesCollectedForSum == gCCDParams.m_bKeepSumOfPrevNFrames) ){
      	printf("Frame %d, analysing sum of previous %d frames\n", ccd_pipeline.GetFrameIndex(),gCCDParams.m_bKeepSumOfPrevNFrames);fflush(0);
      	AnalyseSumOfPrevNFrames( ccd_pipeline );
   	}
	}

	if( gCCDParams.m_nCompareToOldFreqInSec>0 ){
		if( ( ccd_pipeline.m_FrameUnixTime - ccd_pipeline.m_PrevSumOfFramesTime)>gCCDParams.m_nCompareToOldFreqInSec && 
			 ( ccd_pipeline.m_nFramesCollectedForSum==gCCDParams.m_bKeepSumOfPrevNFrames ) && 
			 strlen( ccd_pipeline.m_szPrevSumOfNFileName.c_str() ) ){
			printf("OLD_ANAL : analysing average frame to frame %s sec old ...",gCCDParams.m_nCompareToOldFreqInSec);
			AnalyseCompareToOldFrame( ccd_pipeline );
		}
	}

	// for newly found events equatorial coordinates are determined now :
	CalcAdditionalInfoForEvents( ccd_pipeline );
	// CalcAstroCoorinatesOfNewEvents( ccd_pipeline );


	if(gCCDParams.GetMC()){
		// it must be here - after PrepareEventReport is called 
		// not identified events are removed from list :
		ccd_pipeline.UpdateGenEventsList();
	}


	if(bRet || gCCDParams.GetMC()){
		// only in case events found or Monte Carlo (here event is put on picture)
		PROFILER_START
		PrepareEventReport( ccd_pipeline, idx );
		PROFILER_END("Preparation of event report took :")
	}		

	// adding events to global list :
	// here adds found events to global list 
	ccd_pipeline.AddFoundEventsToList();

	// now verify tracks
	// tracks verification is performed here because I wanted to skip generated 
	// events in tracks fit - so I must be after Report Compilation :
	// ( in real analysis can be uncommented in function : StartOptimizedAnalyse
	// and removed here :
	//if( !gCCDParams.m_bCCDDouble ){
	//   VerifyTracks( ccd_pipeline );
	//}


	return bRet;
}

BOOL_T CCD_Analyser::AnalyseNewFrame(CCDPipeline& ccd_pipeline, BOOL_T bReport, LONG_T idx){
	if(ccd_pipeline.GetCount()!=ccd_pipeline.GetPipelineSize())
		return FALSE;
	
	m_AllEvents.clear();

	BOOL_T bRetMax = TRUE;
	if(gCCDParams.m_bAnalyseMaxFromPrevFrames){
		bRetMax = AnalysePrevMaxNew( ccd_pipeline );
		if(bReport){
			PrepareEventReport( ccd_pipeline, idx );
		}		
	}

	BOOL_T bRetMaxOnPrevInPipeline=TRUE;
	if(gCCDParams.m_bAnalyseMaxFromPrevInPipeline){
		bRetMaxOnPrevInPipeline = AnalysePrevMaxInBuffer( ccd_pipeline );
		if(bReport){
			UpdateEventReport( ccd_pipeline, idx );
		}
	}

	BOOL_T bRetSum = TRUE;
	if(gCCDParams.m_bAnalyseSumAround){
		bRetSum = AnalyseSumOfNeighbours(ccd_pipeline);
	}


	return (bRetMax && bRetSum && bRetMaxOnPrevInPipeline);
}



void CCD_Analyser::CalcCenterOfHit( ELEM_TYPE* p_data, CLongList& cluster,
                   				      long xSize, long& x0, long& y0 )
{
	double x_sum=0,y_sum=0,hit_sum=0;	
	CLongList::iterator pPoint;
#ifdef _DEBUG
	CLongList used;
#endif
   for(pPoint=cluster.begin();pPoint!=cluster.end();pPoint++){	
		long pos = (*pPoint);

#ifdef _DEBUG
		Assert(!used.FindPoint(pos),"Point pos=%d already used",pos);
		used.Add(pos);		
#endif

		long x = (pos % xSize);
		long y = (pos / xSize);

		x_sum += p_data[pos]*x;
		y_sum += p_data[pos]*y;
		hit_sum += p_data[pos];
	}
	x0 = (long)(x_sum/hit_sum);
	y0 = (long)(y_sum/hit_sum);
}

int CCD_Analyser::CalcCenterOfHitRealOpt( ELEM_TYPE* p_data,
                                  LONG_T* cluster,LONG_T cluster_cnt,
                                  long xSize, double& x0, double& y0 )
{
	double x_sum=0,y_sum=0,hit_sum=0;	
	for(register int i=0;i<cluster_cnt;i++){
		double x = (cluster[i] % xSize);
      double y = (cluster[i] / xSize);
		
		x_sum += p_data[ cluster[i] ]*x;
      y_sum += p_data[ cluster[i] ]*y;
		hit_sum += p_data[ cluster[i] ];
	}

	x0 = 0;
	y0 = 0;
	if(hit_sum!=0){
		x0 = (x_sum/hit_sum);
		y0 = (y_sum/hit_sum);
	}
	return ((int)hit_sum);
}

int CCD_Analyser::CalcCenterOfHitRealOpt( BIG_ELEM_TYPE** p_data,
                                  LONG_T* cluster,LONG_T cluster_cnt,
                                  long xSize, double& x0, double& y0 )
{
	double x_sum=0,y_sum=0,hit_sum=0;	
	for(register int i=0;i<cluster_cnt;i++){
		int x = (cluster[i] % xSize);
      int y = (cluster[i] / xSize);
		
		x_sum += (p_data[y][x])*x;
      y_sum += p_data[y][x]*y;
		hit_sum += p_data[y][x];
	}

	x0 = 0;
	y0 = 0;
	if(hit_sum!=0){
		x0 = (x_sum/hit_sum);
		y0 = (y_sum/hit_sum);
	}
	return (int)hit_sum;
}

int CCD_Analyser::CalcCenterOfHitRealOpt( ELEM_TYPE** p_data,
                                  LONG_T* cluster,LONG_T cluster_cnt,
                                  long xSize, double& x0, double& y0 )
{
	double x_sum=0,y_sum=0,hit_sum=0;	
	for(register int i=0;i<cluster_cnt;i++){
		int x = (cluster[i] % xSize);
      int y = (cluster[i] / xSize);
		
		x_sum += (p_data[y][x])*x;
      y_sum += p_data[y][x]*y;
		hit_sum += p_data[y][x];
	}

	x0 = 0;
	y0 = 0;
	if(hit_sum!=0){
		x0 = (x_sum/hit_sum);
		y0 = (y_sum/hit_sum);
	}
	return (int)hit_sum;
}


LONGLONG_T CCD_Analyser::CalcCenterOfHitOpt( ELEM_TYPE* p_data, 
            							   LONG_T* cluster,LONG_T cluster_cnt,                       
                   				      long xSize, long& x0, long& y0 )
{
	double x_sum=0,y_sum=0,hit_sum=0;	
	for(int i=0;i<cluster_cnt;i++){
		long x = (cluster[i] % xSize);
      long y = (cluster[i] / xSize);
		
		x_sum += p_data[ cluster[i] ]*x;
      y_sum += p_data[ cluster[i] ]*y;
		hit_sum += p_data[ cluster[i] ];
	}

	x0 = 0;
   y0 = 0;
   if(hit_sum!=0){
		x0 = (long)(x_sum/hit_sum);
		y0 = (long)(y_sum/hit_sum);
	}
	return ((LONGLONG_T)hit_sum);
}

void CCD_Analyser::CalcCenterOfHit( ELEM_TYPE* p_data, CPointList& cluster,
												long xSize, long& x0, long& y0 )
{
	double x_sum=0,y_sum=0,hit_sum=0;	
	CPointList::iterator pPoint;
   for(pPoint=cluster.begin();pPoint!=cluster.end();pPoint++){	
		LONG_T x_long = (LONG_T)pPoint->x;
		LONG_T y_long = (LONG_T)pPoint->y;
		LONG_T pos = x_long + y_long*xSize;
		x_sum += p_data[pos]*(pPoint->x);
		y_sum += p_data[pos]*(pPoint->y);
		hit_sum += p_data[pos];
	}
	x0 = (long)(x_sum/hit_sum);
	y0 = (long)(y_sum/hit_sum);
}

long CCD_Analyser::GetNeighbours(long x,long y,long xSize,long ySize,
                                 CPointList& Neighbours)
{
	Neighbours.clear();
	
	if(y-1>=0)
		Neighbours.Add(CPoint(x,y-1));

	if(x-1>=0)
		Neighbours.Add(CPoint(x-1,y));

	if(x+1<xSize)
		Neighbours.Add(CPoint(x+1,y));

	if(y+1<ySize)
		Neighbours.Add(CPoint(x,y+1));

	return Neighbours.size(); 
}

void CCD_Analyser::FindSumCluster(CCDMatrix& newMatrix,CCDPipeline& ccd_pipeline,
		                     		 long x0,long y0,CLongList& cluster,CPixelList& pixels)
{
	ELEM_TYPE *p_data,*p_prev;
	CLongList newPoints,Neighbours,allUsed;		
	long xSize = newMatrix.GetXSize();
	long ySize = newMatrix.GetYSize();
	LONG_T pos = x0+y0*xSize;
	LONG_T nPixelTreshold = gCCDParams.m_TresholdPerPixel;
	LONG_T MatrixIdx = newMatrix.GetIndex();

	
	cluster.RemoveAll();
	newPoints.Add( pos );
	allUsed.Add( pos );
	p_data = newMatrix.get_data_buffer();

	while(newPoints.GetCount()){
		long& last = newPoints[0];
		cluster.Add(last);
      newPoints.erase(newPoints.begin());
		long x = ( last % xSize );
		long y =  ( last / xSize );
      GetNeighbPositions(x,y,xSize,ySize,Neighbours);

		CLongList::iterator pPoint;
	   for(pPoint=Neighbours.begin();pPoint!=Neighbours.end();pPoint++){
			LONG_T _pos = (*pPoint);
			if(!allUsed.FindPoint( _pos )){		
				cCCD* pFrame;
				LONG_T num=0;
				BOOL_T bAccepted=TRUE;
				ccd_pipeline.GetCurr();
            for(pFrame = ccd_pipeline.GetPrev();pFrame!=NULL && num<gCCDParams.m_FramesBack
                ;pFrame= ccd_pipeline.GetPrev()){
					CCDMatrix& prevMatrix = (*pFrame)[MatrixIdx];
					p_prev = prevMatrix.get_data_buffer();
					if((p_data[_pos]-p_prev[_pos])<=nPixelTreshold){
						bAccepted=FALSE;
						break;
					}
					num++;
				}
				if(bAccepted){
					newPoints.push_back( _pos );
					allUsed.push_back( _pos );
					pixels.HitPixel( _pos );
				}
			}
		}		 
	}
}


void CCD_Analyser::FindMaxCluster(CCDMatrix& ccd_matrix,CCDMatrix& ccd_max,
                                  long x0,long y0,
                                  CPointList& cluster)
{	
	long x_n,y_n,sum,sum_max;
	CPointList newPoints,Neighbours,allUsed;
	ELEM_TYPE *p_data,*p_max_data;
	long xSize = ccd_matrix.GetXSize();
	long ySize = ccd_matrix.GetYSize();

	cluster.clear();
	CPoint tmp(x0,y0);
	newPoints.push_back(tmp);	
	allUsed.push_back(tmp);
	p_data = ccd_matrix.get_data_buffer();
	p_max_data = ccd_max.get_data_buffer();

	while(newPoints.size()){
		 CPoint& Last = newPoints[0];
		 cluster.push_back(Last);
		 LONG_T x = (LONG_T)Last.x;
       LONG_T y = (LONG_T)Last.y;
		 newPoints.erase(newPoints.begin());
		 GetNeighbours(x,y,xSize,ySize,Neighbours);
		 
		 CPointList::iterator pPoint;
	    for(pPoint=Neighbours.begin();pPoint!=Neighbours.end();pPoint++){
			 if(!allUsed.FindPoint((LONG_T)pPoint->x,(LONG_T)pPoint->y)){
				 long sum = CalcSum( p_data, (LONG_T)pPoint->x, (LONG_T)pPoint->y ,xSize ,ySize );
   		    long max_sum = CalcSum( p_max_data , (LONG_T)pPoint->x, (LONG_T)pPoint->y ,xSize, ySize );
				 if(sum==0 && max_sum>0){
                // pixel skiped - probably after coregister value is 0					 
				    continue;
				 }
				 if(sum - max_sum > gCCDParams.m_SumTresholdForNewFrame){
				     newPoints.push_back(CPoint( pPoint->x, pPoint->y ));	     				 	      
					  allUsed.push_back(CPoint( pPoint->x, pPoint->y ));
				 }
			 }
		 }		
	}
}

BOOL_T CCD_Analyser::ConfirmEvent_MaxPrevSums( LONG_T MatrixIdx, CCDPipeline& ccd_pipeline,
                                               CCDMatrix& ccd_matrix, LONG_T CorePos, 
                                               LONG_T xSize, LONG_T ySize, double r0 )
{
	CCDMatrix* pMaxMatrix;
	CCDPipelineIterator i(&ccd_pipeline);
	LONG_T x = ( CorePos % xSize );
	LONG_T y = ( CorePos / xSize );
	double nConfRedial = r0;
	
	LONGLONG_T MaxSum=-1,newSum=-1;
	LONGLONG_T MaxIdx=-1;
	for(;!i.curr();i++){
		CCDMatrix* pMatrix;
		pMatrix = i->GetMatrix( MatrixIdx );
		LONGLONG_T Sum = pMatrix->CalcSumAround( x, y , (LONG_T)nConfRedial );
		if(Sum>MaxSum){
			MaxSum = Sum;
			pMaxMatrix = pMatrix;	
		}
	}
	if(MaxSum>0){
		LONGLONG_T nConfTreshold = (LONGLONG_T)(PI*nConfRedial*nConfRedial)*gCCDParams.m_ConfTresholdPerPixel;
		newSum = ccd_matrix.CalcSumAround( x, y , (LONG_T)nConfRedial );
		if(newSum - MaxSum > nConfTreshold){
			MYTRACE4(gCCDTrace,"Event at (" << x << "," << y << ") accepted MaxSum=" << MaxSum << ", newSum=" << newSum );
			return TRUE;
		}		
	}else{
		// error - claim and trace - reject event
	}
	MYTRACE4(gCCDTrace,"REJECTED by confirmation ConfirmEvent_MaxPrevSums procedure,  Event at (" << x << "," << y << ") MaxSum=" << MaxSum << ", newSum=" << newSum );	
	return FALSE;
}

BOOL_T CCD_Analyser::ConfirmEvent_InCluster( LONG_T MatrixIdx, CCDPipeline& ccd_pipeline,
	                                          CCDMatrix& ccd_matrix, LONG_T CorePos, 
                                             LONG_T xSize, LONG_T ySize, 
	                                          CLongList& cluster )
{
	CCDMatrix* pMaxMatrix;
	CCDPipelineIterator i(&ccd_pipeline);
	LONG_T x = ( CorePos % xSize );
	LONG_T y = ( CorePos / xSize );
	
	LONGLONG_T MaxSum=-1,newSum=-1;
	LONGLONG_T MaxIdx=-1;
	for(;!i.curr();i++){
		CCDMatrix* pMatrix;
		pMatrix = i->GetMatrix( MatrixIdx );
		LONGLONG_T Sum = pMatrix->CalcSum( cluster );
		if(Sum>MaxSum){
			MaxSum = Sum;
			pMaxMatrix = pMatrix;	
		}
	}
	LONGLONG_T nConfTreshold = (cluster.size())*gCCDParams.m_ConfTresholdPerPixel;
	LONGLONG_T nMaxOnPrevAllowed = (cluster.size())*gCCDParams.m_ConfMaxPrevPerPixel;
	if(MaxSum>=0){
		newSum = ccd_matrix.CalcSum( cluster );
		if( CheckConfCondition(newSum,MaxSum,nConfTreshold,nMaxOnPrevAllowed) ){
			MYTRACE4(gCCDTrace,"Event at (" << x << "," << y << ") ACC_InCluster MaxSum=" << MaxSum << ", newSum=" << newSum 
		                  << ", points in cluster = " << cluster.size() << ", req treshold ="<< nConfTreshold );	 
			return TRUE;
		}		
	}else{
		// error - claim and trace - reject event
		mystring szTrace,szPoints;
		CLongList::iterator i;
		for(i=cluster.begin();i!=cluster.end();i++){
			long x1 = ( (*i) % xSize );
			long y1 = ( (*i) / xSize );
	 		szPoints << "(" << x1 << "," << y1 << "),";
		}
		szTrace << "Error in procedure - NEGATIVE MAX_SUM in cluster = " << MaxSum << ", cluster : " << szPoints << " - exiting ...";
		MYTRACE1(gCCDTrace,szTrace.c_str());
		exit(-1);
	}
	MYTRACE4(gCCDTrace,"REJECTED by conf ConfirmEvent_InCluster procedure,  Event at (" << x << "," << y << ") MaxSum=" << MaxSum << ", newSum=" << newSum
                      << ", points in cluster = " << cluster.size() << ", required treshold ="<< nConfTreshold
                      << ", req for prev<" << nMaxOnPrevAllowed);	
	return FALSE;
	
}

void CCD_Analyser::DumpCluster( LONG_T* cluster, LONG_T cluster_cnt, 
                                ELEM_TYPE* p_data, LONG_T xSize, LONG_T ySize,
                                const char* desc )
{
	mystring szTrace,szPoints;
	for(int k=0;k<cluster_cnt;k++){
		long x1 = ( cluster[k] % xSize );
		long y1 = ( cluster[k] / xSize );
		szPoints << "(" << x1 << "," << y1 << ")=" << p_data[cluster[k]]  << ",";
	}	
	MYTRACE4(gCCDTrace,"CLUSTER_DUMP " << desc << ": " << szPoints);
}

BOOL_T CCD_Analyser::ConfirmEvent_InClusterOpt_RotCorrected( 
							  const CPixelAnalyseIn& in,
                       LONG_T* cluster, LONG_T cluster_cnt,
                       LONG_T& newClusterSum,LONG_T* prevClusterSum )
{
	Table2D<ELEM_TYPE>* pMaxMatrix=NULL;
	
	LONGLONG_T MaxSum=-1;
	LONGLONG_T MaxIdx=-1;
	newClusterSum = -1;

	// starting from the newest frame !	
	for(register long i=0;i<in.PrevMatrixPtrCnt;i++){	
		double frame_dx = -((i+1)*gCCDParams.m_FrameDX);
      double frame_dy = -((i+1)*gCCDParams.m_FrameDY);		

		LONGLONG_T Sum = (in.PrevMatrixPtr[i])->CalcSumRotCorrected( cluster,
									cluster_cnt, (long)frame_dx, (long)frame_dy ); 
		prevClusterSum[i] = Sum;
		if(Sum>MaxSum){
			MaxSum = Sum;
			pMaxMatrix = in.PrevMatrixPtr[i];
		}
	}

	LONGLONG_T nConfTreshold = cluster_cnt*gCCDParams.m_ConfTresholdPerPixel;
	LONGLONG_T nMaxOnPrevAllowed = cluster_cnt*gCCDParams.m_ConfMaxPrevPerPixel;
	if(MaxSum>=0){
		newClusterSum = CalcSumOpt( cluster, cluster_cnt, in.p_data );
		prevClusterSum[in.PrevMatrixPtrCnt] = newClusterSum;
		MYTRACE4(gCCDTrace,"Cluster check condition ( for " << cluster_cnt     
           << " pixels : " << GetClusterCheckDesc(newClusterSum,MaxSum,nConfTreshold,nMaxOnPrevAllowed) );
		if( CheckConfCondition(newClusterSum,MaxSum,nConfTreshold,nMaxOnPrevAllowed) ){
			MYTRACE4(gCCDTrace,"ACCEPTED Event at (" << in.x << "," << in.y << ")" );
			return TRUE;
		}		
	}else{
		// error - claim and trace - reject event
		DumpCluster( cluster, cluster_cnt , in.p_data, in.xSize, in.ySize, "" );
		exit(-1);
	}
	MYTRACE4(gCCDTrace,"REJECTED Event at (" << in.x << "," << in.y << ")" );
	if( gCCDParams.m_bDumpClusters ){
		DumpCluster( cluster, cluster_cnt , pMaxMatrix->get_data_buffer(), in.xSize, in.ySize, "MaxCluster :" );
	}

	return FALSE;

}


BOOL_T CCD_Analyser::ConfirmEvent_InClusterOpt( 
							   Table2D<ELEM_TYPE>** PrevMatrixPtr, LONG_T frames_back,
	                     ELEM_TYPE* p_data, LONG_T CorePos, 
                        LONG_T xSize, LONG_T ySize, 
                        LONG_T* cluster, LONG_T cluster_cnt )
{
	Table2D<ELEM_TYPE>* pMaxMatrix=NULL;
	LONG_T x = ( CorePos % xSize );
	LONG_T y = ( CorePos / xSize );
	
	LONGLONG_T MaxSum=-1,newSum=-1;
	LONGLONG_T MaxIdx=-1;

	for(int i=0;i<frames_back;i++){	
		LONGLONG_T Sum = (PrevMatrixPtr[i])->CalcSum( cluster, cluster_cnt );
		if(Sum>MaxSum){
			MaxSum = Sum;
			pMaxMatrix = PrevMatrixPtr[i];
		}
	}

	LONGLONG_T nConfTreshold = cluster_cnt*gCCDParams.m_ConfTresholdPerPixel;
	LONGLONG_T nMaxOnPrevAllowed = cluster_cnt*gCCDParams.m_ConfMaxPrevPerPixel;
	if(MaxSum>=0){
		newSum = CalcSumOpt( cluster, cluster_cnt, p_data );
		MYTRACE4(gCCDTrace,"Cluster check condition : " << GetClusterCheckDesc(newSum,MaxSum,nConfTreshold,nMaxOnPrevAllowed) );
		if( CheckConfCondition(newSum,MaxSum,nConfTreshold,nMaxOnPrevAllowed) ){
			MYTRACE4(gCCDTrace,"ACCEPTED Event at (" << x << "," << y << ")" );
			return TRUE;
		}		
	}else{
		// error - claim and trace - reject event
		DumpCluster( cluster, cluster_cnt , p_data, xSize, ySize, "" );
		exit(-1);
	}
	MYTRACE4(gCCDTrace,"REJECTED Event at (" << x << "," << y << ")" );
	if( gCCDParams.m_bDumpClusters ){
		DumpCluster( cluster, cluster_cnt , pMaxMatrix->get_data_buffer(), xSize, ySize, "MaxCluster :" );
	}

	return FALSE;
	
}



BOOL_T CCD_Analyser::ConfirmEvent_MaxAlgorithm( CCDMatrix& ccd_matrix,CCDMatrix& ccd_max,
                                                LONG_T CorePos, LONG_T xSize, LONG_T ySize )
{
	if(!gCCDParams.m_bConfirmReq || !gCCDParams.m_ConfRedial) 
		return TRUE;

	LONG_T x,y;
	ccd_matrix.GetXY_FromPos( CorePos, x ,y );
	LONG_T start_x = MAX(0,x-gCCDParams.m_ConfRedial);
	LONG_T start_y = MAX(0,y-gCCDParams.m_ConfRedial);	
	LONG_T end_x = MIN(xSize-1,x+gCCDParams.m_ConfRedial);
	LONG_T end_y = MIN(ySize-1,y+gCCDParams.m_ConfRedial);	
	LONGLONG_T MaxSum=0,NewSum=0;
	ELEM_TYPE* p_max_data = ccd_max.get_data_buffer();
	ELEM_TYPE* p_data = ccd_matrix.get_data_buffer();
	CLongList neighb_list,outer_neighb;

	for(LONG_T x0=start_x;x0<end_x;x0++){
		for(LONG_T y0=start_y;y0<end_y;y0++){
			LONG_T pos = x0 + y0*xSize;
			double r = sqrt((x0-x)*(x0-x)+(y0-y)*(y0-y));
			if(r <= gCCDParams.m_ConfRedial){
				MaxSum += p_max_data[pos];
				GetAnalNeighb(x0,y0,xSize,ySize,neighb_list,outer_neighb);
				NewSum += CalcSum( neighb_list, p_data );
			}
		}
	}
	
	if(NewSum - MaxSum >= gCCDParams.m_ConfTreshold)
		return TRUE;	

	return FALSE;	
}

void CCD_Analyser::FindClusterAboveTreshold(CCDMatrix& ccd_matrix, 
								 					     LONG_T x, LONG_T y,LONG_T CorePos,
                                            LONG_T& x0,LONG_T& y0, double& r_max,
                                            CPixelList& in_clusters, CLongList& cluster )
{
	CLongList Neighbours,allUsed,startPoints;
	ELEM_TYPE *p_data;
	LONG_T xSize = ccd_matrix.GetXSize();
	LONG_T ySize = ccd_matrix.GetYSize();

	cluster.clear();
	startPoints.Add(CorePos);
	// allUsed.Add(CorePos);
	p_data = ccd_matrix.get_data_buffer();
	while(startPoints.size()){
		 long pos = startPoints[0];
		 long x = ( pos % xSize );
	    long y =  ( pos / xSize );
		 cluster.Add( pos );
		 in_clusters.HitPixel( pos );
		 startPoints.erase(startPoints.begin());
		 GetNeighbours(pos,xSize,ySize,Neighbours,FALSE);
		
		 CLongList::iterator pPoint;
		 for(pPoint=Neighbours.begin();pPoint!=Neighbours.end();pPoint++){
           // if(allUsed.FindPoint(*pPoint))
			  //    continue;
			  if( in_clusters.CheckPixel( *pPoint ) )
               continue;			 
			  if(p_data[*pPoint]>gCCDParams.m_MaxNoiseLevel){
   	            startPoints.Add( *pPoint );						
			  }
           in_clusters.HitPixel( *pPoint );
			  //  allUsed.Add( *pPoint );
		 } 	
	}

	CalcCenterOfHit( ccd_matrix.get_data_buffer(), 	cluster, xSize, x0, y0);
	CLongList::iterator pP;
	r_max = -1;
	for(pP=cluster.begin();pP!=cluster.end();pP++){
		long x = ( (*pP) % xSize );
		long y =  ( (*pP) / xSize );
		double r = sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0) );
		if( r > r_max )
			r_max = r;				 
	}
}

/*void ADD_POINT( LONG_T* table, LONG_T& idx, LONG_T& pos )
{
	if(idx>=MAX_CLUSTER_SIZE){ 
		printf("Cluster size exceeded, exiting\n");
		for(int i=0;i<idx;i++){
			printf("%d,",pos);
		}
		exit(-1); 
	} 
	table[idx] = pos;
	idx++;	
}*/


double CCD_Analyser::GetMaxNoiseLevel( const CPixelAnalyseIn* pIn )
{
	if(gCCDParams.m_ClusterIfNSigmaAboveBackgr>-1000){
		double SigmaB,AverB;
		CCDPipeline::GetBackgroundStat( pIn->pPipeline, pIn->ccd_index, 
												  AverB, SigmaB );
		return SigmaB*gCCDParams.m_ClusterIfNSigmaAboveBackgr+AverB;
	}
	return gCCDParams.m_MaxNoiseLevel;
}

double CCD_Analyser::GetMaxNoiseLevelLaplace( CCDPipeline* pPipeline, int ccd_index,
															 double nSigmaAbove )
{
	if( nSigmaAbove<-100.00 ){
		// default is -1000.00 - so if default :
		nSigmaAbove = gCCDParams.m_ClusterIfNSigmaAboveBackgr;
	}
	if(nSigmaAbove>-1000){
		double SigmaB,AverB;
		CCDPipeline::GetBackgroundStat( pPipeline, ccd_index, AverB, SigmaB, 
										 		  gCCDParams.m_eLaplaceType );
		return SigmaB*nSigmaAbove+AverB;
	}
	return gCCDParams.m_MaxNoiseLevel;
}

double CCD_Analyser::GetMaxNoiseLevelLaplace( const CPixelAnalyseIn* pIn, double nSigmaAbove )
{
	return GetMaxNoiseLevelLaplace( pIn->pPipeline, pIn->ccd_index, nSigmaAbove );
}

BOOL_T CCD_Analyser::CheckEventShape( const CPixelAnalyseIn& in, CPixelAnalyseOut& out )
{
	BOOL_T bRet = TRUE;
	FindClusterAboveTresholdOpt2( in, out.m_PixelOut.x0_real,out.m_PixelOut.y0_real, out.cluster, 
											out.cluster_cnt, out.m_PixelOut.max_noise_level );

	out.m_PixelOut.max_redial = CCDFastCalc::FindMaxRedial( out.cluster, out.cluster_cnt, 
																			  out.m_PixelOut.x0_real, out.m_PixelOut.y0_real, in.xSize );
	out.m_PixelOut.pos0 = out.m_PixelOut.y0*in.xSize + out.m_PixelOut.x0;
	out.m_PixelOut.m_Sphericity = CCDFastCalc::CalcSphericity( out.m_PixelOut.max_redial, out.cluster_cnt );
//	printf("SPHERICITY = %f\n",spher);

	for(register int i=0;i<out.cluster_cnt;i++){
		(in.pPixelList)->ClearPixel( out.cluster[i] );
	}
	
	if(out.m_PixelOut.m_Sphericity < gCCDParams.m_CheckEventShape)
		bRet = FALSE;

	if(gCCDParams.m_bCheckCenterInCluster){
		if(find_value( out.cluster, out.cluster_cnt, out.m_PixelOut.pos0 )==NOT_FOUND){
			// center of cluster outside the cluster - this cannot be a star :
			bRet = FALSE;
		}
	}
	
	out.m_PixelOut.m_bRejectedByEventShape = TRUE;
	MYTRACE2(gCCDTrace,"Sphericity=" << mystring::double2string( out.m_PixelOut.m_Sphericity ) << ", rejected");		
	
	return bRet;
}

LONGLONG_T CCD_Analyser::FindClusterAboveTresholdOpt( const CPixelAnalyseIn* pIn,
														  ELEM_TYPE* p_data, 
								 					     LONG_T x, LONG_T y,LONG_T CorePos,
														  LONG_T xSize, LONG_T ySize,
                                            LONG_T& x0,LONG_T& y0,
                                            CPixelList& in_clusters,
                                            LONG_T* cluster,LONG_T& cluster_cnt)
{
	LONG_T Neighbours[100];
	static LONG_T startPoints[MAX_CLUSTER_SIZE];
	LONG_T* neighb_list = Neighbours;
	LONG_T neighb_cnt = 0;
	LONG_T all_cnt = 0;
	LONG_T startPoints_cnt = 0;
	LONG_T ncnt;
	LONGLONG_T sum=0;
	double maxNoiseLevel = GetMaxNoiseLevel( pIn );

	cluster_cnt = 0;	

	if(p_data[ CorePos ]>maxNoiseLevel){	
		ADD_POINT(startPoints,startPoints_cnt,CorePos);
	}

	while(startPoints_cnt){
		 long last = startPoints_cnt-1;
		 long pos = startPoints[last];			
		 long x = ( pos % xSize );
	    long y =  ( pos / xSize );
		 ADD_POINT(cluster,cluster_cnt,pos);
		 in_clusters.HitPixel( pos );
		 
		 startPoints_cnt--;
       GetAnalNeighbNoSelf(x,y,pos,xSize,ySize,neighb_list,ncnt);

		for(int i=0;i<ncnt && neighb_list[i]!=NEIGHB_LIST_END;i++){
			if( neighb_list[i]==pos || in_clusters.CheckPixel( neighb_list[i] ) )
               continue;
			if(p_data[ neighb_list[i] ]>maxNoiseLevel){
		      ADD_POINT(startPoints,startPoints_cnt,neighb_list[i]);
		   }
			in_clusters.HitPixel( neighb_list[i] );
		}
	}

	sum = CalcCenterOfHitOpt( p_data, cluster, cluster_cnt, xSize, x0, y0);	
	return sum;
}

int CCD_Analyser::FindClusterAboveTresholdOpt2( const CPixelAnalyseIn& in,
 	 	                                        double& x0,double& y0,
     			                                  LONG_T* cluster,LONG_T& cluster_cnt,
															 double& maxNoiseLevel)
{
	LONG_T Neighbours[100];
	static LONG_T startPoints[MAX_CLUSTER_SIZE];
	LONG_T* neighb_list = Neighbours;
	LONG_T neighb_cnt = 0;
	LONG_T all_cnt = 0;
	LONG_T startPoints_cnt = 0;
	LONG_T ncnt;
	int sum=0;
	maxNoiseLevel = GetMaxNoiseLevel( &in );

	cluster_cnt = 0;	

	if(in.p_data[ in.pos ]>maxNoiseLevel){	
		ADD_POINT_func(startPoints,startPoints_cnt,in.pos);
	}
	
   //                +
	// all pixels in +++ used for cluster finding :
   //                +
	if( in.x > 0 ){
		if(in.p_data[ in.pos-1 ]>maxNoiseLevel){
			ADD_POINT_func(startPoints,startPoints_cnt,in.pos-1);
		}
	}
	if( in.x < in.xSize-1 ){
		if(in.p_data[ in.pos+1 ]>maxNoiseLevel){
			ADD_POINT_func(startPoints,startPoints_cnt,in.pos+1);
		}
	}
	if( in.y > 0 ){	
		if(in.p_data[ in.pos-in.xSize ]>maxNoiseLevel){
	      ADD_POINT_func(startPoints,startPoints_cnt,in.pos-in.xSize);
		}
   }
	if( in.y < in.ySize-1 ){
		if(in.p_data[ in.pos+in.xSize ]>maxNoiseLevel){
	      ADD_POINT_func(startPoints,startPoints_cnt,in.pos+in.xSize);
		}
   }

	for(register int i=0;i<startPoints_cnt;i++){
		(in.pPixelList)->HitPixel( startPoints[i] );
	}
		

	while(startPoints_cnt){
		 register int last = startPoints_cnt-1;
		 register int pos = startPoints[last];			
		 register int x = ( pos % in.xSize );
       register int y =  ( pos / in.xSize );			
		 if(!ADD_POINT_func(cluster,cluster_cnt,pos)){
	       cluster_cnt=0;
			 return 0;
		 }
	 	 (in.pPixelList)->HitPixel( pos );			
			
		 
		 startPoints_cnt--;
       // GetAllNeighbNoSelf(x,y,pos,in.xSize,in.ySize,neighb_list,ncnt);
		 GetAnalNeighbNoSelf(x,y,pos,in.xSize,in.ySize,neighb_list,ncnt);

		for(register int i=0;i<ncnt && neighb_list[i]!=NEIGHB_LIST_END;i++){
			// if( neighb_list[i]==pos || find_value( cluster,cluster_cnt,neighb_list[i] )!=NOT_FOUND )
			if( neighb_list[i]==pos || (in.pPixelList)->CheckPixel( neighb_list[i] ) )
               continue;
			if(in.p_data[ neighb_list[i] ]>maxNoiseLevel){
				if(!ADD_POINT_func(startPoints,startPoints_cnt,neighb_list[i])){
					cluster_cnt=0;
					return 0;
				}
		   }
			(in.pPixelList)->HitPixel( neighb_list[i] );
		}
	}

	x0 = in.x;
	y0 = in.y;
	if( cluster_cnt )
		sum = CalcCenterOfHitRealOpt( in.p_data, cluster, cluster_cnt, in.xSize, x0, y0);	
	return sum;
}


int CCD_Analyser::FindMaxInCluster( ELEM_TYPE* p_data, int xSize,
                  	               LONG_T* cluster,LONG_T cluster_cnt,
                                    int& max_pos )
{
	int max_val = -100000;

	max_pos = -1;
	for(register int i=0;i<cluster_cnt;i++){

		if(p_data[cluster[i]]>max_val){
			max_pos = cluster[i];
			max_val = p_data[max_pos];
		}
	}
	return max_val;
}


BOOL_T CCD_Analyser::CheckForSN( const CPixelAnalyseIn& in, 
											CPixelAnalyseOut& out, 
											int TnForSN, int TvForSN, int MinPrevForSN )
{
	int x0=in.x;
	int y0=in.y;
	int max_value=-10000;
	for(int i=0;i<out.cluster_cnt;i++){
		int x = out.cluster[i] % in.xSize;
		int y = out.cluster[i] / in.xSize;
		if(in.p_curr_data_laplace[y][x]>max_value)
		{
			max_value = in.p_curr_data_laplace[y][x];
			x0 = x;
			y0 = y;
		}
	}
	if( max_value < TnForSN )
		return FALSE;

	int prevCount = MIN(((in.pPipeline)->GetCount()-1),gCCDParams.m_nMaxOfAverageOfPrevN);

	for(int x=(x0-1);x<=(x0+1);x++){
		for(int y=(y0-1);y<=(y0+1);y++){
			if (x>=0 && x<in.xSize && y>=0 && in.y<in.ySize ){
				int sum=0;
				for(register int f=1;f<=prevCount;f++){
					BIG_ELEM_TYPE** pLapPrev = ((CCDMatrix*)(in.PrevMatrixPtr[f]))->get_frame_laplace_fast();
					int prev_x,prev_y;
					GetPrevPixelPos( prev_x, prev_y, x, y, f, in );
					if(prev_x>=0 && prev_y>=0 && prev_x<in.xSize && prev_y<in.ySize){
						int m_val=pLapPrev[prev_y][prev_x];
						for(int xx=(prev_x-1);xx<=(prev_x+1);xx++){
							for(int yy=(prev_y-1);yy<=(prev_y+1);yy++){
								if( pLapPrev[yy][xx] > m_val ){
									m_val = pLapPrev[yy][xx];
								}
							}
						}
						sum += m_val;
					}
				}
				sum = ( sum / prevCount );
				if( sum > TvForSN || sum<gCCDParams.m_nMinLaplaceOnOther ){
					return FALSE;
				}
				if( x==x0 && y==y0 ){
					if( sum<MinPrevForSN ){
						return FALSE;
					}
				}
			}
		}
	}
	return TRUE;
}

BOOL_T CCD_Analyser::CheckForSN_OnSum( const CPixelAnalyseIn& in, 
											CPixelAnalyseOut& out, 
											int TnForSN, int TvForSN, int MinPrevForSN )
{
	int x0=in.x;
	int y0=in.y;
	int max_value=-10000;
	for(register int i=0;i<out.cluster_cnt;i++){
		int x = out.cluster[i] % in.xSize;
		int y = out.cluster[i] / in.xSize;
		if(in.p_curr_data_laplace[y][x]>max_value)
		{
			max_value = in.p_curr_data_laplace[y][x];
			x0 = x;
			y0 = y;
		}
	}
	if( max_value < TnForSN )
		return FALSE;


	int x0_1 = (x0-1);
	int x0p1 = (x0+1);
	int y0_1 = (y0-1);
	int y0p1 = (y0+1);
	for(register int x=x0_1;x<=x0p1;x++){
		for(register int y=y0_1;y<=y0p1;y++){
			if (x>=0 && x<in.xSize && y>=0 && in.y<in.ySize ){
				int sum = in.p_prev_lap_fast[y][x];
				if( sum > TvForSN || sum<gCCDParams.m_nMinLaplaceOnOther ){
					return FALSE;
				}
				if( x==x0 && y==y0 ){
					if( sum<MinPrevForSN ){
						return FALSE;
					}
				}
			}
		}
	}
	return TRUE;
}



BOOL_T CCD_Analyser::FindClustersOnNewAndPrev( const CPixelAnalyseIn& in,
															  CPixelAnalyseOut& out,
															  CPixelAnalyseIn& inForPrev,
															  CPixelAnalyseOut& outForPrev,
															  double& max_on_new,double& max_on_prev )
{
	double x0=0,y0=0,out_tresh=0;
																
	CCD_Analyser::FindClusterAboveTresholdOpt3( in, out, x0,y0, out_tresh, 
																((in.pPipeline)->GetPipelineCfg()).m_nNewLaplace, FALSE  );

	max_on_new = -10000.00;
	int max_pos = -1;
						
	max_on_new = CCD_Analyser::FindMaxInClusterBig( in.p_curr_data_laplace_normal, in.xSize,
															out.cluster, out.cluster_cnt, max_pos );
	if(max_pos>=0){						
		int star_x = (max_pos % in.xSize );
		int star_y = (max_pos / in.xSize );
		if( star_x < gCCDParams.m_nIgnoreEdgeLeft || star_x>=(in.xSize-gCCDParams.m_nIgnoreEdgeRight) || 
			 star_y < gCCDParams.m_nIgnoreEdgeBottom || star_y>=(in.ySize-gCCDParams.m_nIgnoreEdgeUp))
			return FALSE;
			
		// now find max val on previous frames :
		inForPrev.x = star_x;
		inForPrev.y = star_y;
		inForPrev.pos = max_pos;
		double x0OnPrev=0,y0OnPrev=0,out_tresh_for_prev;
		double noise_level = (out.m_PixelOut.maxAverageOfPrev*0.9);
		if(noise_level<0)
			noise_level = fabs(in.p_curr_data_laplace[star_y][star_x]*0.9);
		max_on_prev = FindClusterAboveTresholdOnPrevOpt3( inForPrev, 
														outForPrev.m_PixelOut.x0, 
														outForPrev.m_PixelOut.y0,
														outForPrev.cluster,outForPrev.cluster_cnt,
														out_tresh_for_prev, noise_level, FALSE );

		return TRUE;
	}
	return FALSE;
}						


BOOL_T CCD_Analyser::FindClustersOnNewAndPrevNew( const CPixelAnalyseIn& in,
															  CPixelAnalyseOut& out,
															  CPixelAnalyseIn& inForPrev,
															  CPixelAnalyseOut& outForPrev,
															  double& max_on_new,double& max_on_prev )
{
	double x0=0,y0=0,out_tresh=0;
																
	double noise_level = CCDDataResults::GetTreshold( gCCDParams.m_eLaplaceType,
								gCCDParams.m_nSigmaAboveMeanInRawCluster , 
								(in.pPipeline)->GetBackgroundMap() );

	CCD_Analyser::FindClusterAboveTresholdOpt3( in, out, x0,y0, out_tresh, 
															  noise_level, FALSE  );

	max_on_new = -10000.00;
	int max_pos = -1;
				
	max_on_prev = 0;				
	max_on_new = CCD_Analyser::FindMaxInClusterBig( in.p_curr_data_laplace_normal, in.xSize,
															out.cluster, out.cluster_cnt, max_pos );
	if(max_pos>=0){						
		int star_x = (max_pos % in.xSize );
		int star_y = (max_pos / in.xSize );
		if( star_x < gCCDParams.m_nIgnoreEdgeLeft || star_x>=(in.xSize-gCCDParams.m_nIgnoreEdgeRight) || 
			 star_y < gCCDParams.m_nIgnoreEdgeBottom || star_y>=(in.ySize-gCCDParams.m_nIgnoreEdgeUp))
			return FALSE;
			
		// now find max val on previous frames :
		inForPrev.x = star_x;
		inForPrev.y = star_y;
		inForPrev.pos = max_pos;
		double x0OnPrev=0,y0OnPrev=0,out_tresh_for_prev;
		//double noise_level = (out.m_PixelOut.maxAverageOfPrev*0.9);
		//if(noise_level<0)
		//	noise_level = fabs(in.p_curr_data_laplace[star_y][star_x]*0.9);


		int prev_x,prev_y;
		double x0,y0,outTresh;

		int pCount=0;
//		int prevCount = MIN(((in.pPipeline)->GetCount()-1),gCCDParams.m_nPrevFramesToCheckEdge);
		int prevCount = MIN(((in.pPipeline)->GetCount()-1),gCCDParams.m_nMaxOfAverageOfPrevN);

		for(register int f=1;f<=prevCount;f++){
			GetPrevPixelPos( prev_x, prev_y, star_x, star_y, f, in );
			if(prev_x>=0 && prev_y>=0 && prev_x<in.xSize && prev_y<in.ySize){		
				inForPrev.x = prev_x;
				inForPrev.y = prev_y;
				inForPrev.pos = (prev_x+prev_y*in.xSize);
				inForPrev.p_data = (in.PrevMatrixPtr[f])->get_data_buffer();
				inForPrev.p_data_fast = (in.PrevMatrixPtr[f])->get_data_buffer_fast();
				inForPrev.p_curr_data_laplace = ((CCDMatrix*)(in.PrevMatrixPtr[f]))->get_frame_laplace_fast();
				inForPrev.p_curr_data_laplace_normal = ((CCDMatrix*)(in.PrevMatrixPtr[f]))->get_frame_laplace();

				outForPrev.cluster_cnt = 0;
				FindClusterAboveTresholdOpt3( inForPrev, outForPrev, x0,y0, 
														outTresh, noise_level, FALSE );
				double max_prev = FindMaxInClusterBig( inForPrev.p_curr_data_laplace_normal,
																in.xSize, 	
																outForPrev.cluster, 
																outForPrev.cluster_cnt, max_pos );
				max_on_prev += max_prev;
			}
			pCount++;
		}
		max_on_prev = (max_on_prev/pCount);		
	}else{
		return FALSE;
	}

	return TRUE;
}						

double CCD_Analyser::CalcMaxClusterRadius( 
										LONG_T* cluster, LONG_T cluster_cnt,										
										double x0, double y0, int xSize )
{
	double max_radius=0;
	for( int i=0;i<cluster_cnt;i++ ){
		int x=( cluster[i] % xSize );
		int y=( cluster[i] / xSize ); 

		double r = sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0) );						
		if( r > max_radius ){
			max_radius = r;
		}
	}

	return max_radius;
}

int CCD_Analyser::FindClusterAboveTresholdOpt3( int x, int y, 
																int xSize, int ySize,
																ELEM_TYPE** p_data,
																CPixelAnalyseOut& out,
																LONG_T* cluster, 
																LONG_T& cluster_cnt,
																double clusterTreshold/*=-1*/,
																BIG_ELEM_TYPE** p_laplace_data,
																CCDPipeline* pPipeline,
																CPixelList* pixel_list )
{
	CPixelAnalyseIn* in = new CPixelAnalyseIn();
	in->x = x;
	in->y = y;
	in->p_data_fast = p_data;
	in->p_curr_data_laplace = p_laplace_data;	
	in->xSize = xSize;
	in->ySize = ySize;
	in->pos = y*in->xSize + x;
	CPixelList* pixel_list_local = NULL;
	if( pixel_list ){
		in->pPixelList = pixel_list;
	}else{
		pixel_list_local = new CPixelList( in->xSize*in->ySize);
		in->pPixelList = pixel_list_local;
	}
	in->pPipeline = pPipeline;
	

	double maxNoiseLevel = -1;
	double x0,y0;
	int ret = FindClusterAboveTresholdOpt3( *in, out, cluster, cluster_cnt, 
														 x0, y0, 
														 maxNoiseLevel, clusterTreshold,
														 TRUE, (p_laplace_data!=NULL), FALSE );
	out.m_PixelOut.x0 = (int)(x0);
	out.m_PixelOut.y0 = (int)(y0);		

	delete in;
	if( pixel_list_local ){
		delete pixel_list_local;
	}
	return ret;
}

/*int CCD_Analyser::GetLaplaceValue( int x, int y, int xSize,
											  ELEM_TYPE** p_data, BIG_ELEM_TYPE** p_laplace,
											  BOOL_T bFromTab )
{
	return ( (bFromTab) ? p_laplace[y][x] : Table2D<ELEM_TYPE>::CalcLaplaceSum( x, y, xSize,p_data,gCCDParams.m_eLaplaceType )  );
}*/

int CCD_Analyser::FindClusterAboveTresholdOpt3_LWP( const CPixelAnalyseIn& in,
															 CPixelAnalyseOut& out,
															 LONG_T* cluster, LONG_T& cluster_cnt,
 	 	                                        double& x0,double& y0,
															 double& maxNoiseLevel,
															 double clusterTreshold /*=-1*/,
															 BOOL_T bCalcCenter /*=TRUE*/,
															 BOOL_T bFromTab/*=TRUE*/,
															 BOOL_T bResetList/*=TRUE*/)
{
	out.neighb_count = 0;
	out.startPoints_cnt = 0;
	if(clusterTreshold<0)
		maxNoiseLevel = GetMaxNoiseLevelLaplace( &in );
	else
		maxNoiseLevel = clusterTreshold;

	cluster_cnt = 0;	
	if( bResetList ){
		out.hitpixel_count = 0;
	}

	if( in.p_curr_data_laplace[in.y][in.x] > maxNoiseLevel ){
		ADD_POINT_func( cluster, cluster_cnt,in.pos);
	}

	// number of new points added :
	int added=1;
	int r=1;

	while( added>0 ){
		added=0;
		
		int x1=(in.x-r);
		int x2=(in.x+r); 
		for( int y=(in.y-r);y<=(in.y+r);y++){
			if( y>=0 && y<in.ySize ){
				if( x1>=0 ){
					if( in.p_curr_data_laplace[y][x1] > maxNoiseLevel ){
						ADD_POINT_func( cluster, cluster_cnt, (in.xSize*y+x1) );
						added++;			
					}								
				}
				if( x2<in.xSize ){
					if( in.p_curr_data_laplace[y][x2] > maxNoiseLevel ){
						ADD_POINT_func( cluster, cluster_cnt, (in.xSize*y+x2) );
						added++;
					}								
				}
			}
		}		


		int y1=(in.y-r);
		int y2=(in.y+r); 
		for( int x=(in.x-r);x<=(in.x+r);x++){
			if( x>=0 && x<in.xSize ){
				if( y1>=0 ){
					if( in.p_curr_data_laplace[y1][x] > maxNoiseLevel ){
						ADD_POINT_func( cluster, cluster_cnt, (in.xSize*y1+x) );
						added++;
					}								
				}
				if( y2<in.ySize ){
					if( in.p_curr_data_laplace[y2][x] > maxNoiseLevel ){
						ADD_POINT_func( cluster, cluster_cnt, (in.xSize*y2+x) );
						added++;
					}								
				}
			}
		}		


		r++;
	}

	if( cluster_cnt && bCalcCenter ){
		int sum = 0;
		if( in.p_data ){
			sum = CalcCenterOfHitRealOpt( in.p_data, cluster, cluster_cnt, in.xSize, x0, y0);	
		}else{
			if( in.p_curr_data_laplace ){
				sum = CalcCenterOfHitRealOpt( in.p_curr_data_laplace, cluster, cluster_cnt, in.xSize, x0, y0);	
			}else{
				if( in.p_data_fast ){
					sum = CalcCenterOfHitRealOpt( in.p_data_fast, cluster, cluster_cnt, in.xSize, x0, y0);
				}
			}
		}
		return sum;
	}

	return 0;	
}

int CCD_Analyser::FindCluster( int x0, int y0,
										 ELEM_TYPE* data, ELEM_TYPE** data_fast,	
										 int xSize, int ySize,
										 LONG_T* cluster, LONG_T& cluster_cnt,
										 int tresh )
{
	int added=1;
	int r=1;

	ADD_POINT_func( cluster, cluster_cnt, y0*xSize+x0 );	

	while( added>0 ){
		added=0;
		
		int x1=(x0-r);
		int x2=(x0+r); 
		for( int y=(y0-r);y<=(y0+r);y++){
			if( y>=0 && y<ySize ){
				if( x1>=0 ){
					if( data_fast[y][x1] > tresh ){
						ADD_POINT_func( cluster, cluster_cnt, (xSize*y+x1) );
						added++;			
					}								
				}
				if( x2<xSize ){
					if( data_fast[y][x2] > tresh ){
						ADD_POINT_func( cluster, cluster_cnt, (xSize*y+x2) );
						added++;
					}								
				}
			}
		}		


		int y1=(y0-r);
		int y2=(y0+r); 
		for( int x=(x0-r);x<=(x0+r);x++){
			if( x>=0 && x<xSize ){
				if( y1>=0 ){
					if( data_fast[y1][x] > tresh ){
						ADD_POINT_func( cluster, cluster_cnt, (xSize*y1+x) );
						added++;
					}								
				}
				if( y2<ySize ){
					if( data_fast[y2][x] > tresh ){
						ADD_POINT_func( cluster, cluster_cnt, (xSize*y2+x) );
						added++;
					}								
				}
			}
		}		


		r++;
	}

	return cluster_cnt;	
}


int CCD_Analyser::FindClusterAboveTresholdOpt3( const CPixelAnalyseIn& in,
															 CPixelAnalyseOut& out,
															 LONG_T* cluster, LONG_T& cluster_cnt,
 	 	                                        double& x0,double& y0,
															 double& maxNoiseLevel,
															 double clusterTreshold /*=-1*/,
															 BOOL_T bCalcCenter /*=TRUE*/,
															 BOOL_T bFromTab/*=TRUE*/,
															 BOOL_T bResetList/*=TRUE*/)
{
	out.neighb_count = 0;
	out.startPoints_cnt = 0;
	if(clusterTreshold<0)
		maxNoiseLevel = GetMaxNoiseLevelLaplace( &in );
	else
		maxNoiseLevel = clusterTreshold;

	cluster_cnt = 0;	
	if( bResetList ){
		out.hitpixel_count = 0;
	}

	if(GetLaplaceValue(in.x, in.y, in.xSize, in.ySize, in.p_data_fast, in.p_curr_data_laplace, bFromTab )>maxNoiseLevel){	
		ADD_POINT_func(out.startPoints,out.startPoints_cnt,in.pos);
	}
	
   //                +
	// all pixels in +++ used for cluster finding :
   //                +
	if( in.x > 0 ){
		if( GetLaplaceValue(in.x-1, in.y, in.xSize, in.ySize, in.p_data_fast,in.p_curr_data_laplace, bFromTab )>maxNoiseLevel){
			ADD_POINT_func(out.startPoints,out.startPoints_cnt,in.pos-1);
		}
	}
	if( in.x < in.xSize-1 ){
		if(GetLaplaceValue(in.x+1, in.y, in.xSize, in.ySize, in.p_data_fast, in.p_curr_data_laplace, bFromTab )>maxNoiseLevel){
			ADD_POINT_func(out.startPoints,out.startPoints_cnt,in.pos+1);
		}
	}
	if( in.y > 0 ){	
		if( GetLaplaceValue(in.x, in.y-1, in.xSize, in.ySize, in.p_data_fast, in.p_curr_data_laplace, bFromTab )>maxNoiseLevel){
	      ADD_POINT_func(out.startPoints,out.startPoints_cnt,in.pos-in.xSize);
		}
   }
	if( in.y < in.ySize-1 ){
		if(GetLaplaceValue(in.x, in.y+1, in.xSize, in.ySize, in.p_data_fast, in.p_curr_data_laplace, bFromTab )>maxNoiseLevel){
	      ADD_POINT_func(out.startPoints,out.startPoints_cnt,in.pos+in.xSize);
		}
   }

	for(register int i=0;i<out.startPoints_cnt;i++){
      (in.pPixelList)->HitPixel( out.startPoints[i] );
		ADD_POINT_func(out.hitpixel_list,out.hitpixel_count, out.startPoints[i] );
   }


		

	while(out.startPoints_cnt){
		 register int last = out.startPoints_cnt-1;
		 register int pos = out.startPoints[last];			
		 register int x = ( pos % in.xSize );
       register int y =  ( pos / in.xSize );			
		 if(!ADD_POINT_func(cluster,cluster_cnt,pos)){
	       cluster_cnt=0;
			 return 0;
		 }

	    (in.pPixelList)->HitPixel( pos );
		 ADD_POINT_func(out.hitpixel_list,out.hitpixel_count, pos );
				 
		 out.startPoints_cnt--;
       // GetAllNeighbNoSelf(x,y,pos,in.xSize,in.ySize,out.neighbours,out.neighb_count);
		 // GetAnalNeighbNoSelf(x,y,pos,in.xSize,in.ySize,out.neighbours,out.neighb_count);

		out.neighb_count = 0;
		for(register int xx=-1;xx<=1;xx++){
			for(register int yy=-1;yy<=1;yy++){
				register int _xx = x+xx;
				register int _yy = y+yy;

				if(_xx>0 && _xx<in.xSize && _yy>0 && _yy<in.ySize){
					out.neighbours[out.neighb_count] = _yy*in.xSize + _xx;
					out.neighb_count++;	
				}
			}
		}

		for(register int i=0;i<out.neighb_count && out.neighbours[i]!=NEIGHB_LIST_END;i++){
			// if( out.neighbours[i]==pos || find_value( cluster,cluster_cnt,out.neighbours[i] )!=NOT_FOUND )
			if( out.neighbours[i]==pos || (in.pPixelList)->CheckPixel( out.neighbours[i] ) )
               continue;
			int neighb_x = ( out.neighbours[i] % in.xSize );
			int neighb_y = ( out.neighbours[i] / in.xSize );
			
			if(GetLaplaceValue(neighb_x, neighb_y, in.xSize, in.ySize, in.p_data_fast, in.p_curr_data_laplace, bFromTab )>maxNoiseLevel){
				if(!ADD_POINT_func(out.startPoints,out.startPoints_cnt,out.neighbours[i])){
					cluster_cnt=0;
					return 0;
				}			
		   }

			(in.pPixelList)->HitPixel( out.neighbours[i] );
			ADD_POINT_func(out.hitpixel_list,out.hitpixel_count, out.neighbours[i] );
		}
	}

	x0 = in.x;
	y0 = in.y;

	// cleaning - clusters are cleaned :
	// pixels are marked only to avoid same pixels in single cluster not between
	// different clusters on single frame :
	(in.pPixelList)->ClearPart( out.hitpixel_list,out.hitpixel_count );

	if( cluster_cnt && bCalcCenter ){
		int sum = 0;
		if( in.p_data ){
			sum = CalcCenterOfHitRealOpt( in.p_data, cluster, cluster_cnt, in.xSize, x0, y0);	
		}else{
			if( in.p_curr_data_laplace ){
				sum = CalcCenterOfHitRealOpt( in.p_curr_data_laplace, cluster, cluster_cnt, in.xSize, x0, y0);	
			}else{
				if( in.p_data_fast ){
					sum = CalcCenterOfHitRealOpt( in.p_data_fast, cluster, cluster_cnt, in.xSize, x0, y0);
				}
			}
		}
		return sum;
	}

	return 0;
}


int CCD_Analyser::FindClusterAboveTresholdOpt4( const CPixelAnalyseIn& in,
															 CPixelAnalyseOut& out,
															 LONG_T* cluster, LONG_T& cluster_cnt,
 	 	                                        double& x0,double& y0,
															 double& maxNoiseLevel,
															 double clusterTreshold,
															 BOOL_T bFromTab/*=TRUE*/ )
{
	out.neighb_count = 0;
	out.startPoints_cnt = 0;
	maxNoiseLevel = clusterTreshold;

	cluster_cnt = 0;	

	if(GetLaplaceValue(in.x, in.y, in.xSize, in.ySize, in.p_data_fast, in.p_curr_data_laplace, bFromTab )>maxNoiseLevel){	
		ADD_POINT_func(out.startPoints,out.startPoints_cnt,in.pos);
	}
	
   //                +
	// all pixels in +++ used for cluster finding :
   //                +
	if( in.x > 0 ){
		if( GetLaplaceValue(in.x-1, in.y, in.xSize, in.ySize, in.p_data_fast,in.p_curr_data_laplace, bFromTab )>maxNoiseLevel){
			ADD_POINT_func(out.startPoints,out.startPoints_cnt,in.pos-1);
		}
	}
	if( in.x < in.xSize-1 ){
		if(GetLaplaceValue(in.x+1, in.y, in.xSize, in.ySize, in.p_data_fast, in.p_curr_data_laplace, bFromTab )>maxNoiseLevel){
			ADD_POINT_func(out.startPoints,out.startPoints_cnt,in.pos+1);
		}
	}
	if( in.y > 0 ){	
		if( GetLaplaceValue(in.x, in.y-1, in.xSize, in.ySize, in.p_data_fast, in.p_curr_data_laplace, bFromTab )>maxNoiseLevel){
	      ADD_POINT_func(out.startPoints,out.startPoints_cnt,in.pos-in.xSize);
		}
   }
	if( in.y < in.ySize-1 ){
		if(GetLaplaceValue(in.x, in.y+1, in.xSize, in.ySize, in.p_data_fast, in.p_curr_data_laplace, bFromTab )>maxNoiseLevel){
	      ADD_POINT_func(out.startPoints,out.startPoints_cnt,in.pos+in.xSize);
		}
   }

	for(register int i=0;i<out.startPoints_cnt;i++){
      (in.pPixelList)->HitPixel( out.startPoints[i] );
   }


		

	while(out.startPoints_cnt){
		 register int last = out.startPoints_cnt-1;
		 register int pos = out.startPoints[last];			
		 register int x = ( pos % in.xSize );
       register int y =  ( pos / in.xSize );			
		 if(!ADD_POINT_func(cluster,cluster_cnt,pos)){
	       cluster_cnt=0;
			 return 0;
		 }

	    (in.pPixelList)->HitPixel( pos );
				 
		 out.startPoints_cnt--;
       // GetAllNeighbNoSelf(x,y,pos,in.xSize,in.ySize,out.neighbours,out.neighb_count);
		 // GetAnalNeighbNoSelf(x,y,pos,in.xSize,in.ySize,out.neighbours,out.neighb_count);

		out.neighb_count = 0;
		for(register int xx=-1;xx<=1;xx++){
			for(register int yy=-1;yy<=1;yy++){
				register int _xx = x+xx;
				register int _yy = y+yy;

				if(_xx>0 && _xx<in.xSize && _yy>0 && _yy<in.ySize){
					out.neighbours[out.neighb_count] = _yy*in.xSize + _xx;
					out.neighb_count++;	
				}
			}
		}

		for(register int i=0;i<out.neighb_count && out.neighbours[i]!=NEIGHB_LIST_END;i++){
			// if( out.neighbours[i]==pos || find_value( cluster,cluster_cnt,out.neighbours[i] )!=NOT_FOUND )
			if( out.neighbours[i]==pos || (in.pPixelList)->CheckPixel( out.neighbours[i] ) )
               continue;
			int neighb_x = ( out.neighbours[i] % in.xSize );
			int neighb_y = ( out.neighbours[i] / in.xSize );
			
			if(GetLaplaceValue(neighb_x, neighb_y, in.xSize, in.ySize, in.p_data_fast, in.p_curr_data_laplace, bFromTab )>maxNoiseLevel){
				if(!ADD_POINT_func(out.startPoints,out.startPoints_cnt,out.neighbours[i])){
					cluster_cnt=0;
					return 0;
				}			
		   }

			(in.pPixelList)->HitPixel( out.neighbours[i] );
		}
	}
	x0 = in.x;
	y0 = in.y;

	return 0;
}



int CCD_Analyser::CalcLapClusterOfMaxPoints( CPixelAnalyseOut& out, 
												  int x, int y, int xSize, int ySize,
                                      ELEM_TYPE** p_data_fast, ELEM_TYPE* p_data,
													int n_points, 
													int n_break_size, // internal radius of ring ( 3 added below )
													int n_apert_size, // ring width
													float& sky,
													float& plus_sum_out,
													double& x0, double& y0,
													BOOL_T bSubtractSky )
{
	int n = FindClusterOfMaxNPoints( x, y, xSize, ySize, p_data_fast,
											 	out.cluster, out.cluster_cnt, n_points );

	// [NEW] commented 20041204
	/*GetOutsidePixels( out.cluster, out.cluster_cnt,
                   out.neighb_list, out.ncnt,
                   xSize, ySize, n_break_size );
	for(register int i=0;i<out.cluster_cnt;i++){
		if( find_value( out.neighb_list, out.ncnt , out.cluster[i] )<0 ){
			ADD_POINT( out.neighb_list, out.ncnt, out.cluster[i] );
		}
	}	*/
	// [NEW] commented 20041204
	/*GetOutsidePixels( out.neighb_list, out.ncnt,
                   out.ClusterWithMore, out.ClusterWithMoreCnt,
                   xSize, ySize, n_apert_size  );*/


	int r_internal = n_break_size+3; // +2 for 3x3 square
	int r_external = r_internal+n_apert_size;
	GetOutsidePixels( x , y , out.ClusterWithMore, out.ClusterWithMoreCnt,
							xSize, ySize, r_internal , r_external );
							
	// in fact should subtract further values not just 

	int plus_sum = CalcClusterSum( out.cluster, out.cluster_cnt, p_data );

	/*
	int minus_sum = CalcClusterSum( out.ClusterWithMore, out.ClusterWithMoreCnt, p_data );
	int lap = (int)((plus_sum) - double(out.cluster_cnt)/double(out.ClusterWithMoreCnt)*minus_sum);
	*/

	// median :
	for( register int i=0;i<out.ClusterWithMoreCnt;i++){
		out.ClusterWithMore[i] = p_data[ out.ClusterWithMore[i] ];
	}
	int minus_sum = get_no_stars_median( out.ClusterWithMore, out.ClusterWithMoreCnt );
	int lap = plus_sum;

	if( bSubtractSky ){
		lap = lap - out.cluster_cnt*minus_sum;
	}

	sky = minus_sum;	
	plus_sum_out = plus_sum;

	
	if( CalcCenterOfHitRealOpt( p_data, out.cluster, out.cluster_cnt, xSize, x0, y0)==0 ){
		x0 = x;
		y0 = y;			
	}

	return lap;
}


int CCD_Analyser::CalcLapClusterOfMaxPointsAnyShape( CPixelAnalyseOut& out, 
												  int x, int y, int xSize, int ySize,
                                      ELEM_TYPE** p_data_fast, ELEM_TYPE* p_data,
													int n_points, 
													int n_break_size, // internal radius of ring ( 3 added below )
													int n_apert_size, // ring width
													float& sky,
													float& plus_sum_out,
													double& x0, double& y0,
													BOOL_T bSubtractSky )
{
	int n = FindClusterOfMaxNPointsAnyShape( x, y, xSize, ySize, p_data_fast,
											 	out.cluster, out.cluster_cnt, n_points );

	int r_internal = n_break_size+3; // +2 for 3x3 square
	int r_external = r_internal+n_apert_size;
	GetOutsidePixels( x , y , out.ClusterWithMore, out.ClusterWithMoreCnt,
							xSize, ySize, r_internal , r_external );
							
	// in fact should subtract further values not just 

	int plus_sum = CalcClusterSum( out.cluster, out.cluster_cnt, p_data );

	// median :
	for( register int i=0;i<out.ClusterWithMoreCnt;i++){
		out.ClusterWithMore[i] = p_data[ out.ClusterWithMore[i] ];
	}
	int minus_sum = get_no_stars_median( out.ClusterWithMore, out.ClusterWithMoreCnt );
	int lap = plus_sum;

	if( bSubtractSky ){
		lap = lap - out.cluster_cnt*minus_sum;
	}

	sky = minus_sum;	
	plus_sum_out = plus_sum;

	
	if( CalcCenterOfHitRealOpt( p_data, out.cluster, out.cluster_cnt, xSize, x0, y0)==0 ){
		x0 = x;
		y0 = y;			
	}

	return lap;
}


int CCD_Analyser::FindClusterOfMaxNPoints( int x0, int y0, int xSize, int ySize,
														 ELEM_TYPE** p_data_fast,
														 LONG_T* cluster,LONG_T& cluster_cnt,
													    int n_points )
{
	cluster_cnt = 1;
	int pos = y0*xSize	+ x0;
	cluster[0] = pos;
	
	// cluster 5x5 - just for testing :
	if( n_points >= 25 && n_points<36 ){
		cluster_cnt = 0;
		for(register int yy=(y0-2);yy<=(y0+2);yy++){
			int p0 = yy*xSize;
			for(register int xx=(x0-2);xx<=(x0+2);xx++){
				if( xx>=0 && yy>=0 && xx<xSize && yy<ySize ){
					cluster[cluster_cnt] = (p0+xx);
					cluster_cnt++; 
				}
			}
		}
		if( n_points==25 )
			return cluster_cnt;
	}
	
	if( n_points >= 36 && n_points<49 ){
		cluster_cnt = 0;
		for(register int yy=(y0-3);yy<=(y0+2);yy++){
			int p0 = yy*xSize;
			for(register int xx=(x0-3);xx<=(x0+2);xx++){
				if( xx>=0 && yy>=0 && xx<xSize && yy<ySize ){
					cluster[cluster_cnt] = (p0+xx);
					cluster_cnt++; 
				}
			}
		}
		if( n_points==36 )
			return cluster_cnt;
	}

	if( n_points >= 49 && n_points<64 ){
		cluster_cnt = 0;
		for(register int yy=(y0-3);yy<=(y0+3);yy++){
			int p0 = yy*xSize;
			for(register int xx=(x0-3);xx<=(x0+3);xx++){
				if( xx>=0 && yy>=0 && xx<xSize && yy<ySize ){
					cluster[cluster_cnt] = (p0+xx);
					cluster_cnt++; 
				}
			}
		}
		if( n_points==49 )
			return cluster_cnt;
	}
	
	
	// to have 3x3 always UNCOMMENT :
	cluster_cnt = 0;
	for(register int yy=(y0-1);yy<=(y0+1);yy++){
		int p0 = yy*xSize;
		for(register int xx=(x0-1);xx<=(x0+1);xx++){
			if( xx>=0 && yy>=0 && xx<xSize && yy<ySize ){
				cluster[cluster_cnt] = (p0+xx);
				cluster_cnt++; 
			}
		}
	}
	

	while( cluster_cnt<n_points ){
		int max_val=-1000;
		int p_max = -1;
		for(register int i=0;i<cluster_cnt;i++){
			int x = (cluster[i] % xSize );
			int y = (cluster[i] / xSize );

			for(register int yy=(y-1);yy<=(y+1);yy++){
				int p0 = yy*xSize;
				for(register int xx=(x-1);xx<=(x+1);xx++){
					int p = p0+xx;
					if( xx>=0 && yy>=0 && xx<xSize && yy<ySize ){
						if( p_data_fast[yy][xx]>max_val && find_value( cluster, cluster_cnt, p )<0 ){
							p_max = p;
							max_val = p_data_fast[yy][xx];
						}
					}
				}
			}
		}
		if( p_max>=0 ){
			cluster[cluster_cnt] = p_max;
			cluster_cnt++;
		}
	}
	return cluster_cnt;
}

int CCD_Analyser::FindClusterOfMaxNPointsAnyShape( int x0, int y0, int xSize, int ySize,
														 ELEM_TYPE** p_data_fast,
														 LONG_T* cluster,LONG_T& cluster_cnt,
													    int n_points )
{
	cluster_cnt = 1;
	int pos = y0*xSize	+ x0;
	cluster[0] = pos;
	
	while( cluster_cnt<n_points ){
		int max_val=-1000;
		int p_max = -1;
		for(register int i=0;i<cluster_cnt;i++){
			int x = (cluster[i] % xSize );
			int y = (cluster[i] / xSize );

			for(register int yy=(y-1);yy<=(y+1);yy++){
				int p0 = yy*xSize;
				for(register int xx=(x-1);xx<=(x+1);xx++){
					int p = p0+xx;
					if( xx>=0 && yy>=0 && xx<xSize && yy<ySize ){
						if( p_data_fast[yy][xx]>max_val && find_value( cluster, cluster_cnt, p )<0 ){
							p_max = p;
							max_val = p_data_fast[yy][xx];
						}
					}
				}
			}
		}
		if( p_max>=0 ){
			cluster[cluster_cnt] = p_max;
			cluster_cnt++;
		}
	}
	return cluster_cnt;
}



BOOL_T CCD_Analyser::FindAverageClusterOnPrevRawData( CPixelAnalyseIn& in, CPixelAnalyseOut& out,
   	                  	                    LONG_T* cluster,LONG_T& cluster_cnt, int maxNoiseLevel,
															  int curr_x, int curr_y )
{	
	int prev_x,prev_y;
	double x0,y0;

	int x_sav = in.x;
	int y_sav = in.y;	
	int pos_sav = in.pos;
	ELEM_TYPE** p_data_fast_sav = in.p_data_fast;
	ELEM_TYPE*  p_data_sav = in.p_data;


	BIG_ELEM_TYPE** p_out_data_fast = (out.m_pCluster)->get_data_buffer_fast();
	BIG_ELEM_TYPE* p_out_data = (out.m_pCluster)->get_data_buffer();
	(out.m_pCluster)->SetData(0);

	int center_x = ( (out.m_pCluster)->GetXSize() / 2);
	int center_y = ( (out.m_pCluster)->GetYSize() / 2);

	int prevCount = MIN(((in.pPipeline)->GetCount()-1),gCCDParams.m_nPrevFramesToCheckEdge);


	int pCount=0;
	for(register int f=1;f<=prevCount;f++){
		GetPrevPixelPos( prev_x, prev_y, curr_x, curr_y, f, in );
		if(prev_x>=0 && prev_y>=0 && prev_x<in.xSize && prev_y<in.ySize){		
			in.x = prev_x;
			in.y = prev_y;
			in.pos = (in.x+in.y*in.xSize);
			in.p_data = (in.PrevMatrixPtr[f])->get_data_buffer();
			in.p_data_fast = (in.PrevMatrixPtr[f])->get_data_buffer_fast();

			out.ClusterWithMoreCnt = 0;
			FindClusterAboveTresholdRawDataOpt3( in, out, x0,y0, out.ClusterWithMore, out.ClusterWithMoreCnt,
															 maxNoiseLevel, FALSE );
			for(register int i=0;i<out.ClusterWithMoreCnt;i++){
				int xx = (out.ClusterWithMore[i] % in.xSize);
				int yy = (out.ClusterWithMore[i] / in.xSize);

				int xxx = (xx-prev_x)+center_x;
				int yyy = (yy-prev_y)+center_y;	
				if( xxx>=0 && yyy>=0 && xxx<(out.m_pCluster)->GetXSize() && yyy<(out.m_pCluster)->GetYSize() ){
					p_out_data_fast[yyy][xxx] += in.p_data[out.ClusterWithMore[i]];
				}
			}
			pCount++;
		}else{
			break;
		}
	}

	for(register int i=0;i<(out.m_pCluster)->GetSize();i++){		
		p_out_data[i] = (p_out_data[i] / pCount);
	}

	in.p_data = p_data_sav;
	in.p_data_fast = p_data_fast_sav;
	in.x = x_sav;
	in.y = y_sav;
	in.pos = pos_sav;

	return TRUE;
}

BOOL_T CCD_Analyser::FindAverageClusterOnPrevLaplace( CPixelAnalyseIn& in, CPixelAnalyseOut& out,
   	                  	                    LONG_T* cluster,LONG_T& cluster_cnt, int maxNoiseLevel,
															  int curr_x, int curr_y, int prevCount )
{	
	int prev_x,prev_y;
	double x0,y0;

	int x_sav = in.x;
	int y_sav = in.y;	
	int pos_sav = in.pos;
	BIG_ELEM_TYPE** p_lap_fast_sav = in.p_curr_data_laplace;
	BIG_ELEM_TYPE* p_lap_sav = in.p_curr_data_laplace_normal;


	BIG_ELEM_TYPE** p_out_data_fast = (out.m_pCluster)->get_data_buffer_fast();
	BIG_ELEM_TYPE* p_out_data = (out.m_pCluster)->get_data_buffer();
	(out.m_pCluster)->SetData(0);

	int center_x = ( (out.m_pCluster)->GetXSize() / 2);
	int center_y = ( (out.m_pCluster)->GetYSize() / 2);


	int pCount=0;
	for(register int f=1;f<=prevCount;f++){
		GetPrevPixelPos( prev_x, prev_y, curr_x, curr_y, f, in );
		if(prev_x>=0 && prev_y>=0 && prev_x<in.xSize && prev_y<in.ySize){		
			in.x = prev_x;
			in.y = prev_y;
			in.pos = (in.x+in.y*in.xSize);
			in.p_curr_data_laplace = (in.PrevMatrixPtr[f])->get_frame_laplace_fast();							
			in.p_curr_data_laplace_normal = (in.PrevMatrixPtr[f])->get_frame_laplace();

			out.ClusterWithMoreCnt = 0;
			double max_noise_level;
			FindClusterAboveTresholdOpt3( in, out, out.ClusterWithMore, out.ClusterWithMoreCnt, 
													x0,y0, max_noise_level, maxNoiseLevel );

			for(register int i=0;i<out.ClusterWithMoreCnt;i++){
				int xx = (out.ClusterWithMore[i] % in.xSize);
				int yy = (out.ClusterWithMore[i] / in.xSize);

				int xxx = (xx-prev_x)+center_x;
				int yyy = (yy-prev_y)+center_y;	
				if( xxx>=0 && yyy>=0 && xxx<(out.m_pCluster)->GetXSize() && yyy<(out.m_pCluster)->GetYSize() ){
					p_out_data_fast[yyy][xxx] += in.p_data[out.ClusterWithMore[i]];
				}
			}
			pCount++;
		}else{
			break;
		}
	}

	for(register int i=0;i<(out.m_pCluster)->GetSize();i++){		
		p_out_data[i] = (p_out_data[i] / pCount);
	}

	in.x = x_sav;
	in.y = y_sav;
	in.pos = pos_sav;
	in.p_curr_data_laplace = p_lap_fast_sav;
	in.p_curr_data_laplace_normal = p_lap_sav;

	return TRUE;
}


int CCD_Analyser::FindClusterAboveTresholdRawDataOpt3( CPixelAnalyseIn& in,
															 CPixelAnalyseOut& out,
 	 	                                        double& x0,double& y0,
														    LONG_T* cluster,LONG_T& cluster_cnt,
															 double maxNoiseLevel,
															 BOOL_T bCalcCenter /*=TRUE*/,
															 BOOL_T bResetList/*=TRUE*/)
{
	out.neighb_count = 0;
	out.startPoints_cnt = 0;

	cluster_cnt = 0;	
	if(bResetList){
		out.hitpixel_count = 0;
	}

	if(in.p_data_fast[in.y][in.x]>maxNoiseLevel){	
		ADD_POINT_func(out.startPoints,out.startPoints_cnt,in.pos);
	}
	
   //                +
	// all pixels in +++ used for cluster finding :
   //                +
	if( in.x > 0 ){
		if(in.p_data_fast[in.y][in.x-1]>maxNoiseLevel){
			ADD_POINT_func(out.startPoints,out.startPoints_cnt,in.pos-1);
		}
	}
	if( in.x < in.xSize-1 ){
		if(in.p_data_fast[in.y][in.x+1]>maxNoiseLevel){
			ADD_POINT_func(out.startPoints,out.startPoints_cnt,in.pos+1);
		}
	}
	if( in.y > 0 ){	
		if(in.p_data_fast[ in.y-1 ][in.x]>maxNoiseLevel){
	      ADD_POINT_func(out.startPoints,out.startPoints_cnt,in.pos-in.xSize);
		}
   }
	if( in.y < in.ySize-1 ){
		if(in.p_data_fast[ in.y+1][in.x]>maxNoiseLevel){
	      ADD_POINT_func(out.startPoints,out.startPoints_cnt,in.pos+in.xSize);
		}
   }

	for(register int i=0;i<out.startPoints_cnt;i++){
      (in.pPixelList)->HitPixel( out.startPoints[i] );
		ADD_POINT_func(out.hitpixel_list,out.hitpixel_count, out.startPoints[i] );
   }


		

	while(out.startPoints_cnt){
		 register int last = out.startPoints_cnt-1;
		 register int pos = out.startPoints[last];			
		 register int x = ( pos % in.xSize );
       register int y =  ( pos / in.xSize );			
		 if(!ADD_POINT_func(cluster,cluster_cnt,pos)){
	       cluster_cnt=0;
			 return 0;
		 }

	    (in.pPixelList)->HitPixel( pos );
		 ADD_POINT_func(out.hitpixel_list,out.hitpixel_count, pos );
				 
		 out.startPoints_cnt--;
       // GetAllNeighbNoSelf(x,y,pos,in.xSize,in.ySize,out.neighbours,out.neighb_count);
		 // GetAnalNeighbNoSelf(x,y,pos,in.xSize,in.ySize,out.neighbours,out.neighb_count);

		out.neighb_count = 0;
		for(register int xx=-1;xx<=1;xx++){
			for(register int yy=-1;yy<=1;yy++){
				register int _xx = x+xx;
				register int _yy = y+yy;

				if(_xx>0 && _xx<in.xSize && _yy>0 && _yy<in.ySize){
					out.neighbours[out.neighb_count] = _yy*in.xSize + _xx;
					out.neighb_count++;	
				}
			}
		}

		for(register int i=0;i<out.neighb_count && out.neighbours[i]!=NEIGHB_LIST_END;i++){
			// if( out.neighbours[i]==pos || find_value( cluster,cluster_cnt,out.neighbours[i] )!=NOT_FOUND )
			if( out.neighbours[i]==pos || (in.pPixelList)->CheckPixel( out.neighbours[i] ) )
               continue;
			int neighb_x = ( out.neighbours[i] % in.xSize );
			int neighb_y = ( out.neighbours[i] / in.xSize );
			
			if(in.p_data_fast[neighb_y][neighb_x]>maxNoiseLevel){
				if(!ADD_POINT_func(out.startPoints,out.startPoints_cnt,out.neighbours[i])){
					cluster_cnt=0;
					return 0;
				}			
		   }

			(in.pPixelList)->HitPixel( out.neighbours[i] );
			ADD_POINT_func(out.hitpixel_list,out.hitpixel_count, out.neighbours[i] );
		}
	}

	x0 = in.x;
	y0 = in.y;

	int max_val=-10000,max_pos=-1;
	for(int i=0;i<cluster_cnt;i++){
		if( in.p_data[cluster[i]]>max_val ){
			max_val = in.p_data[cluster[i]];
			max_pos = cluster[i];
		}
	}
	if( max_pos>=0 ){
		x0 = (max_pos % in.xSize );
		y0 = (max_pos / in.xSize );
	}


	// cleaning - clusters are cleaned :
	// pixels are marked only to avoid same pixels in single cluster not between
	// different clusters on single frame :
	(in.pPixelList)->ClearPart( out.hitpixel_list,out.hitpixel_count );

	if( cluster_cnt && bCalcCenter ){
		int sum = CalcCenterOfHitRealOpt( in.p_data, cluster, cluster_cnt, in.xSize, x0, y0);	
		return sum;
	}

	return 0;
}

// skips pixels already marked as used :
inline int ADD_POINT_func_unique( LONG_T* table, LONG_T& idx, int pos, 
												 CPixelList* usedList,
												 int max_size=MAX_CLUSTER_SIZE )
{
	if( usedList->CheckPixel( pos ) ){
		return 2;
	}
	if(idx>=max_size){
	   return 0;
	}
	table[idx] = pos;
	idx++;
	return 1;   
}


int CCD_Analyser::FindClusterAboveTresholdRawDataOpt3Unique( CPixelAnalyseIn& in,
															 CPixelAnalyseOut& out,
 	 	                                        double& x0,double& y0,
														    LONG_T* cluster,LONG_T& cluster_cnt,
															 double maxNoiseLevel,
															 CPixelList* usedList,
															 BOOL_T bCalcCenter /*=TRUE*/,
															 BOOL_T bResetList/*=TRUE*/)
{
	out.neighb_count = 0;
	out.startPoints_cnt = 0;

	cluster_cnt = 0;	
	if(bResetList){
		out.hitpixel_count = 0;
	}

	if(in.p_data_fast[in.y][in.x]>maxNoiseLevel){	
		ADD_POINT_func_unique(out.startPoints,out.startPoints_cnt,in.pos,usedList);
	}
	
   //                +
	// all pixels in +++ used for cluster finding :
   //                +
	if( in.x > 0 ){
		if(in.p_data_fast[in.y][in.x-1]>maxNoiseLevel){
			ADD_POINT_func_unique(out.startPoints,out.startPoints_cnt,in.pos-1,usedList);
		}
	}
	if( in.x < in.xSize-1 ){
		if(in.p_data_fast[in.y][in.x+1]>maxNoiseLevel){
			ADD_POINT_func_unique(out.startPoints,out.startPoints_cnt,in.pos+1,usedList);
		}
	}
	if( in.y > 0 ){	
		if(in.p_data_fast[ in.y-1 ][in.x]>maxNoiseLevel){
			ADD_POINT_func_unique(out.startPoints,out.startPoints_cnt,in.pos-in.xSize,usedList);
		}
   }
	if( in.y < in.ySize-1 ){
		if(in.p_data_fast[ in.y+1][in.x]>maxNoiseLevel){
			ADD_POINT_func_unique(out.startPoints,out.startPoints_cnt,in.pos+in.xSize,usedList);
		}
   }

	for(register int i=0;i<out.startPoints_cnt;i++){
      (in.pPixelList)->HitPixel( out.startPoints[i] );
		ADD_POINT_func_unique(out.hitpixel_list,out.hitpixel_count, out.startPoints[i], usedList );
   }


		

	while(out.startPoints_cnt){
		 register int last = out.startPoints_cnt-1;
		 register int pos = out.startPoints[last];			
		 register int x = ( pos % in.xSize );
       register int y =  ( pos / in.xSize );			
		 int tmp_ret = ADD_POINT_func_unique(cluster,cluster_cnt,pos,usedList);
		 if(!tmp_ret){
	       cluster_cnt=0;
			 return 0;
		 }
		 if( tmp_ret==2 )
			continue;

	    (in.pPixelList)->HitPixel( pos );
		 ADD_POINT_func_unique(out.hitpixel_list,out.hitpixel_count, pos, usedList );
				 
		 out.startPoints_cnt--;
       // GetAllNeighbNoSelf(x,y,pos,in.xSize,in.ySize,out.neighbours,out.neighb_count);
		 // GetAnalNeighbNoSelf(x,y,pos,in.xSize,in.ySize,out.neighbours,out.neighb_count);

		out.neighb_count = 0;
		for(register int xx=-1;xx<=1;xx++){
			for(register int yy=-1;yy<=1;yy++){
				register int _xx = x+xx;
				register int _yy = y+yy;

				if(_xx>0 && _xx<in.xSize && _yy>0 && _yy<in.ySize){
					out.neighbours[out.neighb_count] = _yy*in.xSize + _xx;
					out.neighb_count++;	
				}
			}
		}

		for(register int i=0;i<out.neighb_count && out.neighbours[i]!=NEIGHB_LIST_END;i++){
			// if( out.neighbours[i]==pos || find_value( cluster,cluster_cnt,out.neighbours[i] )!=NOT_FOUND )
			if( out.neighbours[i]==pos || (in.pPixelList)->CheckPixel( out.neighbours[i] ) )
               continue;
			int neighb_x = ( out.neighbours[i] % in.xSize );
			int neighb_y = ( out.neighbours[i] / in.xSize );
			
			if(in.p_data_fast[neighb_y][neighb_x]>maxNoiseLevel){
				int tmp_ret = ADD_POINT_func_unique(out.startPoints,out.startPoints_cnt,out.neighbours[i],usedList);
				if( !tmp_ret ){
					cluster_cnt=0;
					return 0;
				}			
				if( tmp_ret == 2 )
					continue;
		   }

			(in.pPixelList)->HitPixel( out.neighbours[i] );
			ADD_POINT_func_unique(out.hitpixel_list,out.hitpixel_count, out.neighbours[i], usedList );
		}
	}

	x0 = in.x;
	y0 = in.y;

	int max_val=-10000,max_pos=-1;
	for(int i=0;i<cluster_cnt;i++){
		if( in.p_data[cluster[i]]>max_val ){
			max_val = in.p_data[cluster[i]];
			max_pos = cluster[i];
		}
	}
	if( max_pos>=0 ){
		x0 = (max_pos % in.xSize );
		y0 = (max_pos / in.xSize );
	}


	// cleaning - clusters are cleaned :
	// pixels are marked only to avoid same pixels in single cluster not between
	// different clusters on single frame :
	(in.pPixelList)->ClearPart( out.hitpixel_list,out.hitpixel_count );

	if( cluster_cnt && bCalcCenter ){
		int sum = CalcCenterOfHitRealOpt( in.p_data, cluster, cluster_cnt, in.xSize, x0, y0);	
		return sum;
	}

	return 0;
}



double CCD_Analyser::FindClusterAboveTresholdOnPrevOpt3( CPixelAnalyseIn& in,
 	 	                                        LONG_T& x0, LONG_T& y0,
     			                                  LONG_T* cluster,LONG_T& cluster_cnt,
															 double& maxNoiseLevel,double clusterTreshold /*=-1*/,
															 BOOL_T bCalcCenter /*=TRUE*/)
{
	LONG_T Neighbours[100];
	LONG_T startPoints[MAX_CLUSTER_SIZE];
	LONG_T* neighb_list = Neighbours;
	LONG_T neighb_cnt = 0;
	LONG_T all_cnt = 0;
	LONG_T startPoints_cnt = 0;
	LONG_T ncnt;

	if(clusterTreshold<0)
		maxNoiseLevel = GetMaxNoiseLevelLaplace( &in );
	else
		maxNoiseLevel = clusterTreshold;

	cluster_cnt = 0;	

	double max_val = -10000.00;
	int max_pos = -1;

	double val = CalcAverageOfPrevN( in.x, in.y, in, gCCDParams.m_nMaxOfAverageOfPrevN, TRUE );
	double val1 = 0;
	if(in.x>0)
		val1 = CalcAverageOfPrevN( in.x-1, in.y, in, gCCDParams.m_nMaxOfAverageOfPrevN, TRUE );
	double val2 = 0;
	if(in.x < in.xSize-1)
		val2 = CalcAverageOfPrevN( in.x+1, in.y, in, gCCDParams.m_nMaxOfAverageOfPrevN, TRUE );
	double val3 = 0;
	if( in.y > 0 )
		val3 = CalcAverageOfPrevN( in.x, in.y-1, in, gCCDParams.m_nMaxOfAverageOfPrevN, TRUE );
	double val4 = 0;
	if( in.y < in.ySize-1 )
		val4 = CalcAverageOfPrevN( in.x, in.y+1, in, gCCDParams.m_nMaxOfAverageOfPrevN, TRUE );

	if(val>maxNoiseLevel){	
		ADD_POINT_func(startPoints,startPoints_cnt,in.pos);
		if(val>max_val){
			max_val = val;
			max_pos = in.pos;
		}
	}
	
   //                +
	// all pixels in +++ used for cluster finding :
   //                +
	if(val1>maxNoiseLevel){
		ADD_POINT_func(startPoints,startPoints_cnt,in.pos-1);
		if(val1>max_val){
         max_val = val1;
			max_pos = in.pos-1;
		}
	}
	if(val2>maxNoiseLevel){
		ADD_POINT_func(startPoints,startPoints_cnt,in.pos+1);
		if(val2>max_val){
         max_val = val2;
			max_pos = in.pos+1;
		}
	}
	if(val3>maxNoiseLevel){
      ADD_POINT_func(startPoints,startPoints_cnt,in.pos-in.xSize);
		if(val3>max_val){
         max_val = val3;
			max_pos = in.pos-in.xSize;
		}
   }
	if(val4>maxNoiseLevel){
      ADD_POINT_func(startPoints,startPoints_cnt,in.pos+in.xSize);
		if(val4>max_val){
         max_val = val4;
			max_pos = in.pos+in.xSize;
		}
   }

	for(register int i=0;i<startPoints_cnt;i++){
      (in.pPixelList)->HitPixel( startPoints[i] );
   }


		

	while(startPoints_cnt){
		 register int last = startPoints_cnt-1;
		 register int pos = startPoints[last];			
		 register int x = ( pos % in.xSize );
       register int y =  ( pos / in.xSize );			
		 if(!ADD_POINT_func(cluster,cluster_cnt,pos)){
	       cluster_cnt=0;
			 return 0;
		 }
	    (in.pPixelList)->HitPixel( pos );
				 
		 startPoints_cnt--;
       // GetAllNeighbNoSelf(x,y,pos,in.xSize,in.ySize,neighb_list,ncnt);
		 // GetAnalNeighbNoSelf(x,y,pos,in.xSize,in.ySize,neighb_list,ncnt);

		ncnt = 0;
		for(register int xx=-1;xx<=1;xx++){
			for(register int yy=-1;yy<=1;yy++){
				register int _xx = x+xx;
				register int _yy = y+yy;

				if(_xx>0 && _xx<in.xSize && _yy>0 && _yy<in.ySize){
					neighb_list[ncnt] = _yy*in.xSize + _xx;
					ncnt++;	
				}
			}
		}

		for(register int i=0;i<ncnt && neighb_list[i]!=NEIGHB_LIST_END;i++){
			// if( neighb_list[i]==pos || find_value( cluster,cluster_cnt,neighb_list[i] )!=NOT_FOUND )
			if( neighb_list[i]==pos || (in.pPixelList)->CheckPixel( neighb_list[i] ) )
               continue;
			int neighb_x = ( neighb_list[i] % in.xSize );
			int neighb_y = ( neighb_list[i] / in.xSize );
			
			double new_val = CalcAverageOfPrevN( neighb_x, neighb_y, in, gCCDParams.m_nMaxOfAverageOfPrevN, TRUE );
			if(new_val>maxNoiseLevel){
				if(!ADD_POINT_func(startPoints,startPoints_cnt,neighb_list[i])){
					cluster_cnt=0;
					return 0;
				}			
				if(new_val>max_val){
		         max_val = new_val;
					max_pos = neighb_list[i];
				}
		   }
			(in.pPixelList)->HitPixel( neighb_list[i] );
		}
	}
	
	if(cluster_cnt>0){
		x0 = max_pos % in.xSize;
		y0 = max_pos / in.xSize;
	}

	return max_val;
}



void CCD_Analyser::FindMaxClusterNew(CCDMatrix& ccd_matrix,CCDMatrix& ccd_max,
                                  CLongList& startPoints,
                                  CLongList& cluster,
                                  long CorePoint,CPixelList& in_clusters)
{	
	long x_n,y_n,sum,sum_max;
	CLongList Neighbours,allUsed,AnalNeighb,Outer;
	ELEM_TYPE *p_data,*p_max_data;
	long xSize = ccd_matrix.GetXSize();
	long ySize = ccd_matrix.GetYSize();

	cluster.clear();
	allUsed += startPoints;
	allUsed.Add(CorePoint);
	p_data = ccd_matrix.get_data_buffer();
	p_max_data = ccd_max.get_data_buffer();
	cluster.Add(CorePoint);
	in_clusters.HitPixel( CorePoint );
	while(startPoints.size()){
		 long pos = startPoints[0];
		 long x = ( pos % xSize );
	    long y =  ( pos / xSize );
		 startPoints.erase(startPoints.begin());
		 GetNeighbours(pos,xSize,ySize,Neighbours,FALSE);
		 GetAnalNeighb(x,y,xSize,ySize,AnalNeighb,Outer);			 

		 long sum = CalcSum( AnalNeighb, p_data );; 
	    long max_sum = p_max_data[pos];
 		 if(sum==0 && max_sum>0){
    	    // pixel skiped - probably after coregister value is 0					 
		    continue;
		 }
		 if(sum - max_sum > gCCDParams.m_SumTresholdForNewFrame){
			  in_clusters.HitPixel( pos );
			  cluster.Add(pos);
			  CLongList::iterator pPoint;
			  for(pPoint=Neighbours.begin();pPoint!=Neighbours.end();pPoint++){	
			     if(!allUsed.FindPoint(*pPoint)){
				      startPoints.Add(*pPoint);
						allUsed.Add(*pPoint);
				  }
			  }
		}
	}
}


LONG_T CCD_Analyser::GetPointsAround( LONG_T pos, LONG_T xSize, LONG_T ySize, 
                                      CLongList& around, LONG_T r0 )
{
	LONG_T x = (pos % xSize);
	LONG_T y = (pos / xSize);
	LONG_T start_x = MAX(0,x-r0);
	LONG_T start_y = MAX(0,y-r0);
	LONG_T end_x = MIN(xSize-1,x+r0);
	LONG_T end_y = MIN(ySize-1,y+r0);

	around.RemoveAll();
	for(LONG_T x0=start_x;x0<end_x;x0++){
		for(LONG_T y0=start_y;y0<end_y;y0++){
			LONG_T pos = x0 + y0*xSize;
			double r = sqrt((x0-x)*(x0-x)+(y0-y)*(y0-y));
			if(r <= r0){
				around.Add( pos );
			}
		}
	}
	return around.size();
}

void CCD_Analyser::GetClusterFromShapeMap( LONG_T x0, LONG_T y0, 
                                           LONG_T xSize, LONG_T ySize,
                   				             LONG_T* ClusterWithMore,
                              			    LONG_T& ClusterWithMoreCnt )
{
	ClusterWithMoreCnt = 0;

	for(int j=0;j<m_ConfShapeCount;j++){
		LONG_T x = x0 + m_ShapeToConfirmMap[j].x;
		LONG_T y = y0 + m_ShapeToConfirmMap[j].y;
		LONG_T pos = x + y*xSize;
		if( x>=0 && x<xSize && y>=0 && y<ySize){
			ADD_POINT( ClusterWithMore, ClusterWithMoreCnt, pos );
		}
	}
	ClusterWithMore[ClusterWithMoreCnt] = NEIGHB_LIST_END;
}

int CCD_Analyser::GetClusterEdge( LONG_T* cluster, LONG_T cluster_cnt,
											 LONG_T xSize, LONG_T ySize,
											 LONG_T* EdgePixelsList, LONG_T& EdgePixelsCount )
{
	EdgePixelsCount=0;
	LONG_T prev_y=-1;	
	LONG_T prev_pos=-1;
	LONG_T y0,lower_y;
	LONG_T x,y;
	LONG_T idx=0;


	my_qsort( cluster, cluster_cnt );

	y0 = (cluster[idx] / xSize);
	lower_y = y0;
	y = y0;
	while(y == y0 && idx<cluster_cnt){
		ADD_POINT( EdgePixelsList, EdgePixelsCount, cluster[idx]  )
		idx++;
		y = (cluster[idx] / xSize);
	}

	for(;idx<cluster_cnt;idx++){
		x = ( cluster[idx] % xSize );			
		y = ( cluster[idx] / xSize );			
		
		if(y!=prev_y){
			if(prev_y>=0){
				ADD_POINT( EdgePixelsList, EdgePixelsCount, prev_pos  )
			}
			ADD_POINT( EdgePixelsList, EdgePixelsCount, cluster[idx] )
		}
		prev_pos = cluster[idx];
		prev_y = y;
	}	

	// now add upper pixels - but only in case it is not a single horizontal line :
	y0 = y;

	if(y0 != lower_y && cluster_cnt>1){
		// only in case not single line - recent y!=lower_y 

		idx--;
		y = (cluster[idx] / xSize);
		while(y == y0){
			ADD_POINT( EdgePixelsList, EdgePixelsCount, cluster[idx] );
			idx--;
			y = (cluster[idx] / xSize);
		}
		EdgePixelsCount--;
	}
	return EdgePixelsCount;
}


void CCD_Analyser::GetOutsidePixels( int x, int y,
												 LONG_T* out_cluster, LONG_T& out_cnt,
												 LONG_T xSize, LONG_T ySize,
												 int r0, int r1 )
{
	out_cnt = 0;
	for(register int yy=(y-r1);yy<=(y+r1);yy++){
		for(register int xx=(x-r1);xx<=(x+r1);xx++){
			double dist = sqrt( (xx-x)*(xx-x) + (yy-y)*(yy-y) );
			if( dist>=r0 && dist<=r1 && yy>=0 && yy<ySize && xx>=0 && xx<xSize){
				out_cluster[ out_cnt ] = yy*xSize+xx;
				out_cnt++;
			}
		}
	}
}
										

void CCD_Analyser::GetOutsidePixels( LONG_T* cluster, LONG_T cluster_cnt,
     	                               LONG_T* out_cluster, LONG_T& out_cnt,  
	                                  LONG_T xSize, LONG_T ySize,
												 int nPixels )
{
	LONG_T EdgePixelsList[MAX_CLUSTER_SIZE];
	LONG_T EdgePixelsCount=0;

	out_cnt = 0;
	int edge_cnt = GetClusterEdge( cluster, cluster_cnt, xSize, ySize,
											 EdgePixelsList, EdgePixelsCount );
	// cluster is sorted in GetClusterEdge


	if( edge_cnt ){
	   CValuesCounter val_counter;

		for(int i=0;i<EdgePixelsCount;i++){
			LONG_T x0 = (EdgePixelsList[i] % xSize);
			LONG_T y0 = (EdgePixelsList[i] / xSize);			
			for(int jy=(y0-nPixels);jy<=(y0+nPixels);jy++){
				for(int jx=(x0-nPixels);jx<=(x0+nPixels);jx++){
					if( jx>=0 && jx<xSize && jy>=0 && jy<ySize){
						int pos = jy*xSize+jx;
						if( find_value( cluster, cluster_cnt, pos)<0 ){
							val_counter.Add( pos );
							/*if( find_value( out_cluster, out_cnt, pos )<0 ){
								ADD_POINT( out_cluster, out_cnt, pos );
							}*/
						}
					}
				}
			}
		}		

		out_cnt = val_counter.cnt;
		for( register int l=0;l<val_counter.cnt;l++){
			out_cluster[l] = ((val_counter.m_pValues)[l]).value;
		}

		my_qsort( out_cluster, out_cnt );		
	}		
}

/*int CCD_Analyser::CalcClusterSum( LONG_T* cluster, LONG_T cluster_cnt,
											 ELEM_TYPE* p_data )
{
	register int sum=0;
	for(register int i=0;i<cluster_cnt;i++){
		sum += p_data[ cluster[i] ];
	}
	return sum;
}*/

int CCD_Analyser::CalcClusterLaplace( int x, int y, LONG_T xSize, LONG_T ySize,
                           			  ELEM_TYPE* p_data, ELEM_TYPE** p_data_fast,
												  int treshold )
{
	CPixelAnalyseOut out;	
	CPixelAnalyseIn* in = new CPixelAnalyseIn();
	in->x = x;
	in->y = y;
	in->p_data_fast = p_data_fast;
	in->p_data = p_data;
	in->pos = (x+y*xSize);
	// in->pPipeline = pPipeline;
	in->xSize = xSize;
	in->ySize = ySize;	
	CPixelList pixel_list(in->xSize*in->ySize);
   in->pPixelList = &pixel_list;


	double x0,y0;
	int sum = FindClusterAboveTresholdRawDataOpt3( *in, out, x0, y0, 
																	out.cluster, out.cluster_cnt,
																	treshold );
	GetOutsidePixels( out.cluster, out.cluster_cnt,
						 out.ClusterWithMore, out.ClusterWithMoreCnt,
						 xSize, ySize, 3 );

	delete in;
	if( out.ClusterWithMoreCnt ){
		int plus_sum = CalcClusterSum( out.cluster, out.cluster_cnt, p_data );
		int minus_sum = CalcClusterSum( out.ClusterWithMore, out.ClusterWithMoreCnt, p_data );
		
		int lap = (int)((plus_sum) - double(out.cluster_cnt)/double(out.ClusterWithMoreCnt)*minus_sum);
		return lap;
	}	
	return 0;
}


int CCD_Analyser::GetRawCluster( int x, int y, LONG_T xSize, LONG_T ySize,
                           	   ELEM_TYPE* p_data, ELEM_TYPE** p_data_fast,
											LONG_T* cluster, LONG_T& cluster_cnt,
										   int treshold )
{
	CPixelAnalyseOut out;	
	CPixelAnalyseIn* in = new CPixelAnalyseIn();
	in->x = x;
	in->y = y;
	in->p_data_fast = p_data_fast;
	in->p_data = p_data;
	in->pos = (x+y*xSize);
	// in->pPipeline = pPipeline;
	in->xSize = xSize;
	in->ySize = ySize;	
	CPixelList pixel_list(in->xSize*in->ySize);
   in->pPixelList = &pixel_list;


	double x0,y0;
	int sum = FindClusterAboveTresholdRawDataOpt3( *in, out, x0, y0, 
																	cluster, cluster_cnt,
																	treshold );

	delete in;


	return cluster_cnt;
}


void CCD_Analyser::GetClusterWithPointsAroundOpt( LONG_T* cluster, LONG_T cluster_cnt,
     	                                          LONG_T* out_cluster, LONG_T& out_cnt,  
	                                          LONG_T xSize, LONG_T ySize, LONG_T nPixels )
{
	LONG_T EdgePixelsList[MAX_CLUSTER_SIZE];
	LONG_T EdgePixelsCount=0;
	LONG_T prev_y=-1;	
	LONG_T prev_pos=-1;
	LONG_T y0,lower_y;
	LONG_T x,y;
	LONG_T idx=0;


	my_qsort( cluster, cluster_cnt );

	y0 = (cluster[idx] / xSize);
	lower_y = y0;
	y = y0;
	while(y == y0 && idx<cluster_cnt){
		ADD_POINT( EdgePixelsList, EdgePixelsCount, cluster[idx]  )
		idx++;
		y = (cluster[idx] / xSize);
	}

	for(;idx<cluster_cnt;idx++){
		x = ( cluster[idx] % xSize );			
		y = ( cluster[idx] / xSize );			
		
		if(y!=prev_y){
			if(prev_y>=0){
				ADD_POINT( EdgePixelsList, EdgePixelsCount, prev_pos  )
			}
			ADD_POINT( EdgePixelsList, EdgePixelsCount, cluster[idx] )
		}
		prev_pos = cluster[idx];
		prev_y = y;
	}	

	// now add upper pixels - but only in case it is not a single horizontal line :
	y0 = y;

	if(y0 != lower_y && cluster_cnt>1){
		// only in case not single line - recent y!=lower_y 

		idx--;
		y = (cluster[idx] / xSize);
		while(y == y0){
			ADD_POINT( EdgePixelsList, EdgePixelsCount, cluster[idx] );
			idx--;
			y = (cluster[idx] / xSize);
		}
		EdgePixelsCount--;
	}


	CValuesCounter val_counter;
	// memcpy(out_cluster,cluster,sizeof(LONG_T)*cluster_cnt);
	for(int k=0;k<cluster_cnt;k++){
		val_counter.AddNew( cluster[k] );
	}

	// now we have edge points around clusters, now add "aureola" around cluster
	for(int i=0;i<EdgePixelsCount;i++){
		LONG_T x0 = (EdgePixelsList[i] % xSize);
		LONG_T y0 = (EdgePixelsList[i] / xSize);
		for(int j=0;j<m_nPixelsInAureola;j++){
			LONG_T x = x0 + m_Aureola[j].x;
			LONG_T y = y0 + m_Aureola[j].y;
			if( x>=0 && x<xSize && y>=0 && y<ySize){
				val_counter.Add( y*xSize+x );
			}
		}
	}		

	out_cnt = val_counter.cnt;
	for( int l=0;l<val_counter.cnt;l++){
		out_cluster[l] = ((val_counter.m_pValues)[l]).value;
	}

	my_qsort( out_cluster, out_cnt );
}


void CCD_Analyser::GetClusterWithPointsAround( CLongList& cluster, 
                                               CLongList& out_cluster, 
	                                            LONG_T xSize,LONG_T ySize,
                                               LONG_T nPixels )
{
	CLongList tmp_list = cluster;
	tmp_list.Sort();
	CLongList EdgePixesList;
	CLongList::iterator pP;
	LONG_T prev_y=-1;	
	LONG_T prev_pos=-1;
	LONG_T y0;
	LONG_T x,y;
	
	pP=tmp_list.begin();
	y0 = ((*pP) / xSize);
	y = y0;
	while(y == y0 && pP!=tmp_list.end()){
		EdgePixesList.Add( *pP );
		pP++;
		y = ((*pP) / xSize);
	}

	for(;pP!=tmp_list.end();pP++){
		x = ( (*pP) % xSize );			
		y = ( (*pP) / xSize );			
		
		if(y!=prev_y){
			if(prev_y>=0)
				EdgePixesList.Add( prev_pos );
			EdgePixesList.Add( *pP );								
		}
		prev_pos = (*pP);
		prev_y = y;
	}	

	// now add lower pixels :
	y0 = y;
	pP--;
	y = ((*pP) / xSize);
	while(y == y0){
		EdgePixesList.Add( *pP );
		pP--;
		y = ((*pP) / xSize);
	}
	EdgePixesList.RemoveLast();


	out_cluster = tmp_list;
	// now we have edge points around clusters, now add "aureola" around cluster
	for(pP=EdgePixesList.begin();pP!=EdgePixesList.end();pP++){
		CLongList PointsAround;
		GetPointsAround( *pP, xSize, ySize, PointsAround, nPixels);
		out_cluster.MergeUnique( PointsAround );
	}	
}

BOOL_T CCD_Analyser::CheckNextFrameCondition( CCDMatrix* pNextMatrix, int nextMatrixX, int nextMatrixY,	
															 CccdReport* pEvent, int nextFrameIndex, CPixelAnalyseIn& in  )
{
	BIG_ELEM_TYPE** p_lapace =  pNextMatrix->get_frame_laplace_fast();

	int radius = gCCDParams.m_ConfirmOnNextRadius;	
	register int low_x = MAX((nextMatrixX-radius),0);
	register int low_y = MAX((nextMatrixY-radius),0);
	register int up_x = MIN((pNextMatrix->GetXSize()),(nextMatrixX+radius));
	register int up_y = MIN((pNextMatrix->GetYSize()),(nextMatrixY+radius));

	
	BOOL_T bRet = FALSE;
	register int maxValue = -10000;
	register int max_pos_x=-1;
	register int max_pos_y=-1;
	

	for(register int x=low_x;x<=up_x;x++){
		for(register int y=low_y;y<=up_y;y++){
			if(p_lapace[y][x]>maxValue){
				maxValue = p_lapace[y][x];
				max_pos_x = x;
				max_pos_y = y;
			}
		}
	}	


	if(maxValue>((in.pPipeline)->GetPipelineCfg()).m_nNewLaplace){		
		_TRACE_PRINTF_5("Event at (%d,%d) ACCEPTED on next frame max_value=$d, at (%d,%d)\n",(int)(pEvent->m_MaxPoint).x,(int)(pEvent->m_MaxPoint).y,
							 maxValue,max_pos_x,max_pos_y);
		bRet = TRUE;
	}else{
		_TRACE_PRINTF_5("Event at (%d,%d) REJECTED on next frame max_value=$d, at (%d,%d)\n",(int)(pEvent->m_MaxPoint).x,(int)(pEvent->m_MaxPoint).y,
							 maxValue,max_pos_x,max_pos_y);
	}
	(pEvent->m_PixelAnalResults).LaplaceOnNextFrames[nextFrameIndex] = maxValue;
	return bRet;
}

int CCD_Analyser::GetNextFrameValue( CCDMatrix* pNextMatrix, int nextMatrixX, int nextMatrixY,	
												 CccdReport* pEvent, int nextFrameIndex, CPixelAnalyseIn& in  )
{
	BIG_ELEM_TYPE** p_lapace =  pNextMatrix->get_frame_laplace_fast();
	
	register int low_x = MAX((nextMatrixX-1),0);
	register int low_y = MAX((nextMatrixY-1),0);
	register int up_x = MIN((pNextMatrix->GetXSize()),(nextMatrixX+1));
	register int up_y = MIN((pNextMatrix->GetYSize()),(nextMatrixY+1));

	// NEW CHANGE on 20050528 - due to problem with SLT 
	// when setting value on next, x,y is reacalculated to full frame values 
	// and core dumped here, so must subtract m_X_On_Big !
	// if not ok, just add parameter m_bDoCalcNextValue and disable it in SLT !
	if( pNextMatrix->m_X_On_Big>0 ){
		low_x -= pNextMatrix->m_X_On_Big;
		up_x -= pNextMatrix->m_X_On_Big;
	}
	if( pNextMatrix->m_Y_On_Big>0 ){
		low_y -= pNextMatrix->m_Y_On_Big;
		up_y -= pNextMatrix->m_Y_On_Big;
	}

	
	BOOL_T bRet = FALSE;
	register int maxValue = -10000;
	register int max_pos_x=-1;
	register int max_pos_y=-1;
	int size_x = pNextMatrix->GetXSize();
	int size_y = pNextMatrix->GetYSize();	

	for(register int x=low_x;x<=up_x;x++){
		for(register int y=low_y;y<=up_y;y++){
			if( x>=0 && x<size_x && y>=0 && y<size_y){
				if(p_lapace[y][x]>maxValue){
					maxValue = p_lapace[y][x];
					max_pos_x = x;
					max_pos_y = y;
				}
			}
		}
	}		

	return maxValue;
}


/*BOOL_T CCD_Analyser::CheckNextFrameCondition( LONGLONG_T SumToVerify,
															 LONGLONG_T NextSum,
															 LONG_T cluster_cnt,
															 LONG_T nextFrameIndex )
{
	if(NextSum < (gCCDParams.m_MaxOnNextPerPixel*cluster_cnt)){
		return TRUE;
	}
	return FALSE;	
}*/

BOOL_T CCD_Analyser::CheckConfCondition( LONGLONG_T sum, LONGLONG_T prev_sum_max,
               	     		         LONGLONG_T Tresh_signal, LONGLONG_T Tresh_noise )
{
	if( gCCDParams.m_bConfCheckDifference ){
		return ( (sum-prev_sum_max)>Tresh_signal );
	}
	if( gCCDParams.m_bConfCheckTreshAndMaxPrev ){
		return ( sum>Tresh_signal && prev_sum_max<Tresh_noise );
	}
	return FALSE;
}


BOOL_T CCD_Analyser::IsLocalMax( long pos0, long xSize, ELEM_TYPE* p_data )
{
	register long value=p_data[pos0];
	register long pos_xSize = pos0-xSize;
	register long pos_p_xSize = pos0+xSize;

	if(p_data[pos0-1]>value)
		return FALSE;
	if(p_data[pos0+1]>value)
		return FALSE;

	for(register long pos=pos_xSize-1;pos<=pos_xSize+1;pos++){
		if(p_data[pos]>value)
			return FALSE;
	}

	for(register long pos=pos_p_xSize-1;pos<=pos_p_xSize+1;pos++){
		if(p_data[pos]>value)
			return FALSE;
	}
		

	return TRUE;
}

BOOL_T CCD_Analyser::CheckLocalMaxCondition( const CPixelAnalyseIn& in, 
														   CPixelAnalyseOut& out,
                                             LONG_T& max_pos )
{
	out.m_PixelOut.m_bLocalMaxRejected=FALSE;
	if(gCCDParams.m_bCheckLaplaceCondition){
		for(register int i=0;i<gCCDParams.m_LaplacePlusCount;i++){
			register long x = in.x + gCCDParams.m_LaplacePlusList[i].x;
			register long y = in.y + gCCDParams.m_LaplacePlusList[i].y;
			register long pos = x+y*in.xSize;
			if( IsLocalMax( pos, in.xSize, in.p_data ) ){
				return TRUE;		
			}
		}
	}else{
		if( IsLocalMax( in.pos , in.xSize, in.p_data ) )
			return TRUE;
	}
	out.m_PixelOut.m_bLocalMaxRejected = TRUE;
	return FALSE;	
}

BOOL_T CCD_Analyser::CheckLaplaceConditionMedianMinus( const CPixelAnalyseIn& in,
                                                       CPixelAnalyseOut& out )
{
	static LONG_T laplacePlus[50];
	static LONG_T laplaceMinus[50];
	static LONG_T laplaceMinusValues[50];
	static LONG_T laplacePlusValues[50];
	long PlusCnt;
	long MinusCnt;
	Table2D<ELEM_TYPE>::GetLaplacePlusMinusList( laplacePlus, PlusCnt, laplaceMinus, MinusCnt,
									in.x, in.y, in.pos, in.xSize, gCCDParams.m_eLaplaceType );

	for(register long i=0;i<MinusCnt;i++){
		laplaceMinusValues[i] = in.p_data[ laplaceMinus[i] ];
	}
	my_qsort( laplaceMinusValues, MinusCnt );

	for(register long i=0;i<PlusCnt;i++){
		laplacePlusValues[i] = in.p_data[ laplacePlus[i] ];
	}
	// my_qsort( laplacePlusValues, PlusCnt );
	
	long plusSum = long_sum( laplacePlusValues, PlusCnt );
	long minusSum = long_sum( laplaceMinusValues, MinusCnt-1, 1 );

	out.m_PixelOut.laplaceSum = plusSum - (long)((double(PlusCnt)/double(MinusCnt-2))*minusSum);
	
	BOOL_T bRet = CheckLaplaceCondition( in, out );
	// printf("Confirming by median minus at (%d,%d) , plus=%d, minus=%d, lap=%d , ret=%d\n",in.x,in.y,plusSum,minusSum,out.m_PixelOut.laplaceSum,bRet);
	return bRet;
}

BOOL_T CCD_Analyser::CheckLaplaceCondition( const CPixelAnalyseIn& in,
				                        		  CPixelAnalyseOut& out )
{
	BOOL_T bRet = TRUE;
   // register int check=0;

	/*if( gCCDParams.m_bCheckDifference ){
   	bRet =  ( (newLaplace-prevLaplace)>gCCDParams.m_DiffTreshold );
	   check++;
   }*/
	/*		if( gCCDParams.m_bCombinedCheck && bRet ){
				if( newLaplace>=gCCDParams.m_UseDoubleCheckAboveADU )
					bRet = ( newLaplace>((in.pPipeline)->GetPipelineCfg()).m_nNewLaplace && 
								prevLaplace<((in.pPipeline)->GetPipelineCfg())m_nMaxLaplaceOnOther );
				else
					bRet =  ( (newLaplace-prevLaplace)>gCCDParams.m_DiffTreshold );
				check++;		
			}*/

	// if(out.laplaceSum>0){
	// 	printf("new = %f > %f ?\n",(double)out.laplaceSum,(double)((in.pPipeline)->GetPipelineCfg()).m_nNewLaplace);
	// }

	if(out.m_PixelOut.laplaceSum>((in.pPipeline)->GetPipelineCfg()).m_nNewLaplace){
		int d=2;		

		// check all neignhbours - nothing there !
		register LONG_T maxLaplace = -100000,minLaplace=100000;
		register LONG_T maxLaplaceX=-1,maxLaplaceY=-1;

		/*register int low_x = (in.x-d);
		register int up_x = (in.x+d);
		register int low_y = (in.y-d);
      register int up_y = (in.y+d);
		for(register int y=low_y;y<=up_y;y++){
			for(register int x=low_x;x<=up_x;x++){
				if(maxLaplace<in.p_laplace_data_fast[y][x])
					maxLaplace = in.p_laplace_data_fast[y][x];
			}		
		}*/

		if(gCCDParams.m_nMaxOfAverageOfPrevN>0){
			if(out.m_PixelOut.maxAverageOfPrev>=((in.pPipeline)->GetPipelineCfg()).m_nMaxLaplaceOnOther){
				return FALSE;
				// bRet = FALSE;
			}
		}

		if(bRet){
			if(in.p_laplace_data_fast){
				for(register int i=0;i<gCCDParams.m_nVetoPointsCount;i++){
					register long x = in.x+(gCCDParams.m_VetoArea)[i].x;
					register long y = in.y+(gCCDParams.m_VetoArea)[i].y;
					// printf("x=%d, y=%d\n",x,y);fflush(0);
					if(maxLaplace<in.p_laplace_data_fast[y][x]){
						maxLaplace = in.p_laplace_data_fast[y][x];
						maxLaplaceX = x;
						maxLaplaceY = y;
					}
					if(minLaplace>in.p_laplace_data_fast[y][x])
						minLaplace = in.p_laplace_data_fast[y][x];
				}
				out.m_PixelOut.prevLaplaceSum = maxLaplace;
				out.m_PixelOut.prevLaplaceSumMin = minLaplace;
				out.m_PixelOut.otherSum = maxLaplace;
				out.m_PixelOut.prevLaplaceSumX = maxLaplaceX;
				out.m_PixelOut.prevLaplaceSumY = maxLaplaceY;
			}
		}

	if(out.m_PixelOut.otherSum>=((in.pPipeline)->GetPipelineCfg()).m_nMaxLaplaceOnOther || gCCDParams.m_nMinLaplaceOnOther)
			bRet = FALSE;
	}else{
		bRet = FALSE;
	}
	return (bRet);
}

// function checks all condition in AND statement 
// in case no condition is imposed - returns FALSE !
BOOL_T CCD_Analyser::CheckCondition( LONGLONG_T sum, LONGLONG_T prev_sum_max )
{
	BOOL_T bRet = TRUE;	
	register int check=0;
	if( gCCDParams.m_bCheckDifference ){
		bRet =  ( (sum-prev_sum_max)>gCCDParams.m_DiffTreshold );
		check++;
	}
	if( gCCDParams.m_bCheckTreshAndMaxPrev && bRet ){
		// printf("tresholds : %f, %f\n",(double)gCCDParams.m_SumTresholdForNewFrame,(double)gCCDParams.m_SumOnPrevLessThen );

		bRet =  ( sum>gCCDParams.m_SumTresholdForNewFrame && prev_sum_max<gCCDParams.m_SumOnPrevLessThen );
		check++;
	}
	if( gCCDParams.m_bCombinedCheck && bRet ){
		if( sum>=gCCDParams.m_UseDoubleCheckAboveADU )
			bRet = ( sum>gCCDParams.m_SumTresholdForNewFrame && prev_sum_max<gCCDParams.m_SumOnPrevLessThen );
		else
			bRet =  ( (sum-prev_sum_max)>gCCDParams.m_DiffTreshold );
		check++;		
	}
	return (check>0 && bRet);
}

mystring CCD_Analyser::GetClusterCheckDesc(LONGLONG_T newSum,LONGLONG_T MaxSum,
                                           LONGLONG_T nConfTreshold,LONGLONG_T nMaxOnPrevAllowed)
{
	mystring szRet;

	if(gCCDParams.m_bCheckDifference ){
		szRet << "(" << newSum << "-" << MaxSum << ") > " << nConfTreshold;
	}
	if( gCCDParams.m_bCheckTreshAndMaxPrev ){
		szRet << newSum << ">" << nConfTreshold 
            << " && " << MaxSum << "<" << nMaxOnPrevAllowed;
	}
	return szRet;

}

mystring CCD_Analyser::GetCondDesc( LONGLONG_T newSum,LONGLONG_T maxSum,
                                    LONGLONG_T homeoSum )
{
	mystring szRet;
	LONGLONG_T prevSum = maxSum;
	if(gCCDParams.m_bKeepHomeopaticSum){
		prevSum = homeoSum;
	}

	if(gCCDParams.m_bCheckDifference ){
		szRet << "(" << newSum << "-" << prevSum << ") > " << gCCDParams.m_DiffTreshold;
	}
	if( gCCDParams.m_bCheckTreshAndMaxPrev ){
		szRet << newSum << ">" << gCCDParams.m_SumTresholdForNewFrame 
            << " && " << prevSum << "<" << gCCDParams.m_SumOnPrevLessThen;
	}
	return szRet;
}



double CCD_Analyser::CalcSpherSigma( int x, int y, CCDMatrix& matrix, CCDPipeline* pPipeline, double sigmaN )
{
	// find cluster and fill :
   CPixelAnalyseIn* in = new CPixelAnalyseIn();
   CPixelAnalyseOut* out = new CPixelAnalyseOut();
   in->x = x;
   in->y = y;
   in->p_data = matrix.get_data_buffer();
	in->p_data_fast = matrix.get_data_buffer_fast();
   in->pos = in->y*matrix.GetXSize()+in->x;
   in->xSize = matrix.GetXSize();
   in->ySize = matrix.GetYSize();
	CPixelList pixel_list( in->xSize*in->ySize );
   in->pPixelList = &pixel_list;
   in->ccd_index = 0;
   in->pPipeline = pPipeline;	


	gCCDParams.m_ClusterIfNSigmaAboveBackgr = sigmaN;
	CCD_Analyser::FindClusterAboveTresholdOpt2( *in, out->m_PixelOut.x0_real,out->m_PixelOut.y0_real,
							                         out->cluster, out->cluster_cnt, out->m_PixelOut.max_noise_level );
	//CCD_Analyser::FindClusterAboveTresholdOpt( &in, in->p_data, in->x, in->y, in->pos,
	//														in->xSize, in->ySize, x0_clust,y0_clust,
	//														*(in->pPixelList), out->cluster, out->cluster_cnt );
	double max_redial = CCDFastCalc::FindMaxRedial( out->cluster, out->cluster_cnt,
                                                   out->m_PixelOut.x0_real, out->m_PixelOut.y0_real, in->xSize );
	double spher = CCDFastCalc::CalcSphericity( max_redial, out->cluster_cnt );


	delete in;
	delete out;
	return spher;
}

BOOL_T CCD_Analyser::IsVerifiedOK( CccdReport& event )
{
	return ( (!event.m_PixelAnalResults.m_bRejectedByNextFrames) &&
				(!event.m_PixelAnalResults.m_bLocalMaxRejected) &&
				(!event.m_PixelAnalResults.m_bNotMoreThenNAboveRejected) &&
				(!event.m_PixelAnalResults.m_bRejectedByEventShape) &&
				(!event.m_PixelAnalResults.m_bRejectedDueToTrack) &&
				(!event.m_PixelAnalResults.m_bRejectedByCoic) );
}

int CCD_Analyser::CombineCamAnalResults( CCDEventList& cam1events, CCDEventList& cam2events )
{
	int ret=0;
	if( cam1events.size() != cam2events.size() ){
		printf("ERROR , different number of events on cameras : %d != %d\n",
					cam1events.size(),cam2events.size() );
	}

	int count = MIN( cam1events.size(),cam2events.size() );
	for(int i=0;i<count;i++){
		CccdReport& evt1 = cam1events[i];
		CccdReport& evt2 = cam2events[i];
		BOOL_T bOK=TRUE;

		if( evt1.m_PixelAnalResults.m_bRejectedDueToTrack || evt2.m_PixelAnalResults.m_bRejectedDueToTrack ){
			// in case at least one was rejected due to track - reject both :
			evt1.m_PixelAnalResults.m_bRejectedDueToTrack = TRUE;
			evt2.m_PixelAnalResults.m_bRejectedDueToTrack = TRUE;
			bOK = FALSE;
		}
		if( !bOK )
			ret++;
	}
	return ret;
}

int CCD_Analyser::FillAntyCoicList( CCDEventList& events, CCDEventList& coicevents, 
												CCDEventList& antycoic )
{
	antycoic.clear();
	CCDEventList::iterator i;
	for(i=events.begin();i!=events.end();i++){
		CccdReport* pCoicEvt = coicevents.Find( i->m_DayFrameIndex, 
															 (int)(i->m_MaxPoint).x,
															 (int)(i->m_MaxPoint).y );
		if( !pCoicEvt ){
			if( !antycoic.Find( i->m_DayFrameIndex, (int)(i->m_MaxPoint).x, (int)(i->m_MaxPoint).y ) ){
				antycoic.push_back( *i );
			}
		}
	}
	return antycoic.size();
}

int CCD_Analyser::VerifySingleCameraEvents( CCDPipeline& pipeline1, CCDEventList& eventsOn1,
														 CCDEventList& verifiedOn1  )
{
	CCD_Analyser* pAnalyser1 = pipeline1.GetAnalPtr();
	CImageCreator* pGenObj1 = pipeline1.GetGenObj();
	CCcdCfg* pParamCCD1 = (pipeline1.GetCamParamSet(0));

	verifiedOn1.clear();

	CCDEventList::iterator pEventOn1;		
	for(pEventOn1=eventsOn1.begin();pEventOn1!=eventsOn1.end();pEventOn1++){
		BOOL_T bOK = TRUE;

		if(pAnalyser1->m_pCurrSatList && gCCDParams.m_bCheckIfSatelite){
			mystring szSatName,szClosest1;
			double min_dist1;
			BOOL_T bIsSat = IsSatelite( *pEventOn1, pAnalyser1->m_pCurrSatList, 
												  szSatName, min_dist1,
												  szClosest1, gCCDParams.m_nSatRejRadius );
			double min_dist1_sec = AstroAngle::rad2arcsec( min_dist1 );

			if( m_pAnalFoundEventFunc ){
				(*m_pAnalFoundEventFunc)( &min_dist1_sec, eMinSatDist, pipeline1.GetPipelineIndex(), NULL );
			}

			BOOL_T bNotVisible=FALSE;
			if(!bIsSat){
				// try to check NOT-VISIBLE satelites now :
				if( gCCDParams.m_nNotVisibleSatRejRadius>0 ){
					mystring szClosest1NotVisi;
					double min_dist1_not_visi;

					bIsSat = IsSatelite( *pEventOn1, pAnalyser1->m_pAllSatList ,
												 szSatName, min_dist1_not_visi, 
												 szClosest1NotVisi, gCCDParams.m_nNotVisibleSatRejRadius );
					if( bIsSat ){
						min_dist1_sec = AstroAngle::rad2arcsec( min_dist1_not_visi );
						bNotVisible = TRUE;
					}else{
						double min_dist1_not_visi_sec = AstroAngle::rad2arcsec( min_dist1_not_visi );
						if( min_dist1_not_visi_sec < min_dist1_sec ){
							min_dist1_sec = min_dist1_not_visi_sec;
							szClosest1 = szClosest1NotVisi;
						}
					}
				}
			}

			if(bIsSat){
				// it is satellite 
				pEventOn1->m_EventType = EVENT_TYPE_SAT;
				pEventOn1->m_szSatName = szSatName;
				pEventOn1->m_SatCoicRadiusInSec = min_dist1_sec;
				if(gCCDParams.m_bRejectSatelite){
					bOK = FALSE;
				}else{
					pEventOn1->m_szSatName = szClosest1;
					pEventOn1->m_SatCoicRadiusInSec = min_dist1_sec;
				}				
			}else{
				pEventOn1->m_szSatName = szClosest1;
            pEventOn1->m_SatCoicRadiusInSec = min_dist1_sec;
			}

			if( gCCDParams.m_bCheckIfStar && 
				 pEventOn1->m_EventType==EVENT_TYPE_FLASH ){
				mystring szStarDesc;
				BOOL_T bIsKnownStar = IsStar( *pEventOn1, NULL, szStarDesc );
				if( bIsKnownStar ){
					pEventOn1->m_EventType = EVENT_TYPE_STAR;
					pEventOn1->m_szSatName = szStarDesc;
					if( gCCDParams.m_bRejectStars ){
						bOK = FALSE;
					}
				}
			}
		}
		
		if(bOK){
			if( verifiedOn1.FindEventByMaxPoint( (int)((pEventOn1->m_MaxPoint).x),
														    (int)((pEventOn1->m_MaxPoint).y),
															 gCCDParams.m_OverlapRedial ) ){
				bOK = FALSE;
			}
		}

		if( bOK ){
      	verifiedOn1.Add( *pEventOn1 );
      }
	}

	return verifiedOn1.size();
}

BOOL_T CCD_Analyser::IsSingleFrameEvent( CccdReport& evt1, CccdReport& evt2,
													  CCDPipeline& pipeline1, CCDPipeline& pipeline2 )
{
	if( pipeline1.ReadSingleFrameEvents() > 0 || pipeline2.ReadSingleFrameEvents()){
		int camid1 = pipeline1.GetCameraID();
		int camid2 = pipeline2.GetCameraID();
		int min_frame,max_frame;
		pipeline1.GetMinMaxAver( min_frame, max_frame );

		CccdReport tmpevt1=evt1,tmpevt2=evt2;


		for(int f=min_frame;f<=max_frame;f++){
			tmpevt1.m_DayFrameIndex = f;

			CCDEventList::iterator frame_it = lower_bound( 
											 pipeline1.m_SingleFramesEvents.begin() ,
         	                      pipeline1.m_SingleFramesEvents.end() ,
            	                   tmpevt1 ,
               	                event_finder() );
			for(CCDEventList::iterator it=frame_it;
				 (it!=pipeline1.m_SingleFramesEvents.end() && (*it).m_DayFrameIndex == f);
				 it++){
				if( camid1 == (*it).m_CameraIdx && 
					fabs( evt1.m_MaxPoint.x - (*it).m_MaxPoint.x )<=1 && 
				   fabs( evt1.m_MaxPoint.y - (*it).m_MaxPoint.y )<=1 ){
					return TRUE;
				}				
			}
		}

		for(int f=min_frame;f<=max_frame;f++){
         tmpevt2.m_DayFrameIndex = f;
			CCDEventList::iterator frame_it2 = lower_bound( 
										 pipeline2.m_SingleFramesEvents.begin() ,
                               pipeline2.m_SingleFramesEvents.end() ,
                               tmpevt2 ,
                               event_finder() );
			for(CCDEventList::iterator it=frame_it2;
				 (it!=pipeline2.m_SingleFramesEvents.end() && (*it).m_DayFrameIndex == f);
				 it++){
				if( camid2 == (*it).m_CameraIdx && 
					fabs( evt2.m_MaxPoint.x - (*it).m_MaxPoint.x )<=1 && 
				   fabs( evt2.m_MaxPoint.y - (*it).m_MaxPoint.y )<=1 ){
					return TRUE;
				}				
			}
		}

	}

	return FALSE;
}

BOOL_T CCD_Analyser::IsSingleFrameTrack( CccdReport& evt1, 
													  CCDPipeline& pipeline1 )
{
	if( pipeline1.ReadSingleFrameTracks() > 0 ){
		int camid1 = pipeline1.GetCameraID();
		int min_frame,max_frame;
		pipeline1.GetMinMaxAver( min_frame, max_frame );


		

		
		for(int f=min_frame;f<=max_frame;f++){
			for(int t=0;t<pipeline1.m_SingleFramesTracks.size();t++){
				CTrackDesc& track = pipeline1.m_SingleFramesTracks[t];

				if( track.frame_index<=f && f<=track.frame_index_end ){
//					if( CheckIfEventBelongsToTrack( evt1, f, &pipeline1,
//                                           FALSE /* no velocity check */,
//														 gCCDParams.m_fVelocityError, track ) ){
					double new_a,new_b;
					BOOL_T bReFitted;
					if( CheckIfEventBelongsToTrack( evt1, track, bReFitted, new_a, new_b,
															  FALSE, gCCDParams.m_fVelocityError, FALSE ) ){
						return TRUE;
					}
				}
			}
		}
	}

	return FALSE;
}


int CCD_Analyser::VerifyCoicydence( CCDPipeline& pipeline1, CCDPipeline& pipeline2,
												CCDEventList& eventsOn1, CCDEventList& eventsOn2,
												CCDEventList& verifiedOn1, CCDEventList& verifiedOn2,
												int& nSampleEvents, int& confirmedEvents, int i, 
												BOOL_T bDoCheckStars, const char* szType )
{
	CCD_Analyser* pAnalyser1 = pipeline1.GetAnalPtr();
	CCD_Analyser* pAnalyser2 = pipeline2.GetAnalPtr();
	CImageCreator* pGenObj1 = pipeline1.GetGenObj();
	CImageCreator* pGenObj2 = pipeline2.GetGenObj();	


		CCcdCfg* pParamCCD1 = (pipeline1.GetCamParamSet(i));
		CCcdCfg* pParamCCD2 = (pipeline2.GetCamParamSet(i));		

	if( gCCDParams.m_bOnSumedFrames>0 && gCCDParams.m_bCheckIfSatelite ){
		time_t minTime = pipeline1.GetCurrent()[0].getMinTime();
		time_t maxTime = pipeline1.GetCurrent()[0].getMaxTime();


		// checking if satellite first :
		int dt = 10;
		int loop_index=0;
		for(time_t t=minTime;t<=maxTime;t+=dt){
			// when analysing sumed frames initialize list with time in range 
			pAnalyser1->InitSatList( pipeline1, t );
			pAnalyser2->InitSatList( pipeline2, t );

			CCDEventList::iterator pEventOn1;		
			for(pEventOn1=eventsOn1.begin();pEventOn1!=eventsOn1.end();pEventOn1++){
					if( pEventOn1->m_EventType == EVENT_TYPE_SAT ){
   	      		continue;
	      	   }

					CccdReport* pEventOn2 = NULL;
					double x_prim, y_prim;

					if( gCCDParams.m_bCoicByRaDec && pipeline1.IsAsasTransformOK() && pipeline2.IsAsasTransformOK() ){
						// coicydence by RA,DEC only in case both cameras have 
						// transformation correctly determined :
						pEventOn2 = eventsOn2.FindEventByRaDec( (pEventOn1->m_AstroCoord).RA, 
																		 (pEventOn1->m_AstroCoord).Dec,
																		 gCCDParams.m_nCoicRadiusInRad );


					}else{
						if( gCCDParams.m_bCoicByRaDec ){
							printf("\n\nERROR !!!!");
							printf("COIC by RA,DEC and astrometry is not defined !!!!\n");													
							printf("All events ignored , no coic events found !\n\n\n");
						}
						pEventOn2 = NULL;
					}

					if( pEventOn2 ){
						// first check if this event was not already used
						// this means they are in fact overlaping so 
						// current pEventOn1 can be skiped :
						// note - event on eventsOn2 can be used once :
						if( verifiedOn2.FindEventByMaxPoint( (int)(pEventOn2->m_MaxPoint).x, (int)(pEventOn2->m_MaxPoint).y ) ){
							// event already on list, skiping :
							pEventOn2 = NULL;
						}
					}
	
					if( pEventOn2 ){
                  if( pEventOn2->m_EventType == EVENT_TYPE_SAT ){
                     continue;
                  }

						// checking satellites database :
 						BOOL_T bOK = CheckEventInSatDB( pAnalyser1, pAnalyser2,
				  										     pipeline1, pipeline2, 
															  &(*pEventOn1), &(*pEventOn2) );
					}
			}					
		}
	}


	CCDEventList::iterator pEventOn1;		
	for(pEventOn1=eventsOn1.begin();pEventOn1!=eventsOn1.end();pEventOn1++){
			CccdReport* pEventOn2 = NULL;
			double x_prim, y_prim;

		   // coicydence by pixels - transformation from cooregister used :
			CCDProcState::TransformCCD1_to_CCD2( (int)pEventOn1->m_Point.x, (int)pEventOn1->m_Point.y,
																 x_prim, y_prim );

			if( gCCDParams.m_bCoicByRaDec && pipeline1.IsAsasTransformOK() && pipeline2.IsAsasTransformOK() ){
				// coicydence by RA,DEC only in case both cameras have 
				// transformation correctly determined :
				pEventOn2 = eventsOn2.FindEventByRaDec( (pEventOn1->m_AstroCoord).RA, 
																	 (pEventOn1->m_AstroCoord).Dec,
																	 gCCDParams.m_nCoicRadiusInRad );


				int x_tmp,y_tmp;
				(pipeline2.GetCamStat()).ad2xy( (pEventOn1->m_AstroCoord).RA, 
															(pEventOn1->m_AstroCoord).Dec, 
															x_tmp, y_tmp );
				x_prim = x_tmp;
				y_prim = y_tmp;
				/*if( !pEventOn2 ){
					if( gCCDParams.m_bLogAntyCoic ){
						pipeline1.m_AntyCoicEvents.push_back( *pEventOn1 );
					}
				}*/
			}else{
				if( gCCDParams.m_bCoicByRaDec ){
					printf("\n\nERROR !!!!");
					printf("COIC by RA,DEC and astrometry is not defined !!!!\n");													
					printf("All events ignored , no coic events found !\n\n\n");
				}

				//pEventOn2 = eventsOn2.FindEvent( (int)x_prim, (int)y_prim, 
 				//											(int)gCCDParams.m_nCoicRedial );
				pEventOn2 = NULL;
			}

		
			// make histogram of events distances :
			if( m_pAnalFoundEventFunc ){					
				HistoDistance( pipeline1, pipeline2, *pEventOn1 , eventsOn2 , 
									x_prim, y_prim, TRUE );
			}


			BOOL_T bSampleEvent=FALSE;
			if( gCCDParams.GetMC() && pGenObj1 && pGenObj1->m_pGenEvent && pGenObj2 && pGenObj2->m_pGenEvent && 
				 pEventOn2 ){
				int genX1 = pGenObj1->m_LastPutObjectX;
				int genY1 = pGenObj1->m_LastPutObjectY;
				int genX2 = pGenObj2->m_LastPutObjectX;
				int genY2 = pGenObj2->m_LastPutObjectY;
				if( fabs(genX1-pEventOn1->m_MaxPoint.x)<=10 &&
					 fabs(genY1-pEventOn1->m_MaxPoint.y)<=10 &&
					 fabs(genX2-pEventOn2->m_MaxPoint.x)<=10 &&
					 fabs(genY2-pEventOn2->m_MaxPoint.y)<=10 ){
					bSampleEvent=TRUE;
				}
			}

			if( pEventOn2 ){
				// first check if this event was not already used
				// this means they are in fact overlaping so 
				// current pEventOn1 can be skiped :
				// note - event on eventsOn2 can be used once :
				if( verifiedOn2.FindEventByMaxPoint( (int)(pEventOn2->m_MaxPoint).x, (int)(pEventOn2->m_MaxPoint).y ) ){
					// event already on list, skiping :
					pEventOn2 = NULL;
				}
			}
	
			if(pEventOn2){
				BOOL_T bOK=TRUE;
            if( bSampleEvent ){
           		nSampleEvents++;
            }
				confirmedEvents++;

				pEventOn1->SetTransformedPoint( (pEventOn2->m_MaxPoint).x,
														  (pEventOn2->m_MaxPoint).y );
				pEventOn2->SetTransformedPoint( x_prim, y_prim, TRUE );
				pEventOn2->CalcCoicRadiusInRad( pEventOn1->m_AstroCoord.RA, 
														  pEventOn1->m_AstroCoord.Dec );
		
				if( gCCDParams.m_bOnSumedFrames>0 ){
					if( gCCDParams.m_bRejectSingleEvents ){
						if( IsSingleFrameEvent( *pEventOn1 , *pEventOn2 ,
							 							 pipeline1 , pipeline2 ) ){
							bOK = FALSE;
						}
					}else{
						if( gCCDParams.m_bCheckIfSatelite ){
							bOK = bOK && ( pEventOn1->m_EventType != EVENT_TYPE_SAT );
							if( pEventOn2 && bOK ){
								bOK = ( pEventOn2->m_EventType != EVENT_TYPE_SAT );
							}
						}
					}
					if( gCCDParams.m_bRejectSingleTracks ){
						if( pAnalyser1->IsSingleFrameTrack( *pEventOn1 , pipeline1 ) ||
                      pAnalyser2->IsSingleFrameTrack( *pEventOn2 , pipeline2 ) ){
							(pEventOn1->m_PixelAnalResults).m_bRejectedDueToTrack = TRUE;
							(pEventOn2->m_PixelAnalResults).m_bRejectedDueToTrack = TRUE;
                     bOK = FALSE;
                  }
					}
				}else{
					// checking satellites database :
 					bOK = bOK && CheckEventInSatDB( pAnalyser1, pAnalyser2,
															  pipeline1, pipeline2, 
															  &(*pEventOn1), &(*pEventOn2) );
				}
								
				if( gCCDParams.m_bCheckIfStar && 
					pEventOn1->m_EventType==EVENT_TYPE_FLASH && 
					bDoCheckStars ){
					mystring szStarDesc;
					BOOL_T bIsKnownStar = IsStar( *pEventOn1, pEventOn2, szStarDesc );
					if( bIsKnownStar ){
						pEventOn1->m_EventType = EVENT_TYPE_STAR;
						pEventOn2->m_EventType = EVENT_TYPE_STAR;
						pEventOn1->m_szSatName = szStarDesc;
						pEventOn2->m_szSatName = szStarDesc;								
						if( gCCDParams.m_bRejectStars ){
							if( gCCDParams.GetMC() && pAnalyser1 && 
								(pAnalyser1->samplesStat).size() && pAnalyser2 && 
								(pAnalyser2->samplesStat).size() && 
								pEventOn1->m_bGenerated && pEventOn2->m_bGenerated ){
								if( pGenObj1 && pGenObj1->m_pGenEvent ){
									(pAnalyser1->samplesStat).back().nStarRej++;
									(pAnalyser2->samplesStat).back().nStarRej++;
								}
							}else{
								if( (pAnalyser1->backgrStat).size() && (pAnalyser2->backgrStat).size() ){
									(pAnalyser1->backgrStat).back().nStarRej++;
   		              	   (pAnalyser2->backgrStat).back().nStarRej++;
								}
							}
							bOK = FALSE;
						}
					}
				}

				if(bOK){
					verifiedOn1.Add( *pEventOn1 );	
					verifiedOn2.Add( *pEventOn2 );
				}
		  }
	}


	printf("Coic result (%s) : %d-%d -> %d-%d\n", szType,
			eventsOn1.size(),eventsOn2.size(),
			verifiedOn1.size(),verifiedOn2.size());
}

BOOL_T CCD_Analyser::CheckEventInSatDB( CCD_Analyser* pAnalyser1, CCD_Analyser* pAnalyser2,
													 CCDPipeline& pipeline1, CCDPipeline& pipeline2,
													 CccdReport* pEventOn1, CccdReport* pEventOn2 )
{	
	BOOL_T bIsSat;
	BOOL_T bOK = CheckEventInSatDB( pAnalyser1, pAnalyser2,
										 pipeline1, pipeline2,
									  	 pEventOn1, pEventOn2, bIsSat );

	return bOK;
}


BOOL_T CCD_Analyser::CheckEventInSatDB( CCD_Analyser* pAnalyser1, CCD_Analyser* pAnalyser2,
													 CCDPipeline& pipeline1, CCDPipeline& pipeline2,
													 CccdReport* pEventOn1, CccdReport* pEventOn2,
													 BOOL_T& bIsSat )
{
	BOOL_T bOK=TRUE;
	CImageCreator* pGenObj1 = pipeline1.GetGenObj();
   CImageCreator* pGenObj2 = pipeline2.GetGenObj();

	bIsSat = FALSE;	

	if(pAnalyser1->m_pCurrSatList && pAnalyser2->m_pCurrSatList && gCCDParams.m_bCheckIfSatelite){		
		mystring szSatName,szClosest1,szClosest2;
		double min_dist1,min_dist2;
		bIsSat = IsSatelite( *pEventOn1, *pEventOn2, 
											  pAnalyser1->m_pCurrSatList , pAnalyser2->m_pCurrSatList, 
											  szSatName, min_dist1, min_dist2,
											  szClosest1, szClosest2,
											  gCCDParams.m_nSatRejRadius );
		double min_dist1_sec = AstroAngle::rad2arcsec( min_dist1 );
      double min_dist2_sec = AstroAngle::rad2arcsec( min_dist2 );

		if( m_pAnalFoundEventFunc ){
			(*m_pAnalFoundEventFunc)( &min_dist1_sec, eMinSatDist, pipeline1.GetPipelineIndex(), NULL );
			(*m_pAnalFoundEventFunc)( &min_dist2_sec, eMinSatDist, pipeline2.GetPipelineIndex(), NULL );
		}

		BOOL_T bNotVisible=FALSE;
		if(!bIsSat){
			// try to check NOT-VISIBLE satelites now :
			if( gCCDParams.m_nNotVisibleSatRejRadius>0 ){
				mystring szClosest1NotVisi,szClosest2NotVisi;
				double min_dist1_not_visi, min_dist2_not_visi;

				bIsSat = IsSatelite( *pEventOn1, *pEventOn2, 
											 pAnalyser1->m_pAllSatList , pAnalyser2->m_pAllSatList, 
											 szSatName, 
											 min_dist1_not_visi, min_dist2_not_visi,
											 szClosest1NotVisi,szClosest2NotVisi,
											 gCCDParams.m_nNotVisibleSatRejRadius );
				if( bIsSat ){
					min_dist1_sec = AstroAngle::rad2arcsec( min_dist1_not_visi );
					min_dist2_sec = AstroAngle::rad2arcsec( min_dist2_not_visi );
					bNotVisible = TRUE;
				}else{
					double min_dist1_not_visi_sec = AstroAngle::rad2arcsec( min_dist1_not_visi );
					double min_dist2_not_visi_sec = AstroAngle::rad2arcsec( min_dist2_not_visi );
					if( min_dist1_not_visi_sec < min_dist1_sec || 
						 min_dist2_not_visi_sec < min_dist2_sec){
						 min_dist1_sec = min_dist1_not_visi_sec;
						 min_dist2_sec = min_dist2_not_visi_sec;
						 szClosest1 = szClosest1NotVisi;
						 szClosest2 = szClosest2NotVisi;
					}
				}
			}
		}

		if(bIsSat){
			// it is satellite 
			pEventOn1->m_EventType = EVENT_TYPE_SAT;
			pEventOn2->m_EventType = EVENT_TYPE_SAT;
			pEventOn1->m_szSatName = szSatName;
			pEventOn2->m_szSatName = szSatName;
			pEventOn1->m_SatCoicRadiusInSec = min_dist1_sec;
			pEventOn2->m_SatCoicRadiusInSec = min_dist2_sec;

			if(gCCDParams.m_bRejectSatelite){
				bOK = FALSE;
			}

			// change on 2060428 :
			// flag always not only when sat reject is enabled 
			if( gCCDParams.GetMC() && pAnalyser1 && 
				(pAnalyser1->samplesStat).size() && pAnalyser2 && 
				(pAnalyser2->samplesStat).size() &&
				pEventOn1->m_bGenerated && pEventOn2->m_bGenerated ){
				if( pGenObj1 && pGenObj1->m_pGenEvent  ){
					(pAnalyser1->samplesStat).back().nSatRej++;
					(pAnalyser2->samplesStat).back().nSatRej++;
				}
			}else{
				if( (pAnalyser1->backgrStat).size() && (pAnalyser2->backgrStat).size() ){
					 (pAnalyser1->backgrStat).back().nSatRej++;
					(pAnalyser2->backgrStat).back().nSatRej++;
				}
			}
		}else{
			if( gCCDParams.m_bOnSumedFrames>0 && gCCDParams.m_bCheckIfSatelite ){
				if( strlen( pEventOn1->m_szSatName.c_str() )==0 || 
					 min_dist1_sec < pEventOn1->m_SatCoicRadiusInSec ){
					pEventOn1->m_szSatName = szClosest1;
					pEventOn1->m_SatCoicRadiusInSec = min_dist1_sec;
				}
				if( strlen( pEventOn2->m_szSatName.c_str() )==0 ||
					 min_dist2_sec < pEventOn2->m_SatCoicRadiusInSec ){
					pEventOn2->m_szSatName = szClosest2;
					pEventOn2->m_SatCoicRadiusInSec = min_dist2_sec;
				}
			}else{
				pEventOn1->m_szSatName = szClosest1;
				pEventOn2->m_szSatName = szClosest2;	
				pEventOn1->m_SatCoicRadiusInSec = min_dist1_sec;
				pEventOn2->m_SatCoicRadiusInSec = min_dist2_sec;							
			}
		}				

	}	

	return bOK;
}



int CCD_Analyser::VerifyCoicydence( CCDPipeline& pipeline1, CCDPipeline& pipeline2 )
{
	// uncomment TOGETHER
	// pipeline1.m_allFoundEvents.back().clear();
	// pipeline2.m_allFoundEvents.back().clear();

	CCD_Analyser* pAnalyser1 = pipeline1.GetAnalPtr();
	CCD_Analyser* pAnalyser2 = pipeline2.GetAnalPtr();
	CImageCreator* pGenObj1 = pipeline1.GetGenObj();
	CImageCreator* pGenObj2 = pipeline2.GetGenObj();	

	register int CamNo = (pipeline1.GetCurrent()).GetCount();

	
	int nSampleEvents=0;
	int confirmedEvents=0;
	for(register int i=0;i<CamNo;i++){
		CCDEventList verifiedOn1,verifiedOn2;


		// uncomment TOGETHER
		//	CCDEventList& eventsOn1 = (pipeline1.GetCurrent())[i].GetFoundEvents();
		// CCDEventList& eventsOn2 = (pipeline2.GetCurrent())[i].GetFoundEvents();
		if(i<pipeline1.m_allFoundEvents.back().size() && i<pipeline2.m_allFoundEvents.back().size() ){

			CCDEventList& eventsOn1 = pipeline1.m_allFoundEvents.back()[i];
			CCDEventList& eventsOn2 = pipeline2.m_allFoundEvents.back()[i];

			CCcdCfg* pParamCCD1 = (pipeline1.GetCamParamSet(i));
			CCcdCfg* pParamCCD2 = (pipeline2.GetCamParamSet(i));		

		
			VerifyCoicydence( pipeline1, pipeline2, 
									eventsOn1, eventsOn2,
									verifiedOn1,verifiedOn2,
									nSampleEvents, confirmedEvents, i );

			if( gCCDParams.m_bLogAntyCoic ){
				FillAntyCoicList( eventsOn1, verifiedOn1, pipeline1.m_AntyCoicEvents );
				FillAntyCoicList( eventsOn2, verifiedOn2, pipeline2.m_AntyCoicEvents );
			}


			if(gCCDParams.m_eCheckIfMorePoint!=eAfterCurrentFrameOnly){
				// check against to many events after coic :
				if(gCCDParams.m_MaxNumberOfEventsOnFrameAfterCoic>0){
					if(verifiedOn1.size()>=gCCDParams.m_MaxNumberOfEventsOnFrameAfterCoic){
						printf("REJECTING AFTER COID ALL DUE TO %d > %d\n",
									verifiedOn1.size(),
									gCCDParams.m_MaxNumberOfEventsOnFrameAfterCoic);

						// all events are rejected :
						verifiedOn1.clear();
						verifiedOn2.clear();
						if( (pAnalyser1->backgrStat).size()>0 ){
							(pAnalyser1->backgrStat).back().bAllRejAfterCoic = TRUE;
						}
						if( gCCDParams.GetMC() && pAnalyser2 && (pAnalyser2->samplesStat).size() ){
							if( pGenObj1 && pGenObj1->m_pGenEvent ){
								(pAnalyser2->samplesStat).back().bAllRejAfterCoic = TRUE;
							}
						}
					}
				}
	
				if( gCCDParams.m_bSkipIfMoreThen>0 ){
					// checking for nearby events and reject them :
					CCDMatrix::RejectIfMoreThen( verifiedOn1, &pipeline1, FALSE, TRUE );
					CCDMatrix::RejectIfMoreThen( verifiedOn2, &pipeline2, FALSE, TRUE );
				}
			}

			eventsOn1 = verifiedOn1;
			eventsOn2 = verifiedOn2;
			for(register int evt=0;evt<eventsOn1.size();evt++){
				eventsOn1[evt].EvtIdx = evt;
				eventsOn2[evt].EvtIdx = evt;
			}
		}
	}


	// pipeline1.AddFoundEventsToList();
   // pipeline2.AddFoundEventsToList();
	pipeline1.GetAnalObj().UpdateEventReport( pipeline1, 0 );
   pipeline2.GetAnalObj().UpdateEventReport( pipeline2, 0 );

	if( (pAnalyser1->backgrStat).size() && (pAnalyser2->backgrStat).size() ){
		(pAnalyser1->backgrStat).back().nAfterCoic = confirmedEvents-nSampleEvents;
		(pAnalyser2->backgrStat).back().nAfterCoic = confirmedEvents-nSampleEvents;
		if( (pAnalyser1->backgrStat).back().nAfterCoic<0 )
			(pAnalyser1->backgrStat).back().nAfterCoic = 0;
		if( (pAnalyser2->backgrStat).back().nAfterCoic<0 )
			(pAnalyser2->backgrStat).back().nAfterCoic = 0;
	}
	if( (pAnalyser1->samplesStat).size() && (pAnalyser2->samplesStat).size() ){
		if( pGenObj1 && pGenObj1->m_pGenEvent ){
			(pAnalyser1->samplesStat).back().nAfterCoic = nSampleEvents;
			(pAnalyser2->samplesStat).back().nAfterCoic = nSampleEvents;
		}
	}
	

	if( gCCDParams.m_bAnalyzeSumOfPrevNFrames ){		
		if( pipeline1.m_EventsOnSumedFrame.size()>0 || pipeline2.m_EventsOnSumedFrame.size()>0 ){
			_TRACE_PRINTF_1("AnalysisOfSumFrame : %d,%d -> ",pipeline1.m_EventsOnSumedFrame.size(),pipeline2.m_EventsOnSumedFrame.size());
			if( pipeline1.m_EventsOnSumedFrame.size()>0 && pipeline2.m_EventsOnSumedFrame.size()>0 ){
				CCDEventList verifiedOnSum1,verifiedOnSum2;
				int nSample=0,nCofirmed=0;
			
				VerifyCoicydence(  pipeline1,  pipeline2, 
										 pipeline1.m_EventsOnSumedFrame, pipeline2.m_EventsOnSumedFrame,
										 verifiedOnSum1,verifiedOnSum2,
										 nSample, nCofirmed );
				((pipeline1.GetAnalPtr())->sumBackgrStat).back().nAfterCoic = verifiedOnSum1.size();
				((pipeline2.GetAnalPtr())->sumBackgrStat).back().nAfterCoic = verifiedOnSum2.size();
					
				VerifySatellitesOnSumFrame( pipeline1,  pipeline2,
													 verifiedOnSum1, verifiedOnSum2 );
	
				pipeline1.m_EventsOnSumedFrame = verifiedOnSum1;
				pipeline2.m_EventsOnSumedFrame = verifiedOnSum2;
				pipeline1.m_EventsOnSumedFrame.SetEvtIdx();
				pipeline2.m_EventsOnSumedFrame.SetEvtIdx();
			}else{
				pipeline1.m_EventsOnSumedFrame.clear();
				pipeline2.m_EventsOnSumedFrame.clear();
			}
			printf(" %d coicicing events on sum frames\n",pipeline1.m_EventsOnSumedFrame.size());
		}

		if( gCCDParams.m_bCheckForSUPERNEW ){
			if( pipeline1.m_BrightenOnSumList.size()>0 && pipeline2.m_BrightenOnSumList.size()>0 ){
				CCDEventList verifiedSN1,verifiedSN2;
				int nSample=0,nCofirmed=0;

				pipeline1.m_SNSumEventsLog.DumpNewEvents( pipeline1.m_BrightenList, "0", "SN Sum-Events:", TRUE );
				pipeline2.m_SNSumEventsLog.DumpNewEvents( pipeline2.m_BrightenList, "0", "SN Sum-Events:", TRUE );

				BOOL_T bSkipStars=FALSE;
				VerifyCoicydence(  pipeline1,  pipeline2,
										 pipeline1.m_BrightenOnSumList, pipeline2.m_BrightenOnSumList,
										 verifiedSN1,verifiedSN2,
										 nSample, nCofirmed, 0, bSkipStars, "SN_ON_SUM" );

				// Apply some specific cuts first :
				pipeline1.m_BrightenOnSumList.clear();
				pipeline2.m_BrightenOnSumList.clear();
				// adding to list :
				int count=MIN( verifiedSN1.size(), verifiedSN2.size() );
				for(int j=0;j<count;j++){
					CccdReport& evt1 = verifiedSN1[j];
					CccdReport& evt2 = verifiedSN2[j];
								
					if( gCCDParams.m_CheckEventShape<=0 ||
						 ( evt1.m_PixelAnalResults.m_Sphericity >= gCCDParams.m_CheckEventShape 
							&& evt2.m_PixelAnalResults.m_Sphericity >= gCCDParams.m_CheckEventShape )){
						pipeline1.m_BrightenOnSumList.push_back( evt1 );
						pipeline2.m_BrightenOnSumList.push_back( evt2 );
					}
				}
				pipeline1.m_BrightenOnSumList.SetEvtIdx();
         	pipeline2.m_BrightenOnSumList.SetEvtIdx();
				pipeline1.m_BrightenOnSumVerifList += pipeline1.m_BrightenOnSumList;
				pipeline2.m_BrightenOnSumVerifList += pipeline2.m_BrightenOnSumList;
			}else{
				printf("No SN_ON_SUM events identified on camers\n");
				pipeline1.m_BrightenOnSumList.clear();
				pipeline2.m_BrightenOnSumList.clear();
			}
		}

	}

	if( gCCDParams.m_nCompareToOldFreqInSec>0 ){		
		if( pipeline1.m_EventsFromCompareToOld.size()>0 || pipeline2.m_EventsFromCompareToOld.size()>0 ){
			printf("AnalysisOfOldFrame : %d,%d -> ",pipeline1.m_EventsFromCompareToOld.size(),pipeline2.m_EventsFromCompareToOld.size());
			if( pipeline1.m_EventsFromCompareToOld.size()>0 && pipeline2.m_EventsFromCompareToOld.size()>0 ){
				CCDEventList verifiedOnSum1,verifiedOnSum2;
				int nSample=0,nCofirmed=0;
			
				VerifyCoicydence(  pipeline1,  pipeline2, 
										 pipeline1.m_EventsFromCompareToOld, pipeline2.m_EventsFromCompareToOld,
										 verifiedOnSum1,verifiedOnSum2,
										 nSample, nCofirmed );
				pipeline1.m_EventsFromCompareToOld = verifiedOnSum1;
				pipeline2.m_EventsFromCompareToOld = verifiedOnSum2;									 
			}else{
				pipeline1.m_EventsFromCompareToOld.clear();
				pipeline2.m_EventsFromCompareToOld.clear();
			}
			printf(" %d coicicing events since OLD frame\n",pipeline1.m_EventsFromCompareToOld.size());
		}
	}

	if( gCCDParams.m_bCheckForSUPERNEW ){
		if( pipeline1.m_BrightenList.size()>0 && pipeline2.m_BrightenList.size()>0 ){
			CCDEventList verifiedSN1,verifiedSN2;
			int nSample=0,nCofirmed=0;

			BOOL_T bSkipStars=FALSE;
			VerifyCoicydence(  pipeline1,  pipeline2,
									 pipeline1.m_BrightenList, pipeline2.m_BrightenList,
									 verifiedSN1,verifiedSN2,
									 nSample, nCofirmed, 0, bSkipStars );

			// Apply some specific cuts first :
			pipeline1.m_BrightenList.clear();
			pipeline2.m_BrightenList.clear();
			// adding to list :
			int count=MIN( verifiedSN1.size(), verifiedSN2.size() );
			for(int j=0;j<count;j++){
				CccdReport& evt1 = verifiedSN1[j];
				CccdReport& evt2 = verifiedSN2[j];
								
				if( gCCDParams.m_CheckEventShape<=0 ||
					 ( evt1.m_PixelAnalResults.m_Sphericity >= gCCDParams.m_CheckEventShape 
						&& evt2.m_PixelAnalResults.m_Sphericity >= gCCDParams.m_CheckEventShape )){
					pipeline1.m_BrightenList.push_back( evt1 );
					pipeline2.m_BrightenList.push_back( evt2 );
				}
			}
			pipeline1.m_BrightenList.SetEvtIdx();
         pipeline2.m_BrightenList.SetEvtIdx();
			pipeline1.m_BrightenVerifList += pipeline1.m_BrightenList;
			pipeline2.m_BrightenVerifList += pipeline2.m_BrightenList;
		}else{
			pipeline1.m_BrightenList.clear();
			pipeline2.m_BrightenList.clear();
		}
	}
	
	return confirmedEvents;	
}

void CCD_Analyser::HistoDistance( CCDPipeline& ccd_pipeline1, CCDPipeline& ccd_pipeline2,
											 CccdReport& evt1, CCDEventList& eventsOn2, 
											 double x_prim, double y_prim, BOOL_T bRADEC )
{
	CCDEventList::iterator i;
	double min_dist=1000000.000;
	for(i=eventsOn2.begin();i!=eventsOn2.end();i++){
		if( bRADEC ){
			double dist_in_rad = CccdReport::CalcDist( evt1, (*i) );						
			double dist_in_sec = AstroAngle::rad2arcsec( dist_in_rad );

			if( dist_in_sec<min_dist ){
				min_dist = dist_in_sec;
			}

			/*if( m_pAnalFoundEventFunc ){
				(*m_pAnalFoundEventFunc)( &dist_in_sec, eHistoCoicRADEC, ccd_pipeline1.GetPipelineIndex() );
			}*/
		}						
	}

	if( m_pAnalFoundEventFunc ){
		(*m_pAnalFoundEventFunc)( &min_dist, eHistoCoicRADEC, ccd_pipeline1.GetPipelineIndex(), NULL );
	}

}

BOOL_T CCD_Analyser::IsStarTycho( CccdReport& evt1, CccdReport* evt2, mystring& szStarDesc )
{
	vector<CCatalogStar> starList;
	int nStars = 0;

	double search_radius = AstroAngle::rad2arcsec( gCCDParams.m_fStarRejectRadius );
	if( m_pAnalFoundEventFunc ){
		search_radius = 900;
	}
	if( gCCDParams.m_bRejectIfBigStarNearBy ){
		if( gCCDParams.m_fBigStarRejectRadiusInArcSec > search_radius ){
			search_radius = gCCDParams.m_fBigStarRejectRadiusInArcSec;
		}
	}

	nStars = gStarCatTYCHO.getStarList( evt1.m_AstroCoord.RA, evt1.m_AstroCoord.Dec, 
													AstroAngle::arcsec2rad( search_radius ),
													starList, TRUE, -1000, gCCDParams.m_fStarCatMaxMag );

	if( nStars>0 ){
		vector<CCatalogStar>::iterator i;
		double minDist=10000.000;
		CCatalogStar star;
		double dist_in_rad;
		int ret = gStarCatDefault.findClosestStar( starList, evt1.m_AstroCoord.RA,
														 evt1.m_AstroCoord.Dec, star, dist_in_rad );
		double dist_in_sec = AstroAngle::rad2arcsec( dist_in_rad ); 

		if( m_pAnalFoundEventFunc ){
			(*m_pAnalFoundEventFunc)( &dist_in_sec, eMinDistStar , 0, NULL );
			if( AstroAngle::arcsec2rad( dist_in_sec ) > gCCDParams.m_fStarRejectRadius ){
				// checking condition here in histograming mode :
				ret = FALSE;
			}
		}

		// checking close stars 
		if(ret && dist_in_rad<=gCCDParams.m_fStarRejectRadius ){
			double ra_in_rad = AstroAngle::hours2rad( star.ra );
			szStarDesc << AstroAngle::toString( ra_in_rad, ANGLE_RA_TYPE ).c_str() 
						  << "," << star.dec << ":" << star.mag << "m:" 
						  << AstroAngle::rad2arcsec(dist_in_rad);
			return TRUE;
		}
	}


	if( gCCDParams.m_bRejectIfBigStarNearBy ){
		for(int i=0;i<starList.size();i++){
			CCatalogStar& star = starList[i];

			if( star.ra!=0 || star.dec!=0 ){		
				if( star.mag <= gCCDParams.m_fBigStarMaxMagnitudo ){
					double dist_deg = AstroCCD::CalcDistInRad( AstroAngle::hours2rad( star.ra ),
																			 AstroAngle::deg2rad(star.dec),
																			 evt1.m_AstroCoord.RA,
																			 evt1.m_AstroCoord.Dec );
					if( dist_deg*3600 < gCCDParams.m_fBigStarRejectRadiusInArcSec ){
						szStarDesc << "BIG_STAR:" 
							  <<AstroAngle::toString( AstroAngle::hours2rad( star.ra ), ANGLE_RA_TYPE ).c_str() 
							  << "," << star.dec << ":" << star.mag << "m:" 
							  << dist_deg*3600;
						return TRUE;					
					}																		 
				}
			}
		}
	}
	return FALSE;
	
}

BOOL_T CCD_Analyser::IsStar( CccdReport& evt1, CccdReport* evt2, mystring& szStarDesc )
{
	if( gCCDParams.m_bCheckStarsInTycho ){
		return IsStarTycho( evt1, evt2, szStarDesc );
	}else{
		return IsStarDefault( evt1, evt2, szStarDesc );
	}
}

BOOL_T CCD_Analyser::IsStarDefault( CccdReport& evt1, CccdReport* evt2, mystring& szStarDesc )
{
	vector<CCatalogStar> starList;
	int nStars = 0;
	if( !m_pAnalFoundEventFunc ){
		// normal - fast mode !
		nStars = gStarCatDefault.getStarList( evt1.m_AstroCoord.RA, evt1.m_AstroCoord.Dec, 
												   gCCDParams.m_fStarRejectRadius, 
													starList, TRUE );
	}else{
		// histograming mode - big radius here :
		nStars = gStarCatDefault.getStarList( evt1.m_AstroCoord.RA, evt1.m_AstroCoord.Dec, 
												  AstroAngle::arcsec2rad( 900 ), 
													starList, TRUE );
	}
	if( nStars>0 ){
		vector<CCatalogStar>::iterator i;
		double minDist=10000.000;
		CCatalogStar star;
		double dist_in_rad;
		int ret = gStarCatDefault.findClosestStar( starList, evt1.m_AstroCoord.RA,
														 evt1.m_AstroCoord.Dec, star, dist_in_rad );
		if( m_pAnalFoundEventFunc ){
			double dist_in_sec = AstroAngle::rad2arcsec( dist_in_rad );
			(*m_pAnalFoundEventFunc)( &dist_in_sec, eMinDistStar , 0, NULL );
			if( AstroAngle::arcsec2rad( dist_in_sec ) > gCCDParams.m_fStarRejectRadius ){
				// checking condition here in histograming mode :
				ret = FALSE;
			}
		}
		if(ret){
			double ra_in_rad = AstroAngle::hours2rad( star.ra );
			szStarDesc << AstroAngle::toString( ra_in_rad, ANGLE_RA_TYPE ).c_str() 
						  << "," << star.dec << ":" << star.mag << "m:" 
						  << AstroAngle::rad2arcsec(dist_in_rad);
			return TRUE;
		}
	}


	if( gCCDParams.m_bRejectIfBigStarNearBy ){
		vector<CCatalogStar> bigStarList;
		int nBigStars = gStarCatDefault.getStarList( evt1.m_AstroCoord.RA, evt1.m_AstroCoord.Dec, 
													AstroAngle::arcsec2rad(gCCDParams.m_fBigStarRejectRadiusInArcSec), 
													bigStarList, TRUE );

		// adding HIPPACOS stars to list :
		nBigStars = gStarCatDefault.getStarListHIP( evt1.m_AstroCoord.RA, evt1.m_AstroCoord.Dec, 
													AstroAngle::arcsec2rad(gCCDParams.m_fBigStarRejectRadiusInArcSec), 
													bigStarList, TRUE, FALSE );


		for(int i=0;i<bigStarList.size();i++){
			CCatalogStar& star = bigStarList[i];

			if( star.ra!=0 || star.dec!=0 ){		
				if( star.mag < gCCDParams.m_fBigStarMaxMagnitudo ){
					double dist_deg = AstroCCD::CalcDistInRad( AstroAngle::hours2rad( star.ra ),
																			 AstroAngle::deg2rad(star.dec),
																			 evt1.m_AstroCoord.RA,
																			 evt1.m_AstroCoord.Dec );
					if( dist_deg*3600 < gCCDParams.m_fBigStarRejectRadiusInArcSec ){
						szStarDesc << "BIG_STAR:" 
							  <<AstroAngle::toString( AstroAngle::hours2rad( star.ra ), ANGLE_RA_TYPE ).c_str() 
							  << "," << star.dec << ":" << star.mag << "m:" 
							  << dist_deg*3600;
						return TRUE;					
					}																		 
				}
			}
		}
	}
	return FALSE;
}

void CCD_Analyser::VerifySatellitesOnSumFrame( CCDPipeline& pipeline1, CCDPipeline& pipeline2,
															  CCDEventList& verifiedOnSum1, CCDEventList& verifiedOnSum2 )
{
	time_t startTime = pipeline1.m_PrevSumOfFramesTime;
	time_t endTime   = pipeline1.m_FrameUnixTime;

	time_t t = startTime;
	int step = 10;

	mystring szStart = get_date_time_string( startTime );
	mystring szEnd   = get_date_time_string( endTime );
	printf("Checking satellites on list of frames in time range : %s - %s\n",szStart.c_str(),szEnd.c_str() );
	while( t < endTime ){
		CSatList satlist;

		int ret = CSatInfo::GetSatInfo( satlist, t, TRUE );
		mystring szDTM = get_date_time_string( t );
		printf("read %d satellites for time : %s\n",ret,szDTM.c_str());

		int count = MIN( verifiedOnSum1.size(), verifiedOnSum2.size() );
		int nIdent = 0;
		for(int i=0;i<count;i++){
			CccdReport& evt1 = verifiedOnSum1[i];
			CccdReport& evt2 = verifiedOnSum2[i];
	
			mystring szSatName,szClose1,szClose2;
			double min_dist1,min_dist2;
			BOOL_T bSat = FALSE;
			if( evt1.IsIdentified() && evt2.IsIdentified() ){
				nIdent++;
				bSat = IsSatelite( evt1, evt2, &satlist, &satlist, szSatName,
											  min_dist1,min_dist2, szClose1,szClose2, 
											  gCCDParams.m_nSatRejRadius );
				double min_dist1_sec = AstroAngle::rad2arcsec( min_dist1 );
				double min_dist2_sec = AstroAngle::rad2arcsec( min_dist2 );

				if(bSat){
					// it is satellite 
					evt1.m_EventType = EVENT_TYPE_SAT;
					evt2.m_EventType = EVENT_TYPE_SAT;
					evt1.m_szSatName = szSatName;
					evt2.m_szSatName = szSatName;
					evt1.m_SatCoicRadiusInSec = min_dist1_sec;
					evt2.m_SatCoicRadiusInSec = min_dist2_sec;
				}else{
					if( min_dist1_sec < evt1.m_SatCoicRadiusInSec ){
						evt1.m_szSatName = szClose1;
						evt2.m_szSatName = szClose2;	
						evt1.m_SatCoicRadiusInSec = min_dist1_sec;
						evt2.m_SatCoicRadiusInSec = min_dist2_sec;							
					}
				}				
			}
		}		
		if( nIdent==0 ){
			// no more events to be checked against satellites 
			break;
		}

		t += step;		
	}
}

BOOL_T CCD_Analyser::IsSatelite( CccdReport& evt1, CccdReport& evt2,
											mystring& szSatName,
											double& min_dist1, double& min_dist2,
                                 mystring& szClosestSat1, mystring& szClosestSat2,
                                 double sat_coic_radius )
{
	return FALSE;
}


BOOL_T CCD_Analyser::IsSatelite( CccdReport& evt1, CccdReport& evt2, 
											CSatList* pList1, CSatList* pList2, 
											mystring& szSatName, 
											double& min_dist1, double& min_dist2,
											mystring& szClosestSat1, mystring& szClosestSat2,
											double sat_coic_radius )
{
	CSatList::iterator i;
	min_dist1 = 10000000.00;
	min_dist2 = 10000000.00;
	szClosestSat1 = "";
	szClosestSat2 = "";
	for(i=pList1->begin();i!=pList1->end();i++){

		// double sat_ra_deg = AstroAngle::rad2deg( i->sat_ra );
		// double sat_dec_deg = AstroAngle::rad2deg( i->sat_dec );
		

		// double evt_ra_deg = AstroAngle::rad2deg( evt1.m_AstroCoord.RA );
		// double evt_dec_deg = AstroAngle::rad2deg( evt1.m_AstroCoord.Dec );

		double dist1 = AstroAngle::getDist( evt1.m_AstroCoord.RA, evt1.m_AstroCoord.Dec, 
														i->sat_ra, i->sat_dec );
		double dist2 = AstroAngle::getDist( evt2.m_AstroCoord.RA, evt2.m_AstroCoord.Dec,
														i->sat_ra, i->sat_dec );

		/*double dist1 = sqrt( (evt1.m_AstroCoord.RA-i->sat_ra)*(evt1.m_AstroCoord.RA-i->sat_ra)+
								   (evt1.m_AstroCoord.Dec-i->sat_dec)*(evt1.m_AstroCoord.Dec-i->sat_dec) );
		double dist2 = sqrt( (evt2.m_AstroCoord.RA-i->sat_ra)*(evt2.m_AstroCoord.RA-i->sat_ra)+
								   (evt2.m_AstroCoord.Dec-i->sat_dec)*(evt2.m_AstroCoord.Dec-i->sat_dec) );
		*/

		if( dist1<min_dist1 ){
			min_dist1 = dist1;
			szClosestSat1 = i->sat_name;
		}
		if( dist2<min_dist2 ){
			min_dist2 = dist2;
			szClosestSat2 = i->sat_name;
		}

		// gCCDParams.m_nSatRejRadius
		if( dist1<sat_coic_radius || dist2<sat_coic_radius){
			szSatName = i->sat_name;
		
			evt1.m_EventType = EVENT_TYPE_SAT;
			evt2.m_EventType = EVENT_TYPE_SAT;
			evt1.m_szSatName = szSatName;
			evt2.m_szSatName = szSatName;
			evt1.m_szSatName << " - " << AstroAngle::rad2arcsec( dist1 );
			evt2.m_szSatName << " - " << AstroAngle::rad2arcsec( dist2 );

			if(i->sat_geo_stat){
				evt1.m_SatType = 'G';
				evt2.m_SatType = 'G';
			}
			evt1.visibility = i->visibility;
			evt2.visibility = i->visibility;
			memcpy( evt1.m_SatInfo, &(*i), sizeof(struct satInfo) );
			memcpy( evt2.m_SatInfo, &(*i), sizeof(struct satInfo) );
		
			return TRUE;
		}
	}
	return FALSE;
}

BOOL_T CCD_Analyser::IsSatelite( CccdReport& evt1, CSatList* pList1,
											mystring& szSatName, double& min_dist1,
											mystring& szClosestSat1, double sat_coic_radius )
{
	CSatList::iterator i;
	min_dist1 = 10000000.00;
	szClosestSat1 = "";
	for(i=pList1->begin();i!=pList1->end();i++){
		double dist1 = AstroAngle::getDist( evt1.m_AstroCoord.RA, evt1.m_AstroCoord.Dec, 
														i->sat_ra, i->sat_dec );

//		if( strcmp( i->sat_name, "Comstar_1_Rk" )==0 ){
//			printf("odo");
//		}

		if( dist1<min_dist1 ){
			min_dist1 = dist1;
			szClosestSat1 = i->sat_name;
		}
		// gCCDParams.m_nSatRejRadius
		if( dist1<sat_coic_radius ){
			szSatName = i->sat_name;
		
			evt1.m_EventType = EVENT_TYPE_SAT;
			evt1.m_szSatName = szSatName;
			evt1.m_szSatName << " - " << AstroAngle::rad2arcsec( dist1 );

			if(i->sat_geo_stat){
				evt1.m_SatType = 'G';
			}
			evt1.visibility = i->visibility;
			memcpy( evt1.m_SatInfo, &(*i), sizeof(struct satInfo) );
		
			return TRUE;
		}
	}
	return FALSE;
}


int CCD_Analyser::VerifyCoicydenceOfEvents( CCDEventList& eventsOn1, CCDEventList& eventsOn2, 
											   CCDEventList& verifiedOn1, CCDEventList& verifiedOn2,
												int max_time_diff, double coic_redial )
{
		
	CCDEventList::iterator pEventOn1;		

	
	int coicEventsCount=0;
	for(pEventOn1=eventsOn1.begin();pEventOn1!=eventsOn1.end();pEventOn1++){
		CCDEventList::iterator pEventOn2;

		//if(pEventOn1->m_FrameIndex==295 && ((int)pEventOn1->m_Point.x)==712 )
		//	printf("odo

		for(pEventOn2=eventsOn2.begin();pEventOn2!=eventsOn2.end();pEventOn2++){
			if(fabs( pEventOn1->m_Time - pEventOn2->m_Time ) < max_time_diff){
				// first finding events corresponding in time :
				double x_prim, y_prim;
				CCDProcState::TransformCCD1_to_CCD2( (int)pEventOn1->m_Point.x, (int)pEventOn1->m_Point.y,
																 x_prim, y_prim );
				pEventOn1->SetTransformedPoint( x_prim, y_prim );										
				if( fabs( x_prim - (pEventOn2->m_Point).x)<=coic_redial && 
					 fabs( y_prim - (pEventOn2->m_Point).y)<=coic_redial){	 

					verifiedOn1.Add( *pEventOn1 );	
					verifiedOn2.Add( *pEventOn2 );
					coicEventsCount++;
				}
			}
		}
	}
	return coicEventsCount;	
}



int CCD_Analyser::GetLaplaceValue( int x, int y, int xSize, int ySize,
								ELEM_TYPE** p_data, BIG_ELEM_TYPE** p_laplace,
								BOOL_T bFromTab)
{
	if(x>=5 && y>=5 && x<(xSize-5) && y<(ySize-5)){
		int lap=0;
		if( bFromTab){
			lap = p_laplace[y][x];
		}else{
			lap = Table2D<ELEM_TYPE>::CalcLaplaceSum( x, y,
						xSize,p_data,gCCDParams.m_eLaplaceType );
		}
		return lap;

		// return ( (bFromTab) ? p_laplace[y][x] : Table2D<ELEM_TYPE>::CalcLaplaceSum( x, y, xSize,p_data,gCCDParams.m_eLaplaceType )  );
	}
	return 0;					
}								



BOOL_T CCD_Analyser::VerifyEvent( CccdReport& event )
{
	if( gCCDParams.m_CheckEventShape>=0){
		if( event.m_PixelAnalResults.m_Sphericity < gCCDParams.m_CheckEventShape ||
         event.m_ClusterCount > gCCDParams.m_MaxPixelsInClusterAllowed)
      {
         return FALSE;
      }
	}
	if( event.m_PixelAnalResults.m_bRejByTrackOnSingleCam ){
		return FALSE;
	}
	return TRUE;
}

BOOL_T CCD_Analyser::CheckAnalRange( int x, int y, int xSize, int ySize )
{
	if( x<gCCDParams.m_nIgnoreEdgeLeft || x>=(xSize-gCCDParams.m_nIgnoreEdgeRight) ||
       y<gCCDParams.m_nIgnoreEdgeBottom || y>=(ySize-gCCDParams.m_nIgnoreEdgeUp) ){
		return FALSE;
	}
	return TRUE;
}


BOOL_T CCD_Analyser::CheckAnalRangeAuto( int x, int y, int xSize, int ySize )
{
	int leftEdge  = gCCDParams.m_nIgnoreEdgeLeft;
	int rightEdge = gCCDParams.m_nIgnoreEdgeRight;	
	int bottomEdge = gCCDParams.m_nIgnoreEdgeBottom;
	int upEdge = gCCDParams.m_nIgnoreEdgeUp;

	if( x<=100 ){
		leftEdge = 30;
	}
	if( x>=(xSize-100) ){
		rightEdge = 30;
	}
	if( y<=100 ){
		bottomEdge = 30;
	}
	if( y>=(ySize-100) ){
		upEdge = 30;
	}

	if( x<leftEdge || x>=(xSize-rightEdge) || y<bottomEdge || y>=(ySize-upEdge) ){
		return FALSE;
	}
	return TRUE;
}


double CCD_Analyser::CalcSum( int x, int y, int radius, double sky, ELEM_TYPE** data )
{
	int sum=0,count=0;
	
	for(register int j=(y-radius);j<=(y+radius);j++){
		for(register int i=(x-radius);i<=(x+radius);i++){
			sum += data[j][i];
			count++;
		}
	}
	
	sum = (sum - sky*count);
	return sum;
}

double CCD_Analyser::CalcSumWithSide( int x, int y, int radius, double sky, ELEM_TYPE** data )
{
	int sum=0,count=0;
	
	for(register int j=(y-radius);j<=(y+radius);j++){
		for(register int i=(x-radius);i<=(x+radius);i++){
			sum += data[j][i];
			count++;
		}
	}

	if( ( 2*radius-1 )>0 ){
		int radius2 = ( radius-1 );

		for(register int j=(y-radius2);j<=(y+radius2);j++){
			sum += data[j][x-radius-1];
			count++;

			sum += data[j][x+radius+1];
			count++;
		}

		for(register int i=(x-radius2);i<=(x+radius2);i++){
			sum += data[y-radius-1][x];
         count++;

         sum += data[y+radius+1][x];
         count++;
		}
	}
	
	sum = (sum - sky*count);
	return sum;
}

BOOL_T CCD_Analyser::AnalyseNewFrameWithAverageOfPrevAlgMC(CCDPipeline& ccd_pipeline)
{	
	Initialize();
	BOOL_T bMC = gCCDParams.GetMC();
	LONG_T pipeline_size = ccd_pipeline.GetPipelineSize();
	int pipeline_index = ccd_pipeline.GetPipelineIndex();
	LONG_T ClusterWithMoreCnt=0;
	LONG_T n_found_events = 0;
	LONG_T n_accepted_on_level_1=0;
	LONG_T n_accepted_on_level_2=0;
   LONG_T star_cluster[MAX_CLUSTER_SIZE];
	LONG_T nFramesBack = gCCDParams.m_FramesBack;
	int nDayFrame = ccd_pipeline.GetDayFrameCounter();
	register LONG_T max_pos=-1;

	if(gCCDParams.m_nMaxOfAverageOfPrevN>nFramesBack)
		nFramesBack = gCCDParams.m_nMaxOfAverageOfPrevN;
	if(gCCDParams.m_nNotMoreThenNExceedsBackFramesCount>nFramesBack)
		nFramesBack = gCCDParams.m_nNotMoreThenNExceedsBackFramesCount;

	// to have also current frame
	nFramesBack++;			


	if(ccd_pipeline.GetCount()!=pipeline_size)
		return FALSE;

	BOOL_T bRet = FALSE;

	PROFILER_START
	cCCD& newFrame = ccd_pipeline.GetCurrent();
	long size = newFrame.GetCount();
	if(size<1){
		return FALSE;
	}
	CPixelAnalyseIn in;
   CPixelAnalyseOut out;

	CPixelAnalyseIn* inForPrev = new CPixelAnalyseIn();
	CPixelAnalyseOut* outForPrev = new CPixelAnalyseOut();

	in.PrevMatrixPtrCnt = nFramesBack;
	in.frame_index = ccd_pipeline.GetFrameIndex();
	in.xSize = newFrame[0].GetXSize();
   in.ySize = newFrame[0].GetYSize();	
	int yUpperLimit = (in.ySize-gCCDParams.m_nIgnoreEdgeUp);
	int xUpperLimit = (in.xSize-gCCDParams.m_nIgnoreEdgeRight);		
	m_FirstDy = ((m_NeighbPoints[0].y)*in.xSize);


	CPixelList pixel_list(in.xSize*in.ySize),allclusters(in.xSize*in.ySize);
	CPixelList star_clusters(in.xSize*in.ySize),super_new_clusters(in.xSize*in.ySize);
	in.pPixelList = &pixel_list;


	in.pPipeline = &ccd_pipeline;
	in.pipeline_size_minus_1 = pipeline_size-1;


	// calculating of tresholds - using (0,0) square of 
	// calculated bacground , in case more precise quatisation
	// is needed - tresholds are overwritten later :
	CalculateTresholds( in );
	out.m_PixelOut.treshold_for_not_more_then_n_exceeds = in.treshold_for_not_more_then_n_exceeds;

	out.ncnt = gCCDParams.m_nNeighbToSumCount;

	CManyTab2D<BIG_ELEM_TYPE>* pHomeoFrame = NULL;
	if(gCCDParams.m_bKeepHomeopaticSum)
		pHomeoFrame = ccd_pipeline.GetHomeopaticFrame();
	CManyTab2D<BIG_ELEM_TYPE>* pLaplaceFrame = NULL;
	if(gCCDParams.m_bKeepLaplaceFrame)
		pLaplaceFrame = ccd_pipeline.GetLaplaceFrame();

	// graphical cut for SUPER-NEW stars area :
	CGraphCut graph_cut;

	// loosy cut :
   // graph_cut.AddLineByPoints( 5000, 1250, 18750, 15000, eRelSmaller );
	// strict cut :
	graph_cut.AddLineByPoints( 15000, 7500, 25000, 17500, eRelSmaller );
   // graph_cut.AddHorizontalLine( 1250, eRelGreater );
	graph_cut.AddHorizontalLine( 500, eRelGreater ); // small stars too 
   graph_cut.AddHorizontalLine( 17500, eRelSmaller );


	inForPrev->pPixelList = &star_clusters;
	(*outForPrev) = out;

	int super_new=0;

	register int nNewLaplaceTreshold = ((in.pPipeline)->GetPipelineCfg()).m_nNewLaplace;
	register int nMaxLaplaceOnOther = ((in.pPipeline)->GetPipelineCfg()).m_nMaxLaplaceOnOther;

	// supernova tresholds :
	int TvForSN = ((in.pPipeline)->GetPipelineCfg()).m_fSuperNovaTvInSigma*in.SigmaLap;
	int TnForSN = ((in.pPipeline)->GetPipelineCfg()).m_fSuperNovaTnInSigma*in.SigmaLap;
	int MinPrevForSN = ((in.pPipeline)->GetPipelineCfg()).m_fSuperNovaMinPrevValue*in.SigmaLap;

	CImageCreator* pGenObj = ccd_pipeline.GetGenObj();

	clock_t total_in_check=0;
	for(register int i=0;i<size;i++){
		in.ccd_index = i;
		in.pCCDInfo = &(ccd_pipeline.GetCCDInfoTab()[in.ccd_index]);
		in.Matrix = &(newFrame[i]);
		in.p_data = (in.Matrix)->get_data_buffer();
		in.p_data_fast = (in.Matrix)->get_data_buffer_fast();
		in.p_curr_data_laplace = (in.Matrix)->get_frame_laplace_fast();
		in.p_curr_data_laplace_normal = (in.Matrix)->get_frame_laplace();			

		in.pCamCfg = (CCcdCfg*)((in.pPipeline)->GetCamCfgTab()[in.ccd_index]);		
		(in.pCCDInfo)->ReCalcAlfaMultTab( ((in.pCamCfg)->m_CCDParams).m_RotValueDAlfa, ccd_pipeline.GetPipelineSize() );
		

		if(pHomeoFrame){
			in.p_homeo_data = (*pHomeoFrame)[i].get_data_buffer();
			in.p_fast_homeo_data = (*pHomeoFrame)[i].get_data_buffer_fast();
		}
		if(pLaplaceFrame){
			in.p_laplace_data = (*pLaplaceFrame)[i].get_data_buffer();
			in.p_laplace_data_fast = (*pLaplaceFrame)[i].get_data_buffer_fast();
		}

		in.PrevMatrixPtrCnt = ccd_pipeline.GetAllMatrixPtrsChronologicalInt( i, in, TRUE, nFramesBack );
		// inForPrev->SetPrevMatrix( in );
		(*inForPrev) = in;
																					

		register int nIgnoreEdgeBottom = ccd_pipeline.GetPipelineCfg().m_nIgnoreEdgeBottom;
		register int nIgnoreEdgeUp = ccd_pipeline.GetPipelineCfg().m_nIgnoreEdgeUp;
		register int nIgnoreEdgeRight = ccd_pipeline.GetPipelineCfg().m_nIgnoreEdgeRight;
		register int nIgnoreEdgeLeft = ccd_pipeline.GetPipelineCfg().m_nIgnoreEdgeLeft;

		register int y_pos = (nIgnoreEdgeBottom-1)*in.xSize + nIgnoreEdgeLeft;
		register int pos = 0;
		register int low_y = nIgnoreEdgeBottom;
		register int low_x = nIgnoreEdgeLeft;

		register int genX = -1;
		register int genY = -1;
		if( pGenObj ){
			genX = pGenObj->m_LastPutObjectX;
			genY = pGenObj->m_LastPutObjectY;
		}

		int bContinue=1;

		// preparing statistics objects :
		CFrameIdentStat stat( nDayFrame ), sstat( nDayFrame );
		stat.nTotal = (yUpperLimit-nIgnoreEdgeBottom)*(xUpperLimit-nIgnoreEdgeLeft);
		if( gCCDParams.GetMC() && gCCDParams.m_bPutSample && pGenObj && genX>0 && genY>0 ){
			sstat.nTotal = 1;
			stat.nTotal = 0;				
		}
		

		BOOL_T bIsGenObject=FALSE;
		BOOL_T bDoBreak=FALSE;

		_TRACE_PRINTF_3("Cam%d Analysing area : (%d,%d)-(%d,%d)\n",pipeline_index,
					nIgnoreEdgeLeft,nIgnoreEdgeBottom,(xUpperLimit-1),(yUpperLimit-1));

		for(register int y=nIgnoreEdgeBottom;y<yUpperLimit && bContinue;y++){
			y_pos += in.xSize;
			pos = y_pos;

			for(register int x=nIgnoreEdgeLeft;x<xUpperLimit && bContinue;x++,pos++){


// only in MC :
#ifndef _OPTIMIZED_VERSION_
				// MONTE CARLO - comment out in real analysis :
				bIsGenObject=FALSE;
				if( bMC && gCCDParams.m_bPutSample && pGenObj && pGenObj->m_pGenEvent){
					if(x==genX && y==genY){	
						(pGenObj->m_pGenEvent)->m_PixelAnalResults.laplaceSum = in.p_curr_data_laplace[y][x];
					}
					if( abs(x-genX)<=gCCDParams.m_bGenEventRedial && abs(y-genY)<=gCCDParams.m_bGenEventRedial ){
						bIsGenObject=TRUE;
					}else{
						stat.nTotal++;
					}
				}else{
					stat.nTotal++;
				}
#endif

				if(in.p_data_fast[y][x]>=gCCDParams.m_MaxAllowedVal){
					continue;
				}


				if( gCCDParams.m_bCalcTresholdsByBackgrMap ){
					CalcTresholdsByMap( in, nNewLaplaceTreshold, nMaxLaplaceOnOther );
				}

#ifndef _OPTIMIZED_VERSION_ 
				if( !bIsGenObject ){
					stat.nMaxAllowedValue++;
				}else{
					sstat.nMaxAllowedValue=1;
				}
#endif

				if(in.p_curr_data_laplace[y][x]<=nNewLaplaceTreshold)
					continue;

#ifndef _OPTIMIZED_VERSION_
				// MONTE CARLO - comment out in real analysis :
				if( bIsGenObject && pGenObj && pGenObj->m_pGenEvent){
					if( in.p_curr_data_laplace[y][x] > (pGenObj->m_pGenEvent)->m_PixelAnalResults.laplaceSum )
					{
						(pGenObj->m_pGenEvent)->m_PixelAnalResults.laplaceSum = in.p_curr_data_laplace[y][x];
					}
				}
#endif

				
				in.x = x;
            in.y = y;
            in.pos = pos;
#ifndef _OPTIMIZED_VERSION_
				if( !bIsGenObject ){
					stat.nTnewCut++;
				}else{
					sstat.nTnewCut=1;
				}
#endif
			
	
				//  out.m_PixelOut.maxAverageOfPrev = GetMaxAverageInVetoArea( in, TRUE );
				out.m_PixelOut.maxAverageOfPrev = -10000;
				for(register int v=0;v<gCCDParams.m_nVetoPointsCount;v++){
			       register long xx = (long)(x+(gCCDParams.m_VetoArea)[v].x);
			       register long yy = (long)(y+(gCCDParams.m_VetoArea)[v].y);
	
					 register long prevAverage = CalcAverageOfPrevN( xx, yy, in,  
																		 gCCDParams.m_nMaxOfAverageOfPrevN,
																		 TRUE );
					 if( prevAverage>out.m_PixelOut.maxAverageOfPrev ){
					     out.m_PixelOut.maxAverageOfPrev = prevAverage;
					 }
				}		
	
#ifndef _OPTIMIZED_VERSION_
				// MONTE CARLO - comment out in real analysis :
				// if(x==genX && y==genY && pGenObj && pGenObj->m_pGenEvent){
				if( bIsGenObject && pGenObj && pGenObj->m_pGenEvent){
					(pGenObj->m_pGenEvent)->m_PixelAnalResults.maxAverageOfPrev = out.m_PixelOut.maxAverageOfPrev;
				}
#endif
		
				out.m_PixelOut.eventType = eFlash;
				if(out.m_PixelOut.maxAverageOfPrev>=nMaxLaplaceOnOther){
					if(!gCCDParams.m_bCheckForSUPERNEW){
						continue;
					}else{
						// checking for SN 
						if( out.m_PixelOut.maxAverageOfPrev >= MinPrevForSN ){
							out.m_PixelOut.eventType = eBrighten;
						}else{
							continue;
						}
					}
				}

				if( m_pAnalFoundEventFunc ){
					HistoVariables( x, y, pos, in, out, stat, sstat, varInfo, in.SigmaLap );
				}


				if(!ApplyCuts( x, y, pos, in, out, stat, sstat, bIsGenObject, 
									bDoBreak, allclusters ) ){
					if( bDoBreak ){
						bContinue=0;
						break;
					}
					continue;
				}

				if( gCCDParams.m_bCheckForSUPERNEW && out.m_PixelOut.eventType==eBrighten ){
					if( out.cluster_cnt ){
						// cluster must have been calculated :
						if( !CheckForSN( in, out, TnForSN, TvForSN, MinPrevForSN ) ){
							continue;
						}						
					}else{
						continue;
					}
				}

				// NEW 20050403 - new change , added due to big background 
				// in ccdsingle analysis :
				if( gCCDParams.m_bCheckPrevOfMaxPixel ){
					int max_x = (int)out.m_PixelOut.m_MaxClusterPixel.x;
					int max_y = (int)out.m_PixelOut.m_MaxClusterPixel.y;

					if( max_x!=x || max_y!=y ){
						// in case MAX_PIXEL different then current one 
						// check if not edge of STAR :

						int maxPrev = -10000;

						for(register int v=0;v<gCCDParams.m_nVetoPointsCount;v++){
					       register long xx = (long)(max_x+(gCCDParams.m_VetoArea)[v].x);
					       register long yy = (long)(max_y+(gCCDParams.m_VetoArea)[v].y);
	
							 register long prevAverage = CalcAverageOfPrevN( xx, yy, in,  
																			 gCCDParams.m_nMaxOfAverageOfPrevN,
																			 TRUE );
							 if( prevAverage>maxPrev ){
							     maxPrev = prevAverage;
							 }
						}		

						if( maxPrev >= nMaxLaplaceOnOther){
							if(!gCCDParams.m_bCheckForSUPERNEW){
								continue;
							}else{
								// checking for SN 
								if( maxPrev >= MinPrevForSN ){
									out.m_PixelOut.eventType = eBrighten;
									out.m_PixelOut.maxAverageOfPrev = maxPrev;
								}else{
									continue;
								}
							}
						}
					}
				}

				out.m_PixelOut.PixelRawValue = in.p_data_fast[y][x];
				if( out.m_PixelOut.eventType != eBrighten ){
					(in.Matrix)->AddFoundEvent( in, out );
					n_found_events++;
				}else{
					// super nova events added to other list :
					(in.Matrix)->AddFoundEvent( (in.pPipeline)->m_BrightenList, in, out, TRUE );
				}

				if( gCCDParams.m_bSkipIfMoreThen>0 && gCCDParams.m_eCheckIfMorePoint!=eAfterCoic ){
					// check limit for total number of events - reject all
					// if exceeded :
					if(gCCDParams.m_MaxNumberOfEventsOnFrame>0){
						if( (in.Matrix)->GetFoundEvents().size() > gCCDParams.m_MaxNumberOfEventsOnFrame){
							printf("REJECTING ALL DUE TO %d > %d\n",
										(in.Matrix)->GetFoundEvents().size(),
										gCCDParams.m_MaxNumberOfEventsOnFrame);
							// in case limit exceeded , skip all and do not continue 
							(in.Matrix)->GetFoundEvents().clear();
							bContinue=0;
							stat.bAllRejected=TRUE;
							sstat.bAllRejected=TRUE;
							break;
						}
					}
				}
			}
		}

		int nFrameEvents = (in.Matrix)->GetFoundEventCount();
		
		if( m_pAnalFoundEventFunc ){
			int afterTv = stat.nMinPrevCut;

			(*m_pAnalFoundEventFunc)( &afterTv, eHistoAfterTv, (in.pPipeline)->GetPipelineIndex(), NULL );
			(*m_pAnalFoundEventFunc)( &nFrameEvents, eHistoFrameEvents, (in.pPipeline)->GetPipelineIndex(), NULL );
		}

		int nGenEvents = 0;
		if( gCCDParams.GetMC() && pGenObj && genX>0 && genY>0 ){
			nGenEvents = (in.Matrix)->CountFoundEvents( genX, genY, 10 );
			if( nGenEvents>1 )
				nGenEvents = 1;
		}

		BOOL_T bToAll=FALSE;
		if( !bToAll ){
			if( gCCDParams.m_bFitLineToSingleFrame ){
				if( FitLineToSingleFrameEvents( (in.Matrix)->GetFoundEvents() ) ){
					printf("line fit SUCCEDED !!!\n");
				}
			}
		}

		// adding additional information to events :
		// before VerifySingleCamTracks - to have event time :
		AddAditionalEventInfo( in, &stat );

		if( gCCDParams.m_bCheckTracksOnSingleCam ){
			VerifySingleCamTracks( in, *(in.Matrix), (in.Matrix)->GetFoundEvents(), 
										  (in.pPipeline)->m_SingleCamEvents );
		}

		stat.nIfMoreThen = nFrameEvents;
		if(gCCDParams.m_bSkipIfMoreThen>0 && gCCDParams.m_eCheckIfMorePoint!=eAfterCoic){
	      // requires rejection of all events in case in the redial of gCCDParams.m
   	   // there is more then m_bSkipIfMoreThen events :
      	int nRejectedIfMore = (in.Matrix)->RejectIfMoreThen( in.pPipeline );

			nGenEvents = (in.Matrix)->CountFoundEvents( genX, genY, 10 );
			if( nGenEvents>1 )
            nGenEvents = 1;
						
			stat.nIfMoreThen = ( (in.Matrix)->GetFoundEventCount() - nGenEvents);
			if( stat.nIfMoreThen < 0 ){
				stat.nIfMoreThen = 0;
			}
			sstat.nIfMoreThen = nGenEvents;
	   }
		stat.nTracks = MAX( ( (in.Matrix)->GetFoundEventCount() - nGenEvents ) , 0 );
		sstat.nTracks = nGenEvents;
		if( gCCDParams.m_bLogFrameStat ){
			backgrStat.push_back( stat );
			if( gCCDParams.GetMC() && pGenObj && pGenObj->m_pGenEvent){
				// add to sample only if really added :
				samplesStat.push_back ( sstat );
			}
		}

		// log event rates :
		if(gCCDParams.m_bLogFrameStat && !gCCDParams.m_bCCDDouble){
			LogEventRates( in.pPipeline );
		}
		


	}
	PROFILER_END("Full analysis of new frame took : ")
	PROFILER_PUT("CheckIfOverlaps took : ",total_in_check);

	if( (in.Matrix)->GetFoundEvents().size()>0 ){
		bRet = TRUE;
		(in.Matrix)->SetEventTime();
		(in.Matrix)->GetFoundEvents().SetEvtIdx();
	}

	if( gCCDParams.m_bCheckForSUPERNEW && (in.pPipeline)->m_BrightenList.size()>0 ){
		(in.pPipeline)->m_BrightenList.SetEventTime( (in.pPipeline)->m_FrameUnixTime );
	}

	if( gCCDParams.m_nSamplesToPutOnFrame>1 ){
		FillSamplesInfo( in, nNewLaplaceTreshold );
	}

	delete inForPrev;
	delete outForPrev;

	_TRACE_PRINTF_3("checked for SUPER_NEW %d times\n",super_new);

	PrintFrameAnalyseInfo( in, (in.Matrix)->GetFoundEvents().size() );

	return bRet;
}


int CCD_Analyser::CountClusterMAX( BIG_ELEM_TYPE** p_data, int xSize,
				   LONG_T* cluster,LONG_T cluster_cnt,
				   CPointList& max_list )
{
	int ret=0;
	for(int i=0;i<cluster_cnt;i++){
		int x = ( cluster[i] % xSize );
		int y = ( cluster[i] / xSize );

		int is_local_max = 1;
		for(int xx=(x-1);(xx<=(x+1) && is_local_max);xx++){
			for(int yy=(y-1);yy<=(y+1);yy++){
				if( xx!=x || yy!=y ){
					if( p_data[yy][xx] >= p_data[y][x] ){
						is_local_max=0;
						break;
					}
				}
			}
		}

		if( is_local_max ){
			max_list.Add( x, y );
			ret++;
		}
	}

	return ret;
}



int CCD_Analyser::CountClusterMAX( ELEM_TYPE** p_data, int xSize,
				   LONG_T* cluster,LONG_T cluster_cnt,
				   CPointList& max_list )
{
	int ret=0;
	for(int i=0;i<cluster_cnt;i++){
		int x = ( cluster[i] % xSize );
		int y = ( cluster[i] / xSize );

		int is_local_max = 1;
		for(int xx=(x-1);(xx<=(x+1) && is_local_max);xx++){
			for(int yy=(y-1);yy<=(y+1);yy++){
				if( xx!=x || yy!=y ){
					if( p_data[yy][xx] >= p_data[y][x] ){
						is_local_max=0;
						break;
					}
				}
			}
		}

		if( is_local_max ){
			max_list.Add( x, y );
			ret++;
		}
	}

	return ret;
}



