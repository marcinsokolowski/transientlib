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
#include <stdlib.h>
#include "ccd_defines.h"
#include <string.h>
#include <mystring.h>
#include <myutil.h>
#include "ccd_excp.h"
#include "ccd_log.h"

mystring GetWorkingModeDesc( eDAQMode_T mode, BOOL_T bScanMode  )
{
	mystring szRet = "UNKNOWN";

	if( bScanMode ){
		szRet = "SCAN";
		return szRet;
	}
		
	if(mode==eDAQSatTriggerMode)
		szRet = "TRIGGER";
	if(mode==eDAQTriggerMovingMode)
		szRet = "MOVING TO TRIGGER";
	if(mode==eDAQNormalMode)
		szRet = "NORMAL";
	if(mode==eDAQCollectMode)
		szRet = "COLLECTION";
	if(mode==eDAQDarkCollectionMode)
		szRet = "DARK COLLECTION";

	return szRet;	
}

eDAQMode_T GetNewDAQMode( eCCDRequestType_T req_mode )
{
	if(req_mode==eDAQReq_SatTriggerMode || req_mode==eDAQReq_OnTriggerPosition)
		return eDAQSatTriggerMode;
	if(req_mode==eDAQReq_MovingMove)
		return eDAQTriggerMovingMode;
	if(req_mode==eDAQReq_NormalMode)
		return eDAQNormalMode;

	return eDAQNormalMode;
}


const char* szEventTypeDescTab[4] = { "NONE", "FLASH", "SN" , NULL };

eEventType GetEventType( const char* szDescType )
{
	register int i=0;
	while( szEventTypeDescTab[i] ){
		if(strcmp( szEventTypeDescTab[i], szDescType)==0){
			return (eEventType)i;
		}
		i++;
	}
	return eNone;
}


CPixelAnalyseIn::CPixelAnalyseIn()
: treshold_for_max(0), x(0),y(0),pos(0),xSize(0),ySize(0), Matrix(NULL),
  p_data(NULL),p_data_fast(NULL),p_curr_data_laplace(NULL), p_homeo_data(NULL),
  p_fast_homeo_data(NULL),p_laplace_data(NULL),p_laplace_data_fast(NULL),
  pipeline_size_minus_1(0),PrevMatrixPtrCnt(0),pMaxLaplacePrev(NULL),
  p_max_prev_lap(NULL),ccd_index(0),frame_index(0),pPipeline(NULL),treshold_for_not_more_then_n_exceeds(0),
  pPixelList(NULL), p_curr_data_laplace_normal(NULL),pCCDInfo(NULL), pCamCfg(NULL),
  treshold_for_cluster(0),treshold_for_black_pixel(0),treshold_for_prev(0),
	day_frame_index(0),treshold_for_hot(0),SigmaLap(0),MeanLap(0),SigmaRaw(0),MeanRaw(0),
  p_prev_lap_fast(NULL)
{
}

void CPixelAnalyseIn::SetPrevMatrix( const CPixelAnalyseIn& in )
{
	PrevMatrixPtrCnt = in.PrevMatrixPtrCnt;
	memcpy( PrevMatrixPtr, in.PrevMatrixPtr, PrevMatrixPtrCnt*sizeof(PrevMatrixPtr[0]) );


	memcpy( PrevFramesShiftTab, in.PrevFramesShiftTab, PrevMatrixPtrCnt*sizeof(PrevFramesShiftTab[0]) );
	memcpy( PrevLaplacePtr, in.PrevLaplacePtr, PrevMatrixPtrCnt*sizeof(PrevLaplacePtr[0]) );

}

CSinglePixelData::CSinglePixelData()	
 : newSum(0),maxSum(0),homeoSum(0),m_bLaplaceMedianRejected(FALSE),
   laplaceSum(0),medianLaplaceSum(0),prevLaplaceSum(0),homeoLaplaceSum(0),
   otherSum(0),x0(0),y0(0),pos0(0),eventType(eNone),
	m_bLocalMaxRejected(FALSE),maxAverageOfPrev(0),newRawSum(0),prevLaplaceSumMin(0),
	prevLaplaceSumX(0),prevLaplaceSumY(0),m_nCountAbove(0),treshold_for_not_more_then_n_exceeds(0),
	m_bNotMoreThenNAboveRejected(FALSE),m_bOutSideAcceptanceRegion(FALSE),
	m_bRejectedByCustomEventAnal(FALSE),m_bRejectedByNextFrames(FALSE),
	m_bRejectedByEventShape(FALSE),m_Sphericity(0),m_bRejectedDueToTrack(FALSE),m_bRejectedByCoic(FALSE),
	m_bInTriggerMode(FALSE), PrevLaplaceValuesCount(0), m_Likehood(0), m_Significance(0),
	m_MaxClusterValue(0),m_fBlackRatio(0),m_PrevLapInPixel(0),PixelRawValue(0),laplaceOnNext(0),
	m_TrackType(eNormalTrack),m_bRejByTrackOnSingleCam(FALSE),minChi2(10000000000.000),
	minChi2On3(100000.000),bestTrackID(-1),minDist(0),veloCheckOK(1),rx(0),ry(0),
	MaxPixelRawValue(0)
{

}

void CSinglePixelData::Init()
{
	memset( this, '\0', sizeof(CSinglePixelData) );
}

CPixelAnalyseOut::CPixelAnalyseOut()
:ncnt(0),pixel_list(NULL),prevFramesCount(0),
 cluster_cnt(0),ClusterWithMoreCnt(0),newClusterSum(0), hitpixel_count(0),
 startPoints_cnt(0), neighb_count(0)
{	
	hitpixel_list = new LONG_T[MAX_CLUSTER_SIZE];
	startPoints   = new LONG_T[MAX_CLUSTER_SIZE];
	neighb_list   = new LONG_T[MAX_CLUSTER_SIZE];
	cluster       = new LONG_T[MAX_CLUSTER_SIZE];
	ClusterWithMore = new LONG_T[MAX_CLUSTER_SIZE];

	memset(neighb_list,'\0',MAX_CLUSTER_SIZE*sizeof(LONG_T));
	memset(ClusterWithMore,'\0',MAX_CLUSTER_SIZE*sizeof(LONG_T));
	memset(cluster,'\0',MAX_CLUSTER_SIZE*sizeof(LONG_T));
	memset(prevClusterSum,'\0',MAX_PREV_SUMS*sizeof(LONG_T));
	memset(hitpixel_list,'\0',MAX_CLUSTER_SIZE*sizeof(LONG_T));
	memset(startPoints,'\0',MAX_CLUSTER_SIZE*sizeof(LONG_T));
	memset(neighbours,'\0',MAX_NEIGHB_SIZE*sizeof(LONG_T));

	m_pCluster = new Table2D<BIG_ELEM_TYPE>(50,50);
}

CPixelAnalyseOut::~CPixelAnalyseOut()
{
	if(hitpixel_list){
		delete [] hitpixel_list;
	}
	if(startPoints){
		delete [] startPoints;
	}
	if(neighb_list){
		delete [] neighb_list;
	}
	if(cluster){
		delete [] cluster;
	}
	if(ClusterWithMore){
		delete [] ClusterWithMore;
	}
	if(m_pCluster)
		delete m_pCluster;
}

CPixelAnalyseOut::CPixelAnalyseOut( const CPixelAnalyseOut& right )
{
	(*this) = right;
}

CPixelAnalyseOut& CPixelAnalyseOut::operator=( const CPixelAnalyseOut& right )
{
	ncnt = right.ncnt;
	memcpy( neighb_list, right.neighb_list, ncnt*sizeof(neighb_list[0]) );

	pixel_list = right.pixel_list;

	newClusterSum = right.newClusterSum;
	memcpy( prevClusterSum, right.prevClusterSum, MAX_PREV_SUMS*sizeof(prevClusterSum[0]) );

	prevFramesCount = right.prevFramesCount;

	memcpy( &m_PixelOut, &(right.m_PixelOut), sizeof(m_PixelOut) );

	cluster_cnt = right.cluster_cnt;
	memcpy( cluster, right.cluster, cluster_cnt*sizeof(cluster[0]) );

	ClusterWithMoreCnt = right.ClusterWithMoreCnt;
	memcpy( ClusterWithMore, right.ClusterWithMore, ClusterWithMoreCnt*sizeof(ClusterWithMore[0]) );		

	hitpixel_count = right.hitpixel_count;
	memcpy( hitpixel_list, right.hitpixel_list, hitpixel_count*sizeof(hitpixel_list[0]) );

	startPoints_cnt = right.startPoints_cnt;
	memcpy( startPoints, right.startPoints, startPoints_cnt*sizeof(startPoints[0]) );

	neighb_count = right.neighb_count;
	memcpy( neighbours, right.neighbours, neighb_count*sizeof(neighbours[0]) );

	return (*this);
}

void CPixelAnalyseOut::Dump()
{
	mystring szOutDesc;
	
	szOutDesc << "NeighbCount = " << ncnt << "\n";
	szOutDesc << "newLaplace = " << m_PixelOut.laplaceSum << "\n";
	szOutDesc << "prevLaplace = " << m_PixelOut.prevLaplaceSum << "\n";
	printf("%s\n",szOutDesc.c_str());
}


/*BOOL_T ADD_POINT_func( LONG_T* table, LONG_T& idx, int pos )
{
	if(idx>=MAX_CLUSTER_SIZE){
		// throw CCDClusterSizeExcp();
		return FALSE;
	}
	table[idx] = pos;
	idx++;
	return TRUE;
}*/




CFrameIdentStat::CFrameIdentStat( int frame_index )
: m_FrameIndex( frame_index )
{
	Init();
}

void CFrameIdentStat::Init()
{
	nTotal=0;
	nMaxAllowedValue=0;
	nTnewCut=0;
   nTprevCut=0;

   // rejected after :
   nShapeCut=0;
   nBlackCut=0;
   nHotCut=0;
   nLocalMaxCut=0;
	nMinPrevCut=0;
	bAllRejected=FALSE;		
	nIfMoreThen=0;
	nOverlap = 0;
	nAfterCoic = 0;
	nTracks = 0;
	bAllRejAfterCoic = FALSE;
	nSatRej = 0;
   nStarRej = 0;
	nFinalEvents = 0;
	nClusterOverlap = 0;
}

void CFrameIdentStat::DecTrackAndCoicEvents()
{
	DecZero( nTracks );
	DecZero( nAfterCoic );
}

void CFrameIdentStat::LogToFile( CCDPipeline* pPipeline, const char* fname,
											const char* szSubDir /*="Events"*/ )
{
	CCDLog log("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t\t%d\t%d\t%d\n",
					"Total\tMaxAll\tT_new\tT_prev\tMinPrev\tLocMax\t"\
					"Overlap\tShape\tBlack\tHot\tIfMore\tAllRej\t"\
					"Tracks\tCoic\tAllRejAftCoic\tSatRej\tStarRej\tFinal");

	int all_rej = (int)bAllRejected;
	int all_rej_after_coic = (int)bAllRejAfterCoic;
   log.DumpToFile2( pPipeline, m_FrameIndex, szSubDir, fname, 
						nTotal, nMaxAllowedValue, nTnewCut, nTprevCut, 
						nMinPrevCut,nLocalMaxCut, nOverlap, nShapeCut, 
						nBlackCut, nHotCut, nIfMoreThen, all_rej,
						nTracks,nAfterCoic, all_rej_after_coic,
						nSatRej , nStarRej, nFinalEvents );
	
}

CFrameEventStatTab::CFrameEventStatTab()
{
	
}

CFrameEventStatTab::~CFrameEventStatTab()
{
	
}

int CFrameEventStatTab::IncFinal( int frame_index )
{
	int ret=0;
	vector<CFrameIdentStat>::iterator i;
	for(i=begin();i!=end();i++){
		if( i->m_FrameIndex == frame_index ){
			i->nFinalEvents++;
			ret++;
		}
	}
	return ret;
}

void CFrameEventStatTab::DecTrackEvents( int frame_index )
{
	vector<CFrameIdentStat>::iterator i;
	for(i=begin();i!=end();i++){
		if( i->m_FrameIndex == frame_index ){
			i->nTracks--;
		}
	}
}

void CFrameEventStatTab::DecTrackAndCoicEvents( int frame_index )
{
   vector<CFrameIdentStat>::iterator i;
   for(i=begin();i!=end();i++){
      if( i->m_FrameIndex == frame_index ){
			DecZero( i->nTracks );
         DecZero( i->nAfterCoic );
      }
   }
}


BOOL_T IsCollectionMode( eDAQMode_T mode ){ 
	return ( mode==eDAQCollectMode || mode==eDAQDarkCollectionMode );
}


CVariableInfo::CVariableInfo()
{
	Init();
}

void CVariableInfo::Init()
{
	l_new=0;
	l_prev=0;
	nAfterTv=0;
	nOverlaps=0;				
	fSpericity=0;
	fBlackRatio=0;
	averInPixel=0;
}

const char* GetCamStatusDesc( eCamState_T stat )
{
	switch( stat ){
		case eCamState_OK:
			return "OK";
		case eDoCleanMeasureErr:
			return "DoCleanMeasure error";
		case eDoMeasureErr:
			return "DoMeasure error";
		case eDoCleanErr:
			return "DoClean error";
		case eGetDataErr:
			return "GetData error";
		case eDoStartFrameErr:
			return "StartFrame error";
		case eDoReadChipToRAMOnCamErr:
			return "DoReadChipToRAMOnCam error";
		case eRefreshFullStatusErr:
			return "RefreshFullStatus error";
		case eCriticalErr:
			return "Critical error - could not continue !!!";
		default :
			return "Undefined";
	}
	return "Undefined";
}
