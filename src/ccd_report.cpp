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
#include "ccd_report.h"
#include "ccd_analyse.h"
#include "ccd_globals.h"
#include "ccd_trace.h"
#include "ccd_pipeline.h"
#include <mystring.h>
#include <myfile.h>
#include <cfg.h>
#include <mymacros.h>
#include <mydate.h>
#include <cexcp.h>
#include <math.h>
#include <myutil.h>
#include <Astroangle.h>
#include <mathfunc.h>
#include "ccd_log.h"
#include "ccd_errors.h"
#include <myfits.h>
#include "dbstructdef.h"
#include <satlibdefs.h>
#include <myparser.h>
#include "satcone.h"

eLogFileFormat CCDEventList::m_DefaultFormat=eVerif;
BOOL_T CccdReport::m_bIgnoreReadError=FALSE;
BOOL_T CCDEventList::m_bInitialized=FALSE;
int CccdReport::m_DayEventCounter=-1;
mystring CccdReport::m_EventDescOutFmt = "%d %d %.2f %.2f %d %d %d %d (%.2f,%.2f)-(%.2f,%.2f) %s %s %s %d %d %d %d (%s,%+.4f) %d %.2f (%.2f,%.4f) %.2f %s %d %.2f %.6f %d (%d-%.2f-%.2f %d,%.2f,%.2f) %d %d %s %d %d %.2f\n";
mystring CccdReport::m_EventDescInFmt = "%d %d %lf %lf %d %d %d %d (%lf,%lf)-(%lf,%lf) %s %s %s %d %d %d %d %s %d %lf (%lf,%lf) %lf %s %d %lf %lf %d (%d-%lf-%lf %d,%lf,%lf) %d %d %s %d %d %lf\n";
mystring CccdReport::m_EventDescHeader = "# DayFrame Cam# X Y Gen Ident RejOnNext Track (x0,y0)-(x1,y1) Mag Type DateTime InTrigger RejByCoic Lap_new Lap_prev (ra,dec) Lap_new_max Signi CCD_prim CoicR Type RunFrame CoicR_sec Spher SingleCamTrack bestTrack EvtIdx Star# External RawValue MaxRawValue ConeDist[km]\n";

// final event report 
mystring CccdReport::m_FinalEventHeader = 
  "# ID# DF# F# Evt# x y L_new L_prev L_max_new S_raw " \
	"Spher Signi ClusterSize CoicR_pix AverInPixel "\
	"BlackRatio (RA,DEC) Type CoicR_sec (RA,DEC) UtTime L_next bestTrack EventUT Star# External ConeDist[km]\n";
mystring CccdReport::m_FinalEventOutFmt = 
	"%s %d %d %d %d %d %d %d %d %d %.4f %.4f %d %.2f %d %.2f (%s,%+.4f) %s %.2f (%.8f,%.8f) %d %d (%d-%.2f-%.2f %d,%.2f,%.2f) %s %d %s %.2f\n";
mystring CccdReport::m_FinalEventInFmt = 
	"%s %d %d %d %lf %lf %d %d %d %d %lf %lf %d %lf %d %lf (%dh%dm%lfs,%lf) %s %lf (%lf,%lf) %d %d (%s %d,%lf,%lf) %s %d %s %lf\n";


// additional grb info :
mystring CccdReport::m_GrbInfoHeader ="GRB_ID GCN_ID SOURCE_ID GRB_NAME (GRB_RA,GRB_DEC) \n";
mystring CccdReport::m_GrbInfoFmt    ="%d %d %d %s (%s,%.8f) ";


// tracks :
mystring CTrackDesc::m_LogHeader="# DF# x y track chi2 vx vy rx ry OK minDist why\n";
mystring CTrackDesc::m_LogFmt="%d %d %d %d %.2f %.2f %.2f %.2f %.2f %d %.2f %s\n";


CEventBaseInfo::CEventBaseInfo( const CEventBaseInfo& right )
{
	(*this) = right;
}

CEventBaseInfo& CEventBaseInfo::operator=( const CEventBaseInfo& right )
{
   m_MaxPoint = right.m_MaxPoint;
   m_Point = right.m_Point;
   m_PointTransformed = right.m_PointTransformed;
   m_CoicRadius = right.m_CoicRadius;
   m_CoicRadiusInRad = right.m_CoicRadiusInRad;
   m_FrameIndex = right.m_FrameIndex;
   m_DayFrameIndex = right.m_DayFrameIndex;
   m_Time = right.m_Time;
	m_AstroCoord = right.m_AstroCoord;
	m_bSavedToDB = right.m_bSavedToDB;

	return (*this);
}

CAstroCoord& CAstroCoord::operator=(const CAstroCoord& right){   
	m_Azimuth = right.m_Azimuth;
	m_Altitude = right.m_Altitude;
	Dec = right.Dec;
	RA = right.RA;
	HA = right.HA;

	return (*this);
} 


CEventBaseInfo::CEventBaseInfo( int max_point_x, int max_point_y, int point_x, 
										  int point_y, int frame_index, int day_frame_index )
: m_MaxPoint( max_point_x, max_point_y ), m_Point( point_x, point_y ),m_FrameIndex( frame_index ), 
  m_PointTransformed( max_point_x, max_point_y ), m_CoicRadius(0),
  m_DayFrameIndex( day_frame_index ), m_CoicRadiusInRad(0),m_Time(0),
  m_bSavedToDB(FALSE)
{
}

void CccdReport::CalcCoicRadiusInRad( double ra, double dec )
{
	m_CoicRadiusInRad = sqrt( (ra-m_AstroCoord.RA)*(ra-m_AstroCoord.RA) + 
								  	  (dec-m_AstroCoord.Dec)*(dec-m_AstroCoord.Dec) );
}

void CEventBaseInfo::SetTransformedPoint( double x, double y, BOOL_T bCalcCoicRadius )
{
	m_PointTransformed.x = x;
	m_PointTransformed.y = y;

	if( bCalcCoicRadius ){
		m_CoicRadius = sqrt ( CMyMathFunc::mysqr( x-m_MaxPoint.x) + CMyMathFunc::mysqr( y-m_MaxPoint.y) );	
	}
}

int CEventBaseInfo::GetMemSize( BOOL_T bShowAll )
{
	return sizeof(CEventBaseInfo);
}

void CccdReport::AddGRBInfo( CGRBInfo* pGrbInfo )
{
	m_EventType = EVENT_TYPE_GRB_IN_DB;
	m_szGRBInfo << pGrbInfo->grb_id << " " 
		    << pGrbInfo->gcn_id << " " 
		    << pGrbInfo->source_id;
}


int CccdReport::GetMemSize( BOOL_T bShowAll )
{
	int ret = sizeof(CccdReport);
	ret += m_ClusterCount*sizeof(LONG_T);
	

	return ret;
}

void CccdReport::GetInternalTriggerID( mystring& szID )
{
	char szTmp[64];
	sprintf(szTmp,"%.5d%.5d",m_DayFrameIndex,EvtIdx);
	szID = szTmp;
}

void CccdReport::Init()
{	
	if( m_SatInfo ){
		memset(m_SatInfo,'\0',sizeof(struct satInfo));
	}
	m_szExternal = "-";
	m_nFrameStarCount = 0;
	m_SatCoicRadiusInSec = 0;
	m_Time = 0;
	EvtIdx = 0;
	m_CameraIdx = 0;
	m_bGenerated = FALSE;
	m_bIdentified = FALSE;
   m_MaxPoint.Set(-1,-1);
	m_Point.Set(0,0);
	m_ClusterCenter.Set(0,0);
	m_LowLeft.Set(NOT_DEFINED,NOT_DEFINED);
	m_TopRight.Set(NOT_DEFINED,NOT_DEFINED);
	m_Magnitude = "";
	m_Cluster = NULL;
	m_ClusterCount = 0;
	prevFramesCount = 0;
	m_FrameIndex = 0;
	m_DayFrameIndex= 0;
	m_bFirstLevelAccepted = FALSE;
	nextFramesCount = 0;
	m_bCheckOnNext = FALSE;
	m_nAddInfo = 0;
	m_MaxValue = 0;
	newClusterSum = 0;
	m_LastDumpedFrame = 0;
	m_bRemoveFromList = FALSE;
	m_EventType = EVENT_TYPE_FLASH;
	m_SatType = 'N';
	visibility = SAT_UNKNOWN;

	m_PixelAnalResults.Init();
	m_MinDistOutCone = 0;
}

CccdReport::CccdReport()
: CEventBaseInfo(-1,-1,0,0,0), m_CameraIdx(0),m_bGenerated(FALSE),m_bIdentified(FALSE),
  m_ClusterCenter(0,0), m_LowLeft(NOT_DEFINED,NOT_DEFINED), m_TopRight(NOT_DEFINED,NOT_DEFINED),
  m_Magnitude(""),m_Cluster(NULL),m_ClusterCount(0),
  prevFramesCount(0),m_bFirstLevelAccepted(FALSE),
  nextFramesCount(0),m_bCheckOnNext(FALSE),
  m_nAddInfo(0),m_MaxValue(0),newClusterSum(0),m_LastDumpedFrame(0),m_bRemoveFromList(FALSE),
  m_EventType(EVENT_TYPE_FLASH),m_SatType('N'),EvtIdx(0),m_PipelineIndex(0),
  visibility(SAT_UNKNOWN),m_nFrameStarCount(0),m_SatInfo(NULL),
  m_MinDistOutCone(0)
{
	m_SatInfo = new satInfo();
	memset(m_SatInfo,'\0',sizeof(struct satInfo));
}

CccdReport::CccdReport( LONG_T camNum, const CPixelAnalyseIn& in,
                        const CPixelAnalyseOut& out,
                        BOOL_T bGen,BOOL_T bIdentified,
                        CPoint* pLowLeft, CPoint* pTopRight,
                        const char* mag, BOOL_T bFirstLevelAccepted )
: CEventBaseInfo(in.x,in.y,in.x,in.y,in.frame_index),m_CameraIdx(camNum),m_bGenerated(bGen),m_bIdentified(bIdentified),
  m_ClusterCenter(out.m_PixelOut.x0,out.m_PixelOut.y0),
	m_LowLeft(NOT_DEFINED,NOT_DEFINED),
  m_TopRight(NOT_DEFINED,NOT_DEFINED),
  m_Magnitude(mag),m_Cluster(NULL),m_ClusterCount(out.cluster_cnt),
  prevFramesCount(out.prevFramesCount), m_bFirstLevelAccepted(bFirstLevelAccepted),
  nextFramesCount(0),m_bCheckOnNext(FALSE),
  m_nAddInfo(0),m_MaxValue(0),newClusterSum(0),m_LastDumpedFrame(0),
  m_bRemoveFromList(FALSE),m_EventType(EVENT_TYPE_FLASH),EvtIdx(0),m_SatType('N'),
  m_PipelineIndex(0),visibility(SAT_UNKNOWN),m_nFrameStarCount(0),m_SatInfo(NULL),
  m_MinDistOutCone(0)
{
	m_SatInfo = new satInfo();
	memset(m_SatInfo,'\0',sizeof(struct satInfo));

	if( in.pPipeline ){
		m_DayFrameIndex = (in.pPipeline)->GetDayFrameCounter();
	}else{
		m_DayFrameIndex = in.frame_index;
	}

	if(pLowLeft)
		m_LowLeft = (*pLowLeft);
	if(pTopRight)
		m_TopRight = (*pTopRight);
	if(out.cluster && out.cluster_cnt){
		m_ClusterCount = out.cluster_cnt;
		m_Cluster = new LONG_T[ m_ClusterCount ];
		memcpy(m_Cluster,out.cluster,m_ClusterCount*sizeof(LONG_T));
	}else{
		m_ClusterCount = 0;
		m_Cluster = NULL;
	}
	SetPrevClusterSum( out.prevFramesCount, out.prevClusterSum );
   memset(nextClusterSum,'\0',MAX_PREV_SUMS*sizeof(LONG_T));	
	// memcpy(m_PixelAnalResults.PrevFramesMaxVal,in.PrevFramesMaxVal,(gCCDParams.m_nMaxOfAverageOfPrevN+1)*sizeof(double));


	memcpy(m_PixelAnalResults.PrevFramesX,in.PrevFramesX,(gCCDParams.m_nMaxOfAverageOfPrevN+1)*sizeof(int));
	memcpy(m_PixelAnalResults.PrevFramesY,in.PrevFramesY,(gCCDParams.m_nMaxOfAverageOfPrevN+1)*sizeof(int));
}

CccdReport::CccdReport(const CccdReport& right)
: m_Cluster(NULL),m_ClusterCount(0),m_SatInfo(NULL)
{	
   (*this) = right;
}

CccdReport::~CccdReport()
{
   if(m_Cluster)
      delete [] m_Cluster;
	if( m_SatInfo )
		delete m_SatInfo;
}

BOOL_T CccdReport::IsIdentified() const
{
	return ( m_bIdentified && !m_PixelAnalResults.m_bRejectedDueToTrack && 
			   !m_PixelAnalResults.m_bRejectedByNextFrames && 
				m_EventType==EVENT_TYPE_FLASH );
}

BOOL_T CccdReport::IsIdentifiedForTrack() const 
{
	return ( m_bIdentified && !m_PixelAnalResults.m_bRejectedDueToTrack && 
			   !m_PixelAnalResults.m_bRejectedByNextFrames && 
				( m_EventType==EVENT_TYPE_FLASH || m_EventType==EVENT_TYPE_SAT || m_EventType==EVENT_TYPE_STAR ) );
}


const char* CccdReport::GetAddInfoDesc( eAddInfoType_T add_info_type )
{
	if( add_info_type == eClusterChar )
		return "Cluster characteristics";

	if( add_info_type == eCrossChar )
		return "Cross characteristics";
	
	if( add_info_type ==  eCrossAboveTreshChar )
		return "Cross characteristics with treshold condition";

	return "UNKNOWN";
}

void CccdReport::SetNextFramesSum( LONG_T next_frames_cnt, const LONG_T* _nextClusterSum )
{
	nextFramesCount = next_frames_cnt;
	if(nextFramesCount)
		memcpy(nextClusterSum,_nextClusterSum,sizeof(LONG_T)*MAX_PREV_SUMS);
}

void CccdReport::RejectByNextFrames()
{
	if(gCCDParams.m_bRejectNotConfirmed)
		m_bIdentified = FALSE;
	m_PixelAnalResults.m_bRejectedByNextFrames=TRUE;
}

void CccdReport::CopyCluster( const CccdReport& right)
{
	SetCluster( right.m_Cluster, right.m_ClusterCount );
}

void CccdReport::SetCluster( const LONG_T* cluster, LONG_T cluster_cnt )
{
	m_ClusterCount = cluster_cnt;
	if(m_Cluster){
		delete [] m_Cluster;
		m_Cluster = NULL;
	}
	if( cluster && cluster_cnt ){
		m_Cluster = new LONG_T[m_ClusterCount];
		memcpy(m_Cluster,cluster,m_ClusterCount*sizeof(LONG_T));
		// }else{
		// printf("odo");
	}
}

CccdReport& CccdReport::operator=(const CccdReport& right)
{
	if( !m_SatInfo ){
		m_SatInfo = new satInfo();
	}
	if( m_SatInfo ){
		memset(m_SatInfo,'\0',sizeof(struct satInfo));
		if( right.m_SatInfo ){
			memcpy( m_SatInfo, right.m_SatInfo, sizeof(struct satInfo) );
		}
	}

	((CEventBaseInfo&)(*this)) = ((CEventBaseInfo&)right);
	m_szExternal = right.m_szExternal;
	visibility = right.visibility;
	m_EvtID = right.m_EvtID;
	m_PipelineIndex = right.m_PipelineIndex;
	m_Time = right.m_Time;
	EvtIdx = right.EvtIdx;
	m_CoicRadius = right.m_CoicRadius;
	m_CoicRadiusInRad = right.m_CoicRadiusInRad;
	m_DayFrameIndex = right.m_DayFrameIndex;
	m_FrameIndex = right.m_FrameIndex;
	m_CameraIdx = right.m_CameraIdx;
	m_MaxPoint = right.m_MaxPoint;
	m_MaxValue = right.m_MaxValue;
	m_Point = right.m_Point;
	m_ClusterCenter = right.m_ClusterCenter;
	m_GenPoint = right.m_GenPoint;
	m_bGenerated = right.m_bGenerated;
	m_bIdentified = right.m_bIdentified;
	m_LowLeft = right.m_LowLeft;
	m_TopRight = right.m_TopRight;
	m_Magnitude = right.m_Magnitude;
	m_bFirstLevelAccepted = right.m_bFirstLevelAccepted;
	m_bCheckOnNext = right.m_bCheckOnNext;
	m_LastDumpedFrame = right.m_LastDumpedFrame;
	// m_AstroCoord = right.m_AstroCoord;
	m_PointTransformed = right.m_PointTransformed;
	

	m_EventType = right.m_EventType;
	m_szSatName = right.m_szSatName;
	m_SatType = right.m_SatType;
	m_SatCoicRadiusInSec = right.m_SatCoicRadiusInSec;

	m_nFrameStarCount = right.m_nFrameStarCount;

	CopyCluster( right );

   newClusterSum = right.newClusterSum;
	SetPrevClusterSum( right.prevFramesCount, right.prevClusterSum );


	// next frames :
	SetNextFramesSum( right.nextFramesCount, right.nextClusterSum );	


	m_nAddInfo = right.m_nAddInfo;
	if(m_nAddInfo){
		memcpy(m_AdditionalInfo,right.m_AdditionalInfo,sizeof(m_AdditionalInfo[0])*m_nAddInfo);		
	}	


	memcpy( &m_PixelAnalResults, &(right.m_PixelAnalResults), sizeof(m_PixelAnalResults) );

	m_MinDistOutCone = right.m_MinDistOutCone;

	return (*this);
}

void CccdReport::SetPrevClusterSum( LONG_T prev_frames_cnt, const LONG_T* _prevClusterSum )
{
	prevFramesCount = prev_frames_cnt;
	if(prev_frames_cnt){
		memcpy(prevClusterSum,_prevClusterSum,sizeof(prevClusterSum[0])*MAX_PREV_SUMS);		
	}else{
		memset(prevClusterSum,'\0',MAX_PREV_SUMS*sizeof(prevClusterSum[0]));
	}
}

void CccdReport::SetIdentificationData( const CccdReport& right)
{
	memcpy( &m_PixelAnalResults, &(right.m_PixelAnalResults), sizeof(m_PixelAnalResults) );
	if(right.m_ClusterCount){
		SetCluster( right.m_Cluster, right.m_ClusterCount );
		newClusterSum = right.newClusterSum;
		SetPrevClusterSum(  right.prevFramesCount , right.prevClusterSum );
		m_bFirstLevelAccepted = TRUE;
	}
}

void CccdReport::SetIdentificationData( const CPixelAnalyseIn& in, const CPixelAnalyseOut& out )
{
	// printf("newLaplace = %d, prevLaplaceSum = %d\n",laplaceSum,prevLaplaceSum);	

	if(out.cluster_cnt){
		SetCluster( out.cluster, out.cluster_cnt );
		newClusterSum = out.newClusterSum;
		SetPrevClusterSum(  out.prevFramesCount , out.prevClusterSum );
		m_bFirstLevelAccepted = TRUE;
	}
	
	if(m_bGenerated){
		if((long)m_Point.x != (long)in.x || (long)m_Point.y != (long)in.y ){
			m_Point.x = in.x;
			m_Point.y = in.y;
		}
	}

	memcpy( &m_PixelAnalResults, &(out.m_PixelOut), sizeof(m_PixelAnalResults) );	

	if( in.PrevValues.counter>0 ){
		m_PixelAnalResults.PrevLaplaceValuesCount = in.PrevValues.counter;
		for(register int i=0;i<in.PrevValues.counter;i++){
			m_PixelAnalResults.PrevLaplaceValues[i] = in.PrevValues.values[i];
		}
	}

}


void CccdReport::GetXYRange( LONG_T& x_size, LONG_T& y_size,
                             LONG_T xSize, LONG_T ySize )
{
	x_size=0;
	y_size=0;

	for(int i=0;i<m_ClusterCount;i++){
		LONG_T pos = m_Cluster[i];
		LONG_T x = (pos % xSize);
		LONG_T y = (pos / xSize);

		if(fabs(x-m_ClusterCenter.x)>x_size)
			x_size = (long)fabs(x-m_Point.x);
		if(fabs(y-m_ClusterCenter.y)>y_size)
         y_size = (long)fabs(y-m_Point.y);
	}
}

void CccdReport::CalcAstroCoordinates( CCDMatrix& frame, CCDProcState* pFrameInfo )
{
	pFrameInfo->CalcAstroCoord( my_round(m_MaxPoint.x), my_round(m_MaxPoint.y), 
										 m_AstroCoord.m_Azimuth, m_AstroCoord.m_Altitude, 
										 m_AstroCoord.Dec, m_AstroCoord.RA, m_AstroCoord.HA );
}

double CccdReport::FindMinDist( vector<CccdReport*>& eventList, CccdReport& evt )
{
	vector<CccdReport*>::iterator i;
	double minDist=100000.000;
	for(i=eventList.begin();i!=eventList.end();i++){
		double dist = CPoint::dist( (*i)->m_Point , evt.m_Point );
		if(dist < minDist){
			minDist = dist;
		}
	}
	return minDist;
}

CCDEventList::CCDEventList(int nAutoAlloc)
// :vector<CccdReport>(1000),format(eUnknownLogFormat)
{
	if(nAutoAlloc>0)
		reserve(nAutoAlloc);
	Initialize();
}

void CCDEventList::SetEventTime( time_t ut_time )
{
	CCDEventList::iterator i;
	for(i=begin();i!=end();i++){
		i->m_Time = ut_time;
	}
}

void CCDEventList::SetEvtIdx()
{
	CCDEventList::iterator i;
	int k=0;
	for(i=begin();i!=end();i++){
		i->EvtIdx = k;
		k++;
	}
}

void CCDEventList::SetEvtNo()
{
	CCDEventList::iterator i;
	int k=1;
	char szTmp[128];
	for(i=begin();i!=end();i++){
		// i->EvtIdx = k;
		mystring szDT;
		get_short_night_date_local( i->m_Time, szDT );
		sprintf(szTmp,"%s%.5d",szDT.c_str(),k);
		i->m_EvtID = szTmp;
		k++;
	}
}

int CCDEventList::GetEventsSinceFrame( CCDEventList& out_list, int start_frame )
{
	out_list.clear();
	CCDEventList::iterator i;

	for(i=begin();i!=end();i++){
		if(i->m_DayFrameIndex>=start_frame){
			out_list.push_back( *i );
		}
		if( start_frame>=1133 ){
			printf("CCDEventList::GetEventsSinceFrame %d-(%d,%d)\n",i->m_DayFrameIndex,
					(int)i->m_MaxPoint.x,(int)i->m_MaxPoint.y);
		}
	}
	return out_list.size();
}

int CCDEventList::GetEvents( CCDEventList& out_list, int start_frame, int end_frame, BOOL_T bOnlyFinal )
{
	out_list.clear();
	CCDEventList::iterator i;

	for(i=begin();i!=end();i++){
		if(i->m_DayFrameIndex>=start_frame && i->m_DayFrameIndex<=end_frame){
			BOOL_T bAdded=FALSE;
			if( !bOnlyFinal || i->IsIdentified() ){
				bAdded=TRUE;
				out_list.push_back( *i );
			}
			printf("CCDEventList::GetEventsSinceFrame %d-(%d,%d) added=%d\n",i->m_DayFrameIndex,
					(int)i->m_MaxPoint.x,(int)i->m_MaxPoint.y,bAdded);
			if( i->m_DayFrameIndex == 1133 ){
				printf("%d %d %d %c\n",i->m_bIdentified,
						(i->m_PixelAnalResults).m_bRejectedDueToTrack,
						(i->m_PixelAnalResults).m_bRejectedByNextFrames,
						i->m_EventType);
			}
		}
	}
	return out_list.size();
}

int CCDEventList::GetEvents( CCDEventList& out_list, int day_frame )
{
	out_list.clear();
	CCDEventList::iterator i;

	for(i=begin();i!=end();i++){
		BOOL_T bAdded=FALSE;
		if(i->m_DayFrameIndex==day_frame){
			out_list.push_back( *i );
			bAdded=TRUE;
		}
		printf("CCDEventList::GetEventsSinceFrame %d-(%d,%d) added=%d\n",i->m_DayFrameIndex,
				(int)i->m_MaxPoint.x,(int)i->m_MaxPoint.y,bAdded);
	}
	return out_list.size();
}

void CCDEventList::Initialize()
{	
	if(!m_bInitialized){
		m_bInitialized = TRUE;
	}
}

CCDEventList::CCDEventList( const CCDEventList& right )
{
	(*this) = right;
}

CCDEventList& CCDEventList::operator=( const CCDEventList& right )
{
	format = right.format;
	clear();
	(*this) += right;
	return (*this);
}

CCDEventList& CCDEventList::operator+=( const CCDEventList& right )
{
	for(register int i=0;i<right.size();i++){
		push_back( right[i] );
	}
	return (*this);
}

void CCDEventList::CalcConeDist()
{
	CCDEventList::iterator i;
   for(i=begin();i!=end();i++){
		double ra_in_deg = AstroAngle::rad2deg( (i->m_AstroCoord).RA );
		double dec_in_deg = AstroAngle::rad2deg( (i->m_AstroCoord).Dec );

		CSatCone::calc_min_dist( ra_in_deg, dec_in_deg, i->m_Time, i->m_MinDistOutCone );
	}
}

BOOL_T CCDEventList::RejectEventDueToTrack( int x, int y )
{
	CCDEventList::iterator i;
	for(i=begin();i!=end();i++){	
		if( ((int)i->m_Point.x)==x && ((int)i->m_Point.y)==y ){
			(i->m_PixelAnalResults).m_bRejectedDueToTrack = TRUE;
			return TRUE;
		}
	}
	return FALSE;
}

CCDEventList& CCDEventList::AddNotFound( const CCDEventList& genEvents )
{
	for(int i=0;i<genEvents.size();i++){
		if(!genEvents[i].m_bIdentified)
			Add( genEvents[i] );
   }
   return (*this);
}

void CCDEventList::Add( const CccdReport& newEvent )
{
	push_back( newEvent );
}

BOOL_T CCDEventList::IsAnyIdentifiedNotGen()
{
	for(int i=0;i<size();i++){
		if((*this)[i].m_bIdentified && !(*this)[i].m_bGenerated)
			return TRUE;
	}
	return FALSE;
}

BOOL_T CCDEventList::WasAnyIdentified()
{

	for(int i=0;i<size();i++){
		if((*this)[i].m_bIdentified)
			return TRUE;
	}
	return FALSE;
}

void CCDEventList::DumpAllEvents()
{
	if(gCCDParams.m_bDumpAllEvents){
		if(strlen(gCCDParams.m_szRunEventsLog.c_str())){
   	   mystring szFname;
      	szFname << gCCDParams.GetOutputDir() << "/" << gCCDParams.m_szRunEventsLog;
	      printf("Dumping all run events to file : %s\n",szFname.c_str());
			MYTRACE2(gCCDTrace,"Dumping all events to file : " << szFname.c_str());
      	DumpEventReport( szFname.c_str(), *this );
	   }
	}
}

void CccdReport::GetDescHeader( mystring& szHeader )
{
	szHeader = m_EventDescHeader;
}

void CccdReport::GetDetailEventDesc( mystring& szEventDesc, long Index, long xSize, CCDPipeline* pPipeline )
{
	BOOL_T bFirstLevelAcc = FALSE;
	mystring szMethod;

	szEventDesc << "\n\nFrame# " << m_FrameIndex << "\n";
	szEventDesc << "Day Frame# " << m_DayFrameIndex << "\n";
	if(m_bGenerated && !m_bIdentified){
		szEventDesc << "GENERATED event was rejected\n";
	}
	szEventDesc << "Event# " << Index << "\n";
	szEventDesc << "At x=" << m_Point.x << ", y=" 
					<< m_Point.y << "\n";

	mystring szRA = AstroAngle::toString( m_AstroCoord.RA, ANGLE_RA_TYPE ).c_str();
	szEventDesc << "(RA,DEC) = (" << szRA.c_str()
					<< " , " << AstroAngle::rad2deg( m_AstroCoord.Dec ) << "\n";
	if(m_ClusterCount){
		szEventDesc << "Cluster center of mass at x_center=" << m_ClusterCenter.x << ", y_center=" 
						<< m_ClusterCenter.y << "\n";
	}

	LONG_T ncnt = gCCDParams.m_nNeighbToSumCount;

	if(gCCDParams.m_bAnalyseMaxFromPrevInPipeline){
		szMethod << " PREV_SUM ";
	}
	if(gCCDParams.m_bCheckHomeoRawCond){
		szMethod << " HOMEOPATIC ";
	}				
	if(gCCDParams.m_bCheckLaplaceCondition){
		szMethod << " LAPLACE ";
	}

	if(m_bIdentified){
		szEventDesc << "Event type is : " << szEventTypeDescTab[(LONG_T)m_PixelAnalResults.eventType] << "\n";
	}
	szEventDesc << "At first level trigger " << szMethod << " method was used, using " << ncnt << " pixels \n";	
	szEventDesc << "RAW : newRawSum=" << m_PixelAnalResults.newRawSum 
               << ", otherSum=" << m_PixelAnalResults.otherSum << "\n";

	if(gCCDParams.m_bCheckLaplaceCondition){
		szEventDesc << "LAPLACE : newLaplace=" << m_PixelAnalResults.laplaceSum << ", prevLaplace=" << m_PixelAnalResults.prevLaplaceSum << "\n";
	}

	szEventDesc << "condition was : ";		
	// mystring szCond;
	// GetCondDesc( m_EventType, newSum, otherSum, ncnt )

	if( m_PixelAnalResults.eventType != eBrighten){
		if(gCCDParams.m_bCheckHomeoRawCond || gCCDParams.m_bAnalyseMaxFromPrevInPipeline){
			szEventDesc << "newSum>" << (ncnt*gCCDParams.m_TresholdPerPixel)
   		            << " AND otherSum<" << (ncnt*gCCDParams.m_MaxPrevPerPixel);
		}
		if(gCCDParams.m_bCheckLaplaceCondition){
			szEventDesc << "newLaplace>" << pPipeline->GetPipelineCfg().m_nNewLaplace 
                     << " AND otherLaplace<" << pPipeline->GetPipelineCfg().m_nMaxLaplaceOnOther
							<< " AND otherLaplace>" << gCCDParams.m_nMinLaplaceOnOther;
			szEventDesc << "\nnewLaplace=" << m_PixelAnalResults.laplaceSum 
                     << ", otherLaplaceSum=" << m_PixelAnalResults.prevLaplaceSum << "(at " << m_PixelAnalResults.prevLaplaceSumX << "," << m_PixelAnalResults.prevLaplaceSumY << ")"
							<< ", minLaplace=" << m_PixelAnalResults.prevLaplaceSumMin << "\n";
		}
		if(gCCDParams.m_nMaxOfAverageOfPrevN){
			szEventDesc << "Average of previous " << gCCDParams.m_nMaxOfAverageOfPrevN 
							<< " frames avPrev=" << m_PixelAnalResults.maxAverageOfPrev << "\n";
		}
	}else{		
		if(gCCDParams.m_eBrigtheningCheckType==ePercentIncrease){
			char buf[10];

         sprintf(buf,"%.2f",(double(m_PixelAnalResults.newSum)-double(m_PixelAnalResults.otherSum))/double(m_PixelAnalResults.otherSum));
			szEventDesc << "(newSum-otherSum)/otherSum = " << buf << " > " 
                     << gCCDParams.m_IncreaseTreshPercent;
		}else{
			szEventDesc << "(newSum-otherSum)>" << gCCDParams.m_IncreaseTreshADU;
		}
	}
	szEventDesc << "\n";	

	if(!m_bFirstLevelAccepted){
		szEventDesc << "Event was rejected by first level trigger\n";
	}else{
		// show cluster and sums - causing rejection :
		if(m_ClusterCount>0){
			szEventDesc << "Found cluster has " << m_ClusterCount << " points\n";
			szEventDesc << "Points in cluster :\n";
			int i=0;
			for(i=0;i<m_ClusterCount;i++){
				long x = ((m_Cluster)[i] % xSize);
				long y = ((m_Cluster)[i] / xSize);
				szEventDesc << "  (" << x << "," << y << ")" << ",";
				if((i+1)%6==0 && i)
					szEventDesc << "\n";
			}
			szEventDesc << "\n";			
		}
		if(gCCDParams.m_bConfirmReq){
			LONGLONG_T maxPrevSum = my_find_max_value_long( prevClusterSum, prevFramesCount );

			szEventDesc << "Event was rejected by 2 level trigger\n";
			szEventDesc << "Found cluster has " << m_ClusterCount << " points\n";
			szEventDesc << "Points in cluster :\n";
			int i=0;
			for(i=0;i<m_ClusterCount;i++){
				long x = ((m_Cluster)[i] % xSize);
				long y = ((m_Cluster)[i] / xSize);
				szEventDesc << "  (" << x << "," << y << ")" << ",";
				if((i+1)%6==0 && i)
					szEventDesc << "\n";
			}
			szEventDesc << "\n";
			for(i=0;i<=prevFramesCount;i++){
				szEventDesc << "Sum in cluster on frame " << (i+1) << " sum=" 
   	   	            << prevClusterSum[i] << "\n";
			}
			szEventDesc << "NEW Cluster Sum = " << newClusterSum << "\n";

			if( m_PixelAnalResults.eventType != eBrighten){
				szEventDesc << "Condition : PREV<" << (m_ClusterCount*(gCCDParams.m_ConfMaxPrevPerPixel)) 
   		   	         << " AND NEW>" << (m_ClusterCount*(gCCDParams.m_ConfTresholdPerPixel)) << "\n";
			}else{
				if(gCCDParams.m_eBrigtheningCheckType==ePercentIncrease){
					char buf[10];
					sprintf(buf,"%.2f",(double(newClusterSum)-double(maxPrevSum))/double(maxPrevSum));
					szEventDesc << "(newSum-otherSum)/otherSum = " << buf << " > " 
      	               << gCCDParams.m_IncreaseTreshPercent << "\n";
				}else{
					szEventDesc << "(newSum - otherSum)=(" 
									<< newClusterSum << " - " << maxPrevSum
                           << ")>" << gCCDParams.m_IncreaseTreshADU << "\n";
				}
			}
			szEventDesc << "\n";
		}else{
			szEventDesc << "No second level check required\n";
		}
	}
	if(gCCDParams.m_ConfirmEventsOnNextNFrames>0 && m_bCheckOnNext){
		szEventDesc << "MAX values on next frames :";
		for(int i=0;i<gCCDParams.m_ConfirmEventsOnNextNFrames && i<MAX_NEXT_FRAMES;i++){
			szEventDesc << m_PixelAnalResults.LaplaceOnNextFrames[i] << ",";
		}
		szEventDesc << "\n";
		if(m_PixelAnalResults.m_bRejectedByNextFrames){
			szEventDesc << "Event was rejected by next frames verification\n";
		}else{
			szEventDesc << "Event was accepted by next frames verification\n";
		}
	}

	if(m_nAddInfo){
		szEventDesc << "Additional event characteristics :\n";
		for(int a=0;a<m_nAddInfo;a++){
			const char* desc = CccdReport::GetAddInfoDesc( (CccdReport::eAddInfoType_T)a );
			szEventDesc << desc << " : " << m_AdditionalInfo[a] << "\n";
		}
	}

	if(m_PixelAnalResults.m_bLocalMaxRejected){
		szEventDesc << "Event rejected by Local Max Requirement\n";
	}

	if(m_PixelAnalResults.m_bOutSideAcceptanceRegion){
		szEventDesc << "Event rejected - outside acceptance region\n";
	}

	if(m_PixelAnalResults.m_nCountAbove>0){
		szEventDesc << "Treshold=" << m_PixelAnalResults.treshold_for_not_more_then_n_exceeds
						<< " exceeded on " << m_PixelAnalResults.m_nCountAbove << " frames ";
		if(m_PixelAnalResults.m_bNotMoreThenNAboveRejected)
			szEventDesc << " - rejected\n";
		else
			szEventDesc << " - accepted\n";
	}

	if(gCCDParams.m_CheckEventShape>=0){
		szEventDesc << "Event sphericity=" << m_PixelAnalResults.m_Sphericity << " required >=" << gCCDParams.m_CheckEventShape << "\n";
		szEventDesc << "Cluster center of mass = (" << m_PixelAnalResults.x0 << "," 
						<< m_PixelAnalResults.y0 << ") max_redial=" << m_PixelAnalResults.max_redial << "\n";
		szEventDesc << "Precisly center = (" << mystring::double2string( m_PixelAnalResults.x0_real )
						<< "," << mystring::double2string( m_PixelAnalResults.y0_real ) << ")\n";
		szEventDesc << "Cluster treshold = " << mystring::double2string(m_PixelAnalResults.max_noise_level) << "\n";
		if(m_PixelAnalResults.m_bRejectedByEventShape){
			szEventDesc << " - REJECTED\n";
		}else{
			szEventDesc << " - ACCEPTED\n";
		}
	}

	if(gCCDParams.m_nMaxOfAverageOfPrevN>0){
		szEventDesc << "Position on previous frames :";
		for(register int i=0;i<(gCCDParams.m_nMaxOfAverageOfPrevN+1);i++){
			szEventDesc << "(" << m_PixelAnalResults.PrevFramesX[i] << "," << m_PixelAnalResults.PrevFramesY[i] << ") , ";
		}
		szEventDesc << "\n";
	}
	
}

void CccdReport::WriteEventDetails( CCDPipeline* pPipeline, int eventNo )
{
	mystring szEventDesc(2048,2048);
	char szTmp[256];

	mystring szEventDT;
   if( m_Time>0 ){
      szEventDT = get_gmtime_string( m_Time );
   }
	
	szEventDesc << "Frame#      :\t" << m_FrameIndex <<"\n";
	szEventDesc << "Day Frame#  :\t" << m_DayFrameIndex <<"\n";

	mystring szNightDate;
	get_night_date_local( m_Time, szNightDate );
	szEventDesc << "Night Date  :\t" << szNightDate << "\n";
	szEventDesc << "Event#      :\t" << eventNo << "\n";
	szEventDesc << "EventType   :\t" << szEventTypeDescTab[(LONG_T)m_PixelAnalResults.eventType] << "\n";
	szEventDesc << "(x,y)       :\t" << m_MaxPoint.x << "," << m_MaxPoint.y << "\n";

	mystring szRA = AstroAngle::toString( m_AstroCoord.RA, ANGLE_RA_TYPE ).c_str();
	szEventDesc << "(RA,DEC)    :\t" << szRA << "," << AstroAngle::rad2deg( m_AstroCoord.Dec ) << "\n";

	szEventDesc << "UT_DATE_TIME:\t" << szEventDT << "\n";   
	szEventDesc << "GEN         :\t" << m_bGenerated << "\n";
	szEventDesc << "Lap_new     :\t" << m_PixelAnalResults.laplaceSum << "\n";
	szEventDesc << "Lap_prev_av :\t" << m_PixelAnalResults.maxAverageOfPrev << "\n";
	// szEventDesc << "Lap_next    :\t" << m_PixelAnalResults.laplaceOnNext << "\n";

	sprintf(szTmp,"%.6f",(float)m_PixelAnalResults.m_Sphericity);
	szEventDesc << "Sphericity  :\t" << szTmp << "\n";

	sprintf(szTmp,"%.6f",(float)m_PixelAnalResults.m_fBlackRatio);
	szEventDesc << "Black ratio  :\t" << szTmp << "\n";

	sprintf(szTmp,"%.6f",(float)m_PixelAnalResults.m_PrevLapInPixel);
	szEventDesc << "Aver in pixel:\t" << szTmp << "\n";

	mystring szClusterList;
	register int i=0;
	register int xSize = pPipeline->GetXSize();	
	for(i=0;i<m_ClusterCount;i++){
		szClusterList << "(" << (m_Cluster[i]%xSize) << "," << (m_Cluster[i]/xSize) << "),";
		if((i+1)%5==0 && i>0){
			szClusterList << "\t\t\t";
		}
	}	
	szEventDesc << "Cluster     :\t" << szClusterList << "\n";

	szEventDesc << "Prev-laplace:\t";
	mystring szPrevValues,szPrevPos;
	for(i=0;i<m_PixelAnalResults.PrevLaplaceValuesCount;i++){
		szPrevValues << m_PixelAnalResults.PrevLaplaceValues[i] << ",";
		szPrevPos << "(" << m_PixelAnalResults.PrevFramesX[i] << "," << m_PixelAnalResults.PrevFramesY[i] << "),";
	}
	szEventDesc << szPrevValues << "\n";
	szEventDesc << "Prev Pos    :\t" << szPrevPos << "\n\n";

	szEventDesc << "PLUS VALUES :";
	for(i=0;i<m_PixelAnalResults.laplacePlusCount;i++){
		szEventDesc << m_PixelAnalResults.laplacePlusValues[i] << ",";
	}	

	szEventDesc << "\nMINUS VALUES :";
	for(i=0;i<m_PixelAnalResults.laplaceMinusCount;i++){
		szEventDesc << m_PixelAnalResults.laplaceMinusValues[i] << ",";
	}
	szEventDesc << "\n";
	

	char szSubDir[128];
   sprintf(szSubDir,"Frame%.5d",m_DayFrameIndex);
	
	mystring szFName;
	szFName << gCCDParams.GetOutputDir() << "/Events/" <<  szSubDir << "/";
	if( pPipeline->GetPipelineCount()>1 )
      szFName << "Cam" << pPipeline->GetPipelineIndex() << "/";
	szFName << "event_details.txt";
	MyOFile out( szFName.c_str() , "a+" );
	out.Printf("%s",szEventDesc.c_str());	
}

void CccdReport::GetEventTypeDesc( char* szOutType )
{
	if( strlen( m_szGRBInfo.c_str() ) ){
		sprintf( szOutType, "%c(%s)",m_EventType,m_szGRBInfo.c_str());
		return;
	}

	if( m_EventType==EVENT_TYPE_SAT ){
		char szTmp[16];
		strncpy(szTmp,m_szSatName.c_str(),SAT_NAME_LEN);
		szTmp[6]='\0';
		sprintf( szOutType, "%c(%s:%c:%.2f:%c:%.2f)", m_EventType, szTmp, 
			m_SatType, m_SatCoicRadiusInSec, visibility, m_SatInfo->sat_alt );
	}else{
		if( m_EventType==EVENT_TYPE_STAR ){
			sprintf( szOutType, "%c(%s)", m_EventType, m_szSatName.c_str());
		}
	}

	if( m_EventType==EVENT_TYPE_FLASH ){
		if( strlen( m_szSatName.c_str() ) ){
			sprintf( szOutType, "%c(%s:%c:%.2f:%c)",m_EventType, m_szSatName.c_str(),
								m_SatType, m_SatCoicRadiusInSec, visibility  );
		}else{
			sprintf( szOutType, "%c", m_EventType);
		}
	}
}


void CccdReport::ParseTrackDesc( const char* szTrack  )
{
	if( szTrack && szTrack[0] ){
		MyParser pars = szTrack;
		
		const char* szTmp = pars.GetItemUpTo( "-" );
		m_PixelAnalResults.bestTrackID = atof( szTmp );
		szTmp = pars.GetItemUpTo( "-" );
		m_PixelAnalResults.minChi2 = atof( szTmp );
		szTmp = pars.GetItemUpTo( "- " );
      m_PixelAnalResults.minDist = atof( szTmp );		
	}
}


void CccdReport::ParseRA_DEC( const char* szRA_DEC )
{
	if( szRA_DEC ){
		MyParser pars = szRA_DEC+1;
		const char* szTmp;

		szTmp = pars.GetItemUpTo(",");
		m_AstroCoord.RA =  AstroAngle::deg2rad( AstroAngle::ra2deg( szTmp ) );
		szTmp = pars.GetItemUpTo(",)");
		m_AstroCoord.Dec = AstroAngle::deg2rad( atof(szTmp ) );
	}
}


void CccdReport::ParseEventType( const char* szEventTypeDesc )
{
	if( szEventTypeDesc && szEventTypeDesc[0] ){
		char szTmp[64];

		m_EventType = szEventTypeDesc[0];					
		if( szEventTypeDesc[0]!=EVENT_TYPE_FLASH ){
			if( szEventTypeDesc[0]==EVENT_TYPE_SAT ){
				// SAT_NAME_LEN : 
				MyParser pars = szEventTypeDesc+2;
				// sscanf( szEventTypeDesc,"%c(%6s:%c:%lf:%c)",&m_EventType,szTmp,&m_SatType,&m_SatCoicRadiusInSec,&visibility);				
				// m_szSatName = szTmp;
				const char* szTmp;
				szTmp = pars.GetItemUpTo( ":" );
				if( szTmp )
					m_szSatName = szTmp;
				szTmp = pars.GetItemUpTo( ":" );
				if( szTmp )
					m_SatType = szTmp[0];
				szTmp = pars.GetItemUpTo( ":" );
				if( szTmp )
					m_SatCoicRadiusInSec = atof( szTmp );
				szTmp = pars.GetItemUpTo( ":" );
            if( szTmp )
               visibility = szTmp[0];
				szTmp = pars.GetItemUpTo( ":)" );
            if( szTmp )
               m_SatInfo->sat_alt = atof( szTmp );
			}
			if( szEventTypeDesc[0]==EVENT_TYPE_STAR ){
				// sscanf( szEventTypeDesc,"%c(%s)",&m_EventType,szTmp);
				// m_szSatName = szTmp;
				MyParser pars = szEventTypeDesc+2;
				const char* szTmp;
            szTmp = pars.GetItemUpTo( ")" );
            if( szTmp )
               m_szSatName = szTmp;
			}
		}else{
			// flash - but can contain also closest sat info :
			if( strlen(szEventTypeDesc)>1 ){
				// sscanf( szEventTypeDesc,"%c(%6s:%c:%lf:%c)",&m_EventType,szTmp,
				//	&m_SatType, &m_SatCoicRadiusInSec, &visibility);
				// m_szSatName = szTmp;				
				MyParser pars = szEventTypeDesc+2;
				const char* szTmp;
            szTmp = pars.GetItemUpTo( ":" );
            if( szTmp )
               m_szSatName = szTmp;
				szTmp = pars.GetItemUpTo( ":" );
				if( szTmp )
					m_SatType = szTmp[0];
				szTmp = pars.GetItemUpTo( ":" );
            if( szTmp )
					m_SatCoicRadiusInSec	= atof( szTmp );
				szTmp = pars.GetItemUpTo( ":)" );
            if( szTmp )
               visibility = szTmp[0];
			}
		}
	}
}

// be carefull - you must allocate enough for szLine buffer !!!
long CccdReport::GetEventDescStr( char* szLine )
{
	if(strlen( m_szExternal.c_str())==0){
		m_szExternal = "-";
	}

	mystring szEventDT;
	if( m_Time>0 ){
		szEventDT = get_gmtime_string( m_Time );
	}else{
		szEventDT = get_gmtime_string( get_dttm() );
	}

	mystring szMag = m_Magnitude;
	if(strlen(szMag.c_str())==0){
		szMag = "0.00";
	}

	char szEventType[128];
	GetEventTypeDesc( szEventType );
	double coic_radius_sec = AstroAngle::rad2arcsec( m_CoicRadiusInRad );


	long ret =sprintf(szLine,m_EventDescOutFmt.c_str(),
           (LONG_T)m_DayFrameIndex, (LONG_T)m_CameraIdx, m_MaxPoint.x, m_MaxPoint.y,
           m_bGenerated, m_bIdentified, m_PixelAnalResults.m_bRejectedByNextFrames,
			  m_PixelAnalResults.m_bRejectedDueToTrack,
			  m_LowLeft.x, m_LowLeft.y, m_TopRight.x, m_TopRight.y, szMag.c_str(),
			  szEventTypeDescTab[(LONG_T)m_PixelAnalResults.eventType],szEventDT.c_str(),
			  m_PixelAnalResults.m_bInTriggerMode,m_PixelAnalResults.m_bRejectedByCoic,
			  m_PixelAnalResults.laplaceSum,m_PixelAnalResults.maxAverageOfPrev,
			  AstroAngle::toString( m_AstroCoord.RA, ANGLE_RA_TYPE ).c_str(), AstroAngle::rad2deg( m_AstroCoord.Dec ),
// 			  m_AstroCoord.RA, m_AstroCoord.Dec,
				m_PixelAnalResults.m_MaxClusterValue,m_PixelAnalResults.m_Significance,
				m_PointTransformed.x, m_PointTransformed.y, m_CoicRadius, 
				szEventType, m_FrameIndex, coic_radius_sec,
				m_PixelAnalResults.m_Sphericity, m_PixelAnalResults.m_bRejByTrackOnSingleCam,
				m_PixelAnalResults.bestTrackID,m_PixelAnalResults.minChi2,
				m_PixelAnalResults.minDist, m_PixelAnalResults.veloCheckOK,
				m_PixelAnalResults.rx, m_PixelAnalResults.ry, EvtIdx,
				m_nFrameStarCount, m_szExternal.c_str(),
				m_PixelAnalResults.PixelRawValue,
				m_PixelAnalResults.MaxPixelRawValue, m_MinDistOutCone );

	return ret;
}

BOOL_T CccdReport::ParseOutputLine( const char* szLine )
{
	m_Cluster = NULL;
	m_ClusterCount = 0;
	char szRA[64],szMag[20],szEventType[20],szEventDT[20];
	char szEventTypeDesc[64],szRA_DEC[128],szExt[128];
	double dec;
	int h,m;
	double s;

	int index;
	szExt[0]='\0';
	int ret = sscanf( szLine, m_EventDescInFmt.c_str(), &index, &m_CameraIdx,
							&(m_Point.x), &(m_Point.y), &m_bGenerated, &m_bIdentified,
							&(m_PixelAnalResults.m_bRejectedByNextFrames), &(m_PixelAnalResults.m_bRejectedDueToTrack),
							&(m_LowLeft.x), &(m_LowLeft.y), &(m_TopRight.x), &(m_TopRight.y), szMag,  
							szEventType, 
							szEventDT, 
							&(m_PixelAnalResults.m_bInTriggerMode), 
							&(m_PixelAnalResults.m_bRejectedByCoic),
							&(m_PixelAnalResults.laplaceSum),
							&(m_PixelAnalResults.maxAverageOfPrev),
							szRA_DEC, &(m_PixelAnalResults.m_MaxClusterValue),
							&(m_PixelAnalResults.m_Significance),
					      &(m_PointTransformed.x), &(m_PointTransformed.y), 
							&m_CoicRadius, szEventTypeDesc, &(m_FrameIndex),
							&m_CoicRadiusInRad, // converted in lines below 
							&(m_PixelAnalResults.m_Sphericity),
                     &(m_PixelAnalResults.m_bRejByTrackOnSingleCam),
							&(m_PixelAnalResults.bestTrackID),
							&(m_PixelAnalResults.minChi2),
							&(m_PixelAnalResults.minDist),
							&(m_PixelAnalResults.veloCheckOK),
							&(m_PixelAnalResults.rx), &(m_PixelAnalResults.ry),
							&EvtIdx, &m_nFrameStarCount, szExt,
							&(m_PixelAnalResults.PixelRawValue),
							&(m_PixelAnalResults.MaxPixelRawValue),
							&(m_MinDistOutCone) );

	if( ret<40 ){
		printf("ERROR !\n");
		printf("ONLY %d items found !\n",ret);
		printf("INCORRECT FORMAT OF LOG FILE LINE : \n%s\n",szLine);
		if( !CccdReport::m_bIgnoreReadError ){
			return FALSE;
		}
	}

	m_szExternal = szExt;
	if(strlen(m_szExternal.c_str())==0)
		m_szExternal = "-";
	m_CoicRadiusInRad = AstroAngle::arcsec2rad( m_CoicRadiusInRad );
	m_Magnitude = szMag;
	m_MaxPoint = m_Point;
	m_DayFrameIndex = index;
	m_Time = get_gmtime_from_string( szEventDT );

	// m_EventType = szEventTypeDesc[0];
	ParseEventType( szEventTypeDesc );
		
	// m_AstroCoord.Dec = AstroAngle::deg2rad( dec );
	// m_AstroCoord.RA = mysign(h)*AstroAngle::timeangle2rad( abs(h), m , s );	
	ParseRA_DEC( szRA_DEC );
	
	mystring szTmp = get_gmtime_string( m_Time );

	m_PixelAnalResults.eventType = GetEventType( szEventType );
	return (ret>0);							
}




void CccdReport::GetEventDesc( mystring& desc )
{
	char szLine[1000];
	
	GetEventDescStr( szLine );
	desc = szLine;

/*	desc = "";
	desc << m_FrameIndex << " " << m_CameraIdx << " " << m_Point.x << " " << m_Point.y << " " 
        << m_bGenerated << " " << m_bIdentified << " " 
		  << m_PixelAnalResults.m_bRejectedByNextFrames << " "
        << "(" << m_LowLeft.x << "," << m_LowLeft.y << ")-("
        << m_TopRight.x << "," << m_TopRight.y << ") " << m_Magnitude << "\n";*/
}

int CCDEventList::GetNotIndentifiedCount()
{
	int ret=0;
	for(int i=0;i<size();i++){
		if(!((*this)[i].m_bIdentified))
			ret++;
	}
	return ret;
}

mystring CCDEventList::GetOutputDir()
{
	return gCCDParams.GetOutputDir();
}

mystring CCDEventList::GetEventFileName( CccdReport& event, int FrameIndex, int EventNo )
{
	mystring szRet;
	szRet << "eventframe" << event.m_FrameIndex 
			<< "_camera" << event.m_CameraIdx << "_eventno" << EventNo
			<< "_currframe" << FrameIndex;
	return szRet;
}

CccdReport* CCDEventList::Find( int frameno , int x, int y )
{
	CCDEventList::iterator i;
	for(i=begin();i!=end();i++){
		if( i->m_DayFrameIndex == frameno && ((int)i->m_MaxPoint.x == x)  && ((int)i->m_MaxPoint.y == y) ){
			return &(*i);
		}else{
			// assuming sort by frames 
			if( i->m_DayFrameIndex > frameno )
				return NULL;
		}
	}
	return NULL;
}

CccdReport* CCDEventList::FindFromBack( int frameno , int x, int y )
{
	int s=size();
	for(int i=(s-1);i>=0;i--){
		CccdReport& evt = (*this)[i];

		if( evt.m_DayFrameIndex == frameno && ((int)evt.m_MaxPoint.x == x)  && ((int)evt.m_MaxPoint.y == y) ){
			return &((*this)[i]);
		}else{
			// assuming sort by frames 
			if( evt.m_DayFrameIndex < frameno )
				return NULL;
		}
	}
	return NULL;
}

BOOL_T CCDEventList::CheckSortByFrame()
{
	int prev_frame=-1000;
	CCDEventList::iterator i;

   for(i=begin();i!=end();i++){
		if( prev_frame > i->m_DayFrameIndex ){
			printf("prev_frame=%d > current_frame=%d\n",prev_frame,i->m_DayFrameIndex);
			return FALSE;
		}
	}
	return TRUE;
}

CccdReport* CCDEventList::FindEvent( int x, int y, int frameno, int redial )
{
	CCDEventList::iterator i;
	for(i=begin();i!=end();i++){
		if( i->m_DayFrameIndex == frameno ){
			double r = sqrt( ((i->m_Point).x-x)*((i->m_Point).x-x) + ((i->m_Point).y-y)*((i->m_Point).y-y) );
			if( r<=redial ){
				return &(*i);
			}
		}
	}
	return NULL;

}

CccdReport* CCDEventList::FindEvent( int x, int y, int redial )
{
	CCDEventList::iterator i;
	for(i=begin();i!=end();i++){
		double r = sqrt( ((i->m_Point).x-x)*((i->m_Point).x-x) + ((i->m_Point).y-y)*((i->m_Point).y-y) );
		if( r<=redial ){
			return &(*i);
		}
	}
	return NULL;
}

CccdReport* CCDEventList::FindEventByMaxPoint( int x, int y, int redial )
{
	CCDEventList::iterator i;
	for(i=begin();i!=end();i++){
		double r = sqrt( ((i->m_MaxPoint).x-x)*((i->m_MaxPoint).x-x) + ((i->m_MaxPoint).y-y)*((i->m_MaxPoint).y-y) );
		if( r<=redial ){
			return &(*i);
		}
	}
	return NULL;
}


CccdReport* CCDEventList::FindEventByRaDec( double ra, double dec, double radius )
{
	CCDEventList::iterator i;
	for(i=begin();i!=end();i++){
		double evt_ra = (i->m_AstroCoord).RA;
      double evt_dec = (i->m_AstroCoord).Dec;
      double r = sqrt( (evt_ra-ra)*(evt_ra-ra) + (evt_dec-dec)*(evt_dec-dec) );
//		printf("r = %.8f rad , r = %.8f sec\n",r,AstroAngle::rad2arcsec( r ) );

		if( r <= radius ){
			return &(*i);
		}
	}
	return NULL;	
}

void CCDEventList::Dump( CccdReport& evt, eLogFileFormat _format, const char* szFileName )
{
	char tmp[1028];
	if( _format == eFinal ){
		SprintfEvent( tmp , evt );
	}else{
		evt.GetEventDescStr( tmp );
	}
	if( gCCDParams.m_bDumpAllLogToStdout ){
		_TRACE_PRINTF_0("%s",tmp);
	}
	if( szFileName && szFileName[0] ){
		BOOL_T bHeader=FALSE;
	   if( MyFile::GetFileSize( szFileName )<=0 ){
			bHeader=TRUE;
		}

		MyOFile out( szFileName, "a+" );
		
		if( bHeader ){
			if( _format == eFinal ){
				out.Printf("%s",CccdReport::m_FinalEventHeader.c_str());
			}else{
				out.Printf("%s",CccdReport::m_EventDescHeader.c_str());
			}
		}
      out.Printf("%s",tmp);
	}
}

void CCDEventList::Dump( const char* szFileName ){
	if( size()>0 ){
		if( format == eFinal ){
			_TRACE_PRINTF_0("%s",CccdReport::m_FinalEventHeader.c_str());
		}else{
			mystring szHeader;
			CccdReport::GetDescHeader(szHeader);
			_TRACE_PRINTF_0("%s",szHeader.c_str());
			/*if( szFileName && szFileName[0] ){
				if( !MyFile::DoesFileExist( szFileName )){
					MyOFile out( szFileName );
					out.Printf("%s",szHeader.c_str());
				}
			}*/
		}
		for(int i=0;i<size();i++){
			Dump( (*this)[i], format, szFileName );
		}
		if( !gCCDParams.m_bDumpAllLogToStdout ){
			printf("EVENTS# = %d ( dumped to file %s )\n",size(),szFileName);
		}
	}
}

int CCDEventList::Read( const char* fname, int frameno )
{
	clear();
	if( MyFile::DoesFileExist( fname )){
		// read header line :
		_TRACE_PRINTF_5("checking format of file : %s ...",fname);
		MyIFile infile;
		const char* pLine;
		mystring szHeader;
		if(infile.Open(fname,"r", FALSE)){
			if( pLine = infile.GetLine()	){
				szHeader = pLine;
			}
		}
		infile.Close();
		
		if( strlen( szHeader.c_str() ) ){
			format=eUnknownLogFormat;
			int len_to_check = strlen( szHeader.c_str() )-2; // just to avoid newline character 

			if( strncmp( szHeader.c_str() , CccdReport::m_EventDescHeader.c_str(), len_to_check )==0 ){
				_TRACE_PRINTF_5("verif\n");
				format = eVerif;
			}
			if( strncmp( szHeader.c_str() , CccdReport::m_FinalEventHeader.c_str(), len_to_check )==0 ){
            _TRACE_PRINTF_5("final\n");
				format = eFinal;
         }

			if( format==eUnknownLogFormat ){
				printf("unknown format of log file , using default format\n");
				printf("but this may cause a problem ...\n");
				format = m_DefaultFormat;
			}

			if( format!=eUnknownLogFormat ){
				if( format==eVerif ){
					return ReadEvents( fname, TRUE, -1, frameno );
				}
				if( format==eFinal ){
					return ReadFinalEvents( fname, (*this), 0, 0, frameno );
				}
			}			
		}
	}
	return size();
}

int CCDEventList::ReadEvents( const char* fname, BOOL_T bExcp, int min_frame,
										int frameno )
{
	clear();
	MyIFile infile;
	BOOL_T bAllZero=TRUE;

	if(infile.Open(fname,"r", bExcp)){
		const char* pLine;
		while( pLine = infile.GetLine() ){
			if(strlen(pLine)==0 || pLine[0]=='#')
				continue;
			CccdReport tmp;
			tmp.ParseOutputLine( pLine );
			if( strlen( tmp.m_EvtID.c_str() )>0 ){
				bAllZero=FALSE;
			}
			if( min_frame<0 || tmp.m_DayFrameIndex>=min_frame ){
				if( frameno<=0 || frameno==tmp.m_DayFrameIndex ){
					push_back( tmp );
				}
			}
		}
	}

	if( bAllZero ){
		SetEvtNo();
	}

	return size();
}

int CCDEventList::SprintfEvent( char* ptr, CccdReport& evt )
{
	char szEventType[128];
	mystring szRA = AstroAngle::toString( (evt.m_AstroCoord).RA, ANGLE_RA_TYPE ).c_str();
	evt.GetEventTypeDesc( szEventType );
	int zero_value = 0;
	double coic_radius_sec = AstroAngle::rad2arcsec( evt.m_CoicRadiusInRad );

	mystring szEventDT;
   if( evt.m_Time>0 ){
      szEventDT = get_gmtime_string( evt.m_Time );
   }

	int len = sprintf( ptr, CccdReport::m_FinalEventOutFmt.c_str(), 
							   evt.m_EvtID.c_str(),
								(int)evt.m_DayFrameIndex,
							  (int)evt.m_FrameIndex,
							  (int)evt.EvtIdx, (int)(evt.m_MaxPoint).x, 
							(int)(evt.m_MaxPoint).y,
							  (evt.m_PixelAnalResults).laplaceSum,
							  (evt.m_PixelAnalResults).maxAverageOfPrev,	 
							  (evt.m_PixelAnalResults).m_MaxClusterValue,
							  (evt.m_PixelAnalResults).PixelRawValue , 
								(evt.m_PixelAnalResults).m_Sphericity,
							  (evt.m_PixelAnalResults).m_Significance,
							  evt.m_ClusterCount, evt.m_CoicRadius, 
							  (evt.m_PixelAnalResults).m_PrevLapInPixel,
							  (evt.m_PixelAnalResults).m_fBlackRatio,
							  szRA.c_str(), AstroAngle::rad2deg( (evt.m_AstroCoord).Dec ),
							  szEventType, coic_radius_sec,
							  (evt.m_AstroCoord).RA,(evt.m_AstroCoord).Dec,
							  evt.m_Time, (evt.m_PixelAnalResults).laplaceOnNext,
							  evt.m_PixelAnalResults.bestTrackID,evt.m_PixelAnalResults.minChi2,
							  evt.m_PixelAnalResults.minDist,
							  evt.m_PixelAnalResults.veloCheckOK, evt.m_PixelAnalResults.rx,
							  evt.m_PixelAnalResults.ry, szEventDT.c_str(),
							  evt.m_nFrameStarCount, evt.m_szExternal.c_str(),
							  evt.m_MinDistOutCone );
	return len;
}


int CCDEventList::SprintfGrb( char* ptr, CGRBInfo* pGrbInfo )
{
	double ra_in_rad = AstroAngle::deg2rad( pGrbInfo->ra );
	mystring szRA = AstroAngle::toString( ra_in_rad , ANGLE_RA_TYPE ).c_str();

	int ret = sprintf(ptr,CccdReport::m_GrbInfoFmt.c_str(),
			  pGrbInfo->grb_id,pGrbInfo->gcn_id,pGrbInfo->source_id,
			  pGrbInfo->foundgrbname,szRA.c_str(),pGrbInfo->dec,pGrbInfo->grbTime );
	return ret;	
}


int CCDEventList::DumpFinalEvents( const char* fname , CCDEventList& ccdReport,
											  BOOL_T bUseOutDir /* =TRUE */)
{
	int ret=0;	
	if( ccdReport.size() ){
		mystring szOutName;
		const char* szOutDir = gCCDParams.GetOutputDir();
		if( szOutDir && szOutDir[0] && bUseOutDir)
			szOutName << gCCDParams.GetOutputDir() << "/";
		szOutName << fname;
		szOutName.env2str();
		BOOL_T bHeader=FALSE;
		int file_size = MyFile::GetFileSize( szOutName.c_str() );
		if( file_size <= 0 ){
			bHeader = TRUE;
		}
		MyOFile out( szOutName.c_str(), "a" );
		if(bHeader){
			out.Printf("%s",CccdReport::m_FinalEventHeader.c_str());
		}

		int size=0;
		int string_size = 1000*(ccdReport.size()+2);
      char* pEventsStr = new char[ string_size ];
		char* ptr = pEventsStr;
		
		CCDEventList::iterator i;
		char szEventType[64];
		for(i=ccdReport.begin();i!=ccdReport.end();i++){
			if(size<string_size){
				int len = SprintfEvent( ptr, *i );
				ptr += len;
				size += len;
			}else{
				mystring szErr;
				szErr << "ERROR : could not write all events !!! buffer to small for : " 
						<< (int)ccdReport.size() << " events";
				printf("%s",szErr.c_str());
				PrintError( NULL, i->m_FrameIndex, EVENT_BUFFER_TO_SMALL, szErr.c_str() );
				break;
			}
		}
		out.Printf("%s",pEventsStr);
		out.Close();
		
		delete [] pEventsStr;		
	}
	return ccdReport.size();
}


int CCDEventList::DumpFinalEventLine( const char* fname, CccdReport& evt, CGRBInfo* pGrbInfo )
{
	mystring szOutName=fname;
	szOutName.env2str();
	BOOL_T bHeader=FALSE;
	if( MyFile::GetFileSize( szOutName.c_str() )<=0 ){
		bHeader = TRUE;
	}
	MyOFile out( szOutName.c_str(), "a" );
	if(bHeader){
		mystring szHeader=CccdReport::m_FinalEventHeader;
		szHeader << CccdReport::m_GrbInfoHeader;
		out.Printf("%s",szHeader.c_str());
	}

	int size=0;
	char szEvt[4096],szGRB[4096];
	szGRB[0]='\0';
		
	SprintfEvent( szEvt, evt );
	if( pGrbInfo ){
		SprintfGrb( szGRB, pGrbInfo );	
	}

	mystring szLine;
	szLine << szEvt;
	if( pGrbInfo ){
		szLine << " " << szGRB;
	}
	out.Printf("%s\n",szLine.c_str());
	out.Close();	
	return 1;
}

int CCDEventList::ReadFinalEvents( const char* fname , CCDEventList& ccdReport,
											  time_t startTime, time_t endTime,
											  int frameno )
{
	ccdReport.clear();
	if(!MyFile::DoesFileExist( fname )){
		return 0;
	}
	char szRA[64],szEventType[128],szID[64],szEvtDTM[64],szTrackDesc[128],szExt[128];
	double dec,coic_radius_sec;
	double ra_rad,dec_rad;
	int zero_value,h,m;
	double s;
	BOOL_T bAllZero=TRUE;

	MyIFile in( fname );
	const char* pLine;
	while( pLine = in.GetLine() ){
		if( mystring::get_first_non_white( pLine )=='#' )
			continue;

		CccdReport evt;
		szEvtDTM[0]='\0';
		szExt[0]='\0';
		sscanf( pLine, CccdReport::m_FinalEventInFmt.c_str(), 
				  szID,
				  &(evt.m_DayFrameIndex),
				  &(evt.m_FrameIndex), 
				  &(evt.EvtIdx), &((evt.m_MaxPoint).x), 
				  &((evt.m_MaxPoint).y),
				  &((evt.m_PixelAnalResults).laplaceSum),
				  &((evt.m_PixelAnalResults).maxAverageOfPrev),	 
				  &((evt.m_PixelAnalResults).m_MaxClusterValue),
				  &zero_value, &((evt.m_PixelAnalResults).m_Sphericity),
				  &((evt.m_PixelAnalResults).m_Significance),
				  &(evt.m_ClusterCount), &(evt.m_CoicRadius), 
				  &((evt.m_PixelAnalResults).m_PrevLapInPixel),
				  &((evt.m_PixelAnalResults).m_fBlackRatio),
				  &h,&m,&s, &dec,
				  szEventType, &coic_radius_sec,
				  &(evt.m_AstroCoord.RA),&(evt.m_AstroCoord.Dec),
				  &(evt.m_Time), &((evt.m_PixelAnalResults).laplaceOnNext),
				  szTrackDesc,
				  &(evt.m_PixelAnalResults.veloCheckOK),&(evt.m_PixelAnalResults.rx),
				  &(evt.m_PixelAnalResults.ry), szEvtDTM,
				  &(evt.m_nFrameStarCount), szExt, &(evt.m_MinDistOutCone) );  

		evt.m_bIdentified = TRUE;
		evt.m_szExternal = szExt;
		if(strlen(evt.m_szExternal.c_str())==0)
	      evt.m_szExternal = "-";
		evt.m_EvtID = szID;
		evt.ParseEventType( szEventType );	
		evt.m_CoicRadiusInRad = AstroAngle::arcsec2rad( coic_radius_sec );
		evt.m_Point = evt.m_MaxPoint;
		(evt.m_PixelAnalResults).PixelRawValue = zero_value;

		evt.ParseTrackDesc( szTrackDesc );

		if( (startTime==0 && endTime==0) || 
			 (evt.m_Time>=startTime && evt.m_Time<=endTime) ){

			if( frameno<=0 || frameno==evt.m_DayFrameIndex){
				ccdReport.push_back( evt );
			}
		}
		if( strlen( evt.m_EvtID.c_str() )>0 ){
			bAllZero = FALSE;
		}
	}

	if( bAllZero ){
      ccdReport.SetEvtNo();
   }

	return ccdReport.size();
}

int CCDEventList::DumpEventReport( const char* fname , CCDEventList& ccdReport,
												LONG_T idx, BOOL_T bStdout, BOOL_T bToFile,
												BOOL_T bNew, BOOL_T bOnlyIdentified/*=TRUE*/,
												const char* szComment )
{
	// bToFile = (bToFile && gCCDParams.m_bDumpAllEvents);	
	// printf("bToFile=%d\n",bToFile);
	if( gCCDParams.m_bSameLogFiles ){
		// in case re-usage of log files enabled - ignore bNew
		bNew = FALSE;
	}

	mystring szMode = "w";
	if( !bNew )
		szMode = "a";
	Initialize();

	if(!bStdout && !bToFile)
		return ccdReport.size();

	if(ccdReport.size()){
		MyOFile fEvents;
		mystring fname_evt;
		if( gCCDParams.GetOutputDir() && strlen( gCCDParams.GetOutputDir() )>0 ){
			fname_evt << gCCDParams.GetOutputDir() << "/";
		}
		if(fname && strlen(fname))
			fname_evt << fname;
		else
			fname_evt << gCCDParams.m_szRunEventsLog;
		
		BOOL_T bNewFile=FALSE;
		if( MyFile::GetFileSize( fname_evt.c_str() )<=0 ){
			bNewFile = TRUE;
		}
			

		// printf("file = %s\n",fname_evt.c_str());
		if(bToFile)
			fEvents.Open( fname_evt.c_str(), szMode.c_str() );
		CCDEventList::iterator pEvt;
		mystring szHeader;
		long num=0;

		long string_size = 1000*(ccdReport.size()+10);
		char* pEventsStr = new char[ string_size ];
		// szEvents.AllocBuffer( (ccdReport.size()+10)*500 );
		if( bNewFile ){
			CccdReport::GetDescHeader(szHeader);
			if(bToFile)
				fEvents.Printf("%s",szHeader.c_str());
		}
		BOOL_T bIdentified=FALSE;

		char szLine[1000];
		PROFILER_START

		/*for(pEvt=ccdReport.begin();pEvt!=ccdReport.end();pEvt++){				  				
			if(pEvt->m_bIdentified || !bOnlyIdentified){
				pEvt->GetEventDescStr( szLine );
				szEvents << szLine;
			}
			if(pEvt->m_bIdentified){
				bIdentified=TRUE;
			}
		}*/

		long pos=0;
		long rep_count=ccdReport.size();
		for(register long e=0;e<rep_count;e++){				  				
			if(ccdReport[e].m_bIdentified || ccdReport[e].m_bGenerated){
				pos += ccdReport[e].GetEventDescStr( pEventsStr+pos );
				// szEvents << szLine;
			}
			if(ccdReport[e].m_bIdentified){
				bIdentified=TRUE;
			}
		}
		Assert(pos<=string_size,"CCDEventList::DumpEventReport writing outside buffer !!!");

		mystring szTrace;
		szTrace << "Loop over " << (int)ccdReport.size() << " events in CCDEventList::DumpEventReport took:";
		PROFILER_END(szTrace.c_str())

		if(bToFile){
			fEvents.Printf("%s",pEventsStr);
			fEvents.Close();
		}
		if(bStdout){
			if( gCCDParams.m_bDumpAllLogToStdout ){
				if(bIdentified){
					_TRACE_PRINTF_0("%s",szComment);
				}else{
					_TRACE_PRINTF_0("NO EVENTS\n");
				}
				if(ccdReport.size()){
					_TRACE_PRINTF_0("%s\n",szHeader.c_str());
					_TRACE_PRINTF_0("%s",pEventsStr);
				}
				_TRACE_PRINTF_0("##############################################\n");
			}else{
				printf("Dumping of events to stdout is disabled\n");
				printf("NEW EVENTS# = %d ( dumped to file %s )\n",rep_count,fname_evt.c_str());
			}
		}		
		delete [] pEventsStr;	
	}


	return ccdReport.size();
}


int CFrameEvents::GetMemSize( BOOL_T bShowAll )
{
	int ret = sizeof(CFrameEvents);
	
	/*if(bShowAll){
		printf("Number of tracks = %d\n",m_TracksOnFrame.size());
	}
	for(int i=0;i<m_TracksOnFrame.size();i++){
		int total = m_TracksOnFrame[i].GetMemSize( FALSE );
		if(bShowAll){
			printf("\tTrack %d uses %d kB :\n",i,total/1000);
			m_TracksOnFrame[i].GetMemSize( bShowAll );
		}
		ret += total;
	}*/
	return ret;
}


CFrameEvents::CFrameEvents( int FrameCounter )
: m_FrameCounter(FrameCounter)
{
}

CFrameEvents::~CFrameEvents()
{
	//	m_TracksOnFrame.clear();
}

int CFrameEvents::GetEventsCount( int cam_idx )
{
	Assert(cam_idx<size(),"No camera %d , max count = %d",cam_idx,size());
	return (*this)[cam_idx].size();
}

int CFrameEvents::HasFrameEvents( int cam_idx )
{
	Assert(cam_idx<size(),"No camera %d , max count = %d",cam_idx,size());
	return ( (*this)[cam_idx].size()>0 );
}


int CFrameEvents::GetNotMatchedEventsCountSkipGen( int cam_idx )
{
	Assert(cam_idx<size(),"No camera %d , max count = %d",cam_idx,size());

	int count=0;
	CCDEventList::iterator i;
	for(i=(*this)[cam_idx].begin();i!=(*this)[cam_idx].end();i++){
		if(!(i->m_PixelAnalResults).m_bRejectedDueToTrack && !i->m_bGenerated &&
			!(i->m_PixelAnalResults).m_bRejectedByNextFrames)
			count++;
	}
	return count;	
}

int CCDEventList::GetNotMatchedEventsCountSkipGen()
{
	int count=0;
	CCDEventList::iterator i;
	for(i=begin();i!=end();i++){
		if(!(i->m_PixelAnalResults).m_bRejectedDueToTrack && !i->m_bGenerated &&
			!(i->m_PixelAnalResults).m_bRejectedByNextFrames)
			count++;
	}
	return count;	
}


int CFrameEvents::GetNotMatchedEventsCount( int cam_idx )
{
	Assert(cam_idx<size(),"No camera %d , max count = %d",cam_idx,size());

	int count=0;
	CCDEventList::iterator i;
	for(i=(*this)[cam_idx].begin();i!=(*this)[cam_idx].end();i++){
		if(!i->m_PixelAnalResults.m_bRejectedDueToTrack)
			count++;
	}
	return count;
}

//CFrameEvents::CFrameEvents( const CFrameEvents& right )
//{
//	(*this) = right;
//}


//CFrameEvents& CFrameEvents::operator=( const CFrameEvents& right )
//{
//	m_FrameEvents = right.m_FrameEvents;
//	return (*this);
//}


// 

CMagStat::CMagStat()
: m_Magnitude(""),m_MagGenerated(0),m_MagIdentified(0)
{
}


CIdentStat::CIdentStat()
{
	m_pMagnitTab = new vector<CMagStat>();
	m_pParamValues = new vector<CEnvVar>();
	m_TotalGenerated = 0;
   m_TotalIdentified = 0;
   m_nBackground = 0;
}

CIdentStat::~CIdentStat()
{
	if(m_pMagnitTab)
		delete m_pMagnitTab;
	if(m_pParamValues) 
		delete m_pParamValues;
}

CIdentStat::CIdentStat( const CIdentStat& right )
{
	m_pMagnitTab = new vector<CMagStat>();
	m_pParamValues = new vector<CEnvVar>();
	(*this) = right;
}

CIdentStat& CIdentStat::operator=( const CIdentStat& right )
{
	m_TotalGenerated = right.m_TotalGenerated;
	m_TotalIdentified = right.m_TotalIdentified;	
	m_nBackground = right.m_nBackground;

	if(right.m_pMagnitTab){
		vector<CMagStat>::iterator i;
		for(i=(right.m_pMagnitTab)->begin();i!=(right.m_pMagnitTab)->end();i++){
			m_pMagnitTab->push_back( *i );
		}
	}
	if(right.m_pParamValues){
		vector<CEnvVar>::iterator i;
		for(i=(right.m_pParamValues)->begin();i!=(right.m_pParamValues)->end();i++){
			m_pParamValues->push_back( *i ); 
		}
	}	
	return (*this);	
}


CMagStat* CIdentStat::FindMagInfo( const char* szMag )
{
	vector<CMagStat>::iterator i;

	for(i=m_pMagnitTab->begin();i!=m_pMagnitTab->end();i++){
		if( strcmp( szMag, i->m_Magnitude.c_str() )==0 )
			return &(*i);		
	}
	return NULL;
}

void CIdentStat::AddOrUpdateMagInfo( CccdReport& evt )
{
	CMagStat* pMagInfo = FindMagInfo( evt.m_Magnitude.c_str() );
	if (pMagInfo){
		// found - update
		pMagInfo->m_MagGenerated++;
		if( evt.IsIdentified() )
			pMagInfo->m_MagIdentified++;
	}else{
		CMagStat tmp;
		tmp.m_Magnitude = evt.m_Magnitude.c_str();
		tmp.m_MagGenerated++;
		if( evt.IsIdentified() )
			tmp.m_MagIdentified++;
		m_pMagnitTab->push_back( tmp );
	}	

}

void CIdentStat::AddOrUpdateMagInfo( const char* szMag, BOOL_T bIdent)
{
	CMagStat* pMagInfo = FindMagInfo( szMag );
	if (pMagInfo){
		// found - update
		pMagInfo->m_MagGenerated++;
		if(bIdent)
			pMagInfo->m_MagIdentified++;
	}else{
		CMagStat tmp;
		tmp.m_Magnitude = szMag;
		tmp.m_MagGenerated++;
		if(bIdent)
			tmp.m_MagIdentified++;
		m_pMagnitTab->push_back( tmp );
	}	
}


void CIdentStat::Reset()
{
	m_pMagnitTab->clear();
	if(m_pParamValues)
		m_pParamValues->clear();
}

vector<CEnvVar>& CIdentStat::GetParamTab()
{
	return (*m_pParamValues);
}



// insertion sort :
void CTrackDesc::SortEventListByTime()
{
	register int size = m_EventsOnTrack.size();
	for( register int i=0;i<size;i++){
		time_t minValue = m_EventsOnTrack[i].m_Time;
		int minPos=i;
		for( register int j=(i+1);j<size;j++){
			if( m_EventsOnTrack[j].m_Time < minValue ){
				minValue = m_EventsOnTrack[j].m_Time;
				minPos = j;
			}
		}
		if( minPos != i ){
			// replace min found with i element :
			CEventBaseInfo evt = m_EventsOnTrack[i];
			m_EventsOnTrack[i] = m_EventsOnTrack[minPos];
			m_EventsOnTrack[minPos] = evt;
		}
	}
}

void CTrackDesc::FlagEventsSaved()
{
	register int size = m_EventsOnTrack.size();
   for( register int i=0;i<size;i++){
		m_EventsOnTrack[i].m_bSavedToDB = TRUE;
	}
}

void CTrackDesc::GetEventsList( mystring& szList )
{
	szList = "";
	vector<CEventBaseInfo>::iterator i;
	for(i=m_EventsOnTrack.begin();i!=m_EventsOnTrack.end();i++){
		szList << (int)i->m_DayFrameIndex << "-(" 
				 << (int)i->m_MaxPoint.x << ","
         	 << (int)i->m_MaxPoint.y << "),";
	}	
}

void CTrackDesc::GetEventsList_RADEC( mystring& szList, BOOL_T bDoConvert )
{
	szList = "";
	vector<CEventBaseInfo>::iterator i;
	char szRADEC[128];
	

	for(i=m_EventsOnTrack.begin();i!=m_EventsOnTrack.end();i++){
		if( bDoConvert ){
			if(sprintf(szRADEC,"%.4f,%.4f",AstroAngle::rad2deg(i->m_AstroCoord.RA),
					AstroAngle::rad2deg(i->m_AstroCoord.Dec) )>=128 ){
				printf("ERROR in CTrackDesc::GetEventsList_RADEC\n");
				exit(0);
			}
		}else{
			if(sprintf(szRADEC,"%.4f,%.4f",i->m_AstroCoord.RA,i->m_AstroCoord.Dec )>=128 ){
				printf("ERROR in CTrackDesc::GetEventsList_RADEC - NO CONVERT\n");
				exit(0);
			}
		}


		szList << (int)i->m_DayFrameIndex << "-(" 
				 << (int)i->m_MaxPoint.x << ","
         	 << (int)i->m_MaxPoint.y << ","
				 << szRADEC << "," 
				 << i->m_Time
				 << "),";
	}	
}



void CTrackDesc::get_frame_range( int& min_track, int& max_track )
{
	min_track = 0;
	max_track = 0;
	if( m_EventsOnTrack.size() ){
		min_track = m_EventsOnTrack[0].m_DayFrameIndex;
		max_track = m_EventsOnTrack.back().m_DayFrameIndex;
	}
}

void CTrackDesc::GetSortedEventsList( mystring& szList )
{
	szList="";
	LONG_T* pFramesList = new LONG_T[ m_EventsOnTrack.size() ];
	vector<CEventBaseInfo>::iterator i;
	int count=0;
	for(i=m_EventsOnTrack.begin();i!=m_EventsOnTrack.end();i++){
		if( find_value( pFramesList, count, i->m_FrameIndex )<0 ){
			pFramesList[count] = i->m_FrameIndex;
			count++;
		}
	}	
	my_qsort( pFramesList , m_EventsOnTrack.size() );

	for(int pos=0;pos<count;pos++){
		for(i=m_EventsOnTrack.begin();i!=m_EventsOnTrack.end();i++){
			if( i->m_FrameIndex == pFramesList[pos] ){
				szList << (int)i->m_FrameIndex << "-(" 
						 << (int)i->m_MaxPoint.x << ","
            		 << (int)i->m_MaxPoint.y << "),";
			}
		}
	}

	delete [] pFramesList;
}

int CTrackDesc::GetMemSize( BOOL_T bShowAll )
{
	int ret = sizeof(CTrackDesc);
	for(int i=0;i<m_EventsOnTrack.size();i++){
		int evt_size =  0; //m_EventsOnTrack[i].GetMemSize( bShowAll );

		if(bShowAll){
			printf("\t\tevent %d takes %d bytes\n",i,evt_size);
		}
		ret += evt_size;
	}
	return ret;
}


BOOL_T CTrackDesc::operator==( const CTrackDesc& right )
{
	int a_int_this = (int)(a*100);
   int b_int_this = (int)(b*100);

	int a_int_right = (int)(right.a*100);
   int b_int_right = (int)(right.b*100);

	return (( track_id==right.track_id ) || (	a_int_this==a_int_right && b_int_this==b_int_right && m_cam_idx==right.m_cam_idx ));
}

BOOL_T CTrackDesc::CheckEvent( int x, int y, double& chi2, double chi2_limit )
{
	chi2 = CMyFit::CalcDist2FromLine2Par( a, b, x, y );
	if( chi2 < chi2_limit )
		return TRUE;
	return FALSE;
}

BOOL_T CTrackDesc::CheckEvent_RADEC( double ra, double dec, 
												 double& chi2, double chi2_limit )
{
	if( m_bRADEC ){
		chi2 = CMyFit::CalcDist2FromLine2Par( radec_a, radec_b, ra, dec );
		if( chi2 < chi2_limit )
			return TRUE;
		return FALSE;
	}
	return TRUE;
}


BOOL_T CTrackDesc::CheckMove( int x, int y, int frame )
{
	if( DoesEventBelongs( x, y, frame ) )
		return TRUE;

	if( m_EventsOnTrack.size()>0 ){
		if( frame <= m_EventsOnTrack[0].m_DayFrameIndex ){
			// before any of events on track :
			if( frame==m_EventsOnTrack[0].m_DayFrameIndex ){
				return FALSE; // does not belong to track - see above 
			}else{
				double _vx = ( m_EventsOnTrack[0].m_MaxPoint.x - x );
				if( vx*_vx < 0 )
					return FALSE;
				double _vy = ( m_EventsOnTrack[0].m_MaxPoint.y - y );
				if( vy*_vy < 0 )
					return FALSE;
			}
		}
		if( frame >= m_EventsOnTrack.back().m_DayFrameIndex ){
			// after all :
			if( frame==m_EventsOnTrack.back().m_DayFrameIndex ){
				return FALSE; // does not belong to track - see above
			}else{
				double _vx = ( x - m_EventsOnTrack.back().m_MaxPoint.x );
				if( vx*_vx < 0 )
					return FALSE;
				double _vy = ( y - m_EventsOnTrack.back().m_MaxPoint.y );
				if( vy*_vy < 0 )
					return FALSE;
			}
		} 
		if( frame > m_EventsOnTrack[0].m_DayFrameIndex && frame < m_EventsOnTrack.back().m_DayFrameIndex ){
			double min_dist=5;
			// here we are sure m_EventsOnTrack[0].size()>=2 
			// inside :

			//double stepX=(m_EventsOnTrack[1].m_MaxPoint.x-m_EventsOnTrack[0].m_MaxPoint.x);
			//double stepY=(m_EventsOnTrack[1].m_MaxPoint.y-m_EventsOnTrack[0].m_MaxPoint.y);

			BOOL_T bOK=FALSE;
			int evt_cnt = m_EventsOnTrack.size();
			for(int i=0;i<(evt_cnt-1);i++){
				if( frame>m_EventsOnTrack[i].m_DayFrameIndex && frame<m_EventsOnTrack[i+1].m_DayFrameIndex ){
					if( CPoint::IsBetween( m_EventsOnTrack[i].m_MaxPoint.x, m_EventsOnTrack[i+1].m_MaxPoint.x, x ) &&
 					    CPoint::IsBetween( m_EventsOnTrack[i].m_MaxPoint.y, m_EventsOnTrack[i+1].m_MaxPoint.y, y ) 
					){
						bOK=TRUE;
						break;
					}
				}else{
					if( frame==m_EventsOnTrack[i].m_DayFrameIndex ){
						if( i>0 ){
							if(   CPoint::IsBetween( m_EventsOnTrack[i-1].m_MaxPoint.x, m_EventsOnTrack[i+1].m_MaxPoint.x, x ) &&
								   CPoint::IsBetween( m_EventsOnTrack[i-1].m_MaxPoint.y, m_EventsOnTrack[i+1].m_MaxPoint.y, y ) 
								  ){
								bOK=TRUE;						
								break;
							}
						}						
					}
				}
			}
			return bOK;
		}
	}
	return TRUE;
}

BOOL_T CTrackDesc::CheckMove_RADEC(  int x, int y,
												 double ra, double dec, int frame )
{
	if( DoesEventBelongs( x, y, frame ) )
		return TRUE;

	if( m_EventsOnTrack.size()>0 ){
		if( frame <= m_EventsOnTrack[0].m_DayFrameIndex ){
			// before any of events on track :
			if( frame==m_EventsOnTrack[0].m_DayFrameIndex ){
				return FALSE; // does not belong to track - see above 
			}else{
				double _v_ra = ( m_EventsOnTrack[0].m_AstroCoord.RA - ra );
				if( v_ra*_v_ra < 0 )
					return FALSE;
				double _v_dec = ( m_EventsOnTrack[0].m_AstroCoord.Dec - dec );
				if( v_dec*_v_dec < 0 )
					return FALSE;
			}
		}
		if( frame >= m_EventsOnTrack.back().m_DayFrameIndex ){
			// after all :
			if( frame==m_EventsOnTrack.back().m_DayFrameIndex ){
				return FALSE; // does not belong to track - see above
			}else{
				double _v_ra = ( ra - m_EventsOnTrack.back().m_AstroCoord.RA );
				if( v_ra*_v_ra < 0 )
					return FALSE;
				double _v_dec = ( dec - m_EventsOnTrack.back().m_AstroCoord.Dec );
				if( v_dec*_v_dec < 0 )
					return FALSE;
			}
		} 
		if( frame > m_EventsOnTrack[0].m_DayFrameIndex && frame < m_EventsOnTrack.back().m_DayFrameIndex ){
			double min_dist=5;
			// here we are sure m_EventsOnTrack[0].size()>=2 
			// inside :


			BOOL_T bOK=FALSE;
			int evt_cnt = m_EventsOnTrack.size();
			for(int i=0;i<(evt_cnt-1);i++){
				if( frame>m_EventsOnTrack[i].m_DayFrameIndex && frame<m_EventsOnTrack[i+1].m_DayFrameIndex ){
					if( CPoint::IsBetween( m_EventsOnTrack[i].m_AstroCoord.RA, m_EventsOnTrack[i+1].m_AstroCoord.RA, ra ) &&
 					    CPoint::IsBetween( m_EventsOnTrack[i].m_AstroCoord.Dec, m_EventsOnTrack[i+1].m_AstroCoord.Dec, dec ) 
					){
						bOK=TRUE;
						break;
					}
				}else{
					if( frame==m_EventsOnTrack[i].m_DayFrameIndex ){
						if( i>0 ){
							if(   CPoint::IsBetween( m_EventsOnTrack[i-1].m_AstroCoord.RA, m_EventsOnTrack[i+1].m_AstroCoord.RA, ra ) &&
								   CPoint::IsBetween( m_EventsOnTrack[i-1].m_AstroCoord.Dec, m_EventsOnTrack[i+1].m_AstroCoord.Dec, dec ) 
								  ){
								bOK=TRUE;						
								break;
							}
						}						
					}
				}
			}
			return bOK;
		}
	}
	return TRUE;
}


BOOL_T CTrackDesc::DoesEventBelongs( int x, int y, int frame, double r )
{
	vector<CEventBaseInfo>::iterator i;
	for(i=m_EventsOnTrack.begin();i!=m_EventsOnTrack.end();i++){
		if( i->m_DayFrameIndex == frame ){
			double d = sqrt( (i->m_MaxPoint.x-x)*(i->m_MaxPoint.x-x) + (i->m_MaxPoint.y-y)*(i->m_MaxPoint.y-y) );
		
			if( d<=r ){
				return TRUE;
			}
		}
	}		
	return FALSE;
}

BOOL_T CTrackDesc::Parse( const char* szTrack, int len )
{
	// 20041020_235616 0       41      -0.79 2110.65 28 -0.27 0.22 N 36-(845,1445),
	int camno;
	char szDTM[32],type;
	char* p_events_list = new char[ len ];
	int ret = sscanf( szTrack, "%s %d %d %lf %lf %d %lf %lf %c %s",
				szDTM, &camno, &frame_index, &a, &b, &track_id,
				&vx, &vy, &type, p_events_list );
	MyParser pars=p_events_list;
	const char* szTmp = NULL;
	while( (szTmp = pars.GetItemUpToString( ")," )) ){
		int x,y,f;
		int ret1 = sscanf( szTmp,"%d-(%d,%d)",&f,&x,&y);
		if( ret1 == 3 ){
			CEventBaseInfo evt( x,y,x,y,f,f);
			m_EventsOnTrack.push_back( evt );
		}
	}	
	delete [] p_events_list;
	return (ret==10);
}

BOOL_T CTrackDesc::ParseWithRADEC( const char* szTrack, int len )
{
	// # DTM   CCD#    FRAME#  a b ID vx vy radec_a radec_b v_ra v_dec Oper
	// Events_on_Track

	// 20050107_144434 0       0       -0.0400 1501.6100 49 -0.2500 0.0100 -0.0527 
	// 0.2177 0.0001 -0.0000 N     42-(1133,1459,88.2431,7.8100)
	int camno;
	char szDTM[32],type;
	char* p_events_list = new char[ len ];
	int ret = sscanf( szTrack, "%s %d %d %lf %lf %d %lf %lf %lf %lf %lf %lf %c %s",
				szDTM, &camno, &frame_index, &a, &b, &track_id,
				&vx, &vy, &radec_a, &radec_b, &v_ra, &v_dec, 
				&type, p_events_list );
	//if( track_id==3075 ){
	//	printf("odo");
	//}
	MyParser pars=p_events_list;
	const char* szTmp = NULL;
	while( (szTmp = pars.GetItemUpToString( ")," )) ){
		int x,y,f;
		double ra,dec;
		time_t unix_time=0;
		int ret1 = sscanf( szTmp,"%d-(%d,%d,%lf,%lf,%d)",
							&f,&x,&y,&ra,&dec,&unix_time);
		if( ret1<6 ){
			ret1 = sscanf( szTmp,"%d-(%d,%d,%lf,%lf)",&f,&x,&y,&ra,&dec);
		}

		if( ret1 >=6  ){
			CEventBaseInfo evt( x,y,x,y,f,f);
			evt.m_AstroCoord.RA = ra;
			evt.m_AstroCoord.Dec = dec;
			evt.m_Time = unix_time;
			m_EventsOnTrack.push_back( evt );
		}
	}	
	delete [] p_events_list;
	return (ret==14);
}


CTrackDesc* CTrackList::CheckEvent( int x, int y, double& min_chi2, 
												double chi2_limit, BOOL_T bDoCheck )
{
	min_chi2=1000000.00;
	double chi2;
	CTrackDesc* pBestTrack=NULL;
	CTrackList::iterator i;
	for(i=begin();i!=end();i++){
		i->CheckEvent( x, y, chi2, chi2_limit );
		if( chi2<min_chi2 ){
			min_chi2 = chi2;
			pBestTrack = &(*i);
		}
	}

	if( bDoCheck ){
		if( min_chi2 > chi2_limit ){
			pBestTrack=NULL;
		}
	}

	return pBestTrack;
}

CTrackDesc* CTrackList::CheckEvent_RADEC( int x, int y, double ra, double dec,
												double& min_chi2, 
												double chi2_limit, BOOL_T bDoCheck )
{
	min_chi2=1000000.00;
	double chi2;
	CTrackDesc* pBestTrack=NULL;
	CTrackList::iterator i;
	for(i=begin();i!=end();i++){
		i->CheckEvent_RADEC( ra, dec, chi2, chi2_limit );
		if( chi2<min_chi2 ){
			min_chi2 = chi2;
			pBestTrack = &(*i);
		}
	}

	if( bDoCheck ){
		if( min_chi2 > chi2_limit ){
			pBestTrack=NULL;
		}
	}

	return pBestTrack;
}


int CTrackList::GetBestTracks( int x, int y, int frame,
										 int n, CTrackList& bestlist )
{
	bestlist.clear();
	for(int i=0;i<n;i++){
		double min_chi2=1000000.00;
		CTrackDesc* pBestTrack=NULL;
		CTrackList::iterator it;
		for(it=begin();it!=end();it++){
			double chi2;
			it->CheckEvent( x, y, chi2, 1000000.00 );	
			if( chi2<min_chi2 && !bestlist.Find( it->track_id ) ){
 	      	min_chi2 = chi2;
         	pBestTrack = &(*it);				
	      }
		}
		if( pBestTrack ){
			bestlist.push_back( *pBestTrack );
			bestlist.back().chi2 = min_chi2;
			bestlist.back().bMoveOK = bestlist.back().CheckMove( x, y ,frame );
		}		
	}		

	for( int i=0;i<bestlist.size();i++ ){
		bestlist[i].CheckMove( x, y, frame );
	}

	return bestlist.size();
}

int CTrackList::GetBestTracks_RADEC( double ra, double dec,
										 int x, int y, int frame,
										 int n, CTrackList& bestlist )
{
	bestlist.clear();
	for(int i=0;i<n;i++){
		double min_chi2=1000000.00;
		CTrackDesc* pBestTrack=NULL;
		CTrackList::iterator it;
		for(it=begin();it!=end();it++){
			double chi2;
			it->CheckEvent_RADEC( ra, dec, chi2, 1000000.00 );	
			if( chi2<min_chi2 && !bestlist.Find( it->track_id ) ){
 	      	min_chi2 = chi2;
         	pBestTrack = &(*it);				
	      }
		}
		if( pBestTrack ){
			bestlist.push_back( *pBestTrack );
			// bestlist.back().chi2 = CMyMathFunc::mysqr((sqrt(min_chi2)*3600)/gCCDParams.m_fPixScale);
			bestlist.back().chi2 = min_chi2;
			bestlist.back().bMoveOK = bestlist.back().CheckMove_RADEC( x, y ,ra, dec, frame );
		}		
	}		

	for( int i=0;i<bestlist.size();i++ ){
		bestlist[i].CheckMove_RADEC( x, y, ra, dec, frame );
	}

	return bestlist.size();
}


int CTrackList::MergeRADECInfo( CCDEventList& event_list, BOOL_T bSaveRADEC )
{
	int ret=0;
	for( int i=0;i<size();i++){
		CTrackDesc& track = (*this)[i];

		int cnt=0;
		double* ra_tab = new double[track.m_EventsOnTrack.size()+1];
		double* dec_tab = new double[track.m_EventsOnTrack.size()+1];

		for(int j=0;j<track.m_EventsOnTrack.size();j++){
			CEventBaseInfo& evt = (track.m_EventsOnTrack)[j];

			CccdReport* pEvt = event_list.FindEvent( (int)evt.m_MaxPoint.x, (int)evt.m_MaxPoint.y, evt.m_DayFrameIndex, 1 );
			if( !pEvt ){
				pEvt = event_list.FindEvent( (int)evt.m_Point.x, (int)evt.m_Point.y, evt.m_DayFrameIndex, 1 );
			}
			if ( pEvt ){
				evt.m_AstroCoord = pEvt->m_AstroCoord;
				evt.m_Time = pEvt->m_Time;
				// evt.m_AstroCoord.RA = evt.m_AstroCoord.RA;
				// evt.m_AstroCoord.Dec = evt.m_AstroCoord.Dec;
				
				// NEW 20050120 :
				evt.m_AstroCoord.RA = AstroAngle::rad2deg( evt.m_AstroCoord.RA );
				evt.m_AstroCoord.Dec = AstroAngle::rad2deg( evt.m_AstroCoord.Dec );
				ret++;

				// fit must be in DEGREES :
				// NEW 20050120 :
				// ra_tab[cnt] = AstroAngle::rad2deg( evt.m_AstroCoord.RA );
				// dec_tab[cnt] = AstroAngle::rad2deg( evt.m_AstroCoord.Dec );
				ra_tab[cnt] = evt.m_AstroCoord.RA;
				dec_tab[cnt] = evt.m_AstroCoord.Dec;

				cnt++;
			}
		}
		if( cnt>=3 ){
			double chi2 = CMyFit::FitLineChi2( ra_tab, dec_tab, cnt, track.radec_a, track.radec_b );
			double chi2_deg = sqrt(chi2);
			double chi2_pixel = (chi2_deg*3600.00)/gCCDParams.m_fPixScale;
			
			// if( chi2_pixel < gCCDParams.m_MaxChi2InTrack ){
			printf("RADEC line (track=%d): a=%.2f , b=%.2f , chi2=%.2f (chi2_pixel=%.5f, limit=%.5f), cnt=%d\n",
						track.track_id,track.radec_a,track.radec_b,chi2,
						chi2_pixel, gCCDParams.m_MaxChi2InTrack,
						cnt);
			track.m_bRADEC = TRUE;
			track.v_ra = ( track.m_EventsOnTrack.back().m_AstroCoord.RA - track.m_EventsOnTrack.front().m_AstroCoord.RA ) / (track.m_EventsOnTrack.back().m_Time-track.m_EventsOnTrack.front().m_Time);
			track.v_dec = ( track.m_EventsOnTrack.back().m_AstroCoord.Dec - track.m_EventsOnTrack.front().m_AstroCoord.Dec ) / (track.m_EventsOnTrack.back().m_Time-track.m_EventsOnTrack.front().m_Time);
		}
		delete [] ra_tab;
		delete [] dec_tab;

		if( bSaveRADEC ){
			CCD_Analyser::LogNewTrack_RADEC( track, eNormalTrack, "newtrackslog_radec_ccd0.txt", 
														TRUE, NULL, TRUE, NULL, FALSE, FALSE );
		}
	}
	return ret;
}

int CTrackList::ReadWithRADEC( const char* szFileName )
{
	clear();
	if( MyFile::DoesFileExist( szFileName )){
		MyIFile in( szFileName );

		const char* pLine;
		while( pLine = in.GetLine() ){
      	if( mystring::get_first_non_white( pLine )=='#' )
         	continue;
			CTrackDesc track;
			BOOL_T bParse =  track.ParseWithRADEC( pLine, strlen( pLine) );
			track.m_bRADEC = TRUE;			

			if( bParse ){
				CTrackDesc* pTR = Find( track.track_id );
				if( pTR ){
					// update
					(*pTR) = track;
				}else{
					// add new
					push_back( track );
				}
			}
		}		
	}	
	return size();
	
}

int CTrackList::Read( const char* szFileName, BOOL_T bRADEC, 
							 const char* szVerifFile, BOOL_T bSaveRADEC )
{
	clear();
	if( MyFile::DoesFileExist( szFileName )){
		MyIFile in( szFileName );

		const char* pLine;
		while( pLine = in.GetLine() ){
      	if( mystring::get_first_non_white( pLine )=='#' )
         	continue;
			CTrackDesc track;
			BOOL_T bParse =  track.Parse( pLine, strlen( pLine) );
			if( bParse ){
				CTrackDesc* pTR = Find( track.track_id );
				if( pTR ){
					// update
					(*pTR) = track;
				}else{
					// add new
					push_back( track );
				}
			}
		}		

		if( size()>0 ){
			printf("Read %d tracks\n",size());
			if( bRADEC ){
				printf("reading RA,DEC info using file : %s\n",szVerifFile);
				if( MyFile::DoesFileExist( szVerifFile ) ){
					CCDEventList verifEvents;
					verifEvents.Read( szVerifFile );					
					MergeRADECInfo( verifEvents, bSaveRADEC );
				}else{
					printf("File %s NOT FOUND cannot read RA,DEC info !!!\n");
					exit(0);
				}
			}
		}
	}	
	return size();
}

void CTrackList::remove( int min_pos, int max_pos )
{
	
}


int CTrackList::GetMemSize( BOOL_T bShowAll )
{
	int ret = 0;
	for(int i=0;i<size();i++){
		ret += (*this)[i].GetMemSize( bShowAll );
	}
	return ret;
}

void CTrackList::Add( double a, double b, int cam_idx )
{
	CTrackList::iterator i;
	CTrackDesc tmp( a,b,cam_idx );
	Add( tmp );
}

CTrackDesc* CTrackList::Find( vector<CTrackDesc*>& trackList, int trackID )
{
	vector<CTrackDesc*>::iterator i;
	for(i=trackList.begin();i!=trackList.end();i++){
		if( (*i)->track_id == trackID )
			return ((CTrackDesc*)(*i));
	}
	return NULL;
}

CTrackDesc* CTrackList::Find( int trackID )
{
	CTrackList::iterator i;
	for(i=begin();i!=end();i++){
		if(i->track_id == trackID )
			return &(*i);
	}
	return NULL;
}

CTrackDesc* CTrackList::Find( const CTrackDesc& newelem )
{
	CTrackList::iterator i;
	for(i=begin();i!=end();i++){
		if( (*i)==newelem )
			return &(*i);
	}
	return NULL;
}

void CTrackList::Add( const CTrackDesc& newelem )
{
	if(Find( newelem ))
		return;
	push_back(newelem);
}

CTrackList::~CTrackList()
{
	clear();
}


                                                                                
int CTrackDesc::m_TrackIDGenerator=1;

CTrackDesc::CTrackDesc( double aa, double bb, int cam_idx )
 : a(aa), b(bb), m_cam_idx(cam_idx), track_id(0), vx(0), vy(0),
   bMoveOK(TRUE),chi2(0),m_bRADEC(FALSE),radec_a(0),radec_b(0),
   v_ra(0),v_dec(0),tr_night(0),frame_index_end(0)
{
	track_id = m_TrackIDGenerator;
	m_TrackIDGenerator++;
}

CTrackDesc::~CTrackDesc()
{
	m_EventsOnTrack.clear();
}

CTrackDesc::CTrackDesc( const CTrackDesc& right )
{
	(*this) = right;
}

double CTrackDesc::FindMinDistByMaxPixel( double x, double y )
{
	vector<CEventBaseInfo>::iterator i;
	double minDist=100000.00;
	for(i=m_EventsOnTrack.begin();i!=m_EventsOnTrack.end();i++){
		double dist = CPoint::dist( i->m_MaxPoint.x, i->m_MaxPoint.y, x, y );
		if(dist<minDist){
			minDist=dist;
		}
	}
	return minDist;
}

double CTrackDesc::FindMinDist( double x, double y )
{
	vector<CEventBaseInfo>::iterator i;
	double minDist=100000.00;
	for(i=m_EventsOnTrack.begin();i!=m_EventsOnTrack.end();i++){
		double dist = CPoint::dist( i->m_Point.x, i->m_Point.y, x, y );
		if(dist<minDist){
			minDist=dist;
		}
	}
	return minDist;
}

void CTrackDesc::UpdateTrack( double new_a, double new_b, CccdReport& new_evt )
{
	a = new_a;
	b = new_b;
	m_EventsOnTrack.push_back( new_evt );
}

void CTrackDesc::calcAverageVelocity(){
	calcAverageVelocity( vx, vy );
}


void CTrackDesc::calcAverageVelocity( double& _vx, double& _vy )
{
	_vx = 0;
	_vy = 0;
	if( m_EventsOnTrack.size()<=1 )
		return;
	time_t startTime=m_EventsOnTrack[0].m_Time;
	time_t endTime=m_EventsOnTrack.back().m_Time;
	vector<CEventBaseInfo>::iterator i;	
	CEventBaseInfo* pStart=(&(m_EventsOnTrack[0]));
	CEventBaseInfo* pEnd=(&(m_EventsOnTrack.back()));
	for(i=m_EventsOnTrack.begin();i!=m_EventsOnTrack.end();i++){
		if(i->m_Time>endTime){
			endTime = i->m_Time;
			pEnd = &(*i);
		}
		if(i->m_Time<startTime){
			startTime =  i->m_Time;
			pStart = &(*i);
		}
	}	
	time_t dt = (endTime-startTime);
	if( dt==0 ){
		_vx = 0;
		_vy = 0;
	}else{
		// dt<>0
		if( pEnd && pStart ){
			_vx = ( (pEnd->m_MaxPoint).x - (pStart->m_MaxPoint).x )/dt;
			_vy = ( (pEnd->m_MaxPoint).y - (pStart->m_MaxPoint).y )/dt;
		}
	}
}

BOOL_T CTrackDesc::CheckVelocity( CEventBaseInfo* pEvt1, CEventBaseInfo* pEvt2, 
											 double fVeloError,
										    int idx )
{
	if( pEvt1 && pEvt2 ){
		time_t dt = (pEvt2->m_Time-pEvt1->m_Time);
		if( dt==0 )
			return FALSE;
		double vx_new = ( (pEvt2->m_MaxPoint).x	- (pEvt1->m_MaxPoint).x );
		double vy_new = ( (pEvt2->m_MaxPoint).y - (pEvt1->m_MaxPoint).y );

		if( CCD_Analyser::m_pAnalFoundEventFunc ){
			CCD_Analyser::FillVelocityHisto( vx, vy, vx_new, vy_new, idx );
		}

		if( !CMyFit::CheckVelocityCondition( vx, vy, vx_new, vy_new, fVeloError ) ){
			return FALSE;
		}
	}
	return TRUE;
}

BOOL_T CTrackDesc::CheckVelocity( CEventBaseInfo* pNewEvent, double velError,
											 int idx )
{
	if( m_EventsOnTrack.size()<2 )
      return TRUE;
   CEventBaseInfo* pBefore=NULL;
   CEventBaseInfo* pAfter=NULL;
   time_t beforeDT=1000000;
   time_t  afterDT=1000000;
	vector<CEventBaseInfo>::iterator i;
	for(i=m_EventsOnTrack.begin();i!=m_EventsOnTrack.end();i++){
		int dt = (pNewEvent->m_Time - i->m_Time);
      if( dt>0 && dt<beforeDT ){
         beforeDT = dt;
         pBefore = &(*i);
      }else{
         dt = -dt;
         if( dt>0 && dt<afterDT ){
            afterDT = dt;
				pAfter = &(*i);
         }
      }
	}
	if( pBefore || pAfter ){
		if( pBefore && !CheckVelocity( pBefore, pNewEvent, velError, idx ) ){
			return FALSE;
		}
		if( pAfter && !CheckVelocity( pNewEvent, pAfter, velError, idx ) ){
			return FALSE;
		}
	}

	return TRUE;
}



void CTrackDesc::SetTrackDesc( double aa, double bb, int cam_idx )
{
	a = aa;
	b = bb;
	m_cam_idx = cam_idx;
}

CTrackDesc& CTrackDesc::operator=( const CTrackDesc& right )
{
	a = right.a;
	b = right.b;
	m_cam_idx = right.m_cam_idx;
	track_id = right.track_id;
	tr_night =  right.tr_night;
	frame_index_end = right.frame_index_end;

	m_EventsOnTrack.clear();
	for(register int i=0;i<right.m_EventsOnTrack.size();i++){
		m_EventsOnTrack.push_back( right.m_EventsOnTrack[i] );
	}
	vx = right.vx;
	vy = right.vy;
	frame_index = right.frame_index;
	chi2 = right.chi2;
	bMoveOK = right.bMoveOK;

	radec_a = right.radec_a;
	radec_b = right.radec_b;
	v_ra = right.v_ra;
	v_dec = right.v_dec;
	m_bRADEC = right.m_bRADEC;
	

	return (*this);
}

int CTrackDesc::get_xy_values( double* x_values, double* y_values ) const
{
	for(register int i=0;i<m_EventsOnTrack.size();i++){
		x_values[i] = (m_EventsOnTrack[i]).m_Point.x;
		y_values[i] = (m_EventsOnTrack[i]).m_Point.y;
	}
	return m_EventsOnTrack.size();
}

int CTrackDesc::get_radec_values( double* ra_values, double* dec_values ) const
{
   for(register int i=0;i<m_EventsOnTrack.size();i++){
      ra_values[i] = (m_EventsOnTrack[i]).m_AstroCoord.RA;
      dec_values[i] = (m_EventsOnTrack[i]).m_AstroCoord.Dec;
   }
   return m_EventsOnTrack.size();
}


void CTrackDesc::LogEventCheck( CccdReport& evt, int cam_idx, double chi2,
										  double vx, double vy,
										  double rx, double ry, BOOL_T bOK,
										  double minDist,
										  const char* szReason )
{
	BOOL_T bHeader=FALSE;
	mystring szFileName;
	szFileName << gCCDParams.GetOutputDir() << "/" << TRACKS_SUBDIR << "/Cam" 
				  << cam_idx << "/addevt.log";		
	if( !MyFile::DoesFileExist( szFileName.c_str() ) ){
		bHeader = TRUE;
	}
	MyOFile out( szFileName.c_str(), "a+" );
	if( bHeader ){
		out.Printf("%s",CTrackDesc::m_LogHeader.c_str());
	}
	out.Printf(CTrackDesc::m_LogFmt.c_str(),evt.m_DayFrameIndex,
		(int)evt.m_MaxPoint.x,(int)evt.m_MaxPoint.y,track_id,chi2,
		vx,vy,rx,ry,(int)bOK,minDist,szReason);
}

BOOL_T CccdReport::SaveGenEvent( const char* fname,
									 int x, int y,
	                         const char* mag,
	                         Table2D<ELEM_SAMPLE_TYPE>& OtherImage, 
									 int frame_index, int max_lap )
{
	LONG_T x_s = x;
	LONG_T y_s = y;
	int start_x = (x_s - OtherImage.m_MaxX);
	int start_y = (y_s - OtherImage.m_MaxY);
	CPoint LowLeft( start_x, start_y );
	CPoint TopRight( (start_x + OtherImage.GetXSize()) , (start_y + OtherImage.GetYSize()));

	CPixelAnalyseIn in;
	CPixelAnalyseOut out;

	in.x = x_s;
	in.y = y_s;
	in.frame_index = frame_index;

	out.m_PixelOut.x0 = x_s;
	out.m_PixelOut.y0 = y_s;

	CccdReport tmp(0,	in, out, TRUE, FALSE, &LowLeft,&TopRight,mag);
	tmp.m_PixelAnalResults.laplaceSum = max_lap;

	long maxX,maxY;
	// tmp.m_MaxValue = GetMaxValueAndPos( start_x,start_y, end_x, end_y,
	//				  							   maxX, maxY );
	_TRACE_PRINTF_5("Other Image : (%d,%d)\n",OtherImage.m_MaxX,OtherImage.m_MaxY);

	tmp.m_MaxPoint.x = x_s;
   tmp.m_MaxPoint.y = y_s;
	tmp.m_Time = get_dttm();
	// tmp.m_MaxValue = m_pFastData[(long)tmp.m_MaxPoint.y][(long)tmp.m_MaxPoint.x];

	_TRACE_PRINTF_4("\n\n PUT IMAGE at (%d,%d) MAX:(max_x,max_y) = (%f,%f)\n\n",x,y,tmp.m_MaxPoint.x,tmp.m_MaxPoint.y);

	BOOL_T bNewFile=FALSE;
	if( MyFile::GetFileSize( fname )<=0 ){
		bNewFile = TRUE;
	}

	MyOFile outfile(fname,"a");
	if( bNewFile ){
		mystring szHeader;
		CccdReport::GetDescHeader(szHeader);
		outfile.Printf("%s",szHeader.c_str());
	}
	
	char line[2048];
	tmp.GetEventDescStr( line );
	outfile.Printf("%s",line);
	outfile.Close();

	return TRUE;		
}

// assums sorted by frame_index :	
CccdReport* CCDEventList::FindByFrameIndex( int frame_index, int& pos )
{
	pos=-1;
	for(register int i=0;i<size();i++){
		if(frame_index== (*this)[i].m_FrameIndex){
			pos = i;
			return &((*this)[i]);
		}else{
			if((*this)[i].m_FrameIndex>frame_index)
				break;
		}
	}

	return NULL;
}

int CCDEventList::GetEventsByFrameIndex( CCDEventList& out, int frame_index, int end_frame )
{
	out.clear();
	

	if( end_frame<=0 ){
		int pos;
		CccdReport* pFirst = FindByFrameIndex( frame_index, pos );

		if(pFirst){
			for(int i=pos;i<size();i++){
				if( (*this)[i].m_FrameIndex == frame_index ){
					out.push_back( (*this)[i] );
				}else{
					break;			
				}
			}
		}
	}else{
		for(int i=0;i<=size();i++){
			if( (*this)[i].m_FrameIndex>=frame_index && (*this)[i].m_FrameIndex<=end_frame){
				out.push_back( (*this)[i] );
			}
		}
	}
	
	return out.size();
}

double CccdReport::CalcDist( CccdReport& evt1, CccdReport& evt2 )
{
	double dist_in_rad = sqrt( (evt1.m_AstroCoord.RA-evt2.m_AstroCoord.RA)*(evt1.m_AstroCoord.RA-evt2.m_AstroCoord.RA)+
										(evt1.m_AstroCoord.Dec-evt2.m_AstroCoord.Dec)*(evt1.m_AstroCoord.Dec-evt2.m_AstroCoord.Dec) );
	return dist_in_rad;
}

BOOL_T CccdReport::IsGoodToSave( CccdReport& evt )
{
	if( !evt.m_PixelAnalResults.m_bRejectedByNextFrames && !evt.m_PixelAnalResults.m_bRejectedDueToTrack &&
		 !evt.m_PixelAnalResults.m_bRejectedByCoic && 
		 ( !gCCDParams.m_bSaveOnlyGood || CCD_Analyser::VerifyEvent( evt ) ) && 
       (!gCCDParams.m_bSaveSuperNovaOnly || evt.m_PixelAnalResults.eventType==eBrighten ) ){
		return TRUE;
	}
	return FALSE;		 
}

int CccdReport::getDayEventNum()
{
	if( m_DayEventCounter<0 ){
		MyFile::init_day_file_counter( EVENT_COUNTER_FILE , m_DayEventCounter );
	}
	m_DayEventCounter++;
	MyFile::save_day_file_counter( EVENT_COUNTER_FILE , m_DayEventCounter );
	return m_DayEventCounter;
}

int CccdReport::GetEvtNo( const char* szEvtID )
{
	if( szEvtID && szEvtID[0] ){
		int tmp = atol( szEvtID+6 ); // skip 040802 , 6 bytes, means 
		return tmp;
	}
	return 0;
}

void CccdReport::GetEvtID( mystring& szEvtID )
{
	int day_num = getDayEventNum();
	mystring szNight;
	get_short_night_date_local( m_Time, szNight );
	char szID[20];
	sprintf(szID,"%s%.5d",szNight.c_str(),day_num);
	szEvtID = szID;
}


int CCDEventList::remove( int pos_start, int pos_end )
{
	int tab_size = size();
	if( pos_start==0 && pos_end==(tab_size-1) ){
		clear();		
		return tab_size;
	}
	

	int rm_pos=pos_start;
	for(int i=(pos_end+1);i<tab_size;i++){
		(*this)[rm_pos] = (*this)[i];
		rm_pos++;
	}
		
	int p=0;
	CCDEventList::iterator it;
	for(it=begin();it!=end();it++,p++){
		if( p==rm_pos ){
			break;
		}	
	}
	if( p==rm_pos ){
		erase( it, end() );
	}
	return (tab_size-rm_pos-1);
}

/*int CCDEventList::remove_older( int day_frame )
{
	CCDEventList::iterator i,pStart,pEnd;

	BOOL_T bStart=FALSE;
	BOOL_T bEnd=FALSE;
	for(i=begin();i!=end();i++){
		if( i->m_DayFrameIndex<=day_frame ){
			if( !bStart ){
				pStart = i;
				bStart=TRUE;
			}
			pEnd = i;
			bEnd = TRUE;
		}		
	}
	if( bStart && bEnd ){
		printf_now2("Removing events older then frame : %d ...",day_frame);
		erase( pStart, pEnd+1 );
		my_printf_now("OK\n");
	}

	return 1;
}*/

int CCDEventList::remove_older( int day_frame, CCDEventList* pRemoved )
{
	CCDEventList tmp_list;

	for(int i=0;i<size();i++){
		if( (*this)[i].m_DayFrameIndex>day_frame ){
			tmp_list.push_back( (*this)[i] );
		}else{
			if( pRemoved ){
				pRemoved->push_back( (*this)[i] );
			}
		}
	}
	(*this) = tmp_list;

	return size();
}
