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
#ifndef _CCD_DEFINES_H__
#define _CCD_DEFINES_H__

#include <mytypes.h>
#include <ccd_common_defines.h>
#include <tab2D.h>
#include <myfastlongtab.h>
#include <ccd_piman_interf_defs.h>
#include <mystring.h>
#include <time.h>
#include <mypoints.h>

class CCDMatrix;
class CCDProcState;


#define CURRRENT_EVENTS_DESC "!!! NEW EVENTS IDENTIFIED !!!\n"
#define SAT_NAME_LEN 6

// event types :
#define EVENT_TYPE_FLASH       'F'
#define EVENT_TYPE_SAT         'S'
#define EVENT_TYPE_HOT_PIXEL   'H'
#define EVENT_TYPE_BLACK_PIXEL 'B'
#define EVENT_TYPE_STAR        'G'
#define EVENT_TYPE_GRB_IN_DB   'D'


//exit codes :
#define EXIT_DUE_TO_STAMP_FILE 5


// output subdirectories names :

// cfg files :
#define DYNAMIC_CFG_FILE         "dynamic.cfg"

// fits files :
#define FITS_FILE_FORMAT         "k2%c_%d_%d.fitc"
#define FITS_FILE_BASE_FORMAT    "k2%c_%d_%d"
#define FITS_FILE_FULL_PATH_LIST "full_frames_list_ccd"
#define NEW_LIST_FILE            "new_frames_list_ccd"
#define NORMAL_FRAMES_LIST_FILE  "normal_frames_list_ccd"

// log file names :
#define TRIGGER_EVENTS_LOG     "triggerevents%d_id%d_sid%d.log"
#define FINAL_EVENTS_LOG       "finalevents_%d.log"
#define SN_ALLEVENTS_LOG          "sn_allevents_%d.log"
#define SN_EVENTS_LOG          "sn_final_%d.log"
#define SN_COIC_LOG            "sn_coic_%d.log"
#define SN_SUM_ALLEVENTS_LOG   "sn_sum_allevents_%d.log"
#define SN_SUM_EVENTS_LOG      "sn_sum_final_%d.log"
#define SN_SUM_COIC_LOG            "sn_sum_coic_%d.log"
#define TRIGGERS_SUBDIR        "Triggers"
#define TRIGGERS_PAST_SUBDIR   "TriggersPast"
#define DAQ_EXIT_FILE_STAMP    "daqExitOK.txt"
#define AVER_FRAME_EVENTS      "averframeevents_%d.log"
#define AVER_VERIF_EVENTS      "aververifevents_%d.log"
#define AVER_FINAL_EVENTS      "averfinalevents_%d.log"
#define ANTY_COIC_LOG          "antycoic_%d.log"
#define COMPARE_TO_OLD_EVENTS  "comptooldevents_%d.log"
#define TRACKS_SUBDIR          "Tracks"
#define TRACKS_LOGFILE         "newtrackslog_ccd0.txt"
#define TRACKS_LOGFILE_RADEC   "newtrackslog_radec_ccd0.txt"

// SLT files names :
#define SLT_LAST_EVTNO         "sltlastevtno.txt"
#define SLT_LOG_FILE_NAME      "slt_"
#define SLT_REJ_LOG_FILE_NAME  "slt_rejected_"
#define SLT_ACC_LOG_FILE_NAME  "slt_accepted_"

// event parts file names :
#define EVENT_PART_SUBDIR    "%s/Frame%.5d/Cam%d"
#define EVENT_PART_LIST_FILE "%s/Frame%.5d/Cam%d/list%d"
#define EVENT_PART_FILE_NAME "eventframe%d_camera%d_eventno%d_currframe%d"
#define SINGLE_PART_FILE_NAME "eventframe%d_camera%d_eventno%d_currframe%d_%d"
#define EVENT_AVER_FILE_NAME "avg_evt%d_fr%dto%d"


// system state log :
#define DAQ_STATE_LOG_FILE "daqstat.cfg"

// defines :
#define MAX_NEXT_FRAMES 10
#define MAX_PREV_SUMS 10
#define MAX_ADDITIONAL_CHAR 50
#define MAX_PIPELINE_SIZE 20
#define CCD_ASTROMETRY_FILE "ccd_astrometry_pipeline%d.cfg"

// macros :
#ifdef _CHECK_RATES_
#define INC_RATE_COUNTER( bGen, nBackgr, nSample ) { if( bGen ){ nSample++; }else{ nBackgr++; } }
#else
#define INC_RATE_COUNTER( bGen, nBackgr, nSample )
#endif

// enum eTrackCheckType_T { eNormalTrack=0, ePlaneTrack, eSingleCamTrack, 
//								 eRejIfMoreTrack, eTrackOnSumedFrame };

// enumerators :
enum eDAQMode_T { eDAQNormalMode=0,      // normal nalysing mode 
						eDAQSatTriggerMode,    // trigger mode
						eDAQTriggerMovingMode, // moving to trigger  
					   eDAQWaitingMode,       // analysis stoped 
					   eDAQCollectMode,       // only collecting frames
					   eDAQDarkCollectionMode // dark collection mode 
					 };
					 
enum eCamState_T {
							eCamState_OK=0, // OK
							eDoCleanMeasureErr,
							eDoMeasureErr,  // DoMeasureError
							eDoCleanErr,	 // DoCleanError
							eGetDataErr,    // GetDataError
							eDoStartFrameErr,
							eDoReadChipToRAMOnCamErr,
							eRefreshFullStatusErr,
							eCriticalErr, // cannot continue 
							eOtherError
					  };					 
					  
const char* GetCamStatusDesc( eCamState_T stat );					  

BOOL_T IsCollectionMode( eDAQMode_T mode );
					  
/*                nothing  final   brighten   background */					   
enum eEventType { eNone=0, eFlash, eBrighten, eBkg };
extern const char* szEventTypeDescTab[4];


eEventType GetEventType( const char* szDescType );


// functions :
mystring GetWorkingModeDesc( eDAQMode_T mode, BOOL_T bScanMode=FALSE );


enum eDiffCheckType_T { DiffPerPixel=0, DiffTotal };
enum eBRIGHT_CHECK_TYPE { eADUTotalIncrease=0, ePercentIncrease };
enum eBackgrSubtrType_T { eNoneBackgrSubtr=0, eMostPopularValue, eMedianValue };
enum ePROCSTATE_T { eNormalState=0, eTriggerState };

//enum eLaplaceType_T { eSinglePoint=0, eTwoPoints, eFourPoints, eFivePoints, eFivePlusFourMin,
//                      eNineEightFive, eNineEightSevenVeryBig,
//                      eEightFour, eEightTen, eFiveEight };


eDAQMode_T GetNewDAQMode( eCCDRequestType_T req_mode );

// defines 
#define GENERATED "GENERATED"
#define IDENTIFIED "IDENTIFIED"
#define BACKGROUND "BACKGROUND"
#define MAGNITUDE  "MAGNITUDE"

// output file items :
#define SEP "\t"
#define RECORD_SEPARATOR "-----\n"
#define TOTAL_SEPARATOR  "*****\n"
#define END_SEPARATOR    "=====\n"
#define PARAMS_SEPARATOR "+++++\n"
#define MAG_SEPARATOR    "NEXT_MAGNITUDE\n"
#define VARY_PARAM_SECTION "PARAM_VARIATION_DESCRIPTION"
#define VARY_PARAM_DESC_COUNT 5


#define MAX_CLUSTER_SIZE 100000
// #define MAX_NEIGHB_ALLOWED 100
#define MAX_EVENTS_TO_DUMP 100
#define MIN_FRAME_SIZE 5
#define MAX_NEIGHB_SIZE 100


inline BOOL_T ADD_POINT_func( LONG_T* table, LONG_T& idx, int pos, int max_size=MAX_CLUSTER_SIZE )
{
	if(idx>=max_size){
		// throw CCDClusterSizeExcp();
	   return FALSE;
	}
	table[idx] = pos;
	idx++;
	return TRUE;	                        
}


// in optimized version it will be macro :
#define ADD_POINT(table,idx,pos) { \
if(idx>=MAX_CLUSTER_SIZE){ \
	printf("Cluster size exceeded (idx=%d), exiting\n",idx);\
	exit(-1); \
} table[idx] = pos;idx++;}

	// for(int i=0;i<idx;i++){ \
	// 	printf("(%d),",table[i]); \
	//} \

#define _CCDLIB_PRINTF_ if(gCCDParams.m_bStdoutOn)printf


class CCDPipeline;
class CPixelList;
class CCcdCfg;

struct DxDyFromCurrentFrame 
{
	double frameDX;
	double frameDY;
	DxDyFromCurrentFrame():frameDX(0),frameDY(0){}	
};

struct CPixelDistrValues
{
	double mean;
	double sigma;
};

class CFrameIdentStat
{
public :
	CFrameIdentStat( int frame_index );
	~CFrameIdentStat(){};
	void Init();
	void LogToFile( CCDPipeline* pPipeline, const char* fname="event_rates",
						 const char* szSubDir="Events" );
	void DecTrackAndCoicEvents();
	
	// accepted after :
	int m_FrameIndex;
	
	int nTotal;
	int nMaxAllowedValue;
	int nTnewCut;
	int nTprevCut;
	
	// rejected after :
	int nShapeCut;
	int nBlackCut;
	int nHotCut;
	int nLocalMaxCut;
	int nMinPrevCut;	
	BOOL_T bAllRejected;
	int nIfMoreThen;
	int nOverlap;
	int nAfterCoic;
	int nTracks;
	BOOL_T bAllRejAfterCoic;
	int nSatRej;
	int nStarRej;
	int nFinalEvents;
	int nClusterOverlap;
};

class CFrameEventStatTab : public vector<CFrameIdentStat>
{
public :
	CFrameEventStatTab();
	~CFrameEventStatTab();

	int IncFinal( int frame_index );	
	void DecTrackEvents( int frame_index );
	void DecTrackAndCoicEvents( int frame_index );
};

class CPixelAnalyseIn
{
public :
	CPixelAnalyseIn();
	void SetPrevMatrix( const CPixelAnalyseIn& in );		

	LONG_T x;
	LONG_T y; 
	LONG_T pos;
   LONG_T xSize;
   LONG_T ySize;
   LONG_T Size;
	CCDMatrix* Matrix;
	ELEM_TYPE* p_data;
	ELEM_TYPE** p_data_fast;		

	BIG_ELEM_TYPE*  p_curr_data_laplace_normal;
	BIG_ELEM_TYPE** p_curr_data_laplace;

	BIG_ELEM_TYPE* p_homeo_data;
	BIG_ELEM_TYPE** p_fast_homeo_data;

	BIG_ELEM_TYPE* p_laplace_data;
	BIG_ELEM_TYPE** p_laplace_data_fast;		
	
	BIG_ELEM_TYPE** p_prev_lap_fast;

	CCDPipeline* pPipeline;
	LONG_T pipeline_size_minus_1;

//	Table2D<ELEM_TYPE>* PrevMatrixPtr[MAX_PIPELINE_SIZE];		
	CCDMatrix* PrevMatrixPtr[MAX_PIPELINE_SIZE];
	int PrevMatrixPtrCnt;
	DxDyFromCurrentFrame PrevFramesShiftTab[MAX_PIPELINE_SIZE];	
	double PrevFramesTime[MAX_PIPELINE_SIZE];
	CFastLongTable PrevValues;

	BIG_ELEM_TYPE** PrevLaplacePtr[MAX_PIPELINE_SIZE];

	Table2D<BIG_ELEM_TYPE>* pMaxLaplacePrev;
	BIG_ELEM_TYPE** p_max_prev_lap;

	LONG_T ccd_index;
	LONG_T frame_index;
	LONG_T day_frame_index;

	CPixelList* pPixelList;
	
	int PrevFramesX[MAX_PIPELINE_SIZE];
	int PrevFramesY[MAX_PIPELINE_SIZE];		
	
	CCDProcState* pCCDInfo;
	CCcdCfg* pCamCfg;
	

	//---------------- ALGO - SECTION -------------------------------
	// some important values are stored for whole frame analysis	
	int nAver; // numer of averaged frames 

	//-----------------------------------------------------------
	// additional parameters passed - to avoid re-calculations :
	double SigmaLap;
	double MeanLap;
	double SigmaRaw;
	double MeanRaw;
	long treshold_for_max;
	long treshold_for_not_more_then_n_exceeds;

	// treshold to calculate raw cluster :
	int treshold_for_raw_cluster;
	int treshold_for_raw_cluster_aver;
	                                                                                               
	double treshold_for_cluster;
	double treshold_for_black_pixel;
	double treshold_for_prev;

	int treshold_for_hot;		
	//---------------------------------------------------------------
	
	
	CPixelDistrValues m_pDistrValues[MAX_LAPLACE_DEFINED];
	
	//---------------------- END OF ALGO - SECTION -------------------------
};
	
//----------------------------------------------
// Adding new value :
//
// add fields in CSinglePixelData and CPixelAnalyseOut, fill value in CPixelAnalyseOut
// and program copies it to CCcdReport structure also 
//----------------------------------------------
class CSinglePixelData
{
public :
	CSinglePixelData();
	void Init();

	LONG_T PixelRawValue;
	int    MaxPixelRawValue;
	
	LONG_T newSum;

	// raw data :
	LONG_T newRawSum;
	LONG_T maxSum;
	LONG_T homeoSum;
	
	// laplace data :
	BOOL_T m_bLaplaceMedianRejected;
	
	LONG_T laplaceSum;
	LONG_T medianLaplaceSum;
	LONG_T prevLaplaceSum;
	LONG_T prevLaplaceSumX;		
	LONG_T prevLaplaceSumY;		
	LONG_T prevLaplaceSumMin;
	LONG_T homeoLaplaceSum;
	int m_PrevLapInPixel;
	int laplaceOnNext;

	// average of previous :
	LONG_T maxAverageOfPrev;
	
	LONG_T otherSum;
	
	double x0_real;
	double y0_real;
	LONG_T x0;
	LONG_T y0;
	LONG_T pos0;
	
	eEventType eventType;
	
	// other checks results :
	BOOL_T m_bLocalMaxRejected;		

	// outside analysed region :
	BOOL_T m_bOutSideAcceptanceRegion;

	// check if not more then Tresh on more then N frames :
	BOOL_T m_bNotMoreThenNAboveRejected;
	LONG_T m_nCountAbove;
	LONG_T treshold_for_not_more_then_n_exceeds;
	
	// confirming on next frames :
	LONG_T LaplaceOnNextFrames[MAX_NEXT_FRAMES];
	BOOL_T m_bRejectedByNextFrames;		
	
	// custom event analysis (performed outside class CCD_Analyser and not 
	// affecting code of this class ) :
	BOOL_T m_bRejectedByCustomEventAnal;


	// event shape analysis :
	BOOL_T m_bRejectedByEventShape;
	double m_Sphericity;	
	double max_redial;
	double max_noise_level;

	// black pixel ratio :
	double m_fBlackRatio;	
	
	// track cut :
	BOOL_T m_bRejectedDueToTrack;
	eTrackCheckType_T m_TrackType;		
	BOOL_T m_bRejByTrackOnSingleCam;

	// coicydence rejected :
	BOOL_T m_bRejectedByCoic;
	
	// found in trigger mode :
	BOOL_T m_bInTriggerMode;
	
	// time_t :
	// time_t m_Time;		

	// average of :
	// double PrevFramesMaxVal[MAX_PIPELINE_SIZE];
	int PrevFramesX[MAX_PIPELINE_SIZE];
	int PrevFramesY[MAX_PIPELINE_SIZE];

	int PrevLaplaceValuesCount;
	int PrevLaplaceValues[MAX_PIPELINE_SIZE];

	// event significance and likehood :
	double m_Likehood;
	double m_Significance;

	CPoint m_MaxClusterPixel;
	int m_MaxClusterValue;

	int laplacePlusValues[MAX_PIXELS_IN_LAPLACE];
	int laplaceMinusValues[MAX_PIXELS_IN_LAPLACE];
	int laplacePlusCount;
	int laplaceMinusCount;
	
	// 
	double minChi2; // minimal distance to existing track
	double minChi2On3; // minimal value of chi2 on 3 points including current one
 	double minDist; // minimal distance to track
	int bestTrackID;
	int veloCheckOK;
	double rx;
	double ry;
};

class CPixelAnalyseOut
{
public:
	CPixelAnalyseOut();
	CPixelAnalyseOut( const CPixelAnalyseOut& right );
	~CPixelAnalyseOut();
	CPixelAnalyseOut& operator=( const CPixelAnalyseOut& right );
	void Dump();
	
	LONG_T neighbours[MAX_NEIGHB_SIZE];
	LONG_T neighb_count;
	
	LONG_T* neighb_list;
   LONG_T ncnt;
   CPixelList* pixel_list;

	LONG_T newClusterSum;
   LONG_T prevClusterSum[MAX_PREV_SUMS];
	LONG_T prevFramesCount;

	// cluster :
	LONG_T* cluster;
	LONG_T cluster_cnt;

	// cluster with surrounding :
	LONG_T* ClusterWithMore;
	LONG_T ClusterWithMoreCnt;

	// hit pixel list - to perform cleaning :
	LONG_T hitpixel_count;
	LONG_T* hitpixel_list;

	// list of start points for cluster - in fact stack in algorithm of
	// finding cluster :
	LONG_T startPoints_cnt;
	LONG_T* startPoints;

	Table2D<BIG_ELEM_TYPE>* m_pCluster;	

	CSinglePixelData m_PixelOut;
};   


class CCD_Analyser;

class CVariableInfo
{
public:
	CVariableInfo();
	void Init();
	
	double l_new;
	double l_prev;
	int nAfterTv;
	int nOverlaps;				
	double fSpericity;
	double fBlackRatio;
	double averInPixel;
};

#endif



