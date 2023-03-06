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
#ifndef _CCD_REPORT_H__
#define _CCD_REPORT_H__

#include <mytypes.h>
#include <mypoints.h>
#include <vector>
#include <basestructs.h>
#include "ccd_defines.h"

using namespace std;

class CGRBInfo;
class mystring;
class CCDPipeline;
class CCDEventList;
class CccdReport;
class CCDMatrix;
class CCDProcState;
struct satInfo;

class CAstroCoord
{
public :
	CAstroCoord(): m_Azimuth(0), m_Altitude(0), Dec(0), RA(0), HA(0) {};
	
	CAstroCoord& operator=(const CAstroCoord& right);

	double m_Azimuth;
	double m_Altitude;
	
	double Dec; // in radians 
	double RA;  // in radians 
	double HA;
};


class CEventBaseInfo
{
public:
	CEventBaseInfo( int max_point_x=0, int max_point_y=0, int point_x=0, 
						 int point_y=0, int frame_index=0, 
						 int day_frame_index=0 );
	CEventBaseInfo( const CEventBaseInfo& right );
	CEventBaseInfo& operator=( const CEventBaseInfo& right );						 
	
	int GetMemSize( BOOL_T bShowAll=FALSE );		
	
	void SetTransformedPoint( double x, double y, BOOL_T bCalcCoicRadius=FALSE );

   CPoint m_MaxPoint;
   CPoint m_Point;
   CPoint m_PointTransformed;
   double m_CoicRadius;
   double m_CoicRadiusInRad;
	int m_FrameIndex;   	
	LONG_T m_DayFrameIndex;
	int m_Time;
	
	// official coordinates :
	CAstroCoord m_AstroCoord; // coordinates of current event

	BOOL_T m_bSavedToDB;
};

class CTrackDesc
{
public:
	static mystring m_LogHeader;
	static mystring m_LogFmt;

   double a;
   double b;
	int m_cam_idx;
	int track_id;
	double vx;
	double vy;
	int frame_index;
	int tr_night;
	int frame_index_end;

	// additional info
	double chi2;
	BOOL_T bMoveOK;

	// RA,DEC track parameters :
	double radec_a;
	double radec_b;
	double v_ra;
	double v_dec;
	BOOL_T m_bRADEC;
	
	static int m_TrackIDGenerator;
	vector<CEventBaseInfo> m_EventsOnTrack;	


	int GetMemSize( BOOL_T bShowAll=FALSE );

	CTrackDesc( double aa=0, double bb=0, int cam_idx=0 );
	CTrackDesc( const CTrackDesc& right );
	virtual ~CTrackDesc();

	// returns TRUE if event belongs to track :
	BOOL_T CheckEvent_RADEC( double ra, double dec,
	                         double& chi2, double chi2_limit );
	
	BOOL_T CheckEvent( int x, int y, double& chi2, 
							 double chi2_limit=100000.00 );
	BOOL_T CheckMove( int x, int y, int frame );
	BOOL_T CheckMove_RADEC(  int x, int y,
	                         double ra, double dec, int frame );
	BOOL_T DoesEventBelongs( int x, int y, int frame, double r=2.00 );
	
	void SetTrackDesc( double aa, double bb, int cam_idx );

	BOOL_T operator==( const CTrackDesc& right );      
	CTrackDesc& operator=( const CTrackDesc& right );
	
	int get_xy_values( double* x_values, double* y_values ) const;
	int get_radec_values( double* ra_values, double* dec_values ) const;		
	
	double FindMinDist( double x, double y );
	double FindMinDistByMaxPixel( double x, double y );
	
	void GetEventsList( mystring& szList );		
	void GetEventsList_RADEC( mystring& szList, BOOL_T bDoConvert=TRUE );
	void GetSortedEventsList( mystring& szList );
	void SortEventListByTime();
	void calcAverageVelocity( double& vx, double& vy );
	void calcAverageVelocity();
	void UpdateTrack( double new_a, double new_b, CccdReport& new_evt );

/*	inline void calcAverageVelocity(){
		calcAverageVelocity( vx, vy );			
	}*/
	
	BOOL_T CheckVelocity( CEventBaseInfo* pNewEvent, double velError, 
								 int idx );
	BOOL_T CheckVelocity( CEventBaseInfo* pEvt1, CEventBaseInfo* pEvt2, double fVeloError,
								 int idx	);

	void LogEventCheck( CccdReport& evt, int cam_idx, double chi2, 
							  double vx, double vy,
							  double rx, double ry, BOOL_T bOK,
							  double minDist,
							  const char* szReason );

	void FlagEventsSaved();							  							  
							  
	void get_frame_range( int& min_track, int& max_track );

	// parsing track log file :
	BOOL_T Parse( const char* szTrack, int len );
	BOOL_T ParseWithRADEC( const char* szTrack, int len );
};

class CTrackList : public vector<CTrackDesc>
{
public :
	int GetMemSize( BOOL_T bShowAll=FALSE );

	CTrackList(){};
	virtual ~CTrackList();

	// returns pointer to best track 
	CTrackDesc* CheckEvent( int x, int y, double& min_chi2, 
									double chi2_limit, BOOL_T bDoCheck=FALSE );	
	CTrackDesc* CheckEvent_RADEC( int x, int y, double ra, double dec,
												double& min_chi2, 
												double chi2_limit, BOOL_T bDoCheck );
									
	int GetBestTracks( int x, int y, int frame,
							 int n, CTrackList& bestlist );		
	int GetBestTracks_RADEC( double ra, double dec,
										 int x, int y, int frame,
										 int n, CTrackList& bestlist );
							 
		
	void Add( double a, double b, int cam_idx );
	CTrackDesc* Find( const CTrackDesc& newelem );
	void Add( const CTrackDesc& newelem );
	CTrackDesc* Find( int trackID );
	
	static CTrackDesc* Find( vector<CTrackDesc*>& trackList, int trackID );

	// removes several 
	void remove( int min_pos, int max_pos );
	
	// reading log file :
	int Read( const char* szFileName, BOOL_T bRADEC=FALSE, 
				 const char* szVerifFile="", BOOL_T bSaveRADEC=FALSE );
	int ReadWithRADEC( const char* szFileName );				 
	int MergeRADECInfo( CCDEventList& event_list, BOOL_T bSaveRADEC=FALSE );
};



class CccdReport : public CEventBaseInfo
{
public :
	// error handling when reading files :
	static BOOL_T m_bIgnoreReadError;
	
	
	static int m_DayEventCounter;
	static int getDayEventNum();
	void GetEvtID( mystring& szEvtID );
	static int GetEvtNo( const char* szEvtID );
	int GetEvtNo()
		{ return GetEvtNo( m_EvtID.c_str() ); }

	void AddGRBInfo( CGRBInfo* pGrbInfo );
	static double CalcDist( CccdReport& evt1, CccdReport& evt2 );

	int GetMemSize( BOOL_T bShowAll=FALSE );

	enum eAddInfoType_T { eClusterChar=0, eCrossChar, eCrossAboveTreshChar };
	static const char* GetAddInfoDesc( eAddInfoType_T add_info_type );

	// internal trigger ID
	void GetInternalTriggerID( mystring& szID );

	static mystring m_EventDescOutFmt;
	static mystring m_EventDescInFmt;
	static mystring m_EventDescHeader;
	static mystring m_FinalEventHeader;
	static mystring m_FinalEventOutFmt;
	static mystring m_FinalEventInFmt;
	
	// grb info concatenated :
	static mystring m_GrbInfoHeader;
	static mystring m_GrbInfoFmt;
	

	static BOOL_T IsGoodToSave( CccdReport& evt );	

	static BOOL_T SaveGenEvent( const char* fname,
										 int x, int y,
	                            const char* mag,
 	                            Table2D<ELEM_SAMPLE_TYPE>& OtherImage, 
 	                            int frame_index, int max_lap=0);
	
	void Init();
	
	CccdReport();	
	CccdReport( LONG_T camNum, const CPixelAnalyseIn& in, 
	            const CPixelAnalyseOut& out, 
	            BOOL_T bGen,BOOL_T bIdentified,	             
	            CPoint* pLowLeft=NULL,CPoint* pTopRight=NULL,
	            const char* mag=NULL,BOOL_T bFirstLevelAccepted=FALSE );

	CccdReport(const CccdReport& right);
	
	~CccdReport();
	
	CccdReport& operator=(const CccdReport& right);

	BOOL_T IsIdentified() const;
	BOOL_T IsIdentifiedForTrack() const;

	void CopyCluster( const CccdReport& right);
	
	void SetCluster( const LONG_T* cluster, LONG_T cluster_cnt );

	void CalcCoicRadiusInRad( double ra, double dec );

	void SetIdentificationData( const CccdReport& right);
		
	void SetIdentificationData( const CPixelAnalyseIn& in, const CPixelAnalyseOut& out );

	void SetPrevClusterSum( LONG_T prev_frames_cnt, const LONG_T* _prevClusterSum );     	

	void SetNextFramesSum( LONG_T next_frames_cnt, const LONG_T* _nextClusterSum );

	static void GetDescHeader( mystring& szHeader );
	void GetEventTypeDesc( char* szOutType );

	// parsing output lines :
	void ParseEventType( const char* szEventTypeDesc );
	void ParseRA_DEC( const char* szRA_DEC );
	void ParseTrackDesc( const char* szTrack );

	void GetEventDesc( mystring& desc );	
	long GetEventDescStr( char* szLine );
	BOOL_T ParseOutputLine( const char* szLine );		

	// detailed desc :
	void GetDetailEventDesc( mystring& szEventDesc, long Index, long xSize, CCDPipeline* pPipeline );

	void GetXYRange( LONG_T& x_size, LONG_T& y_size,
	                 LONG_T xSize, LONG_T ySize );

	void RejectByNextFrames();
	
	void WriteEventDetails( CCDPipeline* pPipeline, int eventNo );

	// calculation of astro coordinates :
	void CalcAstroCoordinates( CCDMatrix& frame, CCDProcState* pFrameInfo );			

	//
	static double FindMinDist( vector<CccdReport*>& eventList, CccdReport& evt );
	
	// fields :

	// satellite information :
	struct satInfo* m_SatInfo;
	double m_MinDistOutCone; // [km]

	int m_PipelineIndex;
	
	LONG_T m_CameraIdx; // obsolate do not use 
	double m_MaxValue;
	CPoint m_ClusterCenter;
	CPoint m_GenPoint;
	BOOL_T m_bGenerated;
	BOOL_T m_bIdentified;
	CPoint m_LowLeft;
	CPoint m_TopRight;
	mystring m_Magnitude;
	mystring m_szGRBInfo;
	LONG_T* m_Cluster;
	LONG_T m_ClusterCount;
	double m_AdditionalInfo[MAX_ADDITIONAL_CHAR];
	LONG_T m_nAddInfo;
	
	// identification description :
	BOOL_T m_bFirstLevelAccepted;

	LONG_T newClusterSum;
	LONG_T prevClusterSum[MAX_PREV_SUMS];
	LONG_T     prevFramesCount;

	LONG_T nextClusterSum[MAX_PREV_SUMS];
	LONG_T     nextFramesCount;
	BOOL_T     m_bCheckOnNext;

	// dumpng frames with events to FITS file :
	LONG_T m_LastDumpedFrame;

	CSinglePixelData m_PixelAnalResults;
	
	
	// processing flags :
	BOOL_T m_bRemoveFromList;

	// official coordinates :
	// CAstroCoord m_AstroCoord; 
	// coordinates of current event
	
	
	// program suggestion what event can be :
	char m_EventType;
	mystring m_szSatName;
	char m_SatType;
	int EvtIdx;
	double m_SatCoicRadiusInSec;	
	char visibility;
	
	// ID of internal trigger :
	mystring m_EvtID;
	
	// clouds info :
	int m_nFrameStarCount;
	
	// external triggers info :
	mystring m_szExternal;

};


// struct event_finder : public unary_function<double, void>
struct event_finder
    {
          event_finder(){}
          bool operator()( const CccdReport& left, const CccdReport& right ){
          	return (left.m_DayFrameIndex < right.m_DayFrameIndex);
          }
    };
                                   

class CCDEventList : public vector<CccdReport>
{
protected:
	static BOOL_T m_bInitialized;
	static void Initialize();
public:
	static eLogFileFormat m_DefaultFormat;
	eLogFileFormat format;

	CCDEventList(int nAutoAlloc=0);
	CCDEventList( const CCDEventList& right );		
	void Add( const CccdReport& newEvent );
	CCDEventList& operator+=( const CCDEventList& right );	
	CCDEventList& operator=( const CCDEventList& right );		
	CCDEventList& AddNotFound( const CCDEventList& genEvents );

	BOOL_T WasAnyIdentified();
	BOOL_T IsAnyIdentifiedNotGen();
	int GetNotIndentifiedCount();

	static mystring GetOutputDir();
	static mystring GetEventFileName( CccdReport& event, int FrameIndex, int EventNo );

	void DumpAllEvents();	
	static int DumpEventReport( const char* fname , CCDEventList& ccdReport,
	                            LONG_T idx=-1, BOOL_T bStdout=FALSE, 
	                            BOOL_T bToFile=TRUE, BOOL_T bNew=TRUE,
	                            BOOL_T bOnlyIdentified=FALSE, 
	                            const char* szComment=CURRRENT_EVENTS_DESC );
	                    
	static int SprintfEvent( char* ptr, CccdReport& evt );  
	static int SprintfGrb( char* ptr, CGRBInfo* pGrbInfo );
	static int DumpFinalEvents( const char* fname , CCDEventList& ccdReport,
										 BOOL_T bUseOutDir=TRUE );
	static int ReadFinalEvents( const char* fname , CCDEventList& ccdReport,
										 time_t startTime=0, time_t endTime=0,
										 int frameno=-1 );
										 
	static int DumpFinalEventLine( const char* fname, CccdReport& evt, CGRBInfo* pGrbInfo=NULL );
											 
	                
	BOOL_T CheckSortByFrame();	                
	CccdReport* Find( int frameno , int x, int y );	                	             
	CccdReport* FindFromBack( int frameno , int x, int y );
	CccdReport* FindEvent( int x, int y, int redial=0 );
	CccdReport* FindEvent( int x, int y, int frameno, int redial=0 );
	CccdReport* FindEventByMaxPoint( int x, int y, int redial=0 );
	
	CccdReport* FindEventByRaDec( double ra, double dec, double radius );
	
	BOOL_T RejectEventDueToTrack( int x, int y );

	// new function - automatically checks which format should be used
	// when reading and calls ReadEvents or ReadFinalEvents
	int Read( const char* fname, int frameno=-1 );
	
	static void Dump( CccdReport& evt, eLogFileFormat _format, const char* szFileName=NULL );
	void Dump( const char* szFileName=NULL );
	
	int ReadEvents( const char* fname, BOOL_T bExcp=TRUE , int min_frame=-1,
						 int frameno=-1 );
	
	int GetEventsByFrameIndex( CCDEventList& out, int frame_index, int end_frame=-1 );

	CccdReport* FindByFrameIndex( int frame_index, int& pos );
	
	int GetNotMatchedEventsCountSkipGen();
	
	void SetEvtIdx();
	void SetEvtNo();
	void SetEventTime( time_t ut_time );

	int GetEvents( CCDEventList& out_list, int day_frame );
	int GetEventsSinceFrame( CCDEventList& out_list, int start_frame );
	int GetEvents( CCDEventList& out_list, int start_frame, int end_frame, BOOL_T bOnlyFinal=FALSE );
	
	int remove( int pos_start, int pos_end );
	int remove_older( int day_frame, CCDEventList* pRemoved=NULL );

	void CalcConeDist();
};

class CFrameEvents : public vector<CCDEventList>
{
public:
	int m_FrameCounter;

	CFrameEvents( int FrameCounter );
	//CFrameEvents( const CFrameEvents& right );
	//CFrameEvents& operator=( const CFrameEvents& right );
	virtual ~CFrameEvents();
	
	int GetEventsCount( int cam_idx );
	int HasFrameEvents( int cam_idx );

	int GetNotMatchedEventsCount( int cam_idx );
	int GetNotMatchedEventsCountSkipGen( int cam_idx );		


	int GetMemSize( BOOL_T bShowAll=FALSE );

	// CTrackList m_TracksOnFrame;	
};


class CMagStat
{
public :
	CMagStat();

	mystring m_Magnitude;
	LONGLONG_T m_MagGenerated;
	LONGLONG_T m_MagIdentified;		
};

class CIdentStat
{
public:
	CIdentStat();
	~CIdentStat();	
	CIdentStat( const CIdentStat& right );
	CIdentStat& operator=( const CIdentStat& right );
	
	
	CMagStat* FindMagInfo( const char* szMag );

	void AddOrUpdateMagInfo( const char* szMag, BOOL_T bIdent);
	void AddOrUpdateMagInfo( CccdReport& evt );
	
	vector<CMagStat>* GetTab() { return m_pMagnitTab; }
	void Reset();
	vector<CEnvVar>& GetParamTab();	

	vector<CMagStat>* m_pMagnitTab;
	vector<CEnvVar>*  m_pParamValues;
	LONGLONG_T m_TotalGenerated;
	LONGLONG_T m_TotalIdentified;		
	LONGLONG_T m_nBackground;
};


#endif
