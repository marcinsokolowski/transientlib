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
#ifndef _CCD_ANALYSE_H__
#define _CCD_ANALYSE_H__
#include <mytypes.h>
#include <mypoints.h>
#include <ccdcfg.h>
#include "ccd_defines.h"
#include "ccd_report.h"
#include <tab2D.h>
#include <baseanal.h>
#include <myutil.h>
#include "ccd_datares.h"
#include "ccd_procstate.h"
#include "ccd_globals.h"

class CVariableInfo;
class CCDPipeline;
class CCDMatrix;
class CPointList;
class CLongList;
class CPixelList;
class CDescTab2D;
class ShiftInfo;
class CSatList;
class CFrameIdentStat;
struct sEventDesc;


class CCD_Analyser : public CBaseAnal
{
public :
	static CLongPoint* GetNeighbPoints(){ return m_NeighbPoints; }
	
	// customization analysis - to call from other rutines :
	static BOOL_T (*m_pAnalFoundEventFunc)( void* event_info, int type, int ccd_idx, void* event_info2  );

	static void SetCustomEventAnalFunc( BOOL_T (*pAnalFoundEventFunc)( void* event_info, int type, int ccd_idx, void* event_info2 ) );

	static BOOL_T CheckAnalRange( int x, int y, int xSize, int ySize );
	static BOOL_T CheckAnalRangeAuto( int x, int y, int xSize, int ySize );
	
		
	// ACCESSORS :
	void SetPipelineObj( CCDPipeline* pPipeline ){ m_pPipeline=pPipeline; }
	
public:
	static BOOL_T m_bInitialized;

	// full list of events found with single frame analysis 
	// - sorted in order :
	// 1/ generated events 
	// 2/ fake events 
	CCDEventList m_AllEvents;	 		
	static CCDEventList m_AllRunEvents;

	// satelites list for current frame :
	CSatList* m_pCurrSatList;	// only visible sat are stored in this list
	CSatList* m_pAllSatList;   // all satelites from DB are stored here 

	// optional optimization stuctures 
	CDescTab2D* m_pNeighMap;		
	
	// map of points around to prepare aureola and avoid to much calculations :
	static CLongPoint m_Aureola[MAX_CLUSTER_SIZE];
	static LONG_T m_nPixelsInAureola;

	// map of shape to confirm :
	static CLongPoint m_singlePoint[1];

	static CLongPoint m_ShapeToConfirmMap[MAX_CLUSTER_SIZE];
	static LONG_T m_ConfShapeCount;	

	static CLongPoint m_NeighbPoints[MAX_CLUSTER_SIZE];
	static LONG_T m_SelfPos;
	static LONG_T m_FirstDy;
	
	// rejected / accepted stat :
	CFrameEventStatTab backgrStat;
	CFrameEventStatTab samplesStat;
	
	// sum analysis background stat :
	CFrameEventStatTab sumBackgrStat;

	// info about variables to be histogramed :
	CVariableInfo varInfo;
	   		

	Table2D<BIG_ELEM_TYPE>* m_pMaxLaplacePrev;
	BIG_ELEM_TYPE* m_pWrkTableInX;

	// pipeline object for this :
	CCDPipeline* m_pPipeline;
	
	// just for histograming purposes :
	eTrackCheckType_T m_eCurrFitTrackType;	
	
	// functions for creating optimization maps :
	static void DumpCluster( LONG_T* cluster, LONG_T cluster_cnt, ELEM_TYPE* p_data,
									 LONG_T xSize, LONG_T ySize, const char* desc );
	
	static void CalcAurola();

	void GetClusterFromShapeMap( LONG_T x, LONG_T y, 
										  LONG_T xSize, LONG_T ySize,
										  LONG_T* ClusterWithMore,
										  LONG_T& ClusterWithMoreCnt );

	static void GetOutsidePixels( int x, int y,
												 LONG_T* out_cluster, LONG_T& out_cnt,
												 LONG_T xSize, LONG_T ySize,
											 int r0, int r1 );

	static void GetOutsidePixels( LONG_T* cluster, LONG_T cluster_cnt,
     	                               LONG_T* out_cluster, LONG_T& out_cnt,  
	                                  LONG_T xSize, LONG_T ySize,
	                                  	 int nPixels );

	static double CalcSum( int x, int y, int radius, double sky, ELEM_TYPE** data );
	static double CalcSumWithSide( int x, int y, int radius, double sky, ELEM_TYPE** data );

	static inline int CalcClusterSum( LONG_T* cluster, LONG_T cluster_cnt,
	                           ELEM_TYPE* p_data )
	{
		register int sum=0;
		for(register int i=0;i<cluster_cnt;i++){
			sum += p_data[ cluster[i] ];
		}
		return sum;
	}
	                           
	                                  
	static int GetRawCluster( int x, int y, LONG_T xSize, LONG_T ySize,
	                          ELEM_TYPE* p_data, ELEM_TYPE** p_data_fast,
	                          LONG_T* cluster, LONG_T& cluster_cnt,
	                          int treshold );	                                                                                                   

	static int CalcClusterLaplace( int x, int y, LONG_T xSize, LONG_T ySize,
									ELEM_TYPE* p_data, ELEM_TYPE** p_data_fast,
									int treshold=-1000 );	

	static int GetClusterEdge( LONG_T* cluster, LONG_T cluster_cnt,
							LONG_T xSize, LONG_T ySize,
							LONG_T* EdgePixelsList, LONG_T& EdgePixelsCount );


	void GetClusterWithPointsAround( CLongList& cluster, CLongList& out_cluster, 
	                                 LONG_T xSize, LONG_T ySize, LONG_T nPixels );

	void GetClusterWithPointsAroundOpt( LONG_T* cluster, LONG_T cluster_cnt,
	                                    LONG_T* out_cluster, LONG_T& out_cnt,  
	                                    LONG_T xSize, LONG_T ySize, LONG_T nPixels );

	LONG_T GetPointsAround( LONG_T pos, LONG_T xSize, LONG_T ySize, 
                                CLongList& around, LONG_T r0 );

	// prepares table of matrix pointers - to have it ready to use 
	// for pixel by pixel calculations :
	Table2D<ELEM_TYPE>** GetPrevMatrixPtrs( LONG_T index, CCDPipeline& ccd_pipeline );	                                 		
public :
	// calculate astro coordinates of new events :	
	void CalcCoordForEvents( CCDPipeline& ccd_pipeline, CCDProcState* pFrameInfo,
									 CCDEventList& events, CCDMatrix& matrix, double sigmaB );
	void CalcAdditionalInfoForEvents( CCDPipeline& ccd_pipeline );
	void CalcAdditionalInfoForEvents( CCDPipeline& ccd_pipeline, CCDEventList& events );
	void CalcAstroCoorinatesOfNewEvents( CCDPipeline& ccd_pipeline );

	// coicydence of 2 ccd cameras :
	static int FillAntyCoicList( CCDEventList& events, CCDEventList& coicevents, 
				     CCDEventList& antycoic );

	static int VerifySingleCameraEvents( CCDPipeline& pipeline1, CCDEventList& eventsOn1,
					     CCDEventList& verifiedOn1  );

	static BOOL_T IsSingleFrameEvent( CccdReport& evt1, CccdReport& evt2,
	              					       CCDPipeline& pipeline1, CCDPipeline& pipeline2 );

	BOOL_T IsSingleFrameTrack( CccdReport& evt1, CCDPipeline& pipeline1 );

	static int VerifyCoicydence( CCDPipeline& pipeline1, CCDPipeline& pipeline2 );

	static int VerifyCoicydence( CCDPipeline& pipeline1, CCDPipeline& pipeline2,
				CCDEventList& eventsOn1, CCDEventList& eventsOn2,
				CCDEventList& verifiedOn1, CCDEventList& verifiedOn2,
				int& nSampleEvents, int& confirmedEvents, int i=0,
				BOOL_T bDoCheckStars=TRUE, const char* szType="NORMAL" );


	static int VerifyCoicydenceOfEvents( CCDEventList& eventsOn1, CCDEventList& eventsOn2,
	                             CCDEventList& verifiedOn1, CCDEventList& verifiedOn2,
	                             int max_time_diff, double coic_redial );
	                                                                        

	static int CombineCamAnalResults( CCDEventList& cam1events, CCDEventList& cam2events );

	// function to state if verified event was accepted :
	static BOOL_T IsVerifiedOK( CccdReport& event );
	static BOOL_T VerifyEvent( CccdReport& event );

	BOOL_T CheckBrightening( const CPixelAnalyseIn& in,
            CPixelList& pixel_list,CCDMatrix& Matrix,
				LONG_T newSum, LONG_T maxSum,
				LONG_T homeoSum,LONG_T otherSum,
				LONG_T& newClusterSum,LONG_T* prevClusterSum,
				LONG_T x0, LONG_T y0,LONG_T pos0,
				LONG_T* cluster, LONG_T& cluster_cnt,
				LONG_T* ClusterWithMore, LONG_T& ClusterWithMoreCnt );


	void CalcTresholdsByMap( CPixelAnalyseIn& in, int& nNewLaplaceTreshold, 
									 int& nMaxLaplaceOnOther );
	void CalculateTresholds( CPixelAnalyseIn& in );
	void CalcMaxLaplaceTab( const CPixelAnalyseIn& in );
	long GetMaxLaplaceOnPrevRotCorrect( const CPixelAnalyseIn& in );

	BOOL_T FullAnalysePixel( CPixelAnalyseIn& in,
	                         CPixelAnalyseOut& out,
					  			    BOOL_T bTrace=FALSE);


	BOOL_T AnalysePixel( CPixelAnalyseIn& in, CPixelAnalyseOut& out, BOOL_T bTrace=FALSE );

	BOOL_T AnalysePixelOpt3( const CPixelAnalyseIn& in, CPixelAnalyseOut& out, BOOL_T bTrace=FALSE );



	BOOL_T AnalysePixel( LONG_T x, LONG_T y, CCDMatrix& Matrix, CCDPipeline* pPipeline, LONG_T ccd_index,
								mystring& szOutput );



	// static analysing functions go here :
	static double CalcClusterRatio( ELEM_TYPE* p_data, LONG_T pos,
	                                LONG_T* pixel_list, LONG_T list_count, LONG_T& pixel_cnt,
	                                BOOL_T bAboveTresholdOnly=FALSE, double Treshold=0, BOOL_T bUseMax=FALSE);
	
	
	static LONG_T CalcShapePoints( CLongPoint* shapePoints,
	       		                   eConfShape_T shapeType,
                                  double Redial );
	
	static void CalcShapeToConfirm();


	CCD_Analyser( int sizeX, int sizeY );
	~CCD_Analyser();

	// parameters initialization/refreshing 
	static void Initialize();
	static void RefreshParams();
	   
	

	static void VerifyPrevValue(LONG_T& val);
	static void VerifyValue(LONG_T& val,BOOL_T bMax=FALSE);

	void Init();
	void ClearState();

	CDescTab2D* GetNeighbMap(){ return m_pNeighMap; }

	void CalcCenterOfHit( ELEM_TYPE* p_data, CPointList& cluster,
	                      long xSize, long& x0, long& y0 );

	void CalcCenterOfHit( ELEM_TYPE* p_data, CLongList& cluster,
	                      long xSize, long& x0, long& y0 );

	static LONGLONG_T CalcCenterOfHitOpt( ELEM_TYPE* p_data, 
	          		                LONG_T* cluster,LONG_T cluster_cnt,
	                               long xSize, long& x0, long& y0 );
	                      
	static int CalcCenterOfHitRealOpt( ELEM_TYPE* p_data, 
	          		                LONG_T* cluster,LONG_T cluster_cnt,
	                               long xSize, double& x0, double& y0 );
	                      
	static int CalcCenterOfHitRealOpt( BIG_ELEM_TYPE** p_data, 
	          		                LONG_T* cluster,LONG_T cluster_cnt,
	                               long xSize, double& x0, double& y0 );
	                      
	static int CalcCenterOfHitRealOpt( ELEM_TYPE** p_data, 
	          		                LONG_T* cluster,LONG_T cluster_cnt,
	                               long xSize, double& x0, double& y0 );
	                      

	LONG_T* GetAnalNeighbFromMap( long pos, long& ncnt );

	static long GetAnalNeighbWithSelf( long x,long y,long pos, long xSize,long ySize,
		                     		      LONG_T* neighb_list, LONG_T& ncnt );

	static long GetAllNeighbNoSelf( long x0,long y0,long pos, long xSize,long ySize,
                                   LONG_T* neighb_list, LONG_T& ncnt );

	static long GetAnalNeighbNoSelf( long x,long y,long pos, long xSize,long ySize,
		                     	      LONG_T* neighb_list, LONG_T& ncnt );

	static long GetAnalNeighbOpt( long x,long y,long xSize,long ySize,
                     		      LONG_T* neighb_list, LONG_T& ncnt,
	                                 LONG_T* out_list, LONG_T& ocnt );

	static long GetAnalNeighbOpt( long x,long y,long pos, long xSize,long ySize,
                     		      LONG_T* neighb_list, LONG_T& ncnt, BOOL_T bAddSelf=TRUE );

	static long GetAnalNeighb( long x, long y, long xSize,long ySize, 
	                           CLongList& Neighbours, CLongList& OuterNeighbours );

	static long GetNeighbours( long pos, long xSize,long ySize, 
	                           CLongList& Neighbours, BOOL_T bAddSelf=TRUE );

	static long GetNeighbours( long x,long y,long xSize,long ySize,
	                           CPointList& Neighbours );

	static long GetNeighbPositions( long x,long y,long xSize,long ySize,
	                                CLongList& Neighbours, BOOL_T bAddSelf=TRUE );

	// checking if satelite :
	static void VerifySatellitesOnSumFrame( CCDPipeline& pipeline1, CCDPipeline& pipeline2,
						CCDEventList& verifiedOnSum1, CCDEventList& verifiedOnSum2 );

	static BOOL_T CheckEventInSatDB( CCD_Analyser* pAnalyser1, CCD_Analyser* pAnalyser2,
					 CCDPipeline& pipeline1, CCDPipeline& pipeline2,
					 CccdReport* pEventOn1, CccdReport* pEventOn2,
					 BOOL_T& bIsSat );
					 
	static BOOL_T CheckEventInSatDB( CCD_Analyser* pAnalyser1, CCD_Analyser* pAnalyser2,
					 CCDPipeline& pipeline1, CCDPipeline& pipeline2,
					 CccdReport* pEventOn1, CccdReport* pEventOn2 );
						 

	static BOOL_T IsSatelite( CccdReport& evt1, CSatList* pList1,
			mystring& szSatName, double& min_dist1,
			mystring& szClosestSat1, double sat_coic_radius );


	static BOOL_T IsSatelite( CccdReport& evt1, CccdReport& evt2,
				  mystring& szSatName,
				  double& min_dist1, double& min_dist2,
                                  mystring& szClosestSat1, mystring& szClosestSat2,
                                  double sat_coic_radius );

	static BOOL_T IsSatelite( CccdReport& evt1, CccdReport& evt2, 
							 CSatList* pList1, CSatList* pList2, 
							 mystring& szSatName, 
							 double& min_dist1, double& min_dist2,
							 mystring& szClosestSat1, mystring& szClosestSat2,
							 double sat_coic_radius );
							 
	static BOOL_T IsStarTycho( CccdReport& evt1, CccdReport* evt2, mystring& szStarDesc );
	static BOOL_T IsStar( CccdReport& evt1, CccdReport* evt2, mystring& szStarDesc );
	static BOOL_T IsStarDefault( CccdReport& evt1, CccdReport* evt2, mystring& szStarDesc );


	// Function for finding clusters :

	BOOL_T CheckEventShape( const CPixelAnalyseIn& in, CPixelAnalyseOut& out );
	
	BOOL_T CheckIfHotPixel( CPixelAnalyseIn& in, CPointList& hotList );
	
	BOOL_T CheckBlackPixels( int* plusValues, int* minusValues,
	                         int plusValCount,int minusValCount,
	                         const CPixelAnalyseIn& in,
	                         CPixelAnalyseOut& out );
	static double GetBlackRatio( int* plusValues, int* minusValues,
	                      int plusValCount,int minusValCount );


	static double CalcSpherSigma( int x, int y, CCDMatrix& matrix, CCDPipeline* pPipeline, double sigmaN );
	
	static double CalcMaxClusterRadius( LONG_T* cluster, LONG_T cluster_cnt,
	                                    double x0, double y0, int xSize );

	static int FindCluster( int x0, int y0,
	                        ELEM_TYPE* data, ELEM_TYPE** data_fast,
	                        int xSize, int ySize,
	                        LONG_T* cluster, LONG_T& cluster_cnt,
	                        int tresh=150 );
	
	static int FindClusterAboveTresholdOpt2( const CPixelAnalyseIn& in,
	                                        double& x0,double& y0,
	                                        LONG_T* cluster,LONG_T& cluster_cnt,
	                                        double& maxNoiseLevel);

	static int CalcLapClusterOfMaxPoints( CPixelAnalyseOut& out,
													   int x, int y, int xSize, int ySize,
														ELEM_TYPE** p_data_fast, ELEM_TYPE* p_data,
														int n_points, 
														int n_break_size,
														int n_apert_size,
														float& sky, 
														float& plus_sum_out,
														double& x0, double& y0, BOOL_T bSubtractSky=TRUE );

	static int CalcLapClusterOfMaxPointsAnyShape( CPixelAnalyseOut& out,
													   int x, int y, int xSize, int ySize,
														ELEM_TYPE** p_data_fast, ELEM_TYPE* p_data,
														int n_points, 
														int n_break_size,
														int n_apert_size,
														float& sky, 
														float& plus_sum_out,
														double& x0, double& y0, BOOL_T bSubtractSky=TRUE );

	static int FindClusterOfMaxNPoints( int x0, int y0, int xSize, int ySize,
														 ELEM_TYPE** p_data_fast,
														 LONG_T* cluster,LONG_T& cluster_cnt,
													    int n_points );

	static int FindClusterOfMaxNPointsAnyShape( int x0, int y0, int xSize, int ySize,
														 ELEM_TYPE** p_data_fast,
														 LONG_T* cluster,LONG_T& cluster_cnt,
													    int n_points );


	BOOL_T CheckIfEdgeOfBigStar( CPixelAnalyseIn& in, CPixelAnalyseOut& out, 
										  BOOL_T bClusterFound, int nMaxLaplaceOnOther  );

	static int FindClusterAboveTresholdRawDataOpt3( CPixelAnalyseIn& in,
														  CPixelAnalyseOut& out,
 	 	                                      double& x0,double& y0,
 	 	                                      LONG_T* cluster,LONG_T& cluster_cnt,
														  double maxNoiseLevel,
														  BOOL_T bCalcCenter=TRUE,
														  BOOL_T bResetList=TRUE);

	static int FindClusterAboveTresholdRawDataOpt3Unique( CPixelAnalyseIn& in,
														  CPixelAnalyseOut& out,
 	 	                                      double& x0,double& y0,
 	 	                                      LONG_T* cluster,LONG_T& cluster_cnt,
														  double maxNoiseLevel,
														  CPixelList* usedList,
														  BOOL_T bCalcCenter=TRUE,
														  BOOL_T bResetList=TRUE);


	static  int GetLaplaceValue( int x, int y, int xSize, int ySize,
								ELEM_TYPE** p_data, BIG_ELEM_TYPE** p_laplace,
								BOOL_T bFromTab=TRUE );
/*	{
		if(x>=5 && y>=5 && x<(xSize-5) && y<(ySize-5)){
			return ( (bFromTab) ? p_laplace[y][x] : Table2D<ELEM_TYPE>::CalcLaplaceSum( x, y, xSize,p_data,gCCDParams.m_eLaplaceType )  );
		}
		return 0;					
	}								*/

	static int FindClusterAboveTresholdOpt3( int x, int y, int xSize, int ySize,
															ELEM_TYPE** p_data,
																CPixelAnalyseOut& out,
																LONG_T* cluster, LONG_T& cluster_cnt,
																double clusterTreshold=-1,
																BIG_ELEM_TYPE** p_laplace_data=NULL,
																CCDPipeline* pPipeline=NULL,
																CPixelList* pixel_list=NULL );

	static int FindClusterAboveTresholdOpt3( const CPixelAnalyseIn& in,
															 CPixelAnalyseOut& out,
															 LONG_T* cluster, LONG_T& cluster_cnt,
 	 	                                        double& x0,double& y0,
															 double& maxNoiseLevel,double clusterTreshold=-1,
															 BOOL_T bCalcCenter=TRUE, 
															 BOOL_T bFromTab=TRUE,
															 BOOL_T bResetList=TRUE );

	static int FindClusterAboveTresholdOpt3_LWP( const CPixelAnalyseIn& in,
															 CPixelAnalyseOut& out,
															 LONG_T* cluster, LONG_T& cluster_cnt,
 	 	                                        double& x0,double& y0,
															 double& maxNoiseLevel,double clusterTreshold=-1,
															 BOOL_T bCalcCenter=TRUE, 
															 BOOL_T bFromTab=TRUE,
															 BOOL_T bResetList=TRUE );

	static int FindClusterAboveTresholdOpt4( const CPixelAnalyseIn& in,
															 CPixelAnalyseOut& out,
															 LONG_T* cluster, LONG_T& cluster_cnt,
 	 	                                        double& x0,double& y0,
															 double& maxNoiseLevel,
															 double clusterTreshold,
															 BOOL_T bFromTab=TRUE );

	static inline int FindClusterAboveTresholdOpt3( const CPixelAnalyseIn& in,
														 CPixelAnalyseOut& out,
	                                        double& x0,double& y0,
	                                        double& maxNoiseLevel,double clusterTreshold=-1,
	                                        BOOL_T bCalcCenter=TRUE)
	{
		return FindClusterAboveTresholdOpt3( in, out, out.cluster, out.cluster_cnt,
		                                     x0,y0, maxNoiseLevel, clusterTreshold , bCalcCenter );	                                        	
	}

	static double FindClusterAboveTresholdOnPrevOpt3( CPixelAnalyseIn& in,
	                                        LONG_T& x0, LONG_T& y0,
	                                        LONG_T* cluster,LONG_T& cluster_cnt,
	                                        double& maxNoiseLevel,double clusterTreshold=-1,
	                                        BOOL_T bCalcCenter=TRUE);

	static BOOL_T FindClustersOnNewAndPrev( const CPixelAnalyseIn& in,
														 CPixelAnalyseOut& out,
														 CPixelAnalyseIn& inForPrev,
                                           CPixelAnalyseOut& outForPrev,
                                           double& max_on_new, double& max_on_prev  );

	static BOOL_T CheckForSN( const CPixelAnalyseIn& in,
                             CPixelAnalyseOut& out,
                             int TnForSN, int TvForSN, int MinPrevForSN );

	static BOOL_T CheckForSN_OnSum( const CPixelAnalyseIn& in,
                             CPixelAnalyseOut& out,
                             int TnForSN, int TvForSN, int MinPrevForSN );
                                                                                                             
                                           
	static BOOL_T FindClustersOnNewAndPrevNew( const CPixelAnalyseIn& in,
															  CPixelAnalyseOut& out,
															  CPixelAnalyseIn& inForPrev,
															  CPixelAnalyseOut& outForPrev,
															  double& max_on_new,double& max_on_prev );
                                           
                                           
	static BOOL_T FindAverageClusterOnPrevRawData( CPixelAnalyseIn& in, CPixelAnalyseOut& out,
														 LONG_T* cluster,LONG_T& cluster_cnt, int maxNoiseLevel,
														 int curr_x, int curr_y );   
														 
	static BOOL_T FindAverageClusterOnPrevLaplace( CPixelAnalyseIn& in, CPixelAnalyseOut& out,
																  LONG_T* cluster,LONG_T& cluster_cnt, int maxNoiseLevel,
																  int curr_x, int curr_y, int prevCount );
													

	static int FindMaxInCluster( ELEM_TYPE* p_data, int xSize,
										  LONG_T* cluster,LONG_T cluster_cnt,
										  int& max_pos );				

	static int CountClusterMAX( BIG_ELEM_TYPE** p_data, int xSize,
										 LONG_T* cluster,LONG_T cluster_cnt,
										 CPointList& max_list );

	static int CountClusterMAX( ELEM_TYPE** p_data, int xSize,
										 LONG_T* cluster,LONG_T cluster_cnt,
										 CPointList& max_list );

	static inline int FindMaxInClusterBig( BIG_ELEM_TYPE* p_data, int xSize,
										  LONG_T* cluster,LONG_T cluster_cnt,
										  int& max_pos )
	{
		register int max_val = -100000;
		max_pos = -1;
		for(register int i=0;i<cluster_cnt;i++){
			if(p_data[cluster[i]]>max_val){
				max_pos = cluster[i];
				max_val = p_data[max_pos];
			}
		}
		return max_val;
	}
										  
											  	
	static LONGLONG_T FindClusterAboveTresholdOpt(const CPixelAnalyseIn* pIn,ELEM_TYPE* p_data,
	          		              LONG_T x, LONG_T y,LONG_T CorePos,
	               		        LONG_T xSize, LONG_T ySize,
                         		  LONG_T& x0,LONG_T& y0,
                                CPixelList& in_clusters,                         		         
      		                    LONG_T* cluster,LONG_T& cluster_cnt);
	
	void FindClusterAboveTreshold( CCDMatrix& ccd_matrix,
				       LONG_T x, LONG_T y,LONG_T CorePos,
                                       LONG_T& x0,LONG_T& y0, double& r_max,
                                       CPixelList& in_clusters,CLongList& cluster );
	
	void FindMaxCluster(CCDMatrix& ccd_matrix,CCDMatrix& ccd_max,
	                    long x,long y,CPointList& cluster);

	void FindMaxClusterNew(CCDMatrix& ccd_matrix,CCDMatrix& ccd_max,
                          CLongList& startPoints,CLongList& cluster,
                          long CorePoint, CPixelList& in_clusters);


	void FindSumCluster(CCDMatrix& newMatrix,CCDPipeline& ccd_pipeline,
	                    long x,long y,CLongList& cluster,CPixelList& pixels);

	BOOL_T ConfirmEvent_MaxAlgorithm( CCDMatrix& ccd_matrix,CCDMatrix& ccd_max,
	                                  LONG_T CorePos, LONG_T xSize, LONG_T ySize );

	BOOL_T ConfirmEvent_MaxPrevSums( LONG_T MatrixIdx, CCDPipeline& ccd_pipeline,
	                                 CCDMatrix& ccd_matrix, LONG_T CorePos, 
	                                 LONG_T xSize, LONG_T ySize, double r0 );

	BOOL_T ConfirmEvent_InCluster( LONG_T MatrixIdx, CCDPipeline& ccd_pipeline,
	                               CCDMatrix& ccd_matrix, LONG_T CorePos, 
	                               LONG_T xSize, LONG_T ySize, 
	                               CLongList& cluster );

	// the 2 functions below MUST be unified into one !!!!
	BOOL_T ConfirmEvent_InClusterOpt_RotCorrected( 
							  const CPixelAnalyseIn& in,
	                    LONG_T* cluster, LONG_T cluster_cnt,
                       LONG_T& newClusterSum,LONG_T* prevClusterSum );

	BOOL_T ConfirmEvent_InClusterOpt( Table2D<ELEM_TYPE>** PrevMatrixPtr, LONG_T frames_back,
	 	                          ELEM_TYPE* p_data, LONG_T CorePos, 
                                          LONG_T xSize, LONG_T ySize, 
                                          LONG_T* cluster, LONG_T cluster_cnt );

	BOOL_T CheckBiggerArea(const CPixelAnalyseIn& in,
								  LONG_T* cluster, LONG_T& cluster_cnt,
                          LONG_T* ClusterWithMore, LONG_T& ClusterWithMoreCnt,
                          CPixelList& pixel_list,BOOL_T bTrace=TRUE);


	BOOL_T AnalyseSumOfNeighbours(CCDPipeline& ccd_pipeline);
	
	BOOL_T AnalysePrevMaxNew(CCDPipeline& ccd_pipeline);		
	
	BOOL_T AnalysePrevMaxInBuffer(CCDPipeline& ccd_pipeline);

	// full analysing function - to be used now 
	// it will use several methods on single point - depending on what 
	// confuiguration is
	BOOL_T FullNewFrameAnalyse(CCDPipeline& ccd_pipeline);

	BOOL_T ApplyCuts( int x, int y, int pos,
		               CPixelAnalyseIn& in, CPixelAnalyseOut& out,
	                  CFrameIdentStat& stat, CFrameIdentStat& sstat,
	                  BOOL_T bIsGenObject, BOOL_T& bDoBreak,
	                  CPixelList& allclusters  );

	static void FillVelocityHisto( double vx, double vy,
	                        double vx_new, double vy_new,
	                        int ccd_idx );
	                                                                            

	BOOL_T HistoVariables( int x, int y, int pos,
								  CPixelAnalyseIn& in, CPixelAnalyseOut& out,
								  CFrameIdentStat& stat, CFrameIdentStat& sstat,
								  CVariableInfo& eventVariables, double SigmaB );

	static void HistoDistance( CCDPipeline& ccd_pipeline1, CCDPipeline& ccd_pipeline2,
										CccdReport& evt1, CCDEventList& eventsOn2, 
										double x_prim, double y_prim, BOOL_T bRADEC );
	                  

	BOOL_T AnalyseCompareToOldFrame( CCDPipeline& ccd_pipeline );
	BOOL_T AnalyseSumOfPrevNFrames( CCDPipeline& ccd_pipeline );
	BOOL_T AnalyseSumOfPrevNFrames( CCDPipeline& ccd_pipeline, CCDEventList& event_list,
															 Table2D<BIG_ELEM_TYPE>* pCurrSumFrame,
															 Table2D<BIG_ELEM_TYPE>* pPrevSumFrame,
															 Table2D<BIG_ELEM_TYPE>* pCurrSumFrameLap,
															 Table2D<BIG_ELEM_TYPE>* pPrevSumFrameLap  );


	BOOL_T AnalyseNewFrameWithAverageOfPrevAlg_FromLogFiles( CCDPipeline& ccd_pipeline );
	BOOL_T AnalyseNewFrameWithAverageOfPrevAlg( CCDPipeline& ccd_pipeline );
	BOOL_T AnalyseNewFrameWithAverageOfPrevAlgMC( CCDPipeline& ccd_pipeline );
	
	// filling samples information :
	int FillSamplesInfo( CPixelAnalyseIn& in, int nNewLaplaceTreshold );
	
	BOOL_T CheckGenerated( CCDPipeline& ccd_pipeline );
	
	BOOL_T StartOptimizedAnalyse( CCDPipeline& ccd_pipeline );
	
	// only comparison with previous frames here :
	BOOL_T FullNewFrameAnalysePrev(CCDPipeline& ccd_pipeline);
	
	// verification of frame - by comparing with next frames :
	// returns number of verified 
	LONG_T VerifyEventsOnImage( CCDEventList& eventToVerify, CPixelAnalyseIn& in );

	// set next frame value :
	void SetNextFrameValues( CCDEventList& eventToVerify, CPixelAnalyseIn& in );

	// function adds additional information to each events from in.Matrix.FoundEvents list :
	static void AddAditionalEventInfo( CPixelAnalyseIn& in, CFrameIdentStat* stat=NULL );

	// analysing list of visible satelites :
	virtual int InitSatList( CCDPipeline& ccd_pipeline, time_t ut_time = 0 );

	// new version - trying gather all algorithms in single loop	
	virtual BOOL_T AnalyseNewFrameOpt(CCDPipeline& ccd_pipeline, BOOL_T bReport, LONG_T idx); 
	

	// obsolate will not be used now :
	virtual BOOL_T AnalyseNewFrame(CCDPipeline& ccd_pipeline,BOOL_T bReport=FALSE,LONG_T idx=0);


	// functions for calculating sums around given pixel, or sum of list of pixels	
	static LONG_T GetMaxAverageInVetoArea( CPixelAnalyseIn& in, BOOL_T bLaplace=FALSE );


	static inline void GetPrevPixelPosFromFormula( int& prev_x, int& prev_y, long x, long y,
	                                    int stepsBack, const CPixelAnalyseIn& in, int sec )	                                    
	{
		double prev_x_d,prev_y_d;
		(in.pCCDInfo)->xyAfterTimeNEW( (double)x, (double)y, sec, prev_x_d,prev_y_d );
		prev_x = my_round(prev_x_d);
		prev_y = my_round(prev_y_d);												
	}	                                    	
	                                    
	static void GetPrevPixelPos( int& prev_x, int& prev_y, int x, int y,
                 		           int stepsBack, const CPixelAnalyseIn& in );
                 		            
	// optimized - gets all positions in one call :
	static void GetPrevPixelPos( int x, int y, int* PrevFramesX, int* PrevFramesY,
	                             int backStart, int backTo, const CPixelAnalyseIn& in );
                                                                        
	static void CalcPrevPositionWithRot( long& prev_x, long& prev_y, long x0, long y0, 
													 int stepsBack, const CPixelAnalyseIn& in );

	static int CalcAverageOfPrevN( int x, int y, 
												  int xSize, int ySize,
												  int prevCount, 
												  BOOL_T bLaplace, CCDPipeline* pPipeline );

	static LONG_T CalcAverageOfPrevNInPixel( const long x, const long y,
	                        CPixelAnalyseIn& in, long prevCount,
	                        BOOL_T bLaplac=FALSE );
	                                                

	static LONG_T CalcAverageOfPrevN( const long x, const long y,
	  			          CPixelAnalyseIn& in, long prevCount,
	  			          BOOL_T bLaplace=FALSE );
	

	static LONG_T CalcAverageOfPrevNRecalc( const long x, const long y,
	  			          CPixelAnalyseIn& in, long prevCount,
	  			          BOOL_T bLaplace=FALSE,
	  			          CLongPoint* shapePoints=NULL, long pointCnt=0,
	  			          BOOL_T bPrintMax=FALSE );  	

	static LONG_T CalcAverageOfPrevN( CPixelAnalyseIn& in, long prevCount, 
						  			          BOOL_T bLaplace=FALSE,
	  			   					       CLongPoint* shapePoints=NULL, long pointCnt=0 );
	  			   					       
	static LONG_T CalcWeightedAverageOfPrevN( const CPixelAnalyseIn& in, long prevCount,
	                            InfoTable2D* pShiftInfo, BOOL_T bLaplace=FALSE );
	                            
	
	static LONGLONG_T CalcMaxSumOnPrev( LONG_T* neighb_list, LONG_T ncnt, 
						Table2D<ELEM_TYPE>** PrevMatrixPtr, LONG_T frames_back );

	static LONGLONG_T CalcMaxSumOnPrevRotCorrect( LONG_T* neighb_list, LONG_T ncnt, 
						const CPixelAnalyseIn& in );

	static LONGLONG_T CalcSum( ELEM_TYPE* p_data, long x, long y, long xSize, 
	                           long ySize );
	              
	static LONGLONG_T CalcSum( CLongList& List, const ELEM_TYPE* pData,  
	                           CPixelList* pUsedList=NULL );

	static LONGLONG_T CalcSumOpt( LONG_T* neighb_list, LONG_T ncnt, const ELEM_TYPE* pData );
	
														 	
	// laplacjan of Dorota and Bogumil (Gauss matrix) :
	static double CalcGaussMatrix( Table2D<double>& gaussMatrix, double fwhm, long radius, int mode );														 	

	// static LONGLONG_T CalcSumOptRotCorrected( LONG_T* neighb_list, LONG_T ncnt, const ELEM_TYPE* pData );
	//------------------------------------------------------------------------


	// function defining specific conditions and requirements :
	// static BOOL_T CheckLaplaceCondition( LONGLONG_T newLaplace, LONGLONG_T prevLaplace );
	static BOOL_T CheckLaplaceCondition( const CPixelAnalyseIn& in,
	                                     CPixelAnalyseOut& out );
	static BOOL_T CheckLaplaceConditionMedianMinus( const CPixelAnalyseIn& in,
					                                    CPixelAnalyseOut& out );
	
	static BOOL_T CheckBrighteningCond( LONGLONG_T newSum, LONGLONG_T otherSum );
	
	/*static BOOL_T CheckNextFrameCondition( LONGLONG_T SumToVerify,
	                                       LONGLONG_T NextSum,
	                                       LONG_T cluster_cnt,
                                          LONG_T nextFrameIndex );*/
                                          
	static BOOL_T CheckNextFrameCondition( CCDMatrix* pNextMatrix, int nextMatrixX, int nextMatrixY,
														CccdReport* pEvent, int nextFrameIndex, CPixelAnalyseIn& in );

	static int GetNextFrameValue( CCDMatrix* pNextMatrix, int nextMatrixX, int nextMatrixY,
											CccdReport* pEvent, int nextFrameIndex, CPixelAnalyseIn& in );

	static BOOL_T CheckLocalMaxCondition( const CPixelAnalyseIn& in, 
													  CPixelAnalyseOut& out,
	                                      LONG_T& max_pos );
	static BOOL_T IsLocalMax( long pos0, long xSize, ELEM_TYPE* p_data );
		                      
	static BOOL_T VerifyIfNotMoreThenNExceedsTreshold( const CPixelAnalyseIn& in,
																		CPixelAnalyseOut& out, 
																		LONG_T prevCount,
								                              LONG_T nAllowedToExceed,
			                                             BOOL_T bLaplace,
		                                                CLongPoint* shapePoints=NULL, long pointCnt=0);
		                                                

	static BOOL_T CheckCondition( LONGLONG_T sum, LONGLONG_T prev_sum_max );

	static BOOL_T CheckConfCondition( LONGLONG_T sum, LONGLONG_T prev_sum_max,
	       		          	          LONGLONG_T Tresh_signal, LONGLONG_T Tresh_noise );

	mystring GetClusterCheckDesc(LONGLONG_T newSum,LONGLONG_T MaxSum,
                  LONGLONG_T nConfTreshold,LONGLONG_T nMaxOnPrevAllowed);

	mystring GetCondDesc( LONGLONG_T newSum,LONGLONG_T maxSum,
	                      LONGLONG_T homeoSum );
	                      
	static void AutoCalculateIgnoreEdges( CCDPipeline* pPipeline );
	static void AutoCalculateIgnoreEdges();
	static double GetMaxNoiseLevel( const CPixelAnalyseIn* pIn );
	static double GetMaxNoiseLevelLaplace( const CPixelAnalyseIn* pIn, double nSigmaAbove=-1000.000 );
	static double GetMaxNoiseLevelLaplace( CCDPipeline* pPipeline, int ccd_index, double nSigmaAbove=-1000.000 );
	                                    	

	// tracks :
	void UpdateBestFit( const CPixelAnalyseIn& in, int startFrame, int EndFrame,
	                    CTrackDesc& old_track, CTrackDesc& new_track, const CccdReport& newEvent,
	                    CTrackList& tracks,
	                    BOOL_T bAddNewEventOnly=FALSE );		

	BOOL_T CheckSumEventsToNormalTracks( CCDPipeline* pPipeline );

	BOOL_T CheckNormalTracks( CCDPipeline* pPipeline, CCDEventList& events,
									  BOOL_T bDoCheckVelocity=TRUE );

	BOOL_T CheckIfEventBelongsToTrack( CccdReport& event,
												int confirmedFrameIndex,
												CCDPipeline* pPipeline,
									         BOOL_T bCheckVelocity, double fVelocityError,
									         CTrackDesc& track );

	BOOL_T CheckIfEventBelongsToTrack( CccdReport& event, CTrackDesc& track, BOOL_T& bReFitted,
												  double& new_a, double& new_b, 
												  BOOL_T bCheckVelocity, double velError,
												  BOOL_T bDoFit=TRUE );	
	BOOL_T VerifyIfEventsNotOnTrack( const CPixelAnalyseIn& in, CCDMatrix& Matrix );
	BOOL_T VerifyIfEventsNotOnTrackBase( const CPixelAnalyseIn& in, CCDMatrix& Matrix,
													 CTrackList& trackList, 
													 int MaxEventRate, int nNumBackFramesForTracks,
													 BOOL_T bCheckVelocity, double fVelocityError,
		                                  int nMaxEventsToFitTrack,
						                      int nNumBackFramesHasEventsForTracks,
						                      double MaxChi2InTrack, int nMinEventNoInTrack,
						                      double MaxChi2ForPointToMatchLine,
						                      int nCheckFramesIfNotOlderThenNFrames,
						                      eTrackCheckType_T track_type=eNormalTrack,
						                      BOOL_T bDoCheckForNewTracks=TRUE,
						                      int nMinFramesWithEvent=1 );
	static void AddNewNormalTracks( CCDPipeline* pPipeline );							                      
							                      
	BOOL_T VerifyPlaneTracks( const CPixelAnalyseIn& in, CCDMatrix& Matrix,
									  BOOL_T bDoCheckForNewTracks=FALSE,
									  BOOL_T bDoCheckVelocity=FALSE,
									  double fVelError=0.00 );

	BOOL_T VerifySingleCamTracks( const CPixelAnalyseIn& in, CCDMatrix& Matrix,
	                              CCDEventList& newEvents, CCDEventList& oldEvents );
	BOOL_T VerifySumFrameTracks( CCDPipeline* pPipeline, CCDMatrix& Matrix );
	                              
	                                                                                              
	BOOL_T VerifySingleCamTracks( const CPixelAnalyseIn& in, CCDMatrix& Matrix, 
																CTrackList& oldTracks,
															 CCDEventList& newEvents, CCDEventList& oldEvents,
																int MaxEventRate, int nNumBackFramesForTracks,
																BOOL_T bCheckVelocity, double fVelocityError, 
																int nMaxEventsToFitTrack, 
                                                int nNumBackFramesHasEventsForTracks,
																double MaxChi2InTrack, int nMinEventNoInTrack,
																double MaxChi2ForPointToMatchLine,
																eTrackCheckType_T track_type=eSingleCamTrack );
							                      
	
	
	static BOOL_T IsSingleFrameTrack( eTrackCheckType_T track_type );
	static BOOL_T DoAddFromSameFrame( eTrackCheckType_T track_type );
	
	BOOL_T CheckIfPlaneTrack(  vector<CccdReport*>& recentEvents,
	                           CTrackDesc& track,
	                           CCDPipeline* pPipeline,
	                           int minAddedForPlane );
	
	BOOL_T FlagEventsOnTrack( vector<CccdReport*>& recentEvents, 
									  CTrackDesc& track, CCDPipeline* pPipeline,
									  BOOL_T bCheckVelocity, double fVelocityError,
									  double MaxChi2InTrack, int nMinEventNoInTrack,
									  double MaxChi2ForPointToMatchLine,
									  eTrackCheckType_T track_type=eNormalTrack );

	int VerifyTracks( CCDPipeline& ccd_pipeline, int nConfirmOnNext=0, BOOL_T bForceCheck=FALSE );
	int VerifyPlaneTracks( CCDPipeline& ccd_pipeline, 
								  BOOL_T bDoCheckForNewTracks=FALSE,
								  BOOL_T bDoCheckVelocity=FALSE,
								  double fVelError=0.00  );
	
	
	BOOL_T FitLine( vector<CccdReport*>& rejectedList, CCDPipeline* pPipeline );
	BOOL_T CheckLineForEventsOnCurrFrame( CCDEventList& frameEvents );
	static BOOL_T FitLineToSingleFrameEvents( CCDEventList& recentEvents, int ccd_index,
															double fChi2ForLineOnSingleFrame,
															double& a_line, double& b_line,
															BOOL_T bCheckVelo=FALSE, double fVeloError=0.00 );
	BOOL_T FitLineToSingleFrameEvents( CCDEventList& recentEvents );

	// loging :
	void GetTrackLogName( eTrackCheckType_T track_type,
	                      mystring& szNewTrackLog, mystring& szTrackLog,
	                      mystring& szType, mystring& szSubDir );
	
	void LogNewTrack( CTrackDesc& track, BOOL_T bNew, eTrackCheckType_T track_type );
	void LogNewTrack_RADEC( CTrackDesc& track, BOOL_T bNew, 
									eTrackCheckType_T track_type, BOOL_T bDoConvert=TRUE );
	static int SaveTrack( CTrackDesc& track, int night,
	                BOOL_T bNew, eTrackCheckType_T track_type );
			
	
	static void LogNewTrack_RADEC( CTrackDesc& track, eTrackCheckType_T track_type,
					const char* szFileName, BOOL_T bNew=TRUE,
					CCDPipeline* pPipeline=NULL, BOOL_T bLocal=TRUE,
					CCD_Analyser* pAnal=NULL, BOOL_T bDoFit=FALSE,
					BOOL_T bDoConvert=TRUE );

	void LogEventRates( CCDPipeline* pPipeline );	
	void LogAllEventRates( CCDPipeline* pPipeline );
	void LogSumEventsRates( CCDPipeline* pPipeline );
	
	// checking event significance and if qualifies as internal trigger :
	static BOOL_T CheckIfInternalTrigger( CccdReport& newEvent );
	static BOOL_T CheckIfInternalTrigger( CccdReport& coicEvent1, CccdReport& coicEvent2 );
	
	// reporting :
	void PrepareEventReport( CCDPipeline& ccd_pipeline, LONG_T idx );		
	void UpdateEventReport( CCDPipeline& ccd_pipeline, LONG_T idx );
	CCDEventList& GetNewEvents(){ return m_AllEvents; }
	static CCDEventList& GetAllEvents(){ return m_AllRunEvents; }
	static void ClearRunEvents(){ m_AllRunEvents.clear(); }

	// repoting about frame analyse :
	void PrintFrameAnalyseInfo( CPixelAnalyseIn& in, int found_cnt );

	//       WARNING !!!
	// for fast calculations - statics can cause problems 
	// in multi-threaded applications :
};

#endif
