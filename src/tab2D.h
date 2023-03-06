#ifndef _TAB2D_H__
#define _TAB2D_H__

#include <stdlib.h>
#include "mytypes.h"
#include "mysafekeytab.h"
#include "basedefines.h"
#include <limits.h>

#define DEF_ELEM 0
#define MAX_VALUE 32000
#define MAX_VALUE_ALLOWED 50000

class CMyHisto;
class ImageStat;
class InfoTable2D;
class CLongPoint;
class ShiftInfo;
class Area2DInfo;
class mystring;
class CRowColList;
class CPointList;
class CWindowList;

enum eTableType  { eTable2D=0, eCCDMatrix };


template<class ARG_TYPE>
class Table2D {
public:
	virtual eTableType GetType(){ return eTable2D; }

	BOOL_T m_bAllocHere;
	ARG_TYPE* m_pData;
	
	long m_SizeX;
	long m_SizeY;
	long m_Size;
	LONG_T m_Index;		
	ARG_TYPE** m_pFastData;
	
	
	// for storing FITS header data :
	CSafeKeyTab m_KeyTab;

	// for fast statistics :
	LONG_T  m_MaxValue;
	LONG_T* m_pDistribTab;

	// statistics :
	LONG_T* m_pValuesTable;
	
	LONG_T* m_pWrkTable;

	void InitDistribTab( LONG_T MaxValue );
	void GetDistrib( LONG_T* pDistribTab, LONG_T maxValue );
	void ClearDistribTab();
	
	// rotation calculations :
	static void CalcOverlapParts( double dx, double dy, double& s1, 
								  double& s2, double& s3, double& s4,
								  CLongPoint* overlapedPixels );
	static void CalcOverlapParts2( double dx, double dy, ShiftInfo& pShiftInfo );
																		  
									  	
public:
	// stat :
	// stat :
	LONG_T m_MedianValue;	
   LONG_T m_MaxVal;
   LONG_T m_MaxX;
   LONG_T m_MaxY;
   double m_MeanValue;

	// for geting image from bigger one :
	long m_X_On_Big; // it is used for astrometric calculations, it may be disabled in case image
	long m_Y_On_Big; // has its own astrometry ( aver20 ) pipeline ( ASTTHIS=1 ) , in such a case it is ignored
	int  m_StartOnOrigX; // filled always to know where image started on original image 
	int  m_StartOnOrigY; // it maybe used for bad rows/cols 

	// transalation/rotation corrections :
	LONG_T m_PrevShiftX;
	LONG_T m_PrevShiftY;
	
	double m_totalShiftX;
	double m_totalShiftY;
	      

	Table2D();
	Table2D(long x_size,long y_size,long idx=0,BOOL_T bAllocHere=TRUE);
	Table2D(const Table2D& right);
	
	// must do the same as contructor Table2D(long x_size,long y_size,long idx=0,BOOL_T bAllocHere=TRUE)
	// added due to gcc4.0 requirements :
	BOOL_T InitConstructor(long x_size,long y_size,long idx=0,BOOL_T bAllocHere=TRUE);		

	void SetDataPtr( ARG_TYPE* pData );
	
	virtual ~Table2D();	
	inline long GetXSize() const { return m_SizeX; }
	inline long GetYSize() const { return m_SizeY; }
	inline int GetSize() const { return m_Size; }
	ARG_TYPE* GetPos(long x,long y);
	ARG_TYPE& GetPixel(long x,long y);
	inline ARG_TYPE* get_data_buffer(){ return m_pData; }
	inline ARG_TYPE* get_data_buffer() const{ return m_pData; }
	inline ARG_TYPE** get_data_buffer_fast(){ return m_pFastData; }		
	inline ARG_TYPE** get_data_buffer_fast() const{ return m_pFastData; }		
	inline void SetIndex( LONG_T idx ){ m_Index = idx; }
   inline LONG_T GetIndex() { return m_Index; }
   inline CSafeKeyTab& GetKeyTab(){ return m_KeyTab; }

   BOOL_T Alloc(long x_size,long y_size);

   // making smaller image from larger :
   void CutBorderBase( int left, int up, int right, int bottom, Table2D<ARG_TYPE>& out_image );
   
   // do not change buffer ( it is still large ), only change sizes :
   // WARNING : enlarging will not work - it must be implemented !!!
   BOOL_T ReSize(long x_size,long y_size);

   void ReCalcIndex();
   void ReCreateIndex();
   
   ARG_TYPE getval(long x,long y);
	void Init(ARG_TYPE val=DEF_ELEM);
	void SetData(ARG_TYPE val);
	ARG_TYPE& get(long x,long y);
	void InitRandom(long min_value=0,long max_value=255);	

	void copy( ARG_TYPE* src );
	Table2D<ARG_TYPE>& Assign(const Table2D<ARG_TYPE>& right);	

	Table2D<ARG_TYPE>& OperatorEq(const Table2D<ARG_TYPE>& right);
	Table2D<ARG_TYPE>& operator=(const Table2D<ARG_TYPE>& right);
	Table2D<ARG_TYPE>& operator-=(const Table2D<ARG_TYPE>& right);		
	
	void ClearState();
	
	
	void GetXY_FromPos( LONG_T pos, LONG_T& x ,LONG_T& y );
	LONG_T GetPos_FromXY( LONG_T& x ,LONG_T& y );
	
	const char* getKey(const char* keyname);
	
	virtual BIG_ELEM_TYPE* get_frame_laplace(){ return NULL; }
	virtual BIG_ELEM_TYPE** get_frame_laplace_fast(){ return NULL; }

	int GetPixelsAbove( double treshold, CPointList& pointList, int start_x=0, int end_x=0, int start_y=100000, int end_y=100000 );

	// simple statistics functions :
	void GetSafeStatistics( LONG_T start_x, LONG_T start_y,
	                   		LONG_T end_x, LONG_T end_y,
	                    		double& mean, double& RMS );


	int GetBackgroundValue( LONG_T start_x, LONG_T start_y, 
									LONG_T end_x, LONG_T end_y,
									LONG_T* buffer, int nSize );

	void GetStatistics( LONG_T start_x, LONG_T start_y,
	                    LONG_T end_x, LONG_T end_y,
	                    double& mean, double& RMS );

   void GetSafeStatistics( ImageStat& statInfo, 
   							  BOOL_T bDistrib=FALSE,
                  	     LONG_T start_x=0, LONG_T start_y=0,
                     	  LONG_T end_x=-1, LONG_T end_y=-1,
	                       BOOL_T bDoFit=FALSE );

	
   void GetStatistics( ImageStat& statInfo, 
   						  BOOL_T bDistrib=FALSE,
                       LONG_T start_x=0, LONG_T start_y=0,
                       LONG_T end_x=-1, LONG_T end_y=-1,
                       BOOL_T bDoFit=FALSE,
                       double min_value=-1000000, double max_value=1000000 );

   void GetPartStat( ImageStat& statInfo, 
   					   BOOL_T bDistrib=FALSE,
                     LONG_T start_x=0, LONG_T start_y=0,
                     LONG_T end_x=-1, LONG_T end_y=-1 );

   void GetPartDistrib(  LONG_T* pDistribTab, LONG_T maxValue=MAX_VALUE,
                         LONG_T start_x=0, LONG_T start_y=0,
                         LONG_T end_x=-1, LONG_T end_y=-1 );

	void FillHisto( eLaplaceType_T laplaceType, int low_x, int up_x,
                   int low_y, int up_y,
						 double& mean, double& rms, double& max_val,
						 CMyHisto* pHisto );

	void FillHistoFromData( eLaplaceType_T laplaceType, int low_x, int up_x,
                   int low_y, int up_y,
						 double& mean, double& rms, double& max_val,
						 CMyHisto* pHisto, BIG_ELEM_TYPE** p_data );



	void GetVariableSigma( eLaplaceType_T laplaceType,double* typicalSigmaTable,
				double mean, double& sigma, Area2DInfo& info,
				LONG_T start_x=0, LONG_T start_y=0,
                                LONG_T end_x=-1, LONG_T end_y=-1 );
                   
	void GetVariableMeanAndSigma( eLaplaceType_T laplaceType,
	                              double& mean, double& sigma,
	                              LONG_T start_x=0, LONG_T start_y=0,
	                              LONG_T end_x=-1, LONG_T end_y=-1 );
                         
	BOOL_T GetVariableMeanAndSigma( eLaplaceType_T laplaceType,
											double& mean, double& sigma, Area2DInfo& info,
											LONG_T start_x=0, LONG_T start_y=0,
	  	                           LONG_T end_x=-1, LONG_T end_y=-1, 
	  	                           int frame_index=0,
	  	                           BIG_ELEM_TYPE** p_lap_data=NULL,
	  	                           BOOL_T bUseMeanOnFitFailed=FALSE );
	  	                           
	void GetMeanAndRMS( double& mean, double& rms );   
	void GetMean( double& mean, int border=30 );
	void GetMean( int border=30 );
	
	                     
	LONG_T GetMostPopular( LONG_T MaxValue=MAX_VALUE,
	                       LONG_T start_x=0, LONG_T start_y=0,
                          LONG_T end_x=-1, LONG_T end_y=-1 );

	LONG_T GetMedianValue(){ return m_MedianValue; }

	LONG_T GetMedian( LONG_T MaxValue=MAX_VALUE,
	                  LONG_T start_x=0, LONG_T start_y=0,
                     LONG_T end_x=-1, LONG_T end_y=-1 );

	LONG_T GetFastMedian( LONG_T MaxValue=MAX_VALUE,
	                  LONG_T start_x=0, LONG_T start_y=0,
                     LONG_T end_x=-1, LONG_T end_y=-1 );

	LONG_T GetMaxValueAndPos( long start_x, long start_y, 
									  long end_x, long end_y,
									  long& max_x, long& max_y );
									  
	int GetMaxValue( int& max_x, int& max_y );  
									  
	void CalcMaxAndPos();									  	

	void SubtractBacgroundInFrame( long nFrameSize );
	
	void SubtractConst( long value, BOOL_T bZero=FALSE );
	
	void Subtract( Table2D<ARG_TYPE>& right,Table2D<ARG_TYPE>& result, BOOL_T bZero);
	
	void MultiplyByConst( double mult );

	void Multiply(Table2D<ARG_TYPE>& right,Table2D<ARG_TYPE>& result,
			int nNorm, int max_val=USHRT_MAX);

	void FlipImage( eDriverReverseImage_T flipType );			

	static void FlipImage( ARG_TYPE* data, int sizeX, int sizeY,
          			        eDriverReverseImage_T flipType );
          			        
	void FH();
	
	void FV();          			        
          			        
	void ClearBorder( int border );          			        
		
	void GetDistrib( ImageStat& statInfo, LONG_T MaxValue=MAX_VALUE );

	void GetBackgrDescStr( InfoTable2D& info, mystring& szDesc );

	void CalcTableMapInfo( InfoTable2D& info,
								  LONG_T maxValue, 
								  BOOL_T* pBackgrFlagTable,
								  eLaplaceType_T currLaplace,
								  BIG_ELEM_TYPE** p_lap_data,
								  BOOL_T bCalcMostPopular=FALSE, 
	                       BOOL_T bCalcMedian=FALSE,
	                       int frame_index=0,
	                       int border=0 );
	                       
	BOOL_T CalcStat( CRowColList& rowcollist, double& rms, 
						  double& mean, BOOL_T bFitGauss=FALSE, 
						  int frame_index=0 );

	BOOL_T CalcStat( CWindowList& rowcollist, double& rms, 
						  double& mean, BOOL_T bFitGauss=FALSE,
						  int frame_index=0  );
						  
	static BOOL_T CalcStat( ARG_TYPE** data, int sizeX, int sizeY,
									CWindowList& rowcollist, double& rms,
									double& mean, BOOL_T bFitGauss=FALSE,
									int frame_index=0  );						  

	static BOOL_T CalcStat( ARG_TYPE* data, int sizeX, int sizeY,
									CWindowList& rowcollist, double& rms,
									double& mean, BOOL_T bFitGauss=FALSE,
									int frame_index=0  );						  

	LONGLONG_T CalcSum( LONG_T* pixels, LONG_T cnt );	   		

	BOOL_T MarkStarWithCircle( int x0, int y0, int radius=5, int width=2,
										int value=50000 );

	BOOL_T GetImage( long x0, long y0,long len_x, long len_y,
	                 Table2D<ARG_TYPE>& Image );	

	BOOL_T GetImageELEMTYPE( long x0, long y0,long len_x, long len_y,
	                 Table2D<ELEM_TYPE>& Image );	

	BOOL_T GetImageBIGELEMTYPE( long x0, long y0,long len_x, long len_y,
	                 Table2D<BIG_ELEM_TYPE>& Image );	

	BOOL_T GetImageWithZeros( long x0, long y0,long len_x, long len_y,
	                 Table2D<ARG_TYPE>& Image );	

	// rotation/translation :
	void RotateFrame( LONG_T FrameCounter, double dAlpha, double x_center, double y_center );
	
	void ShiftFrame( LONG_T FrameCounter, double frame_dx, double frame_dy, BOOL_T bUseTotal=TRUE );

	void ShiftFrameWeighted( LONG_T FrameCounter, double frame_dx, double frame_dy );
	
	
	void CalcTotalShift( LONG_T FrameCounter, double frame_dx, double frame_dy );		
	void GetTotalShift( LONG_T FrameCounter, double frame_dx, double frame_dy,
	                    double& totalShiftX, double& totalShiftY, LONG_T& nSteps );
	// 
	
	
	
   LONGLONG_T CalcSumRotCorrected( LONG_T* pixels, LONG_T cnt,
                                   double dx, double dy );

   void GetXYRotCorrected( LONG_T x, LONG_T y, LONG_T dx, LONG_T dy, 
                           LONG_T& x_corr, LONG_T& y_corr );

	// write FITS :
	virtual void WriteToFITSFile( const char* fname="$(DATADIR)/image.fit", BOOL_T bDumpEvents=FALSE,
											CKeyTab* pHduTab=NULL ){}										
                                   

	// arithmetic operations :
	BOOL_T Compare( const Table2D<ARG_TYPE>& right );
		
	void Dump();
	
	void Dump( long _x0, long _y0, long _x1, long _y1 );

// laplace:
	static inline int CalcG54( LONG_T x, LONG_T y, LONG_T xSize, ARG_TYPE** p_data ){
		return (int)( ( p_data[y][x]+p_data[y-1][x]+p_data[y+1][x]+p_data[y][x-1]+p_data[y][x+1])
					-1.25*( p_data[y-2][x-2]+p_data[y-2][x+2]+p_data[y+2][x-2]+p_data[y+2][x+2] ) );
	}

	static inline int CalcG5( int x, int y, int xSize, ARG_TYPE** p_data ){
		return (int)( ( p_data[y][x]+p_data[y-1][x]+p_data[y+1][x]+p_data[y][x-1]+p_data[y][x+1])
					-1.25*( p_data[y-1][x-1]+p_data[y-1][x+1]+p_data[y+1][x-1]+p_data[y+1][x+1] ) );
	}


	static inline LONG_T CalcLaplaceSum( LONG_T x, LONG_T y, LONG_T xSize, ARG_TYPE** p_data,
										eLaplaceType_T laplaceType, BOOL_T bMedianBkg=FALSE );

	static inline LONG_T CalcLaplaceSumBig( LONG_T x, LONG_T y, LONG_T xSize, BIG_ELEM_TYPE** p_data,
										eLaplaceType_T laplaceType );

	static inline LONG_T CalcLaplaceSum( LONG_T x, LONG_T y, LONG_T xSize, ARG_TYPE** p_data,
										eLaplaceType_T laplaceType, LONG_T& plus_sum, LONG_T& minus_sum,
										BOOL_T bMedianBkg=FALSE );

	static inline LONG_T CalcLaplaceSumBig( LONG_T x, LONG_T y, LONG_T xSize, BIG_ELEM_TYPE** p_data,
										eLaplaceType_T laplaceType, LONG_T& plus_sum, LONG_T& minus_sum );
										

	// OBSOLATE - do not USE !
	//static LONG_T CalcLaplaceSumOLD( LONG_T x, LONG_T y, LONG_T xSize, ARG_TYPE** p_data,
	//									eLaplaceType_T laplaceType, LONG_T& plus_sum, LONG_T& minus_sum );
	//
	//static LONG_T CalcLaplaceSumBigOLD( LONG_T x, LONG_T y, LONG_T xSize, BIG_ELEM_TYPE** p_data,
	//									eLaplaceType_T laplaceType, LONG_T& plus_sum, LONG_T& minus_sum );

	static LONG_T CalcLaplaceSumEstimateOnEdge( LONG_T x, LONG_T y, 
															  LONG_T xSize, LONG_T ySize,
   								               		ARG_TYPE** p_data,
   								               		CLongPoint* plus_list, LONG_T plus_cnt,
   								               		CLongPoint* minus_list, LONG_T minus_cnt );

	static void GetPlusMinusList( eLaplaceType_T laplaceType,
											CLongPoint* plus_list, LONG_T& plus_cnt,
					  					   CLongPoint* minus_list, LONG_T& minus_cnt );

	static void GetLaplacePlusMinusValues( ARG_TYPE** p_data,
														int* plusValues, int* minusValues,
														int& plusValCount,int& minusValCount,
														LONG_T x, LONG_T y,
														LONG_T xSize, LONG_T ySize,
														eLaplaceType_T laplaceType );

	static LONG_T GetLaplacePlusMinusList( LONG_T* plusTab, LONG_T& plusCnt,
	                                       LONG_T* minusTab, LONG_T& minusCnt,
	                                       LONG_T x, LONG_T y, LONG_T pos,
                                          LONG_T xSize, eLaplaceType_T laplaceType );
					  					   
					  					   
	static const char* GetLaplaceName( eLaplaceType_T laplaceType );					  					   
	

	// calculating laplace of frame :	
	void Laplace( Table2D<BIG_ELEM_TYPE>* pFrameLaplace, 
					  eLaplaceType_T laplaceType,
					  BOOL_T bEstimateOnEdge=TRUE );
	
	void Laplace( BIG_ELEM_TYPE** p_laplace_fast, eLaplaceType_T laplaceType,
					  long x0, long y0, long x1, long y1 );

	void LaplaceSafe( BIG_ELEM_TYPE** p_laplace_fast, 
							CLongPoint* LaplacePlusList, int LaplacePlusCount,
							CLongPoint* LaplaceMinusList, int LaplaceMinusCount,							
							long x0, long y0, long x1, long y1 );
	
	void Laplace( Table2D<ARG_TYPE>& out, eLaplaceType_T laplaceType,
					  int ignore_edge,
					  long x0=-1, long y0=-1, long x1=-1, long y1=-1 );



	// generator :
	void GenEmptyImage( double sigma, double mean );
};

#endif
