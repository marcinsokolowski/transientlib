#ifndef _TAB2D_DESC_H__
#define _TAB2D_DESC_H__

#include "mytypes.h"
#include "mypoints.h"
#include "myptrtable.h"

// neigbours of pixel will be kept in memory in continuos 
// memory area in order to get neighbours of pixel on position pos
// you have to call :
//
//        GetDesc(pos) 
//
// which is :
//
//       m_pDescTab[ pos*(m_Step+1) ]
//
// where m_Step is number of neighbours for example 5
// in case pixel has lees then m_Step neighbours ending is filled with 
// value NEIGHB_LIST_END (see below)
// first position stores number of neighbours !

#define NEIGHB_LIST_END -10000

// for quicker exectuion :
#define GET_DESC

class CMyHisto;

class CDescTab2D 
{
protected:
	long* m_pDescTab;
	LONG_T m_SizeX;
	LONG_T m_SizeY;		
	LONG_T m_Step;
	
	void Alloc();
public :
	void InitNeigb2D();
	CDescTab2D( LONG_T xSize, LONG_T ySize, LONG_T step );
	LONG_T* GetDesc( LONG_T pos, LONG_T& ncnt );
};

class  ShiftInfo
{
public:
	ShiftInfo():m_LocalDY(0),m_LocalDX(0),m_LocalS0(0),m_LocalS1(0),
				   m_LocalS2(0),m_LocalS3(0) {}
	
	double m_LocalDY;
	double m_LocalDX;
	double m_LocalS0;
	double m_LocalS1;
	double m_LocalS2;
	double m_LocalS3;
	CLongPoint overlapledPoints[4];
};

enum eFitValueType_T { eFitResUnknown=0, eFitResGauss=1, eFitResAVG=2 };

struct CDataInfo
{
	double m_Average;
	double m_Sigma;
	BOOL_T m_bFitOK;
	eFitValueType_T m_eFitResType;

	CDataInfo() : m_Average(0), m_Sigma(0), m_bFitOK(FALSE), m_eFitResType(eFitResUnknown) {};
};

class Area2DInfo
{
public:
	CMyPtrTable<CMyHisto> m_HistoList;

	Area2DInfo( LONG_T start_x, LONG_T start_y, LONG_T end_x, LONG_T end_y );	
	Area2DInfo();
	~Area2DInfo();
	void Init( LONG_T start_x, LONG_T start_y, LONG_T end_x, LONG_T end_y );		

	CPoint m_LowLeft;
	CPoint m_UpRight;
		
	// distributions :
	// LONG_T m_MostPopularValue;
	LONG_T m_MedianValue;

	CDataInfo m_DataInfo[MAX_LAPLACE_DEFINED];
	
	ShiftInfo m_LocalShiftInfo;
};

class InfoTable2D
{
public:
	LONG_T m_SizeX;
	LONG_T m_SizeY;
	LONG_T m_X_count;
	LONG_T m_Y_count;
	LONG_T m_dX;
	LONG_T m_dY;

	InfoTable2D();
	InfoTable2D( LONG_T SizeX, LONG_T SizeY, LONG_T dX, LONG_T dY, 
					 BOOL_T bInPixels=TRUE, int border=30 );	

	void InfoTable2D_InitConstructor( LONG_T SizeX, LONG_T SizeY, LONG_T dX, LONG_T dY,
				          BOOL_T bInPixels=TRUE, int border=30 );

	void Init( LONG_T SizeX, LONG_T SizeY, LONG_T dX, LONG_T dY, 
				  BOOL_T bInPixels=TRUE, int border=30 );	
	void Init2DMap( int border );
	~InfoTable2D();
	void Clear();
	
	inline Area2DInfo& GetAreaDesc( LONG_T x, LONG_T y )
	{
		register int x_elem = (x/m_dX);
		register int y_elem = (y/m_dY);		
		return m_pTable2DMap[y_elem][x_elem];
	}

	inline Area2DInfo& GetElem( LONG_T X_elem, LONG_T Y_elem )
	{
		return m_pTable2DMap[Y_elem][X_elem];
	}
	
	inline void GetAreaDescPos( LONG_T x, LONG_T y, int& x_elem, int& y_elem )
	{
		x_elem = (x/m_dX);
		y_elem = (y/m_dY);	
	}

	double GetSigmaBackground( long x, long y, eLaplaceType_T laptype );
	double GetSigmaBackground( long x, long y, eLaplaceType_T laptype, 
										BOOL_T& bFitOK );
	double GetMeanBackground( long x, long y, eLaplaceType_T laptype );
	BOOL_T GetMeanAndSigmaBackgr( long x, long y, double& mean, double& sigma, eLaplaceType_T laptype );

	void Dump();
	
	Area2DInfo** m_pTable2DMap;
};
#endif
