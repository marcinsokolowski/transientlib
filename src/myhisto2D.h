#ifndef _MY_HISTO_2D_H__
#define _MY_HISTO_2D_H__

#include <math.h>
#include <stdlib.h>
#include "mytypes.h"
#include "mystring.h"
#include "cexcp.h"

class CMyHisto2D
{
public:
	double m_MinValueX;
   double m_MaxValueX;
	double m_MinValueY;
	double m_MaxValueY;	   	
  	int m_BinNoX;
  	int m_BinNoY;
   double m_BinWidthX;
   double m_BinWidthY;

   mystring m_szHistoName;
   
   int* m_pCountTab;

// functions :
	CMyHisto2D( const char* szName, double minX, double maxX, int bin_noX,
					 double minY, double maxY, int bin_noY  );
	~CMyHisto2D();		
	
	void Init( const char* szName, double minX, double maxX, int bin_noX,
	           double minY, double maxY, int bin_noY );

	int GetBinY( double y );
	int GetBinX( double x );
	int& GetCell( int binX, int binY );		
	void FillCell( int binX, int binY );		
	
	void Fill( double x, double y );            	
};


#endif
