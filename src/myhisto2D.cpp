#include "myhisto2D.h"

CMyHisto2D::CMyHisto2D( const char* szName, double minX, double maxX, int bin_noX,
								double minY, double maxY, int bin_noY  )
: m_pCountTab(NULL)
{
	Init( szName, minX, maxX, bin_noX, minY, maxY, bin_noY );
}

CMyHisto2D::~CMyHisto2D()
{
	if( m_pCountTab ){
		delete [] m_pCountTab;
	}
}
	
void CMyHisto2D::Init( const char* szName, double minX, double maxX, int bin_noX,
				           double minY, double maxY, int bin_noY )
{
	m_MinValueX = minX;
   m_MaxValueX = maxX;
  	m_MinValueY = minY;
   m_MaxValueY = maxY;
  	m_BinNoX    = bin_noX;
   m_BinNoY    = bin_noY;
   m_BinWidthX = (maxX-minX)/bin_noX;
   m_BinWidthY = (maxY-minY)/bin_noY;
                                                                               
  	m_szHistoName = szName;                                                                               

	int count = bin_noX*bin_noY;
   m_pCountTab = new int[count];		
}

int CMyHisto2D::GetBinX( double x )
{
	double from_min = (x-m_MinValueX);
	int bin_noX = (int)(from_min/m_BinWidthX);
	return bin_noX;
}

int CMyHisto2D::GetBinY( double y )
{
	double from_min = (y-m_MinValueY);
	int bin_noY = (int)(from_min/m_BinWidthY);
	return bin_noY;
}

int& CMyHisto2D::GetCell( int binX, int binY )
{
	int pos = binY*m_BinNoY + binX;
	return m_pCountTab[pos];
}

void CMyHisto2D::FillCell( int binX, int binY )
{
	int& cell = GetCell( binX , binY );
	cell++; 
}
	
void CMyHisto2D::Fill( double x, double y )
{
	if(x>=m_MinValueX && x<m_MaxValueX && y>=m_MinValueY && y<m_MaxValueY){
		int bin_noX = GetBinX( x );
		int bin_noY = GetBinY( x );
		FillCell( bin_noX, bin_noY );
	}
}

