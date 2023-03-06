#include "mymatrix.h"
#include "cexcp.h"
#include "myfile.h"
#include <stdlib.h>
#include "myparser.h"
#include "mystrtable.h"
#include <math.h>
#include "mymacros.h"

CMatrix::CMatrix()
: m_pElements(NULL), m_SizeN(0), m_SizeM(0)
{

}


CMatrix::CMatrix( int n, int m )
: m_pElements(NULL), m_SizeN(n), m_SizeM(m)
{
	Init( n, m );
}

CMatrix::CMatrix( int n )
: m_pElements(NULL), m_SizeN(n), m_SizeM(n)
{
	Init( n, n );
}

CMatrix::~CMatrix()
{
	if(m_pElements)
      delete [] m_pElements;
}


void CMatrix::Init( int n, int m )
{
	m_SizeN = n;
	m_SizeM = m;
	m_Size = (m_SizeN*m_SizeN);

	if(m_pElements)
		delete [] m_pElements;
	m_pElements = new double[m_Size];

	if( m_SizeN==m_SizeM )
		InitUnitMatrix();
	else
		InitZeroMatrix();
}

void CMatrix::InitZeroMatrix()
{
	for(register int i=0;i<m_Size;i++){
      m_pElements[i]=0;
   }
}


void CMatrix::InitUnitMatrix()
{
	for(register int i=0;i<m_Size;i++){
		m_pElements[i]=0;		
	}
	if(m_SizeN == m_SizeM){
		for(register int i=0;i<m_SizeN;i++){
			m_pElements[i+m_SizeN*i] = 1;
		}
	}
}

void CMatrix::SetElem( int row , int col, double val )
{
	Assert(row>=0 && row<m_SizeN && col>=0 && col<m_SizeM,"Cannot access %d,%d, matrix has sizes %d,%d",row,col,m_SizeM,m_SizeN);

	int pos = row*m_SizeN+col;	
	m_pElements[pos] = val;	
}

double& CMatrix::GetElem( int row , int col )
{
	Assert(row>=0 && row<m_SizeN && col>=0 && col<m_SizeM,"Cannot access %d,%d, matrix has sizes %d,%d",row,col,m_SizeM,m_SizeN);

	int pos = row*m_SizeN+col;	
	return m_pElements[pos];
}


void CMatrix::Dump()
{
	for(int row=0;row<m_SizeM;row++){
		for(int col=0;col<m_SizeN;col++){
			printf("%.2f ",GetElem( row, col ));
		}
		printf("\n");
	}
}

BOOL_T CMatrix::ReadFromFile( const char* fname )
{
	MyIFile in( fname, FALSE );

	if(!in.IsOpened())
		return FALSE;

	int row=0,col=0;
	const char* pLine;
	
	if(m_SizeN==0 || m_SizeM==0)
		return FALSE;	


	BOOL_T bRet = TRUE;
	while(pLine = in.GetLine()){
		if(pLine[0]=='\0' || pLine[0]=='\n')
			break;

		MyParser pars = pLine;
		CMyStrTable items;
		int newCol = pars.GetItems( items );

		if(row>=m_SizeM)
			break;

		for(int i=0;i<m_SizeN;i++){
			m_pElements[ row*m_SizeN + i  ] = atof( items[i].c_str() );
		}
		row++;
	}	

	return (row>0);
}


BOOL_T CMatrix::TimesVec( double* pVec, int nSize, double* outVec )
{
	if(nSize != m_SizeN)
		return FALSE;

	for(register int row=0;row<m_SizeM;row++){
		double sum=0;
		register int rowPos = row*m_SizeN;
		for(register int i=0;i<m_SizeN;i++){		
			sum += m_pElements[ rowPos + i ]*pVec[ i ];
		}
		outVec[row] = sum;
	}	

	return TRUE;
}


CMatrix3X3::CMatrix3X3()
: CMatrix( 3 )
{

}

void CMatrix3X3::TimesVec( double x1, double x2, double x3,
                  double& x1p, double& x2p, double& x3p )
{
	double vec[3],out[3];
	vec[0] = x1;
	vec[1] = x2;
	vec[2] = x3;
	CMatrix::TimesVec( vec, 3, out );
	x1p = out[0];
	x2p = out[1];
	x3p = out[2];
}

 void CMatrix3X3::TimesVec( double x, double y, double& xp, double& yp )
{
	double tmp;
	TimesVec( x, y, 1, xp, yp, tmp );
}


CTransformMatrixInTime::CTransformMatrixInTime( double alfa_per_sec, double dx_per_sec, double dy_per_sec,
																double center_x, double center_y )
: m_AlfaPerSec( alfa_per_sec ), m_DxPerSec( dx_per_sec ), m_DyPerSec( dy_per_sec ), m_CenterX(center_x), m_CenterY(center_y)
{

}

CTransformMatrixInTime::CTransformMatrixInTime()
: m_AlfaPerSec( 0 ), m_DxPerSec( 0 ), m_DyPerSec( 0 )
{

}

void CTransformMatrixInTime::Init( double alfa_per_sec, double dx_per_sec, double dy_per_sec, double center_x, double center_y )
{
	m_AlfaPerSec = alfa_per_sec;
	m_DxPerSec = dx_per_sec;
	m_DyPerSec = dy_per_sec;
	m_CenterX = center_x;
	m_CenterY = center_y;	
}

void CTransformMatrixInTime::PrintMatrix( const char* szFileName )
{
	mystring szDesc;
	char buff[128];
	sprintf(buff,"%.5f\t%.5f\t%.5f\n",cos(m_AlfaPerSec),-sin(m_AlfaPerSec),m_DxPerSec);
	szDesc << buff;
	sprintf(buff,"%.5f\t%.5f\t%.5f\n",sin(m_AlfaPerSec),cos(m_AlfaPerSec),m_DyPerSec);
	szDesc << buff;

	if(szFileName && szFileName[0]){
		MyOFile out( szFileName, "a+" );
   	out.Printf("%s",szDesc.c_str());
	}
	printf("%s\n",szDesc.c_str());
}

void CTransformMatrixInTime::Dump( const char* szFileName, int frame_index, const char* szAddInfo )
{
	mystring szDesc;
	char buff[128];

	szDesc << "Frame#" << frame_index << "\n";
	szDesc << "Transformation matrix with respect to point (" << (int)m_CenterX << "," << (int)m_CenterY << ")\n";
	sprintf(buff,"cos( %.7f *t )\t-sin( %.7f *t )\t %.5f*t\n",m_AlfaPerSec,m_AlfaPerSec,m_DxPerSec);
	szDesc << buff;		
	sprintf(buff,"sin( %.7f *t )\t+cos( %.7f *t )\t %.5f*t\n",m_AlfaPerSec,m_AlfaPerSec,m_DyPerSec);
	szDesc << buff << "\n";
	if(szAddInfo && szAddInfo[0])
		szDesc << "Additional info " << szAddInfo << "\n";
	
	_TRACE_PRINTF_2("%s",szDesc.c_str());
	if(szFileName && szFileName[0]){
		MyOFile out( szFileName, "a+" );
		out.Printf("%s",szDesc.c_str());
	}
}

void CTransformMatrixInTime::TimesVec( double dt_time, double x, double y, double& xp, double& yp )
{
	double dAlfa = dt_time*m_AlfaPerSec;
	double dx = dt_time*m_DxPerSec;
	double dy = dt_time*m_DyPerSec;

	double cos_alfa = cos(dAlfa);
	double sin_alfa = sin(dAlfa);
	
	xp = cos_alfa*x - sin_alfa*y + dx;
	yp = sin_alfa*x + cos_alfa*y + dy;
}

/*void CTransformMatrixInTime::TransformPoint( double dt_time, double x, double y, double& xp, double& yp )
{
	double x_prim = ( x-m_CenterX );
	double y_prim = ( y-m_CenterY );
	TimesVec( dt_time, x_prim, y_prim, xp, yp );
	xp += m_CenterX;
	yp += m_CenterY;
}*/

void CTransformMatrixInTime::TransformWithRespectTo( double dt_time, double x, double y, double rot_center_x, double rot_center_y,
                             							     double& xp, double& yp )
{
	double x_prim = ( x-rot_center_x );
	double y_prim = ( y-rot_center_y );
	TimesVec( dt_time, x_prim, y_prim, xp, yp );
	xp += rot_center_x;
	yp += rot_center_y;
}

