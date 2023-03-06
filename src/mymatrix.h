#ifndef _MY_MATRIX_H__
#define _MY_MATRIX_H__

#include "mytypes.h"
#include <math.h>

// matrix N x M :

class CMatrix
{
protected:
	//  number of columns :
	int m_SizeN;
	
	// numer of rows :
	int m_SizeM;
	
	// numer of elements (m_SizeN*m_SizeM)
	int m_Size;
	
	double* m_pElements;		


public:
	CMatrix();
	CMatrix( int n, int m );
	CMatrix( int n );
	~CMatrix();
	
	double& GetElem( int row , int col );
	void SetElem( int row , int col, double val );

	void Init( int n, int m );
	void InitUnitMatrix();
	void InitZeroMatrix();
	
	// reading from file:
	BOOL_T ReadFromFile( const char* fname );

	// output :
	void Dump();


	// operations :
	BOOL_T TimesVec( double* pVec, int nSize, double* outVec );
};


class CMatrix3X3 : public CMatrix
{
public:
	CMatrix3X3();
	
	void TimesVec( double x1, double x2, double x3,
						double& x1p, double& x2p, double& x3p );
						
	void TimesVec( double x, double y, double& xp, double& yp );						
};

class CTransformMatrixInTime 
{
public:
	CTransformMatrixInTime();
	CTransformMatrixInTime( double alfa_per_sec, double dx_per_sec, double dy_per_sec, double center_x, double center_y );

	void Init( double alfa_per_sec, double dx_per_sec, double dy_per_sec, double center_x, double center_y );
	
	double m_AlfaPerSec;
	double m_DxPerSec;
	double m_DyPerSec;
	double m_CenterX;
	double m_CenterY;

	void Dump( const char* szFileName, int frame_index, const char* szAddInfo="" );
	void PrintMatrix( const char* szFileName );		

	inline void TransformPoint( double dt_time, double x, double y, double& xp, double& yp )
	{
		double x_prim = (x-m_CenterX);	
		double y_prim = (y-m_CenterY);
		// TimesVec( dt_time, , ( y-m_CenterY ), xp, yp );

		double dAlfa = dt_time*m_AlfaPerSec;
		// double dx = dt_time*m_DxPerSec;
		// double dy = dt_time*m_DyPerSec;

		double cos_alfa = cos(dAlfa);
		double sin_alfa = sin(dAlfa);
	
		xp = cos_alfa*x_prim - sin_alfa*y_prim + dt_time*m_DxPerSec + m_CenterX;
		yp = sin_alfa*x_prim + cos_alfa*y_prim + dt_time*m_DyPerSec + m_CenterY;
	}
	
	void TimesVec( double dt_time, double x, double y, double& xp, double& yp );
	void TransformWithRespectTo( double dt_time, double x, double y, double rot_center_x, double rot_center_y,	
										  double& xp, double& yp );
};

#endif

