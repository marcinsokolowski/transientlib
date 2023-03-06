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

