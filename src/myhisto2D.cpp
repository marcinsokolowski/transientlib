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

