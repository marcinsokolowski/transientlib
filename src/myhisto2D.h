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
