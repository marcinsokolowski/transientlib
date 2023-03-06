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
#ifndef _CCD_ROW_COL_LIST_H__
#define _CCD_ROW_COL_LIST_H__

#include "mytypes.h"

#include <vector>
using namespace std;

enum eCCDRowColType_T { eCCDUnknownType=0, eCCDRow, eCCDCol };

class CRowColDesc
{
public :
	CRowColDesc() : num(0),type(eCCDRow) {}
	CRowColDesc( int _num, eCCDRowColType_T _type ) :
		num(_num),type(_type) {}
	CRowColDesc( const CRowColDesc& right ){
		num = right.num;
		type = right.type;
	}		
	static eCCDRowColType_T GetType( const char* type_desc );
	
	int num;
	eCCDRowColType_T type;		
};

// list of row and columns for some specific action 
// can be read from file of syntax :
// COLUMN 1,23,54,900
// ROW    0,12,2047

#define ROW_INDICATOR "ROW"
#define COL_INDICATOR "COLUMN"

class CRowColList : public vector<CRowColDesc>
{
public:
	CRowColList(){};
	CRowColList( const char* filename );

	CRowColDesc* Find( eCCDRowColType_T type, int num );
	void AddUnique( CRowColDesc& elem );
	BOOL_T ReadFromFile( const char* filename );
};


class CCDWindow
{
public:
	int m_LowX;
	int m_LowY;
	int m_UpX;
	int m_UpY;

	CCDWindow() : m_LowX(0),m_LowY(0), m_UpX(0), m_UpY(0) {}
	CCDWindow( int low_x, int low_y, int up_x, int up_y )
		: m_LowX(low_x),m_LowY(low_y), m_UpX(up_x), m_UpY(up_y)		
	{}			
};


class CWindowList : public vector<CCDWindow>
{
public:
	CWindowList(){};
	CWindowList( const char* filename );

	BOOL_T ReadFromFile( const char* filename );	
};

#endif

