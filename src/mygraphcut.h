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
#ifndef _MY_GRAPH_CUT_H__
#define _MY_GRAPH_CUT_H__

#include <vector>
#include <math.h>
#include <mytypes.h>

using namespace std;

enum eRelation { eRelSmaller=0, eRelGreater };

class CCutLineDef
{
public :
	CCutLineDef( double _a, double _b, double _c, eRelation _rel );
	CCutLineDef( double _a, double _b, eRelation _rel );	
	
	double a;
	double b;
	double c;
	eRelation rel;
};

class CGraphCut : public vector<CCutLineDef>
{
public:
	CGraphCut();


	void AddLine( double a, double b, double c, eRelation rel );
	void AddLine( double a, double b, eRelation rel );
	void AddLineByPoints( double x0, double y0, double x1, double y1,
								 eRelation rel);
	void AddVerticalLine( double x_const, eRelation rel );
	void AddHorizontalLine( double y_const, eRelation rel );								 
								 
	BOOL_T CheckPoint( double x, double y );								 	
	
	static const char* GetRelDesc( eRelation rel );	
	void Dump();
};


#endif
