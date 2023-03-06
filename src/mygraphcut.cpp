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
#include "mygraphcut.h"
#include "calcrot.h"
#include <stdio.h>


CCutLineDef::CCutLineDef( double _a, double _b, double _c, eRelation _rel )
: a(_a), b(_b), c(_c), rel(_rel)
{

}

CCutLineDef::CCutLineDef( double _a, double _b, eRelation _rel )
: a(_a), b(-1.00), c(b), rel(_rel)
{

}


CGraphCut::CGraphCut()
{

}


void CGraphCut::AddLine( double a, double b, eRelation rel )
{
	CCutLineDef newline( a, b , rel );
	push_back( newline );
}

void CGraphCut::AddLine( double a, double b, double c, eRelation rel )
{
	CCutLineDef newline( a, b , c, rel );
	push_back( newline );
}

void CGraphCut::AddVerticalLine( double x_const, eRelation rel )
{
	AddLine( 1.00, 0.00, -x_const, rel );
}

void CGraphCut::AddHorizontalLine( double y_const, eRelation rel )
{
	AddLine( 0.00, 1.00, -y_const, rel );
}

void CGraphCut::AddLineByPoints( double x0, double y0, double x1, double y1,
		  								   eRelation rel )
{
	double a,b,c;
	CMyCalcRot::CalcLine( x0,y0,x1,y1,a,b,c  );
	
//	printf("LINE : %f,%f,%f\n",a,b,c);

	AddLine( a, b, c, rel );	
}


BOOL_T CGraphCut::CheckPoint( double x, double y )
{
	vector<CCutLineDef>::iterator i;

	for(i=begin();i!=end();i++){
		double a = i->a;
		double b = i->b;
		double c = i->c;

		if(b==0){
			// x = CONST
			// in this case smaller means - x<CONST , greater x>CONST
			double const_x = (-c/a);
			if( i->rel == eRelSmaller ){
				if( x>const_x )
					return FALSE;
			}
			if( i->rel == eRelGreater ){
				if( x<const_x )
					return FALSE;
			}
		}else{
			double y_on_line = (-a*x-c)/b;
			if( i->rel == eRelSmaller ){
				if(y_on_line<y)
					return FALSE;
			}
			if( i->rel == eRelGreater ){
				if(y<y_on_line)
					return FALSE;
			}
		}			
	}
	return TRUE;
}


const char* CGraphCut::GetRelDesc( eRelation rel )
{
	if(rel==eRelSmaller)
		return "SMALLER";
	if(rel==eRelGreater)
		return "GREATER";
	
	return "UNKNOWN";
}

void CGraphCut::Dump()
{
	vector<CCutLineDef>::iterator i;

	int j=0;
   for(i=begin();i!=end();i++){
		printf("Line %d : %.5f*x + %.5f*y + %.5f = 0 , relation=%s\n",j,i->a,i->b,i->c,GetRelDesc(i->rel));
		j++;
	}
}
