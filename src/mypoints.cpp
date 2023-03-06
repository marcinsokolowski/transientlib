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
#include "mypoints.h"
#include <algorithm>
#ifdef _DEBUG
#include <stdio.h>
#endif
#include "mathfunc.h"
#include "myfile.h"
#include "myparser.h"
#include "mystrtable.h"
#include "mymatrix.h"
#include "calcrot.h"

void CRectangle::Init(long _s_x,long _s_y,long _e_x, long _e_y, long col )
{
	start_x = _s_x;
	start_y = _s_y;
	end_x   = _e_x;
	end_y   = _e_y;
	color   = col;
}


CLongPoint::CLongPoint(LONG_T _x,LONG_T _y)
: x(_x),y(_y)
{

}

CLongPoint::CLongPoint(const CLongPoint& right)
{
	x = right.x;
	y = right.y;
}

CLongPoint::CLongPoint()
: x(0),y(0)
{
}


CPoint::CPoint(): x(0),y(0), m_FrameIndex(0), frame_time(0), m_Value(0), m_IdentType( eStarIdent_OK )
{
}

CPoint::CPoint(double _x,double _y,int frame_index,time_t _time, double val,eStarIdentType_T identType)
: x(_x), y(_y), m_FrameIndex(frame_index),frame_time(_time), m_Value(val), m_IdentType( identType )
{
	
}

CPoint::CPoint(const CPoint& right)
{
	x = right.x;
	y = right.y;
	m_FrameIndex = right.m_FrameIndex;
	frame_time = right.frame_time;
	m_Value = right.m_Value;
	m_IdentType = right.m_IdentType;
}


void CPoint::Set( LONG_T _x, LONG_T _y )
{
	x = (LONG_T)_x;
	y = (LONG_T)_y;
}

CPoint& CPoint::Subtract( CPoint& elem1, CPoint& elem2 )
{
	m_FrameIndex = elem1.m_FrameIndex;
	x = (elem1.x - elem2.x);
	y = (elem1.y - elem2.y);

	return (*this);
}

CPoint& CPoint::Add( CPoint& elem1, CPoint& elem2 )
{
	m_FrameIndex = elem1.m_FrameIndex;
	x = (elem1.x + elem2.x);
	y = (elem1.y + elem2.y);

	return (*this);
}

const char* CPoint::get_ident_desc( eStarIdentType_T identType, BOOL_T bNotEnd )
{
	if(!bNotEnd){
		return "AT END";
	}

	if(identType==eStarIdent_OK){
		return "OK";
	}
	if(identType==eStarIdent_ReUsedOld){
		return "PREDICTED";
	}

	return "OK";	
}

double CPoint::dist( double x1, double y1, double x2, double y2 )
{
	double ret = sqrt( CMyMathFunc::mysqr(x1-x2) + CMyMathFunc::mysqr(y1-y2) );
	return ret;
}

double CPoint::dist( CPoint& elem1, CPoint& elem2 )
{
	double ret = sqrt( CMyMathFunc::mysqr(elem1.x-elem2.x) + CMyMathFunc::mysqr(elem1.y-elem2.y) );
	return ret;
}

int CPoint::find_transform( CPoint& f1_s1, CPoint& f1_s2, CPoint& f2_s1, CPoint& f2_s2,
									 int size_x, int size_y,
									 CTransformMatrixInTime& matrix, int nLevel )
{
	CPoint shift_vec;
	shift_vec.Subtract( f1_s1, f2_s1 );
	
	CPoint f2_s2_shifted;
	f2_s2_shifted.Add( f2_s2, shift_vec );

	CPoint center = f1_s1;
	double alfa = CMyCalcRot::CalcRotAngle( f1_s2.x , f1_s2.y, 
														 f2_s2_shifted.x, f2_s2_shifted.y,
														 center.x, center.y );
	matrix.Init( alfa , -shift_vec.x, -shift_vec.y, center.x, center.y );		
	
	printf("transform center = (%.5f,%.5f)\n",center.x, center.y );
	int ret=1;
/*	if(nLevel==0){
		// CPoint frame_center( size_x/2 , size_y/2 );
		CPoint frame_center( 0, 0 );
		CPoint frame_center_prim;
		//matrix.TransformPoint( 1, frame_center.x, frame_center.y, 
		//							  frame_center_prim.x, frame_center_prim.y );
		CTransformMatrixInTime tmp;
		tmp.Init( -alfa, shift_vec.x, shift_vec.y, center.x, center.y );
		matrix.TransformPoint( 1, frame_center.x, frame_center.y,
                             frame_center_prim.x, frame_center_prim.y );
		//frame_center_prim.x = frame_center.x + shift_vec.x;
		//frame_center_prim.y = frame_center.y + shift_vec.y;
		ret = find_transform( frame_center, f1_s2, frame_center_prim, f2_s2, 
									 size_x, size_y, matrix, 1 );
	}*/

	return ret;
}

BOOL_T CPoint::IsBetween( double x0, double x1, double x )
{
	double min_x = MIN(x0,x1);
	double max_x = MAX(x0,x1);
	if( min_x<=x && x<=max_x )
		return TRUE;
	return FALSE;
}

BOOL_T CPoint::IsBetweenBig( double x0, double x1, double x, double min_dist )
{
	double min_x = MIN(x0,x1);
	double max_x = MAX(x0,x1);
	if( fabs( x-min_x )>min_dist || fabs( x-max_x )>min_dist || fabs( min_x-max_x )>min_dist ){
		if( x<min_x || max_x<x )
			return FALSE;
	}
	return TRUE;
}



CSector::CSector( double x1, double y1, double x2, double y2, int frame_index )
: m_Begin( x1, y1 ), m_End( x2, y2 ), m_FrameIndex(frame_index)
{

}

CSector::CSector( const CPoint& begin, const CPoint& end, int frame_index )
: m_Begin(begin), m_End(end), m_FrameIndex(frame_index)
{

}

void CSectorList::Add( double x1, double y1, double x2, double y2 )
{
	CSector sec( x1, y1, x2, y2 );
	push_back(sec);
}

void CSectorList::Add( const CPoint& begin, const CPoint& end )
{
	CSector sec( begin, end );
	push_back(sec);
}

void CSectorList::Add( const CSector& sec )
{
	push_back(sec);
}

void CSectorList::Clear()
{
	clear();
}

void CPointList::Add( int x, int y )
{
	CPoint tmp( x,y );
	push_back( tmp );
}

void CPointList::Add( const CPoint& new_elem )
{
	push_back( new_elem );
}

void CPointList::Clear()
{
	clear();
}


CPointList& CPointList::operator+=(CPointList& newElems)
{
	CPointList::iterator i;
	for(i=newElems.begin();i!=newElems.end();i++){
		Add( *i );
	}
	return (*this);
}


CPoint* CPointList::FindPoint( long x, long y )
{
	CPointList::iterator i;
	for(i=begin();i!=end();i++){
		if(i->x==x && i->y==y)
			return &(*i); 
	}
	return NULL;
}

CPoint* CPointList::FindByFrame( int frame_index )
{
	CPointList::iterator i;
   for(i=begin();i!=end();i++){
      if(i->m_FrameIndex == frame_index)
         return &(*i);
   }
	return NULL;
}

int CPointList::ReadFromFile( const char* filename )
{
	mystring szFile = filename;
	szFile.env2str();
	if(!MyFile::DoesFileExist( szFile.c_str() ))
		return 0;
	MyIFile in( szFile.c_str() );
	const char* pLine;
	while(pLine=in.GetLine()){
		if( mystring::get_first_non_white( pLine )!='#' ){
			MyParser pars=pLine;
			CMyStrTable items;
			pars.GetItems( items );
			if(items.size()>=2){
				CPoint tmp;
				tmp.x = atol( items[0].c_str() );
				tmp.y = atol( items[1].c_str() );

				if( items.size()>=3){
					tmp.m_FrameIndex = atol( items[2].c_str() );
				}
				push_back( tmp );
			}
		}
	}

	return size();
}


void CPointDescList::Add( const CPointDesc& new_elem )
{
	push_back( new_elem );
}


CPointDesc* CPointDescList::FindPoint( long x, long y )
{
	CPointDescList::iterator i;
   for(i=begin();i!=end();i++){
      if(i->x==x && i->y==y)
         return &(*i);
   }
   return NULL;
}

CLongList::CLongList()
{}

CLongList::CLongList(CLongList& newElems)
{
	(*this) = newElems;
}

void CLongList::Add(const LONG_T new_elem)
{
	push_back(new_elem);
}

void CLongList::AddUnique(const LONG_T new_elem)
{
	if(!FindPoint(new_elem))
		push_back(new_elem);
}

BOOL_T CLongList::FindPoint(LONG_T pos)
{
	CLongList::iterator i;
	for(i=begin();i!=end();i++){
		if((*i) == pos)
			return TRUE;
	}
	return FALSE;
}

void CLongList::Clear()
{
	clear();
}

CLongList& CLongList::operator+=(CLongList& newElems)
{
	CLongList::iterator i;
	for(i=newElems.begin();i!=newElems.end();i++){
		Add( *i );
	}
	return (*this);

}

CLongList& CLongList::operator=(CLongList& newElems)
{
	Clear();
	CLongList::iterator i;
	(*this) += newElems;
	return (*this);
}

CLongList& CLongList::operator=(const CLongList& newElems)
{
	Clear();
	for(int i=0;i<newElems.size();i++){
		push_back( newElems[i] );	
	}
	return (*this);
}

void CLongList::MergeUnique(CLongList& newElems)
{
	CLongList::iterator i;
	for(i=newElems.begin();i!=newElems.end();i++){
		if(!FindPoint(*i))
			Add( *i );
	}
}


void CLongList::Sort()
{
	sort(begin(),end());
}

void CLongList::RemoveLast()
{
	erase( end()-1 );
}

#ifdef _DEBUG
void CLongList::dump()
{
	CLongList::iterator i;
   for(i=begin();i!=end();i++){
		printf("%d\n",(*i));
	}
}
#endif
