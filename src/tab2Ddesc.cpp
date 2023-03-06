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
#include "tab2Ddesc.h"
#include <stdio.h>
#include <math.h>
#include "myhisto.h"


void CDescTab2D::Alloc()
{
	LONGLONG_T total_size = ( m_SizeX*m_SizeY*m_Step );
	m_pDescTab = new LONG_T[ total_size ];
}


// m_Step = (step+1) - to keep number of neighbours on first position !
CDescTab2D::CDescTab2D( LONG_T xSize, LONG_T ySize, LONG_T step )
: m_SizeX(xSize), m_SizeY(ySize), m_Step(step+1)
{
	Alloc();
}

void CDescTab2D::InitNeigb2D()
{
	LONG_T size = ( m_SizeX*m_SizeY );
	LONG_T pos=0;
	for(pos=0;pos<size;pos++){
		LONG_T start_pos = pos*m_Step;
		LONG_T x = ( pos % m_SizeX );
		LONG_T y = ( pos / m_SizeX );
		LONGLONG_T new_pos = start_pos+1;

      // add neighbours sorted !
      if(y>0){
	      m_pDescTab[ new_pos ] = pos-m_SizeX;
         new_pos++;
      }
      if(x>0){
	      m_pDescTab[ new_pos ] = pos-1;
         new_pos++;
      }

      // curent point - only to neighb - not outer
      m_pDescTab[ new_pos ] = pos;
      new_pos++;

      if(x<(m_SizeX-1)){
 	       m_pDescTab[ new_pos ] = pos+1;
          new_pos++;
      }
      if(y<(m_SizeY-1)){
     		 m_pDescTab[ new_pos ] = pos+m_SizeX;
          new_pos++;
      }

		if(new_pos != ((pos+1)*m_Step))
			m_pDescTab[ new_pos ] = NEIGHB_LIST_END;

		LONG_T cnt = (new_pos - (start_pos+1));
		m_pDescTab[ start_pos ] = cnt;
	}
}

LONG_T* CDescTab2D::GetDesc( LONG_T pos, LONG_T& ncnt )
{
	// skip first possition - where numer of neighbours is stroed 
	LONG_T* start_pos = m_pDescTab + (pos*m_Step);
	ncnt = start_pos[0];		// number of neighbours is stored in first position
	return ( start_pos + 1 );	
}



Area2DInfo::Area2DInfo( LONG_T start_x, LONG_T start_y, LONG_T end_x, LONG_T end_y )
: m_MedianValue(0)
{
	Init( start_x,start_y,end_x,end_y);
}

Area2DInfo::Area2DInfo()
: m_MedianValue(0)
{
	Init( 0, 0, 0, 0);
}

void Area2DInfo::Init( LONG_T start_x, LONG_T start_y, LONG_T end_x, LONG_T end_y )
{
	m_LowLeft.Set( start_x, start_y );
	m_UpRight.Set( end_x, end_y );	
}

Area2DInfo::~Area2DInfo()
{
	// printf("~Area2DInfo::~Area2DInfo() !!!\n");
}

InfoTable2D::InfoTable2D()
: m_SizeX(0), m_SizeY(0), m_dX(0), m_dY(0), m_X_count(0), m_Y_count(0), m_pTable2DMap(NULL)
{

}

void InfoTable2D::InfoTable2D_InitConstructor( LONG_T SizeX, LONG_T SizeY, LONG_T dX, LONG_T dY,
									BOOL_T bInPixels, int border )
{
	m_SizeX = SizeX;
	m_SizeY = SizeY;
	m_pTable2DMap = NULL;
		
	if(bInPixels){
		m_dX = dX;
		m_dY = dY;
		m_X_count = (m_SizeX/m_dX) + (m_SizeX % m_dX ? 1 : 0);
		m_Y_count = (m_SizeY/m_dY) + (m_SizeY % m_dY ? 1 : 0);
	
	}else{
		m_X_count = dX;
		m_Y_count = dY;
		m_dX = (m_SizeX/m_X_count);
		m_dY = (m_SizeY/m_Y_count);
		if((m_dX*m_X_count)<m_SizeX){
			m_dX++;
		}
		if((m_dY*m_Y_count)<m_SizeY){
			m_dY++;
		}
	}

	Init2DMap( border );

}

InfoTable2D::InfoTable2D( LONG_T SizeX, LONG_T SizeY, LONG_T dX, LONG_T dY, 
								  BOOL_T bInPixels, int border )
: m_SizeX(SizeX), m_SizeY(SizeY), m_pTable2DMap(NULL)
{		
	if(bInPixels){
		m_dX = dX;
		m_dY = dY;
		m_X_count = (m_SizeX/m_dX) + (m_SizeX % m_dX ? 1 : 0);
		m_Y_count = (m_SizeY/m_dY) + (m_SizeY % m_dY ? 1 : 0);
	
	}else{
		m_X_count = dX;
		m_Y_count = dY;
		m_dX = (m_SizeX/m_X_count);
		m_dY = (m_SizeY/m_Y_count);
		if((m_dX*m_X_count)<m_SizeX){
			m_dX++;
		}
		if((m_dY*m_Y_count)<m_SizeY){
			m_dY++;
		}
	}

	Init2DMap( border );
}

void InfoTable2D::Clear()
{
	if(m_pTable2DMap){
		for(int i=0;i<m_Y_count;i++){
			delete [] m_pTable2DMap[i];
		}
		delete [] m_pTable2DMap;
		m_pTable2DMap = NULL;
	}
}

void InfoTable2D::Init2DMap( int border )
{
   m_pTable2DMap = new Area2DInfo*[m_Y_count];
   for(int i=0;i<m_Y_count;i++){
      m_pTable2DMap[i] = new Area2DInfo[m_X_count];
		for(int j=0;j<m_X_count;j++){
			long end_x = (j+1)*m_dX;
			if(end_x>m_SizeX)
				end_x = m_SizeX;

			// do not put outside j<m_X_count loop !!!
			long end_y = (i+1)*m_dY;
			if(end_y>m_SizeY)
      	      end_y = m_SizeY;


			int start_x = j*m_dX;
			int start_y = i*m_dY;
			if( j==0 ){
				start_x += border;
			}
			if( i==0 ){
				start_y += border;
			}
			if( j==(m_X_count-1) ){
				end_x -= border;
			}
			if( i==(m_Y_count-1) ){
				end_y -= border;
			}
			
			m_pTable2DMap[i][j].Init( start_x, start_y, end_x, end_y );
		}
   }
}

void InfoTable2D::Init( LONG_T SizeX, LONG_T SizeY, LONG_T dX, LONG_T dY, 
								BOOL_T bInPixels, int border )
{
	m_SizeX = SizeX;
	m_SizeY = SizeY;

	if(bInPixels){
		m_dX = dX;
		m_dY = dY;
		m_X_count = (m_SizeX/m_dX) + (m_SizeX % m_dX ? 1 : 0);
		m_Y_count = (m_SizeY/m_dY) + (m_SizeY % m_dY ? 1 : 0);
	
	}else{
		m_X_count = dX;
		m_Y_count = dY;
		m_dX = (m_SizeX/m_X_count);
		m_dY = (m_SizeY/m_Y_count);
		if((m_dX*m_X_count)<m_SizeX){
			m_dX++;
		}
		if((m_dY*m_Y_count)<m_SizeY){
			// m_dY += (LONG_T)ceil(double(m_SizeY-m_dY*m_Y_count)/double(m_Y_count));
			m_dY++;
		}
	}


	Clear();
	Init2DMap( border );
}


InfoTable2D::~InfoTable2D()
{
	Clear();
}


/*Area2DInfo& InfoTable2D::GetAreaDesc( LONG_T x_pixel, LONG_T y_pixel )
{
	LONG_T x_elem = (x_pixel/m_dX);
	LONG_T y_elem = (y_pixel/m_dY);		

	return m_pTable2DMap[y_elem][x_elem];
}*/

/*Area2DInfo& InfoTable2D::GetElem( LONG_T X_elem, LONG_T Y_elem )
{
	return m_pTable2DMap[Y_elem][X_elem];
}*/


void InfoTable2D::Dump()
{
	for(register int y=0;y<m_Y_count;y++){
		for(register int x=0;x<m_X_count;x++){
			int count = (int)m_pTable2DMap[y][x].m_DataInfo[0].m_Average;
			double aver_mag=0;
			if( count > 0 )
				aver_mag = (m_pTable2DMap[y][x].m_DataInfo[0].m_Sigma/m_pTable2DMap[y][x].m_DataInfo[0].m_Average);

			printf("%d/%.2f ",count,aver_mag);
		}
		printf("\n");
	}
}



double InfoTable2D::GetMeanBackground( long x, long y, eLaplaceType_T laptype )
{
	int laptype_int = (int)laptype;
	if( laptype_int>=0 && laptype_int<(int)eLastEnumElem)
		return m_pTable2DMap[0][0].m_DataInfo[laptype_int].m_Average;
	
	return 0.00;	
}

double InfoTable2D::GetSigmaBackground( long x, long y, eLaplaceType_T laptype,
                    				          BOOL_T& bFitOK )
{
	int laptype_int = (int)laptype;
	if( laptype_int>=0 && laptype_int<(int)eLastEnumElem){
		bFitOK = m_pTable2DMap[0][0].m_DataInfo[laptype_int].m_bFitOK;
		return m_pTable2DMap[0][0].m_DataInfo[laptype_int].m_Sigma;		
	}

	return 0.00;	
}

double InfoTable2D::GetSigmaBackground( long x, long y, eLaplaceType_T laptype )
{
	int laptype_int = (int)laptype;
	if( laptype_int>=0 && laptype_int<(int)eLastEnumElem)
		return m_pTable2DMap[0][0].m_DataInfo[laptype_int].m_Sigma;

	return 0.00;	
}

BOOL_T InfoTable2D::GetMeanAndSigmaBackgr( long x, long y, double& mean, double& sigma, eLaplaceType_T laptype )
{
	int laptype_int = (int)laptype;

	mean = 0;
	sigma = 0;
   if( laptype_int>=0 && laptype_int<(int)eLastEnumElem){
		sigma = m_pTable2DMap[0][0].m_DataInfo[laptype_int].m_Sigma;
		mean = m_pTable2DMap[0][0].m_DataInfo[laptype_int].m_Average;
		return TRUE;
	}
	return FALSE;
}