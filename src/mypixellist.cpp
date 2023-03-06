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
#include "mypixellist.h"
#include "mybits.h"
#include <string.h>
#include "cexcp.h"

// #include <stdio.h>

// BOOL_T CPixelList::m_bActive=TRUE;

CPixelList::CPixelList(long nSize):m_nPixels(nSize)
{
	m_Size = (m_nPixels/8) + 1;
	m_pList = new char[m_Size];
	Clear();	
}

CPixelList::~CPixelList()
{
	delete [] m_pList;
}

void CPixelList::Clear()
{
	memset(m_pList,0,m_Size);
}

void CPixelList::ClearPart( LONG_T* list_to_clear, LONG_T cnt )
{
	for(register int i=0;i<cnt;i++){
		// printf("%d,%d\n",i,list_to_clear[i]);fflush(0);
		int pos = (list_to_clear[i] >> 3);
		Assert(pos>=0 && pos<m_Size,"CPixelList::ClearPart Range of pixel list exceeded %d , size=%d",pos,m_Size);
		m_pList[pos]=0;
	}
}

BOOL_T CPixelList::CheckPixel(long nPos)
{
	int pos = (nPos >> 3); // >>3 ==  /8
	if(pos<0 || pos>=m_Size){
		Assert(pos>=0 && pos<m_Size,"CPixelList::CheckPixel Range of pixel list exceeded %d , size=%d",pos,m_Size);
	}
	char& EightPixels = m_pList[pos];
	long nBit = (nPos-(pos<<3));
	return CheckBit( EightPixels, nBit );
}

void CPixelList::HitPixel(long nPos)
{
	int pos = (nPos >> 3);
	if(pos<0 || pos>=m_Size){
		Assert(pos>=0 && pos<m_Size,"CPixelList::HitPixel Range of pixel list exceeded %d , size=%d",pos,m_Size);
	}
	char& EightPixel = m_pList[pos];
	long nBit = (nPos-(pos << 3));
	HitBit( EightPixel, nBit );
}


void CPixelList::ClearPixel(long nPos)
{
	int pos = (nPos >> 3);
	if(pos<0 || pos>=m_Size){
		Assert(pos>=0 && pos<m_Size,"CPixelList::ClearPixel Range of pixel list exceeded %d , size=%d",pos,m_Size);
	}
   char& EightPixel = m_pList[pos];
   long nBit = (nPos-(pos << 3));
	ClearBit( EightPixel , nBit );
}

