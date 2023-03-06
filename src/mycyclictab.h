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
#ifndef _MY_CYCLIC_TAB_H
#define _MY_CYCLIC_TAB_H


// the idea of re-puting sample is that I keep number N of 
// samples to be re-put in table , adding sample 7 overwirtes sample 0 
// - but GetFirst/GetNext functon is incorrect I suppose !

template<class ARG_TYPE>
class CCycleTab2D {
protected:
	int m_Size;
	int m_Count;
	
	
	ARG_TYPE* m_pTable;


	int m_Newest;		
	int m_Oldest;
	int m_CurrPos;
public:
	CCycleTab2D( int size );
	~CCycleTab2D();
	
	void Add( const ARG_TYPE& elem );


	ARG_TYPE* GetFirst();
	ARG_TYPE* GetNext();	
	
	int GetSize() { return m_Size; }
	int GetCount() { return m_Count; }

};






#endif
