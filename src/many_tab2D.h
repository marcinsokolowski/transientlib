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
#ifndef _MANY_TAB_2D_H__
#define _MANY_TAB_2D_H__

#include "tab2D.h"

template<class ARG_TYPE>
class CManyTab2D {
protected:
	Table2D<ARG_TYPE>* m_pCCD;
	long m_SizeX;
	long m_SizeY;
	long m_Count;
	double m_Average;
public :
	CManyTab2D();
	~CManyTab2D();
	CManyTab2D( long xSize, long ySize, long count);
	CManyTab2D(const CManyTab2D<ARG_TYPE>& right);
	BOOL_T Init( long xSize, long ySize, long nCCD );
	void InitFrame();
	BOOL_T ClearState();
	void SetData(ARG_TYPE val);
	
	long GetCount() const { return m_Count; }
	Table2D<ARG_TYPE>& operator[](const int pos);
	Table2D<ARG_TYPE>* GetMatrix(const int pos);
	CManyTab2D<ARG_TYPE>& operator=(const CManyTab2D<ARG_TYPE>& right);
};


#endif
