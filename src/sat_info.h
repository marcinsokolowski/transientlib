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
#ifndef _CCD_SAT_H__
#define _CCD_SAT_H__

#include <time.h>
extern "C" {
#include "predict.h"
}

#include <mytypes.h>
#include <vector>


using namespace std;

class CSatList : public vector<satInfo>
{
public:
	CSatList(){}
	~CSatList(){}
};

class CSatInfo
{
public :
	CSatInfo();

	static void InitSateliteLibrary( const char* qth_file, const char* tle_file, const char* db_file=NULL );
	static int GetSatInfo( vector<satInfo>& satList, time_t ut_time, BOOL_T bVisibleOnly=TRUE );
	static int GetVisibleOnly( vector<satInfo>& allSatList, vector<satInfo>& visibleSatList );
};

#endif

