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
#ifndef _CCDSTAR_CATALOG_H__
#define _CCDSTAR_CATALOG_H__

#include <mytypes.h>
#include <vector>
#include <starcat_base_defines.h>
#include "ccd_starcat_cache.h"

using namespace std;


// structure according to asas catalog structure :
// 14 RECORD SIZE
// RA DEC VT BT BV
//  4 4 2 2 2
/*struct CCatalogStar
{
	float ra;
	float dec;
	float mag;
	short int vt;
	short int bt;
	short int bv;
};*/


// usage :
// OLD :
//      just use functions : getStarList, findClosestStar
// NEW :
//      call :
//         CStarCat::Initialize();
//      use : getStarListFast  , findClosestStarFast

class CStarCat 
{
public :
	CStarCatCache m_StarCatCache;

//	char m_sStarCat[1024];
	void Initialize();
	
	BOOL_T m_bCacheStarCatalog;

	CStarCat( const char* szCatBase=NULL );
	~CStarCat();

	// takes radian input :	
	int getStarList( double ra, double dec, double radius,
									vector<CCatalogStar>& starList, 
									BOOL_T bCheck=FALSE, 
									double min_mag=-1000, double max_mag=1000 );	

	int getStarListFAST( double ra, double dec, double radius,
	                           vector<CCatalogStar>& starList,
                              BOOL_T bCheck=FALSE,
                              double min_mag=-1000, double max_mag=1000 );

	int getStarList_NoCache( double ra, double dec, double radius,
									vector<CCatalogStar>& starList, 
									BOOL_T bCheck=FALSE );	

	//
	int findClosestStar( vector<CCatalogStar>& starList, 
										 double ra, double dec, CCatalogStar& star,
										 double& distInRad );


	int getStarListHIP( double ra, double dec, double radius,
									vector<CCatalogStar>& starList, 
									BOOL_T bCheck=FALSE, BOOL_T bClearList=TRUE );	
										  

	// other functions :
	
	// calculating radius of bright star - should be updated for 
	// different cameras ( same as in perl script : do_flareevents.pl )
	static double get_bigstar_radius( double mag, double max_radius=600 );
};

extern CStarCat gStarCatDefault;
extern CStarCat gStarCatTYCHO;

#endif
