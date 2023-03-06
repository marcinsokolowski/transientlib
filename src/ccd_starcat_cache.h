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
#ifndef _STAR_CAT_CACHE_H__
#define _STAR_CAT_CACHE_H__

#include <mystring.h>
#include <vector>

using namespace std;

// read :
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


class CCatalogStar;

class CStarCatCache
{
public:
	CStarCatCache( const char* szCatBase=NULL );

	int get_cat( double ra_min,double ra_max,double dec_min,double dec_max,
					 vector<CCatalogStar>& star_list );
	int get_cat( double ra, double dec, double radius_arcsec,
					 vector<CCatalogStar>& star_list );
	int get_cat_fast( double ra, double dec, double radius_arcsec,
                vector<CCatalogStar>& star_list );
					 
					 
	int count_cat_stars( double ra_min,double ra_max,double dec_min,double dec_max );					 

	BOOL_T m_bInitialized;	
	mystring m_szCatBase;
	int offset[120][60];
	int cnt[120][60];
	int record_size;
	static int m_MagFieldIndex; // 0-visual (vt),1-bt
	vector<CCatalogStar> star_catalog;

	int Init();
};

#endif
