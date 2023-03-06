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
#ifndef _STAR_LIST_FILE_H__
#define _STAR_LIST_FILE_H__

#include "myfitsfile.h"
#include <ccd_common_struct.h>
#include <vector>

using namespace std;

// #define MAX_ADDITIONAL_MAG_COUNT 10

/*struct cStarCat
{

	cStarCat():
		star_ident_count(0),id_star(0),id_frm(0),hjd(0),m_pMeasureTab(NULL)
	{}

	double ra;
	double dec;
	double xc;
	double yc;
	double shape;
	double sharp;
	double sky;
	char starname[16];
	double mag;
	double mag_add[MAX_ADDITIONAL_MAG_COUNT];
	
	double mag_cat;    
	double mag_piphoto;
	int    cat_star_number;
	double cat_star_dist;
	
	// star ID :
	double hjd;
	int id_star;
	int id_frm;
	int star_ident_count;

	vector<cStarCat>* m_pMeasureTab;
};*/


class CStarListFile : public CMyFITSFile<float>
{
public :
	CStarListFile();
	
	CStarListFile( const CStarListFile& right );

	double maglimit;	
	int m_SizeX;
	int m_SizeY;
	time_t unix_time;
	double ra0;
	double dec0;
	double ra2000;
	double dec2000;
	double hjd;
	double jd;
	double xmass;	
	char   object[24];
	int nmag;
	int iDayNight;
	int dImage;
	double ast_err;
	int nApert;
	int iCamid;
	vector<cStarCat> m_StarList;		
	
};


#endif
