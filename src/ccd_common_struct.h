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
#ifndef _CCD_COMMON_STRUCT_CFGLIB_H__
#define _CCD_COMMON_STRUCT_CFGLIB_H__

#include <stdlib.h>
#include <vector>

using namespace std;

#define MAX_ADDITIONAL_MAG_COUNT 10
#define MAX_MAG_COUNT 5

// flags :
#define MANY_MEASUREMENTS_ON_FRAME 0x01

struct cStarMeasure
{
	cStarMeasure():star(0),hjd(0),mag(0),error(0),id_frm(0),ra(0),dec(0),
		mag_piphoto(0),x(0),y(0),new_star(0), flag(0) {
		for(int i=0;i<MAX_MAG_COUNT;i++)mag_ap[i]=0;		
	}
	
	int star;
	double hjd;
	double mag;
	double error;
	int id_frm;
	double ra;
	double dec;
	double mag_piphoto;
	double x;
	double y;
	double mag_ap[MAX_MAG_COUNT];
	int    new_star;
	char   flag;
};

struct cStarCat
{

	cStarCat():
		star_ident_count(0),id_star(0),id_frm(0),hjd(0),m_pMeasureTab(NULL),
		no_measurements(0),min_mag(0),max_mag(0),sigma_mag(0),magsum(0),
		mag2sum(0),bUpdated(0),x_prim(0),y_prim(0),cat_star_number(0)
	{}
	cStarCat( const cStarCat& right );
	~cStarCat();
	
	cStarCat& operator=( const cStarCat& right );

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
	
	double min_mag;
	double max_mag;
	double sigma_mag;
	int no_measurements;
	double magsum;
	double mag2sum;
				
	double mag_cat;    
	double mag_piphoto;
	int    cat_star_number;
	double cat_star_dist;

	// temporary for astrometry purposes :
	double x_prim;
	double y_prim;
	
	// star ID :
	double hjd;
	int id_star;
	int id_frm;
	int star_ident_count;

	int bUpdated;

	vector<cStarMeasure>* m_pMeasureTab;
};




#endif
