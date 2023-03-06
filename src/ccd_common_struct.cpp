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
#include "ccd_common_struct.h"
#include <string.h>


cStarCat::~cStarCat()
{
	if( m_pMeasureTab ){
		delete m_pMeasureTab;
	}
}

cStarCat::cStarCat( const cStarCat& right )
: m_pMeasureTab(NULL)
{
	(*this) = right;
}

cStarCat& cStarCat::operator=( const cStarCat& right )
{
	bUpdated = right.bUpdated;
	ra = right.ra;
	dec = right.dec;
	xc = right.xc;
	yc = right.yc;
	shape = right.shape;
	sharp = right.sharp;
	sky = right.sky;
	strcpy( starname, right.starname );
	mag = right.mag;
	memcpy( mag_add , right.mag_add, MAX_ADDITIONAL_MAG_COUNT*sizeof(double) );
	
	min_mag = right.min_mag;
	max_mag = right.max_mag;
	sigma_mag = right.sigma_mag;
	no_measurements = right.no_measurements;
	magsum = right.magsum;
	mag2sum = right.mag2sum;

			
	mag_cat = right.mag_cat;    
	mag_piphoto = right.mag_piphoto;
	cat_star_number = right.cat_star_number;
	cat_star_dist = right.cat_star_dist;
	
	hjd = right.hjd;
	id_star = right.id_star;
	id_frm = right.id_frm;
	star_ident_count = right.star_ident_count;

	if( m_pMeasureTab ){
		delete m_pMeasureTab;
		m_pMeasureTab = NULL;
	}
	if( right.m_pMeasureTab ){
		m_pMeasureTab = new vector<cStarMeasure>();
		for( vector<cStarMeasure>::iterator i=(right.m_pMeasureTab)->begin();
			  i!=(right.m_pMeasureTab)->end();i++){
			m_pMeasureTab->push_back( *i );
		}		
	}

	// temporary for astrometry purposes :
	x_prim = right.x_prim;
	y_prim = right.y_prim;

	return (*this);
}

