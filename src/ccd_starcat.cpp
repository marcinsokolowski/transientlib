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
#include "ccd_starcat.h"
extern "C" {
#include "asas_astrometry.h"
#include <hip.h>
#include <cat_new.h>
}
#include <Astroangle.h>
#include <AstroCCD.h>
#include "ccd_globals.h"
#include <cexcp.h>

// char CStarCat::m_sStarCat[1024];
// BOOL_T CStarCat::m_bCacheStarCatalog=TRUE;

// default star catalog ( usually ASAS astrometric ):
CStarCat gStarCatDefault;
CStarCat gStarCatTYCHO( "/opt/pi/dev/pisys/daq/ndir/data/cat/tycho" ); // TYCHO catalog 

void CStarCat::Initialize()
{
//	strcpy( m_sStarCat, gCCDParams.m_szASASStarCatalog.c_str() );
}

CStarCat::CStarCat( const char* szCatBase )
: m_bCacheStarCatalog(TRUE), m_StarCatCache(szCatBase)
{}

CStarCat::~CStarCat()
{}

int CStarCat::getStarList( double ra, double dec, double radius,
              	            vector<CCatalogStar>& starList,
									BOOL_T bCheck,
									double min_mag, double max_mag )
{	
	vector<CCatalogStar> tmp_list;

	if( CStarCat::m_bCacheStarCatalog ){
		starList.clear();
	
		double ra_h = AstroAngle::rad2hours( ra );
		double dec_deg = AstroAngle::rad2deg( dec );
		double radius_deg = AstroAngle::rad2deg( radius );
		double radius_h = AstroAngle::rad2hours( radius );
		double radius_arcsec = radius_deg*3600.00;

		m_StarCatCache.get_cat( ra_h, dec_deg, radius_arcsec, tmp_list );
	}else{
		getStarList_NoCache( ra, dec, radius, tmp_list, bCheck );
	}

	
	int ret=0;
	for(vector<CCatalogStar>::iterator i=tmp_list.begin();i!=tmp_list.end();i++){	
		if( i->mag >= min_mag && i->mag <= max_mag ){
			ret++;
			starList.push_back( *i );
		}
	}

	return ret;
}

// FAST VERSION - BUT NOT CORRECT NEAR POLE !!!
// GOOD FOR SMALL RADIUSES - to find single star radius<1'
int CStarCat::getStarListFAST( double ra, double dec, double radius,
              	            vector<CCatalogStar>& starList,
									BOOL_T bCheck,
									double min_mag, double max_mag )
{	
	vector<CCatalogStar> tmp_list;

	if( CStarCat::m_bCacheStarCatalog ){
		starList.clear();
	
		double ra_h = AstroAngle::rad2hours( ra );
		double dec_deg = AstroAngle::rad2deg( dec );
		double radius_deg = AstroAngle::rad2deg( radius );
		double radius_h = AstroAngle::rad2hours( radius );
		double radius_arcsec = radius_deg*3600.00;

		m_StarCatCache.get_cat_fast( ra_h, dec_deg, radius_arcsec, tmp_list );
	}else{
		getStarList_NoCache( ra, dec, radius, tmp_list, bCheck );
	}

	
	int ret=0;
	for(vector<CCatalogStar>::iterator i=tmp_list.begin();i!=tmp_list.end();i++){	
		if( i->mag >= min_mag && i->mag <= max_mag ){
			ret++;
			starList.push_back( *i );
		}
	}

	return ret;
}


int CStarCat::getStarList_NoCache( double ra, double dec, double radius,
              	            vector<CCatalogStar>& starList,
									BOOL_T bCheck )
{
	starList.clear();
	
	int nGscStars=0;
	double* pListRA=NULL;
	double* pListDEC=NULL;
	double* pListMAG=NULL;
	
	double ra_h = AstroAngle::rad2hours( ra );
	double dec_deg = AstroAngle::rad2deg( dec );
	double radius_deg = AstroAngle::rad2deg( radius );
	double radius_h = AstroAngle::rad2hours( radius );

	double ra_min = ra_h - radius_h;
	double ra_max = ra_h + radius_h;
	double dec_min = dec_deg - radius_deg;
	double dec_max = dec_deg + radius_deg;
	int field=0;
	char sStarCat[1024];
	strcpy( sStarCat, gCCDParams.m_szASASStarCatalog.c_str() );

	int ret = get_cat( ra_min, ra_max, dec_min, dec_max, 
							 &pListRA, &pListDEC, &pListMAG, &nGscStars,
							 sStarCat, field );
	if(ret==0){
		for(register int i=0;i<nGscStars;i++){
			double ra_rad = AstroAngle::hours2rad( pListRA[i] );
			double dec_rad = AstroAngle::deg2rad( pListDEC[i] );

			double dist = AstroAngle::getDist( ra, dec, ra_rad, dec_rad );
			// double dist = sqrt( (ra_rad-ra)*(ra_rad-ra) + (dec_rad-dec)*(dec_rad-dec) );

			if(dist<=radius || !bCheck){
				CCatalogStar star;
				star.ra = pListRA[i];
				star.dec = pListDEC[i];
  				star.mag = pListMAG[i];
				star.vt = 0;
		   	star.bt = 0;
  				star.bv = 0;
		
				starList.push_back( star );
			}
		}
	}else{
		printf("Could not read star catalog from file : %s\n",sStarCat);
		printf("Check parameter CCD_ASAS_STAR_CAT\n");
		printf("Error : %d\n",ret);
		Assert(FALSE,"get_cat error : %d",ret);
	}

	if( pListRA ){
		free( pListRA );
	}
	if( pListDEC ){
		free( pListDEC );
	}
	if( pListMAG ){
		free( pListMAG );
	}

	return starList.size();							
}




int CStarCat::findClosestStar( vector<CCatalogStar>& starList, 
                               double ra, double dec, CCatalogStar& star,
										 double& distInRad )
{
	double minDist=1000000.00;
	vector<CCatalogStar>::iterator i;
	CCatalogStar* pStar=NULL;
	for(i=starList.begin();i!=starList.end();i++){
		double ra_in_rad = AstroAngle::hours2rad( i->ra );
      double dec_in_rad = AstroAngle::deg2rad( i->dec );

		
		double d = AstroAngle::getDist( ra , dec , ra_in_rad , dec_in_rad );
		//double d = sqrt( (ra-ra_in_rad)*(ra-ra_in_rad)+  
		//					  (dec-dec_in_rad)*(dec-dec_in_rad) );

		if( d<minDist ){
			minDist = d;
			pStar = &(*i);
		}
	}
	distInRad = minDist;
	if( pStar ){
		star = (*pStar);
		return 1;
	}
	return 0;
}


int CStarCat::getStarListHIP( double ra, double dec, double radius,
              	            vector<CCatalogStar>& starList,
									BOOL_T bCheck, BOOL_T bClearList )
{
/*	if( bClearList )
		starList.clear();
	
	int nGscStars=0;
	double* pListRA=NULL;
	double* pListDEC=NULL;
	double* pListMAG=NULL;
	
	double ra_h = AstroAngle::rad2hours( ra );
	double ra_deg = AstroAngle::rad2deg( ra );
	double dec_deg = AstroAngle::rad2deg( dec );
	double radius_deg = AstroAngle::rad2deg( radius );
	double radius_h = AstroAngle::rad2hours( radius );

	double ra_min = ra_h - radius_h;
	double ra_max = ra_h + radius_h;
	double dec_min = dec_deg - radius_deg;
	double dec_max = dec_deg + radius_deg;
	int field=0;
	char sStarCat[1024];
	strcpy( sStarCat, gCCDParams.m_szHIPStarCatalog.c_str() );

	double equinox=2000.0;
	int filter = 0;

	struct HIP *list;
	
	int ret = SelectHipparcosStars( ra_h, dec_deg, equinox, 2*radius_deg*3600.00,
											  2*radius_deg*3600.00, &list, sStarCat, filter );
	if(ret>0){
		for(register int i=0;i<ret;i++){
			double ra_rad = AstroAngle::hours2rad( list[i].ra );
			double dec_rad = AstroAngle::deg2rad( list[i].dec );

			double dist = AstroCCD::CalcDistInDeg( ra_deg, dec_deg,
									AstroAngle::hours2deg(list[i].ra), list[i].dec );
			// double dist = sqrt( (ra_rad-ra)*(ra_rad-ra) + (dec_rad-dec)*(dec_rad-dec) );

			if(dist<=radius_deg || !bCheck){
				CCatalogStar star;
				star.ra = list[i].ra;
				star.dec = list[i].dec;
  				star.mag = list[i].mag;
				star.vt = 0;
		   	star.bt = 0;
  				star.bv = 0;
		
				starList.push_back( star );
			}
		}
	}

	free( list );

	return starList.size();							 */

	printf("ERROR : not available in this version !!!\n");
	return 0;
}


double CStarCat::get_bigstar_radius( double mag, double max_radius )
{
	double ret = -107.26 * mag + 1150.00;
	return ret;
}


