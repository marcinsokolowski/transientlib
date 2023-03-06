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
#include "sat_info.h"
#include <AstroCCD.h>
#include <mylock.h>

// due to globals in predict.c this is currently needed :
static CMyMutex gSatLibLock;

void CSatInfo::InitSateliteLibrary( const char* qth_file, const char* tle_file, const char* db_file )
{
	gSatLibLock.Lock();
	
	int ret = ::InitSatLib( qth_file, tle_file, db_file );
	if(ret<=0){
		printf("could not initialize satelite library !!!\n");
		printf("check TLE file : %s\n",tle_file);
		printf("check QTH file : %s\n",qth_file);
		printf("trying to use predict.tle, predict.qth files\n");
		ret = ::InitSatLib( "predict.qth", "predict.tle", "" );
		if( ret<=0 ){
			printf("CSatInfo::InitSateliteLibrary failed\n");
			exit(0);
		}
	}

	gSatLibLock.UnLock();
}

int CSatInfo::GetSatInfo( vector<satInfo>& satList, time_t ut_time, BOOL_T bVisibleOnly )
{
	satList.clear();
	gSatLibLock.Lock();

	InitSatLib( NULL, NULL, NULL );
	
	AstroCCD astro( 0.050, 0.000015, 2062, 2048, 
						 obs_geodetic.lat, obs_geodetic.lon, 
						 0, 0, 0 );

	for(int x=0;x<gSatNumber;x++){
		struct satInfo sat_info;
		CalcCurrentPositionOfSingle( x, ut_time, &sat_info );

		if(!bVisibleOnly || sat_info.visibility==SAT_VISIBLE){	
			double azim_rad_geo = AstroAngle::deg2rad( sat_info.sat_azim );
			double alt_rad = AstroAngle::deg2rad( sat_info.sat_ele );

//			double azim_rad_astro = azim_rad_geo;
			double azim_rad_astro = (PI_VALUE+azim_rad_geo);
			if(azim_rad_astro>TWO_PI_VALUE){
				azim_rad_astro = (azim_rad_astro-TWO_PI_VALUE);
			}

			astro.calculateEquatorialCoordinates( azim_rad_astro, alt_rad,		
															  ut_time, sat_info.sat_ha,
															  sat_info.sat_dec,
															  sat_info.sat_ra );	

			mystring szSatName = sat_info.sat_name;
			szSatName.replace_char(' ','_');
			strcpy( sat_info.sat_name, szSatName.c_str() );
			satList.push_back( sat_info );		
		}
	}		
	gSatLibLock.UnLock();

	return satList.size();
}


int CSatInfo::GetVisibleOnly( vector<satInfo>& allSatList, vector<satInfo>& visibleSatList )
{
	visibleSatList.clear();
		
	vector<satInfo>::iterator i;
	for(i=allSatList.begin();i!=allSatList.end();i++){
		if( i->visibility==SAT_VISIBLE ){
			visibleSatList.push_back( *i );
		}
	}
	
	return visibleSatList.size();
}
