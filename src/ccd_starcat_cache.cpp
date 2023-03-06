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
#include <myutil.h>
#include "ccd_globals.h"
#include <unistd.h>
#include <myfile.h>
#include <mystrtable.h>
#include "ccd_starcat.h"
#include <AstroCCD.h>
#include <Astroangle.h>
#include <vector>
#include "ccd_starcat_cache.h"

using namespace std;

// read :
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


extern "C" {
void  swab4(char * cp,int n);
void  swab2(char * cp,int n);
int  needswab();
// void coo2xy(double ra,double dec,int* ia,int* id);
}

int CStarCatCache::m_MagFieldIndex=0;

void coo2xy(double ra,double dec,int* ia,int* id)
{
  *ia = (ra*15)/3;
  if (*ia < 0) *ia = 0;
  if (*ia > 119) *ia =119;
  *id = (dec+90)/3;
  if (*id < 0) *id = 0;
  if (*id > 59) *id = 59;
}


CStarCatCache::CStarCatCache( const char* szCatBase )
{
	m_szCatBase="/opt/pi/dev/pisys/daq/ndir/data/cat/act";

	if( szCatBase && szCatBase[0] ){
		m_szCatBase = szCatBase;
	}
	
	m_bInitialized=FALSE;
	record_size=14;
}


int CStarCatCache::Init()
{
	printf("Initializing starcat cache from %s\n",m_szCatBase.c_str());fflush(stdout);

	mystring szCatFile,szDefFile,szIdxFile;
	szCatFile << m_szCatBase.c_str() << ".bin";
	szDefFile << m_szCatBase.c_str() << ".def";
	szIdxFile << m_szCatBase.c_str() << ".idx";

	FILE *fp;
                                                                                
   fp = fopen(szDefFile.c_str(), "r");
   if(fp==NULL){
		printf("could not open def file : %s\n",szDefFile.c_str());
		exit(-1);
	}
   fscanf(fp,"%d",&record_size);
   fclose(fp);

	int fd = open(szIdxFile.c_str(), O_RDONLY);
   if (fd <= 0){
		printf("could not open idx file : %s\n",szIdxFile.c_str());
		exit(-1);
	}
   read(fd, offset, 60*120*4);
   read(fd, cnt, 60*120*4);
   if(needswab()){
   	swab4((char*)offset,60*120);
   	swab4((char*)cnt,60*120);
   }
   close(fd);

	
	int file_size = MyFile::GetFileSize( szCatFile.c_str() );
	int stars_count = file_size/record_size;	
	
	char* buffer = new char[file_size];

	int fdc = open( szCatFile.c_str(), O_RDONLY);
   if (fdc <= 0) {
		printf("Could not open star catalog file : %s\n",szCatFile.c_str());
		exit(-1);
//		return 0;
   }

	int read_ret = read( fdc, buffer, file_size );
	if( read_ret < file_size ){
		printf("error reading file : %s\n",szCatFile.c_str());
		printf("read %d bytes of expected %d\n",read_ret,file_size);
		close(fdc);
		return 0;
	}
	close( fdc );

	int field = CStarCatCache::m_MagFieldIndex; // VISUAL MAGNITUDO
	int ira,idec;
   short imag;

	CCatalogStar star;
	star.vt = 0;
	star.bt = 0;
	star.bv = 0;
	for(register int i=0;i<stars_count;i++){
		int pos = i*record_size;
		
		memcpy(&ira, &buffer[pos], 4);
      memcpy(&idec, &buffer[pos+4], 4);
      memcpy(&imag, &buffer[pos+8+field*2], 2);		
		if(needswab()){
      	swab4( (char*)(&ira),1);
         swab4( (char*)(&idec),1);
         swab2( (char*)(&imag),1);
      }
	
		star.ra = (double)ira/15./3600./1000.;
		star.dec = (double)idec/3600./1000.;
		star.mag = (double)imag/1000.;

		star_catalog.push_back( star );
	}	

	delete [] buffer;

	m_bInitialized=TRUE;
	return star_catalog.size();
}


int CStarCatCache::get_cat( double ra_min,double ra_max,double dec_min,double dec_max,
                				 vector<CCatalogStar>& star_list )
{
  int id0,id1,ia0,ia1,ia0_,ia1_;
  int i,j,k,n,ira,idec;

	star_list.clear();

	if(!m_bInitialized ){
		int nStars = Init();
		printf("Cache catalog initialized %d stars\n",nStars);
		if( nStars<=0 ){
			return 0;
		}
	}

  if(ra_min < 0) ra_min += 24.;
  if(ra_min >24) ra_min -= 24.;
  if(ra_max < 0) ra_max += 24.;
  if(ra_max >24) ra_max -= 24.;

  coo2xy(ra_min,dec_min,&ia0,&id0);
  coo2xy(ra_max,dec_max,&ia1,&id1);

  printf("%d %d %d %d\n",ia0,id0,ia1,id1);

  if ( ia1 < ia0 ){
    ia0_ = 0;
    ia1_ = ia1;
    ia1  = 119;
  }else{
    ia0_ = 0;
    ia1_ = -1;
  }
  if ( id1 < id0 ){
    i=id0;
    id0=id1;
    id1=i;
  }

  n=0;
  for(j=id0; j<=id1; j++){
    for(i=ia0; i<=ia1; i++){
/*
      printf("%d %d\n",offset[i][j],cnt[i][j]);
*/
       n+=cnt[i][j];
    }
    for(i=ia0_; i<=ia1_; i++) n+=cnt[i][j];
  }

	printf("n=%d\n",n);

  n=0;
  for(j=id0; j<=id1; j++){
    for(i=ia0; i<=ia1; i++){
      for(k=0; k<cnt[i][j]; k++){		  
		  star_list.push_back( star_catalog[offset[i][j]+k] );
        n++;
      }
    }
    for(i=ia0_; i<=ia1_; i++){
      for(k=0; k<cnt[i][j]; k++){
		  star_list.push_back( star_catalog[offset[i][j]+k] );
        n++;
      }
    }
  }


	return n;
}

int CStarCatCache::get_cat( double ra, double dec, double radius_arcsec,
				                vector<CCatalogStar>& star_list )
{
  int id0,id1,ia0,ia1,ia0_,ia1_;
  int i,j,k,n,ira,idec;

	star_list.clear();

	if(!m_bInitialized ){
		int nStars = Init();
		printf("Cache catalog initialized %d stars\n",nStars);
		if( nStars<=0 ){
			return 0;
		}
	}

	double radius_in_h = (radius_arcsec/3600.00)/15.00;
   double radius_in_deg = (radius_arcsec/3600.00);

//ms change - 20071109 , due to fact that : genimage 15.85 -72.03 -pixscale=36.6
// did not generate good image , but only part of it ( a triangle ) 
// now I changed it sa all stars in DEC slice are selected 
// and distance is changed here for every star :
	double ra_min = 0;
	double ra_max = 24;

//   double ra_min = ra - radius_in_h;
//   double ra_max = ra + radius_in_h;

// THIS IS NOT IDEAL AND SHOULD BE IMPROVED BY PROPER FORMULA 
/*	if( (fabs(dec)+radius_in_deg)>45.00 ){
		printf("Far from equator ra_min=%.2f,ra_max=%.2f changed to 0,24\n",ra_min,ra_max);
		ra_min = 0;
		ra_max = 24;
	}*/
   double dec_min = dec - radius_in_deg;
   double dec_max = dec + radius_in_deg;


  if(ra_min < 0) ra_min += 24.;
  if(ra_min >24) ra_min -= 24.;
  if(ra_max < 0) ra_max += 24.;
  if(ra_max >24) ra_max -= 24.;

  coo2xy(ra_min,dec_min,&ia0,&id0);
  coo2xy(ra_max,dec_max,&ia1,&id1);

//  printf("%d %d %d %d\n",ia0,id0,ia1,id1);

  if ( ia1 < ia0 ){
    ia0_ = 0;
    ia1_ = ia1;
    ia1  = 119;
  }else{
    ia0_ = 0;
    ia1_ = -1;
  }
  if ( id1 < id0 ){
    i=id0;
    id0=id1;
    id1=i;
  }

  n=0;
  for(j=id0; j<=id1; j++){
    for(i=ia0; i<=ia1; i++){
/*
      printf("%d %d\n",offset[i][j],cnt[i][j]);
*/
       n+=cnt[i][j];
    }
    for(i=ia0_; i<=ia1_; i++) n+=cnt[i][j];
  }

//	printf("n=%d\n",n);

  n=0;
  for(j=id0; j<=id1; j++){
    for(i=ia0; i<=ia1; i++){
      for(k=0; k<cnt[i][j]; k++){		  
		  CCatalogStar& cache_star = star_catalog[offset[i][j]+k];
		  if( AstroCCD::CalcDistInDeg( AstroAngle::hours2deg( cache_star.ra ),
												 cache_star.dec,
												 AstroAngle::hours2deg( ra ),
												 dec ) <= radius_in_deg ){ 	
			  star_list.push_back( star_catalog[offset[i][j]+k] );
   	     n++;
		  }
      }
    }
    for(i=ia0_; i<=ia1_; i++){
      for(k=0; k<cnt[i][j]; k++){
		  CCatalogStar& cache_star = star_catalog[offset[i][j]+k];
		  if( AstroCCD::CalcDistInDeg( AstroAngle::hours2deg( cache_star.ra ),
												 cache_star.dec,
												 AstroAngle::hours2deg( ra ),
												 dec ) <= radius_in_deg ){ 	
			  star_list.push_back( star_catalog[offset[i][j]+k] );
   	     n++;
		  }
      }
    }
  }


	return n;

}

int CStarCatCache::get_cat_fast( double ra, double dec, double radius_arcsec,
				                vector<CCatalogStar>& star_list )
{
  int id0,id1,ia0,ia1,ia0_,ia1_;
  int i,j,k,n,ira,idec;

	star_list.clear();

	if(!m_bInitialized ){
		int nStars = Init();
		printf("Cache catalog initialized %d stars\n",nStars);
		if( nStars<=0 ){
			return 0;
		}
	}

	double radius_in_h = (radius_arcsec/3600.00)/15.00;
   double radius_in_deg = (radius_arcsec/3600.00);

	// this function is special for fast queries - until proper formula is 
	// calulated - it works only for small radiuses !
	// IT DOES NOT WORK PERFECTLY NEAR POLE !!!
   double ra_min = ra - radius_in_h;
   double ra_max = ra + radius_in_h;


   double dec_min = dec - radius_in_deg;
   double dec_max = dec + radius_in_deg;


  if(ra_min < 0) ra_min += 24.;
  if(ra_min >24) ra_min -= 24.;
  if(ra_max < 0) ra_max += 24.;
  if(ra_max >24) ra_max -= 24.;

  coo2xy(ra_min,dec_min,&ia0,&id0);
  coo2xy(ra_max,dec_max,&ia1,&id1);

//  printf("%d %d %d %d\n",ia0,id0,ia1,id1);

  if ( ia1 < ia0 ){
    ia0_ = 0;
    ia1_ = ia1;
    ia1  = 119;
  }else{
    ia0_ = 0;
    ia1_ = -1;
  }
  if ( id1 < id0 ){
    i=id0;
    id0=id1;
    id1=i;
  }

  n=0;
  for(j=id0; j<=id1; j++){
    for(i=ia0; i<=ia1; i++){
/*
      printf("%d %d\n",offset[i][j],cnt[i][j]);
*/
       n+=cnt[i][j];
    }
    for(i=ia0_; i<=ia1_; i++) n+=cnt[i][j];
  }

//	printf("n=%d\n",n);

  n=0;
  for(j=id0; j<=id1; j++){
    for(i=ia0; i<=ia1; i++){
      for(k=0; k<cnt[i][j]; k++){		  
		  CCatalogStar& cache_star = star_catalog[offset[i][j]+k];
		  if( AstroCCD::CalcDistInDeg( AstroAngle::hours2deg( cache_star.ra ),
												 cache_star.dec,
												 AstroAngle::hours2deg( ra ),
												 dec ) <= radius_in_deg ){ 	
			  star_list.push_back( star_catalog[offset[i][j]+k] );
   	     n++;
		  }
      }
    }
    for(i=ia0_; i<=ia1_; i++){
      for(k=0; k<cnt[i][j]; k++){
		  CCatalogStar& cache_star = star_catalog[offset[i][j]+k];
		  if( AstroCCD::CalcDistInDeg( AstroAngle::hours2deg( cache_star.ra ),
												 cache_star.dec,
												 AstroAngle::hours2deg( ra ),
												 dec ) <= radius_in_deg ){ 	
			  star_list.push_back( star_catalog[offset[i][j]+k] );
   	     n++;
		  }
      }
    }
  }


	return n;

}


int CStarCatCache::count_cat_stars( double ra_min,double ra_max,double dec_min,double dec_max )
{
  int id0,id1,ia0,ia1,ia0_,ia1_;
  int i,j,k,n,ira,idec;

	if(!m_bInitialized ){
		int nStars = Init();
		printf("Cache catalog initialized %d stars\n",nStars);
		if( nStars<=0 ){
			return 0;
		}
	}

  if(ra_min < 0) ra_min += 24.;
  if(ra_min >24) ra_min -= 24.;
  if(ra_max < 0) ra_max += 24.;
  if(ra_max >24) ra_max -= 24.;

  coo2xy(ra_min,dec_min,&ia0,&id0);
  coo2xy(ra_max,dec_max,&ia1,&id1);

  printf("%d %d %d %d\n",ia0,id0,ia1,id1);

  if ( ia1 < ia0 ){
    ia0_ = 0;
    ia1_ = ia1;
    ia1  = 119;
  }else{
    ia0_ = 0;
    ia1_ = -1;
  }
  if ( id1 < id0 ){
    i=id0;
    id0=id1;
    id1=i;
  }

  n=0;
  for(j=id0; j<=id1; j++){
    for(i=ia0; i<=ia1; i++){
/*
      printf("%d %d\n",offset[i][j],cnt[i][j]);
*/
       n+=cnt[i][j];
    }
    for(i=ia0_; i<=ia1_; i++) n+=cnt[i][j];
  }

	printf("n=%d\n",n);

  n=0;
  for(j=id0; j<=id1; j++){
    for(i=ia0; i<=ia1; i++){
      for(k=0; k<cnt[i][j]; k++){		  
		  if( star_catalog[offset[i][j]+k].ra>=ra_min && star_catalog[offset[i][j]+k].ra<=ra_max &&
				star_catalog[offset[i][j]+k].dec>=dec_min && star_catalog[offset[i][j]+k].dec<=dec_max ){
	        n++;
			}
      }
    }
    for(i=ia0_; i<=ia1_; i++){
      for(k=0; k<cnt[i][j]; k++){
		   if( star_catalog[offset[i][j]+k].ra>=ra_min && star_catalog[offset[i][j]+k].ra<=ra_max &&
			 	 star_catalog[offset[i][j]+k].dec>=dec_min && star_catalog[offset[i][j]+k].dec<=dec_max ){
	        n++;
			}
      }
    }
  }


	return n;
}
