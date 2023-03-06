#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include "asas_errors.h"
#include "starcat_base_defines.h"

int needswab();
void swab4(), swab2();

int get_cat( double ra_min,double ra_max,double dec_min,double dec_max,
             double** ra,double** dec,double** mag,int* ngsc,char* path,
             int field );

static int record_size=10;

static int offset_tab[120][60];
static int cnt_tab[120][60];
static int gCacheInitialized=0;
int gAsasCatalogUseCache=0;
struct CCatalogStar* gStarCatalogCache=NULL;

static int get_file_size( const char* filename ){
	struct stat buf;
	stat( filename, &buf );
	return (int)buf.st_size;
}

int initialize_starcat_cache( const char* cat_path )
{
	if( !gCacheInitialized ){
		char szCatFile[1024],szDefFile[1024],szIdxFile[1024];
		int field=0; // VISUAL MAGNITUDO
		int ira,idec;
	   short imag;
	   struct CCatalogStar star;
	   int file_size = 0;
	   int stars_count = 0;
		FILE *fp;
		int fd = -1;
		int fdc = -1;
		int read_ret=0;
		int i=0;


		printf("Initializing star catalog cache in memory ...\n");fflush(0);
		sprintf(szCatFile,"%s.bin",cat_path );
		sprintf(szDefFile,"%s.def",cat_path );
		sprintf(szIdxFile,"%s.idx",cat_path );

	   fp = fopen(szDefFile, "r");
   	if(fp==NULL){
			printf("could not open def file : %s\n",szDefFile);
		}
	   fscanf(fp,"%d",&record_size);
   	fclose(fp);

		fd = open(szIdxFile, O_RDONLY);
	   if (fd <= 0){
			printf("could not open idx file : %s\n",szIdxFile);
		}
   	read(fd, offset_tab, 60*120*4);
	   read(fd, cnt_tab, 60*120*4);
   	if(needswab()){
   		swab4((char*)offset_tab,60*120);
	   	swab4((char*)cnt_tab,60*120);
	   }
   	close(fd);

	
		file_size = get_file_size( szCatFile );
		stars_count = file_size/record_size;	
	
		char* buffer = malloc( file_size );

		fdc = open( szCatFile, O_RDONLY);
   	if (fdc <= 0) {
			printf("Could not open star catalog file : %s\n",szCatFile);
			return 0;
   	}

		read_ret = read( fdc, buffer, file_size );
		if( read_ret < file_size ){
			printf("error reading file : %s\n",szCatFile);
			printf("read %d bytes of expected %d\n",read_ret,file_size);
			close(fdc);
			return 0;
		}
		close( fdc );

		gStarCatalogCache = malloc( stars_count*sizeof(struct CCatalogStar) );
		star.vt = 0;
		star.bt = 0;
		star.bv = 0;
		for(i=0;i<stars_count;i++){
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


			memcpy( &(gStarCatalogCache[i]) , &star, sizeof( struct CCatalogStar ));
		}	

		free( buffer );

		
		gCacheInitialized = 1;
		printf("star catalog cache initilized\n");
		return stars_count;
	}
	return 0;
}

int get_cat_new( double ra_min,double ra_max,double dec_min,double dec_max,
             double** ra,double** dec,double** mag,int* ngsc,char* path,
             int field )
{
  if( !gAsasCatalogUseCache ){
  		return get_cat( ra_min, ra_max, dec_min, dec_max, ra, dec, mag, ngsc, path, field );
  }    	
  
  if( !gCacheInitialized ){
  		initialize_starcat_cache( path );
  }

  int id0,id1,ia0,ia1,ia0_,ia1_;
  int i,j,k,n,ira,idec;
  short imag;
  char buf[256];


  if(ra_min < 0) ra_min += 24.;
  if(ra_min >24) ra_min -= 24.;
  if(ra_max < 0) ra_max += 24.;
  if(ra_max >24) ra_max -= 24.;

  coo2xy(ra_min,dec_min,&ia0,&id0);
  coo2xy(ra_max,dec_max,&ia1,&id1);
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
      printf("%d %d\n",offset_tab[i][j],cnt_tab[i][j]);
*/
       n+=cnt_tab[i][j];
    }
    for(i=ia0_; i<=ia1_; i++) n+=cnt_tab[i][j];
  }
  *ngsc = n;
  *ra = (double *)malloc(n*sizeof(double));
  *dec = (double *)malloc(n*sizeof(double));
  *mag = (double *)malloc(n*sizeof(double));
  n=0;
  for(j=id0; j<=id1; j++){
    for(i=ia0; i<=ia1; i++){
      for(k=0; k<cnt_tab[i][j]; k++){
        *(*ra  + n) = gStarCatalogCache[offset_tab[i][j]+k].ra;
        *(*dec + n) = gStarCatalogCache[offset_tab[i][j]+k].dec;
        *(*mag + n) = gStarCatalogCache[offset_tab[i][j]+k].mag;
        n++;
      }
    }
    for(i=ia0_; i<=ia1_; i++){
      for(k=0; k<cnt_tab[i][j]; k++){
        *(*ra  + n) = gStarCatalogCache[offset_tab[i][j]+k].ra;
        *(*dec + n) = gStarCatalogCache[offset_tab[i][j]+k].dec;
        *(*mag + n) = gStarCatalogCache[offset_tab[i][j]+k].mag;
        n++;
      }
    }
  }


  return(0);
}

int read_offset_cnt( const char* cat_path, int offset_tab[120][60], int cnt_tab[120][60] )
{
	char szCatFile[512],szDefFile[512],szIdxFile[512];
	sprintf(szCatFile,"%s.bin",cat_path );
	sprintf(szIdxFile,"%s.idx",cat_path );

	int fd = open(szIdxFile, O_RDONLY);
	if (fd <= 0){
		printf("could not open idx file : %s\n",szIdxFile);
	}
   int ret = read(fd, offset_tab, 60*120*4);
	ret += read(fd, cnt_tab, 60*120*4);
   if(needswab()){
   	swab4((char*)offset_tab,60*120);
	   swab4((char*)cnt_tab,60*120);
	}
   close(fd);

	return ret;
}
