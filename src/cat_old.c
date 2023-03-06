#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include "asas_errors.h"

int needswab();
void swab4(), swab2();

int rec_size=10;

void coo2xy(ra,dec,ia,id)
double ra,dec;
int *ia, *id;
{
  *ia = (ra*15)/3;
  if (*ia < 0) *ia = 0;
  if (*ia > 119) *ia =119;
  *id = (dec+90)/3;
  if (*id < 0) *id = 0;
  if (*id > 59) *id = 59;
}
/*int get_cat(ra_min,ra_max,dec_min,dec_max,ra,dec,mag,ngsc,path,field)
double ra_min,ra_max,dec_min,dec_max,**ra,**dec,**mag;
int *ngsc,field;
char *path;*/

int get_cat_old( double ra_min,double ra_max,double dec_min,double dec_max,
             double** ra,double** dec,double** mag,int* ngsc,char* path,
             int field )
{
  static int offset[120][60],cnt[120][60];
  static int fdc;
  static int opened = 0;
  int id0,id1,ia0,ia1,ia0_,ia1_;
  int i,j,k,n,ira,idec;
  short imag;
  char buf[256];

  if ( !opened ){
    int fd;
    FILE *fp;

    sprintf(buf,"%s.def",path);
    fp = fopen(buf, "r");
    if(fp==NULL){
    	printf("could not open file : %s\n",buf);
    	return(E_no_cat_stars);
    }
    fscanf(fp,"%d",&rec_size);
    fclose(fp);
//
    sprintf(buf,"%s.idx",path);
    fd = open(buf, O_RDONLY);
    if (fd <= 0){
    	printf("could not open file : %s\n",buf);
    	return(E_no_cat_stars);
    }
    read(fd, offset, 60*120*4);
    read(fd, cnt, 60*120*4);
    if(needswab()){
      swab4(offset,60*120);
      swab4(cnt,60*120);
    }
    close(fd);
    sprintf(buf,"%s.bin",path);
    fdc = open(buf, O_RDONLY);
    if (fdc <= 0) {
    	printf("could not open file : %s\n",buf);
      return(E_no_cat_stars);
    }
    opened =1;
  }

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
      printf("%d %d\n",offset[i][j],cnt[i][j]);
*/
       n+=cnt[i][j];
    }
    for(i=ia0_; i<=ia1_; i++) n+=cnt[i][j];
  }
  *ngsc = n;
  *ra = (double *)malloc(n*sizeof(double));
  *dec = (double *)malloc(n*sizeof(double));
  *mag = (double *)malloc(n*sizeof(double));
  n=0;
  for(j=id0; j<=id1; j++){
    for(i=ia0; i<=ia1; i++){
      lseek(fdc, offset[i][j]*rec_size,SEEK_SET);
      for(k=0; k<cnt[i][j]; k++){
        read(fdc, &buf, rec_size);
        memcpy(&ira, &buf[0], 4);
        memcpy(&idec, &buf[4], 4);
        memcpy(&imag, &buf[8+field*2], 2);
        if(needswab()){
          swab4(&ira,1);
          swab4(&idec,1);
          swab2(&imag,1);
        }
        *(*ra  + n) = (double)ira/15./3600./1000.;
        *(*dec + n) = (double)idec/3600./1000.;
        *(*mag + n) = (double)imag/1000.;
        n++;
      }
    }
    for(i=ia0_; i<=ia1_; i++){
      lseek(fdc, offset[i][j]*rec_size,SEEK_SET);
      for(k=0; k<cnt[i][j]; k++){
        read(fdc, &buf, rec_size);
        memcpy(&ira, &buf[0], 4);
        memcpy(&idec, &buf[4], 4);
        memcpy(&imag, &buf[8+field*2], 2);
        if(needswab()){
          swab4(&ira,1);
          swab4(&idec,1);
          swab2(&imag,1);
        }
        *(*ra  + n) = (double)ira/15./3600./1000.;
        *(*dec + n) = (double)idec/3600./1000.;
        *(*mag + n) = (double)imag/1000.;
        n++;
      }
    }
  }
  return(0);
}
int get_next_gsc(fd,ra,dec,mag,field)
int fd, field;
double *ra,*dec,*mag;
{
  int ira,idec;
  short imag;
  char buf[256];

  read(fd, &buf, rec_size);
  memcpy(&ira, &buf[0], 4);
  memcpy(&idec, &buf[4], 4);
  memcpy(&imag, &buf[8+field*2], 2);

  if(needswab()){
    swab4(&ira,1);
    swab4(&idec,1);
    swab2(&imag,1);
  }
  *ra  = (double)ira/15./3600./1000.;
  *dec = (double)idec/3600./1000.;
  *mag = (double)imag/1000.;
  return(0);
}

