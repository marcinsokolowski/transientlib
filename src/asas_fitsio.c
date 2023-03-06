/****************************************************
* FITS - IO routines for ASAS software
*
*
load_gparam()
*
*****************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>                                /* floor()               */
#include <unistd.h>
#include <time.h>
#include <sys/types.h>                           /* time_t time()         */
#include <sys/time.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "asas_fitsio.h"
#include "asas_errors.h"
#include "asas_gparam_def.h"
/* #include "apocam.h" */
/* #include "astutil.h" */


char *shrink();
char *d2hms();
int fits_rdecomp ();

extern int verb;
/* int verb = 0; */

/*struct FITS {
  int fd;
  char *head;
  int hlines;
  float bs,bz;
  int bp;
  int naxis,nx,ny,nz;
  int compressed,blocksize,divisor,left,more;
  char *obuf;
};*/

struct FITS_HDR {
  int bitpix;
  int naxis;
  int naxis1;
  int naxis2;
  int naxis3;
  double bscale;
  double bzero;
  char object[32];
  double alpha;
  double delta;
  double itime;
  ulong apos;
  ulong dpos;
  char origin[32];
  char site[32];
  char observer[32];
  double tellong;
  double tellat;
  double telalt;
  char instrument[32];
  char camera[32];
  double gain;
  double rnoise;
  char camoptic[32];
  char filter[32];
  double focus;
  double pixscale;
  char savearea[32];
  short abinn;
  short sbinn;
  double case_temp;
  double chip_temp;
  time_t time_ut;
  double JD;
  double HJD;
  double equinox;
  double epoch;
  char date_obs[32];
  char time_obs[32];
  double ha;
  double ST;
  double airmass;
  double zenith_d;
  int badline;
  int rotate;                   /* image rotated - East(?) of mount */
  int nimage;                   /* sequential image number */
  char filename[32];
  char software[32];
  int tracking_pulse;
  double ha_corr;
  double dec_corr;
  double steps_per_sec;
  double steps_per_arcsec;
  char comment[64];
};
                

/****************************************************
* load transformation parameters from the ast-file
*****************************************************/


void set_fits_times(fh)
struct FITS_HDR *fh;
{
  struct tm *gt;
  time_t tt;
  double st, zd, ut, exptime, dtime, ra, dec, lon, lat;

  ra = fh->alpha;
  dec = fh->delta;
  lon = fh->tellong;
  lat = fh->tellat;
  exptime = fh->itime;
  tt=fh->time_ut;
  gt = gmtime(&tt);

/* set TIME_OBS & DATE_OBS */
  strftime(fh->date_obs,32,"%d/%m/%y",gt);
  strftime(fh->time_obs,32,"%H:%M:%S",gt);

/* set EQUINOX & EPOCH */
  fh->equinox = equinox(tt);
  fh->epoch = fh->equinox;

/* set AIRMASS, ZENITH_DIST, ST */
  fh->airmass = xmass(ra, dec, tt, &st, &zd, lon, lat);
  fh->ST = st;
  fh->zenith_d = zd;

/* set JD, HJD */
  ut = (double)gt->tm_hour + (double)(double)gt->tm_min/60. +
    (double)gt->tm_sec/3600.;
  dtime = exptime/2./3600.;

  fh->JD=juldate(gt->tm_year,gt->tm_mon+1,gt->tm_mday,ut+dtime);
  fh->HJD=heljuldate(gt->tm_year,gt->tm_mon+1,gt->tm_mday,ut+dtime,ra,dec);
}


/****************************************
** create simple FITS header and file
****************************************/
int create_simple_fits(filename,title,nx,ny,bp,header,hlines)
char *filename, *title;
int nx,ny,bp;
char **header;
int *hlines;
{
  int fd;
  int naxis, line;
  double bs,bz;

  if ( (fd = open(filename, O_WRONLY|O_TRUNC|O_CREAT, 0644)) <= 0){
    return(-1);
  }
 
  *header = (char *)malloc(36*80);
  *hlines = 36;
  naxis = 2;
  bs= 1.;
  bz = 0.;

  line=0;
  write_keyword(header,line++,"SIMPLE",FITS_EXACT,"T","");
  write_keyword(header,line++,"BITPIX",FITS_INT,&bp,"");
  write_keyword(header,line++,"NAXIS",FITS_INT,&naxis,"");
  write_keyword(header,line++,"NAXIS1",FITS_INT,&nx,"");
  write_keyword(header,line++,"NAXIS2",FITS_INT,&ny,"");
  write_keyword(header,line++,"BSCALE",FITS_DOUBLE,&bs,"");
  write_keyword(header,line++,"BZERO",FITS_DOUBLE,&bz,
    "real = bzero + bscale*value");
  write_keyword(header,line++,"OBJECT",FITS_STRING,title,"");
  for(; line < *hlines-1; line++){
    write_keyword(header,line,"COMMENT",FITS_STRING,"","");
  }
  write_keyword(header,line++,"END",FITS_STRING,"","");
  write_fits_header(fd, header, hlines, bp);
 
  return(fd);
}

/****************************************
** create simple FITS file using floating data
****************************************/
int create_double_fits(filename,title,data,nx,ny)
char *filename, *title;
double *data;
int nx,ny;
{
  int i, j, fd, hlines, err, bp = -32;
  char *header;
  float *fl;

  fl = (float *)malloc(nx*sizeof(float));
 
  fd = create_simple_fits(filename,title,nx,ny,bp,&header,&hlines);
  if(fd <= 0)return(E_open_fits);
  for(j=0; j<ny; j++){
    for(i=0; i<nx; i++)fl[i]=data[i+j*nx];
    err = write_fits_data(fd,fl,nx,1,bp,0);
  }
  free(fl);
  if(err != 0){
    close(fd);
    free(header);
    return(err);
  }
  close_fits(fd,nx,ny,bp);
  free(header);
  return(E_OK);
}
/****************************************
** create simple FITS file using floating data
****************************************/
int create_float_fits(filename,title,data,nx,ny)
char *filename, *title;
float *data;
int nx,ny;
{
  int fd, hlines, err, bp = -32;
  char *header;
 
  fd = create_simple_fits(filename,title,nx,ny,bp,&header,&hlines);
  if(fd <= 0)return(E_open_fits);
  err = write_fits_data(fd,data,nx,ny,bp,0);
  if(err != 0){
    close(fd);
    free(header);
    return(err);
  }
  close_fits(fd,nx,ny,bp);
  free(header);
  return(E_OK);
}
/****************************************
** create simple FITS file using floating data
****************************************/
int create_short_fits(filename,title,data,nx,ny)
char *filename, *title;
short *data;
int nx,ny;
{
  int fd, hlines, err, bp = 16;
  char *header;
 
  fd = create_simple_fits(filename,title,nx,ny,bp,&header,&hlines);
  if(fd <= 0)return(E_open_fits);
  err = write_fits_data(fd,data,nx,ny,bp,0);
  if(err != 0){
    close(fd);
    free(header);
    return(err);
  }
  close_fits(fd,nx,ny,bp);
  free(header);
  return(E_OK);
}

/****************************************
** create simple FITS file using floating data
****************************************/
int create_ushort_fits(filename,title,data,nx,ny)
char *filename, *title;
u_short *data;
int nx,ny;
{
  int fd, hlines, err, bp = -16;
  char *header;
 
  fd = create_simple_fits(filename,title,nx,ny,bp,&header,&hlines);
  if(fd <= 0)return(E_open_fits);
  err = write_fits_data(fd,data,nx,ny,bp,0);
  if(err != 0){
    close(fd);
    free(header);
    return(err);
  }
  close_fits(fd,nx,ny,bp);
  free(header);
  return(E_OK);
}

/****************************************
** open fits file for reading or writing
****************************************/
int open_fits(filename,mode)
char *filename;
int mode;
{
  char buf[256];
  int fd;

  strcpy(buf, filename);
  if (mode&(O_WRONLY|O_RDWR)){
    mode = mode | O_CREAT;
    if(strchr(filename,'.') == NULL){
      sprintf(buf,"%s.fits",filename);
    }
  }else{
    strcpy(buf, filename);
  }
  if ((fd = open(buf,mode, 0644)) == -1) {
    if(strstr(filename,".fits") == NULL){
      sprintf(buf,"%s.fitsc",filename);
      if ((fd = open(buf,mode, 0644)) == -1) {
        sprintf(buf,"%s.fits",filename);
        if ((fd = open(buf,mode, 0644)) == -1) {
          if(strstr(filename,".fts") == NULL){
            sprintf(buf,"%s.fts",filename);
            fd = open(buf,mode, 0644);
          }
        }
      }
    }
    if(fd < 0 ) return(-1);
  }
  return(fd);
}
/****************************************
** open fits file for reading or writing
****************************************/
int open_fitsc(filename,mode,pf)
char *filename;
int mode;
struct FITS *pf;
{
  char buf[256];

  strcpy(buf, filename);
  if (mode&(O_WRONLY|O_RDWR)){
    mode = mode | O_CREAT;
    if(strchr(filename,'.') == NULL){
      sprintf(buf,"%s.fits",filename);
    }
  }else{
    strcpy(buf, filename);
  }
  if ((pf->fd = open(buf,mode, 0644)) == -1) {
    if(strstr(filename,".fits") == NULL){
      sprintf(buf,"%s.fitsc",filename);
      if ((pf->fd = open(buf,mode, 0644)) == -1) {
        sprintf(buf,"%s.fits",filename);
        if ((pf->fd = open(buf,mode, 0644)) == -1) {
          if(strstr(filename,".fts") == NULL){
            sprintf(buf,"%s.fts",filename);
            pf->fd = open(buf,mode, 0644);
          }
        }
      }
    }
    if(pf->fd < 0 ) return(-1);
  }
  pf->head=NULL;
  pf->hlines=0;
  pf->bs=pf->bz=0.;
  pf->bp=0;
  pf->naxis=pf->nx=pf->ny=pf->nz=0;
  pf->compressed=pf->blocksize=pf->divisor=pf->left=pf->more=0;
  pf->obuf = NULL;

  return(pf->fd);
}
/****************************************
** open fits file for reading or writing
****************************************/
int open_fitsc_keep(filename,mode,pf)
char *filename;
int mode;
struct FITS *pf;
{
  char buf[256];

  strcpy(buf, filename);
  if (mode&(O_WRONLY|O_RDWR)){
    mode = mode | O_CREAT;
    if(strchr(filename,'.') == NULL){
      sprintf(buf,"%s.fits",filename);
    }
  }else{
    strcpy(buf, filename);
  }
  if ((pf->fd = open(buf,mode, 0644)) == -1) {
    if(strstr(filename,".fits") == NULL){
      sprintf(buf,"%s.fitsc",filename);
      if ((pf->fd = open(buf,mode, 0644)) == -1) {
        sprintf(buf,"%s.fits",filename);
        if ((pf->fd = open(buf,mode, 0644)) == -1) {
          if(strstr(filename,".fts") == NULL){
            sprintf(buf,"%s.fts",filename);
            pf->fd = open(buf,mode, 0644);
          }
        }
      }
    }
    if(pf->fd < 0 ) return(-1);
  }
//  pf->head=NULL;
//  pf->hlines=0;
//  pf->bs=pf->bz=0.;
//  pf->bp=0;
//  pf->naxis=pf->nx=pf->ny=pf->nz=0;
  pf->compressed=pf->blocksize=pf->divisor=pf->left=pf->more=0;
  pf->obuf = NULL;

  return(pf->fd);
}
/***************************************
** duplicate fitsc description file contents (header, sizes, etc)
****************************************/
int dup_fitsc(pfd, pfs)
struct FITS *pfd, *pfs;
{
  if(pfd == NULL || pfs == NULL)return(-1);
  pfd->hlines = pfs->hlines;
  if((pfd->head = (char *)calloc(pfd->hlines*80,1)) == NULL)return(-1);
  memcpy(pfd->head,pfs->head,pfd->hlines*80);
  pfd->naxis = pfs->naxis;
  pfd->bs = pfs->bs;
  pfd->bp = pfs->bp;
  pfd->bz = pfs->bz;
  pfd->nx = pfs->nx;
  pfd->ny = pfs->ny;
  pfd->nz = pfs->nz;
  pfd->compressed = 0;
  pfd->blocksize = 32;
  pfd->divisor = 1;
  pfd->left = 0;
  pfd->more = 1;
  if(pfd->obuf != NULL) free(pfd->obuf);
  pfd->obuf = NULL;
  return(0);
} 
/****************************************
** read header
****************************************/
int read_fitsc_header(pf)
struct FITS *pf;

{
  char buf[256];
  int   endflag = 0;
  int   j;
 
  pf->head = NULL;
  pf->hlines=0;
  pf->bs = 1.;
  pf->bz = 0.;
  pf->bp = 16;
  pf->naxis = 2;
  pf->nx = 0;
  pf->ny = 0;
  pf->nz = 1;
  pf->compressed = 0;
  pf->blocksize = 32;
  pf->divisor = 1;
  pf->left = 0;
  pf->more = 1;
  if(pf->obuf != NULL) free(pf->obuf);
  pf->obuf = NULL;

  while (endflag == 0) {
    if((  pf->head = (char *)realloc(  pf->head,(  pf->hlines+36)*80)) == NULL){
      return(E_malloc);
    }
    for (j=0; j<36; j++) {
      if(read(pf->fd,buf,80) != 80){
        close(pf->fd);
        free(  pf->head);
          pf->head = NULL;
          pf->hlines =0;
        return(E_read_fits);
      }
      memcpy((char *)(  pf->head+80*((  pf->hlines)++)),buf,80);
      buf[80]='\0';
      if (!strncmp(buf,"BSCALE",6)) {          /* get BSCALE    */
        sscanf(buf+10,"%f",&  pf->bs);
      } else if (!strncmp(buf,"BZERO",5)) {    /* get BZERO     */
        sscanf(buf+10,"%f",&  pf->bz);
      } else if (!strncmp(buf,"BITPIX",6)) {   /* get BITPIX    */
          pf->bp = atoi(buf+10);
      } else if (!strncmp(buf,"NAXIS1",6)) {   /* get NAXIS1    */
          pf->nx = atoi(buf+10);
      } else if (!strncmp(buf,"NAXIS2",6)) {   /* get NAXIS2    */
          pf->ny = atoi(buf+10);
      } else if (!strncmp(buf,"NAXIS3",6)) {   /* get NAXIS3    */
          pf->nz = atoi(buf+10);
      } else if (!strncmp(buf,"NAXIS",5)) {   /* get NAXIS    */
          pf->naxis = atoi(buf+10);
      } else if (!strncmp(buf,"COMPRES",7)) {   /* get COMPRESS  */
          pf->compressed = 1;
      } else if (!strncmp(buf,"BLOCKSZ",7)) {   /* get BLOCKSIZE  */
          pf->blocksize = atoi(buf+10);
      } else if (!strncmp(buf,"DIVISOR",7)) {   /* get DIVISOR  */
          pf->divisor = atoi(buf+10);
      } else if (!strncmp(buf,"END",3)) {      /* end of header */
        endflag = 1;
      }
    }
  }
  return(0);
} 
/****************************************
** open fits file and read header; return file descriptor
****************************************/
int read_fits_header(fd,header,hlines,naxis1,naxis2,bitpix,bscale,bzero)
int   fd;
char **header;
int *hlines;
int *naxis1,*naxis2;
int *bitpix;
float *bscale,*bzero;

{
  char buf[256];
  int   endflag = 0;
  int   j;
 
  *header=NULL;
  *hlines=0;
  *bscale = 1.;
  *bzero = 0.;
  *bitpix = 16;
  *naxis1 = 0;
  *naxis2 = 0;
  while (endflag == 0) {
    if((*header = (char *)realloc(*header,(*hlines+36)*80)) == NULL){
      return(E_malloc);
    }
    for (j=0; j<36; j++) {
      if(read(fd,buf,80) != 80){
        close(fd);
        free(*header);
        *header = NULL;
        *hlines =0;
        return(E_read_fits);
      }
      memcpy((char *)(*header+80*((*hlines)++)),buf,80);
      buf[80]='\0';
      if (!strncmp(buf,"BSCALE",6)) {          /* get BSCALE    */
        sscanf(buf+10,"%f",bscale);
      } else if (!strncmp(buf,"BZERO",5)) {    /* get BZERO     */
        sscanf(buf+10,"%f",bzero);
      } else if (!strncmp(buf,"BITPIX",6)) {   /* get BITPIX    */
        *bitpix = atoi(buf+10);
      } else if (!strncmp(buf,"NAXIS1",6)) {   /* get NAXIS1    */
        *naxis1 = atoi(buf+10);
      } else if (!strncmp(buf,"NAXIS2",6)) {   /* get NAXIS2    */
        *naxis2 = atoi(buf+10);
      } else if (!strncmp(buf,"END",3)) {      /* end of header */
        endflag = 1;
      }
    }
  }
  return(0);
} 
/****************************************
** get fits data (file must have been opened using open_fits)
****************************************/
int read_fits_data(fd,nx,ny,bitpix,data)
int fd,nx,ny,bitpix;
void **data;
{
  int bpp=4;

  switch(bitpix){
  case 16:
  case -16:
    bpp=2;
    break;
  case 32:
  case -32:
    bpp=4;
    break;
  }

  if( (*data = (void *)malloc(nx*ny*bpp)) == NULL)return(E_malloc);
 
  read(fd,(char *)(*data),nx*ny*bpp);
  if(needswab()){
    switch(bpp){
    case 2:
      swab2(*data,nx*ny);
      break;
    case 4:
      swab4(*data,nx*ny);
      break;
    }
  }
  return(0);
}


int read_fits_data_new(fd,nx,ny,bitpix,data)
int fd,nx,ny,bitpix;
void *data;
{
  int bpp=4;
  int ret=0;

  switch(bitpix){
  case 16:
  case -16:
    bpp=2;
    break;
  case 32:
  case -32:
    bpp=4;
    break;
  }

 
  ret = read(fd,(char *)data,nx*ny*bpp);
  if(needswab()){
    switch(bpp){
    case 2:
      swab2(data,nx*ny);
      break;
    case 4:
      swab4(data,nx*ny);
      break;
    }
  }
  return ret;
}

/****************************************
** get fitsc data (file must have been opened using open_fitsc)
****************************************/
int read_fitsc_data(pf,nx,ny,data)
struct FITS *pf;
int nx,ny;
void **data;
{
  int bpp=4;

  bpp = fabs(pf->bp)/8;

  if( (*data = (void *)malloc(nx*ny*bpp)) == NULL)return(E_malloc);
 
  if(pf->compressed == 0){
    read(pf->fd,(char *)(*data),nx*ny*bpp);
    if(needswab()){
      switch(bpp){
      case 2:
        swab2(*data,nx*ny);
        break;
      case 4:
        swab4(*data,nx*ny);
        break;
      case 8:
        swab8(*data,nx*ny);
        break;
      }
    }
  } else {
    int n,j, *idata = NULL;
    void *dat;
    float bs,bz;

    if(pf->obuf == NULL){
      if ( (pf->obuf = (char *)malloc(4*nx)) == NULL){
        free(*data);
        return(-1);
      }
    }

    bs = pf->bs*pf->divisor;
    bz = pf->bz*pf->divisor;
    ny = nx*ny/pf->nx;
    nx = pf->nx;

    for(j=0; j<ny; j++){
      if(pf->more){
        n = read (pf->fd, pf->obuf+pf->left,4*nx-pf->left);
        if(n <=0 ) {
          pf->more = 0;
          n = pf->left;
        }else n += pf->left;
      }else n = pf->left;
      idata = (int *)malloc(nx*sizeof(int));
      pf->left=fits_rdecomp(pf->obuf,n,idata,nx,pf->blocksize,2);
      if(pf->left)memcpy(pf->obuf, pf->obuf+n-pf->left, pf->left);
//      dat = convert_data(idata,nx,32,bs,bz,pf->bp);
      dat = convert_data(idata,nx,32,1.,0.,pf->bp);
      memcpy(*data+j*nx*bpp,dat,nx*bpp);
      free(dat);
    }
  }
  return(0);
}

/****************************************
** read fits file
****************************************/
int read_fits(filename,header,hlines,data,nx,ny,bitpix,bscale,bzero)
char *filename;
char **header;
int *hlines;
void *data;
int *nx,*ny,*bitpix;
float *bscale, *bzero;
{
  int fd,err;

  if( (fd = open_fits(filename,O_RDWR)) < 0) return(E_open_fits);
  err = read_fits_header(fd,header,hlines,nx,ny,bitpix,bscale,bzero);
  if(!err)err = read_fits_data(fd,*nx*(*ny),*bitpix,data);
  return(err);
}


/****************************************
** convert fits data (boupix_inp -> bitpix_out)
****************************************/
void *convert_data(void* inp,int npixel,int bpi,float bscale,float bzero,int bpo)
{
  int i,bppo=16;
  void *out, *tmp;
  u_short *up, *up0;
  short   *sp, *sp0;
  long    *lp, *lp0;
  float   *fp, *fp0;
  double  *dp;
  long long   *xp;
  u_char  *cp;
 
  
  switch(bpo){
  case 16:
  case -16:
    bppo=2;
    break;
  case 32:
  case -32:
    bppo=4;
    break;
  }


  out = (void *)malloc(npixel*bppo);
 
  switch(bpi){
  case 8:
    switch(bpo){
    case 16:
      for (i=0, cp=(u_char *)inp, sp=(short *)out; i<npixel; i++)
        *sp++ = (short)rintf(bzero+(float)(*cp++)*bscale);
      break;
    case -16:
      for (i=0, cp=(u_char *)inp, up=(u_short *)out; i<npixel; i++)
        *up++ = (u_short)rintf(bzero+(float)(*cp++)*bscale);
      break;
    case 32:
      for (i=0, cp=(u_char *)inp, lp=(long *)out; i<npixel; i++)
        *lp++ = (long)rintf(bzero+(float)(*cp++)*bscale);
      break;
    case -32:
      for (i=0, cp=(u_char *)inp, fp=(float *)out; i<npixel; i++)
        *fp++ = (float)(bzero+(float)(*cp++)*bscale);
      break;
    }
    break;
  case -64:
    switch(bpo){
    case 16:
      for (i=0, dp=(double *)inp, sp=(short *)out; i<npixel; i++)
        *sp++ = (short)rint(bzero+(float)(*dp++)*bscale);
      break;
      break;
    case -16:
      for (i=0, dp=(double *)inp, up=(u_short *)out; i<npixel; i++)
        *up++ = (u_short)rint(bzero+(float)(*dp++)*bscale);
      break;
    case 32:
      for (i=0, dp=(double *)inp, lp=(long *)out; i<npixel; i++)
        *lp++ = (long)rint(bzero+(float)(*dp++)*bscale);
      break;
    case -32:
      for (i=0, dp=(double *)inp, fp=(float *)out; i<npixel; i++)
        *fp++ = (float)(bzero+(float)(*dp++)*bscale);
      break;
    }
    break;
  case 64:
    switch(bpo){
    case 16:
      for (i=0, xp=(long long *)inp, sp=(short *)out; i<npixel; i++)
        *sp++ = (short)rintf(bzero+(float)(*xp++)*bscale);
      break;
      break;
    case -16:
      for (i=0, xp=(long long *)inp, up=(u_short *)out; i<npixel; i++)
        *up++ = (u_short)rintf(bzero+(float)(*xp++)*bscale);
      break;
    case 32:
      for (i=0, xp=(long long *)inp, lp=(long *)out; i<npixel; i++)
        *lp++ = (long)rintf(bzero+(float)(*xp++)*bscale);
      break;
    case -32:
      for (i=0, xp=(long long *)inp, fp=(float *)out; i<npixel; i++)
        *fp++ = (float)(bzero+(float)(*xp++)*bscale);
      break;
    }
    break;
  case 16:
    switch(bpo){
    case 16:
      for (i=0, sp0=(short *)inp, sp=(short *)out; i<npixel; i++)
        *sp++ = (short)rintf(bzero+(float)(*sp0++)*bscale);
//      tmp = inp;
//      inp = out;
//      out = tmp;
      break;
    case -16:
      for (i=0, sp=(short *)inp, up=(u_short *)out; i<npixel; i++)
        *up++ = (u_short)rintf(bzero+(float)(*sp++)*bscale);
      break;
    case 32:
      for (i=0, sp=(short *)inp, lp=(long *)out; i<npixel; i++)
        *lp++ = (long)rintf(bzero+(float)(*sp++)*bscale);
      break;
    case -32:
      for (i=0, sp=(short *)inp, fp=(float *)out; i<npixel; i++)
        *fp++ = (float)(bzero+(float)(*sp++)*bscale);
      break;
    }
    break;
  case -16:
    switch(bpo){
    case 16:
      for (i=0, up=(u_short *)inp, sp=(short *)out; i<npixel; i++)
        *sp++ = (short)rintf(bzero+(float)(*up++)*bscale);
      break;
    case -16:
      for (i=0, up0=(u_short *)inp, up=(u_short *)out; i<npixel; i++)
        *up++ = (u_short)rintf(bzero+(float)(*up0++)*bscale);
//      tmp = inp;
//      inp = out;
//      out = tmp;
      break;
    case 32:
      for (i=0, up=(u_short *)inp, lp=(long *)out; i<npixel; i++)
        *lp++ = (long)rintf(bzero+(float)(*up++)*bscale);
      break;
    case -32:
      for (i=0, up=(u_short *)inp, fp=(float *)out; i<npixel; i++)
        *fp++ = (float)(bzero+(float)(*up++)*bscale);
      break;
    }
    break;
  case 32:
    switch(bpo){
    case 16:
      for (i=0, lp=(long *)inp, sp=(short *)out; i<npixel; i++)
        *sp++ = (short)rintf(bzero+(float)(*lp++)*bscale);
      break;
    case -16:
      for (i=0, lp=(long *)inp, up=(u_short *)out; i<npixel; i++)
        *up++ = (u_short)rintf(bzero+(float)(*lp++)*bscale);
      break;
    case 32:
      for (i=0, lp0=(long *)inp, lp=(long *)out; i<npixel; i++)
        *lp++ = (long)rintf(bzero+(float)(*lp0++)*bscale);
//      tmp = inp;
//      inp = out;
//      out = tmp;
      break;
    case -32:
      for (i=0, lp=(long *)inp, fp=(float *)out; i<npixel; i++)
        *fp++ = (float)(bzero+(float)(*lp++)*bscale);
      break;
    }
    break;
  case -32:
    switch(bpo){
    case 16:
      for (i=0, fp=(float *)inp, sp=(short *)out; i<npixel; i++)
        *sp++ = (short)rintf(bzero+(float)(*fp++)*bscale);
      break;
    case -16:
      for (i=0, fp=(float *)inp, up=(u_short *)out; i<npixel; i++)
        *up++ = (u_short)rintf(bzero+(float)(*fp++)*bscale);
      break;
    case 32:
      for (i=0, fp=(float *)inp, lp=(long *)out; i<npixel; i++)
        *lp++ = (long)rintf(bzero+(float)(*fp++)*bscale);
      break;
    case -32:
      for (i=0, fp0=(float *)inp, fp=(float *)out; i<npixel; i++)
        *fp++ = (float)(bzero+(float)(*fp0++)*bscale);
//      tmp = inp;
//      inp = out;
//      out = tmp;
      break;
    }
    break;
  }
  // free(inp);
  return(out);
}

/****************************************
** create fits file (header must be prepared)
****************************************/
int write_fits_header(fd,header,hlines,bitpix)
int fd;
char **header;
int *hlines;
int bitpix;
{
  if(bitpix != 0)replace_keyword(header,hlines,"BITPIX",FITS_INT,&bitpix,"");
  write(fd,*header,*hlines*80);

  return(0);
}
/****************************************
** create fits file (header must be prepared)
****************************************/
int write_fitsc_header(pf)
struct FITS *pf;
{
  replace_keywordc(pf,"BITPIX",FITS_INT,&pf->bp,"");
  replace_keywordc(pf,"BZERO",FITS_FLOAT,&pf->bz,"");
  replace_keywordc(pf,"BSCALE",FITS_FLOAT,&pf->bs,"");
  if(pf->naxis > 0)replace_keywordc(pf,"NAXIS",FITS_INT,&pf->naxis,"");
  if(pf->nx > 0)replace_keywordc(pf,"NAXIS1",FITS_INT,&pf->nx,"");
  if(pf->ny > 0)replace_keywordc(pf,"NAXIS2",FITS_INT,&pf->ny,"");
  if(pf->nz > 0)replace_keywordc(pf,"NAXIS3",FITS_INT,&pf->nz,"");
  write(pf->fd,pf->head,pf->hlines*80);

  return(0);
}
/****************************************
** put fits data (file must have been opened using cretae_fits)
****************************************/
int write_fits_data(int fd,void* data,int nx,int ny,int bitpix,int flip)
{
  int bpp=4, linelen;

  switch(bitpix){
  case -32:
  case  32:
    bpp=4;
    break;
  case -16:
  case  16:
    bpp=2;
    break;
  }
  linelen = nx*bpp;

  switch(flip){
  case 0:                               /* do not flip */
  default:
    if(needswab()){
      char *buf;
      int j;
 
      buf=(char *)malloc(linelen);
      for(j=0; j<ny; j++){
        memcpy(buf,(char *)data+j*linelen,linelen);
        if(bpp == 2)swab2(buf,linelen/bpp);
          else swab4(buf,linelen/bpp);
        write(fd,(char *)buf,linelen);
      }
      free(buf);
    }else{
      write(fd,(char *)data,ny*linelen);
    }
    break;
  case 1:                               /* flip LEFT<->RIGHT */
    {
      char *buf;
      long *lp0,*lp1;
      short *sp0,*sp1;
      int i,j;
 
      buf=(char *)malloc(linelen);
 
      for(j=0; j<ny; j++){
        if(bpp == 2){
          sp0=(short *)data+(j+1)*nx-1;
          sp1=(short  *)buf;
          for(i=0; i<nx; i++)   *sp1++=*sp0--;
        }else{
          lp0=(long *)data+(j+1)*nx-1;
          lp1=(long *)buf;
          for(i=0; i<nx; i++)   *lp1++=*lp0--;
        }
        if(needswab()){
          if(bpp == 2)swab2(buf,linelen/bpp);
            else swab4(buf,linelen/bpp);
        }
        write(fd,(char *)buf,linelen);
      }
      free(buf);
    }
    break;
  case 2:                               /* fliup UP<->DOWN */
    if(needswab()){
      char *buf;
      int j;
 
      buf=(char *)malloc(linelen);
      for(j=ny-1; j>-0; j--){
        memcpy(buf,(char *)data+j*linelen,linelen);
        if(bpp == 2)swab2(buf,linelen/bpp);
          else swab4(buf,linelen/bpp);
        write(fd,(char *)buf,linelen);
      }
      free(buf);
    }else{
      int i;

      for(i=ny-1; i>=0; i--){
        write(fd,(char *)data+i*linelen,linelen);
      }
    }
    break;
  case 3:                               /* flip UP<->DOWN & LEFT<->RIGHT */
    {
      void *buf;
      long *lp0,*lp1;
      short *sp0,*sp1;
      int i,j;
 
      buf=(char *)malloc(linelen);
 
      for(j=ny-1; j>0; j--){
        if(bpp == 2){
          sp0=(short *)data+(j+1)*nx-1;
          sp1=(short *)buf;
          for(i=0; i<nx; i++)   *sp1++=*sp0--;
        }else{
          lp0=(long *)data+(j+1)*nx-1;
          lp1=(long *)buf;
          for(i=0; i<nx; i++)   *lp1++=*lp0--;
        }
        if(needswab()){
          if(bpp == 2)swab2(buf,linelen/bpp);
            else swab4(buf,linelen/bpp);
        }
        write(fd,(char *)buf,linelen);
      }
      free(buf);
    }
    break;
  }
  return(0);
}
/******************************************/
/* pad data to the multiple of 2880 words */
/******************************************/
int close_fits(fd,nx,ny,bitpix)
int fd,nx,ny,bitpix;
{
  int npad,bpp;
  char *cp;

  if(bitpix == 32 || bitpix == -32)bpp=4;
    else bpp =2;
  npad = (nx*ny*bpp) % 2880;
  if(npad != 0){
    npad = 2880 - npad;
    cp = (char *)calloc(npad,1);
    write(fd,cp,npad);
    free(cp);
  }
  close(fd);
  return(0);
}

/****************************************
** write image image into fits_file
****************************************/
/*int write_fits_image(im,filename,bitpix,flip,mframe,mframes)
struct IMAGE *im;
char *filename;
int bitpix,flip,mframe,mframes;
{
  int err;
  char *header=NULL;
  int hlines=0;
  static int fd;
  struct FITS_HDR *fh;

  fh = im->fits_hdr;
 

  fh->bitpix=bitpix;
  fh->naxis=2;
  if(mframes > 1) fh->naxis=3;
  fh->naxis1=(int)im->sx/(int)fh->abinn/(int)fh->sbinn;
  fh->naxis2=(int)im->sy/(int)fh->abinn/(int)fh->sbinn;
  fh->naxis3=mframes;
  fh->bscale=1.;
  fh->bzero=0.;
  
  fill_header(&header,&hlines,im,filename);

  if(mframe == 1){
    if( (fd = open_fits(filename,O_RDWR)) < 0) return(E_open_fits);
    err = write_fits_header(fd,&header,&hlines,bitpix);
    free(header);
    if(err) return(err);
  }
  err = write_fits_data(fd,(void *)im->data,im->wx,im->wy,bitpix,1);
  if(err) return(err);
  if(mframe == mframes) err = close_fits(fd,im->wx,im->wy,bitpix);

  return(err);
}*/

/************************************************************************
* get keyword from the fits header
************************************************************************/
 
char* get_keyword(header,hlines,name)
char **header;
int *hlines;
char *name;
{
  int j,found=-1;
  char *lp = *header;
  static char content[81];
  for (j=0; j<*hlines; j++) {
      strncpy(content, lp,70);
      content[10] = '\0';
    if (!strncmp(lp,name,strlen(name))) {          /* test for "name" */
      strncpy(content, lp+10,70);
      content[20] = '\0';
      found=j;
      break;
    } else if (!strncmp(lp,"END",3)) {      /* end of header */
      found=-1;
      break;
    }
    lp += 80;
  }
  if(found < 0){
    return(NULL);
  }else{
    return(shrink(content));
  }
}
/************************************************************************
* delete keyword in the fits header
************************************************************************/
 
void delete_keyword(header,hlines,name)
char **header;
int *hlines;
char *name;
{
  int j;
  char *lp = *header;

  for (j=0; j<*hlines; j++) {
    if (!strncmp(lp,name,strlen(name))) {          /* test for "name" */
      write_keyword(header,j,"COMMENT",FITS_STRING,"","");
    } else if (!strncmp(lp,"END",3)) {      /* end of header */
      break;
    }
    lp += 80;
  }
  return;
}

/************************************************************************
* replace keyword in the fits header
************************************************************************/
 
void replace_keyword(header,hlines,name,type,value,comment)
char **header;
int *hlines;
char *name;
int type;
void *value;
char *comment;
{
  int j,found=-1;
  char *lp = *header, content[80];

  for (j=0; j<*hlines; j++) {
    if (!strncmp(lp,name,strlen(name))) {          /* test for "name" */
      strncpy(content, lp,80);
      content[79] = '\0';
      found=j;
      break;
    } else if (!strncmp(lp,"END",3)) {      /* end of header */
      found=-1;
      break;
    }
    lp += 80;
  }
  if(found < 0){
    add_keyword(header,hlines,name,type,value,comment);
    return;
  }
  write_keyword(header,found,name,type,value,comment);
}
/************************************************************************
* replace keyword in the fits header
************************************************************************/
 
void add_keyword(header,hlines,name,type,value,comment)
char **header;
int *hlines;
char *name;
int type;
void *value;
char *comment;
{
  int i,j,found=0;
  char *lp = *header;

  for (j=0; j<*hlines; j++) {
    if (!strncmp(lp,"COMMENT",7) ||	/* search for COMMENT  */
       !strncmp(lp,"        ",8) ){	/* or empty space */
      found = j;
      for(i=12; i<30; i++){
        if(*(lp+i) != ' ' ){
          found = 0;
          break;
        }
      }
    }
    else if (!strncmp(lp,"END",3)) {      /* end of header */
      found = -j;
      break;
    }
    if(found != 0)break;
    lp += 80;
  }
  if(found == 0){
    fprintf(stderr,"cannot find END keyword\n");
    return;
  } else if(found < 0){
    found = -found;
    if(found > *hlines-2){
      *header = (char *)realloc(*header,(*hlines+36)*80);
      *hlines += 36;
    }
    write_keyword(header,found,name,type,value,comment);
    for(i=found+1; i<*hlines-1; i++){
      write_keyword(header,i,"COMMENT",FITS_STRING,"","");
    }
    write_keyword(header,i,"END",FITS_EXACT,"","");
    return;
  }else{
    write_keyword(header,found,name,type,value,comment);
  }
}
/*******************************************************
** write keyword at given line
*******************************************************/
void  write_keyword(header,fline,name,type,value,comment)
char **header;
int fline;
char *name;
int type;
void *value;
char *comment;
{
  char line[80],cval[80];


  switch(type){
  case FITS_SHORT:
    sprintf(cval,"%d",*(short *)value);
    break;
  case FITS_INT:
    sprintf(cval,"%d",*(int *)value);
    break;
  case FITS_LONG:
    sprintf(cval,"%ld",*(long *)value);
    break;
  case FITS_ULONG:
    sprintf(cval,"%lu",*(ulong *)value);
    break;
  case FITS_STRING:
    if(strlen((char *)value) ) sprintf(cval,"'%s'",(char *)value);
      else strcpy(cval,"");
    break;
  case FITS_EXACT:
    if(strlen((char *)value) ) sprintf(cval,"%s",(char *)value);
      else strcpy(cval,"");
    break;
  case FITS_FLOAT:
    if(fabs(*(float *)value) < 0.001 || *(float *)value > 9999999.){
      sprintf(cval,"%12.8e",*(float *)value);
    }else{
      sprintf(cval,"%10.8f",*(float *)value);
    }
    break;
  case FITS_DOUBLE:
    if(fabs(*(double *)value) < 0.001 || *(double *)value > 9999999.){
      sprintf(cval,"%12.8e",*(double *)value);
    }else{
      sprintf(cval,"%10.8f",*(double *)value);
    }

    break;

  }
  strcpy(line,name);
  if(strcmp(name,"END") != 0){
    extend(line,' ',8);
    strcat(line,"= ");
    extend(line,' ',30-strlen(cval));
    strcat(line,cval);
    extend(line,' ',31);
    strcat(line,"/ ");
    strcat(line,comment);
  }
  extend(line,' ',80);
  strncpy(*header+fline*80,line,80);
}
void extend(s,c,l)             /* fill string 's' to length 'l'  with 'c' */
char s[],c;
int l;
{
  int i;
 
  i = strlen(s);
  while (i < l) s[i++] = c;
 
  if(i < 80) s[i] = '\0';                      /* s[l] = '\0';   */
}
char    *shrink(line)
char    *line;
{
  int i,i1=0,i2;
 
  i2 = strlen(line) - 1;
 
  while ((line[i1]==' '  || line[i1]==0x27) && i1<i2) i1++;
  while ((line[i2]==' '  || line[i2]==0x27) && i2>i1) i2--;
 
  if (i2 == i1 && (line[i1]==' '  || line[i1]==0x27) ){
    line[0]='\0';
  }else{
    for (i = 0; i <= i2-i1; i++)
    line[i] = line[i1+i];
    line[i] = '\0';
  }
  return line;
}


/************************************************************************
* get keyword from the fits header
************************************************************************/
 
char* get_keywordc(ft,name)
struct FITS *ft;
char *name;
{
  int j,found=-1;
  char *lp = ft->head;
  static char content[81];
  for (j=0; j<ft->hlines; j++) {
      strncpy(content, lp,70);
      content[10] = '\0';
    if (!strncmp(lp,name,strlen(name))) {          /* test for "name" */
      strncpy(content, lp+10,70);
      content[20] = '\0';
      found=j;
      break;
    } else if (!strncmp(lp,"END",3)) {      /* end of header */
      found=-1;
      break;
    }
    lp += 80;
  }
  if(found < 0){
    return(NULL);
  }else{
    return(shrink(content));
  }
}
/************************************************************************
* delete keyword in the fits header
************************************************************************/
 
void delete_keywordc(ft,name)
struct FITS *ft;
char *name;
{
  int j;
  char *lp = ft->head;

  for (j=0; j<ft->hlines; j++) {
    if (!strncmp(lp,name,strlen(name))) {          /* test for "name" */
      write_keywordc(ft,j,"COMMENT",FITS_STRING,"","");
    } else if (!strncmp(lp,"END",3)) {      /* end of header */
      break;
    }
    lp += 80;
  }
  return;
}

/************************************************************************
* replace keyword in the fits header
************************************************************************/
 
void replace_keywordc(ft,name,type,value,comment)
struct FITS *ft;
char *name;
int type;
void *value;
char *comment;
{
  int j,found=-1;
  char *lp = ft->head, content[80];

  for (j=0; j<ft->hlines; j++) {
    if (!strncmp(lp,name,strlen(name))) {          /* test for "name" */
      strncpy(content, lp,80);
      content[79] = '\0';
      found=j;
      break;
    } else if (!strncmp(lp,"END",3)) {      /* end of header */
      found=-1;
      break;
    }
    lp += 80;
  }
  if(found < 0){
    add_keywordc(ft,name,type,value,comment);
    return;
  }
  write_keywordc(ft,found,name,type,value,comment);
}
/************************************************************************
* replace keyword in the fits header
************************************************************************/
 
void add_keywordc(ft,name,type,value,comment)
struct FITS *ft;
char *name;
int type;
void *value;
char *comment;
{
  int i,j,found=0;
  char *lp = ft->head;

  for (j=0; j<ft->hlines; j++) {
    if (!strncmp(lp,"COMMENT",7) ||	/* search for COMMENT  */
       !strncmp(lp,"        ",8) ){	/* or empty space */
      found = j;
      for(i=12; i<30; i++){
        if(*(lp+i) != ' ' ){
          found = 0;
          break;
        }
      }
    }
    else if (!strncmp(lp,"END",3)) {      /* end of header */
      found = -j;
      break;
    }
    if(found != 0)break;
    lp += 80;
  }
  if(found == 0){
    fprintf(stderr,"cannot find END keyword\n");
    return;
  } else if(found < 0){
    found = -found;
    if(found > ft->hlines-2){
      ft->head = (char *)realloc(ft->head,(ft->hlines+36)*80);
      ft->hlines += 36;
    }
    write_keywordc(ft,found,name,type,value,comment);
    for(i=found+1; i<ft->hlines-1; i++){
      write_keywordc(ft,i,"COMMENT",FITS_STRING,"","");
    }
    write_keywordc(ft,i,"END",FITS_EXACT,"","");
    return;
  }else{
    write_keywordc(ft,found,name,type,value,comment);
  }
}
/*******************************************************
** write keyword at given line
*******************************************************/
void  write_keywordc(ft,fline,name,type,value,comment)
struct FITS *ft;
int fline;
char *name;
int type;
void *value;
char *comment;
{
  char line[80],cval[80];

  switch(type){
  case FITS_SHORT:
    sprintf(cval,"%d",*(short *)value);
    break;
  case FITS_INT:
    sprintf(cval,"%d",*(int *)value);
    break;
  case FITS_LONG:
    sprintf(cval,"%ld",*(long *)value);
    break;
  case FITS_ULONG:
    sprintf(cval,"%lu",*(ulong *)value);
    break;
  case FITS_STRING:
    if(strlen((char *)value) ) sprintf(cval,"'%s'",(char *)value);
      else strcpy(cval,"");
    break;
  case FITS_EXACT:
    if(strlen((char *)value) ) sprintf(cval,"%s",(char *)value);
      else strcpy(cval,"");
    break;
  case FITS_FLOAT:
    if(fabs(*(float *)value) < 0.001 || *(float *)value > 9999999.){
      sprintf(cval,"%12.8e",*(float *)value);
    }else{
      sprintf(cval,"%10.8f",*(float *)value);
    }
    break;
  case FITS_DOUBLE:
    if(fabs(*(double *)value) < 0.001 || *(double *)value > 9999999.){
      sprintf(cval,"%12.8e",*(double *)value);
    }else{
      sprintf(cval,"%10.8f",*(double *)value);
    }
    break;

  }
  strcpy(line,name);
  if(strcmp(name,"END") != 0){
    extend(line,' ',8);
    strcat(line,"= ");
    extend(line,' ',30-strlen(cval));
    strcat(line,cval);
    extend(line,' ',31);
    strcat(line,"/ ");
    strcat(line,comment);
  }
  extend(line,' ',80);
  strncpy(ft->head+fline*80,line,80);
}
/******************************************/
/* pad data to the multiple of 2880 words */
/******************************************/
int close_fitsc(fts)
struct FITS *fts;
{
  int npad;
  char *cp;

  npad = lseek(fts->fd,0L,SEEK_END) % 2880;
//  if(fts->bp == 32 || fts->bp == -32)bpp=4;
//    else bpp =2;
//  npad = (fts->nx*fts->ny*bpp) % 2880;
  if(npad != 0){
    npad = 2880 - npad;
    cp = (char *)calloc(npad,1);
    write(fts->fd,cp,npad);
    free(cp);
  }
  close(fts->fd);
  fts->fd = 0;
  if(fts->head != NULL)free(fts->head);
  if(fts->obuf != NULL)free(fts->obuf);
  fts->head = NULL;
  fts->obuf = NULL;
  fts->hlines = 0;
  fts->left = 0;
  fts->more = 1;
  return(0);
}
/****************************************
** put fits data (file must have been opened using cretae_fits)
****************************************/
int write_fitsc_data(struct FITS * pt, void *data,int nx,int ny,int flip)
/*struct FITS *pt;
void *data;
int nx,ny;
int flip;*/
{
  int bpp=4, linelen;
  int bitpix = pt->bp;
  int fd = pt->fd;

  switch(bitpix){
  case -32:
  case  32:
    bpp=4;
    break;
  case -16:
  case  16:
    bpp=2;
    break;
  }
  linelen = nx*bpp;

  switch(flip){
  case 0:                               /* do not flip */
  default:
    if(needswab()){
      char *buf;
      int j;
 
      buf=(char *)malloc(linelen);
      for(j=0; j<ny; j++){
        memcpy(buf,(char *)data+j*linelen,linelen);
        if(bpp == 2)swab2(buf,linelen/bpp);
          else swab4(buf,linelen/bpp);
        write(fd,(char *)buf,linelen);
      }
      free(buf);
    }else{
      write(fd,(char *)data,ny*linelen);
    }
    break;
  case 1:                               /* flip LEFT<->RIGHT */
    {
      char *buf;
      long *lp0,*lp1;
      short *sp0,*sp1;
      int i,j;
 
      buf=(char *)malloc(linelen);
 
      for(j=0; j<ny; j++){
        if(bpp == 2){
          sp0=(short *)data+(j+1)*nx-1;
          sp1=(short  *)buf;
          for(i=0; i<nx; i++)   *sp1++=*sp0--;
        }else{
          lp0=(long *)data+(j+1)*nx-1;
          lp1=(long *)buf;
          for(i=0; i<nx; i++)   *lp1++=*lp0--;
        }
        if(needswab()){
          if(bpp == 2)swab2(buf,linelen/bpp);
            else swab4(buf,linelen/bpp);
        }
        write(fd,(char *)buf,linelen);
      }
      free(buf);
    }
    break;
  case 2:                               /* fliup UP<->DOWN */
    if(needswab()){
      char *buf;
      int j;
 
      buf=(char *)malloc(linelen);
      for(j=ny-1; j>-0; j--){
        memcpy(buf,(char *)data+j*linelen,linelen);
        if(bpp == 2)swab2(buf,linelen/bpp);
          else swab4(buf,linelen/bpp);
        write(fd,(char *)buf,linelen);
      }
      free(buf);
    }else{
      int i;

      for(i=ny-1; i>=0; i--){
        write(fd,(char *)data+i*linelen,linelen);
      }
    }
    break;
  case 3:                               /* flip UP<->DOWN & LEFT<->RIGHT */
    {
      void *buf;
      long *lp0,*lp1;
      short *sp0,*sp1;
      int i,j;
 
      buf=(char *)malloc(linelen);
 
      for(j=ny-1; j>0; j--){
        if(bpp == 2){
          sp0=(short *)data+(j+1)*nx-1;
          sp1=(short *)buf;
          for(i=0; i<nx; i++)   *sp1++=*sp0--;
        }else{
          lp0=(long *)data+(j+1)*nx-1;
          lp1=(long *)buf;
          for(i=0; i<nx; i++)   *lp1++=*lp0--;
        }
        if(needswab()){
          if(bpp == 2)swab2(buf,linelen/bpp);
            else swab4(buf,linelen/bpp);
        }
        write(fd,(char *)buf,linelen);
      }
      free(buf);
    }
    break;
  }
  return(0);
}

/****************************************************
* load transformation parameters from the ast-file
*****************************************************/

int load_gparam(fts,param)
struct FITS *fts;
struct GPARAM *param;
{
  char *keyword=NULL;
  int i, nord[6]={1,3,6,10,15,21};
  char buf[256];
  double rd;

  param->xc = fts->nx/2;
  param->yc = fts->ny/2;
  for(i=0; i<21; i++) param->px[i]=param->py[i]=0.;
/* load transformation parameters */
  keyword = get_keywordc(fts,"RA");
  if(keyword == NULL){
    fprintf(stderr,"load_gparam: cannot find RA in header.\n");
    return(-1);
  }
  param->ra = strtod(keyword, NULL);

  keyword = get_keywordc(fts,"DEC");
  if(keyword == NULL){
    fprintf(stderr,"load_gparam: cannot find DEC in header.\n");
    return(-1);
  }
  param->dec = strtod(keyword, NULL);

  keyword = get_keywordc(fts,"PIXSCALE");
  if(keyword == NULL){
    fprintf(stderr,"load_gparam: cannot find PIXSCALE in header.\n");
    return(-1);
  }
  param->pixscale = strtod(keyword, NULL);

  param->fi = 0.;
  keyword = get_keywordc(fts,"POSANGLE");
  if(keyword != NULL) param->fi = strtod(keyword, NULL);
    else if(verb)printf("load_gparam: cannot find POSANGLE in header.\n");


  keyword = get_keywordc(fts,"NX");
  if(keyword != NULL) param->xc = atoi(keyword)/2;
    else if(verb) fprintf(stderr,"load_gparam: cannot find NX in header.\n");

  keyword = get_keywordc(fts,"NY");
  if(keyword != NULL) param->yc = atoi(keyword)/2;
    else if(verb)fprintf(stderr,"load_gparam: cannot find NY in header.\n");

  keyword = get_keywordc(fts,"AST_ORD");
  if(keyword == NULL){
    if(verb)printf("load_gparam: cannot find AST_ORD in header.\n");
    goto missing;
  }
  param->order = atoi(keyword);

  for(i=0; i<nord[param->order]; i++){
/* x par */
    sprintf(buf,"PAR_X_%d",i);
    keyword = get_keywordc(fts,buf);
    if(keyword == NULL){
      if(verb) printf("load_gparam: cannot find %s in header.\n", buf);
      goto missing;
    }
    param->px[i] = strtod(keyword,NULL);

/* y par */
    sprintf(buf,"PAR_Y_%d",i);
    keyword = get_keywordc(fts,buf);
    if(keyword == NULL){
      if(verb) printf("load_gparam: cannot find %s in header.\n", buf);
      goto missing;
    }
    param->py[i] = strtod(keyword,NULL);
  }
  return(0);
missing:
  param->order = 1;
  rd = param->fi * 0.017453292;
  param->px[0] = param->py[0] = 0.;
  param->px[1] = param->pixscale * cos(rd);
  param->px[2] = -param->pixscale * sin(rd);
  param->py[1] = param->pixscale * sin(rd);
  param->py[2] = param->pixscale * cos(rd);
  for(i=3; i<21; i++) param->px[i]=param->py[i]=0.;
  return(1);
}


const char* get_keyword_value( const char* fname, const char* key )
{
	const char* keyword=NULL;
	struct FITS fts;	
	if(open_fitsc(fname,O_RDONLY, &fts) <= 0){
		printf("Cannot open %s\n",fname);
		return NULL;
	}
	read_fitsc_header(&fts);
	keyword = get_keywordc(fts,key);  	
	return keyword;
}