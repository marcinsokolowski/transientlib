#ifndef _ASTRO_ASAS_H_
#define _ASTRO_ASAS_H_

/***************************************************************************\
 **  transientlib : library for identification of short optical flashes 
 **						with the wide field cameras.
 **  This software was mostly written by Marcin Sokolowski ( msok@fuw.edu.pl ) 
 **	it was a substantial part of analysis performed for PHD thesis 
 **  ( Investigation of astrophysical phenomena in short time scales with "Pi of the Sky" apparatus )
 **	it can be used under the terms of GNU General Public License.
 **	In case this software is used for scientific purposes and results are
 **	published, please refer to the PHD thesis submited to astro-ph :
 **
 **		http://arxiv.org/abs/0810.1179
 ** 
 ** This file was created by Grzegorz Pojmanski from ASAS experiment, please 
 ** refer also to : 
 **      Pojmanski, G., 1997, Acta Astronomica, 47, 467. 
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


#define AST_VER "V 3.11 25/07/2001"

#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <time.h>
/* #include "asas_fitsio.h" */
#include "asas_errors.h"
#include "astutil.h"


/* Defines */
#define ASTROMETRIC_CAT "/scratch/pi1/asas/cat/act"
#define maxim(A,B) ((A)>(B)?(A):(B))
#define minim(A,B) ((A)<(B)?(A):(B))
#define SQR(A) ((A)*(A))
#define RD 0.0174532925
#define POSANG 0.0
#define NTEST 30
#define MAX_RADIUS 8.	/* maximum radius for the matching routine */
#define DFI 0.20	/* angle step for angle itarations */
#define MAX_TRY 20	/* maximum number of trying different angles */
#define MAX_ITER 30	/* number of iteration of the fitting procedure */
#define MIN_COMMON 10	/* minimum number of comon stars to start fitting*/
#define XDIM 2048
#define YDIM 2048
#define XSCALE 1024
#define MAX_LIST 1200	/* maximum number of stars for get_offset */
#define _MAX_AST_ERR_FATAL 0.6 /* maximum allowable astromeric error(pixels) */
#define _MAX_AST_ERR 0.3 /* maximum allowable astromeric error (in pixels) */
#define MAX_GET_RADIUS 1.5
#define _SHARP_MIN 0.2
#define _SHARP_MAX 4.5
#define _SHAPE_MAX 1.0
#define _SHAPE_MIN 0.
/*
#define BINN1 16
#define BINN2 2
*/
#define BINN1 8
#define BINN2 2
#define MAX_GSC 11.
#define DMAG 3.

extern double MAX_AST_ERR;
extern double MAX_AST_ERR_FATAL;
extern double SHARP_MIN;
extern double SHARP_MAX;
extern double SHAPE_MAX;
extern double SHAPE_MIN;
/* struct GPARAM param; */


/*Function definitions */
int get_cat( double ra_min,double ra_max,double dec_min,double dec_max, 
				double** ra,double** dec,double** mag,int* ngsc,char* path,
				int field );

/* Structure defines */
struct Coo{
  double ra,dec,x,y,mag,sky;
  float *mags;
  float sharp,shape;
  int link;
};
struct Data{
  double ra,dec;
  double xc,yc;
  double xoff, yoff;
  double pixscale;
  double rotation;
  int nmag;
  int n;
  int type;
  struct Coo *data;
  double xmin,xmax,ymin,ymax,vmin,vmax;
};

extern int verb;
extern int enable_mag;
extern double gOutRA2000,gOutDEC2000;
extern int gMaxTimeForAstrometryInSec;
extern double gInstrumFOV;
// extern int gFixRACloseZero;
extern int gMinNumberOfMatchesToAccept;
extern int gResetGPARAM;

// flag for stoping running astrometry loop :
extern int gForceAstrometryBreak;

#endif
