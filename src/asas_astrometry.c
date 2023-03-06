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
#include <stdio.h>
#include <string.h>
#include "asas_astrometry.h"
#include "asas_fitsio.h"


int verb=0;

double calc_diff_ra( double ra, double ra0 );
/*
{
  double diff_ra=(ra-ra0);
  if( ra > 12.00 && ra0<6.00 ){
    // ra=23.43 ra0=0.23
    diff_ra = (ra-(ra0+24.00));
  }      
  if( ra<6.00 && ra0>12.00 ){
    // ra=0.23 ra0=23.4
    diff_ra = ( ra+24.00 - ra0 );
  }
  return diff_ra;
}*/


double MAX_AST_ERR = _MAX_AST_ERR/2.;		//.i.e. divided by FWHM=2.
double MAX_AST_ERR_FATAL = _MAX_AST_ERR_FATAL/2.; // --""--
double SHARP_MIN = _SHARP_MIN;
double SHARP_MAX = _SHARP_MAX;
double SHAPE_MAX = _SHAPE_MAX;
double SHAPE_MIN = _SHAPE_MIN;

double *sao_x, *sao_y, *star_x, *star_y;

double gOutRA2000,gOutDEC2000;
int gMaxTimeForAstrometryInSec=-1;
double gInstrumFOV=22.00;
// int gFixRACloseZero=0;
int gMinNumberOfMatchesToAccept=150;
int gResetGPARAM=1;

// NEW - 20080414 ( due to problem with old field astrometry overwritting
// new coordinates ) - 2008-04-12/13 - see LM e-mail
// flag to force stopping of running astrometry 
// for example due to field change ( no point to continue astrometry
// for old field ) :
int gForceAstrometryBreak=0;

// int verb=0;
int enable_mag=0;


