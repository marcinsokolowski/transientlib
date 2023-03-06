#ifndef _AST_UTIL_H__
#define _AST_UTIL_H__

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


#include "asas_gparam_def.h"

/******************************************************
* astutil definition file
******************************************************/
#define RD 0.0174532925
#define RADS    57.29577951             /* Degrees per radian */

int needswab();
void swab2();
void swab4();
void swab8();
void eq2ecl();
void eq2eq1();
double gmst();
double  xmass();
void dsortindx();
void fsortindx();
void isortindx();
double equinox();
char* d2hms();
double hms2d();
char *designate();
void undesignate();
/* void Precession(); */
double  Distance();
void Transform();
void Offset();
void geo_trans();
/*void ad2xy();
void xy2ad();*/
/*
*
*/
/*struct GPARAM{
  int order;            // number of parameters
  double px[21];
  double py[21];
  double pixscale;      // first order pixscale
  double fi;            // first order rotation
  double ra,dec,xc,yc;
};*/

double juldate(int y,int m,int d,double ut);

double heljuldate(int y,int m,int d,double ut,double ra,double dec);

double jd2hjd(double jd,double ra,double dec);

void Precession( double alpha0, double delta0, double eqnx0,
					  double u_alpha, double u_delta,
					  double* alpha, double* delta, double eqnx);

double refraction(double z,double pres,double temp);


#endif
