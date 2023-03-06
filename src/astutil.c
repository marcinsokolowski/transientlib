/************************************************************
* astro utilities for ASAS software
*
*
* 
* xy2ad()
* ad2xy()
* geo_trans();
* Transform()
* Offset()
* needswab()
* swab2()
* swab4()
* swab8()
************************************************************/

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
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "astutil.h"

#define maxim(A,B) ((A)>(B)?(A):(B))
#define minim(A,B) ((A)<(B)?(A):(B))

void Transform();
void Offset();
void geo_trans();


double jd0jan73=2441682.5, jd1jan70=2440587.5;
double d2r=3.14159265358979323846/180.0;
int doytab[]={0,0,31,59,90,120,151,181,212,243,273,304,334};
 
 
/*******************************************
                  !!!!!!!!!!!!!!!!
* xy2ad ( returns ra - in hours, dec - in deg ... !!!!!!!!!!!!
********************************************/
void xy2ad(param,x,y,ra,dec)
struct GPARAM *param;
double *ra,*dec,x,y;
{
  double da,dd;

  geo_trans(x-param->xc, y-param->yc, param, &da, &dd);
  Transform(param->ra,param->dec,da/3600.,dd/3600.,ra,dec);
}

/*******************************************
* ad2xy
********************************************/
void ad2xy( struct GPARAM * param,double ra,double dec,double* x,double* y )
/*struct GPARAM *param;
double ra,dec,*x,*y;*/
{
  double da,dd,x0,y0,ra1,dec1;
  int i;

  Offset(param->ra,param->dec,ra,dec,&da,&dd);

  x0 = param->xc - da * 3600 / param->px[1];
  y0 = dd * 3600 / param->px[1] + param->yc;
  for(i=0; i<5; i++){
    geo_trans(x0-param->xc,y0-param->yc, param, &da, &dd);
    Transform(param->ra,param->dec,da/3600.,dd/3600.,&ra1,&dec1);
    Offset(ra,dec,ra1,dec1,&da,&dd);
    x0 += da * 3600 / param->px[1];
    y0 += -dd * 3600 / param->px[1];
  }


  *x = x0;
  *y = y0;
}
/***********************************************************
* geometric transformation
***********************************************************/
void geo_trans(x,y,param,xt,yt)
double x,y,*xt,*yt;
struct GPARAM *param;
{
  double u,v;
  double powx[6],powy[6];
  int ix[21]={0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,4,3,2,1,0};
  int iy[21]={0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,1,2,3,4,5};
  int np[6]={1,3,6,10,15,21};
  int i, mp;

  mp = minim(param->order,5);
  powx[0]=powy[0]=1.;
  for(i=1; i<mp+1; i++){
    powx[i] = powx[i-1]*x;
    powy[i] = powy[i-1]*y;
  }
  for(i=0, u=0; i<np[mp]; i++)u += param->px[i]*powx[ix[i]]*powy[iy[i]];
  for(i=0, v=0; i<np[mp]; i++)v += param->py[i]*powx[ix[i]]*powy[iy[i]];
  *xt = u;
  *yt = v;
}

/***********************************************************
* geometric transformation
***********************************************************/
void geo_trans_check(double x,double y,struct GPARAM* param,double* xt,double* yt,
							double* linx,double* nonlinx,double* liny,double* nonliny)
{
  double u,v;
  double powx[6],powy[6];
  int ix[21]={0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,4,3,2,1,0};
  int iy[21]={0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,1,2,3,4,5};
  int np[6]={1,3,6,10,15,21};
  int i, mp;
  double norm;

  (*linx) = 0;
  (*liny) = 0;
  (*nonlinx) = 0;
  (*nonliny) = 0;

  mp = minim(param->order,5);
  powx[0]=powy[0]=1.;
  for(i=1; i<mp+1; i++){
    powx[i] = powx[i-1]*x;
    powy[i] = powy[i-1]*y;
  }
  for(i=0, u=0; i<np[mp]; i++){
  		double val = param->px[i]*powx[ix[i]]*powy[iy[i]];
  		u += val;
  		if( i<3 ){
  			(*linx) = (*linx) + val;
  		}else{
  			(*nonlinx) = (*nonlinx) + val;
  		}
  }		
  for(i=0, v=0; i<np[mp]; i++){
  		double val = param->py[i]*powx[ix[i]]*powy[iy[i]];
  		v += val;
  		if( i<3 ){
  			(*liny) = (*liny) + val;
  		}else{
  			(*nonliny) = (*nonliny) + val;
  		}
  }
  
  norm = ( fabs(*linx) + fabs(*nonlinx) );
  (*linx) = fabs(*linx) / norm;
  (*liny) = fabs(*liny) / norm;
  (*nonlinx) = fabs(*nonlinx) / norm;
  (*nonliny) = fabs(*nonliny) / norm;
  
  *xt = u;
  *yt = v;
}


/*************************************************************
* Transform
* Function calculates coordinates of (RA,DEC) knowing the 
* position (RA0,DEC0) and angular distance (dRA,dDEC) 
* trigonometric formulas are used 
**************************************************************/
void Transform(al0,dl0,da,dd,al,dl)
double al0,dl0,da,dd,*al,*dl;
{
  double d,sin_g,cos_g,sin_a,cos_a,sin_dl,a;

  da = -da*RD;
  dd = dd*RD;

  d = sqrt(da*da + dd*dd);
  if(d ==0){
    *al = al0;
    *dl = dl0;
    return;
  }
  sin_g = da/d;
  cos_g = dd/d;

  sin_dl=sin(dl0*RD)*cos(d) + cos(dl0*RD)*sin(d)*cos_g;
  *dl=asin(sin_dl);
  sin_a=sin(d)*sin_g/cos(*dl);
  cos_a=(cos(d) - sin(dl0*RD)*sin_dl)/cos(dl0*RD)*sqrt(1-sin_dl*sin_dl);
  a=asin(sin_a)/RD;
  if(cos_a < 0) a = 180.-a;
  *al=al0 + a/15.;
  if(*al < 0) *al += 24.;
  if(*al > 24.) *al -= 24.;
  *dl /= RD;
}

/**************************************************************
*  calc. x,y offset (arc degrees) for two alpha,delta points **
* Calculates position difference of 2 points (RA,DEC) and (RA0,DEC0)
**************************************************************/
void  Offset(al0,dl0,al,dl,da,dd)
double  al0,dl0,al,dl,*da,*dd;
{
  double cos_d,sin_d,cos_g,sin_g;

/* correction 200608 */  
  if( al0 > 12.00 && al<6.00){
    // situation on border of 24<->0 :
    // example al0=23.65 and al=0.1
    al = al + 24.00;
/*      al0 = al0 - 24.00;*/
  }
  if( al0 < 6.00 && al>12.00 ){
    // situation on border of 24<->0 :
    // example al0=0.1 and al=23.43
    al0 = al0 + 24.00;
/*      al = al - 24.00;*/
  }

  al0=al0*RD*15.;
  dl0=dl0*RD;
  al=al*RD*15.;
  dl=dl*RD;
  cos_d=sin(dl0)*sin(dl)+cos(dl0)*cos(dl)*cos(al-al0);

	
  if(cos_d >= 0.99999999){ // !!!!!!! in fact >=1.00 - but in O3 and double sometimes does not work fine !!!
    *da = 0;
    *dd = 0;
  }else{	
    sin_d=sqrt(1.-cos_d*cos_d);
    if( fabs(sin_d) > 0.0000 ){
    	sin_g=sin(al-al0)*cos(dl)/sin_d;
    	cos_g=(sin(dl)*cos(dl0)-cos(dl)*sin(dl0)*cos(al-al0))/sin_d;
    	
    	*da = acos(cos_d)*sin_g/RD;
    	*dd = acos(cos_d)*cos_g/RD;
    }else{
    	*da = 0;
    	*dd = 0;
    }
  }
}

/************************************************************
* needswab
************************************************************/
int needswab()   /* 1 on Intel, 0 on Suns  - MSz*/
{
  union { short s; char c[2]; } u;
  u.s=1;
  return((int)(u.c[0]));
}

void  swab2(cp,n)
char *cp;
int n;
{
  char c;
  int i;
  for(i=0; i<n; i++){
    c = *cp;
    *cp = *(cp+1);
    *(++cp) = c;
    cp++;
  }
}
void  swab4(cp,n)
char *cp;
int n;
{
  char c;
  int i;
  for(i=0; i<n; i++){
    c = *cp;
    *cp = *(cp+3);
    *(cp+3) = c;
    c=*(++cp);
    *cp = *(cp+1);
    *(++cp) = c;
    cp+=2;
  }
}
void swab8(cp,n)
char *cp;
int n;
{
  char c,*cp1,*cp2;
  int i;
  cp1=cp;
  for(i=0; i<n; i++) {
    cp2=cp+7;
    while(cp2>cp1){
      c = *cp1;
      *cp1++ = *cp2;
      *cp2-- = c;
    }
    cp1=(cp += 8);
  }
}

void eq2ecl(a,d,l,b)
double a,d,*l,*b;		/* all in degrees */
{
  static double cose=0.917475, sine=0.397793;

  a *= d2r;
  d *= d2r;
  *b=asin(sin(d)*cose - cos(d)*sin(a)*sine)/d2r;
  *l=atan2(cos(d)*sin(a)*cose+sin(d)*sine,cos(d)*cos(a))/d2r;
  if(*l < 0.0) *l += 360.0;
}
 
double gmst(y,m,d,ut)
int y,m,d;
double ut;
{
  double jd,t, juldate();
 
  jd=juldate(y,m,d,0.0);
  t=(jd-2451545.0)/36525.0;
  t=6.6973745583333333+t*(8640184.812866+t*(0.093104-6.2e-6*t))/3600.0;
  t+=(1.00273791*ut);
  t=fmod(t,24.0);
  if ( t < 0.0 ) t+=24.0;
  return(t);
}
void jul2date(jd,y,m,d,ut)
int *y,*m,*d;
double jd, *ut;
{
  double v; 
  int i,*dt;
  static int dt1[12]={31,28,31,30,31,30,31,31,30,31,30,31};
  static int dt2[12]={31,29,31,30,31,30,31,31,30,31,30,31};

  *y = (jd - jd0jan73)/365. + 1973;
  v = jd - (jd0jan73+(*y-1973)*365.0+(*y-1973)/4);
  if (v < 1){
    (*y)--;
    v = jd - (jd0jan73+(*y-1973)*365.0+(*y-1973)/4);
  }
  dt = dt1;
  if ( *y % 4 == 0) dt = dt2;
  *m = 1;
  for (i=0; i<12; i++)  {
    if ((int)v <= dt[i]) break;
    v -= dt[i];
    (*m)++;
  }
  *d = v;
  *ut = 24.*(v - *d);
}

double juldate(int y,int m,int d,double ut)
{
  double jd;

  if(m > 12 || m < 1) return((double)-1.);
  if(y < 70) y += 2000;
    else  y += 1900;
 
  jd=jd0jan73+(y-1973)*365.0+(y-1973)/4+doytab[m]+d+ut/24.0;
  if ( ( y % 4 == 0 ) && ( m > 2 ) ) jd++;
  return(jd);
}

double heljuldate(int y,int m,int d,double ut,double ra,double dec)
 /* ra in hours, dec in degrees */
{
  double jd,hjd,jd2000,anom,longsol,lamsol,rsol,lamb,beta;

  if((jd=juldate(y,m,d,ut)) < 0)return(jd);

  jd2000=jd-2451545.0;
  longsol=280.46+0.9856474*jd2000;
  anom=357.528+0.9856003*jd2000;
  longsol=fmod(longsol,360.0);
  anom=fmod(anom,360.0);
  if ( longsol < 0.0 ) longsol+=360.0;
  if ( anom < 0.0 ) anom+=360.0;
  anom*=d2r;
  lamsol=longsol+1.915*sin(anom)+0.02*sin(anom+anom);
  rsol=1.00014-0.01671*cos(anom)-0.00014*cos(anom+anom);
  eq2ecl(15.*ra,dec,&lamb,&beta);
  hjd=jd-0.0057755*rsol*cos(beta*d2r)*cos(d2r*(lamsol-lamb));

  return(hjd);
}

double jd2hjd(double jd,double ra,double dec)
 /* ra in hours, dec in degrees */
{
  double hjd,jd2000,anom,longsol,lamsol,rsol,lamb,beta;

  jd2000=jd-2451545.0;
  longsol=280.46+0.9856474*jd2000;
  anom=357.528+0.9856003*jd2000;
  longsol=fmod(longsol,360.0);
  anom=fmod(anom,360.0);
  if ( longsol < 0.0 ) longsol+=360.0;
  if ( anom < 0.0 ) anom+=360.0;
  anom*=d2r;
  lamsol=longsol+1.915*sin(anom)+0.02*sin(anom+anom);
  rsol=1.00014-0.01671*cos(anom)-0.00014*cos(anom+anom);
  eq2ecl(15.*ra,dec,&lamb,&beta);
  hjd=jd-0.0057755*rsol*cos(beta*d2r)*cos(d2r*(lamsol-lamb));

  return(hjd);
}

#define SCALE   750.0                   /* Atmospheric scale height */
#define STG0    24055.3   /* Greenwich ST (in sec) on 1970.01.01.00:00:00 GMT */

void  eq2eq1(ra, dec, ut, ha, longitude)
/* ra in (double) hours; dec in (double) degrees */
/* longitude, latitude in degrees */
double  ra, dec;
time_t  ut;                             /* time from 1970.01.01.00:00:00 */
double  *ha;
double longitude;
{
  double  st, stG;
 
  stG=(double)((int)(STG0+366.2422/365.2422*ut)%(int)86400)*24./86400.;
  st=stG+longitude/15.;
  if(st<0)st=st+24.;
  if(st>24)st=st-24.;
  *ha=st-ra;
}

double  xmass(ra, dec, ut, st_out, z_out, longitude, latitude)
/* ra in (double) hours; dec in (double) degrees */
/* longitude, latitude in degrees */
double  ra, dec;
time_t  ut;                             /* time from 1970.01.01.00:00:00 */
double  *st_out,*z_out;
double latitude, longitude;
{
  double  ha,cos_zd, x, st, stG;
 
  stG=(double)((int)(STG0+366.2422/365.2422*ut)%(int)86400)*24./86400.;
  st=stG+longitude/15.;
  if(st<0)st=st+24.;
  if(st>24)st=st-24.;
  ha=st-ra;
  *st_out=st;
 
  cos_zd = sin (latitude*d2r) * sin (dec*d2r) +
           cos (latitude*d2r) * cos (dec*d2r) * cos (15*ha*d2r);
  if(cos_zd>1.0000001)return (double)9.99;
  else if(cos_zd<=1.) *z_out=acos(cos_zd)/d2r;
  else  *z_out=0.;
  if(cos_zd < 0) return (double)99.9;
 
  x  = SCALE * cos_zd;
  x=(sqrt (x*x + 2*SCALE + 1) - x);
  if(x>99.9)x=(double)99.9;
 
  return(x);
}
void dsortindx(n,arrin,indx)
int n, *indx;
double *arrin;
{
  int l,ir,i,j,indxt;
  double q;
 
  if(n<=0)return;
  for ( i=0; i<n; i++ ) indx[i]=i;
  if ( n <= 1 ) return;
  l=n/2;
  ir=n-1;
  while ( 10 ) {
    if ( l > 0 ) {
      l--;
      indxt=indx[l];
      q=arrin[indxt];
    } else {
      indxt=indx[ir];
      q=arrin[indxt];
      indx[ir]=indx[0];
      ir--;
      if ( ir == 0 ) {
        indx[0]=indxt;
        return;
      }
    }
    i=l;
    j=l+l+1;
    while ( j <= ir ) {
      if ( j < ir )
        if ( arrin[indx[j]] < arrin[indx[j+1]] ) j++;
      if ( q < arrin[indx[j]] ) {
        indx[i]=indx[j];
        i=j;
        j+=j+1;
      } else {
        j=ir+1;
      }
    }
    indx[i]=indxt;
  }
}

void fsortindx(n,arrin,indx)
int n, *indx;
float *arrin;
{
  int l,ir,i,j,indxt;
  float q;
 
  if(n<=0)return;
  for ( i=0; i<n; i++ ) indx[i]=i;
  if ( n == 1 ) return;
  l=n/2;
  ir=n-1;
  while ( 10 ) {
    if ( l > 0 ) {
      l--;
      indxt=indx[l];
      q=arrin[indxt];
    } else {
      indxt=indx[ir];
      q=arrin[indxt];
      indx[ir]=indx[0];
      ir--;
      if ( ir == 0 ) {
        indx[0]=indxt;
        return;
      }
    }
    i=l;
    j=l+l+1;
    while ( j <= ir ) {
      if ( j < ir )
        if ( arrin[indx[j]] < arrin[indx[j+1]] ) j++;
      if ( q < arrin[indx[j]] ) {
        indx[i]=indx[j];
        i=j;
        j+=j+1;
      } else {
        j=ir+1;
      }
    }
    indx[i]=indxt;
  }
}
void isortindx(n,arrin,indx)
int n, *indx;
int *arrin;
{
  int l,ir,i,j,indxt;
  int q;
 
  if(n<=0)return;
  for ( i=0; i<n; i++ ) indx[i]=i;
  if ( n == 1 ) return;
  l=n/2;
  ir=n-1;
  while ( 10 ) {
    if ( l > 0 ) {
      l--;
      indxt=indx[l];
      q=arrin[indxt];
    } else {
      indxt=indx[ir];
      q=arrin[indxt];
      indx[ir]=indx[0];
      ir--;
      if ( ir == 0 ) {
        indx[0]=indxt;
        return;
      }
    }
    i=l;
    j=l+l+1;
    while ( j <= ir ) {
      if ( j < ir )
        if ( arrin[indx[j]] < arrin[indx[j+1]] ) j++;
      if ( q < arrin[indx[j]] ) {
        indx[i]=indx[j];
        i=j;
        j+=j+1;
      } else {
        j=ir+1;
      }
    }
    indx[i]=indxt;
  }
}
void ssortindx(n,arrin,indx)
int n, *indx;
short *arrin;
{
  int l,ir,i,j,indxt;
  int q;
 
  if(n<=0)return;
  for ( i=0; i<n; i++ ) indx[i]=i;
  if ( n == 1 ) return;
  l=n/2;
  ir=n-1;
  while ( 10 ) {
    if ( l > 0 ) {
      l--;
      indxt=indx[l];
      q=arrin[indxt];
    } else {
      indxt=indx[ir];
      q=arrin[indxt];
      indx[ir]=indx[0];
      ir--;
      if ( ir == 0 ) {
        indx[0]=indxt;
        return;
      }
    }
    i=l;
    j=l+l+1;
    while ( j <= ir ) {
      if ( j < ir )
        if ( arrin[indx[j]] < arrin[indx[j+1]] ) j++;
      if ( q < arrin[indx[j]] ) {
        indx[i]=indx[j];
        i=j;
        j+=j+1;
      } else {
        j=ir+1;
      }
    }
    indx[i]=indxt;
  }
}
double equinox(ut)
time_t ut;
{
  int    yr;
  double doy,hr,mi,sc,dpy,ep;
  char   buffer[128];
  struct tm *gtime;
 
  gtime = gmtime(&ut);                            /* get tm-structure (UT)*/
  strftime(buffer,80,"%Y",gtime);                 /* calculate epoch      */
  yr  = atoi(buffer);                             /* year                 */
  if ((yr % 4) == 0) {
    if ((yr % 100) == 0) {
      if ((yr % 400) == 0)
        dpy = 366.;
      else
        dpy = 355.;
    } else {
      dpy = 366.;
    }
  } else {
    dpy = 355.;
  }
  strftime(buffer,80,"%j",gtime);                 /* day                  */
  doy = atof(buffer) - 1.;
  strftime(buffer,80,"%H",gtime);                 /* hour                 */
  hr = atof(buffer);
  strftime(buffer,80,"%M",gtime);                 /* minute               */
  mi = atof(buffer);
  strftime(buffer,80,"%S",gtime);                 /* second               */
  sc = atof(buffer);
  ep = (double)yr + doy/dpy + hr/(24.*dpy) + mi/(1440.*dpy) + sc/(86400.*dpy);
 
  return(ep);
}
void undesignate(desig, ra, dec)
char *desig;
double  *ra, *dec;
{
  int hh,mm,ss,ad;
  double ff;

  sscanf(desig,"%02d%02d%02d%03d%4lf",&hh,&mm,&ss,&ad,&ff);

  *ra = (double)hh +(double)mm/60. + (double)ss/3600;
  *dec = (double)abs(ad) + ff/60.;
  if(desig[6] == '-')*dec *= -1;
}
char *designate(ra,dec)
double ra,dec;
{
  static char desig[32]; 
  int hh,mm,ss,ad;
  char sg;
  double ff;
  ff=ra;
  hh=ff;
  ff = (ff-hh)*60;
  mm=ff;
  ff = (ff-mm)*60;
  ss=ff+0.5;
  if(ss == 60)ss=0,mm++;
  if(mm == 60)mm=0,hh++;
  if(hh == 24)hh=0;
  ff=fabs(dec);
  ad=ff;
  ff = (ff-ad)*60;
  ad = abs(ad);
  if(dec < 0)sg = '-';
    else sg ='+';
  sprintf(desig,"%02d%02d%02d%c%02d%04.1f",hh,mm,ss,sg,ad,ff);
  return(desig);
}

char *d2hms(x,mode)
double x;
int mode;
{
  static char buf[32];
  double v,ss,mm;
  int vh,vm,vs;
/* modes:
   4 -  2:02
   2 -  02:02.22
   1 -  02:02:02
  oth-  02:02:02.2
*/
  v=fabs(x);
  vh = v;
  mm=(v-vh)*60.;
  vm = mm;
  ss = (v-vh - vm/60.)*3600;
  vs = (int)(ss+0.5);
  if(x < 0) vh = -vh;
  if(mode == 4){
    if(vs > 60)vm++;
    if(vm == 60)vm=0., vh++;
    sprintf(buf,"%2d:%02d",vh,vm);
  }else if(mode == 2){
    sprintf(buf,"%02d:%05.2f",vh,mm);
  }else if(mode == 1){
    if(vs == 60)vs=0., vm++;
    if(vm == 60)vm=0., vh++;
    sprintf(buf,"%02d:%02d:%02d",vh,vm,vs);
  }else{
    if(vm == 60)vm=0., vh++;
    sprintf(buf,"%02d:%02d:%04.1f",vh,vm,ss);
  }
  return(buf);
}
double hms2d(cval)
char *cval;
{
  double h,m,s;
  static double dval;
 
  h=m=s=0.;
  if(strchr(cval,':') != NULL){
    sscanf(cval,"%lf:%lf:%lf",&h,&m,&s);
  }else{
    sscanf(cval,"%lf",&h);
    m=s=0.;
  }
  dval=(fabs(h)+fabs(m)/60.+fabs(s)/3600.);
  if(h<0) dval = -dval;
  return(dval);
}
#define RADARCS (double)(0.0174532925/3600.)    /* radians per arcsecond */
#define RADDEG (double)0.0174532925     /* radians per arcdeg */
#define RADHOUR (double)0.2617993878    /* radians per hour */
/**************************************************************
*  calc. distance in arcseconds between  two alpha,delta points **
**************************************************************/
double  Distance(al0,dl0,al,dl)
double  al0,dl0,al,dl;
{
  double cos_d;

  al0=al0*RADHOUR;
  dl0=dl0*RADDEG;
  al=al*RADHOUR;
  dl=dl*RADDEG;

  cos_d=sin(dl0)*sin(dl)+cos(dl0)*cos(dl)*cos(al-al0);
  if(cos_d >= 1.){
    return((double)0.);
  }else{
    return(acos(cos_d)/RADARCS);
  }
}
/******************************
*** precession calculation  ***
******************************/
void Precession(alpha0,delta0,eqnx0,u_alpha,u_delta,alpha,delta,eqnx)
double  alpha0,delta0,eqnx0;            /* coo in eqnx0 */
double  u_alpha,u_delta;                /* proper motions */
double  *alpha,*delta,eqnx;             /* coo in eqnx */
{
  double  c1=2.062648e5;
  double  pi=3.141592653;

  double al,del;
  double  tau0,tau;
  double  zeta,theta,z,ac,bc,cc,zx,a1,d1;

  if(eqnx0==eqnx){
    *alpha=alpha0;
    *delta=delta0;
    return;
  }
  al=(alpha0+u_alpha/3600.*(eqnx-eqnx0))*pi/12.;  /*  coord in radians */
  del=(delta0+u_delta/3600.*(eqnx-eqnx0))*pi/180.;

/* NEWCOMB: */
  tau0=(eqnx0-1900.)/100.;
  tau=(eqnx-1900.)/100.-tau0;

  zeta=(2304.25+1.396*tau0)*tau+.302*tau*tau+0.018*tau*tau*tau;
  z= zeta+0.791*tau*tau+0.001*tau*tau*tau;
  theta=(2004.682-0.853*tau0)*tau-0.426*tau*tau-0.042*tau*tau*tau;
  zeta=zeta/c1;
  z=z/c1;
  theta=theta/c1;

  ac=cos(del)*sin(al+zeta);
  bc=cos(theta)*cos(del)*cos(al+zeta)-sin(theta)*sin(del);
  cc=sin(theta)*cos(del)*cos(al+zeta)+cos(theta)*sin(del);
  d1=asin(cc);
  zx=sqrt(1.-cc*cc);
  a1=bc/zx;
  if(a1 > 1.) a1 =  1.;
  if(a1 <-1.) a1 = -1.;
  a1=acos(a1);
  if((ac/zx)<0.) a1=2.*pi-a1;
  a1=a1+z;

  *alpha=a1*12./pi;
  if(*alpha<0)*alpha+=24;
  if(*alpha>=24)*alpha-=24;
  *delta=d1*180./pi;
}


double refraction(z,pres,temp)
double z;
double pres,temp;
{
  double R,a;
 
  a=90.-z;
 
  if(z>15.){
    R=0.00452*pres/(273.+temp)/tan(a/RADS);
  }else{
    R=pres*(0.1594+0.0196*a+0.00002*a*a)/(273.+temp)/(1.+0.505*a+0.0845*a*a);
  }
  z=z+R;
  return(z);
}

double calc_diff_ra( double ra, double ra0 )
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
}
