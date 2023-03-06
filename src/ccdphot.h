#ifndef _CCDPHOT_ASAS_H__
#define _CCDPHOT_ASAS_H__


/* 	All Sky Automated Survey software
*
*	Copyright 1998+ - Grzegorz Pojmanski gp@sirius.astrouw.edu.pl
*
*	ccdphot.h 
*/
#define R_APERT_1 2.0
#define R_APERT_2 3.0
#define R_APERT_3 4.0
#define R_APERT_4 5.0
#define R_APERT_5 6.0
#define THRESHOLD 2.0
#define NAPERTS 1

#define VERBOSE		0x0001
#define PRINT_VER	0x0002
#define NOCHECK		0x0004

int phot(char* finp,char* fout,char* fraw,double thre,double* ap,
         int nap,int flag,char* use_ker);

#endif
