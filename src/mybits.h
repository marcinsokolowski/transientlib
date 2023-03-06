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
#ifndef _MYBITS_H_INCLUDED_
#define _MYBITS_H_INCLUDED_

#include "mytypes.h"

// bits start from bit 0 = 2^0
#define BIT_0 0x01
#define BIT_1 0x02
#define BIT_2 0x04
#define BIT_3 0x08
#define BIT_4 0x10
#define BIT_5 0x20
#define BIT_6 0x40
#define BIT_7 0x80


#define GET_BIT_N(n) (0x01 << n)
#define GET_ZERO_N(n) (~(0x01 << n))

inline int CheckBit(const char& val,char n)
{
	return ((val & GET_BIT_N(n))!=0);	
}

inline void HitBit(char& val,char n)
{
	val = (val | GET_BIT_N(n));	
}

inline void ClearBit(char& val,char n)
{
	(val & GET_ZERO_N(n));
}


char DecodeU2( char u2_value );

unsigned short GetValueFromBytes( char lsb, char msb );

int GetValueOfRightBits_int( int val, char nBits );
int GetValueOfLeftBits_int( int val, char nBits );

char GetValueOfRightBits_char( char val, char nBits );
char GetValueOfLeftBits_char( char val, char nBits );

long GetValueOfRightBits_long( long val, char nBits );
long GetValueOfLeftBits_long( long val, char nBits );

double GetValueOfRightBits_double( double val, char nBits );
double GetValueOfLeftBits_double( double val, char nBits );


#endif
