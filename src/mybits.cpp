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
#include "mybits.h"

/*char GetBitN(char n){
	if(n==1)
		return BIT_1;
	if(n==2)
		return BIT_2;
	if(n==3)
		return BIT_3;
	if(n==4)
		return BIT_4;
	if(n==5)
		return BIT_5;
	if(n==6)
		return BIT_6;
	if(n==7)
		return BIT_7;
	if(n==8)
		return BIT_8;
	return 0;
}*/


/*BOOL_T CheckBit(const char& val,char n)
{
	char tmp = (val & GET_BIT_N(n));
	return (tmp!=0);	 	
}

void HitBit(char& val,char n)
{
	val = (val | GET_BIT_N(n));
}


void ClearBit(char& val,char n)
{
	if(CheckBit(val,n)){
		val = (val ^ GET_BIT_N(n));
	}
}*/



int GetValueOfRightBits_int( int val, char nBits )
{
	int tmp = (val << nBits);
   int ret = (tmp >> nBits);
   return ret;
}

int GetValueOfLeftBits_int( int val, char nBits )
{
	char nBits2 = (sizeof(int)-nBits);
	int tmp = (val >> nBits2);
   int ret = (tmp << nBits2);
   return ret;
}


char GetValueOfRightBits_char( char val, char nBits )
{
	char tmp = (val << nBits);
	char ret = (tmp >> nBits);
	return ret;	
}

char GetValueOfLeftBits_char( char val, char nBits )
{
	char nBits2 = (8-nBits);
	char tmp = (val >> nBits2);
	char ret = (tmp << nBits2);
	return ret;	
}

long GetValueOfRightBits_long( long val, char nBits )
{
	long tmp = (val << nBits);
	long ret = (tmp >> nBits);
	return ret;	
}


long GetValueOfLeftBits_long( long val, char nBits )
{
	char nBits2 = (sizeof(long)-nBits);
	long tmp = (val >> nBits2);
	long ret = (tmp << nBits2);
	return ret;	
}


double GetValueOfRightBits_double( double val, char nBits )
{
	double tmp = val;
//	tmp << nBits;
//	tmp >> nBits;
	return tmp;	
}

double GetValueOfLeftBits_double( double val, char nBits )
{
	double tmp = val;
//	tmp >> nBits;
//	tmp << nBits;
	return tmp;
}

char DecodeU2( char u2_value )
{
	char ret=0;

	ret = - (CheckBit( u2_value, 7 )*128);

	char ret2 =0;
	for(register int i=6;i>=0;i--){
		ret2 *= 2;
		ret2 += CheckBit( u2_value, i );
	}	
	return (ret+ret2);
}

unsigned short GetValueFromBytes( char lsb, char msb )
{
	unsigned short ret=0;
	((char*)(&ret))[0] = lsb; // LSB -> 0 
	((char*)(&ret))[1] = msb; // MSB -> 1
	return ret;
}

