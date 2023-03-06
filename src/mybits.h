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
	val = (val & GET_ZERO_N(n));
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
