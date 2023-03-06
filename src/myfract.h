#ifndef _MY_FRACTION_H__
#define _MY_FRACTION_H__

#include "mytypes.h"
#include <math.h>

class CFraction
{
public :
	LONG_T up;
	LONG_T down;
	
	CFraction( LONG_T _up, LONG_T _down );
	CFraction( const char* szFraction );
	
	void Set( LONG_T _up, LONG_T _down );
	double Val();
};


#endif

