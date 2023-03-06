#include "cexcp.h"
#include "myfract.h"
#include "myparser.h"

CFraction::CFraction( LONG_T _up, LONG_T _down )
{
	Set( _up, _down );
}

CFraction::CFraction( const char* szFraction )
{
	MyParser pars;
	pars.GetFraction( szFraction, (*this) );
}

void CFraction::Set( LONG_T _up, LONG_T _down )
{
	up = _up;
	down = _down;
	Assert(down!=0,"Fraction cannot have 0 on the bottom");
}

double CFraction::Val()
{
	return (double(up)/double(down));
}
