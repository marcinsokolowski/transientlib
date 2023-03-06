#ifndef _MY_VAL_COUNTER_H__
#define _MY_VAL_COUNTER_H__

#include "mytypes.h"

class CValueCounter
{
public:
	CValueCounter():value(0),counter(0){};

	LONG_T value;
	LONG_T counter;
};

#define INITIAL_SIZE 500

class CValuesCounter 
{
public :	
	CValuesCounter();
	~CValuesCounter();
	void Add(LONG_T value);
	void AddNew(LONG_T value);
	CValueCounter* FindValue( LONG_T value );

	CValueCounter* m_pValues;
	LONG_T cnt;		
	LONG_T m_Limit;	
};



#endif
