#include "myvalcounter.h"
#include <string.h>


CValuesCounter::CValuesCounter()
:m_Limit(INITIAL_SIZE),cnt(0)
{	
	m_pValues = new CValueCounter[m_Limit]();
}

CValuesCounter::~CValuesCounter()
{
	delete [] m_pValues;
}

void CValuesCounter::AddNew(LONG_T value)
{
		if(cnt>=m_Limit){
			CValueCounter* pOld = m_pValues;
			LONG_T prev_size=m_Limit;
			m_Limit += INITIAL_SIZE;
			m_pValues = new CValueCounter[m_Limit]();
			memcpy(m_pValues,pOld,sizeof(CValueCounter)*prev_size);
			delete [] pOld;
			
		}	
		m_pValues[cnt].value = value;
      m_pValues[cnt].counter++;
		cnt++;
}

void CValuesCounter::Add(LONG_T value)
{
	CValueCounter* pVal = FindValue(value);
	if(pVal){
		pVal->counter++;
	}else{
		AddNew( value );
	}
}

CValueCounter* CValuesCounter::FindValue( LONG_T value )
{
	for(int i=0;i<cnt;i++){
		if(m_pValues[i].value == value )
			return (&(m_pValues[i]));
	}
	return NULL;
}
