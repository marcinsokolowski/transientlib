#ifndef _ASAS_INTERFACE_H__
#define _ASAS_INTERFACE_H__

#include <mytypes.h>
#include <mystring.h>


enum eDomeStatus_T { eDomeOpened=0, eDomeClosed, eDomeUnknown };

class CAsasInterface 
{
public :
	CAsasInterface();
	~CAsasInterface();
	
	
	static eDomeStatus_T IsDomeOpened();
	static mystring GetDomeStatusDesc( eDomeStatus_T domeStatus );
	
};


#endif
