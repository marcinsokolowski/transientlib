#ifndef _PARAMTAB_H__
#define _PARAMTAB_H__

#include "basestructs.h"
#include <vector>

using namespace std;



class CParamTab : public vector<CEnvVar>
{
public :
	CParamTab(){};
	~CParamTab(){};
	
	BOOL_T IsParamDefined( const char* name );
	const char* GetParamValue( const char* name );
};

#endif
