#include "paramtab.h"


BOOL_T CParamTab::IsParamDefined( const char* name )
{
	CParamTab::iterator i;
	for(i=begin();i!=end();i++){
		if(strcmp(i->szName.c_str(),name)==0)
			return TRUE;		
	}
	return FALSE;
}


const char* CParamTab::GetParamValue( const char* name )
{
	CParamTab::iterator i;
	for(i=begin();i!=end();i++){
		if(strcmp(i->szName.c_str(),name)==0)
			return i->szValue.c_str();		
	}
	return NULL;
}

