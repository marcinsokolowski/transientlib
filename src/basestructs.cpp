#include "basestructs.h"

CEnvVar::CEnvVar( const char* sName, const char* sValue, const char* sComment )
: szName(sName),szValue(sValue)
{
	if(sComment && sComment[0])
		szComment = sComment;
}

CEnvVar::CEnvVar(const CEnvVar& right){
	(*this) = right; 
}				

CEnvVar& CEnvVar::operator=(const CEnvVar& right){
	szName = right.szName; 
	szValue = right.szValue;	
	szComment = right.szComment;
	return (*this);
}				



