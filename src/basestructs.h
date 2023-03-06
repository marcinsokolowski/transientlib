#ifndef _BASESTRUCTS_H__
#define _BASESTRUCTS_H__ 

#include "mystring.h"

class  CEnvVar {
public:
	mystring szName;
	mystring szValue;
	mystring szComment;

	CEnvVar(){};		
	CEnvVar( const char* sName, const char* sValue, const char* sComment="" );
	CEnvVar(const CEnvVar& right);
	CEnvVar& operator=(const CEnvVar& right);


	friend bool operator > ( const CEnvVar& left, const CEnvVar& right ){
		return (strcmp(left.szName.c_str(),right.szName.c_str())>0);
	}
	friend bool operator < ( const CEnvVar& left, const CEnvVar& right ){
		return (strcmp(left.szName.c_str(),right.szName.c_str())<0);
	}
	friend bool operator == ( const CEnvVar& left, const CEnvVar& right ){
		return (strcmp(left.szName.c_str(),right.szName.c_str())==0 &&
		        strcmp(left.szValue.c_str(),right.szValue.c_str())==0);
	}
};

struct sCommand
{
   int command;
   mystring szParam1;
	mystring szParam2;
};
         

/*bool CompareEnvVar( const CEnvVar& left, const CEnvVar& right )
{
	return (strcmp(left.szName.c_str(),right.szName.c_str())==1);
}*/


#endif
