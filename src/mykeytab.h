#ifndef _MYKEYTAB_H__
#define _MYKEYTAB_H__

#include "mystring.h"
#include "basestructs.h"
#include <vector>


using namespace std;



class CKeyTab : public vector<CEnvVar>
{
public:
	CKeyTab(){};
	~CKeyTab(){};
	void Add( const CEnvVar& elem );		
	void Add(const char* key,const char* val, const char* comment="");
	void Add(const char* key,int val,const char* comment="");
	void Add(const char* key,double val, const char* comment="");
	void Delete( const char* key);
	CEnvVar* Find(const char* key);
	CEnvVar* Find(const char* key, int pos);

	void Set( const char* key,const char* val, const char* comment="" );
	void Set( const char* key,double val, const char* comment="" );
	void Set( const char* key,int val, const char* comment="" );
};


#endif
