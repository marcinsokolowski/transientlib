#ifndef _MY_STRING_TABLE_HH
#define _MY_STRING_TABLE_HH

#include "mytypes.h"
#include "basedefines.h"
#include "mystring.h"
#include <vector>

using namespace std;

BASELIB_EI const char* mygetenv(const char* varname);

class BASELIB_EI  CMyStrTable : public vector<mystring>
{
public:
   CMyStrTable(){};
	CMyStrTable(const CMyStrTable& right);   
   
   CMyStrTable& operator=(const CMyStrTable& right);
   
	void Add(const char* str);
	void AddInt( int value );
	const char* Find(const char* str);
	int FindPos(const char* str);
	int Remove(const char* str);

	mystring& get(LONG_T pos);
	LONG_T GetCount(){ return size(); };
	void Sort();

	void SortByDateFormat( const char* fmt );
	
	void DumpToFile( const char* szFName );
};


#endif
