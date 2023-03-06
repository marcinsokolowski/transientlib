#ifndef _MY_PTR_TABLE_H__
#define _MY_PTR_TABLE_H__


#include <vector>

using namespace std;

template<class ARG_TYPE>
class CMyPtrTable : public vector<void*>
{
public:	
	CMyPtrTable();
	~CMyPtrTable();
	void Add(ARG_TYPE* ptr);

	ARG_TYPE* Find( const char* szObjectName ); 	
};


#endif
