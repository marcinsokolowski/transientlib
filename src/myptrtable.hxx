#ifndef _MY_PTR_TABLE_HXX__
#define _MY_PTR_TABLE_HXX__

#include "myptrtable.h"

template<class ARG_TYPE>
CMyPtrTable<ARG_TYPE>::CMyPtrTable()
{
}


template<class ARG_TYPE>
CMyPtrTable<ARG_TYPE>::~CMyPtrTable()
{
	vector<void*>::iterator i;
	for(i=begin();i!=end();i++){
		delete (ARG_TYPE*)(*i);
	}
}

template<class ARG_TYPE>
void CMyPtrTable<ARG_TYPE>::Add(ARG_TYPE* ptr)
{
	vector<void*>::push_back( (void*)ptr );
}


template<class ARG_TYPE>
ARG_TYPE* CMyPtrTable<ARG_TYPE>::Find( const char* szObjectName )
{
	vector<void*>::iterator i;
	for(i=begin();i!=end();i++){
		if(strcmp( ((ARG_TYPE*)(*i))->GetName(), szObjectName )==0)
			return (ARG_TYPE*)(*i);
	}
	return NULL;
}

#endif
