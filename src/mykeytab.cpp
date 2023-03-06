#include "mykeytab.h"


void CKeyTab::Add(const char* key,const char* val, const char* comment)
{
	CEnvVar tmp;
	tmp.szName = key;
	tmp.szValue = val;
	tmp.szComment = comment;
	push_back(tmp);
}

void CKeyTab::Add(const char* key,int val, const char* comment)
{
	CEnvVar tmp;
	tmp.szName = key;
	tmp.szValue << val;
	tmp.szComment = comment;
	push_back(tmp);
}

void CKeyTab::Add( const char* key,double val, const char* comment )
{
	CEnvVar tmp;
	tmp.szName = key;
	tmp.szComment = comment;
	char szTmp[64];
	sprintf(szTmp,"%.8f",val);
	tmp.szValue = szTmp;
	push_back(tmp);
}

void CKeyTab::Add( const CEnvVar& elem )
{
	push_back( elem );
}

void CKeyTab::Delete( const char* key )
{
	vector<CEnvVar>::iterator i;
	for(i=begin();i!=end();i++){
		if(strcmp(i->szName.c_str(),key)==0){
			erase( i );	
			return;
		}
	}
}

CEnvVar*  CKeyTab::Find(const char* key)
{
	vector<CEnvVar>::iterator i;
	for(i=begin();i!=end();i++){
		if(strcmp(i->szName.c_str(),key)==0)
			return &(*i);
	}	
	return NULL;
}

CEnvVar*  CKeyTab::Find(const char* key, int pos)
{
	int count=size();
	if( pos>=0 && pos<count ){
		for(int i=pos;i<count;i++){
			if( strcmp( (*this)[i].szName.c_str(),key)==0 ){
				return &( (*this)[i] );
			}			
		}		
	}
	return CKeyTab::Find( key );
}


void CKeyTab::Set( const char* key,double val, const char* comment )
{
	char szTmp[64];
   sprintf(szTmp,"%.8f",val);

	Set( key, szTmp, comment );
}

void CKeyTab::Set( const char* key,int val, const char* comment )
{
	char szTmp[64];
   sprintf(szTmp,"%d",val);

	Set( key, szTmp, comment );
}

void CKeyTab::Set( const char* key,const char* val, const char* comment )
{
	CEnvVar* pKey = Find( key );
	if(pKey){
		pKey->szValue = val;
		if(comment && comment[0]){
			pKey->szComment = comment;
		}
	}else{
		Add( key, val, comment );
	}
}