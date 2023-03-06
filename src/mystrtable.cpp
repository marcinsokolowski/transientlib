#include "mystrtable.h"
#include "mydate.h"
/*#if __GNUC__ >= 4 && __GNUC_MINOR__ >=6 && __GNUC_PATCHLEVEL__>=1
	#include <algorithm>
	using namespace std;
#else
	#include <algo.h>
#endif*/
#include <algorithm>
using namespace std;
   
#include "myfile.h"

CMyStrTable::CMyStrTable(const CMyStrTable& right)
{
	(*this) = right;
}

void CMyStrTable::Add(const char* str)
{
	mystring tmp(str);
	push_back(tmp);
}


CMyStrTable& CMyStrTable::operator=(const CMyStrTable& right)
{
	clear();
	for(register int i=0;i<right.size();i++){
		push_back( right[i] );
	}
	return (*this);
}


int CMyStrTable::Remove(const char* str)
{
	CMyStrTable::iterator i;

	int pos=0;
	for(i=begin();i!=end();i++){
		if(strcmp(i->c_str(),str)==0){
			erase( i );
			return pos;
		}
		pos++;
	}
	return -1;
}

int CMyStrTable::FindPos(const char* str)
{
	CMyStrTable::iterator i;

	int pos=0;
	for(i=begin();i!=end();i++){
		if(strcmp(i->c_str(),str)==0)
			return pos;
		pos++;	
	}
	return -1;
}

const char* CMyStrTable::Find(const char* str)
{
	CMyStrTable::iterator i;

   for(i=begin();i!=end();i++){
		if(strcmp(i->c_str(),str)==0)
			return i->c_str();
	}
	return NULL;	
}

mystring&  CMyStrTable::get(LONG_T pos)
{
	return (*this)[pos];
}


void CMyStrTable::Sort()
{
	sort( begin(), end() );	
}


int CompareByDate( const char* left, const char* right, const char* fmt )
{
	double leftVal = CMyDate::getTime( left, fmt );
	double rightVal = CMyDate::getTime( right, fmt );


	if(leftVal>rightVal)
		return 1;
	if(leftVal<rightVal)
		return -1;

	return 0;
}

// insertion sort :
void CMyStrTable::SortByDateFormat( const char* fmt )
{
	CMyStrTable restab;


	for(int i=0;i<(size()-1);i++){
		mystring firstVal = (*this)[i];
		mystring szMinVal = firstVal.c_str();
		int min_pos = i;

		for(int j=i+1;j<size();j++){
			if( CompareByDate( (*this)[j].c_str(), szMinVal.c_str(), fmt ) < 0 ){
				// means that minVal > currVal :
				szMinVal = (*this)[j];
				min_pos = j;
			}			
		}
		if(min_pos!=i){
			(*this)[i] = (*this)[min_pos];
			(*this)[min_pos] = firstVal;
		}
	}
}

void CMyStrTable::DumpToFile( const char* szFName )
{
	MyOFile out( szFName , "a" );
	for(int i=0;i<size();i++){
		out.Printf("%s\n",(*this)[i].c_str());		
	}
}
