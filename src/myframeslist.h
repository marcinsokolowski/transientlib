#ifndef _MY_FRAMES_LIST_H__
#define _MY_FRAMES_LIST_H__

#include <mystring.h>
#include <vector>

using namespace std;

struct CFrameInfo
{	
	CFrameInfo() : m_FrameUT(0) {}

	mystring m_szFileName;
	time_t   m_FrameUT;
};

class CFramesList : public vector<CFrameInfo>
{
public :
	CFramesList();
	~CFramesList();	

	int ReadList( const char* listname );
	int UpdateList( const char* listname, BOOL_T bCheck=TRUE );
	
	void AddFrame( const char* szFileName, time_t ut_time );
	CFrameInfo* FindClosest( time_t ut_time );
	CFrameInfo* Find( const char* szFName );
};


#endif