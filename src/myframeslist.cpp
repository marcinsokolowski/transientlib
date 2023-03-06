#include "myframeslist.h"
#include "myfile.h"

CFramesList::CFramesList()
{
}

CFramesList::~CFramesList()
{
}

int CFramesList::UpdateList( const char* listname, BOOL_T bCheck )
{
	if( MyFile::DoesFileExist( listname )){
		MyIFile in( listname );

		const char* pLine = NULL;
		while(pLine = in.GetLine(TRUE)){
      	if(strlen(pLine) && pLine[0]!='#'){
				if( !bCheck || !Find( pLine ) ){
	         	AddFrame( pLine , 0 );
				}
			}
	   }
	}else{
		printf("File : %s does not exist\n",listname);
	}

}
                                                                                
int CFramesList::ReadList( const char* listname )
{
	clear();
	UpdateList( listname, FALSE );

	return size();
}


void CFramesList::AddFrame( const char* szFileName, time_t ut_time )
{
	CFrameInfo tmp;
	tmp.m_szFileName = szFileName;
	tmp.m_FrameUT = ut_time;

	push_back( tmp );	
}


CFrameInfo* CFramesList::FindClosest( time_t ut_time )
{
	for(int i=1;i<size();i++){
		if( (*this)[i-1].m_FrameUT<ut_time && (*this)[i].m_FrameUT>=ut_time ){
			return &((*this)[i]);
		}
	}
	return NULL;
}



CFrameInfo* CFramesList::Find( const char* szFName )
{
	for(int i=0;i<size();i++){
		if( strcmp( (*this)[i].m_szFileName.c_str(), szFName )==0 ){
			return &((*this)[i]);
		}
	}
	return NULL;
}