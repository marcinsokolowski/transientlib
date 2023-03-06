/***************************************************************************\
 **  transientlib : library for identification of short optical flashes 
 **						with the wide field cameras.
 **  This software was written by Marcin Sokolowski ( msok@fuw.edu.pl ) 
 **	it was a substantial part of analysis performed for PHD thesis 
 **  ( Investigation of astrophysical phenomena in short time scales with "Pi of the Sky" apparatus )
 **	it can be used under the terms of GNU General Public License.
 **	In case this software is used for scientific purposes and results are
 **	published, please refer to the PHD thesis submited to astro-ph :
 **
 **		http://arxiv.org/abs/0810.1179
 **
 ** Public distribution was started on 2008-10-31
 **
 ** 
 ** NOTE : some of the files (C files) were created by other developers and 
 **        they maybe distributed under different conditions.
 ** 

 ******************************************************************************
 ** This program is free software; you can redistribute it and/or modify it
 ** under the terms of the GNU General Public License as published by the
 ** Free Software Foundation; either version 2 of the License or any later
 ** version. 
 **
 ** This program is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY; without even the implied warranty of
 ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 ** General Public License for more details. 
 **
 *\**************************************************************************

*/           
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