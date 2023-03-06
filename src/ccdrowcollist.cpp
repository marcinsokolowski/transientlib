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
#include "ccdrowcollist.h"
#include "myfile.h"
#include "mystrtable.h"
#include "myparser.h"


eCCDRowColType_T CRowColDesc::GetType( const char* type_desc )
{
	if(strncmp(type_desc,ROW_INDICATOR,strlen(ROW_INDICATOR))==0){
		return eCCDRow;
	}
	if(strncmp(type_desc,COL_INDICATOR,strlen(COL_INDICATOR))==0){
		return eCCDCol;
	}
	return eCCDUnknownType;
}

CRowColList::CRowColList( const char* filename )
{
	ReadFromFile( filename );
}

BOOL_T CRowColList::ReadFromFile( const char* filename )
{
	if(!MyFile::DoesFileExist( filename )){
		return FALSE;
	}

	MyIFile in(filename);
	const char* line=NULL;
	CMyStrTable items;	

	clear();	
	while(line = in.GetLine(TRUE)){
		MyParser pars = line;
		if( mystring::get_first_non_white( pars.c_str() )=='#' )
			continue;
		pars.GetItems( items, "= " ); // white space separated 
		if(items.size()>=2){
			CMyStrTable list;
			MyParser pars2=items[1].c_str();
			pars2.GetItems( list, "," ); // coma separated 
			eCCDRowColType_T type = CRowColDesc::GetType( items[0].c_str() );
			if(type!=eCCDUnknownType){
				for(int i=0;i<list.size();i++){
					if(atol(list[i].c_str())>=0){
						CRowColDesc tmp;
						tmp.type = type;
						tmp.num = atol( list[i].c_str() );
						AddUnique( tmp );
					}
				}
			}
		}
	}

	return TRUE;
}


CRowColDesc* CRowColList::Find( eCCDRowColType_T type, int num )
{
	CRowColList::iterator i;
	for(i=begin();i!=end();i++){
		if( i->type==type && i->num==num ){
			return (&(*i));
		}
	}
	return NULL;
}

void CRowColList::AddUnique( CRowColDesc& elem )
{
	if(!Find( elem.type, elem.num )){
		push_back( elem );
	}
}


CWindowList::CWindowList( const char* filename )
{
	ReadFromFile( filename );	
}
                                                                                
BOOL_T CWindowList::ReadFromFile( const char* filename )
{
	if(!MyFile::DoesFileExist( filename )){
		return FALSE;
	}

	MyIFile in(filename);
	const char* line=NULL;
	CMyStrTable items;	

	clear();	
	while(line = in.GetLine(TRUE)){
		MyParser pars = line;
		if( mystring::get_first_non_white( pars.c_str() )=='#' )
			continue;
		pars.GetItems( items ); // white space separated 
		if(items.size()>=4){
			CCDWindow tmp;
			tmp.m_LowX = atol( items[0].c_str() );			
			tmp.m_LowY = atol( items[1].c_str() );			
			tmp.m_UpX = atol( items[2].c_str() );			
			tmp.m_UpY = atol( items[3].c_str() );			
			push_back( tmp );
		}
	}
	return TRUE;	
}
