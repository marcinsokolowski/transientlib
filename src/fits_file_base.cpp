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
#include "fits_file_base.h"
#include "ccd_fits_header_defs.h"
#include "fitslib_globals.h"
#include "mymacros.h"
#include "cfg.h"
#include "myparser.h"
#include "mystrtable.h"
#include <fitsio.h>


#include <cmath>
#include <vector>
#include <string>
// The library is enclosed in a namespace.

using namespace std;

void fits_file_extend(char* s,char c,int l,int line_size)             /* fill string 's' to length 'l'  with 'c' */
{
	register int i=strlen(s);
	while(i < l){
		s[i] = c;
		i++;
	}
	if(i < line_size){
		s[i] = '\0';                      /* s[l] = '\0';   */
	}
}

int get_fits_format(int type)
{
	return LONG_IMG;
}

int get_fits_format(unsigned short type)
{
	return USHORT_IMG; 
}


int get_fits_format(unsigned char type)
{
	return BYTE_IMG; 
}

int get_fits_format(short type)
{
	return SHORT_IMG; 
}

int get_fits_format(long type)
{
	return LONG_IMG; 
}

int get_fits_format(float type)
{
	return FLOAT_IMG; 
}

int get_fits_format(double type)
{
	return DOUBLE_IMG; 
}

const char* get_image_type(int type)
{
	 switch (type) {
      case BYTE_IMG:
         return "BYTE_IMG";
         break;
      case SHORT_IMG:
         return "SHORT_IMG";
         break;
      case LONG_IMG:
         return "LONG_IMG";
         break;
      case FLOAT_IMG:
         return "FLOAT_IMG";
         break;
      case DOUBLE_IMG:
         return "DOUBLE_IMG";
         break;
   }
}


// globals to avoid incuding <vector> by every file including fits_file.h
int GetKnownHeaderList( std::vector<string>& keyNames )
{
	keyNames.clear();
	int i=0;
	while(fitsHeaderKeyDefTab[i].keyType!=eEND && fitsHeaderKeyDefTab[i].szKeyName){
		// string newElem = all_header_list[i];
		string newElem = fitsHeaderKeyDefTab[i].szKeyName;
		// printf("KEY = %s %s\n",newElem.c_str(),fitsHeaderKeyDefTab[i].szKeyComment);
		keyNames.push_back(newElem);
		i++;
	}
	return keyNames.size();	
}


CFITSFileBase::CFITSFileBase()
{
}

void CFITSFileBase::Init()
{
	if(!m_bStaticInitialized){
		GetKnownHeaderList( m_HeaderKeys );
		const char* szAdd = GetGlobalParam("CCD_ADDITIONAL_HEADER_KEYS",TRUE);
		if(szAdd && szAdd[0]){
			MyParser addit = szAdd;
			CMyStrTable table;
			addit.GetItems( table );
			for(int i=0;i<table.size();i++){
				string newElem = table[i].c_str();
				m_HeaderKeys.push_back(newElem);
			}
		}
		m_bStaticInitialized = TRUE;
	}
}


void CFITSFileBase::AddHdu(const char* key,const char* val)
{
	m_HduList.Add(key,val);
}

void CFITSFileBase::SetHdu(const char* key,const char* val)
{
	m_HduList.Set(key,val);
}


void CFITSFileBase::Clear()
{
	m_HduList.Clear();
}



void CFITSFileBase::SortHeaderBySections()
{
	int nSectionCount=eSecOther+1;
	CSafeKeyTab tmpSecList[eSecOther+1];
	CSafeKeyTab tmpList;
	int sec_len = strlen( SECTION );
	int pos=0;
	for(register int h=0;h<m_HduList.GetCount();h++){
		sFITSHeaderDesc* keyDesc = GetKeyDesc( m_HduList[h].szName.c_str() , pos  );
		if( keyDesc ){
			if( keyDesc->eSection>=0 && keyDesc->eSection<nSectionCount 
				 && strncmp( m_HduList[h].szName.c_str(), SECTION, sec_len ) ){
				if( strlen( m_HduList[h].szComment.c_str() )==0 && keyDesc->szKeyComment ){
					m_HduList[h].szComment = keyDesc->szKeyComment;
				}
				tmpSecList[ keyDesc->eSection ].Add( m_HduList[h] );
			}
		}
	}

	for(int i=1;i<=eSecOther;i++){
		mystring szSecName;
		szSecName << SECTION << i;		
		mystring szSecVal;
		szSecVal << "---------- " << GetSectionName((eHeaderSectionName)i) << " ----------";

		if( tmpSecList[i].GetCount() ){
			tmpList.Add( szSecName.c_str(), szSecVal.c_str() );
			for(register int h=0;h<tmpSecList[i].GetCount();h++){
				tmpList.Add( tmpSecList[i][h] );			
			}
		}
	}


	/*CSafeKeyTab tmpList;
	
	for(int i=1;i<=eSecOther;i++){
		mystring szSecName;
		szSecName << SECTION << i;		
		mystring szSecVal;
		szSecVal << "---------- " << GetSectionName((eHeaderSectionName)i) << " ----------";

		int nInSection=0;
		for(register int h=0;h<m_HduList.GetCount();h++){
			sFITSHeaderDesc* keyDesc = GetKeyDesc( m_HduList[h].szName.c_str() );

			if( keyDesc && keyDesc->eSection == (eHeaderSectionName)i ){
				if( nInSection==0 ){
					tmpList.Add( szSecName.c_str(), szSecVal.c_str() );
				}
				tmpList.Add( m_HduList[h] );
				nInSection++;
			}
		}
	}*/	

	/*mystring szSecName;
   szSecName << SECTION << (int)eSecUnknown;
   mystring szSecVal;
   szSecVal << "---------- " << GetSectionName(eSecUnknown) << " ----------";
   tmpList.Add( szSecName.c_str(), szSecVal.c_str() );*/

	for(register int h=0;h<m_HduList.GetCount();h++){
		sFITSHeaderDesc* keyDesc = GetKeyDesc( m_HduList[h].szName.c_str() );
		if(!keyDesc || keyDesc->eSection==eSecUnknown){
			if(!tmpList.Find( m_HduList[h].szName.c_str() )){
				tmpList.Add( m_HduList[h] );
			}
		}	
	}

	const char* szFILE = tmpList.getKeyVal( FILENAME );
	if( szFILE && strlen(szFILE)>30 && strstr(szFILE,"/") ){
		mystring tmp=szFILE;
		mystring newFile = getfname( tmp );
		tmpList.Set( FILENAME, newFile.c_str() );
	}

	m_HduList.Clear();
	m_HduList = tmpList;
}

const char* CFITSFileBase::getKeyVal( const char* szKey )
{
	return m_HduList.getKeyVal( szKey );
}
