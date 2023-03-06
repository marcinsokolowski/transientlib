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
#include "myinifile.h"


CIniFile::CIniFile( const char* inifile ) : CCfgFile("") 
{
	if( MyFile::DoesFileExist(inifile) ){
		Init(inifile);
		m_szRequstedIniFile = inifile;
	}else{
/*		mystring szInHome = "$(HOME)/";
		szInHome << "." << inifile;*/

		// NEW 20050228 :
		mystring szInHome = inifile;
		szInHome.env2str();
		Init(szInHome.c_str());
		m_szRequstedIniFile = szInHome;
	}
}

CIniFile::~CIniFile()
{
}


BOOL_T CIniFile::Save( BOOL_T bForceNew )
{
	if(m_CfgTab.size() && ( strlen(m_CfgFile.GetFileName()) || bForceNew  ) ){
		mystring szFileName = m_szRequstedIniFile;
		if(strlen(m_CfgFile.GetFileName()))
			szFileName = m_CfgFile.GetFileName();

		MyOFile out( szFileName.c_str() );
		vector<CEnvVar>::iterator i;

	   for(i=m_CfgTab.begin();i!=m_CfgTab.end();i++){			
			out.Printf("%s=%s\n",i->szName.c_str(),i->szValue.c_str());
		}
		out.Flush();
		out.Close();			
		return TRUE;	
	}
	return FALSE;
}
