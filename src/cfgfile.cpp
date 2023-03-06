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
#include <algorithm>
#include <functional>
#include "cfgfile.h"
#include "cexcp.h"
#include "myparser.h"
#include "mytrace.h"
#include "mycmnglobals.h"
#include "mymacros.h"

CCfgFile::CCfgFile(const char* filename)
: m_bInitialized(FALSE)
{
	Init( filename );
}

void CCfgFile::Init(const char* filename)
{
	if(filename && filename[0]){
		mystring szFileName;
		szFileName = filename;
		szFileName.env2str();
		if (MyFile::DoesFileExist(szFileName.c_str())){
			// szFileName = DEFAULT_CFG_FILE;
			m_CfgFile.Open(szFileName.c_str());
			InitCfgTab();
		}
	}

}

const char* CCfgFile::GetFileName()
{
	return m_CfgFile.GetFileName();
}

BOOL_T CCfgFile::FindParam(const char* cfgcode)
{
	vector<CEnvVar>::iterator i;
	for(i=m_CfgTab.begin();i!=m_CfgTab.end();i++){
		if (strcmp(i->szName.c_str(),cfgcode)==0){
			return TRUE;
		}
	}
	return FALSE;
}

const char* CCfgFile::GetParam(const char* cfgcode,BOOL_T bAllowNull)
{
	vector<CEnvVar>::iterator i;

	InitCfgTab();
	for(i=m_CfgTab.begin();i!=m_CfgTab.end();i++){
		if (strcmp(i->szName.c_str(),cfgcode)==0){
			MYTRACE3(gCmnTrace,"Requested parameter " << cfgcode << "=" << i->szValue.c_str());
			return i->szValue.c_str();
		}
	}
	if (!bAllowNull){
		Assert(FALSE,"Parameter %s not defined",cfgcode);
	}
	return NULL;
}

const char* CCfgFile::GetParamNoInit(const char* cfgcode,BOOL_T bAllowNull)
{
	vector<CEnvVar>::iterator i;

	for(i=m_CfgTab.begin();i!=m_CfgTab.end();i++){
		if (strcmp(i->szName.c_str(),cfgcode)==0){
			MYTRACE3(gCmnTrace,"Requested parameter " << cfgcode << "=" << i->szValue.c_str());
			return i->szValue.c_str();
		}
	}
	if (!bAllowNull){
		Assert(FALSE,"Parameter %s not defined",cfgcode);
	}
	return NULL;
}

BOOL_T CCfgFile::GetValue( const char* cfgcode, BOOL_T& boolVal )
{
	const char* szVal = GetParam( cfgcode , TRUE );
	if(szVal && szVal[0]){
		boolVal = (atol(szVal)>0);
		return TRUE;
	}
	return FALSE;
}

BOOL_T CCfgFile::GetValue( const char* cfgcode, LONG_T& longVal )
{
	const char* szVal = GetParam( cfgcode , TRUE );
	if(szVal && szVal[0]){
		longVal = atol(szVal);
		return TRUE;
	}
	return FALSE;
}


BOOL_T CCfgFile::GetValue( const char* cfgcode, double& doubleVal )
{
	const char* szVal = GetParam( cfgcode , TRUE );
	if(szVal && szVal[0]){
		doubleVal = atof(szVal);
		return TRUE;
	}
	return FALSE;
}

BOOL_T CCfgFile::GetValue( const char* cfgcode, mystring& stringVal )
{
	const char* szVal = GetParam( cfgcode , TRUE );
	if(szVal && szVal[0]){
		stringVal = szVal;
		return TRUE;
	}
	return FALSE;
}


void CCfgFile::SetValue( const char* cfgcode, BOOL_T boolVal )
{
	mystring szTmp;
	szTmp << boolVal;
	SetParam(cfgcode,szTmp.c_str());
}

void CCfgFile::SetValue( const char* cfgcode, LONG_T longVal )
{
	mystring szTmp;
	szTmp << longVal;
	SetParam(cfgcode,szTmp.c_str());
}

void CCfgFile::SetValue( const char* cfgcode, double doubleVal )
{
	char szTmp[100];
	sprintf(szTmp,"%.4f",doubleVal);
	SetParam(cfgcode,szTmp);
}
void CCfgFile::SetValue( const char* cfgcode, const char* stringVal )
{
	mystring szTmp(stringVal);
	SetParam(cfgcode,szTmp.c_str());
}

void CCfgFile::SetParam( const char* cfgcode, const char* cfgval )
{
	vector<CEnvVar>::iterator i;

	BOOL_T bFound=FALSE;
	for(i=m_CfgTab.begin();i!=m_CfgTab.end();i++){
		if (strcmp(i->szName.c_str(),cfgcode)==0){
			i->szValue = cfgval;
			bFound = TRUE;
			break;
		}
	}
	if(!bFound){
		CEnvVar tmp;
		tmp.szName = cfgcode;
		tmp.szValue = cfgval;
		m_CfgTab.push_back(tmp);
	}
}

BOOL_T CCfgFile::IsInitialized()
{
	return m_bInitialized;
}

BOOL_T CCfgFile::ReadCfgFile( MyIFile* pFile )
{
	int count=0;
	if (pFile->IsOpened()){
		_TRACE_PRINTF_0("Reading config file : %s ...\n",pFile->GetFileName());

		const char* pLine;
		mystring szParams;
		while(pLine = pFile->GetLine(TRUE)){
			if (strlen(pLine)==0 || pLine[0]=='\n')
				continue;
			// skip comments :
			if (pLine[0]=='#' || strncmp(pLine,"//",2)==0)
				continue;
			if(strncmp(pLine,"%LOAD%",6)==0 || strncmp(pLine,"%SAFELOAD%",10)==0 ){
				// loading include file :
				mystring szTmp;
				if( strncmp(pLine,"%LOAD%",6)==0 )
			      szTmp = pLine+7;
				else
					szTmp = pLine+11;
				mystring szFileName = szTmp.TrimApostro('"');
		      szFileName.env2str();
		      if (MyFile::DoesFileExist(szFileName.c_str())){
      		   // szFileName = DEFAULT_CFG_FILE;
					MyIFile includeFile;
		         includeFile.Open(szFileName.c_str());
		         ReadCfgFile( &includeFile );
	      	}else{
					if( strncmp( pLine,"%LOAD%", 6 )==0 ){
						Assert(FALSE,"Could not find include file : %s\n",szFileName.c_str());
					}else{
						if( gSafeLoadWarning ){
							_TRACE_PRINTF_0("could not load file : %s\n",szFileName.c_str());
							_TRACE_PRINTF_0("WARN : please verify if this can be accepted !!!\n");
							_TRACE_PRINTF_0("IGNORED command : %s\n",pLine);
						}
					}
				}
				continue;
			}

			MyParser parser = pLine;
			CEnvVar tmp;

			/*if( strstr(pLine,"dark") ){
				printf("odo");
			}*/

			if(!parser.GetVarAndValueFromLine(tmp.szName,tmp.szValue))
				continue;
			tmp.szValue.env2str();

			if(strlen(tmp.szName.c_str())==0 && strlen(tmp.szValue.c_str())==0)
				continue;
			
			// if(!FindParam(tmp.szName.c_str())){
				// adding only once same param - first in the file
				// takes over :
				m_CfgTab.push_back(tmp);
			// }
			szParams << pLine << "\n";			
			count++;
		}
		pFile->Close();
		MYTRACE2(gCmnTrace,"System parameters read :\n" << szParams);
		m_bInitialized = TRUE;
	}
	return count;
	
}

BOOL_T CCfgFile::InitCfgTab()
{
	if (m_CfgTab.size()==0 || m_CfgFile.IsOpened()){
		ReadCfgFile( &m_CfgFile );
	}	
	return m_CfgTab.size();
}

void CCfgFile::GetParamValues( vector<CEnvVar>& params )
{
	params.clear();
	vector<CEnvVar>::iterator i;

   for(i=m_CfgTab.begin();i!=m_CfgTab.end();i++){
		params.push_back( *i );
	}
}


void CCfgFile::GetParams( mystring& szParams )
{
	vector<CEnvVar>::iterator i;
	szParams = "";

	for(i=m_CfgTab.begin();i!=m_CfgTab.end();i++){
		szParams << i->szName << "=" << i->szValue << "\n";
	}
}

void CCfgFile::SaveToFile( const char* filename )
{
	vector<CEnvVar>::iterator i;
	mystring szMsg(2000,2000);
	for(i=m_CfgTab.begin();i!=m_CfgTab.end();i++){
		szMsg << i->szName.c_str() << "=" << i->szValue.c_str() << "\n";
	}
	mystring szName=filename;
	szName.env2str();
	MyOFile out( szName.c_str() );
	out.Printf("%s",szMsg.c_str());
}

void CCfgFile::Dump()
{
	vector<CEnvVar>::iterator i;
	mystring szMsg;
	for(i=m_CfgTab.begin();i!=m_CfgTab.end();i++){
		szMsg << i->szName.c_str() << "=" << i->szValue.c_str() << "\n";
	}
	printf("\n\nPARAMETER VALUES :\n%s",szMsg.c_str());
	printf("#############################\n");
}

void CCfgFile::Sort()
{
	if(m_CfgTab.size()){
		sort(m_CfgTab.begin(),m_CfgTab.end(),less<CEnvVar>());		
	}
}

BOOL_T CCfgFile::CompareEnvVarTables( vector<CEnvVar>& left, vector<CEnvVar>& right, 
												  mystring& szDifferent, mystring& szOnlyInLeft,
												  mystring& szOnlyInRight )
{	
	BOOL_T bRet = TRUE;

	szDifferent = "";
	szOnlyInLeft = "";
	szOnlyInRight = "";

	int l = 0;
	int r = 0;

	if(left.size()!=right.size() )
		bRet = FALSE;

	while( l<left.size() && r<right.size() ){
		if(left[l] == right[r]){
			l++;
			r++;
		}else{
			bRet = FALSE;
			if(left[l].szName==right[r].szName){
				szDifferent << left[l].szName << " : " << left[l].szValue 
								<< "     " << right[r].szValue << "\n";				
				l++;
				r++;
			}else{
				if(left[l].szName<right[r].szName){
					while(left[l].szName<right[r].szName && l<left.size()){
						szOnlyInLeft << left[l].szName << " : " << left[l].szValue << "\n";
						l++;
					}
				}else{
					if(left[l].szName>right[r].szName){
						while(left[l].szName>right[r].szName && r<right.size()){
							szOnlyInRight << right[r].szName << "=" << right[r].szValue << "\n";
							r++;
						}
					}
				}
			}
		}			
	}

	while(l<left.size()){
		szOnlyInLeft << left[l].szName << "=" << left[l].szValue << "\n";
		l++;
	}

	while(r<right.size()){
		szOnlyInRight << right[r].szName << "=" << right[r].szValue << "\n";
		r++;
	}


	return bRet;
}

void CCfgFile::CopyParamsTab( vector<CEnvVar>& dest, vector<CEnvVar>& src )
{
	dest.clear();
	vector<CEnvVar>::iterator i;
	for(i=src.begin();i!=src.end();i++){
		dest.push_back( *i );
	}
}


