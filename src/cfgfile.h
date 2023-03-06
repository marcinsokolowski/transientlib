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
#ifndef _CFG_FILE_H__
#define _CFG_FILE_H__

#include "basestructs.h"
#include <vector>
#include "myfile.h"

using namespace std;

class BASELIB_EI CCfgFile 
{
public:
	CCfgFile(const char* filename);
	void Init(const char* filename);
	BOOL_T FindParam(const char* cfgcode);
	const char* GetParam(const char* cfgcode,BOOL_T bAllowNull=FALSE);
	const char* GetParamNoInit(const char* cfgcode,BOOL_T bAllowNull=FALSE);
	void SetParam( const char* cfgcode, const char* cfgval );
	void GetParamValues( vector<CEnvVar>& params );
	void GetParams( mystring& szParams );
	void Dump();
	void SaveToFile( const char* filename );
	vector<CEnvVar>& GetParamTable(){ return m_CfgTab; }		
	void Sort();
	static BOOL_T CompareEnvVarTables( vector<CEnvVar>& left, vector<CEnvVar>& right,
										 mystring& szDifferent, mystring& szOnlyInLeft,
                               mystring& szOnlyInRight );
	BOOL_T ReadCfgFile( MyIFile* pFile );
	BOOL_T InitCfgTab();
	BOOL_T IsOpened(){ return m_CfgFile.IsOpened(); }
	BOOL_T IsInitialized();
	const char* GetFileName();

	BOOL_T GetValue( const char* cfgcode, BOOL_T& boolVal );
	BOOL_T GetValue( const char* cfgcode, LONG_T& longVal );		
	BOOL_T GetValue( const char* cfgcode, double& doubleVal );		
	BOOL_T GetValue( const char* cfgcode, mystring& stringVal );		

	void SetValue( const char* cfgcode, BOOL_T boolVal );
	void SetValue( const char* cfgcode, LONG_T longVal );		
	void SetValue( const char* cfgcode, double doubleVal );		
	void SetValue( const char* cfgcode, const char* stringVal );		
	
	static void CopyParamsTab( vector<CEnvVar>& dest, vector<CEnvVar>& src );
protected:
	vector<CEnvVar> m_CfgTab;
	MyIFile m_CfgFile;
	BOOL_T m_bInitialized;
};



#endif
