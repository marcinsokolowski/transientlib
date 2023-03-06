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
#ifndef _MC_INPUT_FILE_H__
#define _MC_INPUT_FILE_H__

#include <mytypes.h>
#include <mystring.h>
#include <mystrtable.h>
#include <stdlib.h>
#include <vector>

using namespace std;

#define PARAM_LONG   "LONG"
#define PARAM_DOUBLE "DOUBLE"
#define PARAM_STRING "STRING"
#define DEFAULT_INPUT "mcinput.txt"
#define SEPARATOR "\t"

enum eParamType  { ParamUnknown=-1, ParamLong, ParamDouble, ParamString };

enum eCheckType { NormalCheck=0, OnlyListed };

class CVaryParamDesc
{
public :

	// This function pointers will store pointers to functions
	// from other liraries to be called when refreshing parameters !
	static void (*ParamRefreshFunc1)(void);
	static void (*ParamRefreshFunc2)(void);
	static void (*ParamRefreshFunc3)(void);
	static void (*ParamRefreshFunc4)(void);


	static const char* GetParamTypeStr( eParamType type );
	static eParamType GetParamType( const char* szType );


	CVaryParamDesc();
	
	CVaryParamDesc( LONG_T low, LONG_T up, LONG_T step, mystring Name );
	CVaryParamDesc( double low, double up, double step, mystring Name );
	CVaryParamDesc( mystring low, mystring up, mystring step, mystring Name );

	void InitParams( mystring low, mystring up, mystring step, mystring Name ,eParamType type );
	void Init( CMyStrTable& tab );	
	void Reset();
		
	BOOL_T CheckIfIntegral( const char* szValue );
	BOOL_T ValidateParam();			 

	void UpdateParams();
	

	mystring m_ParamName;
	eParamType m_Type;
	mystring m_LowEnd;
	mystring m_UpEnd;
	mystring m_Step;
	mystring m_CurrVal;

	// member for OnlyListed type - values to test :
	CMyStrTable m_ValuesToCheck;
	int m_Index;

	// type of check :
	eCheckType m_CheckType;

	CVaryParamDesc& operator++(int);
	BOOL_T NotEnd();
	
	operator LONG_T(){ return atol(m_CurrVal.c_str()); }
	operator double(){ return atof(m_CurrVal.c_str()); }
	
};

class CVaryParamDescTab : public vector<CVaryParamDesc>
{
public:
	CVaryParamDescTab();
	void UpdateParams();	
};

//-----------------------------------------------------------------
// Imput File Syntax :
//  1 line - list of magnitudes to be used as sample objects
//  parameters to vary in form :
//  NAME  LOWER  UPPER   STEP    TYPE
//  Example :
//
//  1,2,3,4
//  CCD_PIXELS_AROUND_TO_CONFIRM  1   10   1   LONG
//
//
//  special values : 
//
//   STEP=-1 - means that only values listed in LOWER column will be examined 
//
//-----------------------------------------------------------------
class CInputMC
{
protected:
	CVaryParamDescTab m_ParamsTab;
	CVaryParamDescTab m_ConstParamsTab;
	CMyStrTable m_Magnitudes;
	BOOL_T m_bOK;
public:

	
	CInputMC(const char* fname= DEFAULT_INPUT, BOOL_T bAutoRead=TRUE);
	
	BOOL_T ReadInputFile( const char* fname = DEFAULT_INPUT  );
	BOOL_T IsOK(){ return m_bOK; }
	CVaryParamDesc& GetParam( int i ){ return m_ParamsTab[i]; }
	CMyStrTable& GetMagTab() { return m_Magnitudes; }

	int GetParamCount(){ return m_ParamsTab.size(); }
	void GetVaryParams( CMyStrTable& ParamTab );
	void GetVaryParamsValues( CMyStrTable& ParamTab, mystring& szParams );
	CVaryParamDescTab& GetVaryParamsDesc() { return m_ParamsTab; }		
	CVaryParamDescTab& GetConstParamsDesc() { return m_ConstParamsTab; }
};



#endif

