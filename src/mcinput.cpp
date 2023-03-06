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
#include "mcinput.h"
#include <myparser.h>
#include <myfile.h>
#include <cfg.h>
#include <mytrace.h>
#include "ccd_globals.h"
// #include "ccd_image_creator.h"
#include <stdlib.h>

void (*CVaryParamDesc::ParamRefreshFunc1)(void)=NULL;
void (*CVaryParamDesc::ParamRefreshFunc2)(void)=NULL;
void (*CVaryParamDesc::ParamRefreshFunc3)(void)=NULL;
void (*CVaryParamDesc::ParamRefreshFunc4)(void)=NULL;

const char* CVaryParamDesc::GetParamTypeStr( eParamType type )
{
	if(type == ParamLong)
		return PARAM_LONG;
	if(type == ParamDouble)
		return PARAM_DOUBLE;
	if(type == ParamString)
		return PARAM_STRING;
	return NULL;
}

eParamType CVaryParamDesc::GetParamType( const char* szType )
{
	if(strcmp(szType,PARAM_LONG)==0)
		return ParamLong;
	if(strcmp(szType,PARAM_DOUBLE)==0)
		return ParamDouble;
	if(strcmp(szType,PARAM_STRING)==0)
		return ParamString;
	return ParamUnknown;	
}



CVaryParamDesc::CVaryParamDesc()
:m_Type(ParamLong),m_CheckType(NormalCheck),m_Index(0)
{}

void CVaryParamDesc::Reset()
{
	m_Index = 0;
	if(m_CheckType != OnlyListed){
		m_CurrVal = m_LowEnd;
	}else{
		m_CurrVal = m_ValuesToCheck[0];
	}
}

CVaryParamDesc::CVaryParamDesc( LONG_T low, LONG_T up, LONG_T step, mystring Name )
{
	mystring szLow,szUp,szStep;
	szLow << low;
	szUp << up;
	szStep << step;
	InitParams( szLow, szUp, szStep, Name, ParamLong ); 
}

BOOL_T CVaryParamDesc::NotEnd()
{
	BOOL_T bRet = TRUE;
	if(m_CheckType == NormalCheck){
		if(m_Type==ParamLong){
			LONG_T cur = atol( m_CurrVal.c_str() );
			LONG_T low = atol(m_LowEnd.c_str());
			LONG_T up  = atol(m_UpEnd.c_str());

			if( low <= up )
				bRet = ( low<=cur && cur<=up );
			else
				bRet = ( low>=cur && cur>=up );
				
			return bRet;
		}
		if(m_Type==ParamDouble){
			double cur = atof( m_CurrVal.c_str() );
			double low = atof(m_LowEnd.c_str());
			double up  = atof(m_UpEnd.c_str());

			if( low <= up )
				bRet = ( low<=cur && cur<=up );
			else
				bRet = ( low>=cur && cur>=up );

			return bRet;
		}
	}
	if(m_CheckType == OnlyListed){
		return (m_Index<m_ValuesToCheck.size());
	}
	return bRet;
}

CVaryParamDesc::CVaryParamDesc( double low, double up, double step, mystring Name )
{
	mystring szLow,szUp,szStep;
   szLow << low;
   szUp << up;
   szStep << step;
	InitParams( szLow, szUp, szStep, Name, ParamDouble );
} 

CVaryParamDesc::CVaryParamDesc( mystring low, mystring up, mystring step, mystring Name )
{
	InitParams( low, up, step, Name, ParamString );
}

void CVaryParamDesc::InitParams( mystring low, mystring up, mystring step, 
                                mystring Name, eParamType type )
{
	m_LowEnd = low;
	m_UpEnd = up;
	m_Step = step;	
	m_Type = type;
	m_CurrVal = m_LowEnd;
	m_ParamName = Name;

	if(strcmp(m_Step.c_str(),"CONST")==0){
		// only listed in m_LowEnd will be checked !!!
		MyParser pars = m_LowEnd.c_str();				
		pars.GetItems( m_ValuesToCheck );
		m_CheckType = OnlyListed;
		m_Index = 0;
		m_CurrVal = m_ValuesToCheck[0];
	}
	if(strcmp(m_Step.c_str(),"LIST")==0){
		m_CheckType = OnlyListed;
		m_Index = 0;
		m_CurrVal = m_ValuesToCheck[0];
		m_LowEnd = m_ValuesToCheck[0];
		m_UpEnd = m_ValuesToCheck[m_ValuesToCheck.size()-1];
	}

	// Validatation of parametrs :
	ValidateParam();
}

BOOL_T CVaryParamDesc::CheckIfIntegral( const char* szValue )
{
	char asFloat[100];
	char asLong[100];

	sprintf(asFloat,"%f.10",atof(szValue));
	sprintf(asLong,"%f.10",(double)atol(szValue));

	if ( strcmp(asFloat,asLong) ){
		return FALSE;
	}
	return TRUE;
}

BOOL_T CVaryParamDesc::ValidateParam()
{
	// validate type :
	if(m_Type == ParamLong){
		const char* pStr="";
		if(!CheckIfIntegral( (pStr=m_LowEnd.c_str()) ) || 
			!CheckIfIntegral( (pStr=m_UpEnd.c_str()) ) ||
			!CheckIfIntegral( (pStr=m_Step.c_str()) ) ){
			mystring szErrMsg;
			szErrMsg << "Parameter " << m_ParamName << " of type long has floating value : " 
						<< pStr << " please verify mcinput.txt file, exiting";			
			MYTRACE1(gCmnTrace,szErrMsg.c_str());
			printf("%s\n",szErrMsg.c_str());
			exit(-1);
		}
	}
	return TRUE;
}
	
void CVaryParamDesc::Init( CMyStrTable& tab )
{
	eParamType type = ParamLong;
	if(tab.size()>=4){
		if(tab.size()>=5){
			if(strcmp(tab[4].c_str(),PARAM_DOUBLE)==0){
				type = ParamDouble;
			}
			if(strcmp(tab[4].c_str(),PARAM_STRING)==0){
				type = ParamString;
			}
		}
		InitParams( tab[1], tab[2], tab[3], tab[0], type );
	}
}

CVaryParamDesc& CVaryParamDesc::operator++(int i)
{
	if(m_CheckType == NormalCheck){
		if(m_Type==ParamLong){
			LONG_T tmp;
			tmp = atol(m_CurrVal.c_str());
			tmp += atol(m_Step.c_str());
			m_CurrVal = tmp;
		}
		if(m_Type==ParamDouble){
			double tmp;
			tmp = atof(m_CurrVal.c_str());
			tmp += atof(m_Step.c_str());
			
			// only two digits of precision here 
			// not default as in = operator :
			// ((BaseString&)m_CurrVal) = tmp;
			char szCurrVal[100];
			sprintf(szCurrVal,"%.2f",tmp);
			m_CurrVal = szCurrVal;
		}
		if(m_Type == ParamString){
			// not handled - yet remains unchanged		
		}
	}
	if(m_CheckType == OnlyListed){
		m_Index++;		
		if(m_Index<m_ValuesToCheck.size())		
			m_CurrVal = m_ValuesToCheck[m_Index];		
	}
	return (*this);
}

void CVaryParamDesc::UpdateParams()
{
	mystring szValue = m_CurrVal;
	gCCDParams.SetParam( m_ParamName.c_str() , szValue.c_str() );

	// this two below will be included in single global function - RefreshParams
	RefreshStaticParamsInCcdlib();

	if(ParamRefreshFunc1)
		(*ParamRefreshFunc1)();	
	if(ParamRefreshFunc2)
		(*ParamRefreshFunc2)();	
	if(ParamRefreshFunc3)
		(*ParamRefreshFunc3)();	
	if(ParamRefreshFunc4)
		(*ParamRefreshFunc4)();	
}


CVaryParamDescTab::CVaryParamDescTab()
{}

void CVaryParamDescTab::UpdateParams()
{
	CVaryParamDescTab::iterator i;
	for(i=begin();i!=end();i++){
		i->UpdateParams();
	}
}

	
CInputMC::CInputMC(const char* fname, BOOL_T bAutoRead)
{
	if(bAutoRead)
		ReadInputFile( fname );
}


BOOL_T CInputMC::ReadInputFile( const char* fname )
{
	char* param;
	MyIFile input;

	m_bOK = FALSE;

	if(!input.Open(fname,"r",FALSE)){
		return FALSE;
	}

	while( (param = input.GetLine()) && (param[0]=='#' || param[0]=='\n') ){}
	MyParser Mags = param;
	Mags.GetItems( m_Magnitudes );
	m_ParamsTab.clear();	
	m_ConstParamsTab.clear();
	while( (param = input.GetLine()) && (param[0]=='#' || param[0]=='\n') ){}
	
	// reading parameters to vary :
	while( param ){
		if(param[0]!='#' && param[0]!='\n'){
			MyParser pars = param;
			CMyStrTable tab;
			if(pars.GetItems( tab, SEPARATOR ) >= 4){
				CVaryParamDesc tmp;
				tmp.Init( tab );
				
				if(strcmp( tmp.m_LowEnd.c_str(), tmp.m_UpEnd.c_str() ))
					m_ParamsTab.push_back(tmp);						
				else
					m_ConstParamsTab.push_back(tmp);
			}else{
				if(tab.size()==3){
					// there is parameter :
					// NAME VAL1,VAL2,VAL3 TYPE
					mystring szParName = tab[0];
					MyParser szValues = tab[1].c_str();
					mystring szType = tab[2];
					CVaryParamDesc tmp;
					szValues.GetItems(tmp.m_ValuesToCheck,",");
					//if(tmp.m_ValuesToCheck.size()>1){
						tmp.InitParams( "0", "0", "LIST", szParName, CVaryParamDesc::GetParamType(szType.c_str()) );					
						m_ParamsTab.push_back(tmp);					
					//}else{
						// single value listed - tread as const :
						
					//}
				}
			}
		}		
		param = input.GetLine();
	}
	m_bOK = TRUE;

	gCCDParams.m_szMcInputFile = fname;

	return TRUE;
}

void CInputMC::GetVaryParamsValues( CMyStrTable& ParamTab, mystring& szParams )
{
	vector<CVaryParamDesc>::iterator i;

	ParamTab.clear();
	szParams = "";
	for(i=m_ParamsTab.begin();i!=m_ParamsTab.end();i++){
		mystring szValue;
		szValue << i->m_ParamName << "=" << gCCDParams.GetParam( i->m_ParamName.c_str() );
		ParamTab.Add( szValue );
		szParams << szValue << "\n";
	}	

	for(i=m_ConstParamsTab.begin();i!=m_ConstParamsTab.end();i++){
		mystring szValue;
		szValue << i->m_ParamName << "=" << gCCDParams.GetParam( i->m_ParamName.c_str() );
		ParamTab.Add( szValue );
		szParams << szValue << "\n";
	}

}

void CInputMC::GetVaryParams( CMyStrTable& ParamTab )
{
	vector<CVaryParamDesc>::iterator i;

	ParamTab.clear();
	for(i=m_ParamsTab.begin();i!=m_ParamsTab.end();i++){
		ParamTab.Add( i->m_ParamName );
	}	
}
