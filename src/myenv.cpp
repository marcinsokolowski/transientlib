// myenv.cpp: implementation of the CMyEnv class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _UNIX
#include "stdafx.h"
#endif
#include "myenv.h"
#include "myfile.h"
#include "myparser.h"
#include "basestructs.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

BOOL_T CMyEnv::m_bInitialized=FALSE;
CKeyTab CMyEnv::m_EnvList;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CMyEnv::CMyEnv()
{

}

CMyEnv::~CMyEnv()
{

}

void CMyEnv::InitEnv()
{
	if(!m_bInitialized){
		MyIFile envfile(ENV_FILE_NAME,FALSE);
		if (envfile.IsOpened()){
			const char* pLine;
			while(pLine = envfile.GetLine(TRUE)){
				MyParser parser = pLine;
				CEnvVar tmp;
				parser.GetVarAndValue(tmp.szName,tmp.szValue);
				m_EnvList.push_back(tmp);
			}
		}
		m_bInitialized = TRUE;
		vector<CEnvVar>::iterator i;
		for(i=m_EnvList.begin();i!=m_EnvList.end();i++){
			mystring szTmp = i->szValue.c_str();
			szTmp.env2str();
			i->szValue = szTmp.c_str();
		}
	}
}


const char* CMyEnv::mygetenv(const char* name)
{
	vector<CEnvVar>::iterator i;
	InitEnv();
	for(i=m_EnvList.begin();i!=m_EnvList.end();i++){
		if (STRCASECMP(i->szName.c_str(),name)==0){
			return i->szValue.c_str();
		}
	}
	const char* szRet = ::getenv(name);
	return szRet;
}

void CMyEnv::GetEnvList(mystring& szEnv) 
{ 
	vector<CEnvVar>::iterator i;
	InitEnv();
	szEnv = "";
	for(i=m_EnvList.begin();i!=m_EnvList.end();i++){
		szEnv << i->szName << "=" << i->szValue << "\n";
	}
}
