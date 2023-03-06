// mytrace.cpp: implementation of the CMyTrace class.
//
//////////////////////////////////////////////////////////////////////
#ifndef _UNIX
#include "stdafx.h"
#endif
#include "mytrace.h"
#include "cmncfg.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CMyTrace::CMyTrace(const char* filename,long TrLevel)
{
	m_TrLevel = TrLevel;
	m_FileName = filename;
	m_Attr = "a+";
#ifdef _DEBUG
	printf("Constructing trace file %s, trace level=%d\n",filename,TrLevel);
#endif
}

CMyTrace::~CMyTrace()
{

}



CCmnTrace gCmnTrace;

CCmnTrace::CCmnTrace()
: CMyTrace( "$(TRACEDIR)/$(USER)/cmn.trace", 
            atol(mystring("$(CMN_TRACE_LEVEL)").env2str()))
{
}		
