#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "cexcp.h"
#include "myfile.h"
#include "mydate.h"
#include "mymacros.h"
#include "mycmnglobals.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void PrintToTraceFile( const char* fmt,...)
{
		mystring szFmt;
		szFmt << CMyDate::getdate() << " : " << fmt << "\n";
		MyOFile o_file(ASSERT_FILE_NAME,"a+");
		va_list plist;
		va_start(plist,fmt);
		o_file.VPrintf(szFmt.c_str(),plist);
		char buff[1000];
		o_file.Close();
		if (vsprintf(buff,szFmt.c_str(),plist)>=1000){
			o_file.Printf("Cirtical error buffer in Assert to small exiting");
			exit(0);
		}
		printf("%s\n",buff);
}

void AssertNULLFunc(const void* param,BaseString str) 
{
	if (!param){
		BaseString szTmp;
		szTmp << CMyDate::getdate() << " : Unexpected NULL pointer at : " << str.c_str() << "\n";
		printf("%s\n",szTmp.c_str());
		MyOFile o_file(ASSERT_FILE_NAME,"a+");
		o_file.Printf(szTmp.c_str());
		throw CExcBasic(szTmp.c_str());
	}
}


void AssertFunc(BOOL_T bAssert,const char* fmt,...)
{
	if(!bAssert){
		mystring szFmt;
		szFmt << CMyDate::getdate() << " : " << fmt << "\n";
		va_list plist;
      va_start(plist,fmt);

		// pritnf to STDOUT
		char buff[1000];
		if (vsprintf(buff,szFmt.c_str(),plist)>=1000){
			printf("Cirtical error buffer in Assert to small exiting");
			exit(-1);
		}
		printf("%s\n",buff);

		// to FILE :
		MyOFile o_file(ASSERT_FILE_NAME,"a+");
		o_file.VPrintf(szFmt.c_str(),plist);
		o_file.Close();
		exit(-1);
		// throw CExcBasic(buff);
	}
}

CExcBasic::CExcBasic(const char* pMsg, const char* szCode)
: m_Errcode( ERR_UNKNOWN ), m_errno(0)
{
	if(gPrintfLevel>=0){
		_TRACE_PRINTF_0("Exception %s occured due to : %s\n",szCode,pMsg);
		fflush(0);
	}

	m_Msg = pMsg;
	m_szCode = szCode;
}

CExcBasic::CExcBasic(int _errno,const char* pMsg,eERRCODE_T errcode)
{
	if(gPrintfLevel>=0){
		_TRACE_PRINTF_0("Exception errno=%d occured due to : %s\n",_errno,pMsg);
		_TRACE_PRINTF_0("Errno desc : %s\n",strerror(_errno));
		fflush(0);
	}

	m_Msg = pMsg;
	m_Errcode = ERR_UNKNOWN;
	m_errno = _errno;
}
	

CExcBasic::CExcBasic(const char* pMsg,eERRCODE_T errcode)
{
	m_Msg = pMsg;
	m_Errcode = ERR_UNKNOWN;
}


const char* CExcBasic::GetMsg() const
{
	return m_Msg.c_str();
}

eERRCODE_T CExcBasic::GetErrcode() const
{
	return m_Errcode;
}

mystring CMemAlloc::GetReport()
{
	mystring szMsg;
	szMsg << "Allocating : " << (int)m_size << " bytes\n";
	szMsg << "Errno      : " << m_errno;
	// szMsg << "Message    : " << m_Msg;
	return szMsg;
}
