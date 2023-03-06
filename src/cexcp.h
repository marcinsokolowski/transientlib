#ifndef _C_EXCP_H__
#define _C_EXCP_H__
#include "mytypes.h"
#include "basedefines.h"
#include "mystring.h"

#define ASSERT_FILE_NAME "$(TRACEDIR)/lib.trace"

class BASELIB_EI CExcBasic 
{
protected:
	BaseString m_Msg;
	BaseString m_szCode;
	eERRCODE_T m_Errcode;
	int m_errno;
public:	
	CExcBasic(const char* pMsg, const char* szCode);
	CExcBasic(const char* pMsg=NULL,eERRCODE_T errcode=ERR_UNKNOWN);
	CExcBasic(int _errno,const char* pMsg=NULL,eERRCODE_T errcode=ERR_UNKNOWN);	
	const char* GetMsg() const;
	eERRCODE_T GetErrcode() const;
	int GetErrno() const{ return m_errno;}
};

class BASELIB_EI CExcFile : public CExcBasic
{
public:
	CExcFile(int _errno,const char* pMsg=NULL,eERRCODE_T errcode=ERR_UNKNOWN)
		: CExcBasic(_errno,pMsg,errcode) {};

};

class BASELIB_EI CMemAlloc
{
   long m_size;	
   int m_errno;
	public :
		CMemAlloc(int _errno,int size) : m_errno(_errno),m_size(size){}
	   mystring GetReport();
};

void PrintToTraceFile( const char* fmt,...);
void BASELIB_EI AssertNULLFunc(const void* param,BaseString str);
void BASELIB_EI AssertFunc(BOOL_T bAssert,const char* fmt,...);

#define AssertNULL(param) AssertNULLFunc(param,BaseString()<<__FILE__<<":"<<__LINE__)
#define Assert AssertFunc


#endif
