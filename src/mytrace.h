// mytrace.h: interface for the CMyTrace class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MYTRACE_H__5212D424_3F44_11D6_B636_B42D07C10000__INCLUDED_)
#define AFX_MYTRACE_H__5212D424_3F44_11D6_B636_B42D07C10000__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "myfile.h"
#include "mydate.h"

class BASELIB_EI CMyTrace : public MyOFile  
{
protected:
	long m_TrLevel;
public:
	CMyTrace(const char* filename,long TrLevel);
	virtual ~CMyTrace();
	long GetTrLevel(){ return m_TrLevel;}
	void SetTrLevel( int TrLevel ){ m_TrLevel = TrLevel; }
};

#ifdef _ENABLE_TRACING_

#define MYTRACE(TrObj,Cmd)  { TrObj << "[" << CMyDate::getdate() << "] " << Cmd <<"\n"; TrObj.Flush(); TrObj.Close();}
#define MYTRACE0(TrObj,Cmd) if(TrObj.GetTrLevel()>=0){ TrObj << "[" << CMyDate::getdate() << "] " << Cmd <<"\n"; TrObj.Flush(); TrObj.Close();}
#define MYTRACE1(TrObj,Cmd) if(TrObj.GetTrLevel()>=1){ TrObj << "[" << CMyDate::getdate() << "] " << Cmd <<"\n"; TrObj.Flush(); TrObj.Close();}
#define MYTRACE2(TrObj,Cmd) if(TrObj.GetTrLevel()>=2){ TrObj << "[" << CMyDate::getdate() << "] " << Cmd <<"\n"; TrObj.Flush(); TrObj.Close();}
#define MYTRACE3(TrObj,Cmd) if(TrObj.GetTrLevel()>=3){ TrObj << "[" << CMyDate::getdate() << "] " << Cmd <<"\n"; TrObj.Flush(); TrObj.Close();}
#define MYTRACE4(TrObj,Cmd) if(TrObj.GetTrLevel()>=4){ TrObj << "[" << CMyDate::getdate() << "] " << Cmd <<"\n"; TrObj.Flush(); TrObj.Close();}
#define MYTRACE5(TrObj,Cmd) if(TrObj.GetTrLevel()>=5){ TrObj << "[" << CMyDate::getdate() << "] " << Cmd <<"\n"; TrObj.Flush(); TrObj.Close();}
#else
#define MYTRACE0(TrObj,Cmd)
#define MYTRACE1(TrObj,Cmd)
#define MYTRACE2(TrObj,Cmd)
#define MYTRACE3(TrObj,Cmd)
#define MYTRACE4(TrObj,Cmd)
#define MYTRACE5(TrObj,Cmd)
#endif
#define ALWAYS_TRACE(TrObj,Cmd) { TrObj << "[" << CMyDate::getdate() << "] " << Cmd <<"\n"; TrObj.Flush(); TrObj.Close();}


class CCmnTrace : public CMyTrace
{
public :
   CCmnTrace();
};
      
extern CCmnTrace gCmnTrace;
   

#endif // !defined(AFX_MYTRACE_H__5212D424_3F44_11D6_B636_B42D07C10000__INCLUDED_)
