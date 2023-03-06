// myenv.h: interface for the CMyEnv class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MYENV_H__9D2077ED_B738_11D5_B636_382E07C10000__INCLUDED_)
#define AFX_MYENV_H__9D2077ED_B738_11D5_B636_382E07C10000__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include <vector>
#include "mystring.h"
#include "basedefines.h"
#include "mykeytab.h"

using namespace std;

#define ENV_FILE_NAME "env.cfg"


class BASELIB_EI CMyEnv  
{
public:
	CMyEnv();
	virtual ~CMyEnv();
	static const char* mygetenv(const char* name);
	static void GetEnvList(mystring& szEnv);
protected :
	static void InitEnv();
	static BOOL_T m_bInitialized;
	static CKeyTab  m_EnvList;
};

#endif // !defined(AFX_MYENV_H__9D2077ED_B738_11D5_B636_382E07C10000__INCLUDED_)
