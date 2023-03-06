//
//////////////////////////////////////////////////////////////////////

#ifndef _CMNCFG_H_INCLUDED_
#define _CMNCFG_H_INCLUDED_

#include "cfgfile.h"

class CCmnCfg
{
public:
	CCmnCfg();
	virtual ~CCmnCfg();
	static const char* GetParam(const char* cfgcode,BOOL_T bAllowNull=FALSE);
protected :
	static CCfgFile* m_CfgFile;
};

#endif
