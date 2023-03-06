#include "cmncfg.h"


CCfgFile* CCmnCfg::m_CfgFile=NULL;


CCmnCfg::CCmnCfg()
{}

CCmnCfg::~CCmnCfg()
{}


const char* CCmnCfg::GetParam(const char* cfgcode,BOOL_T bAllowNull)
{
	if(!m_CfgFile){
		m_CfgFile = new CCfgFile("$(CMNCFG)");
	}
	return m_CfgFile->GetParam(cfgcode,bAllowNull);
}

