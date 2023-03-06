#ifndef _MY_INIFILE_H__
#define _MY_INIFILE_H__

#include "cfgfile.h"

class CIniFile : public CCfgFile
{

public :
	CIniFile( const char* initfile );
	~CIniFile();
	BOOL_T Save( BOOL_T bForceNew=TRUE );	
	
	mystring m_szRequstedIniFile;
};


#endif
