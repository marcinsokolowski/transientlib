#include "myinifile.h"


CIniFile::CIniFile( const char* inifile ) : CCfgFile("") 
{
	if( MyFile::DoesFileExist(inifile) ){
		Init(inifile);
		m_szRequstedIniFile = inifile;
	}else{
/*		mystring szInHome = "$(HOME)/";
		szInHome << "." << inifile;*/

		// NEW 20050228 :
		mystring szInHome = inifile;
		szInHome.env2str();
		Init(szInHome.c_str());
		m_szRequstedIniFile = szInHome;
	}
}

CIniFile::~CIniFile()
{
}


BOOL_T CIniFile::Save( BOOL_T bForceNew )
{
	if(m_CfgTab.size() && ( strlen(m_CfgFile.GetFileName()) || bForceNew  ) ){
		mystring szFileName = m_szRequstedIniFile;
		if(strlen(m_CfgFile.GetFileName()))
			szFileName = m_CfgFile.GetFileName();

		MyOFile out( szFileName.c_str() );
		vector<CEnvVar>::iterator i;

	   for(i=m_CfgTab.begin();i!=m_CfgTab.end();i++){			
			out.Printf("%s=%s\n",i->szName.c_str(),i->szValue.c_str());
		}
		out.Flush();
		out.Close();			
		return TRUE;	
	}
	return FALSE;
}
