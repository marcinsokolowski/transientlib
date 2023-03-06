#ifndef _CFG_FILE_H__
#define _CFG_FILE_H__

#include "basestructs.h"
#include <vector>
#include "myfile.h"

using namespace std;

class BASELIB_EI CCfgFile 
{
public:
	CCfgFile(const char* filename);
	CCfgFile( const CCfgFile& right );
	CCfgFile& operator=( const CCfgFile& right );
	
	void Init(const char* filename);
	BOOL_T FindParam(const char* cfgcode);
	const char* GetParam(const char* cfgcode,BOOL_T bAllowNull=FALSE);
	const char* GetParamNoInit(const char* cfgcode,BOOL_T bAllowNull=FALSE);
	void SetParam( const char* cfgcode, const char* cfgval );
	void GetParamValues( vector<CEnvVar>& params );
	void GetParams( mystring& szParams );
	void Dump();
	void SaveToFile( const char* filename );
	vector<CEnvVar>& GetParamTable(){ return m_CfgTab; }		
	void Sort();
	static BOOL_T CompareEnvVarTables( vector<CEnvVar>& left, vector<CEnvVar>& right,
										 mystring& szDifferent, mystring& szOnlyInLeft,
                               mystring& szOnlyInRight );

	static BOOL_T CompareEnvVarTablesNotSorted( vector<CEnvVar>& left, vector<CEnvVar>& right,
										 mystring& szDifferent, mystring& szOnlyInLeft,
                               mystring& szOnlyInRight );


	BOOL_T ReadCfgFile( MyIFile* pFile );
	BOOL_T InitCfgTab();
	BOOL_T IsOpened(){ return m_CfgFile.IsOpened(); }
	BOOL_T IsInitialized();
	const char* GetFileName();

	BOOL_T GetValue( const char* cfgcode, BOOL_T& boolVal );
	BOOL_T GetValue( const char* cfgcode, LONG_T& longVal );		
	BOOL_T GetValue( const char* cfgcode, double& doubleVal );		
	BOOL_T GetValue( const char* cfgcode, mystring& stringVal );		

	void SetValue( const char* cfgcode, BOOL_T boolVal );
	void SetValue( const char* cfgcode, LONG_T longVal );		
	void SetValue( const char* cfgcode, double doubleVal );		
	void SetValue( const char* cfgcode, const char* stringVal );		
	
	static void CopyParamsTab( vector<CEnvVar>& dest, vector<CEnvVar>& src );
protected:
	vector<CEnvVar> m_CfgTab;
	MyIFile m_CfgFile;
	BOOL_T m_bInitialized;
};



#endif
