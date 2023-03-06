#include <algorithm>
#include <functional>
#include "cfgfile.h"
#include "cexcp.h"
#include "myparser.h"
#include "mytrace.h"
#include "mycmnglobals.h"
#include "mymacros.h"

CCfgFile::CCfgFile(const char* filename)
: m_bInitialized(FALSE)
{
	Init( filename );
}

CCfgFile::CCfgFile( const CCfgFile& right )
{
	(*this) = right;
}


CCfgFile& CCfgFile::operator=( const CCfgFile& right )
{
	m_CfgTab.clear();
	for(int i=0;i<right.m_CfgTab.size();i++){
		m_CfgTab.push_back( right.m_CfgTab[i] );
	}
	m_bInitialized = TRUE;

	return (*this);
}

void CCfgFile::Init(const char* filename)
{
	if(filename && filename[0]){
		mystring szFileName;
		szFileName = filename;
		szFileName.env2str();
		if (MyFile::DoesFileExist(szFileName.c_str())){
			// szFileName = DEFAULT_CFG_FILE;
			m_CfgFile.Open(szFileName.c_str());
			InitCfgTab();
		}
	}

}

const char* CCfgFile::GetFileName()
{
	return m_CfgFile.GetFileName();
}

BOOL_T CCfgFile::FindParam(const char* cfgcode)
{
	vector<CEnvVar>::iterator i;
	for(i=m_CfgTab.begin();i!=m_CfgTab.end();i++){
		if (strcmp(i->szName.c_str(),cfgcode)==0){
			return TRUE;
		}
	}
	return FALSE;
}

const char* CCfgFile::GetParam(const char* cfgcode,BOOL_T bAllowNull)
{
	vector<CEnvVar>::iterator i;

	InitCfgTab();
	for(i=m_CfgTab.begin();i!=m_CfgTab.end();i++){
		if (strcmp(i->szName.c_str(),cfgcode)==0){
			MYTRACE3(gCmnTrace,"Requested parameter " << cfgcode << "=" << i->szValue.c_str());
			return i->szValue.c_str();
		}
	}
	if (!bAllowNull){
		Assert(FALSE,"Parameter %s not defined",cfgcode);
	}
	return NULL;
}

const char* CCfgFile::GetParamNoInit(const char* cfgcode,BOOL_T bAllowNull)
{
	vector<CEnvVar>::iterator i;

	for(i=m_CfgTab.begin();i!=m_CfgTab.end();i++){
		if (strcmp(i->szName.c_str(),cfgcode)==0){
			MYTRACE3(gCmnTrace,"Requested parameter " << cfgcode << "=" << i->szValue.c_str());
			return i->szValue.c_str();
		}
	}
	if (!bAllowNull){
		Assert(FALSE,"Parameter %s not defined",cfgcode);
	}
	return NULL;
}

BOOL_T CCfgFile::GetValue( const char* cfgcode, BOOL_T& boolVal )
{
	const char* szVal = GetParam( cfgcode , TRUE );
	if(szVal && szVal[0]){
		boolVal = (atol(szVal)>0);
		return TRUE;
	}
	return FALSE;
}

BOOL_T CCfgFile::GetValue( const char* cfgcode, LONG_T& longVal )
{
	const char* szVal = GetParam( cfgcode , TRUE );
	if(szVal && szVal[0]){
		longVal = atol(szVal);
		return TRUE;
	}
	return FALSE;
}


BOOL_T CCfgFile::GetValue( const char* cfgcode, double& doubleVal )
{
	const char* szVal = GetParam( cfgcode , TRUE );
	if(szVal && szVal[0]){
		doubleVal = atof(szVal);
		return TRUE;
	}
	return FALSE;
}

BOOL_T CCfgFile::GetValue( const char* cfgcode, mystring& stringVal )
{
	const char* szVal = GetParam( cfgcode , TRUE );
	if(szVal && szVal[0]){
		stringVal = szVal;
		return TRUE;
	}
	return FALSE;
}


void CCfgFile::SetValue( const char* cfgcode, BOOL_T boolVal )
{
	mystring szTmp;
	szTmp << boolVal;
	SetParam(cfgcode,szTmp.c_str());
}

void CCfgFile::SetValue( const char* cfgcode, LONG_T longVal )
{
	mystring szTmp;
	szTmp << longVal;
	SetParam(cfgcode,szTmp.c_str());
}

void CCfgFile::SetValue( const char* cfgcode, double doubleVal )
{
	char szTmp[100];
	sprintf(szTmp,"%.4f",doubleVal);
	SetParam(cfgcode,szTmp);
}
void CCfgFile::SetValue( const char* cfgcode, const char* stringVal )
{
	mystring szTmp(stringVal);
	SetParam(cfgcode,szTmp.c_str());
}

void CCfgFile::SetParam( const char* cfgcode, const char* cfgval )
{
	vector<CEnvVar>::iterator i;

	BOOL_T bFound=FALSE;
	for(i=m_CfgTab.begin();i!=m_CfgTab.end();i++){
		if (strcmp(i->szName.c_str(),cfgcode)==0){
			i->szValue = cfgval;
			bFound = TRUE;
			break;
		}
	}
	if(!bFound){
		CEnvVar tmp;
		tmp.szName = cfgcode;
		tmp.szValue = cfgval;
		m_CfgTab.push_back(tmp);
	}
}

BOOL_T CCfgFile::IsInitialized()
{
	return m_bInitialized;
}

BOOL_T CCfgFile::ReadCfgFile( MyIFile* pFile )
{
	int count=0;
	if (pFile->IsOpened()){
		_TRACE_PRINTF_0("Reading config file : %s ...\n",pFile->GetFileName());

		const char* pLine;
		mystring szParams;
		while(pLine = pFile->GetLine(TRUE)){
			if (strlen(pLine)==0 || pLine[0]=='\n')
				continue;
			// skip comments :
			if (pLine[0]=='#' || strncmp(pLine,"//",2)==0)
				continue;
			if(strncmp(pLine,"%LOAD%",6)==0 || strncmp(pLine,"%SAFELOAD%",10)==0 ){
				// loading include file :
				mystring szTmp;
				if( strncmp(pLine,"%LOAD%",6)==0 )
			      szTmp = pLine+7;
				else
					szTmp = pLine+11;
				mystring szFileName = szTmp.TrimApostro('"');
		      szFileName.env2str();
		      if (MyFile::DoesFileExist(szFileName.c_str())){
      		   // szFileName = DEFAULT_CFG_FILE;
					MyIFile includeFile;
		         includeFile.Open(szFileName.c_str());
		         ReadCfgFile( &includeFile );
	      	}else{
					if( strncmp( pLine,"%LOAD%", 6 )==0 ){
						Assert(FALSE,"Could not find include file : %s\n",szFileName.c_str());
					}else{
						if( gSafeLoadWarning ){
							_TRACE_PRINTF_0("could not load file : %s\n",szFileName.c_str());
							_TRACE_PRINTF_0("WARN : please verify if this can be accepted !!!\n");
							_TRACE_PRINTF_0("IGNORED command : %s\n",pLine);
						}
					}
				}
				continue;
			}

			MyParser parser = pLine;
			CEnvVar tmp;

			/*if( strstr(pLine,"dark") ){
				printf("odo");
			}*/

			if(!parser.GetVarAndValueFromLine(tmp.szName,tmp.szValue))
				continue;
			tmp.szValue.env2str();

			if(strlen(tmp.szName.c_str())==0 && strlen(tmp.szValue.c_str())==0)
				continue;
			
			// if(!FindParam(tmp.szName.c_str())){
				// adding only once same param - first in the file
				// takes over :
				m_CfgTab.push_back(tmp);
			// }
			szParams << pLine << "\n";			
			count++;
		}
		pFile->Close();
		MYTRACE2(gCmnTrace,"System parameters read :\n" << szParams);
		m_bInitialized = TRUE;
	}
	return count;
	
}

BOOL_T CCfgFile::InitCfgTab()
{
	if (m_CfgTab.size()==0 || m_CfgFile.IsOpened()){
		ReadCfgFile( &m_CfgFile );
	}	
	return m_CfgTab.size();
}

void CCfgFile::GetParamValues( vector<CEnvVar>& params )
{
	params.clear();
	vector<CEnvVar>::iterator i;

   for(i=m_CfgTab.begin();i!=m_CfgTab.end();i++){
		params.push_back( *i );
	}
}


void CCfgFile::GetParams( mystring& szParams )
{
	vector<CEnvVar>::iterator i;
	szParams = "";

	for(i=m_CfgTab.begin();i!=m_CfgTab.end();i++){
		szParams << i->szName << "=" << i->szValue << "\n";
	}
}

void CCfgFile::SaveToFile( const char* filename )
{
	vector<CEnvVar>::iterator i;
	mystring szMsg(2000,2000);
	for(i=m_CfgTab.begin();i!=m_CfgTab.end();i++){
		szMsg << i->szName.c_str() << "=" << i->szValue.c_str() << "\n";
	}
	mystring szName=filename;
	szName.env2str();
	MyOFile out( szName.c_str() );
	out.Printf("%s",szMsg.c_str());
}

void CCfgFile::Dump()
{
	vector<CEnvVar>::iterator i;
	mystring szMsg;
	for(i=m_CfgTab.begin();i!=m_CfgTab.end();i++){
		szMsg << i->szName.c_str() << "=" << i->szValue.c_str() << "\n";
	}
	printf("\n\nPARAMETER VALUES :\n%s",szMsg.c_str());
	printf("#############################\n");
}

void CCfgFile::Sort()
{
	if(m_CfgTab.size()){
		sort(m_CfgTab.begin(),m_CfgTab.end(),less<CEnvVar>());		
	}
}

BOOL_T CCfgFile::CompareEnvVarTables( vector<CEnvVar>& left, vector<CEnvVar>& right, 
												  mystring& szDifferent, mystring& szOnlyInLeft,
												  mystring& szOnlyInRight )
{	
	BOOL_T bRet = TRUE;

	szDifferent = "";
	szOnlyInLeft = "";
	szOnlyInRight = "";

	int l = 0;
	int r = 0;

	if(left.size()!=right.size() )
		bRet = FALSE;

	while( l<left.size() && r<right.size() ){
		if(left[l] == right[r]){
			l++;
			r++;
		}else{
			bRet = FALSE;
			if(left[l].szName==right[r].szName){
				szDifferent << left[l].szName << " : " << left[l].szValue 
								<< "     " << right[r].szValue << "\n";				
				l++;
				r++;
			}else{
				if(left[l].szName<right[r].szName){
					while(left[l].szName<right[r].szName && l<left.size()){
						szOnlyInLeft << left[l].szName << " : " << left[l].szValue << "\n";
						l++;
					}
				}else{
					if(left[l].szName>right[r].szName){
						while(left[l].szName>right[r].szName && r<right.size()){
							szOnlyInRight << right[r].szName << "=" << right[r].szValue << "\n";
							r++;
						}
					}
				}
			}
		}			
	}

	while(l<left.size()){
		szOnlyInLeft << left[l].szName << "=" << left[l].szValue << "\n";
		l++;
	}

	while(r<right.size()){
		szOnlyInRight << right[r].szName << "=" << right[r].szValue << "\n";
		r++;
	}


	return bRet;
}

BOOL_T CCfgFile::CompareEnvVarTablesNotSorted( vector<CEnvVar>& left, vector<CEnvVar>& right, 
												  mystring& szDifferent, mystring& szOnlyInLeft,
												  mystring& szOnlyInRight )
{	
	BOOL_T bRet = TRUE;

	szDifferent = "";
	szOnlyInLeft = "";
	szOnlyInRight = "";

	for(int i=0;i<left.size();i++){
		CEnvVar& var_left = left[i];
		CEnvVar* pVarRight = NULL;

		BOOL_T bRepeated=FALSE;
		for(int k=0;k<i;k++){
			if( strcmp( var_left.szName.c_str(), left[k].szName.c_str() ) == 0 ){
				bRepeated = TRUE;
				break;
			}
		}
		if( bRepeated ){
			continue;
		}
		
		for(int j=0;j<right.size();j++){
			if( strcmp( right[j].szName.c_str(), var_left.szName.c_str() ) == 0 ){
				// found :
				pVarRight = &(right[j]);
				break;
			}
		}
		
		if( pVarRight ){
			// found :
			if( strcmp( var_left.szValue.c_str() , pVarRight->szValue.c_str() ) ){
				szDifferent << var_left.szName << " : " << var_left.szValue 
				                      << "     " << pVarRight->szValue 
				                      << " ( line : " << (i+1) << " ) "
				                      << "\n";
				bRet = FALSE;				                      
			}		 	
		}else{
			szOnlyInLeft << var_left.szName << "=" << var_left.szValue << "\n";
			bRet = FALSE;
		}
	}

	for(int i=0;i<right.size();i++){
		CEnvVar& var_right = right[i];
		CEnvVar* pVarLeft = NULL;
		
		for(int j=0;j<left.size();j++){
			if( strcmp( left[j].szName.c_str(), var_right.szName.c_str() ) == 0 ){
				// found :
				pVarLeft = &(left[j]);
				break;
			}
		}
		
		if( !pVarLeft ){
			szOnlyInRight << var_right.szName << "=" << var_right.szValue << "\n";
			bRet = FALSE;
		}
	}



	return bRet;
}


void CCfgFile::CopyParamsTab( vector<CEnvVar>& dest, vector<CEnvVar>& src )
{
	dest.clear();
	vector<CEnvVar>::iterator i;
	for(i=src.begin();i!=src.end();i++){
		dest.push_back( *i );
	}
}


