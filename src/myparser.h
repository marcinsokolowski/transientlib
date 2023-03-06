#ifndef _MYPARSER_H__
#define _MYPARSER_H__

#include "mystring.h"
#include "basedefines.h"
#include "mytypes.h"

class CMyStrTable;
class CFraction;
class CEnvVar;

class BASELIB_EI MyParser : public mystring
{
public:
	MyParser();
	MyParser(const char* szString);
	const char* GetNextItem(const char* upto=",");
	const char* GetNextItemNew(const char* upto=",");
	const char* GetNextLine( BOOL_T bAddEndLine=TRUE );
	const char* GetBetween( char b,char e);
	const char* GetItem(int pos);

	BOOL_T GetVarAndValue(mystring& varname,mystring& value, const char* sep="=");
	BOOL_T GetVarAndValue( CEnvVar& paramDesc, const char* sep="=" );
	BOOL_T GetVarAndValueFromLine( mystring& varname,mystring& value );		
	
	void SkipWhite();
	void SkipWhiteAndSeps(const char* seps=",");
	void SkipNonNumbers();
	const char* GetItemUpTo(const char* upto,BOOL_T bReadNull=FALSE);
	const char* GetItemUpToString(const char* upto,BOOL_T bReadNull=FALSE);
	const char* GetItemUpToWithNull(const char* upto);
	int GetItems( CMyStrTable& tab, const char* sep="," );
	int GetItemsNew( CMyStrTable& tab, const char* sep="," );
	void Reset();
	void trimleft();
	void operator++();
	void PlusPlus();
	const char* GetCurrent(){ return &( (*this)[m_pos] ); }

	static BOOL_T IsWhite(char c);
	static BOOL_T IsNumber(char c);
	static BOOL_T IsNumerical(const char* str);
	
	static BOOL_T GetFraction( const char* szFraction, CFraction& fract  );
	static BOOL_T ParseSaveArea( const char* szVal, int& x0, int& y0, int& x1, int& y1 );
protected:
	int m_pos;
	int m_linepos;
	mystring m_szItem;
	mystring m_szLine;
};

#endif
