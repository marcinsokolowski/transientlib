#ifndef _MY_STRING_HH
#define _MY_STRING_HH

#include "mytypes.h"
#include "basedefines.h"
#include "basestring.h"


BASELIB_EI const char* mygetenv(const char* varname);

class BASELIB_EI mystring : public BaseString
{
public:
	mystring(){};
	virtual ~mystring();
	mystring(const char* str):BaseString(str){};
	mystring(int n);
	mystring(const mystring& str);
	mystring(int size,int size2);
	// mystring(double db);

	/*BOOL_T operator==(const mystring& right);
	BOOL_T operator>(const mystring& right);		
	BOOL_T operator<(const mystring& right);		*/

	static mystring double2string( double x, const char* fmt="%.2f" );

	static mystring get_number( mystring& szStr );
	static int get_item( const char* pLine, int start_pos, int end_pos, mystring& szItem );
	static int get_item( const char* pLine, int start_pos, int end_pos, mystring& szItem, 
								int& end_out, int& last_space );
	static int get_item_safe( const char* pLine, int start_pos, int end_pos, mystring& szItem );
	static const char* skip_white( const char* line );

	mystring getupper() const;
	mystring getlower() const;
	mystring fromgreat();
	mystring& replace_char(char old,char nnew);
	mystring& extend(int min_length, char fill_char=' ');
	static void replace_char( char* ptr, char old,char nnew);
	
	
	void splitpath(BaseString& drv,BaseString& dir,
		 	       BaseString& fname,BaseString& ext) const;

	mystring& SkipApostrophs();
	mystring& SkipChar( char z );

	inline int Strcmp( const char* str ){ return strcmp( m_pchData, str ); }

	const char* env2str();
	
	static char get_first_non_white( const char* szString );
	static char get_last_non_white( const char* szString );

	const char* Fgets();	
	
	mystring TrimApostro( char z='\'' );
	mystring&  TrimApostrophs( char z='\'' );
	static BOOL_T IsInApostrophs( const char* szString );

	void SkipTrailSpaces();
	static void SkipTrailSpaces( char* ptr );

	// saving string to file :
	void DumpToFile( const char* szFileName, const char* mode="a+" );
};



#endif
