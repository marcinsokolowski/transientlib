#ifndef _MYSAFEKEYTAB_H__
#define _MYSAFEKEYTAB_H__


class CEnvVar;
class CKeyTab;

class CSafeKeyTab 
{
protected:
	CKeyTab* m_Table;
public :
	CSafeKeyTab();
	CSafeKeyTab(const CSafeKeyTab& right );		
	~CSafeKeyTab();
	CSafeKeyTab& operator=(const CSafeKeyTab& right );			
	void Add( const char* key, const char* val, const char* comment="");
	void Add( const char* key, int val, const char* comment="" );
	void Add( const char* key, double val, const char* comment="" );		
	void Add( const CEnvVar& newelem );

	// removs key :
	void Delete( const char* key );
	
	// in case not found adds , changes otherwise :
	void Set( const char* key, int val, const char* comment="" );
	void Set( const char* key, const char* val, const char* comment="" );
	void Set( const char* key, double val, const char* comment="" );
	
	BOOL_T Get( const char* key, int& val );
	BOOL_T Get( const char* key, double& val );
	
	
	CEnvVar* Find(const char* key);	
	void Clear();
	long GetCount();
	CEnvVar& operator[](long i);
	
	CKeyTab& getKeyTab();		
	
	
	const char* getKeyVal( const char* key );
	
	// reading from file :
	BOOL_T ReadFromFile( const char* fname );
};


#endif
