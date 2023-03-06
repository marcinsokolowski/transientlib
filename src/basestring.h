#ifndef _BaseString_H__
#define _BaseString_H__

/* class BaseString */
#include "basedefines.h"
#include "mytypes.h"
#include <stdio.h>
#include <string.h>

int my_strlen(const char* str);

class BASELIB_EI BaseString {
protected:
	// member variables
	char* m_pchData;
	int m_Size;
	int m_Length;

	// member functions 
	void ReAllocBuffer(int size=DEFAULT_SIZE);
	void Assign(const char* pszString,int len);
	void MakeReplacement(const char* pOld,const char* old,
					     const char* nnew);
public:
	BaseString();
	~BaseString();
	BaseString(long number);
	BaseString(long size,long size2);
	BaseString(const char z);
	BaseString(const char* szString);
	BaseString(const BaseString& szString);
	BaseString operator +(const char* pszString);
	BaseString operator +(const BaseString& szString);
	BaseString& operator +=(const char* pszString);
	BaseString& operator +=(const BaseString& szString);
	BaseString& operator =(int number);
	BaseString& operator =(const char* pszString);
	BaseString& operator =(const BaseString& szString);
	BaseString& operator =(double);
	void AllocBuffer(int size=DEFAULT_SIZE);
	
	int operator ==(const BaseString& szString) const;
	int operator ==(const char* szString) const;
	int operator <(const BaseString& szString) const;
	int operator <(const char* szString) const;
	int operator >(const BaseString& szString) const;
	int operator >(const char* szString) const;

	BaseString& operator<<(const char* str);
	BaseString& operator<<(const BaseString& str);
	BaseString& operator<<(const int n);
	BaseString& operator<<(const long n);
	BaseString& operator<<(INTN_T n);
	BaseString& operator<<(const char z);
	BaseString& operator<<(LONGLONG_T);
	BaseString& operator<<(double);

	// same as << but clears buffer before << 
	BaseString& operator<<=(const char* str);


	char& operator[](int pos);
	const char& operator[](int pos) const;
	inline operator const char*(){	return m_pchData; }
	inline int length() const { return m_Length; } 
	inline const char* c_str() const { return m_pchData; }
	inline const char* c_str(){ return m_pchData; }
	inline char* str(){ return m_pchData; }
	inline void clear() { m_pchData[0]='\0'; }
	int replace(const char* old,const char* nnew,BOOL_T bAll);
	int replace(const char* old,const char* nnew);
	void replace(int pos,int len,const char* nnew);
	void Reset();
	long find(const char* txt,long pos=0);
	long find(char z,long pos=0);
	BaseString substr(long pos,long len=0);
	BOOL_T empty() const { return ((!m_pchData) || (!m_pchData[0]));}
};


#endif
