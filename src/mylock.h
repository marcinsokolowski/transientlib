#ifndef _MY_LOCK_H__
#define _MY_LOCK_H__

#include <pthread.h>
#include "mytypes.h"


class CMyMutex
{
public :
	CMyMutex();
	~CMyMutex();

	int Lock();
	int UnLock();
	int TryLock();
	int IsLocked(){ return is_locked; }
	
protected :	
	pthread_mutex_t mutex;	
	int id;
	int is_locked;
};


class CMutexLock
{
protected:
	CMyMutex* m_pLock;
public:
	CMutexLock( CMyMutex* pLock ): m_pLock(pLock) { 
		m_pLock->Lock(); 
	}
	~CMutexLock(){ m_pLock->UnLock(); }
};


class CMutexUnLock
{
protected:	
   CMyMutex* m_pLock;
   BOOL_T m_bDoUnLock;
public:
	CMutexUnLock( CMyMutex* pLock, BOOL_T bDoUnlock=TRUE ): m_pLock(pLock), m_bDoUnLock(bDoUnlock) {}
	~CMutexUnLock(){ 
		if( m_bDoUnLock && m_pLock ){
			m_pLock->UnLock(); 
		}			
	}
};

#endif
