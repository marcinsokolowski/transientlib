#include "mylock.h"
#include <stdio.h>
#include "mymacros.h"

CMyMutex::CMyMutex()
: is_locked(0)
{
	id = pthread_mutex_init( &mutex, NULL );
}

CMyMutex::~CMyMutex()
{
	UnLock();
	pthread_mutex_destroy( &mutex );
}

int CMyMutex::Lock()
{
	int ret = pthread_mutex_lock( &mutex );
	is_locked = 1;
	_TRACE_PRINTF_6("Mutex %ld : LOCKED\n",(long int)(this));
	
	return ret;
}


int CMyMutex::UnLock()
{
	int ret = pthread_mutex_unlock( &mutex );
	is_locked = 0;
	_TRACE_PRINTF_6("Mutex %ld : UN-LOCKED\n",(long int)(this));
	
	return ret;
}

int CMyMutex::TryLock()
{
   return pthread_mutex_trylock( &mutex );
}

