#ifndef _MY_PIPE_H__
#define _MY_PIPE_H__

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#define RELEASE_TYPE_NORMAL       'Z'
#define RELEASE_TYPE_REPEAT_SAME  'R'

class CMyPipe
{
public:
	int fd[2];
	FILE* pipe_write;
	FILE* pipe_read;
	
	int m_Reader;
	int m_Writer;
	
	CMyPipe();
	~CMyPipe();
	
	int WaitForData( char& data );
	int WriteData( char c='Z' );
};


class CPipeLock : public CMyPipe
{
public:
	CPipeLock(){};
	~CPipeLock(){};
	
	char Wait();
	int Release( char c = RELEASE_TYPE_NORMAL );
};

#endif

