#ifndef _MYPIXELLIST_H_INCLUDED_H__
#define _MYPIXELLIST_H_INCLUDED_H__

#include "mytypes.h"

class CPixelList
{
protected :
	char* m_pList;
	LONGLONG_T m_nPixels;
	LONG_T m_Size;	
	
	// static BOOL_T m_bActive;
public :
	CPixelList(long nSize);	
	~CPixelList();

	// void Activate(){ m_bActive=TRUE; }
	// void Deactivate(){ m_bActive=FALSE; }

	void ClearPixel(long nPos);	
	BOOL_T CheckPixel(long nPos);
	void HitPixel(long nPos);		
	void Clear();
	
	void ClearPart( LONG_T* list_to_clear, LONG_T cnt );
};


#endif
