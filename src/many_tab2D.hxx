#ifndef _MANY_TAB2D_HXX_
#define _MANY_TAB2D_HXX_

#include "many_tab2D.h"

template<class ARG_TYPE>
CManyTab2D<ARG_TYPE>::CManyTab2D()
: m_SizeX(0), m_SizeY(0), m_Count(0), m_pCCD(NULL)
{
}

template<class ARG_TYPE>
CManyTab2D<ARG_TYPE>::CManyTab2D(long xSize, long ySize, long count)
:m_SizeX(xSize),m_SizeY(ySize),m_Count(count),m_pCCD(NULL)
{
	Init( m_SizeX, m_SizeY, m_Count );
}

template<class ARG_TYPE>
CManyTab2D<ARG_TYPE>::CManyTab2D(const CManyTab2D<ARG_TYPE>& right)
:m_pCCD(NULL)
{
	(*this) = right;
}

template<class ARG_TYPE>
CManyTab2D<ARG_TYPE>::~CManyTab2D()
{
	delete [] m_pCCD;
}

template<class ARG_TYPE>
BOOL_T CManyTab2D<ARG_TYPE>::Init( long xSize, long ySize,long nCCD )
{
	m_SizeX = xSize;
	m_SizeY = ySize;
	m_Count = nCCD;
	if(m_pCCD)
		delete [] m_pCCD;
	if(nCCD){ 
		m_pCCD = new Table2D<ARG_TYPE>[nCCD];
		for(int i=0;i<nCCD;i++){
//   NEW change for gcc4.0 :
//			BOOL_T bOK = m_pCCD[i].Alloc(m_SizeX,m_SizeY);
			BOOL_T bOK = m_pCCD[i].InitConstructor( m_SizeX , m_SizeY );
		
			m_pCCD[i].SetIndex(i); 				
			if (!bOK){
				printf("Memory alloction problem occured, when allocating memory for CCD number %d\n",i);
				printf("Errno = %d\n",errno);
				printf("Size  = %dMB\n",(m_SizeX*m_SizeY)*sizeof(ARG_TYPE)/1000000);
				printf("Desc  = %s\n",strerror(errno));
				exit(-1);
			}
		}
	}	
	return TRUE;
}

template<class ARG_TYPE>
void CManyTab2D<ARG_TYPE>::InitFrame()
{
	SetData( 0 );
}

template<class ARG_TYPE>
CManyTab2D<ARG_TYPE>& CManyTab2D<ARG_TYPE>::operator=(const CManyTab2D<ARG_TYPE>& right)
{
	if(right.m_SizeX!=m_SizeX || right.m_SizeY!=m_SizeY || 
      right.m_Count!=m_Count)
		Init( right.m_SizeX, right.m_SizeY, right.m_Count );	
	for(int i=0;i<m_Count;i++){
		m_pCCD[i] = right.m_pCCD[i];
	}
	return (*this);
}

template<class ARG_TYPE>
Table2D<ARG_TYPE>* CManyTab2D<ARG_TYPE>::GetMatrix(const int pos)
{
	Assert(pos>=0 && pos<m_Count,"CCD aparatus has only %d cameras, you want to access camers number %d",m_Count,pos);
	return &(m_pCCD[pos]);
}


template<class ARG_TYPE>
Table2D<ARG_TYPE>& CManyTab2D<ARG_TYPE>::operator[](const int pos)
{
	Assert(pos>=0 && pos<m_Count,"CCD aparatus has only %d cameras, you want to access camers number %d",m_Count,pos);
	return m_pCCD[pos];
}


template<class ARG_TYPE>
void CManyTab2D<ARG_TYPE>::SetData(ARG_TYPE val)
{
   for(int i=0;i<m_Count;i++){
      m_pCCD[i].SetData(val);
   }
}

template<class ARG_TYPE>
BOOL_T CManyTab2D<ARG_TYPE>::ClearState()
{
	for(int i=0;i<m_Count;i++){
		m_pCCD[i].ClearState();
	}
   return TRUE;
}


#endif
