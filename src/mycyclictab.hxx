#include "mycyclictab.h"

template<class ARG_TYPE>
CCycleTab2D<ARG_TYPE>::CCycleTab2D( int size )
: m_Size(size),m_Count(0), m_Newest(-1), m_CurrPos(-1), m_Oldest(0),
  m_pTable(NULL)
{
	if( m_Size > 0 ){
		m_pTable = new ARG_TYPE[m_Size];
	}
}

template<class ARG_TYPE>
CCycleTab2D<ARG_TYPE>::~CCycleTab2D()
{
	if( m_pTable ){
		delete [] m_pTable;
	}
}

template<class ARG_TYPE>
void CCycleTab2D<ARG_TYPE>::Add( const ARG_TYPE& elem )
{
	m_Newest++;
	if( m_Count>=m_Size ){
		// inc oldest only in case cycle is full :
		m_Oldest++;
	}
	if(m_Newest>=m_Size)
		m_Newest = 0;
	if(m_Oldest>=m_Size)
		m_Oldest = 0;
	if( m_pTable && m_Newest<m_Size ){
		m_pTable[m_Newest] = elem;
	}
	if(m_Count<m_Size){
		m_Count++;
		if(m_Count>m_Size)
			m_Count=m_Size;
	}
}

template<class ARG_TYPE>
ARG_TYPE* CCycleTab2D<ARG_TYPE>::GetFirst()
{
	m_CurrPos=m_Oldest;
	if(m_CurrPos>=0 && m_Count>0 && m_pTable && m_CurrPos<m_Size)
		return &(m_pTable[m_CurrPos]);
	else
		return NULL;
}


template<class ARG_TYPE>
ARG_TYPE* CCycleTab2D<ARG_TYPE>::GetNext()
{
	if(!m_pTable || m_CurrPos==m_Newest){
		return NULL;
	}
	m_CurrPos++;
	if(m_CurrPos>=m_Size)
		m_CurrPos = 0;
	return &(m_pTable[m_CurrPos]);
}

