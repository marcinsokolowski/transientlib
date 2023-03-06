#ifndef _MY_CYCLIC_TAB_H
#define _MY_CYCLIC_TAB_H


// the idea of re-puting sample is that I keep number N of 
// samples to be re-put in table , adding sample 7 overwirtes sample 0 
// - but GetFirst/GetNext functon is incorrect I suppose !

template<class ARG_TYPE>
class CCycleTab2D {
protected:
	int m_Size;
	int m_Count;
	
	
	ARG_TYPE* m_pTable;


	int m_Newest;		
	int m_Oldest;
	int m_CurrPos;
public:
	CCycleTab2D( int size );
	~CCycleTab2D();
	
	void Add( const ARG_TYPE& elem );


	ARG_TYPE* GetFirst();
	ARG_TYPE* GetNext();	
	
	int GetSize() { return m_Size; }
	int GetCount() { return m_Count; }

};






#endif
