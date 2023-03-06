#ifndef _MANY_TAB_2D_H__
#define _MANY_TAB_2D_H__

#include "tab2D.h"

template<class ARG_TYPE>
class CManyTab2D {
protected:
	Table2D<ARG_TYPE>* m_pCCD;
	long m_SizeX;
	long m_SizeY;
	long m_Count;
	double m_Average;
public :
	CManyTab2D();
	~CManyTab2D();
	CManyTab2D( long xSize, long ySize, long count);
	CManyTab2D(const CManyTab2D<ARG_TYPE>& right);
	BOOL_T Init( long xSize, long ySize, long nCCD );
	void InitFrame();
	BOOL_T ClearState();
	void SetData(ARG_TYPE val);
	
	long GetCount() const { return m_Count; }
	Table2D<ARG_TYPE>& operator[](const int pos);
	Table2D<ARG_TYPE>* GetMatrix(const int pos);
	CManyTab2D<ARG_TYPE>& operator=(const CManyTab2D<ARG_TYPE>& right);
};


#endif
