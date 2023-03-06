#ifndef _MY_GRAPH_CUT_H__
#define _MY_GRAPH_CUT_H__

#include <vector>
#include <math.h>
#include <mytypes.h>

using namespace std;

enum eRelation { eRelSmaller=0, eRelGreater };

class CCutLineDef
{
public :
	CCutLineDef( double _a, double _b, double _c, eRelation _rel=eRelGreater );
	CCutLineDef( double _a, double _b, eRelation _rel );	

	CCutLineDef( double minX, double maxX, double _a, double _b, eRelation _rel=eRelGreater );

	BOOL_T IsOK( double x, double y );

	double m_MinX;
	double m_MaxX;	
	
	double a;
	double b;
	double c;
	eRelation rel;
};

class CGraphCut : public vector<CCutLineDef>
{
public:
	CGraphCut();


	void AddLine( double a, double b, double c, eRelation rel );
	void AddLine( double a, double b, eRelation rel );
	void AddLineByPoints( double x0, double y0, double x1, double y1,
								 eRelation rel);
	void AddVerticalLine( double x_const, eRelation rel );
	void AddHorizontalLine( double y_const, eRelation rel );								 
								 
	BOOL_T CheckPoint( double x, double y );								 	
	
	static const char* GetRelDesc( eRelation rel );	
	void Dump();
};


#endif
