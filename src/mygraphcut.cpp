#include "mygraphcut.h"
#include "calcrot.h"
#include <stdio.h>


CCutLineDef::CCutLineDef( double _a, double _b, double _c, eRelation _rel )
: a(_a), b(_b), c(_c), rel(_rel), m_MinX(-10000000.00), m_MaxX(+10000000.00)
{

}

CCutLineDef::CCutLineDef( double _a, double _b, eRelation _rel )
: a(_a), b(-1.00), c(_b), rel(_rel), m_MinX(-10000000.00), m_MaxX(+10000000.00)
{

}

CCutLineDef::CCutLineDef( double minX, double maxX, double _a, double _b, eRelation _rel )
: a(_a), b(-1.00), c(_b), rel(_rel), m_MinX(minX), m_MaxX(maxX)
{	
}

BOOL_T CCutLineDef::IsOK( double x, double y )
{
	if( x>=m_MinX && x<=m_MaxX ){
		double border_y = (-a*x-c)/b;
		
		if( rel == eRelGreater ){
			if( y > border_y ){
				return TRUE;
			}			
		}else{
			if( y < border_y ){
				return TRUE;
			}						
		}

		return FALSE;
	}
	
	return TRUE;
}
   

CGraphCut::CGraphCut()
{

}


void CGraphCut::AddLine( double a, double b, eRelation rel )
{
	CCutLineDef newline( a, b , rel );
	push_back( newline );
}

void CGraphCut::AddLine( double a, double b, double c, eRelation rel )
{
	CCutLineDef newline( a, b , c, rel );
	push_back( newline );
}

void CGraphCut::AddVerticalLine( double x_const, eRelation rel )
{
	AddLine( 1.00, 0.00, -x_const, rel );
}

void CGraphCut::AddHorizontalLine( double y_const, eRelation rel )
{
	AddLine( 0.00, 1.00, -y_const, rel );
}

void CGraphCut::AddLineByPoints( double x0, double y0, double x1, double y1,
		  								   eRelation rel )
{
	double a,b,c;
	CMyCalcRot::CalcLine( x0,y0,x1,y1,a,b,c  );
	
//	printf("LINE : %f,%f,%f\n",a,b,c);

	AddLine( a, b, c, rel );	
}


BOOL_T CGraphCut::CheckPoint( double x, double y )
{
	vector<CCutLineDef>::iterator i;

	for(i=begin();i!=end();i++){
		double a = i->a;
		double b = i->b;
		double c = i->c;

		if(b==0){
			// x = CONST
			// in this case smaller means - x<CONST , greater x>CONST
			double const_x = (-c/a);
			if( i->rel == eRelSmaller ){
				if( x>const_x )
					return FALSE;
			}
			if( i->rel == eRelGreater ){
				if( x<const_x )
					return FALSE;
			}
		}else{
			double y_on_line = (-a*x-c)/b;
			if( i->rel == eRelSmaller ){
				if(y_on_line<y)
					return FALSE;
			}
			if( i->rel == eRelGreater ){
				if(y<y_on_line)
					return FALSE;
			}
		}			
	}
	return TRUE;
}


const char* CGraphCut::GetRelDesc( eRelation rel )
{
	if(rel==eRelSmaller)
		return "SMALLER";
	if(rel==eRelGreater)
		return "GREATER";
	
	return "UNKNOWN";
}

void CGraphCut::Dump()
{
	vector<CCutLineDef>::iterator i;

	int j=0;
   for(i=begin();i!=end();i++){
		printf("Line %d : %.5f*x + %.5f*y + %.5f = 0 , relation=%s\n",j,i->a,i->b,i->c,GetRelDesc(i->rel));
		j++;
	}
}
