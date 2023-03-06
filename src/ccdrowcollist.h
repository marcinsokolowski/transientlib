#ifndef _CCD_ROW_COL_LIST_H__
#define _CCD_ROW_COL_LIST_H__

#include "mytypes.h"

#include <vector>
using namespace std;

enum eCCDRowColType_T { eCCDUnknownType=0, eCCDRow, eCCDCol };

class CRowColDesc
{
public :
	CRowColDesc() : num(0),type(eCCDRow) {}
	CRowColDesc( int _num, eCCDRowColType_T _type ) :
		num(_num),type(_type) {}
	CRowColDesc( const CRowColDesc& right ){
		num = right.num;
		type = right.type;
	}		
	static eCCDRowColType_T GetType( const char* type_desc );
	
	int num;
	eCCDRowColType_T type;		
};

// list of row and columns for some specific action 
// can be read from file of syntax :
// COLUMN 1,23,54,900
// ROW    0,12,2047

#define ROW_INDICATOR "ROW"
#define COL_INDICATOR "COLUMN"

class CRowColList : public vector<CRowColDesc>
{
public:
	CRowColList(){};
	CRowColList( const char* filename );

	CRowColDesc* Find( eCCDRowColType_T type, int num );
	void AddUnique( CRowColDesc& elem );
	BOOL_T ReadFromFile( const char* filename );
};


class CCDWindow
{
public:
	int m_LowX;
	int m_LowY;
	int m_UpX;
	int m_UpY;

	CCDWindow() : m_LowX(0),m_LowY(0), m_UpX(0), m_UpY(0) {}
	CCDWindow( int low_x, int low_y, int up_x, int up_y )
		: m_LowX(low_x),m_LowY(low_y), m_UpX(up_x), m_UpY(up_y)		
	{}			
};


class CWindowList : public vector<CCDWindow>
{
public:
	CWindowList(){};
	CWindowList( const char* filename );
	void Add( int low_x, int low_y, int up_x, int up_y );

	BOOL_T ReadFromFile( const char* filename );	
};

// Class for storing CCD defects or other artifacts which should be flagged and excluded from analysis :
enum eCCDDefectType  { eCCDDefectSinglePixel=0, eCCDDefectRectangle=1, eCCDDefectCircle=2 };
// 
class CCDDefect
{
public :
	// hotpixel or other bad pixel :
	CCDDefect(double _x, double _y )
		: x(_x),y(_y),dx(0),dy(0),m_eDefectType(eCCDDefectSinglePixel)
	{
	}
	
	// rectangle ( ice crystals or other bad parts of CCD )
	CCDDefect(double _x, double _y, double _dx, double _dy )
		:  x(_x),y(_y),dx(_dx),dy(_dy),m_eDefectType(eCCDDefectRectangle)
	{
	}	

	// circle like object on CCD - rather for future of for ICE CRYSTALS as well 
	CCDDefect(double _x, double _y, double _radius )
		: x(_x),y(_y),dx(_radius),dy(_radius),radius(_radius)
	{
	}

	BOOL_T IsOK( double test_x, double test_y  );
	BOOL_T IsBAD( double test_x, double test_y );
	
	BOOL_T IsBadPixel( double test_x, double test_y, double check_radius=1.00  );
	
	void Dump();
	
	double x; // if rectangle - X of lower-left corner of defect 
	double y; // if rectangle - Y of lower-left corner of defect
	double dx; // for rectangle X size 
	double dy; // for rectangle Y size 
	double radius; // for cricle like defect 
	eCCDDefectType m_eDefectType;	   	
};


class CCDDefectList : public vector<CCDDefect>
{   
public:
	CCDDefectList(){};
	int ReadFromFile(const char* filename, int shift_dx=0, int shift_dy=0 );
	
	BOOL_T IsOK( double test_x, double test_y  );
	BOOL_T IsBAD( double test_x, double test_y );	   

	BOOL_T IsBadPixel( double test_x, double test_y, double check_radius=1.00  );	
};

#endif

