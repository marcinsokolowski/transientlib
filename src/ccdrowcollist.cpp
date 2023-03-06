#include "ccdrowcollist.h"
#include "myfile.h"
#include "mystrtable.h"
#include "myparser.h"
#include <math.h>


eCCDRowColType_T CRowColDesc::GetType( const char* type_desc )
{
	if(strncmp(type_desc,ROW_INDICATOR,strlen(ROW_INDICATOR))==0){
		return eCCDRow;
	}
	if(strncmp(type_desc,COL_INDICATOR,strlen(COL_INDICATOR))==0){
		return eCCDCol;
	}
	return eCCDUnknownType;
}

CRowColList::CRowColList( const char* filename )
{
	ReadFromFile( filename );
}

BOOL_T CRowColList::ReadFromFile( const char* filename )
{
	if(!MyFile::DoesFileExist( filename )){
		return FALSE;
	}

	MyIFile in(filename);
	const char* line=NULL;
	CMyStrTable items;	

	clear();	
	while(line = in.GetLine(TRUE)){
		MyParser pars = line;
		if( mystring::get_first_non_white( pars.c_str() )=='#' )
			continue;
		pars.GetItems( items, "= " ); // white space separated 
		if(items.size()>=2){
			CMyStrTable list;
			MyParser pars2=items[1].c_str();
			pars2.GetItems( list, "," ); // coma separated 
			eCCDRowColType_T type = CRowColDesc::GetType( items[0].c_str() );
			if(type!=eCCDUnknownType){
				for(int i=0;i<list.size();i++){
					if(atol(list[i].c_str())>=0){
						CRowColDesc tmp;
						tmp.type = type;
						tmp.num = atol( list[i].c_str() );
						AddUnique( tmp );
					}
				}
			}
		}
	}

	return TRUE;
}


CRowColDesc* CRowColList::Find( eCCDRowColType_T type, int num )
{
	CRowColList::iterator i;
	for(i=begin();i!=end();i++){
		if( i->type==type && i->num==num ){
			return (&(*i));
		}
	}
	return NULL;
}

void CRowColList::AddUnique( CRowColDesc& elem )
{
	if(!Find( elem.type, elem.num )){
		push_back( elem );
	}
}


CWindowList::CWindowList( const char* filename )
{
	ReadFromFile( filename );	
}
                                                                                
BOOL_T CWindowList::ReadFromFile( const char* filename )
{
	if(!MyFile::DoesFileExist( filename )){
		return FALSE;
	}

	MyIFile in(filename);
	const char* line=NULL;
	CMyStrTable items;	

	clear();	
	while(line = in.GetLine(TRUE)){
		MyParser pars = line;
		if( mystring::get_first_non_white( pars.c_str() )=='#' )
			continue;
		pars.GetItems( items ); // white space separated 
		if(items.size()>=4){
			CCDWindow tmp;
			tmp.m_LowX = atol( items[0].c_str() );			
			tmp.m_LowY = atol( items[1].c_str() );			
			tmp.m_UpX = atol( items[2].c_str() );			
			tmp.m_UpY = atol( items[3].c_str() );			
			push_back( tmp );
		}
	}
	return TRUE;	
}

void CWindowList::Add( int low_x, int low_y, int up_x, int up_y )
{	
	CCDWindow win( low_x, low_y, up_x, up_y );
	push_back( win );	
}

int CCDDefectList::ReadFromFile(const char* filename, int shift_dx, int shift_dy)
{
	if(!MyFile::DoesFileExist( filename )){
		return -1;
	}

	MyIFile in(filename);
	const char* line=NULL;
	CMyStrTable items;	
	int n_badarea=0;
	int n_hot=0;
	int n_circle=0;

	clear();	
	while(line = in.GetLine(TRUE)){
		MyParser pars = line;
		if( mystring::get_first_non_white( pars.c_str() )=='#' )
			continue;
		pars.GetItems( items ); // white space separated 
		
		double x0 = atof( items[0].c_str() ) + shift_dx;
		double y0 = atof( items[1].c_str() ) + shift_dy;
		
		if(items.size()>=4){
			// rectangle :
			CCDDefect tmp( x0, y0, atof( items[2].c_str() ), atof( items[3].c_str() ) );
			push_back( tmp );
		}else{
			if(items.size()>=3){
				// circle :
				CCDDefect tmp( x0, y0, atof( items[2].c_str() ) );
				push_back( tmp );				
			}else{
				if(items.size()>=2){
					// single pixel :
					CCDDefect tmp( x0, y0 );
					push_back( tmp );
				}
			}
		}
		
		if( items.size()>=5 ){
			// if 5th column found it contains information about defect type :
			// enum eCCDDefectType  { eCCDDefectSinglePixel=0, eCCDDefectRectangle=1, eCCDDefectCircle=2 };
			int size_1=(size()-1);
			((*this)[size_1]).m_eDefectType = (eCCDDefectType)atol( items[4].c_str() );
			
			if( ((*this)[size_1]).m_eDefectType == eCCDDefectRectangle ){
				n_badarea++;
			}
			if( ((*this)[size_1]).m_eDefectType == eCCDDefectSinglePixel ){
				n_hot++;
			}
			if( ((*this)[size_1]).m_eDefectType == eCCDDefectCircle ){
				n_circle++;
			}
			
		}		
		
	}

	printf("CCDDefectList::ReadFromFile Statistics : total=%d, badarea=%d, hotpixels=%d, circlelike=%d\n",size(),n_badarea,n_hot,n_circle);
	return size();		
}

BOOL_T CCDDefect::IsOK( double test_x, double test_y )
{
	switch( m_eDefectType )
	{
		case eCCDDefectSinglePixel :
			if( fabs(test_x-x)<=1 && fabs(test_y-y)<=1 ){
				return FALSE;
			}  
			break;
		case eCCDDefectRectangle :
			if( test_x>=x && test_x<=(x+dx) && test_y>=y && test_y<=(y+dy) ){
				return FALSE;
			}
			break;
		case eCCDDefectCircle :
			if( sqrt((test_x-x)*(test_x-x)+(test_y-y)*(test_y-y)) <= radius ){
				return FALSE;
			}
			break;
		default :
			return TRUE;						                
	}
	
	return TRUE;
}

void CCDDefect::Dump()
{
	switch( m_eDefectType )
	{
		case eCCDDefectSinglePixel :
			printf("(%d,%d) - bad/hot pixel\n",(int)x,(int)y);
			break;
		case eCCDDefectRectangle :
			printf("(%d,%d)-(%d,%d) - bad rectangle area\n",(int)x,(int)y,(int)(x+dx),(int)(y+dy));
			break;
		case eCCDDefectCircle :
			printf("(%d,%d) , radius = %d\n",(int)x,(int)y,(int)radius);
			break;
	}
}

BOOL_T CCDDefect::IsBAD( double test_x, double test_y )
{
	return (!IsOK(test_x,test_y));
}

BOOL_T CCDDefect::IsBadPixel( double test_x, double test_y, double check_radius  )
{
	if( m_eDefectType == eCCDDefectSinglePixel ){
		if( fabs(test_x-x)<=check_radius && fabs(test_y-y)<=check_radius ){
			return TRUE;
		}  		
	}
	
	return FALSE;
}

BOOL_T CCDDefectList::IsBAD( double test_x, double test_y )
{
	for( CCDDefectList::iterator it=begin();it!=end();it++){
		if( it->IsBAD(test_x,test_y) ){
			return TRUE;
		}
	}

	return FALSE;
}

BOOL_T CCDDefectList::IsOK( double test_x, double test_y )
{
	return (!IsBAD(test_x,test_y));
}

BOOL_T CCDDefectList::IsBadPixel( double test_x, double test_y, double check_radius )
{
	for( CCDDefectList::iterator it=begin();it!=end();it++){
		if( it->IsBadPixel(test_x,test_y,check_radius) ){
			return TRUE;
		}
	}
	
	return FALSE;	
}