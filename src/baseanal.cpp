#include "baseanal.h"
#include <mathconst.h>
#include "tab2D.h"

CBaseAnal::CBaseAnal( int sizeX, int sizeY )
: m_SizeX(sizeX), m_SizeY(sizeY)
{

}



const char* CBaseAnal::GetLaplaceName( eLaplaceType_T laplaceType )
{
	return Table2D<ELEM_TYPE>::GetLaplaceName( laplaceType );
}
