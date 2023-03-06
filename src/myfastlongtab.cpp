#include <string.h>
#include "myfastlongtab.h"

CFastLongTable::CFastLongTable()
: counter(0)
{
	memset(values,0,sizeof(int)*MAX_LONG_TAB_SIZE);
}

