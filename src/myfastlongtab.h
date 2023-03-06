#ifndef _MY_FASTLONG_TAB_H__
#define _MY_FASTLONG_TAB_H__

#define MAX_LONG_TAB_SIZE 100

class CFastLongTable
{
public:
	int values[MAX_LONG_TAB_SIZE];
	int counter;
	
	CFastLongTable();
	inline void Add( int new_val ){ values[counter]=new_val;counter++; }
};




#endif
