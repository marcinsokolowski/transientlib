#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include "asas_errors.h"

#include "cat_new.h"
#include "cat_old.h"

int get_cat( double ra_min,double ra_max,double dec_min,double dec_max,
		       double** ra,double** dec,double** mag,int* ngsc,char* path,
             int field )
{
	if( !gAsasCatalogUseCache ){	
		return get_cat_old( ra_min, ra_max, dec_min, dec_max,
								  ra, dec, mag, ngsc, path, field );
	}else{
		return get_cat_new( ra_min, ra_max, dec_min, dec_max,
								  ra, dec, mag, ngsc, path, field );
	}
}            

