#ifndef _ASAS_FITSIO_EXP_H__
#define _ASAS_FITSIO_EXP_H__

#define FITS_INT	1
#define FITS_STRING	2
#define FITS_FLOAT	3
#define FITS_DOUBLE	4
#define FITS_LONG	5
#define FITS_SHORT	6
#define FITS_EXACT	7
#define FITS_ULONG	8
#define EMPTY_HDR {0,NULL,0,0.,0.,0,0,0,0,0,0,0,0,0,0,NULL}

/* compression */
int fits_rcomp(int a[],int nx,unsigned char *c,int clen,int nblock,int bsize);
int fits_rdecomp( unsigned char *c,    /* input buffer         */
			         int clen,       /* length of input      */
                	unsigned int array[], /* output array         */
                  int nx,         /* number of output pixels */
                  int nblock,     /* coding block size    */
                  int bsize);
                                        

int write_fits_data( int fd,void* data,int nx,int ny,int bitpix,int flip );
int read_fits_data_new(int fd,int nx,int ny,int bitpix,void* data);

void *convert_data( void* inp,int npixel,int bpi,float bscale,float bzero,int bpo );
int close_fits(int fd,int nx,int ny,int bitpix);
const char* get_keyword_value( const char* fname, const char* key );


#endif

