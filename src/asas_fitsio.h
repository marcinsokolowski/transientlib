#ifndef _ASAS_FITSIO_H_
#define _ASAS_FITSIO_H_

#define FITS_INT	1
#define FITS_STRING	2
#define FITS_FLOAT	3
#define FITS_DOUBLE	4
#define FITS_LONG	5
#define FITS_SHORT	6
#define FITS_EXACT	7
#define FITS_ULONG	8
#define EMPTY_HDR {0,NULL,0,0.,0.,0,0,0,0,0,0,0,0,0,0,NULL}

struct FITS {
  int fd;
  char *head;
  int hlines;
  float bs,bz;
  int bp;
  int naxis,nx,ny,nz;
  int compressed,blocksize,divisor,left,more;
  char *obuf;
};

/* typedef struct FITS FITS; */

/*struct FITS_HDR {
  int bitpix;
  int naxis;
  int naxis1;
  int naxis2;
  int naxis3;
  double bscale;
  double bzero;
  char object[32];
  double alpha;
  double delta;
  double itime;
  ulong apos;
  ulong dpos;
  char origin[32];
  char site[32];
  char observer[32];
  double tellong;
  double tellat;
  double telalt;
  char instrument[32];
  char camera[32];
  double gain;
  double rnoise;
  char camoptic[32];
  char filter[32];
  double focus;
  double pixscale;
  char savearea[32];
  short abinn;
  short sbinn;
  double case_temp;
  double chip_temp;
  time_t time_ut;
  double JD;
  double HJD;
  double equinox;
  double epoch;
  char date_obs[32];
  char time_obs[32];
  double ha;
  double ST;
  double airmass;
  double zenith_d;
  int badline;
  int rotate;                   
  int nimage;                  
  char filename[32];
  char software[32];
  int tracking_pulse;
  double ha_corr;
  double dec_corr;
  double steps_per_sec;
  double steps_per_arcsec;
  char comment[64];
};*/

/* compression */
int fits_rcomp(int a[],int nx,unsigned char *c,int clen,int nblock,int bsize);

int open_fits();
int close_fits();
int create_fits();
int create_float_fits();
int create_ushort_fits();
int create_simple_fits();
int read_fits_header();
int read_fits_data();
int write_fits_header();
int write_fits_data( int fd,void* data,int nx,int ny,int bitpix,int flip );
void *convert_data( void* inp,int npixel,int bpi,float bscale,float bzero,int bpo );
void add_keyword();
void replace_keyword();
void write_keyword();
void delete_keyword();
void extend();
char *get_keyword();
int close_fitsc();

int open_fitsc();
int open_fitsc_keep();
int read_fitsc_header();
int read_fitsc_data();
int write_fitsc_header();
int write_fitsc_data(struct FITS * pt, void *data,int nx,int ny,int flip);
void add_keywordc();
void write_keywordc();
void replace_keywordc();
void delete_keywordc();
char *get_keywordc();
const char* get_keyword_value( const char* fname, const char* key );


#endif

