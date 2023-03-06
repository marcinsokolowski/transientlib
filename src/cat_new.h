#ifndef _CAT_NEW_H__
#define _CAT_NEW_H__

extern int gAsasCatalogUseCache;
extern struct CCatalogStar* gStarCatalogCache;

int get_cat_new( double ra_min,double ra_max,double dec_min,double dec_max,
             double** ra,double** dec,double** mag,int* ngsc,char* path,
             int field );

int initialize_starcat_cache( const char* cat_path );

int read_offset_cnt( const char* cat_path, int offset_tab[120][60], int cnt_tab[120][60] );

#endif
