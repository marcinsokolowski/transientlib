#ifndef _PREDICT_H__
#define _PREDICT_H__

/* This file was generated by the installer program */

/* char *predictpath={"/opt/pi/ext/dload/sat_soft/PREDICT/predict-2.2.2/"}, soundcard=1, *version={"2.2.2"};*/
extern char *predictpath;
extern char soundcard;
extern char *version;

/* trace level */
extern int gSatlibTraceLevel;

#include "satlibdefs.h"


typedef struct {
	double lat, lon, alt, theta;
}  geodetic_t;
               

int CalcCurrentPosition( int x );
int CalcCurrentPositionOfSingle( int x, time_t curr_time, struct satInfo* sat_info );
void CheckAll( time_t ut_time );
int InitSatLib( const char* qth_file, const char* tle_file, const char* db_file );

// 
extern int gSatNumber;
extern geodetic_t obs_geodetic;


#endif
