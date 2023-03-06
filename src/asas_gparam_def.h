#ifndef _ASAS_GPARAM_DEF_H__
#define _ASAS_GPARAM_DEF_H__

#define ASAS_TRANSFORM_PARAM_COUNT 21

//
// WARNING !!!!!!!!!!!!!!!!!!!
// when adding new fields to this structure , they must be initialized in 
// ccd_asastransform.cpp file in constructor : CCDAsasTransform::CCDAsasTransform
// WARNING !!!!!!!!!!!!!!!!!!!
//
struct GPARAM{
  int order;            // number of parameters
  double px[21];
  double py[21];
  double pixscale;      // first order pixscale
  double fi;            // first order rotation
  double ra,dec,xc,yc;
  double ast_err;       // NEW - value of astrometry error 
  int flip;             // flip of image before transformation can be aplied , 0-none,2-horizontal,3-vertical
};


#endif

