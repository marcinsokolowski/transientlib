/***************************************************************************\
 **  transientlib : library for identification of short optical flashes 
 **						with the wide field cameras.
 **  This software was written by Marcin Sokolowski ( msok@fuw.edu.pl ) 
 **	it was a substantial part of analysis performed for PHD thesis 
 **  ( Investigation of astrophysical phenomena in short time scales with "Pi of the Sky" apparatus )
 **	it can be used under the terms of GNU General Public License.
 **	In case this software is used for scientific purposes and results are
 **	published, please refer to the PHD thesis submited to astro-ph :
 **
 **		http://arxiv.org/abs/0810.1179
 **
 ** Public distribution was started on 2008-10-31
 **
 ** 
 ** NOTE : some of the files (C files) were created by other developers and 
 **        they maybe distributed under different conditions.
 ** 

 ******************************************************************************
 ** This program is free software; you can redistribute it and/or modify it
 ** under the terms of the GNU General Public License as published by the
 ** Free Software Foundation; either version 2 of the License or any later
 ** version. 
 **
 ** This program is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY; without even the implied warranty of
 ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 ** General Public License for more details. 
 **
 *\**************************************************************************

*/           
#ifndef _CCD_PHOTOMETRY_H__
#define _CCD_PHOTOMETRY_H__


#include <stdio.h>
#include <vector>
#include <mytypes.h>
#include "ccd_defines.h"
#include "astfile.h"

using namespace std;

class CMyStrTable;

#define MAG_0 21

struct cStarDesc
{
	cStarDesc() : mag_ap0(0),mag_ap1(0),mag_ap2(0),mag_ap3(0),mag_ap4(0),
					  cluster_n_mag(0) {}

	float x;
	float y;
	float mag;     // from out photometry 
	float sky;
	float sky_per_pixel;
	float sharp;
	float shape;
	double ra;      // RA
	double dec;     // DEC
	float mag_cat; // magnitudo from catalog 
	float lap_value; // laplace value
	float mag_lap;   // mag from laplace
	float cluster_lap_value; // cluster laplace value 
	float mag_cluster;  // mag from cluster
	float apert_lap_value;

	float cluster_n_plus;
	float cluster_n;
	float cluster_n_mag;

	float cat_star_dist; // distance to catalog star in arcsec
	// char szStarName[20]; // can be removed for optiomization reasons !

	// used in picalib to show number of catalog stars nearby 
	int cat_star_number;
	
	// value of maximum pixel :	
	int max_pixel_value;
	int max_pixel_two;
	int max_pixel_three;
	int max_pixel_many;
	int max_pixel_5x5;
	
	// apertures to be written to mag file :
	float mag_ap0;
	float mag_ap1;
	float mag_ap2;
	float mag_ap3;
	float mag_ap4;
};

struct cMagSlice
{
/*	cMagSlice() 
		: min_mag(0),max_mag(0),mag_diff(0),mag_count(0), min_cat_mag(1000),
		  max_cat_mag(-1)
	{}			*/
/*	cMagSlice( double _min_mag, double _max_mag, double _min_cat_mag, 
				  double _max_cat_mag, double _mag_diff, double _mag_count )
	{
	}				*/
		 
	double min_mag;
	double max_mag;
	double mag_aver;
	double min_cat_mag;
	double max_cat_mag;		
	double mag_diff;
	int mag_count;			
};

extern int gMagCorrSliceCount;
extern cMagSlice gMagCorrSlices[];

class CSafeKeyTab;
class mystring;
class CCDMatrix;
class InfoTable2D;
class CPixelAnalyseOut;
class CPixelAnalyseIn;
class CPixelList;
class CCDFastPhotoCorrFile;

class CPiPhotometry
{
public:
	static int m_bVerb;
	static int m_DetailX;
	static int m_DetailY;
	static CCDFastPhotoCorrFile m_CorrTable;
	static int gMapSize;
	static BOOL_T m_bFitLine;
	static BOOL_T m_bSubtractBackgr;
	static BOOL_T gSaveMapFITS;
	static int m_TimeLimit;
	static double m_MaxBlackRatio;
	
	// option for testing new code :
	static BOOL_T m_bTestVersion;

	// cataloging parameters 
	static double m_MinNormalizeMag;
	static double m_MaxNormalizeMag;
	static BOOL_T m_bSaveCorrFile;
	static int    m_MinSmoothStars;
	static double m_MaxSmoothRadius;
	static BOOL_T m_bCalcMaxClusterRadius;
	
	
	// small parts cataloging :
	static BOOL_T m_bAddShifts;
	                                                                                
	                                                                                
	
	// testing parameters :
	static BOOL_T m_bFillMaxPixelInfo;
	
	CPiPhotometry();
	~CPiPhotometry();
	
	// determination of image shift :
	static int CalcImageShift( vector<cStarCat>& list1, vector<cStarCat>& list2,
								double& dx, double& dy, double& sigma_dx, double& sigma_dy,
								int star_count=20,
								double mag_min=-1000, double mag_max=1000,
								int min_x=-1, int max_x=10000,
								int min_y=-1, int max_y=10000,
								int min_dist=50, double delta_mag=0.5,
								BOOL_T bVerb=FALSE, int min_cluster_points=-1  );

	static BOOL_T GetFastPhotoList( CCDMatrix& frame, vector<cStarCat>& starList,									
								CPixelAnalyseOut& out,
								CPixelAnalyseIn& in, CPixelList& usedList,
								double tresh_min, double tresh_max,
								double tresh_cluster, double sigma_g,
								int border=20, 
								int x_start=-1, int y_start=-1, int x_end=-1, int y_end=-1 );


	static int ExecPhotometryAndSave( CCDMatrix& frame, const char* mag_name,
													 double tresh, double tresh_cluster,
													 int border, int n_pixels, double r_apert=-1,
													 int x_start=-1, int y_start=-1, int x_end=-1, 
													 int y_end=-1, 
													 eDriverReverseImage_T flip=eReverseImageNone );													  

	static int ExecPhotometryNoSave( CCDMatrix& frame, const char* mag_name,
													 double tresh, double tresh_cluster,
													 int border, int n_pixels, double r_apert=-1,
													 int x_start=-1, int y_start=-1, int x_end=-1, 
													 int y_end=-1, 
													 eDriverReverseImage_T flip=eReverseImageNone );													  

	static BOOL_T ExecPhotometryOpt( CCDMatrix& frame, vector<cStarDesc>& starList,
												  double tresh, double tresh_cluster,
												  int border, int n_pixels, double r_apert=-1,
												  int x_start=-1, int y_start=-1, int x_end=-1, 
												  int y_end=-1, BOOL_T bDoCorr=FALSE );

	static BOOL_T ExecPhotometryBase( CCDMatrix& frame, vector<cStarDesc>& starList,
												  InfoTable2D& info, CPixelAnalyseOut& out,
												  CPixelAnalyseIn& in, CPixelList& usedList,
												  double tresh, double tresh_cluster,
												  int border, int n_pixels, double r_apert=-1,
												  int x_start=-1, int y_start=-1, int x_end=-1, 
												  int y_end=-1, BOOL_T bDoCorr=FALSE  );


	static BOOL_T DumpMagFile( vector<cStarDesc>& starList , CSafeKeyTab& keytab,
										double tresh, int xSize, int ySize,
										const char* mag_file, BOOL_T bDoGzip=FALSE,
										int naperts=5, int flip=-1 );	

	static BOOL_T FillMaxPixelInfo( ELEM_TYPE* data, LONG_T* cluster, int cluster_cnt,
											  int xSize, cStarDesc& star );


	static double calc_mag_aperture( double xc, double yc, ELEM_TYPE** data,
	                       double r_apert, double& field_size );


	static BOOL_T SaveStarMeasure( cStarDesc& tmp, CPixelAnalyseOut& out,
	                               CPixelAnalyseIn&  in, int lap_value,
	                               int raw_value, int tresh );

	static BOOL_T FindCircleCrossing_LineX( double const_x, 
												double x0, double y0, double r,
											   double& x_out, double& y_out,															
												int close_x, int close_y );

	static BOOL_T FindCircleCrossing_LineY( double const_y, 
												double x0, double y0, double r,
											   double& x_out, double& y_out,															
												int close_x, int close_y );


};


#endif

