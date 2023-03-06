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
#include "ccd_photometry.h"
#include "ccd_matrix.h"
#include "ccd_globals.h"
#include <myprogress.h>
#include <mypixellist.h>
#include "ccd_analyse.h"
#include "ccd_asastransform.h"
#include "ccd_starcat.h"
#include <Astroangle.h>
#include <tab2Ddesc.h>
#include <mymacros.h>
#include "ccd_fastcalc.h"
#include <myfile.h>
#include "ccd_pipeline.h"
#include <fits_file.h>
#include "ccd_corr_file.h"
#include <mathfunc.h>
#include <AstroCCD.h>
#include "ccd_util.h"
#include <myfits.h>
#include <mystrtable.h>
#include <calcrot.h>

extern "C" {
int xy2ad_ast(const char* astfile,double x,double y,double* ra,double* dec);
int ad2xy_ast(const char* astfile,double ra,double dec,double* x,double* y);
#include "asas_fitsio.h"
}

#define SCALE   750.0     /* Atmospheric scale height */
#define STG0    24055.3   /* Greenwich ST (in sec) on 1970.01.01.00:00:00 GMT */

int CPiPhotometry::m_bVerb=0;
int CPiPhotometry::m_DetailX=-1;
int CPiPhotometry::m_DetailY=-1;
CCDFastPhotoCorrFile CPiPhotometry::m_CorrTable;
BOOL_T CPiPhotometry::m_bFillMaxPixelInfo=FALSE;
int CPiPhotometry::gMapSize=10;
BOOL_T CPiPhotometry::m_bCalcMaxClusterRadius=FALSE;
int CPiPhotometry::m_TimeLimit=300;
double CPiPhotometry::m_MaxBlackRatio=0.50;

// option for testing new code :
BOOL_T CPiPhotometry::m_bTestVersion=FALSE;

// fiting line 
BOOL_T CPiPhotometry::m_bFitLine=FALSE;

// cataloging :
double CPiPhotometry::m_MinNormalizeMag=5;
double CPiPhotometry::m_MaxNormalizeMag=10;
BOOL_T CPiPhotometry::m_bSaveCorrFile=FALSE;
int    CPiPhotometry::m_MinSmoothStars=30;
double CPiPhotometry::m_MaxSmoothRadius=5;
BOOL_T CPiPhotometry::m_bSubtractBackgr=FALSE;
BOOL_T CPiPhotometry::gSaveMapFITS=FALSE;

// small parts cataloging :
BOOL_T CPiPhotometry::m_bAddShifts = FALSE;


/* double  CPiPhotometry::m_MagCorrFitParams[5] = { -1.30690e-03 , 
3.76745e-02 , -3.89668e-01, 1.67886e+00, -2.38787e+00 };*/

// calculates shift ( list2.x - list1.x ) - > shift from list1 to list2 
int CPiPhotometry::CalcImageShift( vector<cStarCat>& list1, vector<cStarCat>& list2,
												double& dx, double& dy, double& sigma_dx, double& sigma_dy,
												int star_count,
												double mag_min, double mag_max,
												int min_x, int max_x,
												int min_y, int max_y,
												int min_dist, double delta_mag,
												BOOL_T bVerb, int min_cluster_points  )
{
	double prev_x=-1000,prev_y=-1000;
	int checked_stars=0;
	double sum_dx=0,sum_dy=0;
	int count=0;
	double sum2_dx=0,sum2_dy=0;

	vector<cStarCat>::iterator i;
	for(i=list1.begin();i!=list1.end();i++){
		if( checked_stars >= star_count ){
			break;
		}
		if( (*i).mag >= mag_min && (*i).mag<=mag_max && 
			 (*i).xc >= min_x && (*i).xc <= max_x && 
			 (*i).yc >= min_y && (*i).yc <= max_y && 
			 (fabs((*i).xc-prev_x)>=min_dist || fabs((*i).yc-prev_y)>=min_dist) && 
			 (*i).no_measurements>= min_cluster_points ){
			
			double min_dist=10000.00;
			cStarCat* pStar=NULL;
			// find closest star on second list :
			for(vector<cStarCat>::iterator j=list2.begin();j!=list2.end();j++){
				if( (*j).yc >= ( (*i).yc + 10 ) ){
					break;
				}
				if( fabs( (*j).yc - (*i).yc ) <= 10 && fabs((*j).mag-(*i).mag)<=delta_mag && 
					 (*j).no_measurements>= min_cluster_points ){
					double d = sqrt( ((*i).xc-(*j).xc)*((*i).xc-(*j).xc) + ((*i).yc-(*j).yc)*((*i).yc-(*j).yc) );
					if( d < min_dist ){
						min_dist = d;
						pStar = &(*j);
					}
				}
			}
			if( pStar && min_dist<10.00 ){
				if( bVerb ){
					printf("(%.2f,%2f,%d) -> (%.2f,%.2f,%d)  => (dx,dy) = ( %.2f , %.2f )\n",
						(*i).xc,(*i).yc,(*i).no_measurements,
						pStar->xc,pStar->yc,pStar->no_measurements,
						(pStar->xc-(*i).xc),
						(pStar->yc-(*i).yc));
				}
	
				double dx = (pStar->xc-(*i).xc);
				double dy = (pStar->yc-(*i).yc);
				sum_dx = sum_dx + dx;
				sum_dy = sum_dy + dy;
				sum2_dx = sum2_dx + (dx*dx);
				sum2_dy = sum2_dy + (dy*dy);
				count++;
			}
		}
	}

	dx = 0;
	dy = 0;
	printf("Number of matched stars = %d\n",count);
	if( count > 0 ){		
		dx = (sum_dx / count);
		dy = (sum_dy / count);

		sigma_dx = sqrt( (sum2_dx/double(count)) - dx*dx );
		sigma_dy = sqrt( (sum2_dy/double(count)) - dy*dy );
	}
	return count;
}

BOOL_T CPiPhotometry::GetFastPhotoList( CCDMatrix& frame, vector<cStarCat>& starList,
												  CPixelAnalyseOut& out,
												  CPixelAnalyseIn& in, CPixelList& usedList,
												  double tresh_min, double tresh_max,
												  double tresh_cluster, double sigma_g,
												  int border, 
											     int x_start, int y_start, int x_end, int y_end )
{
	int n_pixels=12;
	eLaplaceType_T laplace_type=gCCDParams.m_eLaplaceType;
	starList.clear();

	int xSize = frame.GetXSize();
	int ySize = frame.GetYSize();
	double mean,sigma,mean_raw,sigma_raw;

	double mean_sigma=0,mean_mean=0;
	int map_count=0;

	mean_sigma = ( mean_sigma / map_count );
	mean_mean  = ( mean_mean  / map_count );
	
	// int edge = gCCDParams.m_nIgnoreEdge;
	int edge = border;
	int low_y = edge+3;
	int low_x = edge+3;
	int up_x = ( xSize-edge-3 );
	int up_y = ( ySize-edge-3 );
	
	if( x_start>0 && y_start>0 ){
		low_x = x_start;
		low_y = y_start;
	}

	if( x_end>0 && y_end>0 ){
		up_x = x_end;
		up_y = y_end;
	}

	BIG_ELEM_TYPE** lap = (frame.get_laplace_data_fast());
	BIG_ELEM_TYPE* lap_data = (frame.m_pFrameLaplace)->get_data_buffer();
	CMyProgressBar bar( 0, ySize );

	CPixelAnalyseOut out2;
	in.p_data_fast = frame.get_data_buffer_fast();
	in.p_data = frame.get_data_buffer();
	in.p_curr_data_laplace = lap;
	in.xSize = xSize;
	in.ySize = ySize;

	int nAboveTreshold=0;
	double x0,y0,maxNoiseLevel;
	int max_pos;
	cStarCat tmp;	

	int nBreakWidth = 5; // [NEW] before was 1 
	int nRingWidth = 2;
	int nLastStarTresh=0;

	// temporary check :
	for(register int y=low_y;y<=up_y;y++){
		for(register int x=low_x;x<=up_x;x++){
			int pos = (y*xSize+x);
			if( usedList.CheckPixel( pos ) ){
				printf("CPiPhotometry::GetFastPhotoList : ERROR IN CODE !!!!!!!!!!!\n");	
				exit(-1);
			}
		}
	}

	time_t start_time = get_dttm();

	for(register int y=low_y;y<=up_y;y++){
		register int y_pos = y*xSize;
		for(register int x=low_x;x<=up_x;x++){
			int pos = (y_pos+x);

			if( !usedList.CheckPixel( pos ) ){
				int starTresh = tresh_min*sigma_g;
				int starTreshMax = tresh_max*sigma_g;
				int TreshCluster = tresh_cluster*sigma_g;
// TEST OF TRESHOLD SMOOTHING :
//				starTresh = tresh_map_ptr[y][x];
				nLastStarTresh = starTresh;

				if( lap[y][x] > starTresh && lap[y][x] < starTreshMax ){
					int val = lap[y][x];
					nAboveTreshold++;

					int dump_cluster=0;
					if( abs(x-1360)<=2 && abs(y-195)<=2 ){
						printf("odo");
						dump_cluster=1;
					}

					// Black pixels rejection :
					// black pixels are some kind of eletronical effect 
					// due to charge loss or something
				   Table2D<ELEM_TYPE>::GetLaplacePlusMinusValues( in.p_data_fast,
				      out2.m_PixelOut.laplacePlusValues, out2.m_PixelOut.laplaceMinusValues,
				      out2.m_PixelOut.laplacePlusCount, out2.m_PixelOut.laplaceMinusCount,
				      x,y,in.xSize,in.ySize, laplace_type );
					double black_ratio = CCD_Analyser::GetBlackRatio(
											 out2.m_PixelOut.laplacePlusValues, out2.m_PixelOut.laplaceMinusValues,
								          out2.m_PixelOut.laplacePlusCount, out2.m_PixelOut.laplaceMinusCount );
					if( black_ratio < 0.7 ){
							continue;
					}


					// now find cluster and choose max point :
					in.x = x;
					in.y = y;
					in.pos = pos;

					// to test LWP cluster add _LWP :
					CCD_Analyser::FindClusterAboveTresholdOpt3( in, out, 
																			  out.cluster, out.cluster_cnt,							
																			  x0,y0,
																			  maxNoiseLevel,
																			  TreshCluster, //clusterTreshold,
																			  // FALSE,
																			  FALSE /* , FALSE */ );
					CCD_Analyser::CalcCenterOfHitRealOpt( in.p_curr_data_laplace, out.cluster, out.cluster_cnt, in.xSize, x0, y0);
					
					int ret = CCD_Analyser::FindMaxInClusterBig( lap_data, xSize, out.cluster,
															 out.cluster_cnt, max_pos );

					// CPointList max_list;
					// int max_count = CCD_Analyser::CountClusterMAX( lap,  xSize, out.cluster,
					//										 out.cluster_cnt, max_list );


					// check if any of pixels in cluster already in used clusters:
					BOOL_T bSkipPixel=FALSE;
					for( register int cc=0;cc<out.cluster_cnt;cc++){
						if( usedList.CheckPixel( out.cluster[cc] ) ){
							bSkipPixel=TRUE;
							break;
						}
					}
					if( bSkipPixel )
						continue;



					// if(max_pos>=0 && !alreadyUsed.CheckPixel( max_pos ) ){
					if(max_pos>=0 && !usedList.CheckPixel( max_pos ) ){
						int xc = (max_pos%xSize);
						int yc = (max_pos/xSize);

						int lap_max_value = lap[yc][xc];


						// now calculate magnitudo :
						tmp.xc = xc;
						tmp.yc = yc;



						float plus_sum_out=0;
						double x0_npixel,y0_npixel;
						double cluster_n = 0;
						float sky=0;
						cluster_n = CCD_Analyser::CalcLapClusterOfMaxPoints( out2,
												xc, yc, xSize, ySize, in.p_data_fast, 
												in.p_data, 
												n_pixels , nBreakWidth, nRingWidth,
												sky, plus_sum_out, x0_npixel, y0_npixel,
												!CPiPhotometry::m_bSubtractBackgr );
						tmp.mag = 0;
						if( cluster_n > 0 ){
							tmp.mag = MAG_0-2.5*log10( cluster_n );
						}

						// x,y from CMS of cluster :
						tmp.xc = (float)x0;
						tmp.yc = (float)y0;

	
						// save number of pixels in cluster :
						tmp.no_measurements = out.cluster_cnt;
		

						// tmp.sky = 400;
						if( xc>low_x && yc>low_y && xc<up_x && yc<up_y && tmp.mag>0 ){
							BOOL_T bFound=FALSE;

							if( !bFound ){
								starList.push_back(tmp);
							}
						}
						for(register int ii=0;ii<out.cluster_cnt;ii++){
							usedList.HitPixel( out.cluster[ii] );
						}
						for(register int ii=0;ii<out2.cluster_cnt;ii++){
							usedList.HitPixel( out2.cluster[ii] );
						}
						usedList.HitPixel( max_pos );
					}
				}
			}
		}
		bar.SetValue( y );
		bar.Update();


		if( m_TimeLimit>0 ){
			time_t curr_time = get_dttm();						
			if( ( curr_time - start_time ) >= m_TimeLimit ){
				printf("Time limit for fast photometry exceeded (t=%d sec , limit=%d sec), stoped at y=%d\n",(curr_time - start_time),m_TimeLimit,y);
				break;
			}
		}
	}


	printf("\n%d pixels above treshold Tn=%d ADU \n",nAboveTreshold,nLastStarTresh);

	return TRUE;
}

int  CPiPhotometry::ExecPhotometryAndSave( CCDMatrix& frame, const char* mag_name,
													 double tresh, double tresh_cluster,
													 int border, int n_pixels, double r_apert,
													 int x_start, int y_start, int x_end, 
													 int y_end, eDriverReverseImage_T flip )
{
	vector<cStarDesc> starList;
	int xSize = frame.GetXSize();
   int ySize = frame.GetYSize();

	printf("Using tresholds %.2f / %.2f\n",tresh,tresh_cluster);

	BOOL_T bRet = ExecPhotometryOpt( frame, starList,
												tresh, tresh_cluster, border, n_pixels,
												r_apert, x_start, y_start, x_end, y_end );

	if( flip!=eReverseImageNone ){
		for(int i=0;i<starList.size();i++){
			if( flip==eReverseImageHor || flip==eReverseImageHorVert ){
				starList[i].x  = (xSize-starList[i].x);
			}
			if( flip==eReverseImageVert || flip==eReverseImageHorVert ){
				starList[i].y  = (ySize-starList[i].y);
			}
		}
	}


	if( bRet ){
		bRet = DumpMagFile( starList, frame.GetKeyTab(), tresh, xSize, ySize, mag_name );
	}

	if( !bRet )
		return 0;
	else
		return starList.size();
}

int  CPiPhotometry::ExecPhotometryNoSave( CCDMatrix& frame, const char* mag_name,
													 double tresh, double tresh_cluster,
													 int border, int n_pixels, double r_apert,
													 int x_start, int y_start, int x_end, 
													 int y_end, eDriverReverseImage_T flip )
{
	vector<cStarDesc> starList;
	int xSize = frame.GetXSize();
   int ySize = frame.GetYSize();

	printf("Using tresholds %.2f / %.2f\n",tresh,tresh_cluster);

	BOOL_T bRet = ExecPhotometryOpt( frame, starList,
												tresh, tresh_cluster, border, n_pixels,
												r_apert, x_start, y_start, x_end, y_end );

	if( !bRet )
		return 0;
	else
		return starList.size();
}

BOOL_T CPiPhotometry::ExecPhotometryOpt( CCDMatrix& frame, vector<cStarDesc>& starList,
												  double tresh, double tresh_cluster,
												  int border, int n_pixels, double r_apert,
												  int x_start, int y_start, int x_end, 
												  int y_end, BOOL_T bDoCorr  )
{
	if( bDoCorr ){
		if( !(CPiPhotometry::m_CorrTable).m_CorrTab ){
			if( !(CPiPhotometry::m_CorrTable).Read( gCCDParams.m_szFastPhotoCorrFile.c_str() ) ){
				printf("could not read correction file %s , exiting\n",gCCDParams.m_szFastPhotoCorrFile.c_str());
				exit(0);
			}
		}
	}

	int xSize = frame.GetXSize();
	int ySize = frame.GetYSize();

	// static removed !
	// due to daq crash on 20070507 night
	// in DAQ it is called in 2 threads - cannot use static variables
	// it should make no difference now :	
	InfoTable2D info( xSize, ySize, gMapSize, gMapSize, FALSE, border );
	CPixelList alreadyUsed( xSize*ySize ),usedList( xSize*ySize );
	CPixelAnalyseIn in;
	CPixelAnalyseOut out;
	alreadyUsed.Clear();
	usedList.Clear();
	
	in.pPixelList = &alreadyUsed;

	BOOL_T bRet=FALSE;

	printf("Running Normal version of photometry\n");
	bRet = ExecPhotometryBase( frame, starList, info, out, in, usedList, 
				 tresh, tresh_cluster, border, n_pixels, 
				 r_apert, x_start, y_start, x_end, y_end,
				 bDoCorr ); 
	return bRet;
}

BOOL_T CPiPhotometry::ExecPhotometryBase( CCDMatrix& frame, vector<cStarDesc>& starList,
												  InfoTable2D& info, CPixelAnalyseOut& out,
												  CPixelAnalyseIn& in, CPixelList& usedList,
												  double tresh, double tresh_cluster,
												  int border, int n_pixels, double r_apert,
											     int x_start, int y_start, int x_end, int y_end, 
												  BOOL_T bDoCorr )
{
	eLaplaceType_T laplace_type=eFivePlusFourMin;

	if( bDoCorr ){
		if( !(CPiPhotometry::m_CorrTable).m_CorrTab ){
			if( !(CPiPhotometry::m_CorrTable).Read( gCCDParams.m_szFastPhotoCorrFile.c_str() ) ){
				printf("could not read correction file %s , exiting\n",gCCDParams.m_szFastPhotoCorrFile.c_str());
				exit(0);
			}
		}
	}

	starList.clear();

	int xSize = frame.GetXSize();
	int ySize = frame.GetYSize();
	double mean,sigma,mean_raw,sigma_raw;

	double mean_sigma=0,mean_mean=0;
	int map_count=0;
	printf_now("calculating bacground map ...\n");

/*	if( CPiPhotometry::m_bSubtractBackgr ){
		for( int j=(gMapSize-1);j>=0;j--){
			for( int i=0;i<gMapSize;i++){
				Area2DInfo& elem = info.GetElem(i,j);

				BOOL_T bFit = frame.GetVariableMeanAndSigma( eRawS, mean_raw, sigma_raw, elem, 
										    (int)elem.m_LowLeft.x, (int)elem.m_LowLeft.y,
											 (int)elem.m_UpRight.x, (int)elem.m_UpRight.y,
											 j*gMapSize+i );
				if( !bFit ){
					printf("Fit of gauss FAILED !\n");
					return FALSE;
				}
				_TRACE_PRINTF_2("%d/%d ",(int)mean_raw,(int)sigma_raw);
			}
			_TRACE_PRINTF_2("\n");
		}
		my_printf_now("OK\n");	

		Table2D<float> map( frame.GetXSize(), frame.GetYSize() );
		if( SubtractBackground( frame, info, map, "bacgr" ) ){
			printf("bacground map subtracted OK\n");
		}	
	}*/

// changed on 20070202 - due to fact that laplace=12 is used on pi2 
// in flash recognition algorithm and it is slower in fast photometry
// maybe correction in cluster finding procedure is nessesary to improve it
// but I decided to always use laplace=4 ( g54 ) in fast photometry 
// for finding stars on image :
//	frame.Laplace();
	frame.Laplace( laplace_type );

	double avgStarTresh=0;
	int n_cell_count=0;
	printf("Background map size = %d x %d\n",gMapSize,gMapSize);
	for( int j=(gMapSize-1);j>=0;j--){
		for( int i=0;i<gMapSize;i++){
			Area2DInfo& elem = info.GetElem(i,j);

			BOOL_T bFit = frame.GetVariableMeanAndSigma( laplace_type , mean, sigma, elem, 
										    (int)elem.m_LowLeft.x, (int)elem.m_LowLeft.y,
											 (int)elem.m_UpRight.x, (int)elem.m_UpRight.y,
											 j*gMapSize+i, frame.get_laplace_data_fast(), 
											 TRUE /* use mean and rms when fit failed */ );
			if( !bFit ){
				printf("Fit of gauss FAILED , frame skiped !\n");
				return FALSE;
			}

			int starTresh = elem.m_DataInfo[laplace_type].m_Average
								+tresh*elem.m_DataInfo[laplace_type].m_Sigma;			
			if( starTresh<100 ){
				printf("WARNING (fit=%d) : treshold in (%d,%d)-(%d,%d) = %d ADU , it is too low - probably bad image (clouds, dome, or system error)!!!\n",
							bFit,
							(int)elem.m_LowLeft.x,(int)elem.m_LowLeft.y,
							(int)elem.m_UpRight.x,(int)elem.m_UpRight.y,starTresh);
				printf("Changing to VERY HIGH VALUE = 1000\n");
				elem.m_DataInfo[laplace_type].m_Average = 1000;
				elem.m_DataInfo[laplace_type].m_Sigma = 1000;
			}
			avgStarTresh = avgStarTresh + starTresh;
			n_cell_count++;


//			_TRACE_PRINTF_2("%d/%d ",(int)mean,(int)sigma);
//			printf("%d ",starTresh);

			mean_sigma = mean_sigma + elem.m_DataInfo[laplace_type].m_Sigma;
			mean_mean  = mean_mean  + elem.m_DataInfo[laplace_type].m_Average;
			map_count++;

		}
//		 _TRACE_PRINTF_2("\n");
//		printf("\n");fflush(stdout);
	}
	fflush(stdout);
	avgStarTresh = ( avgStarTresh / n_cell_count );
	printf("Average star treshold = %.2f\n",avgStarTresh);

/*	Table2D<float> tresh_map(xSize,ySize);
	SmoothTresholds( tresh_map, info, tresh, laplace_type );
	CCDUtil::WriteToFITSFile( tresh_map , "tresh_map.fit" );
	float** tresh_map_ptr = tresh_map.get_data_buffer_fast();*/


	mean_sigma = ( mean_sigma / map_count );
	mean_mean  = ( mean_mean  / map_count );
	
	// int edge = gCCDParams.m_nIgnoreEdge;
	int edge = border;
	int low_y = edge+3;
	int low_x = edge+3;
	int up_x = ( xSize-edge-3 );
	int up_y = ( ySize-edge-3 );
	
	if( x_start>0 && y_start>0 ){
		low_x = x_start;
		low_y = y_start;
	}

	if( x_end>0 && y_end>0 ){
		up_x = x_end;
		up_y = y_end;
	}

	BIG_ELEM_TYPE** lap = (frame.get_laplace_data_fast());
	BIG_ELEM_TYPE* lap_data = (frame.m_pFrameLaplace)->get_data_buffer();
	CMyProgressBar bar( 0, ySize );

	CPixelAnalyseOut out2;
	in.p_data_fast = frame.get_data_buffer_fast();
	in.p_data = frame.get_data_buffer();
	in.p_curr_data_laplace = lap;
	in.xSize = xSize;
	in.ySize = ySize;

	int nAboveTreshold=0;
	double x0,y0,maxNoiseLevel;
	int max_pos;
	cStarDesc tmp;	

	int nBreakWidth = 5; // [NEW] before was 1 
	int nRingWidth = 2;
	int nLastStarTresh=0;

	// temporary check :
	for(register int y=low_y;y<=up_y;y++){
		for(register int x=low_x;x<=up_x;x++){
			int pos = (y*xSize+x);
			if( usedList.CheckPixel( pos ) ){
				printf("CPiPhotometry::ExecPhotometryBase : ERROR IN CODE !!!!!!!!!!!\n");	
				exit(-1);
			}
		}
	}

	time_t start_time = get_dttm();

	for(register int y=low_y;y<=up_y;y++){
		register int y_pos = y*xSize;
		for(register int x=low_x;x<=up_x;x++){
			int pos = (y_pos+x);

//			if( abs( x-1312 )<5 && abs( y-1451 )<5 ){
//				printf("odo");
//			}

			if( !usedList.CheckPixel( pos ) ){
				Area2DInfo& elem = info.GetAreaDesc( x ,y );
				int starTresh = elem.m_DataInfo[laplace_type].m_Average
									+tresh*elem.m_DataInfo[laplace_type].m_Sigma;

// TEST OF TRESHOLD SMOOTHING :
//				starTresh = tresh_map_ptr[y][x];
				nLastStarTresh = starTresh;

//				if( abs( x-1312 )<5 && abs( y-1451 )<5 ){
//					printf("odo");
//				}

				if( lap[y][x] > starTresh ){
					int val = lap[y][x];
					nAboveTreshold++;

					// Black pixels rejection :
					// black pixels are some kind of eletronical effect 
					// due to charge loss or something
				   Table2D<ELEM_TYPE>::GetLaplacePlusMinusValues( in.p_data_fast,
				      out2.m_PixelOut.laplacePlusValues, out2.m_PixelOut.laplaceMinusValues,
				      out2.m_PixelOut.laplacePlusCount, out2.m_PixelOut.laplaceMinusCount,
				      x,y,in.xSize,in.ySize, laplace_type );
					double black_ratio = CCD_Analyser::GetBlackRatio(
											 out2.m_PixelOut.laplacePlusValues, out2.m_PixelOut.laplaceMinusValues,
								          out2.m_PixelOut.laplacePlusCount, out2.m_PixelOut.laplaceMinusCount );
					if( CPiPhotometry::m_MaxBlackRatio>0 ){
						if( black_ratio < CPiPhotometry::m_MaxBlackRatio )
							continue;
					}


					//int clusterTreshRaw = elem.m_DataInfo[eRawS].m_Average + 
					//						 elem.m_DataInfo[eRawS].m_Sigma*tresh_cluster;
					//double clusterTreshold = elem.m_DataInfo[laplace_type].m_Average+
					//							tresh_cluster*elem.m_DataInfo[laplace_type].m_Sigma;



					// now find cluster and choose max point :
					in.x = x;
					in.y = y;
					in.pos = pos;

					// to test LWP cluster add _LWP :
					CCD_Analyser::FindClusterAboveTresholdOpt3( in, out, 
																			  out.cluster, out.cluster_cnt,							
																			  x0,y0,
																			  maxNoiseLevel,
																			  starTresh, //clusterTreshold,
																			  // FALSE,
																			  TRUE /* , FALSE */ );
					
					int ret = CCD_Analyser::FindMaxInClusterBig( lap_data, xSize, out.cluster,
															 out.cluster_cnt, max_pos );

					// CPointList max_list;
					// int max_count = CCD_Analyser::CountClusterMAX( lap,  xSize, out.cluster,
					//										 out.cluster_cnt, max_list );


					// check if any of pixels in cluster already in used clusters:
					BOOL_T bSkipPixel=FALSE;
					for( register int cc=0;cc<out.cluster_cnt;cc++){
						if( usedList.CheckPixel( out.cluster[cc] ) ){
							bSkipPixel=TRUE;
							break;
						}
					}
					if( bSkipPixel )
						continue;


					// shuold be COMMENTED OUT - only for testing reasons !
					if( CPiPhotometry::m_bFillMaxPixelInfo ){
						FillMaxPixelInfo( in.p_data, out.cluster, out.cluster_cnt, xSize, tmp );
					}

					// temporary 11
					/*int max_raw_pos=-1;
					int max_raw_value=-10000;
					for(int pp=0;pp<out.cluster_cnt;pp++){
						if( in.p_data[ out.cluster[pp] ]>max_raw_value ){
							max_raw_value = in.p_data[ out.cluster[pp] ];
							max_raw_pos=out.cluster[pp];
						}
					}*/
					
				
					// if(max_pos>=0 && !alreadyUsed.CheckPixel( max_pos ) ){
					if(max_pos>=0 && !usedList.CheckPixel( max_pos ) ){
						int xc = (max_pos%xSize);
						int yc = (max_pos/xSize);



						// TEMPORARY TEST 11 :
						//int xc = ( max_raw_pos % xSize );
						//int yc = ( max_raw_pos / xSize );

						int lap_max_value = lap[yc][xc];

						/*double x_sum=0;
						double y_sum=0;
						double total_sum=0;
						for( int yyy=(yc-1);yyy<=(yc+1);yyy++){
							for( int xxx=(xc-1);xxx<=(xc+1);xxx++){
								x_sum += (in.p_data_fast[yyy][xxx])*xxx;
								y_sum += (in.p_data_fast[yyy][xxx])*yyy;
								total_sum += in.p_data_fast[yyy][xxx];		
							}
						}
						x0 = ( x_sum / total_sum );
						y0 = ( y_sum / total_sum );*/



						

						// now calculate magnitudo :
						tmp.x = xc;
						tmp.y = yc;

						// double l = (in_sum-(nInPix/nOutPix)*out_sum);
						// tmp.apert_lap_value = l;
						// tmp.mag =   MAG_0-2.5*log10(l);
						tmp.apert_lap_value = 0;
						tmp.mag = 0;

						// laplace of max pixel :
						tmp.lap_value = lap_max_value;
						tmp.mag_lap = MAG_0-2.5*log10( tmp.lap_value );

						// cluster :
						tmp.cluster_lap_value = 0;
						tmp.mag_cluster = 0;

						float plus_sum_out=0;
						double x0_npixel,y0_npixel;
						tmp.cluster_n = CCD_Analyser::CalcLapClusterOfMaxPoints( out2,
												xc, yc, xSize, ySize, in.p_data_fast, 
												in.p_data, 
												n_pixels , nBreakWidth, nRingWidth,
												tmp.sky, plus_sum_out, x0_npixel, y0_npixel,
												!CPiPhotometry::m_bSubtractBackgr );
						tmp.cluster_n_plus = plus_sum_out;
						if( tmp.cluster_n > 0 ){
							tmp.cluster_n_mag = MAG_0-2.5*log10( tmp.cluster_n );

							// temporary test - using max pixel : 11
							// tmp.x = (float)x0_npixel;
							// tmp.y = (float)y0_npixel;

							tmp.cluster_n_plus = plus_sum_out;
							// tmp.mag_lap = plus_sum_out;
						}else{
							_TRACE_PRINTF_2("ERROR at ( %d,%d ) cluster_n = %d\n",x,y,tmp.cluster_n);
							tmp.cluster_n_mag = 0;
						}
						tmp.mag_ap0 = tmp.cluster_n_mag;

						// x,y from CMS of cluster :
						tmp.x = (float)x0;
						tmp.y = (float)y0;
		
						tmp.sky_per_pixel = tmp.sky;

						if( r_apert>0 ){
							double field_size;
							tmp.cluster_n = calc_mag_aperture( tmp.x, tmp.y, in.p_data_fast,
																		  r_apert, field_size );
							tmp.cluster_n_mag = MAG_0-2.5*log10( tmp.cluster_n - field_size*tmp.sky );
							tmp.sky = tmp.sky*field_size;
							tmp.cluster_n_plus = tmp.cluster_n;
						}

						// filing other apertures :
						double tmp_value = CCD_Analyser::CalcSum( xc, yc, 1, tmp.sky, in.p_data_fast );
						tmp.mag_ap1 = MAG_0-2.5*log10( tmp_value );

						tmp_value = CCD_Analyser::CalcSumWithSide( xc, yc, 1, tmp.sky, in.p_data_fast );
						tmp.mag_ap2 = MAG_0-2.5*log10( tmp_value );

						/*tmp_value = CCD_Analyser::CalcSumWithSide( xc, yc, 1, tmp.sky, in.p_data_fast );
						double field_size;
						// tmp_value = calc_mag_aperture( tmp.x, tmp.y, in.p_data_fast, 2, field_size );
						// tmp.mag_ap3 = MAG_0-2.5*log10( tmp_value - field_size*tmp.sky );
						tmp.mag_ap3 = MAG_0-2.5*log10( tmp_value );

						tmp_value = CCD_Analyser::CalcSumWithSide( xc, yc, 2, tmp.sky, in.p_data_fast );
						// tmp_value = calc_mag_aperture( tmp.x, tmp.y, in.p_data_fast, 3, field_size );
						// tmp.mag_ap4 = MAG_0-2.5*log10( tmp_value - field_size*tmp.sky );
						tmp.mag_ap4 = MAG_0-2.5*log10( tmp_value );*/


						// if( bDoCorr && (fabs(tmp.x-100)<=2 && fabs(tmp.y-100)<=2) ){
							if( bDoCorr ){
								if( (CPiPhotometry::m_CorrTable).m_CorrTab ){
									float adu = ( tmp.cluster_n - tmp.sky );
									float adu_corr = (CPiPhotometry::m_CorrTable).Correct( tmp.x+0.5, tmp.y+0.5, adu );
									tmp.cluster_n = adu_corr;
									tmp.cluster_n_plus = tmp.cluster_n;
									tmp.cluster_n_mag = MAG_0-2.5*log10( adu_corr );
								}else{
									printf("ERROR : (CPiPhotometry::m_CorrTable).m_CorrTab=NULL !!!\n");
								}
							}
						// }
																									
	
						//if( bUseLaplace ){
						//tmp.mag = MAG_0-2.5*log10( in.p_data_fast[yc][xc] );
						//}else{
						//}
						tmp.sharp = 1.00;
						tmp.shape = 0.5;

//						if( CPiPhotometry::m_bCalcMaxClusterRadius ){
//							double max_radius = CCD_Analyser::CalcMaxClusterRadius( out.cluster, out.cluster_cnt, x0, y0, xSize );
							double max_radius = CCDFastCalc::FindMaxRedial( out.cluster, out.cluster_cnt, x0, y0, xSize );
							tmp.sharp = CCDFastCalc::CalcSphericity( max_radius, out.cluster_cnt );
//						}

						// tmp.sky = 400;
						if( xc>low_x && yc>low_y && xc<up_x && yc<up_y && tmp.cluster_n_mag>0 ){
							BOOL_T bFound=FALSE;

							/*int count_1=starList.size()-1;
							for(register int f=count_1;f>=0;f--){
								if( fabs( xc-starList[f].x )<=1 && fabs( yc-starList[f].y )<=1 ){
									bFound = TRUE;
									break;
								}
								if( y-starList[f].y > 10 ){
									break;
								}
							}*/

							if( !bFound ){
								if( CPiPhotometry::m_bAddShifts ){
									tmp.x = tmp.x + frame.m_X_On_Big;
									tmp.y = tmp.y + frame.m_Y_On_Big;
								}

/*								if( black_ratio < 0.8 ){
									printf("skiped (x,y)=(%d,%d)\n",x,y);
									printf("black_ratio = %.4f %.2f %.2f\n",black_ratio,tmp.x,tmp.y);
									continue;
								}else{
									printf("black_ratio = %.4f\n",black_ratio);
								}*/

//								if( fabs( tmp.x-1312 )<5 && fabs( tmp.y-1451 )<5 ){
//									printf("odo");
//								}

								starList.push_back(tmp);

// TESTING SPHERICITY :			
//	double max_radius = CCDFastCalc::FindMaxRedial( out.cluster, out.cluster_cnt, x0, y0, in.xSize );
// double spher = CCDFastCalc::CalcSphericity( max_radius, out.cluster_cnt );
//	printf("%.2f %.2f %.2f %d\n",x0,y0,spher,out.cluster_cnt);



								if( m_DetailX>0 && m_DetailY>0 ){
									if( fabs( xc-m_DetailX )<5 && fabs( yc-m_DetailY )<5 ){
										SaveStarMeasure( tmp, out2, in, lap[yc][xc], in.p_data[max_pos], starTresh );
									}
								}
							}
						}
						for(register int ii=0;ii<out.cluster_cnt;ii++){
							if( m_bVerb>=5 ){
								MyOFile out_file("clusters.txt","a+");
								int x = out.cluster[ii] % in.xSize;
								int y = out.cluster[ii] / in.ySize;
								out_file.Printf("%d %d %d\n",out.cluster[ii],x,y);
							}
							usedList.HitPixel( out.cluster[ii] );
						}
						for(register int ii=0;ii<out2.cluster_cnt;ii++){
							if( m_bVerb>=5 ){
                        MyOFile out_file("clusters.txt","a+");
                        int x = out.cluster[ii] % in.xSize;
                        int y = out.cluster[ii] / in.ySize;
                        out_file.Printf("%d %d %d\n",out.cluster[ii],x,y);
                     }
							usedList.HitPixel( out2.cluster[ii] );
						}
						usedList.HitPixel( max_pos );
					}
				}
			}
		}
		bar.SetValue( y );
		bar.Update();


		if( m_TimeLimit>0 ){
			time_t curr_time = get_dttm();						
			if( ( curr_time - start_time ) >= m_TimeLimit ){
				printf("Time limit for fast photometry exceeded (t=%d sec , limit=%d sec), stoped at y=%d\n",(curr_time - start_time),m_TimeLimit,y);
				break;
			}
		}
	}

	if( CPiPhotometry::m_bFitLine ){
		double* lap_tab = new double[starList.size()];
		double* mag_tab = new double[starList.size()];

		for(register int i=0;i<starList.size();i++){
			lap_tab[i] = MAG_0-2.5*log10( starList[i].lap_value );
			mag_tab[i] = starList[i].cluster_n_mag;
			// printf("FIT : %.8f %.8f\n",lap_tab[i],mag_tab[i]);
		}

		double a=0,b=0;
		double chi2 = CMyFit::FitLineChi2( lap_tab , mag_tab , starList.size(), a , b );

		// change on 2008-03-23 , uses average star treshold in background cells 
		// in order to have possibility to determine limit on the whole image
		// not only on image part with local background
		double tresh = avgStarTresh; // was nLastStarTresh		
		double lap_mag_value = MAG_0-2.5*log10( (double)tresh );
		double mag_limit = a*lap_mag_value + b;

		CSafeKeyTab& keys = frame.GetKeyTab();
		keys.Set( "SIGMAL", mean_sigma );
		keys.Set( "MEANL",  mean_mean );
		keys.Set( "TRESHL", tresh ); // change on 2008-03-23
		keys.Set( "MAGLIMIT", mag_limit );
		keys.Set( "FITA" ,  a );
		keys.Set( "FITB" ,  b );
		keys.Set( "FITCHI2", chi2 );
		

		char tmp_str[64];
		sprintf(tmp_str,"%d",starList.size() );
		keys.Set( "FITSTARS" , tmp_str );

		delete [] lap_tab;
		delete [] mag_tab;
	}

	printf("\n%d pixels above treshold Tn=%d ADU \n",nAboveTreshold,nLastStarTresh);

	return TRUE;
}

BOOL_T CPiPhotometry::DumpMagFile( vector<cStarDesc>& starList , CSafeKeyTab& keytab,
											  double tresh, int xSize, int ySize,
											  const char* mag_file,
											  BOOL_T bDoGzip, int naperts /*=5*/, int flip /*=-1*/ )
{
		int nfields = 5+naperts;
		int size = starList.size()*nfields;
		float* obuf = new float[ size ];
	
		for(int i=0;i<starList.size();i++){
			cStarDesc& star = starList[i];
			int pos=i*nfields;
			obuf[pos+0] = star.x+1;	
			obuf[pos+1] = star.y+1;
		
			// normal aperture :
			// obuf[pos+2] = starList[i].mag;

			// Bogumil - cluster :
			obuf[pos+2] = star.cluster_n_mag;
			

			obuf[pos+3] = star.sky;
			obuf[pos+4] = star.sharp;
			obuf[pos+5] = star.shape;
	
			// laplace :
			if( naperts >= 2 )
	         obuf[pos+6] = non_zero( star.mag_ap1 , star.cluster_n_mag );
			if( naperts >= 3 )
				obuf[pos+7] = non_zero( star.mag_ap2 , star.cluster_n_mag );
			if( naperts >= 4 )
				obuf[pos+8] = non_zero( star.mag_ap3 , star.cluster_n_mag );
			if( naperts >= 5 )
				obuf[pos+9] = non_zero( star.mag_ap4 , star.cluster_n_mag );
		}

		mystring szError;
		CMyFITSFile<float> out_fits;
		CSafeKeyTab keys = keytab;

		keys.Add( "PHOTOMTR", "APP_PHOT" );
		keys.Add( "N_APERT", naperts );
		keys.Add( "R_APERT", "2.00000000" );
		keys.Add( "R_APERT1", "2.00000000" );
		keys.Add( "FIELD_01" , "XC" );
		keys.Add( "FIELD_02"	, "YC" );
		keys.Add( "FIELD_03" , "MAG" );
		keys.Add( "FIELD_04" , "SKY" );
		keys.Add( "FIELD_05" , "SHARP" );
		keys.Add( "FIELD_06" , "SHAPE" );
		if( naperts>=2 ){
			keys.Add( "FIELD_07" , "MAG2" );
		}
		keys.Add( "KER_SIZE", 9 );
		keys.Add( "THRESHLD" , tresh );
		keys.Add( "FWHM" , 3.15979362 ); // TODO 
		keys.Add( "ELL" , 0.22567204 );  // TODO 
		keys.Add( "NAXIS3", 1 );
		keys.Set( FH_BZERO, "0.00000000e+00" );
		keys.Add( "FLAT_DIV" , "FLAT" );
		keys.Add( "DARK_SUB" , "DARK" );
		if( !keys.Find( "PIXSCALE" ) ){
			printf("Adding missing pixscale keyword = %.2f\n",gCCDParams.m_fPixScale );
			keys.Add( "PIXSCALE" , gCCDParams.m_fPixScale );
		}
		keys.Add( "NX" , xSize );
		keys.Add( "NY" , ySize );

		// setting flip according to what was done :
		if( flip>=0 ){
			keys.Set( FLIP , flip );
		}		

		keys.Delete( FH_BZERO );

		mystring szOutFile = mag_file;
		if(!out_fits.WriteToFile( NULL, nfields, starList.size(), 0, keys, 
										  szOutFile.c_str(), szError )){
			printf("error while writing mag file header\n");
		}

		FITS ftso;
		ftso.bp =-32;
		ftso.fd = out_fits.getFileDesc();
		int ret = write_fitsc_data( &ftso, obuf, nfields, starList.size(), 0 );
		printf("ret=%d\n",ret);
	
		delete [] obuf;
		out_fits.Close();
		
		if( bDoGzip ){
			mystring szCmd;
			szCmd << "gzip -f " << szOutFile;
			system( szCmd.c_str() );
		}

	return TRUE;
}

BOOL_T CPiPhotometry::FillMaxPixelInfo( ELEM_TYPE* data, LONG_T* cluster, int cluster_cnt,
													 int xSize, cStarDesc& star )
{
	LONG_T* tmp_cluster = new LONG_T[cluster_cnt];
	memcpy( tmp_cluster, cluster, cluster_cnt*sizeof(LONG_T) );
	
	for(int i=0;i<cluster_cnt;i++){
		int max_val=data[ tmp_cluster[i] ];
		int max_pos=-1;
		for( int j=(i+1);j<cluster_cnt;j++ ){
			if( data[ tmp_cluster[j] ] > max_val ){
				max_val = data[ tmp_cluster[j] ];
				max_pos = j;
			}
		}
		if( max_pos>i ){
			int tmp = tmp_cluster[i];
			tmp_cluster[i] = tmp_cluster[max_pos];
			tmp_cluster[max_pos] = tmp;
		}
	}

	if( cluster_cnt>0 )
		star.max_pixel_value = data[ tmp_cluster[0] ];	
	if( cluster_cnt>1 ){
		star.max_pixel_two = data[ tmp_cluster[1] ]+star.max_pixel_value;
	}
	if( cluster_cnt>2 ){
		star.max_pixel_three = data[ tmp_cluster[2] ]+star.max_pixel_two;
	}

	int pos = tmp_cluster[0];
	star.max_pixel_many = data[ pos ] + data[pos-1] + data[pos+1] +
								 data[ pos+xSize-1 ] + data[ pos+xSize ] + data[ pos+xSize+1 ] +
								 data[ pos-xSize-1 ] + data[ pos-xSize ] + data[ pos-xSize+1 ];

	star.max_pixel_5x5=0;
	for(int y=-2;y<=+2;y++){	
		for(int x=-2;x<=+2;x++){
			int p = pos+x+y*xSize;
			star.max_pixel_5x5 += data[p];
		}
	}


	/*for(int i=0;i<cluster_cnt;i++){
		printf("%d ",data[cluster[i]]);
	}
	printf("\n");
	for(int i=0;i<cluster_cnt;i++){
		printf("%d ",data[tmp_cluster[i]]);
	}
	printf("\n");*/
	for(int i=0;i<(cluster_cnt-1);i++){
		if( data[tmp_cluster[i]]<data[tmp_cluster[i+1]] ){
			printf("ERROR in code CPiPhotometry::FillMaxPixelInfo - sort not correct !!!\n");
			printf("%d < %d\n",data[tmp_cluster[i]],data[tmp_cluster[i+1]]);
			exit(0);
		}
	}
	

	delete [] tmp_cluster;		

	return TRUE;
}

double CPiPhotometry::calc_mag_aperture( double xc, double yc,
												  ELEM_TYPE** data,
												  double r_apert, double& field_size )
{
	int xc_int = (int)xc;
	int yc_int = (int)yc;

	// aperture limits :
	int y_low_ap = (int)(yc-r_apert);
	int y_up_ap = (int)(yc+r_apert+1);

	int x_low_ap = (int)(xc-r_apert);
	int x_up_ap = (int)(xc+r_apert+1);

	double x1,y1,x2,y2;
	double sum=0;	

	field_size = 0.00;
	for( int y=y_low_ap;y<=y_up_ap;y++){
		for( int x=x_low_ap;x<=x_up_ap;x++){	
			// if( x==735 && y==39 )
			//	printf("odo");
			
			int x_low=x;
			int x_up=(x+1);

			int y_low=y;	
			int y_up=(y+1);

			// calculating distances of pixel corners :
			double dist1 = CPoint::dist( x_low, y_low, xc, yc );
			double dist2 = CPoint::dist( x_up, y_low, xc, yc );
			double dist3 = CPoint::dist( x_low, y_up, xc, yc );
			double dist4 = CPoint::dist( x_up, y_up, xc, yc );


			// checkingif any of coners is closer to central pixel then requiered
			// aperture radius :
			if( dist1<=r_apert || dist2<=r_apert || dist3<=r_apert || dist4<=r_apert ){
				// at least one of pixel corners must be inside aperture :
				if( dist1<=r_apert && dist2<=r_apert && dist3<=r_apert && dist4<=r_apert ){
					// whole pixel inside 
					sum += data[y][x];
					field_size = field_size + 1.00;
				}else{
					// ONLY PART OF PIXEL INSIDE APERTURE :

					double fract=0.00;
					BOOL_T bHandled=FALSE;

					// 3 pixles inside , one outside  :
					if( ( dist1>r_apert && dist2<=r_apert && dist3<=r_apert && dist4<=r_apert ) || 
						 ( dist1<=r_apert && dist2>r_apert && dist3<=r_apert && dist4<=r_apert ) ||
						 ( dist1<=r_apert && dist2<=r_apert && dist3>r_apert && dist4<=r_apert ) ||
						 ( dist1<=r_apert && dist2<=r_apert && dist3<=r_apert && dist4>r_apert ) 
					  ){
						// ONE OF PIXELS INSIDE :
						if( dist1>r_apert ){
							if( FindCircleCrossing_LineX( x_low, xc, yc, r_apert, x1, y1, x_low, y_low ) && 
	 							FindCircleCrossing_LineY( y_low, xc, yc, r_apert, x2, y2, x_low, y_low ) ){
								// crossing found :
								fract = 0.5*(y1-y_low)*(x2-x_low);
							}
						}
						if( dist2>r_apert ){
							if( FindCircleCrossing_LineX( x_up, xc, yc, r_apert, x1, y1, x_low, y_low ) && 
	 							FindCircleCrossing_LineY( y_low, xc, yc, r_apert, x2, y2, x_low, y_low ) ){
								// crossing found :
								fract = 0.5*(y1-y_low)*(x_up-x2);
							}
						}
						if( dist3>r_apert ){
							if( FindCircleCrossing_LineX( x_low, xc, yc, r_apert, x1, y1, x_low, y_low ) && 
	 							FindCircleCrossing_LineY( y_up, xc, yc, r_apert, x2, y2, x_low, y_low ) ){
								// crossing found :
								fract = 0.5*(y_up-y1)*(x2-x_low);
							}
						}
						if( dist4>r_apert ){
							if( FindCircleCrossing_LineX( x_up, xc, yc, r_apert, x1, y1, x_low, y_low ) && 
	 							FindCircleCrossing_LineY( y_up, xc, yc, r_apert, x2, y2, x_low, y_low ) ){
								// crossing found :
								fract = 0.5*(y_up-y1)*(x_up-x2);
							}
						}

						bHandled=TRUE;						
						fract = 1.00 - fract;
					}
					
					if( !bHandled ){
						// One inside  , 3 outside  :
						if( ( dist1<=r_apert && dist2>r_apert && dist3>r_apert && dist4>r_apert ) || 
							 ( dist1>r_apert && dist2<=r_apert && dist3>r_apert && dist4>r_apert ) ||
							 ( dist1>r_apert && dist2>r_apert && dist3<=r_apert && dist4>r_apert ) ||
							 ( dist1>r_apert && dist2>r_apert && dist3>r_apert && dist4<=r_apert ) 
						  ){
							// ONE OF PIXELS INSIDE :
							if( dist1<=r_apert ){
								if( FindCircleCrossing_LineX( x_low, xc, yc, r_apert, x1, y1, x_low, y_low ) && 
	 								FindCircleCrossing_LineY( y_low, xc, yc, r_apert, x2, y2, x_low, y_low ) ){
									// crossing found :
									fract = 0.5*(y1-y_low)*(x2-x_low);
								}
							}
							if( dist2<=r_apert ){
								if( FindCircleCrossing_LineX( x_up, xc, yc, r_apert, x1, y1, x_low, y_low ) && 
	 								FindCircleCrossing_LineY( y_low, xc, yc, r_apert, x2, y2, x_low, y_low ) ){
									// crossing found :
									fract = 0.5*(y1-y_low)*(x_up-x2);
								}
							}
							if( dist3<=r_apert ){
								if( FindCircleCrossing_LineX( x_low, xc, yc, r_apert, x1, y1, x_low, y_low ) && 
	 								FindCircleCrossing_LineY( y_up, xc, yc, r_apert, x2, y2, x_low, y_low ) ){
									// crossing found :
									fract = 0.5*(y_up-y1)*(x2-x_low);
								}
							}
							if( dist4<=r_apert ){
								if( FindCircleCrossing_LineX( x_up, xc, yc, r_apert, x1, y1, x_low, y_low ) && 
	 								FindCircleCrossing_LineY( y_up, xc, yc, r_apert, x2, y2, x_low, y_low ) ){
									// crossing found :
									fract = 0.5*(y_up-y1)*(x_up-x2);
								}
							}							
							bHandled = TRUE;
						}
					}


					if( !bHandled ){
						if( ( dist1>r_apert && dist2>r_apert && dist3<=r_apert && dist4<=r_apert ) ){
							if( FindCircleCrossing_LineX( x_low, xc, yc, r_apert, x1, y1, x_low, y_low ) &&	
								 FindCircleCrossing_LineX( x_up, xc, yc, r_apert, x2, y2, x_low, y_low ) ){
								fract = 1.0*( y_up - (y1+y2)/2 );
							}
							bHandled = TRUE;
						}
						if( ( dist1<=r_apert && dist2<=r_apert && dist3>r_apert && dist4>r_apert ) ){
							if( FindCircleCrossing_LineX( x_low, xc, yc, r_apert, x1, y1, x_low, y_low ) &&	
								 FindCircleCrossing_LineX( x_up, xc, yc, r_apert, x2, y2, x_low, y_low ) ){
								fract = 1.0*( (y1+y2)/2  - y_low );
							}
							bHandled = TRUE;
						}

						if( ( dist1<=r_apert && dist2>r_apert && dist3<=r_apert && dist4>r_apert ) ){
							if( FindCircleCrossing_LineY( y_low, xc, yc, r_apert, x1, y1, x_low, y_low ) &&	
								 FindCircleCrossing_LineY( y_up, xc, yc, r_apert, x2, y2, x_low, y_low ) ){
								fract = 1.0*( (x1+x2)/2 - x_low );
							}
							bHandled = TRUE;
						}

						if( ( dist1>r_apert && dist2<=r_apert && dist3>r_apert && dist4<=r_apert ) ){
							if( FindCircleCrossing_LineY( y_low, xc, yc, r_apert, x1, y1, x_low, y_low ) &&	
								 FindCircleCrossing_LineY( y_up, xc, yc, r_apert, x2, y2, x_low, y_low ) ){
								fract = 1.0*( x_up - (x1+x2)/2 );
							}
							bHandled = TRUE;
						}
					}

					if( fract < 0.00 || fract>1.00 || !bHandled ){
						printf("ERROR in code !!!!!!!!!!!! fract = %.8f\n",fract);
						exit(0);
					}


					sum += fract*data[y][x];
					field_size = field_size + fract;
				}
			}					
		}
	}

	return sum;	
}

BOOL_T CPiPhotometry::SaveStarMeasure( cStarDesc& tmp, CPixelAnalyseOut& out,
												   CPixelAnalyseIn&  in, int lap_value,
													int raw_value, int tresh )
{
	mystring szFile;
	szFile << "star_" << m_DetailX << "_" << m_DetailY << ".txt";
	MyOFile out_file( szFile.c_str(), "a+" );
	mystring szLine;
	szLine << tmp.x << " " << tmp.y << " " << tmp.cluster_n_mag << " " 
			 << tmp.sky << " " << tmp.mag_lap << " " 
			 << lap_value << " "	<< raw_value << " "
			 << out.cluster_cnt << " " << tresh << " ";
	for( int i=0;i<out.cluster_cnt;i++){
		int x = (out.cluster[i] % in.xSize );
		int y = (out.cluster[i] / in.xSize );

		szLine << "(" << x << "," << y << ")=" << in.p_data[out.cluster[i]] << ",";
	}

	printf("star at (%d,%d) saved to file : %s\n",(int)tmp.x,(int)tmp.y,szFile.c_str());
	out_file.Printf("%s\n",szLine.c_str());
	
	return TRUE;
}

BOOL_T CPiPhotometry::FindCircleCrossing_LineX( double const_x, 
																double x0, double y0, double r,
															   double& x_out, double& y_out,															
																int close_x, int close_y )
{
	double d2 = r*r - (const_x-x0)*(const_x-x0);
	if( d2 >= 0){
		double sqrt_val = sqrt( d2 );	
		double y1 = y0 - sqrt_val;
		double y2 = y0 + sqrt_val;
		x_out = const_x;

		if( y1>=close_y && y1<=(close_y+1) ){
			y_out = y1;
			return TRUE;		
		}
		if( y2>=close_y && y2<=(close_y+1) ){
			y_out = y2;
			return TRUE;		
		}
	}
	return FALSE;
}							


BOOL_T CPiPhotometry::FindCircleCrossing_LineY( double const_y, 
																double x0, double y0, double r,
															   double& x_out, double& y_out,															
																int close_x, int close_y )
{
	double d2 = r*r - (const_y-y0)*(const_y-y0);
	if( d2 >= 0){
		double sqrt_val = sqrt( d2 );	
		double x1 = x0 - sqrt_val;
		double x2 = x0 + sqrt_val;
		y_out = const_y;
		if( x1>=close_x && x1<=(close_x+1) ){
			x_out = x1;
			return TRUE;		
		}
		if( x2>=close_x && x2<=(close_x+1) ){
			x_out = x2;
			return TRUE;		
		}
	}
	return FALSE;
}							


