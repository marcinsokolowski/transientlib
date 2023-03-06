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
#include "myfits.h"
#include <string.h>
#include <stdio.h>
#include <cexcp.h>
#include "mathfunc.h"
// #include <datsrc.h>
// #include <grsrc.h>
#include <myutil.h>
#include <mylock.h>


void* CMyFit::root_func=NULL; 


#ifndef _NO_ROOT_
#include <TGraph.h>
#include <TF1.h>
#include <TF2.h>
#include <TRandom.h>

TF1* pGaussFunc=NULL;
CMyMutex gGlobalLock;
                                                                                
BOOL_T MyFits_FitGauss( double minValue, double& rms, double& center,
                      double &max, int binNo, double binWidth,
                      int* pCountTab, int& nstep, int first_bin )
{
	gGlobalLock.Lock();

	if( !pGaussFunc ){
		pGaussFunc = new TF1("fit_func","gaus",-100000,100000);
	}

	pGaussFunc->SetParameter( 2, rms );
	pGaussFunc->SetParameter( 1, center );
	pGaussFunc->SetParameter( 0, max );

	TGraph graph( binNo );

	double half = (binWidth/2.00);
	for(int i=0;i<binNo;i++){		
		float x_f = minValue+i*binWidth+half;
      float y_f = pCountTab[i];


		graph.SetPoint( i, x_f, y_f );
	}

	Double_t par[3];
	_TRACE_PRINTF_3("Fitting gauss....\n");

	// verbose :
	/*if( gPrintfLevel>=4 ){
		graph.Fit("gaus","V");
	}else{
		// quite :
		graph.Fit("gaus","Q");
	}*/


	int ret=0;
	if( gPrintfLevel>=4 ){
		ret = graph.Fit( "fit_func", "V" );
	}else{
		ret = graph.Fit( "fit_func", "Q" );
	}

	pGaussFunc->GetParameters(par);

	/*rms = (graph.GetFunction("gaus"))->GetParameter(2);
	center = (graph.GetFunction("gaus"))->GetParameter(1);
	max = (graph.GetFunction("gaus"))->GetParameter(0);*/

	rms = par[2];
	center = par[1];
	max = par[0];
	
	_TRACE_PRINTF_2("Fit_gauss (mean,rms,norm) = ( %.4f , %.4f, %.4f )\n",center,rms,max);

	gGlobalLock.UnLock();
	
	return (ret==0);
}

BOOL_T (*CMyFit::m_pFitGauss)( double minValue, double& rms, double& center,
	 			 double &max, int binNo, double binWidth,
	 			 int* pCountTab, int& nstep, int first_bin )=MyFits_FitGauss;

TRandom gRootRND;
CMyMutex gRootRNDLock;

double CMyFit::GetGaussFast( double sigma, double mean )
{
	gRootRNDLock.Lock();
	double ret = gRootRND.Gaus( mean, sigma );
	gRootRNDLock.UnLock();
	return ret;
}

void CMyFit::SetROOT_RandomSeed( int seed ){
	gRootRNDLock.Lock();
	printf("setting root seed = %d\n",seed);
	gRootRND.SetSeed( seed );
	gRootRNDLock.UnLock();
}

double CMyFit::GaussIntegral( double x0, double y0, double x1, double y1 )
{
	if( !root_func ){
		printf("Gauss function not initialized !!!\n");
		return -1;
	}
	TF2* func = (TF2*)root_func;
	double ret = func->Integral( x0, x1, y0, y1, 0.001 );
	return ret;
}

Double_t StarDistrGauss( Double_t* x, Double_t* y )
{
	Double_t valG = CMyMathFunc::GlobalGauss( x[0], x[1] );

	return valG;
}


double CMyFit::CreateGaussFunc( double x, double y, double radius ){
	if( root_func )
	{
		delete ((TF2*)root_func);
	}
	root_func = new TF2( "star_gauss", StarDistrGauss, x-radius, y-radius, x+radius, y+radius, 0 ); 	
}


#else
BOOL_T (*CMyFit::m_pFitGauss)( double minValue, double& rms, double& center,
	 			 double &max, int binNo, double binWidth,
	 			 int* pCountTab, int& nstep, int first_bin )=NULL;

double CMyFit::CreateGaussFunc( double x, double y, double radius ){
	return 0;
}

double CMyFit::GetGaussFast( double sigma, double mean )
{
	return 0;
}

double CMyFit::GaussIntegral( double x0, double y0, double x1, double y1 )
{
	return CMyMathFunc::GlobalGauss( (x0+x1)/2.00, (y0+y1)/2.00 );
}

#endif

BOOL_T (*CMyFit::m_pFillFunc)( void* event_info, int type, int ccd_idx, void* event_info2 )=NULL;



// [MS] - change on 20041003 - before 0.1 , but tracks with 0.06 
// were not checked and velocity was not checked correctly :
// now set to 0.05
// 20041005 - set to 0.02
double gMinVeloToCheck=0.02;

static void get_xy( sEventDesc* list, int cnt, double* x, double* y )
{
	for( register int i=0;i<cnt;i++){
		x[i] = list[i].x;
		y[i] = list[i].y;
	}
}

void sort_event_list( sEventDesc* list, int cnt )
{
	for( register int i=0;i<cnt-1;i++ ){
		register int min_value=list[i].timeUT;
		register int min_pos=i;
		for( register int j=i+1;j<cnt;j++){
			if( list[j].timeUT < min_value ){
				min_value = list[j].timeUT;
				min_pos = j;
			}
		}
		if( i!=min_pos ){
			sEventDesc tmp;
			tmp = list[i];
			list[i] = list[min_pos];
			list[min_pos] = tmp;
		}
	}	
}

static int find_frame( sEventDesc* list, int cnt, int frame )
{
	for( register int i=0;i<cnt;i++){
		if( list[i].frame == frame )
			return i;
	}
	return -1;
}




CMyFit::CMyFit()
{

}

CMyFit::~CMyFit()
{

}

double CMyFit::FitLineChi2( double* x_values, double* y_values, int cnt,
                            double& a, double& b )
{
	FitLine( x_values, y_values, cnt, a, b );
	double chi2 = CalcChi2( x_values, y_values, cnt, a, b );
	return chi2;
}

BOOL_T CMyFit::FitLineHorizontal( double* x_values, double* y_values, int cnt,
                                  double& c )
{
	return FALSE;
}

BOOL_T CMyFit::FitLine( double* x_values, double* y_values, int cnt,
                        double& a, double& b, int exceptPos/*=-1*/ )
{
	register double xy_sum = 0;
	register double x_sum  = 0;
	register double x2_sum = 0;
	register double y_sum  = 0; 	

	int usedCount=0;
	for(register int i=0;i<cnt;i++){
		if(exceptPos<0 || i!=exceptPos){
			xy_sum += x_values[i]*y_values[i];
			x_sum  += x_values[i];
			y_sum  += y_values[i];
			x2_sum += (x_values[i]*x_values[i]);
			usedCount++;
		}
	}

	
	double bottom = usedCount*x2_sum - x_sum*x_sum;
	a = (usedCount*xy_sum-x_sum*y_sum)/bottom;
	b = (x2_sum*y_sum-x_sum*xy_sum)/bottom;

	return TRUE;
}


double CMyFit::getLineChi2( double x, double y, double a, double b )
{
	double line_y = a*x+b;
	
	double chi2 = (y-line_y)*(y-line_y);
	return chi2;
}


double CMyFit::CalcChi2( double* x_values, double* y_values, int cnt,
                         double a, double b, int exceptPos/*=-1*/ )
{
	double sum = 0;
	for(register int i=0;i<cnt;i++){
		if(exceptPos<0 || exceptPos!=i){
			sum += getLineChi2( x_values[i], y_values[i], a, b );
		}
	}
	return sum;
}

double CMyFit::CalcMaxChi2( double* x_values, double* y_values, int cnt,
									 double a, double b, int& pos )
{
	double max_chi2 = -1000;
	
	pos = -1;
	for(register int i=0;i<cnt;i++){
		double chi2 = getLineChi2( x_values[i], y_values[i], a, b );
		if(chi2 > max_chi2){
			max_chi2 = chi2;
			pos = i;
		}
	}

	return max_chi2;
}

void CMyFit::RejectPoint( double* x_values, double* y_values, int cnt, int pos )
{
	int i=pos+1;

	printf("rejected point : (%.5f,%.5f)\n",x_values[pos],y_values[pos]);
	while(i<cnt){
		x_values[i-1] = x_values[i];
		y_values[i-1] = y_values[i];
		i++;
	}
}

BOOL_T CMyFit::FindPointsOnLine( double* x_values, double* y_values, int cnt,
										   double max_chi2, 
										   double* line_x, double* line_y, int& cnt_on_line,
										 	double& a, double& b, double& maxchi2_out )
{
	double aa,bb;
	int pos;


	double* wrk_x = new double[cnt];
	double* wrk_y = new double[cnt];

	memcpy( wrk_x, x_values, cnt*sizeof(double) );
	memcpy( wrk_y, y_values, cnt*sizeof(double) );
	
	
	FitLine( wrk_x, wrk_y, cnt, aa, bb );
	double maxchi2 = CalcMaxChi2( wrk_x, wrk_y, cnt, aa, bb, pos );

	cnt_on_line = cnt;
	while(maxchi2>max_chi2 && pos>=0 && cnt_on_line>2){
		RejectPoint( wrk_x, wrk_y, cnt_on_line, pos );
		cnt_on_line--;											
		maxchi2 = CalcMaxChi2( wrk_x, wrk_y, cnt_on_line, aa, bb, pos );
	}

	memcpy( line_x, wrk_x, cnt_on_line*sizeof(double) );
	memcpy( line_y, wrk_y, cnt_on_line*sizeof(double) );
	
	delete [] wrk_x;
	delete [] wrk_y;

	maxchi2_out = maxchi2;
	return (maxchi2<=max_chi2 && cnt_on_line>2);
}


BOOL_T CMyFit::FindPointsOnLine2( double* x_values, double* y_values, int cnt,
										   double max_chi2, 
										   double* line_x, double* line_y, int& cnt_on_line,
										 	double& a, double& b, double& maxchi2_out )
{
	double aa,bb;
	int pos;


	double* wrk_x = new double[cnt];
	double* wrk_y = new double[cnt];
	double* chi2  = new double[cnt];

	memcpy( wrk_x, x_values, cnt*sizeof(double) );
	memcpy( wrk_y, y_values, cnt*sizeof(double) );
	
	
	FitLine( wrk_x, wrk_y, cnt, aa, bb );
	double maxchi2 = CalcMaxChi2( wrk_x, wrk_y, cnt, aa, bb, pos );
	double allChi2 = CalcChi2( wrk_x, wrk_y, cnt, aa, bb );

	cnt_on_line = cnt;
	while(maxchi2>max_chi2 && pos>=0 && cnt_on_line>2){
		

		double dx1 = (x_values[1]-x_values[0]);
		double dy1 = (y_values[1]-y_values[0]);
		if(maxchi2<=500.00){
			// reject not chrionological :
			
			int bCont=1;
			while(bCont){
				bCont = 0;
				for(int k=2;k<cnt_on_line;k++){
					double dx = (wrk_x[k]-wrk_x[k-1]);		
					double dy = (wrk_y[k]-wrk_y[k-1]);

					if(dx*dx1<0 || dy*dy1<0){
						bCont = 1;						
						RejectPoint( wrk_x, wrk_y, cnt_on_line, k );
						cnt_on_line--;
						break;
					}
				}			
			}
		}

		double minChi2 = 10000000.00;		
		int minPos=-1;
		double bestA,bestB;
		for(register int i=0;i<cnt_on_line;i++){
			FitLine( wrk_x, wrk_y, cnt_on_line, aa, bb, i );
			chi2[i] = CalcChi2( wrk_x, wrk_y, cnt_on_line, aa, bb, i );	
			if(chi2[i]<minChi2){
				minChi2 = chi2[i];
				minPos = i;
				bestA = aa;
				bestB = bb;
			}
		}				


		// it was best to reject point at minPos :
		if(minPos>=0){
			RejectPoint( wrk_x, wrk_y, cnt_on_line, minPos );
			cnt_on_line--;			
			maxchi2 = CalcMaxChi2( wrk_x, wrk_y, cnt_on_line, bestA, bestB, pos );
			allChi2 = minChi2;
		}else{
			Assert(FALSE,"Could not fit line !!!");
		}
	}

	memcpy( line_x, wrk_x, cnt_on_line*sizeof(double) );
	memcpy( line_y, wrk_y, cnt_on_line*sizeof(double) );
	
	delete [] wrk_x;
	delete [] wrk_y;
	delete [] chi2;

	maxchi2_out = maxchi2;
	return (maxchi2<=max_chi2 && cnt_on_line>2);
}

double CMyFit::CalcMaxDist( sEventDesc* events, int cnt )
{
	double maxDist=-1;
	for(register int i=0;i<cnt;i++){
		for(register int j=(i+1);j<cnt;j++){
			double dist2 = ( (events[j].x-events[i].x)*(events[j].x-events[i].x)+
								  (events[j].y-events[i].y)*(events[j].y-events[i].y) );
			double dist = sqrt(dist2);
			if( dist>maxDist ){
				maxDist = dist;
			}
		}
	}
	return maxDist;
}

BOOL_T CMyFit::CheckVelocity( sEventDesc* events, int cnt,
                              double fVelocityError, int maxTimeDiff,
										int ccd_idx, eTrackCheckType_T type,
										eHistoVariableType_T histo_type,
										double* rx_min, double* ry_min )
{
	if(cnt<=2)
		return TRUE;

	// maybe to be sure add sort by timeUT here 
	// this function assumes table events is already sorted in this
	// maner ...
	int dt0 = (events[1].timeUT-events[0].timeUT);

	if( dt0!=0 ){
		double vx0 = (events[1].x-events[0].x)/dt0;
		double vy0 = (events[1].y-events[0].y)/dt0;
		double rx=0,ry=0;
			
		if( rx_min && ry_min ){
			(*rx_min) = 100000.000;
			(*ry_min) = 100000.000;
		}
	

		if( abs(dt0)>maxTimeDiff )
			return FALSE;

		for(register int i=1;i<(cnt-1);i++){	
			int dt = (events[i+1].timeUT-events[i].timeUT);

			if( dt!=0 ){
				// check only if events from different frame :
				double vx = (events[i+1].x-events[i].x)/dt;
				double vy = (events[i+1].y-events[i].y)/dt;			

				if( abs(dt)>maxTimeDiff )
	      		return FALSE;
		
				BOOL_T bCheck= CheckVelocityCondition( vx0, vy0, vx, vy, fVelocityError, rx , ry, ccd_idx, type, histo_type );
		
				if( rx_min && ry_min ){
					if( rx < (*rx_min) ){
						(*rx_min) = rx;
					}
					if( ry < (*ry_min) ){
						(*ry_min) = ry;
					}
				}

				if( !bCheck ){
					return FALSE;
				}
			}
		}	
	}
	return TRUE;
}

BOOL_T CMyFit::CheckVelocity( sEventDesc* events, int cnt,
                              sEventDesc& newEvent, 
										double fVelocityError, 
										int maxTimeDiff, eTrackCheckType_T type,
										eHistoVariableType_T histo_type )
{
	if( cnt<2 )
		return 1;
	time_t startTime=events[0].timeUT;
	time_t endTime=events[cnt-1].timeUT;
	sEventDesc* pStart=(&(events[0]));
	sEventDesc* pEnd=(&(events[cnt-1]));
	sEventDesc* pBefore=NULL;
	sEventDesc* pAfter=NULL;
	time_t beforeDT=1000000;
	time_t  afterDT=1000000;

	for(register int i=0;i<cnt;i++){
		if(events[i].timeUT<startTime){
			startTime = events[i].timeUT;
			pStart = &(events[i]);
		}
		if(events[i].timeUT>endTime){
			endTime = events[i].timeUT;
			pEnd = &(events[i]);
		}
		int dt = (newEvent.timeUT - events[i].timeUT);
		if( dt>0 && dt<beforeDT ){
			beforeDT = dt;
			pBefore = &(events[i]);
		}else{
			dt = -dt;
			if( dt>0 && dt<afterDT ){
				afterDT = dt;
				pAfter = &(events[i]);
			}
		}
	}
	if( pStart && pEnd && (pAfter || pBefore)){
		time_t dt = (pEnd->timeUT - pStart->timeUT);
		double vx = (pEnd->x - pStart->x)/dt;
		double vy = (pEnd->y - pStart->y)/dt;

		
		double vx_new,vy_new;
		if( pAfter ){
			if( abs(afterDT)>maxTimeDiff  ){
				return FALSE;
			}

			vx_new = (pAfter->x-newEvent.x)/afterDT;
			vy_new = (pAfter->y-newEvent.y)/afterDT;
			/*if( !CheckVelocityCondition( vx, vy, vx_new, vy_new, fVelocityError, 
												  DEFAULT_CCD_IDX_FOR_CHECK_COND, type, histo_type) ){
				return FALSE;
			}*/

			// NEW - 20041118 - velocity checked only if >1 pixels difference
			if( fabs(pAfter->y-newEvent.y)>1 && 
				 !CheckVelocityCondition( 1.00, vy, 1.00, vy_new, fVelocityError, 
												  DEFAULT_CCD_IDX_FOR_CHECK_COND, type, histo_type) ){
				return FALSE;
			}
			if( fabs(pAfter->x-newEvent.x)>1 && 
				 !CheckVelocityCondition( vx, 1.00, vx_new, 1.00, fVelocityError, 
												  DEFAULT_CCD_IDX_FOR_CHECK_COND, type, histo_type) ){
				return FALSE;
			}

		}
		if( pBefore ){
			vx_new = (newEvent.x-pBefore->x)/beforeDT;
			vy_new = (newEvent.y-pBefore->y)/beforeDT;
			if( abs(beforeDT)>maxTimeDiff  ){
				return FALSE;
			}

			/*if( !CheckVelocityCondition( vx, vy, vx_new, vy_new, fVelocityError, 
												  DEFAULT_CCD_IDX_FOR_CHECK_COND, type, histo_type ) ){ 
				return FALSE;
			}*/

			// NEW - 20041118 - velocity checked only if >1 pixels difference
			if( fabs(pBefore->y-newEvent.y)>1 && 
				 !CheckVelocityCondition( 1.00, vy, 1.00, vy_new, fVelocityError, 
												  DEFAULT_CCD_IDX_FOR_CHECK_COND, type, histo_type) ){
				return FALSE;
			}
			if( fabs(pBefore->x-newEvent.x)>1 && 
				 !CheckVelocityCondition( vx, 1.00, vx_new, 1.00, fVelocityError, 
												  DEFAULT_CCD_IDX_FOR_CHECK_COND, type, histo_type) ){
				return FALSE;
			}
		}
	}
	return TRUE;
}

BOOL_T CMyFit::CheckVelocityCondition( double vx, double vy,
                                       double vx_new, double vy_new,
                                       double error, int ccd_idx,
													eTrackCheckType_T type,
													eHistoVariableType_T histo_type )
{
/*	if( vx_new*vx<0 || vy_new*vy<0 || fabs( vx_new-vx )>error || fabs( vy_new-vy )>error ){
		return FALSE;
	}*/
	double vx_min=MIN(fabs(vx),fabs(vx_new));
	double vy_min=MIN(fabs(vy),fabs(vy_new));
	double vx_max=MAX(fabs(vx),fabs(vx_new));
	double vy_max=MAX(fabs(vy),fabs(vy_new));	


	if( m_pFillFunc ){
		double rx = ( vx_min/vx_max )*mysign(vx*vx_new);
		if(vx_max<gMinVeloToCheck)
			rx = 1.00;


		double ry = ( vy_min/vy_max )*mysign(vy*vy_new);
		if(vy_max<gMinVeloToCheck)
			ry = 1.00;

		(*m_pFillFunc)( &rx, eRXall, ccd_idx, NULL );
		(*m_pFillFunc)( &ry, eRYall, ccd_idx, NULL );
		(*m_pFillFunc)( &vx, eVXvsVX, ccd_idx, &vx_new );
		(*m_pFillFunc)( &vy, eVYvsVY, ccd_idx, &vy_new );

		if(gDebugTracks>1)printf("(rx,ry)=(%.2f,%.2f)\n",rx,ry);
	}

	if(gDebugTracks>1)printf("VELO : (vx_min,vx_max) = (%.2f,%.2f), (vy_min,vy_max) = (%.2f,%.2f) =>",vx_min,vx_max,vy_min,vy_max);

	// if max is slower then gMinVeloToCheck pixel / sec - do not perform this check :
	if( ( vx_max>=gMinVeloToCheck && ( vx_new*vx<0 || (vx_min/vx_max)<error ) ) || // bad vx 	
		 ( vy_max>=gMinVeloToCheck && ( vy_new*vy<0 || (vy_min/vy_max)<error ) ) ){ // or bad vy
		if(gDebugTracks>1)printf("return FALSE\n");
      return FALSE;
	} 


	/*if( (vx_new*vx<0 && vx_max>=gMinVeloToCheck ) || 
		 (vy_new*vy<0 && vy_max>=gMinVeloToCheck ) || 
		 ( (vx_min/vx_max)<error && vx_max>=gMinVeloToCheck ) ||  // special for very slow 
		 ( (vy_min/vy_max)<error && vy_max>=gMinVeloToCheck ) ){  // special for very slow
		if(gDebugTracks>1)printf("return FALSE\n");
		return FALSE;
	}*/

	if(gDebugTracks>1)printf("return TRUE\n");
	
	return TRUE;	
}


BOOL_T CMyFit::CheckVelocityCondition( double vx, double vy, 
													  double vx_new, double vy_new, 
													  double error,
													  double& rx, double& ry,
													  int ccd_idx, 
													  eTrackCheckType_T type,
													  eHistoVariableType_T histo_type )
{
	double vx_min=MIN(fabs(vx),fabs(vx_new));
	double vy_min=MIN(fabs(vy),fabs(vy_new));
	double vx_max=MAX(fabs(vx),fabs(vx_new));
	double vy_max=MAX(fabs(vy),fabs(vy_new));	

	

	if( m_pFillFunc ){
		double rx = ( vx_min/vx_max )*mysign(vx*vx_new);
		if(vx_max<gMinVeloToCheck)
			rx = 1.00;


		double ry = ( vy_min/vy_max )*mysign(vy*vy_new);
		if(vy_max<gMinVeloToCheck)
			ry = 1.00;

		(*m_pFillFunc)( &rx, eRXall, ccd_idx, NULL );
		(*m_pFillFunc)( &ry, eRYall, ccd_idx, NULL );
		(*m_pFillFunc)( &vx, eVXvsVX, ccd_idx, &vx_new );
		(*m_pFillFunc)( &vy, eVYvsVY, ccd_idx, &vy_new );

		if(gDebugTracks>1)printf("(rx,ry)=(%.2f,%.2f)\n",rx,ry);
	}

	if(gDebugTracks>1)printf("VELO : (vx_min,vx_max) = (%.2f,%.2f), (vy_min,vy_max) = (%.2f,%.2f) =>",vx_min,vx_max,vy_min,vy_max);

	// if max is slower then gMinVeloToCheck pixel / sec - do not perform this check :

	if( vx_max!=0 ){
		rx = ( vx_min/vx_max )*mysign(vx*vx_new); // mysign NEW [20041006]
	}else{
		rx = 1.00;
	}
	if( vy_max!=0 ){
	   ry = ( vy_min/vy_max )*mysign(vy*vy_new); // mysign NEW [20041006]
	}else{
		ry = 1.00;
	}

	if( ( vx_max>=gMinVeloToCheck && ( vx_new*vx<0 || rx<error ) ) || // bad vx 	
		 ( vy_max>=gMinVeloToCheck && ( vy_new*vy<0 || ry<error ) ) ){ // or bad vy
		if(gDebugTracks>1)printf("return FALSE\n");
      return FALSE;
	} 

	if( gDebugTracks>1 ){
		printf("(vx_min,vx_max) = (%.2f,%.2f), (vy_min,vy_max) = (%.2f,%.2f) ==> (rx,ry) = (%.2f,%.2f)\n",vx_min,vx_max,vy_min,vy_max,rx,ry);
	}

	/*if( (vx_new*vx<0 && vx_max>=gMinVeloToCheck ) || 
		 (vy_new*vy<0 && vy_max>=gMinVeloToCheck ) || 
		 ( (vx_min/vx_max)<error && vx_max>=gMinVeloToCheck ) ||  // special for very slow 
		 ( (vy_min/vy_max)<error && vy_max>=gMinVeloToCheck ) ){  // special for very slow
		if(gDebugTracks>1)printf("return FALSE\n");
		return FALSE;
	}*/
	if(gDebugTracks>1)printf("return TRUE\n");
	
	return TRUE;	

}


//	FindBest3Points( wrk_x, wrk_y, frames, cnt, a, b, best3pos, minChi2_for3 );
int CMyFit::FindBest3Points( sEventDesc* events, int cnt,
									  double& a, double& b, 
									  int* best3pos, double& minChi2_for3,
									  BOOL_T bCheckVelocity, double fVelocityError,
									  int ccd_idx, BOOL_T bAddFromSame,
									  eTrackCheckType_T type )
{
	double pointsX[3],pointsY[3];
	int nTry=0;
	int nOK=0;
	double rx_min, ry_min;

	for(register int i=0;i<=cnt-3;i++){
		pointsX[0]=events[i].x;
		pointsY[0]=events[i].y;

		for(register int j=(i+1);j<=(cnt-2);j++){
			if( events[j].frame==events[i].frame && !bAddFromSame ){
				// same frame 
				continue;
			}

			pointsX[1]=events[j].x;
	      pointsY[1]=events[j].y;
			for(register int k=j+1;k<=(cnt-1);k++){
				if( (events[j].frame==events[k].frame ||
					  events[j].frame==events[i].frame ) && !bAddFromSame ){
					// same frame
					continue;
				}


				// have 3 points i,j,k now check them :
				double a3,b3;
				pointsX[2]=events[k].x;
		      pointsY[2]=events[k].y;

				if( bCheckVelocity ){
					sEventDesc tocheck[3];
					tocheck[0]=events[i];
					tocheck[1]=events[j];
					tocheck[2]=events[k];
               BOOL_T bCheck = CheckVelocity( tocheck, 3,fVelocityError, DEFAULT_MAX_DIFF_TIME, ccd_idx, type, eHistoVXRatioTo3, &rx_min, &ry_min );
	
					if( m_pFillFunc ){ 
						(*m_pFillFunc)( &rx_min, eHistoVXRatioTo3, ccd_idx, NULL );
		            (*m_pFillFunc)( &ry_min, eHistoVYRatioTo3, ccd_idx, NULL );	
					}
					
					if( !bCheck ){
						// 3 point do not satisfy velocity criteria 
						continue;
					}
            }

				FitLine( pointsX, pointsY, 3, a3, b3 );
				double chi2_3 = CalcChi2( pointsX, pointsY, 3, a3, b3 );

				if( gShowChi2Points3 && chi2_3<1000.00 ){
					printf("chi2_3points = %.4f\n",chi2_3);fflush(0);
				}

				if(chi2_3<minChi2_for3 || nTry==0){					
					minChi2_for3 = chi2_3;
					best3pos[0]=i;
					best3pos[1]=j;
					best3pos[2]=k;		
					a = a3;
					b = b3;
					nOK++;

					// just for testing reasons :
					/*sEventDesc tocheck[3];
               tocheck[0]=events[i];
               tocheck[1]=events[j];
               tocheck[2]=events[k];
               BOOL_T bVelo = CheckVelocity( tocheck, 3,fVelocityError );*/
				}
				nTry++;
			}
		}
	}
	return nOK;
}

double CMyFit::ReCalcMaxChi2( double maxDist, double maxchi2 )
{
	// return maxchi2;

	int r = (maxDist/100)+1;
	double ret = r*maxchi2;
	if(ret>100){
		ret = 100;
	}
	return ret;
}

BOOL_T CMyFit::FindPointsOnLine3New( sEventDesc* events, int cnt,
                             double max_chi2_per_point,
                             sEventDesc* events_on_line,int& cnt_on_line,
                             double& a, double& b, double& maxchi2_out,
                             sEventDesc* events_rejected, int& rejected_cnt,
									  BOOL_T bCheckVelocity, double fVelocityError,
									  double& minChi2PerEvent3Points, int ccd_idx,
									  BOOL_T bAddFromSame, eTrackCheckType_T type )
{
	double aa,bb;
   int pos;
	maxchi2_out = 0;
                                                                                
   if(cnt<3){
      return FALSE;
   }

	int cnt_original=cnt;
	rejected_cnt=0;                                                                                
   cnt_on_line=0;
	sEventDesc* wrk = new sEventDesc[cnt];
	memcpy( wrk, events, cnt*sizeof(sEventDesc) );	
	int best3[3];
   double minChi2_for3=10000000000.00;
   int nOK = FindBest3Points( wrk, cnt, a, b, best3, minChi2_for3,
										bCheckVelocity, fVelocityError, ccd_idx,
									   bAddFromSame, type );

	/*if(nOK>0){
		sEventDesc best[3];
		best[0] = wrk[best3[0]];
		best[1] = wrk[best3[1]];
		best[2] = wrk[best3[2]];
		double maxDist = CalcMaxDist( best, 3 );
		max_chi2_per_point = ReCalcMaxChi2( maxDist, max_chi2_per_point );
	}*/

	minChi2PerEvent3Points = (minChi2_for3/3);
	if(minChi2PerEvent3Points<=max_chi2_per_point){
		if( gShowChi2OfAdded ){
      	printf("Chi2_of_added %d %.8f\n",3,minChi2PerEvent3Points);
      }

		// now use 3 points as seed and add next :
		for(register int ii=0;ii<3;ii++){
			int bestpos=best3[ii];
			events_on_line[ii] = wrk[ bestpos ];
		}
		cnt_on_line = 3;

		// best3 is sorted by FindBest3Points
		for(int ii=2;ii>=0;ii--){
			// assuming sorted 
			int bestpos=best3[ii];
			for(int jj=bestpos;jj<(cnt-1);jj++){
				wrk[ jj ] = wrk[ jj+1 ];
			}
			cnt--;
		}

		// now try adding points until only very bad left :	
		int new_pos=3;
		double* line_x = new double[cnt_original];
		double* line_y = new double[cnt_original];

		while(cnt>0){
			int checkpos=cnt-1;
			BOOL_T bAdded=FALSE;
			BOOL_T bVelOK=TRUE;
			events_on_line[new_pos] = wrk[checkpos];
			if( (bAddFromSame || find_frame( events_on_line, cnt_on_line, wrk[checkpos].frame )<0) && 
				 (!bCheckVelocity || CheckVelocity( events_on_line, cnt_on_line, wrk[checkpos], fVelocityError, DEFAULT_MAX_DIFF_TIME, type, eHistoVXRatioToOld ) )){

				get_xy( events_on_line, cnt_on_line+1, line_x, line_y );
				FitLine( line_x, line_y, new_pos+1, aa, bb );
				double chi2 = CalcChi2( line_x, line_y, new_pos+1, aa, bb );

				if( gShowChi2OfAdded ){
					printf("Chi2_of_added %d %.8f\n",(cnt_on_line+1),chi2);
				}

				// skip point
				if( (chi2/(cnt_on_line+1)) <= max_chi2_per_point){
					// add to line if good 
					new_pos++;
					cnt_on_line++;					

					// update best line coeficents :
					a = aa;
   	         b = bb;
					bAdded = TRUE;
				}
			}
			if( !bAdded ){
				events_rejected[rejected_cnt] = events_on_line[new_pos];
				rejected_cnt++;
			}
			cnt--;
		}

		delete [] line_x;
		delete [] line_y;
	}

	
	delete [] wrk;
	return (cnt_on_line>2);
}

//
// 1/ finding 3 best points - and adding next :
//
BOOL_T CMyFit::FindPointsOnLine3( double* x_values, double* y_values, 
										 	int cnt,
										   double max_chi2_per_point, 
										   double* line_x, double* line_y, int& cnt_on_line,
										 	double& a, double& b, double& maxchi2_out,
											double* rejected_x, double* rejected_y, int& rejected_cnt )
{
	double aa,bb;
	int pos;

	if(cnt<3){
		return FALSE;
	}

	cnt_on_line=0;

	double* wrk_x = new double[cnt];
	double* wrk_y = new double[cnt];


	memcpy( wrk_x, x_values, cnt*sizeof(double) );
	memcpy( wrk_y, y_values, cnt*sizeof(double) );


	//printf("POTENTIAL TRACK POINTS :\n");
	//for(register int i=0;i<cnt;i++){
	//	printf("%d %d\n",(int)wrk_x[i],(int)wrk_y[i]);
	//}    	
	//printf("\n");

	int best3pos=-1;
   double minChi2_for3=10000000000.00;
	// FindBest3Points( wrk_x, wrk_y, frames, cnt, a, b, best3pos, minChi2_for3 );
		
	for(register int i=0;i<=cnt-3;i++){
		double a3,b3;
		FitLine( wrk_x+i, wrk_y+i, 3, a3, b3 );
		double chi2_3 = CalcChi2( wrk_x+i, wrk_y+i, 3, a3, b3 );
		if(chi2_3<minChi2_for3 || i==0){
			minChi2_for3 = chi2_3;
			best3pos = i;

			// update best coeficients here :
			a = a3;
			b = b3;
		}
	}

	if(best3pos<0){
		delete [] wrk_x;
	   delete [] wrk_y;
		return FALSE;
	}

	Assert(best3pos>=0,"No 3 points found (all=%d) ???",cnt);
	if((minChi2_for3/3)<=max_chi2_per_point){
		// now use 3 points as seed and add next :
		memcpy( line_x, wrk_x+best3pos, 3*sizeof(double) );
		memcpy( line_y, wrk_y+best3pos, 3*sizeof(double) );
		cnt_on_line = 3;
	
		int numToCopy = (cnt-(best3pos+3));
		if(numToCopy>0){
			//memcpy( wrk_x+best3pos, wrk_x+best3pos+3, numToCopy*sizeof(double) );
			//memcpy( wrk_y+best3pos, wrk_y+best3pos+3, numToCopy*sizeof(double) );
			for(register int cc=-0;cc<numToCopy;cc++){
				wrk_x[best3pos+cc] = wrk_x[best3pos+3+cc];
				wrk_y[best3pos+cc] = wrk_y[best3pos+3+cc];
			}
		}
		cnt = cnt - 3;


		// now try adding points until only very bad left :	
		int new_pos=3;
		while(cnt>0){
			int checkpos=cnt-1;
			line_x[new_pos] = wrk_x[checkpos];
			line_y[new_pos] = wrk_y[checkpos];

			FitLine( line_x, line_y, new_pos+1, aa, bb );
			double chi2 = CalcChi2( line_x, line_y, new_pos+1, aa, bb );

			// skip point
			cnt--;
			if( (chi2/(cnt_on_line+1)) <= max_chi2_per_point){
				// add to line if good 
				new_pos++;
				cnt_on_line++;

				// update best line coeficents :
				a = aa;
            b = bb;
			}else{
				rejected_x[rejected_cnt] = line_x[new_pos];
				rejected_y[rejected_cnt] = line_y[new_pos];
				rejected_cnt++;
			}
		}
	}


	
	delete [] wrk_x;
	delete [] wrk_y;


	/*if(cnt_on_line>2){
		printf("POINTS ON TRACK :\n");
		for(register int i=0;i<cnt_on_line;i++){
			printf("%d %d\n",(int)line_x[i],(int)line_y[i]);
		}
		printf("line = %.5f * x + %.5f\n",a,b);
	}*/

	// maxchi2_out = maxchi2;
	return (cnt_on_line>2);
}

double CMyFit::CalcDist2FromLine( double a, double b, double c,
								 			 double x, double y )
{
	double up = CMyMathFunc::mysqr( a*x+b*y+c );
	double bottom = ( a*a+b*b ); 
	double ret = (up/bottom);
	
	return ret;
}

double CMyFit::CalcDist2FromLine2Par( double a, double b,
                                      double* x, double* y, int cnt )
{
	double sum=0;
	for(register int i=0;i<cnt;i++){
		sum += CalcDist2FromLine2Par( a, b, x[i], y[i] );
	}
	return sum;
}

double CMyFit::CalcDist2FromLine2Par( double a, double b,
                                      double x, double y )
{
	return CalcDist2FromLine( a, -1 , b , x, y );
}


BOOL_T CMyFit::TryToAddNewPointsNew( sEventDesc* rejected, int rejected_cnt,
											double max_chi2_per_point,
                                 sEventDesc* line, int& cnt_on_line,
                                 double& a, double& b,  double& maxchi2_out,
											int total_count,
										   BOOL_T (*fill_func)( void*, int, int, void* ),
											int cam_idx, BOOL_T bAddFromSame,
											eTrackCheckType_T type )
{
	double original_a = a;
	double original_b = b;

	int added_total=0;

	maxchi2_out = 0.00;

	if(rejected_cnt>0 && cnt_on_line>2){
		double* line_x = new double[ total_count ];	
		double* line_y = new double[ total_count ];

		// re-trying to fit rest :
		int added=1;
		int cnt_on_line_save = cnt_on_line;
		int iter=0;
		while(added && rejected_cnt && iter<1000){
			added=0;

			double min_chi2=100000.00;			
			int min_pos = -1;
			for(int i=0;i<rejected_cnt;i++){
				if( bAddFromSame || find_frame( line, cnt_on_line, rejected[i].frame )<0 ){
					// chosing points from frames not already belonging to track :
					double point_chi2 = CalcDist2FromLine( a, -1, b, rejected[i].x, rejected[i].y );
					if(point_chi2<min_chi2){
						min_pos = i;
						min_chi2 = point_chi2;
					}				
				}
			}

			if( fill_func ){
				// 2 - is eChi2ToOld - to add to histogram when
				// adding new points to exisiting track :
				if( type==eNormalTrack ){
					(*fill_func)( &min_chi2, eChi2ToOld, cam_idx, NULL );
				}else{
					if( DoAddFromSameFrame( type ) ){
						(*fill_func)( &min_chi2, eChi2ToOldSum, cam_idx, NULL );
					}
				}
			}

			if(min_chi2<max_chi2_per_point && min_pos>=0){
				// in case point satisfying chi2 criteria found - check
				// how new fit line will match :

				if( min_chi2 < max_chi2_per_point ){				
					double aa,bb;

					line[cnt_on_line] = rejected[min_pos];
					get_xy( line, cnt_on_line+1, line_x, line_y ); 

					FitLine( line_x, line_y, cnt_on_line+1, aa, bb );
	   	      double chi2 = CalcChi2( line_x, line_y, cnt_on_line+1, aa, bb );				
					// if( (chi2/(cnt_on_line+1)) <= max_chi2_per_point){
						cnt_on_line++;
						added++;
						added_total++;
						for(int j=min_pos;j<(rejected_cnt-1);j++){
							rejected[j] = rejected[j+1];
						}
						rejected_cnt--;
						a = aa;
						b = bb;
					// }	
				}
			}
			iter++;
		}
		
		delete [] line_x;
		delete [] line_y;
	}	

	return (added_total>0);	
}

BOOL_T CMyFit::TryToAddNewPoints( double* rejected_x, double* rejected_y, int rejected_cnt,
 											 double max_chi2_per_point,
											 double* line_x, double* line_y, int& cnt_on_line,
											 double& a, double& b,  double& maxchi2_out )
{
	double original_a = a;
	double original_b = b;

	int added_total=0;

	maxchi2_out = 0.00;

	if(rejected_cnt>0 && cnt_on_line>2){
		// re-trying to fit rest :
		int added=1;
		int cnt_on_line_save = cnt_on_line;
		int iter=0;
		while(added && rejected_cnt && iter<1000){
			added=0;

			double min_chi2=100000.00;			
			int min_pos = -1;
			for(int i=0;i<rejected_cnt;i++){
				double point_chi2 = CalcDist2FromLine( a, -1, b, rejected_x[i], rejected_y[i] );
				if(point_chi2<min_chi2){
					min_pos = i;
					min_chi2 = point_chi2;
				}				
			}

			if(min_chi2<max_chi2_per_point && min_pos>=0){
				// in case point satisfying chi2 criteria found - check
				// how new fit line will match :

				if( min_chi2 < max_chi2_per_point ){				
					double aa,bb;

					line_x[cnt_on_line] = rejected_x[min_pos];
	            line_y[cnt_on_line] = rejected_y[min_pos];

					FitLine( line_x, line_y, cnt_on_line+1, aa, bb );
	   	      double chi2 = CalcChi2( line_x, line_y, cnt_on_line+1, aa, bb );				
					// if( (chi2/(cnt_on_line+1)) <= max_chi2_per_point){
						cnt_on_line++;
						added++;
						added_total++;
						for(int j=min_pos;j<(rejected_cnt-1);j++){
							rejected_x[j] = rejected_x[j+1];
							rejected_y[j] = rejected_y[j+1];
						}
						rejected_cnt--;
						a = aa;
						b = bb;
					// }	
				}
			}
			iter++;
		}
	}	

	return (added_total>0);	
}




BOOL_T CMyFit::FitGauss( double minValue, double& rms, double& center, 
							 double &max, int binNo, double binWidth, 
							 int* pCountTab, int& nstep, int first_bin)
{
	if( CMyFit::m_pFitGauss ){
		return (*CMyFit::m_pFitGauss)( minValue, rms, center, max, binNo, binWidth, pCountTab, nstep,  first_bin);
	}

	if(binNo>100){		
		 return FALSE;
	}
	return FALSE;
}


void show_points( const char* cmt, sEventDesc* list, int cnt )
{
	printf("%s",cmt);
	for(int ii=0;ii<cnt;ii++){
		double rx=0,ry=0;
		if(ii>0){
			double dx = (list[ii].x-list[ii-1].x);
			double dy = (list[ii].y-list[ii-1].y);
			double dt = (list[ii].timeUT-list[ii-1].timeUT);
			
			if( fabs(dt)>0 ){
				rx = (dx/dt);
				ry = (dy/dt);
			}
		}

 		printf("%d-(%.2f,%.2f)-(%.2f,%.2f),",list[ii].frame,list[ii].x,list[ii].y,rx,ry);
   }
   printf("\n");
}


BOOL_T CMyFit::DoAddFromSameFrame( eTrackCheckType_T track_type )
{
   return ( track_type==ePlaneTrack || track_type==eTrackOnSumedFrame );
}
