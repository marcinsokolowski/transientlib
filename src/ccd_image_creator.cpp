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
#include "ccd_image_creator.h" 
#include <cexcp.h>
#include <math.h>
#include <random.h>
#include <myparser.h>
#include "ccd_globals.h"
#include "ccd_trace.h"
#include "ccd.h"
#include <cfg.h>
#include <fits_file.h>
#include <mystrtable.h>
#include "mcinput.h"
#include "ccd_samples.h"
#include "ccd_pipeline.h"
#include "ccd_datares.h"
#include "ccd_analyse.h"
#include "ccd_util.h"
#include <mypixellist.h>
#include <myutil.h>
#include "gendistr.h"



/*double CImageCreator::m_ObjProb=1;
BOOL_T CImageCreator::m_bInit=FALSE;
BOOL_T CImageCreator::m_bFromSampleFile=TRUE;
mystring CImageCreator::m_MaxBright=1;
mystring CImageCreator::m_MinBright=1;
mystring CImageCreator::m_szListname="object_list";
vector<CListFile*> CImageCreator::m_ListFiles;
CccdSamples* CImageCreator::m_pObjectSamples=NULL;*/



//void gRefreshParams(void){
//	CImageCreator::RefreshParams();
//}

//void (*CVaryParamDesc::ParamRefreshFunc1)(void) = gRefreshParams;


//-----------------------------------------------------------------------
// class for cleaing globa pointers :
class CImageCreatorGlobalCleaner
{
public:
	CImageCreatorGlobalCleaner(){};
	~CImageCreatorGlobalCleaner();
};

CImageCreatorGlobalCleaner::~CImageCreatorGlobalCleaner()
{
	/*if(CImageCreator::m_pObjectSamples){
		delete (CImageCreator::m_pObjectSamples);
		CImageCreator::m_pObjectSamples = NULL;
	}*/
}

static CImageCreatorGlobalCleaner gImageCreatorGlobalCleaner;
//--------------------------------------------------------------------------

CSamplePutInfo::CSamplePutInfo( int x, int y, int frame_index, int cam_idx, 
											const char* szMag, double prevTime, 
											int last_pos, int last_sample_idx)
: m_X(x),m_Y(y),m_FrameIndex(frame_index),m_CamIdx(cam_idx),m_szMag(szMag),
  m_prevTime(prevTime), m_LastPos( last_pos ), m_LastSampleIdx( last_sample_idx )
{
}

CSamplePutInfo::CSamplePutInfo( const CSamplePutInfo& right )
{
	m_X = right.m_X;
	m_Y = right.m_Y;
	m_FrameIndex = right.m_FrameIndex;
	m_CamIdx = right.m_CamIdx;
	m_szMag = right.m_szMag;
	m_prevTime = right.m_prevTime;
	m_LastPos = right.m_LastPos;
	m_LastSampleIdx = right.m_LastSampleIdx;
}

void CImageCreator::InitParams()
{
	gCCDParams.InitParams();
	if(!m_bInit){
		m_MaxBright = GetGlobalParam("CCD_MC_MAX_BRIGHT",TRUE);
		m_MinBright = GetGlobalParam("CCD_MC_MIN_BRIGHT",TRUE);
		m_szListname = (m_pPipeline->GetPipelineCfg()).m_szListName.c_str();
		m_bInit = TRUE;
	}
}

void CImageCreator::RefreshParams()
{
	m_bInit = FALSE;
	gCCDParams.RefreshParams();
	InitParams();
}


CImageCreator::CImageCreator( CCDPipeline* pPipeline )
: m_SampleIndex(0), m_pPrevPutObj(NULL), m_LastPutObjectX(-1), m_LastPutObjectY(-1), 
  m_PrevValue(0), m_pGenEvent(NULL), m_pObjectSamples(NULL),
  m_ObjProb(1),m_bInit(FALSE),m_bFromSampleFile(TRUE),m_MaxBright(1),
  m_MinBright(1),m_szListname("object_list"), m_pPipeline( pPipeline ),
  m_LastPutObjectRA(0),m_LastPutObjectDEC(0), m_LastPos(0), m_LastSampleIdx(0)
{
	InitParams();
	if(gCCDParams.m_nPutObjOnNFrames>0){
		// size of cyclic table is number of frames to re-put sample times
		// number of samples to put on single frame :		

// changed on 20060418 - due to fact that tests without confirmation on 
// next frames also put sample on next frame, currently 
// it is based on value of parameter which sets how many frames is needed
// to confirm event :
	int prev_samples = gCCDParams.m_ConfirmEventsOnNextNFrames * gCCDParams.m_nSamplesToPutOnFrame;
//		int prev_samples = gCCDParams.m_nPutObjOnNFrames * gCCDParams.m_nSamplesToPutOnFrame;

//		m_pPrevPutObj = new CCycleTab2D<CSamplePutInfo>( gCCDParams.m_nPutObjOnNFrames );
		m_pPrevPutObj = new CCycleTab2D<CSamplePutInfo>( prev_samples );
	}
}

CImageCreator::~CImageCreator()
{
	if(m_pPrevPutObj){
		delete m_pPrevPutObj;
	}
}

CListFile* CImageCreator::FindListFile( const char* fname )
{
	vector<CListFile*>::iterator i;
	for(i=m_ListFiles.begin();i!=m_ListFiles.end();i++){
		if(strcmp((*i)->GetFileName(),fname)==0)
			return (*i);
	}
	return NULL;
}

void CImageCreator::ModifyImage(CCDMatrix& Image)
{
	
}

// randomly put Mion signal on the image
void CImageCreator::AddMuon(CCDMatrix& Image)
{}

// adds some noise to image
void CImageCreator::AddNoise(CCDMatrix& Image)
{}

// adds new star to image
void CImageCreator::AddStar(CCDMatrix& Image,Table2D<ELEM_TYPE>& StarImage)
{}

// adds GRB
void CImageCreator::AddGRBBurst(CCDMatrix& Image,Table2D<ELEM_TYPE>& GRBImage)
{}

LONG_T CImageCreator::GetGauss(LONG_T x0,LONG_T y0,LONG_T r0,LONG_T x,LONG_T y,LONG_T max_value)
{
	Assert(r0>0,"r0 must be greater then 0");
	double ret = double(max_value)*exp(-((x-x0)*(x-x0)+(y-y0)*(y-y0))/(r0*r0));
	return (LONG_T)ret;
}

void CImageCreator::InitSimpleStar(Table2D<ELEM_TYPE>& Image,LONG_T nPoints,
                                   LONG_T max_value,LONG_T x0,LONG_T y0,LONG_T r0,
											  LONG_T step)
{
	double r;
	CGenDistr distr;
	LONG_T x_size = Image.GetXSize();
	LONG_T y_size = Image.GetYSize();
	

	if(x0==-1)
		x0 = x_size/2;
	if(y0==-1)
      y0 = y_size/2;
	if(r0==-1)
		r0 = (LONG_T)(sqrt(double(x_size)*double(x_size)/4.00 + 
					        double(y_size)*double(y_size)/4.00)*0.2);

	for(int i=0;i<nPoints;i++){
		double x,y;
		distr.GetGauss( 1.00, r0, x0, y0, x ,y );  		
		LONG_T x_l = (LONG_T)x;
		LONG_T y_l = (LONG_T)y;
		if(x_l>=0 && x_l<Image.GetXSize() &&
         y_l>=0 && y_l<Image.GetYSize() )
			Image.GetPixel(x_l,y_l)+=step;
	}
}

// creates simple simulation star, according to Gauss shape
BOOL_T CImageCreator::AddSimpleStar(Table2D<ELEM_TYPE>& Image,double brightness,       
                                    LONG_T width,LONG_T x0,LONG_T y0)
{
	LONG_T nPoints = GetNumberOfPoints( brightness );
	if(x0<0 || x0>=Image.GetXSize())
		return FALSE;
	if(y0<0 || y0>=Image.GetYSize())
		return FALSE;
	
	double r0;
	CGenDistr distr;
	LONG_T x_size = Image.GetXSize();
	LONG_T y_size = Image.GetYSize();
	

	if(x0==-1)
		x0 = x_size/2;
	if(y0==-1)
      y0 = y_size/2;
	if(r0==-1)
		r0 = (LONG_T)(sqrt(double(x_size)*double(x_size)/4.00 + 
					        double(y_size)*double(y_size)/4.00)*0.2);

	for(int i=0;i<nPoints;i++){
		double x,y;
		distr.GetGauss( 1.00, r0, x0, y0, x ,y );  		
		LONG_T x_l = (LONG_T)x;
		LONG_T y_l = (LONG_T)y;
		if(x_l>=0 && x_l<Image.GetXSize() &&
         y_l>=0 && y_l<Image.GetYSize() )
			Image.GetPixel(x_l,y_l)+=1;
	}	
}


BOOL_T CImageCreator::AddSimpleStar(Table2D<ELEM_TYPE>& Image,double brightness,LONG_T width)
{
	CRandom rnd;	
	LONG_T x0 = (LONG_T)(Image.GetXSize()*rnd.GetRandom());
	LONG_T y0 = (LONG_T)(Image.GetYSize()*rnd.GetRandom());	

	AddSimpleStar( Image, brightness, width, x0, y0 );
	return TRUE;
}

void CImageCreator::AddSampleObject( CCDMatrix& Image, Table2D<ELEM_SAMPLE_TYPE>& SampleObject, const char* szMag )
{
	InitParams();
	LONG_T pos_x = (LONG_T)(Image.GetXSize()*(CRandom::GetRandom()));
	LONG_T pos_y = (LONG_T)(Image.GetYSize()*(CRandom::GetRandom()));
	
	m_pGenEvent = Image.AddSample( pos_x, pos_y, SampleObject, 0,0,-1,-1, TRUE, szMag );
}

void CImageCreator::ScaleSample( const char* szMagnitude, 
											const Table2D<ELEM_SAMPLE_TYPE>& SampleObject,
                    				   Table2D<ELEM_SAMPLE_TYPE>& ScaledSample,
											const char* szInSampleMag )
{
/*	ScaledSample = SampleObject;		


	double maxS1 = ScaledSample.m_MaxVal;
	double theoS1 = theo_s1( atof(szMagnitude) );
	double r = (theoS1/maxS1);
	ScaledSample.MultiplyByConst( r );

	static int num=0;
	printf("Sample of magnitude %s scaled to magnitude %s , by factor %f=(%f/%f)\n",szInSampleMag,szMagnitude, r,theoS1,maxS1);
	mystring szName;	
	szName << szInSampleMag << "_scaled_to_" << szMagnitude << "_" << num <<".fit";
	num++;
	CCDUtil::WriteToFITSFile( ScaledSample, szName.c_str() );		*/
}

BOOL_T CImageCreator::GetSampleObject( Table2D<ELEM_SAMPLE_TYPE>& SampleObject )
{
	InitSampleCache();
	int cnt = m_pObjectSamples->GetCount();
	m_LastPos = CRandom::GetRandomInteger( 0, cnt );
	CMagSamples* pMagSamples = m_pObjectSamples->GetMagSamples( m_LastPos );
	Assert( pMagSamples!=NULL , "Sample at pos : %d not found",m_LastPos);
	SampleObject = pMagSamples->GetRandomSample( &m_LastSampleIdx );
	return TRUE;	
}

BOOL_T CImageCreator::GetSampleObject( Table2D<ELEM_SAMPLE_TYPE>& SampleObject, const char* fname, const char* szMagnitude )
{
	InitSampleCache();
	// getting object from stored in memory	
	if( !gCCDParams.m_bPutScaledSamples ){
		// take only sample of required magnitude - assert if not found 
		CMagSamples* pMagSamples = m_pObjectSamples->GetMagSamples( szMagnitude, &m_LastPos );
		Assert( pMagSamples!=NULL , "Sample of magnitude %s not found",szMagnitude);
		SampleObject = pMagSamples->GetRandomSample( &m_LastSampleIdx );
		return TRUE;
	}else{
		// scaled samples - taking any of samples and scale it to required 
		// magnitude !
		mystring szSampleMag;
		Table2D<ELEM_SAMPLE_TYPE>& sample = m_pObjectSamples->GetRandomSample( szMagnitude, &szSampleMag );
		ScaleSample( szMagnitude, sample, SampleObject, szSampleMag.c_str() );
		return TRUE; 		
	}		
	return FALSE;
}


BOOL_T CImageCreator::GetSampleObject_ByName(  Table2D<ELEM_SAMPLE_TYPE>& SampleIn,
         	                          Table2D<ELEM_SAMPLE_TYPE>& SampleOut,
      	                             const char* fname,
   	                                const char* szMagnitude,
	                                   int x, int y, int nFrame )
{
	InitSampleCache();

	CMagSamples* pMagSamples = m_pObjectSamples->GetMagSamples( szMagnitude, &m_LastPos );
	Assert( pMagSamples!=NULL , "Sample of magnitude %s not found",szMagnitude);
	BOOL_T bRet = pMagSamples->GetRandomSample_ByName( x , y, SampleIn, SampleOut, nFrame );
	return bRet;		
}


BOOL_T CImageCreator::GetSampleObject_CheckXY( Table2D<ELEM_SAMPLE_TYPE>& SampleObject,
                      				              const char* fname,
                                 			     const char* szMagnitude,
															  int x, int y,
															  int nFrame )
{
	InitSampleCache();

	CMagSamples* pMagSamples = m_pObjectSamples->GetMagSamples( szMagnitude, &m_LastPos );
	Assert( pMagSamples!=NULL , "Sample of magnitude %s not found",szMagnitude);
	BOOL_T bRet = pMagSamples->GetRandomSample_CheckXY( x , y, SampleObject, nFrame );
	return bRet;	
}


void CImageCreator::AddSampleObjectAt( CCDMatrix& Image, Table2D<ELEM_SAMPLE_TYPE>& sample_img, 
													long x, long y, 
													const char* szMag, int frame_index, 
													BOOL_T bRePut,
													int cam_idx,
													double ra_in_rad, double dec_in_rad )
{
	m_pGenEvent = Image.AddSample( x, y, sample_img, 0,0,-1,-1, !bRePut, 
											szMag, frame_index, &m_LastPutObjectX, 
											&m_LastPutObjectY );

	if( m_pGenEvent ){
		m_pGenEvent->m_AstroCoord.RA = ra_in_rad;
		m_pGenEvent->m_AstroCoord.Dec = dec_in_rad;
	}

	mystring szLog;
	szLog << gCCDParams.GetOutputDir() << "/" << "sample" << cam_idx << ".log";
	MyOFile log( szLog.c_str(), "a+" );
	log.Printf("%d %d %d\n",frame_index,x,y);

	if(!bRePut){
		if(m_pPrevPutObj)
			m_pPrevPutObj->Add( CSamplePutInfo(x, y, frame_index, 0, szMag, Image.getObsTime(TRUE) ) );
	}
}

void CImageCreator::GetRandomPosition( CCDMatrix& Image, Table2D<ELEM_SAMPLE_TYPE>& sample_img,
													LONG_T& x, LONG_T& y, BOOL_T bOnlyInAnalysiedRange )
{
	GetRandomPosition( Image, sample_img.GetXSize(), sample_img.GetYSize(),
							 x, y, bOnlyInAnalysiedRange );
}

void CImageCreator::GetRandomPosition( CCDMatrix& Image, int sample_x_size, int sample_y_size,
													LONG_T& x, LONG_T& y, BOOL_T bOnlyInAnalysiedRange )
{
	LONG_T pos_x = 0;
   LONG_T pos_y = 0;

	if(bOnlyInAnalysiedRange){
		long allowed_size_x = Image.GetXSize()-2*gCCDParams.m_nIgnoreEdge-sample_x_size;
		long allowed_size_y = Image.GetYSize()-2*gCCDParams.m_nIgnoreEdge-sample_y_size;

		pos_x = (LONG_T)( gCCDParams.m_nIgnoreEdge + (CRandom::GetRandomInteger( 0 , allowed_size_x )));
		pos_y = (LONG_T)( gCCDParams.m_nIgnoreEdge + (CRandom::GetRandomInteger( 0 , allowed_size_y )));
	}else{
		long allowed_size_x = Image.GetXSize()-sample_x_size;
		long allowed_size_y = Image.GetYSize()-sample_y_size;

		pos_x = (LONG_T)( CRandom::GetRandomInteger( 0 , allowed_size_x ) );
		pos_y = (LONG_T)( CRandom::GetRandomInteger( 0 , allowed_size_y ) );
			
	}
	x = pos_x;
   y = pos_y;
}



void CImageCreator::AddSampleObject( CCDMatrix& Image, Table2D<ELEM_SAMPLE_TYPE>& sample_img,
												 long& x, long& y, const char* szMag, 
												 BOOL_T bOnlyInAnalysiedRange, int frame_index )
{
	LONG_T pos_x = 0;
   LONG_T pos_y = 0;
		
	GetRandomPosition( Image, sample_img, pos_x, pos_y, bOnlyInAnalysiedRange );

	x = pos_x;
   y = pos_y;

	m_pGenEvent = Image.AddSample( pos_x, pos_y, sample_img, 0,0,-1,-1, TRUE, szMag, frame_index, 
						 &m_LastPutObjectX, &m_LastPutObjectY );

	if(m_pPrevPutObj)
      m_pPrevPutObj->Add( CSamplePutInfo(x, y, frame_index, 0, szMag, Image.getObsTime(TRUE), m_LastPos, m_LastSampleIdx ) );	
}


void CImageCreator::RePutPrevImages( CCDPipeline& ccd_pipeline )
{
	ccd_pipeline.m_ReGenEventsList.clear();
	if(m_pPrevPutObj){
		cCCD& currFrame = ccd_pipeline.GetCurrent();
	
		CSamplePutInfo* pPrevSample = m_pPrevPutObj->GetFirst();
		while(pPrevSample){
			Table2D<ELEM_SAMPLE_TYPE> sample_img(0,0);

// temporary changed - not taking random sample, but re-put same as originaly 
/*			if(!GetSampleObject( sample_img, NULL, (pPrevSample->m_szMag).c_str() ) ){
		      //
      		// Assert(FALSE,"Cauld not read sample file : %s",szSampleFName.c_str());
		      MYTRACE3(gCCDTrace,"Could not read sample file");
      		return;
		   }*/

			// putting same sample as previously :
			CMagSamples* pMagSamples = m_pObjectSamples->GetMagSamples( pPrevSample->m_LastPos );
			if( pMagSamples ){
				sample_img = (pMagSamples->GetSamples())[ pPrevSample->m_LastSampleIdx ];
			}



			double curr_x,curr_y;

			// change on 20060407 - in pPrevSample->m_FrameIndex there is
			// day frame counter :
			int nSteps = (ccd_pipeline.GetDayFrameCounter() - pPrevSample->m_FrameIndex);			
//			int nSteps = (ccd_pipeline.GetFrameIndex() - pPrevSample->m_FrameIndex);

			CCDConfig& camParams = ccd_pipeline.GetPipelineCfg();

			CCDDataResults::CalcStarPositionAuto( (double)pPrevSample->m_X, (double)pPrevSample->m_Y,
																pPrevSample->m_prevTime, currFrame[0].getObsTime(TRUE), nSteps,
																curr_x,curr_y, 
															 	camParams.m_RotCenterX, camParams.m_RotCenterY,
																camParams.m_FrameDX, camParams.m_FrameDY,
																camParams.m_bUseRotInAverageOfPrevN, 
																camParams.m_RotValueDAlfa,
																&(ccd_pipeline.GetCCDInfoTab()[0]),
																camParams.m_RotValueDAlfaPerSec,
																&(ccd_pipeline.m_PipelineCfg) );

			AddSampleObjectAt( currFrame[0], sample_img, (long)curr_x, (long)curr_y,
									 (pPrevSample->m_szMag).c_str(), ccd_pipeline.GetFrameIndex(), TRUE );


			int max_lap=0;
			if( currFrame[0].m_pFrameLaplace ){
				long max_x,max_y;
				int start_x = curr_x-sample_img.m_MaxX;
				int start_y = curr_y-sample_img.m_MaxY;
				max_lap = (currFrame[0].m_pFrameLaplace)->GetMaxValueAndPos( start_x, start_y,
									 	start_x+sample_img.GetXSize(),
										start_y+sample_img.GetYSize(),
										max_x,max_y );
			}			

			// saving to special file :
			mystring szOutFile = gCCDParams.GetOutputDir();
			szOutFile << "/" << ccd_pipeline.m_ReGenEventsLog.m_szFileName.c_str();

			CccdReport::SaveGenEvent( szOutFile.c_str() , (long)curr_x, (long)curr_y,
											  (pPrevSample->m_szMag).c_str(), sample_img,
											  ccd_pipeline.GetDayFrameCounter(), max_lap );
			
			pPrevSample = m_pPrevPutObj->GetNext();

			ccd_pipeline.m_ReGenEventsList.Add( my_round(curr_x), my_round(curr_y) );
		}		
	}
}

void CImageCreator::AddSampleObject( CCDMatrix& Image, const char* fname, 
												 long& x, long& y,
                                     const char* szMag, BOOL_T bOnlyInAnalysiedRange, int frame_index )
{
	InitParams();
	Table2D<ELEM_SAMPLE_TYPE> sample_img(0,0);

	if(!GetSampleObject( sample_img, fname, szMag ) ){
		//
		// Assert(FALSE,"Cauld not read sample file : %s",szSampleFName.c_str());
		MYTRACE3(gCCDTrace,"Could not read sample file");
		return;
	}
	AddSampleObject( Image, sample_img, x, y, szMag, bOnlyInAnalysiedRange, frame_index );
}

BOOL_T CImageCreator::InitSampleCache()
{
	if(!m_pObjectSamples){
		m_pObjectSamples = new CccdSamples();		
		if(!m_pObjectSamples->ReadSamples( (m_pPipeline->GetPipelineCfg()).m_szSampleDir.c_str() )){
  	      printf("Could not read samples, exiting\n");
			exit(-1);
      }
	}
	return TRUE;
}

void CImageCreator::AddSampleObject( CCDMatrix& Image, LONG_T SampleIdx )
{
	InitParams();
	LONG_T minBright,maxBright;
	minBright = atol( m_MinBright.c_str() );
	maxBright = atol( m_MaxBright.c_str() );

	LONG_T bright = CRandom::GetRandomInteger(maxBright,minBright);		
	char szMag[10];
	sprintf(szMag,"%d",bright);	

	long x,y;
	if(gCCDParams.m_bKeepSamplesInMemory){
		// using cached obiect samples from memory
	   AddSampleObject( Image, NULL, x ,y, szMag );
 	}else{
		mystring szSampleFName = GetRandomSampleFile( szMag, SampleIdx );
		if(strlen(szSampleFName.c_str())==0)
			return;
		AddSampleObject( Image , szSampleFName.c_str(), x, y, szMag );
	}
}

void CImageCreator::AddNextObject( cCCD& CCDFrame )
{
	InitParams();
	double r = CRandom::GetRandom();
	if(r<=m_ObjProb){
		// now - choose camera :
		LONG_T ccd_cnt = CCDFrame.GetCount();
		LONG_T idx = CRandom::GetRandomInteger(0,ccd_cnt);
		AddSampleObject( CCDFrame[idx], m_SampleIndex	);
		m_SampleIndex++;
	}	
}


void CImageCreator::AddObject( CCDMatrix& Image )
{
	InitParams();
	AddSampleObject( Image );	
}

void CImageCreator::AddObjectWithCoord( cCCD& CCDFrame, const char* fname,
                   				 long& x, long& y,
				                   const char* szMag, BOOL_T bOnlyInAnalysied )
{
	InitParams();
	LONG_T idx = CRandom::GetRandomInteger(0,CCDFrame.GetCount());


	AddSampleObject( CCDFrame[idx] , fname, x,y, szMag, bOnlyInAnalysied );
}

void CImageCreator::AddObject( cCCD& CCDFrame, CCDPipeline& ccd_pipeline, const char* fname, 
                               const char* szMag, BOOL_T bOnlyInAnalysied )
{
	InitParams();
	LONG_T idx = CRandom::GetRandomInteger(0,CCDFrame.GetCount());
	long x,y;

	AddSampleObject( CCDFrame[idx] , fname, x, y,
                    szMag, bOnlyInAnalysied, ccd_pipeline.GetDayFrameCounter() );
}



BOOL_T CImageCreator::SimulateSuperNew( CCDMatrix& frame, int frame_index,
                    				        const char* szMinMag, const char* szMaxMag,
												  double increase_ratio,
												  int& star_x, int& star_y,
												  double& initial_mag, double& res_mag )
{
	int xSize = frame.GetXSize();
	int ySize = frame.GetYSize();	
	int ignore_edge = 20;
	int search_size = 40;

	int x0 = ignore_edge + (int)((xSize-2*ignore_edge)*CRandom::GetRandom());
	int y0 = ignore_edge + (int)((ySize-2*ignore_edge)*CRandom::GetRandom());

	int x_min = MAX( (x0 - search_size), 0 );
	int x_max = MIN( (x0 + search_size), (xSize-1) );
	int y_min = MAX( (y0 - search_size), 0);
	int y_max = MIN( (y0 + search_size), (ySize-1) );
		
	double minMag = atof( szMinMag );
	double maxMag = atof( szMaxMag );	

	double g54min = CCDDataResults::GetG54FromMag( maxMag );
	double g54max = CCDDataResults::GetG54FromMag( minMag );

	BIG_ELEM_TYPE** p_data = frame.get_frame_laplace_fast();

	CPixelAnalyseIn in;
	CPixelAnalyseOut out;
	CPixelList list(xSize*ySize);
	
	in.pPipeline = NULL;
	in.ccd_index = 0;
	in.pPixelList = &list;
	in.xSize = xSize;
	in.ySize = ySize;
	in.p_curr_data_laplace = p_data;
	in.p_data = frame.get_data_buffer();

	BOOL_T bFound=FALSE;
	star_x=0;
	star_y=0;
	for(register int y=y_min;y<=y_max && !bFound;y++){
		for(register int x=x_min;x<=x_max && !bFound;x++){
			if( p_data[y][x]>=g54min && p_data[y][x]<=g54max ){
				// find cluster and make sure that it is star :
				in.x = x;
				in.y = y;
				in.pos = x + y*xSize;
				double x0,y0,maxNoise;
				CCD_Analyser::FindClusterAboveTresholdOpt3( in, out, x0, y0,
														maxNoise, -1, FALSE );
				int max_pos = -1;
				int max_val = CCD_Analyser::FindMaxInCluster( in.p_data, xSize,
													out.cluster, out.cluster_cnt,
													max_pos );
				star_x = (max_pos % xSize );
				star_y = (max_pos / xSize );		

				if( p_data[star_y][star_x]>=g54min && p_data[star_y][star_x]<=g54max ){
					bFound = TRUE;
				}
			}
		}
	}

	if(bFound){
		double maxNoiseLevel = CCD_Analyser::GetMaxNoiseLevelLaplace( &in );				
		int count=0;

		int min_x=100000;
		int min_y=100000;
		int max_x=-1;
		int max_y=-1;
		for( register int i=0;i<out.cluster_cnt;i++){
			int x = (out.cluster[i] % xSize);
			int y = (out.cluster[i] / xSize);
			

			if( in.p_data[out.cluster[i]] >= maxNoiseLevel ){
				in.p_data[out.cluster[i]] = (int)(in.p_data[out.cluster[i]]*increase_ratio);
				count++;
			}

			if( x<min_x )
				min_x = x;
			if( y<min_y )
				min_y = y;
			if( x>max_x )
				max_x = x;
			if( y>max_y )
				max_y = y;
		}
		if(count>0){
			frame.AddGeneratedEvent( star_x, star_y, "UNKNOWN", frame_index, eBrighten );
			m_LastPutObjectX = star_x;
			m_LastPutObjectY = star_y;
			
			frame.LaplaceSafe( min_x, min_y, max_x, max_y );
		}

		return (count>0);
	}

	return FALSE;
}

mystring CImageCreator::GetRandomSampleFile( const char* szMag, LONG_T SampleIdx /* =-1 */ )
{
	mystring szSampleFName;
   mystring szSamplesList;
	mystring szSamplesDir = gCCDParams.m_szSampleDir;
	szSamplesList << szSamplesDir << "/" << m_szListname;
	
	BOOL_T bExists =  MyFile::DoesFileExist(szSamplesList.c_str());
	if(!bExists){
		MYTRACE3(gCCDTrace,"Sample directory :  " << szSamplesDir.c_str() << " does not exist" );
		return szSampleFName;
	}
	// Assert(bExists,"Directory with sample object of magnitude %d, does not exist",bright);

	InitSampleCache();
	
	int count = (m_pObjectSamples->m_DescFile).GetCount();
	int idx = CRandom::GetRandomInteger( 0 , count );
	mystring szFileNames = (m_pObjectSamples->m_DescFile).GetDescTab()[idx].szValue;	

	CMyStrTable names;
	MyParser pars = szFileNames.c_str();
	pars.GetItems( names );
	int namesCount = names.GetCount();
	int fileIdx = CRandom::GetRandomInteger( 0 , namesCount );
	szSampleFName = names[ fileIdx ];
	return szSampleFName;		
}

void CImageCreator::AddObjects( cCCD& CCDFrame )
{
	InitParams();
	double r = CRandom::GetRandom();
	if(r<=m_ObjProb){
		// now - choose camera :
		LONG_T ccd_cnt = CCDFrame.GetCount();
		LONG_T idx = CRandom::GetRandomInteger(0,ccd_cnt);
		AddObject( CCDFrame[idx]	);
	}
}


double CImageCreator::get_gauss( double s0, double sigma )
{
	double sum = 0;
	for(int i=0;i<12;i++){
		sum += CRandom::GetRandom();
	}
	sum -= 6.00;
	sum = s0+sum*sigma;
	return sum;
	
}

BOOL_T CImageCreator::IsGenEvent( CPoint& evt )
{
	if( m_LastPutObjectX>0 && m_LastPutObjectY>0 ){
		if ( fabs(m_LastPutObjectX-evt.x)<=gCCDParams.m_bGenEventRedial 
				&& fabs(m_LastPutObjectY-evt.y)<=gCCDParams.m_bGenEventRedial ){
			return TRUE;
		}
	}
	return FALSE;
}


void CImageCreator::LogSamplePut( Table2D<ELEM_SAMPLE_TYPE>* sample, int cam, int run,
											 int frame, int x, int y, const char* szError )
{
	if( gCCDParams.m_bLogSamplePut ){
		char szTmp[16];
		sprintf(szTmp,"run%.5d",run);
	
		mystring szFile;
		szFile << gCCDParams.GetOutputDir() << "/Samples/" << szTmp
				 << "/samples_" << cam << ".log";

		if( !MyFile::DoesFileExist( szFile.c_str() ) ){
			MyOFile out(szFile.c_str(),"a+");
			out.Printf("Frame# MAG Name X Y\n");
		}

		MyOFile out(szFile.c_str(),"a+");
		

		if( !szError ){
			if( sample ){
				const char* szName = (sample->GetKeyTab()).getKeyVal( PARTNAME );
				const char* szMag = (sample->GetKeyTab()).getKeyVal( PARTMAG );
				if( szName && szName[0] && szMag && szMag[0] ){
					out.Printf( "%d %s %s %d %d\n",frame,szMag,szName,x,y);
				}
			}
		}else{
			if ( szError[0] ){
				out.Printf( "%d ERROR:%s %d %d\n",frame,szError,x,y);
			}
		}
	}
}

