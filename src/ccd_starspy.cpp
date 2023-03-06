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
#include "ccd_starspy.h"
#include "ccd_globals.h"
#include "ccd_trace.h"
#include "ccd_analyse.h"
#include <cexcp.h>
#include <baseanal.h>
#include <math.h>
#include <calcrot.h>
#include <myfile.h>
#include <mystrtable.h>
#include <myparser.h>
#include <mymacros.h>
#include <myfits.h>
#include <mathfunc.h>
#include <cfg.h>
#include <myutil.h>
#include "ccd_pipeline.h"
#include "ccd_fastcalc.h"

CStarDesc::CStarDesc( int low_x, int low_y, int up_x, int up_y )
: m_StartX(0), m_StartY(0), m_CurrX(0), m_CurrY(0), m_PrevStepX(0), m_PrevStepY(0), m_AverageStepX(0), 
  m_AverageStepY(0), m_PrevMAX(0),
  m_LowX(low_x),m_LowY(low_y),m_UpX(up_x),m_UpY(up_y),m_OrtoLineA(0),
  m_OrtoLineB(0), m_OrtoLineC(0), m_starDX(0), m_starDY(0), 
  m_starDXPerSec(0), m_starDYPerSec(0),m_bNotEnd(TRUE), m_bGood(TRUE),
  m_nCurrentStep(0), m_MaxDX(0.00), m_MaxDY(0.00), m_MeanValue(0), m_RMS(0),
  m_StarValues("Star"), m_nDrasticChangeCount(0), m_bIdentOK(TRUE),
  m_PredictedX(0), m_PredictedY(0), m_StarRA(0), m_StarDec(0)
{
	m_PrevStepX = gCCDParams.m_FrameDX;
	m_PrevStepY = gCCDParams.m_FrameDY;
}


CStarDesc& CStarDesc::operator=(const CStarDesc& right)
{
	m_StarRA = right.m_StarRA;
	m_StarDec = right.m_StarDec;
	m_LowX = right.m_LowX;
	m_LowY = right.m_LowY;
	m_UpX = right.m_UpX;
	m_UpY = right.m_UpY;	
	m_StartX = right.m_StartX;
	m_StartY = right.m_StartY;	
	m_CurrX = right.m_CurrX;
	m_CurrY = right.m_CurrY;
	m_PrevStepX = right.m_PrevStepX;
	m_PrevStepY = right.m_PrevStepY;
	m_AverageStepX = right.m_AverageStepX;	
	m_AverageStepY = right.m_AverageStepY;
	m_PrevMAX = right.m_PrevMAX;
	m_PosTab = right.m_PosTab;
	m_bNotEnd = right.m_bNotEnd;	
   m_starDX = right.m_starDX;
   m_starDY = right.m_starDY;
   m_starDXPerSec = right.m_starDXPerSec;
   m_starDYPerSec = right.m_starDYPerSec;
	m_OrtoLineA = right.m_OrtoLineA;
	m_OrtoLineB = right.m_OrtoLineB;
	m_OrtoLineC = right.m_OrtoLineC;
	m_nCurrentStep = right.m_nCurrentStep;
	m_MaxDX = right.m_MaxDX;
	m_MaxDY = right.m_MaxDY;
	m_MeanValue = right.m_MeanValue;
	m_RMS = right.m_RMS;
	m_PredictedX = right.m_PredictedX;
	m_PredictedY = right.m_PredictedY;

	return (*this);
}


double CStarDesc::UpdateRMS()
{
	double rms =0;
	double mean = m_MeanValue/m_PosTab.size();
	for(register int i=0;i<m_PosTab.size();i++){
		rms += CMyMathFunc::mysqr(m_PosTab[i].m_Value-mean);
	}
	rms = sqrt(rms/m_PosTab.size());
	m_RMS = rms;
	return m_RMS;
}

CPoint* CStarDesc::FindFrame( int frame_index )
{
	for(int i=0;i<m_PosTab.size();i++){
		if(m_PosTab[i].m_FrameIndex == frame_index ){
			return &(m_PosTab[i]);
		}
	}
	return NULL;
}

int CStarDesc::GetAverageShifts( double& aver_dx, double& aver_dy, int back )
{
	aver_dx = 0;
	aver_dy = 0;
	int cnt=0;
	int downto=(m_PosTab.size()-back);
	for(int i=(m_PosTab.size()-1);(i>=1 && i>=downto);i--){
		aver_dx += ( m_PosTab[i].x - m_PosTab[i-1].x );
		aver_dy += ( m_PosTab[i].y - m_PosTab[i-1].y );
		cnt++;
	}
	aver_dx = ( aver_dx / cnt );
	aver_dy = ( aver_dy / cnt );

	return cnt;
}

CCDStarSpy::CCDStarSpy( int sizeX, int sizeY, int nSize, CCDPipeline* pPipeline )
: m_SizeX(sizeX), m_SizeY(sizeY), m_StarsCount(nSize*nSize), m_nSize(nSize),
  m_RotCenterX(0), m_RotCenterY(0), m_dAlfa(0), m_DX(0), m_DY(0), m_bRotation(FALSE),
  m_DXPerSec(0), m_DYPerSec(0), m_bFromFile(FALSE),
  m_StartTime(0), m_EndTime(0), m_dAlfaPerSec(0), m_pPipeline(pPipeline),
	m_nFramesCount(0)
{
	if(m_SizeX==0){
		m_SizeX = gCCDParams.m_SizeX;
	}
	if(m_SizeY==0){
		m_SizeY = gCCDParams.m_SizeY;
	}

//	Assert( (m_nSize%2)==0 || m_nSize==1,"Only even sizes and 1 are allowed for star tracer );
	m_StarDesc = new CStarDesc[ m_StarsCount ];

	InitAreas( 50 );
}

void CCDStarSpy::InitAreas( int BorderSize )
{
	int PixelSizeX = (m_SizeX-2*BorderSize)/m_nSize;
	int PixelSizeY = (m_SizeY-2*BorderSize)/m_nSize;

	int i=0;
	for(int y=0;y<m_nSize;y++){
		for(int x=0;x<m_nSize;x++){
			m_StarDesc[i].m_LowX = x*PixelSizeX + BorderSize;
			m_StarDesc[i].m_LowY = y*PixelSizeY + BorderSize;
			m_StarDesc[i].m_UpX = (x+1)*PixelSizeX + BorderSize;
			m_StarDesc[i].m_UpY = (y+1)*PixelSizeY + BorderSize;
			i++;						
		}
	}

}


void CCDStarSpy::ResetTracedStar()
{
	for(int i=0;i<m_StarsCount;i++){
		m_StarDesc[i].m_PosTab.clear();
		m_StarDesc[i].m_bNotEnd=FALSE;
	}
}

CStarDesc& CCDStarSpy::GetStarDesc( int x , int y )
{
	Assert( x>=0 && x<m_nSize && y>=0 && y<m_nSize,"Number of columns or row exceeded");
	int pos = y*m_nSize+x;
	Assert( pos>=0 && pos<m_StarsCount,"Not in star numer range");
	return m_StarDesc[pos];
}


CCDStarSpy::CCDStarSpy( int sizeX, int sizeY, double rotX, double rotY, double dAlfa, double dx, 
								double dy, BOOL_T bRot, double dx_per_sec, double dy_per_sec,
								double dAlfaPerSec, CCDPipeline* pPipeline  )
: m_SizeX(sizeX), m_SizeY(sizeY), m_bFromFile(TRUE), m_RotCenterX(rotX),
  m_RotCenterY(rotY), m_dAlfa(dAlfa), m_DX(dx), m_DY(dy), m_DXPerSec(dx_per_sec),
  m_DYPerSec( dy_per_sec ), m_bRotation( bRot ), m_StarsCount(1), m_nSize(1),
  m_dAlfaPerSec( dAlfaPerSec ), m_pPipeline(pPipeline),
  m_StartTime(0), m_EndTime(0), m_nFramesCount(0)
{
	m_StarDesc = new CStarDesc[1];
	InitAreas( 50 );
}


CCDStarSpy::CCDStarSpy( const CCDStarSpy& right )
: m_StarDesc(NULL)
{
	(*this) = right;
}

CCDStarSpy& CCDStarSpy::operator=(const CCDStarSpy& right )
{
	if(m_StarDesc)
		delete [] m_StarDesc;
	
   m_StarsCount = right.m_StarsCount;
   m_nSize = right.m_nSize;
	m_StarDesc = new CStarDesc[m_StarsCount];
	for(int i=0;i<m_StarsCount;i++){
		m_StarDesc[i] = right.m_StarDesc[i];
	}

   m_SizeX = right.m_SizeX;
   m_SizeY = right.m_SizeY;

   m_StartTime = right.m_StartTime;
   m_EndTime = right.m_EndTime;
   m_bFromFile = right.m_bFromFile;
   m_bRotation = right.m_bRotation;
   m_RotCenterX = right.m_RotCenterX;
   m_RotCenterY = right.m_RotCenterY;
   m_dAlfa = right.m_dAlfa;
	m_dAlfaPerSec = right.m_dAlfaPerSec;
   m_DX = right.m_DX;
   m_DY = right.m_DY;
   m_DXPerSec = right.m_DXPerSec;	
   m_DYPerSec = right.m_DYPerSec;	
	m_nFramesCount = right.m_nFramesCount;
	m_pPipeline = right.m_pPipeline;

	return (*this);
}

CCDStarSpy::~CCDStarSpy()
{
	if(m_StarDesc){
		delete [] m_StarDesc;
	}
}

int CCDStarSpy::GetFramesCount()
{
	int nCount=0;
	for(register int i=0;i<m_StarsCount;i++){
		if(m_StarDesc[i].m_PosTab.size()>nCount){
			nCount = m_StarDesc[i].m_PosTab.size();
		}
	}
	return nCount;
}

void CCDStarSpy::InitWithSpecificStar( ELEM_TYPE** p_data, double startTime, int starX0, int starY0, int frame_index )
{
	m_StartTime = startTime;

	m_nFramesCount = 0;
	for(register int i=0;i<m_StarsCount;i++){		
		int x_max = starX0;
		int y_max = starY0;
	
		m_StarDesc[i].m_MeanValue = 0;
		m_StarDesc[i].m_RMS = 0;
		m_StarDesc[i].m_MaxDX = 0.00;
		m_StarDesc[i].m_MaxDY = 0.00;
		m_StarDesc[i].m_nCurrentStep = 0;
		m_StarDesc[i].m_StartX = x_max;
		m_StarDesc[i].m_StartY = y_max;
		m_StarDesc[i].m_CurrX = m_StarDesc[i].m_StartX;	
		m_StarDesc[i].m_CurrY = m_StarDesc[i].m_StartY;
		m_StarDesc[i].m_PrevStepX = gCCDParams.m_FrameDX;
		m_StarDesc[i].m_PrevStepY = gCCDParams.m_FrameDY;
		m_StarDesc[i].m_AverageStepX = 0;
		m_StarDesc[i].m_AverageStepY = 0;
		m_StarDesc[i].m_PrevMAX = p_data[y_max][x_max];
		m_StarDesc[i].m_bIdentOK = TRUE;
		// m_StarValues.Init( 0, 30000, 200 );

		m_StarDesc[i].m_PosTab.Clear();
		CPoint tmp( x_max, y_max, frame_index, (int)startTime, m_StarDesc[i].m_PrevMAX );
		m_StarDesc[i].m_PosTab.Add( tmp );
	}

}

void CCDStarSpy::ReInitMAXStarIfNeeded( ELEM_TYPE** p_data, double startTime, int frame_index, BOOL_T bForceReInit/*=FALSE*/ )
{
	for(int i=0;i<m_StarsCount;i++){
		if( ( !m_StarDesc[i].m_bNotEnd && (frame_index-m_StarDesc[i].m_PosTab.back().m_FrameIndex)>=1 ) || bForceReInit ){
			// for this star re-initialization with new star is required
			if(m_StarDesc[i].InitMAXStar( p_data, startTime, frame_index, m_SizeX, m_SizeY, m_StarDesc, m_StarsCount, m_pPipeline )){
				m_StartTime = startTime;
				printf("\n\nINFO FRAME %d : !!!!! Star %d re-initialized, new position : (%d,%d)\n",frame_index,i,(int)m_StarDesc[i].m_CurrX,(int)m_StarDesc[i].m_CurrY);
			}else{
				printf("\n\nWARNING FRAME %d !!!!! Could not re-initialize star %d\n",frame_index,i);
			}
		}
	}
}

void CCDStarSpyTab::MarkToReInit()
{
	vector<CCDStarSpy>::iterator i;

   for(i=begin();i!=end();i++){
		for(int j=0;j<i->m_StarsCount;j++){
			((i->m_StarDesc)[j]).m_bNotEnd = FALSE;
		}
	}
}

BOOL_T CStarDesc::IsDistFromOthersOK( int x, int y, CStarDesc* pOtherStarsTab, int otherCount )
{
	for(register int i=0;i<otherCount;i++){
		if( fabs(x-pOtherStarsTab[i].m_CurrX)<gCCDParams.m_MinDistBetweenTracedStars &&
			 fabs(y-pOtherStarsTab[i].m_CurrY)<gCCDParams.m_MinDistBetweenTracedStars ){
			return FALSE;
		}
	}
	return TRUE;
}

BOOL_T CStarDesc::InitMAXStar( ELEM_TYPE** p_data, double startTime, int frame_index, int SizeX,
										 int SizeY, CStarDesc* pOtherStarsTab, int otherCount,
										 CCDPipeline* pPipeline ) 
{
		int x_max=-1;
		int y_max=-1;
		double max_val = -100000;

		CPixelAnalyseOut* pOut = new CPixelAnalyseOut();
		BOOL_T bFound=FALSE;

		int retry=10;
		while( !bFound && retry>0){
			x_max=-1;
			y_max=-1;
			max_val = -100000;
			for(register int x=m_LowX;x<=m_UpX;x++){
				for(register int y=m_LowY;y<=m_UpY;y++){
					if( x>=40 && y>=40 && x<(SizeX-40) && y<(SizeY-40)){
						double val = Table2D<ELEM_TYPE>::CalcG54(x,y,SizeX,p_data);
						if( val > max_val && ( !pOtherStarsTab || IsDistFromOthersOK( x, y, pOtherStarsTab,otherCount)) ){
							int pos = SizeX*y + x;
							if( find_value( pOut->hitpixel_list, pOut->hitpixel_count, pos )<0 ){
								max_val = val;
								x_max = x;
								y_max = y;
							}
						}
					}
				}
			} 	
			
			bFound=TRUE;

			double cluster_treshold = CCD_Analyser::GetMaxNoiseLevelLaplace( pPipeline, 0, gCCDParams.m_MAXStarNSigmaAbove );

			int count = CCD_Analyser::FindClusterAboveTresholdOpt3( x_max, y_max, SizeX, SizeY,
								p_data, (*pOut), pOut->cluster, pOut->cluster_cnt, cluster_treshold );
			if(pOut->cluster_cnt>0 && gCCDParams.m_MAXStarShapeLimit>0){
				(pOut->m_PixelOut.max_redial) = CCDFastCalc::FindMaxRedial( pOut->cluster, pOut->cluster_cnt, 
																			  (pOut->m_PixelOut).x0, (pOut->m_PixelOut).y0, SizeX );
				(pOut->m_PixelOut).pos0 = (pOut->m_PixelOut).y0*SizeX + (pOut->m_PixelOut).x0;
				(pOut->m_PixelOut).m_Sphericity = CCDFastCalc::CalcSphericity( (pOut->m_PixelOut).max_redial, pOut->cluster_cnt );
				if( (pOut->m_PixelOut).m_Sphericity< gCCDParams.m_MAXStarShapeLimit){
					printf("MAX star at (%d,%d) rejected, due to shape = %.3f , below limit of %.3f\n",x_max,y_max,(pOut->m_PixelOut).m_Sphericity,gCCDParams.m_MAXStarShapeLimit);
					bFound=FALSE;
				}else{
					_TRACE_PRINTF_0("Initialized star at (%d,%d)\n",x_max,y_max);
				}
			}													
			retry--;
		}
		delete pOut;

		if(y_max<0)
			return FALSE;
	
		m_MeanValue = max_val;
		m_RMS = 0;
		m_MaxDX = 0.00;
      m_MaxDY = 0.00;
		m_StartX = x_max;
		m_StartY = y_max;
		m_CurrX = m_StartX;	
		m_CurrY = m_StartY;
		m_PrevStepX = gCCDParams.m_FrameDX;
		m_PrevStepY = gCCDParams.m_FrameDY;
		m_AverageStepX = 0;
		m_AverageStepY = 0;
		m_PrevMAX = max_val;
		m_bIdentOK = TRUE;
		m_nCurrentStep=0;


		m_PosTab.Clear();
		CPoint tmp( x_max, y_max, frame_index, (int)startTime, m_PrevMAX );
		m_PosTab.Add( tmp );
		m_bNotEnd = TRUE;


		return TRUE;
}

void CCDStarSpy::InitWithMAXStar( ELEM_TYPE** p_data, double startTime, 
											 int frame_index )
{
	m_StartTime = startTime;

	m_nFramesCount = 0;
	for(int i=0;i<m_StarsCount;i++){		
		m_StarDesc[i].m_nCurrentStep = 0;

		int bottomX = m_StarDesc[i].m_LowX;
		int bottomY = m_StarDesc[i].m_LowY;
		int upX = m_StarDesc[i].m_UpX;
		int upY = m_StarDesc[i].m_UpY;

		if( !m_StarDesc[i].InitMAXStar( p_data, startTime, frame_index,
												  m_SizeX, m_SizeY, NULL, 0, m_pPipeline ) ){
			Assert(FALSE,"No positive value found ???");
		}


		/*int x_max=-1;
		int y_max=-1;
		double max_val = -1000;
		for(register int x=bottomX;x<=upX;x++){
			for(register int y=bottomY;y<=upY;y++){
				if( x>=2 && y>=2 && x<(m_SizeX-2) && y<(m_SizeY-2)){
					double val = Table2D<ELEM_TYPE>::CalcG54(x,y,m_SizeX,p_data);
					if( val > max_val ){
						max_val = val;
						x_max = x;
						y_max = y;
					}
				}
			}
		} 	
		Assert(y_max>=0,"No positive value found ???");
	
		m_StarDesc[i].m_MeanValue = max_val;
		m_StarDesc[i].m_RMS = 0;
		m_StarDesc[i].m_MaxDX = 0.00;
      m_StarDesc[i].m_MaxDY = 0.00;
		m_StarDesc[i].m_StartX = x_max;
		m_StarDesc[i].m_StartY = y_max;
		m_StarDesc[i].m_CurrX = m_StarDesc[i].m_StartX;	
		m_StarDesc[i].m_CurrY = m_StarDesc[i].m_StartY;
		m_StarDesc[i].m_PrevStepX = gCCDParams.m_FrameDX;
		m_StarDesc[i].m_PrevStepY = gCCDParams.m_FrameDY;
		m_StarDesc[i].m_AverageStepX = 0;
		m_StarDesc[i].m_AverageStepY = 0;
		m_StarDesc[i].m_PrevMAX = max_val;


		m_StarDesc[i].m_PosTab.Clear();
		CPoint tmp( x_max, y_max, frame_index, (int)startTime, m_StarDesc[i].m_PrevMAX );
		m_StarDesc[i].m_PosTab.Add( tmp );*/
	}
}


double CStarDesc::CalcAverageDX( int new_x , int nPrev )
{
	int startPos = m_PosTab.size()-nPrev;
	if(startPos<0)
		startPos = 0;
	int nFrames = (m_PosTab.size()-startPos); // -1 + 1 = 0 
	double dx = (new_x - m_PosTab[startPos].x)/nFrames;
	return dx;
}

double CStarDesc::CalcAverageDY( int new_y , int nPrev )
{
	int startPos = m_PosTab.size()-nPrev;
	if(startPos<0)
		startPos = 0;
	int nFrames = (m_PosTab.size()-startPos); // -1 + 1 = 0 
	double dy = (new_y - m_PosTab[startPos].y)/nFrames;
	return dy;
}



BOOL_T CCDStarSpy::ReCalcNewPositionByRaDec( int frame_time, CCDPipeline* pPipeline )
{
	return FALSE;
}

void CCDStarSpy::LogMAXStarsPositions( CCDPipeline* pPipeline )
{
	mystring szFileName,szBase;
	szBase << "startracks_" << pPipeline->GetPipelineIndex() << ".log";
	CCDPipeline::GetOutFileName( szFileName, "StarTrackLog", szBase.c_str(),
										  pPipeline, -1, FALSE );

	if( !MyOFile::DoesFileExist( szFileName.c_str() ) ){
		MyOFile out( szBase.c_str() );
		out.Printf("%s","# Frame StartIdx x y\n");
	}

	MyOFile out( szFileName.c_str() , "a" );
	for(int i=0;i<m_StarsCount;i++){
		out.Printf("%d %d %d %d\n",pPipeline->GetFrameIndex(),i,
						(int)(m_StarDesc[i].m_CurrX),(int)(m_StarDesc[i].m_CurrY));
	}
	out.Printf("\n");
}

BOOL_T CCDStarSpy::ReCalcNewPosition( ELEM_TYPE** p_data, int frame_time, int frame_index,
												  int searchSize, double sizeTolerance, BOOL_T bUseHint,
												  int reCalcSearchSize,
												  BOOL_T bUseSteps, double dx, double dy,
												  BOOL_T* reInit )
{
	if( reInit ){
		(*reInit) = FALSE;
	}	

	BOOL_T bRet = FALSE;
	int nNotInTolerance=0;


//	if(frame_index>=227 || frame_index==32){
//		printf("odo");
//	}

	for(int i=0;i<m_StarsCount;i++){
		if(!m_StarDesc[i].m_bNotEnd)
			continue;

		m_StarDesc[i].m_bIdentOK = TRUE;

		register int x_start = (int)MAX((m_StarDesc[i].m_CurrX-searchSize),0);
		register int y_start = (int)MAX((m_StarDesc[i].m_CurrY-searchSize),0);
		register int x_end = (int)MIN((m_StarDesc[i].m_CurrX+searchSize),(m_SizeX-1));
		register int y_end = (int)MIN((m_StarDesc[i].m_CurrY+searchSize),(m_SizeY-1));

		int x_max=-1;
	   int y_max=-1;
   	double max_val = -1000;


		register int predictionX = my_round(m_StarDesc[i].m_CurrX + m_StarDesc[i].m_AverageStepX);
      register int predictionY = my_round(m_StarDesc[i].m_CurrY + m_StarDesc[i].m_AverageStepY);		

		if( bUseHint ){	
			if( !bUseSteps ){
				if( m_StarDesc[i].m_nCurrentStep==0 ){
					predictionX = my_round(m_StarDesc[i].m_CurrX + gCCDParams.m_FrameDX );
					predictionY = my_round(m_StarDesc[i].m_CurrY + gCCDParams.m_FrameDY );
				}	
			}else{
				predictionX = my_round( m_StarDesc[i].m_CurrX + dx );
				predictionY = my_round( m_StarDesc[i].m_CurrY + dy );
			}
		}

		if(m_StarDesc[i].m_nCurrentStep==0 && !bUseHint){
				for(register int x=x_start;x<=x_end && m_StarDesc[i].m_bNotEnd;x++){
					for(register int y=y_start;y<=y_end && m_StarDesc[i].m_bNotEnd;y++){
						if( x>=2 && y>=2 && x<(m_SizeX-2) && y<(m_SizeY-2)){
							double val = Table2D<ELEM_TYPE>::CalcG54(x,y,m_SizeX,p_data);
							if( val > max_val ){
								max_val = val;
								x_max = x;
								y_max = y;
							}			
						}else{
							m_StarDesc[i].m_bNotEnd = FALSE;
						}
					}
				}
		}else{

			for(register int x=(predictionX-reCalcSearchSize);x<(predictionX+reCalcSearchSize) && m_StarDesc[i].m_bNotEnd;x++){
				for(register int y=(predictionY-reCalcSearchSize);y<(predictionY+reCalcSearchSize) && m_StarDesc[i].m_bNotEnd;y++){
					if( x>=2 && y>=2 && x<(m_SizeX-2) && y<(m_SizeY-2)){
						double val = Table2D<ELEM_TYPE>::CalcG54(x,y,m_SizeX,p_data);
						if( val > max_val ){
							max_val = val;
							x_max = x;
							y_max = y;
						}
					}else{
						m_StarDesc[i].m_bNotEnd = FALSE;
					}
				}
			}			
		}

		if(!m_StarDesc[i].m_bNotEnd){
			// during re-calculation of new positin could turned out that already at the edge :
      	continue;
		}
	
		BOOL_T bIsStarOK=TRUE;		
		// fabs((max_val-m_StarDesc[i].m_PrevMAX)/m_StarDesc[i].m_PrevMAX) < sizeTolerance

		if( bIsStarOK ){
			double newStepX = (x_max - m_StarDesc[i].m_CurrX);
			double newStepY = (y_max - m_StarDesc[i].m_CurrY);

			

			
			eStarIdentType_T identType=eStarIdent_OK;
			if( m_StarDesc[i].m_PosTab.size()>10 ){ // was >30 VERIFY !!!!
				double aver_dx,aver_dy;
				int cnt = m_StarDesc[i].GetAverageShifts( aver_dx, aver_dy , 20 );

				double new_aver_dx = ( (aver_dx*cnt) + newStepX ) / (cnt+1);
				double new_aver_dy = ( (aver_dy*cnt) + newStepY ) / (cnt+1);
				// double new_aver_dx = newStepX;
				// double new_aver_dy = newStepY;
				double diff_dx = ( new_aver_dx - aver_dx );
				double diff_dy = ( new_aver_dy - aver_dy );
				double mean = ( m_StarDesc[i].m_MeanValue/m_StarDesc[i].m_PosTab.size() );
				if(mean==0)
					mean = max_val;
				double diff_signal = fabs( (max_val - mean )/mean );
				
				if( ( (fabs(x_max - predictionX)>2.00*fabs(m_StarDesc[i].m_MaxDX) && m_StarDesc[i].m_MaxDX!=0) || (m_StarDesc[i].m_MaxDX==0 && fabs(x_max - predictionX)>1.00) ) ||
					 ( (fabs(y_max - predictionY)>2.00*fabs(m_StarDesc[i].m_MaxDY) && m_StarDesc[i].m_MaxDY!=0) || (m_StarDesc[i].m_MaxDY==0 && fabs(y_max - predictionY)>1.00) ) ||
					 ( (aver_dx!=0 && fabs(diff_dx/aver_dx)>1.0) || (aver_dx==0 && fabs(diff_dx)>1.0 ) ) ||
					 ( (aver_dy!=0 && fabs(diff_dy/aver_dy)>1.0) || (aver_dy==0 && fabs(diff_dy)>1.0 ) ) ||
					 (diff_signal>gCCDParams.m_MaxDiffSignalInSigma)){
					 // (diff_signal>gCCDParams.m_MaxDiffSignalInSigma*m_StarDesc[i].m_RMS)){
					// in case max_dx is exceeeded by factor 2 or 
					// new_average would exceed average from previous by factor 100%
					// claim problem :
					// checking if change not too drastic it would mean something 
					// bad with frame occured :
					// in such case we assume predicted position and wait for next frame :
					printf("WARNING FRAME:%d Drastic change of star:%d movemement ????\n",frame_index,i);
					printf("WARNING FRAME:%d Average shift was (%.2f , %.2f) , now we have (%.2f , %.2f)\n",
							 frame_index,(double)m_StarDesc[i].m_AverageStepX,(double)m_StarDesc[i].m_AverageStepY,(double)newStepX,(double)newStepY);
					printf("WARNING FRAME:%d Ignoring new position (%d,%d), re-using averages\n",frame_index,x_max,y_max);
					printf("WARNING FRAME:%d max shifts : (%.2f,%.2f)\n",frame_index,m_StarDesc[i].m_MaxDX,m_StarDesc[i].m_MaxDY);
					if(diff_signal>gCCDParams.m_MaxDiffSignalInSigma){
						printf("WARNING FRAME:%d BIG value change mean=%.5f new_value=%.5f\n",frame_index,mean,max_val);
					}
	
					if(m_StarDesc[i].m_PosTab.size()>gCCDParams.m_nAutoShiftsMinStepsToIgnore){
						printf("\n\nMINIMUM NUMBER OF STEPS : %d REACHED\n",gCCDParams.m_nAutoShiftsMinStepsToIgnore);
						printf("CURRENT STEPS NUMBER = %d\n",m_StarDesc[i].m_PosTab.size());
						m_StarDesc[i].m_bNotEnd = FALSE;
						continue;
					}

					newStepX = m_StarDesc[i].m_AverageStepX;
					newStepY = m_StarDesc[i].m_AverageStepY;
					x_max = predictionX;
					y_max = predictionY;
					identType=eStarIdent_ReUsedOld;
					m_StarDesc[i].m_nDrasticChangeCount++;
					m_StarDesc[i].m_bIdentOK = FALSE;
				}
			}


			if(fabs(newStepX)>fabs(m_StarDesc[i].m_MaxDX)){
				m_StarDesc[i].m_MaxDX = newStepX;
			}
			if(fabs(newStepY)>fabs(m_StarDesc[i].m_MaxDY)){
				m_StarDesc[i].m_MaxDY = newStepY;
			}
			m_StarDesc[i].m_PrevStepX = newStepX;
			m_StarDesc[i].m_PrevStepY = newStepY;
			m_StarDesc[i].m_AverageStepX = (m_StarDesc[i].m_AverageStepX*m_StarDesc[i].m_nCurrentStep + m_StarDesc[i].m_PrevStepX)/(m_StarDesc[i].m_nCurrentStep+1);
			m_StarDesc[i].m_AverageStepY = (m_StarDesc[i].m_AverageStepY*m_StarDesc[i].m_nCurrentStep + m_StarDesc[i].m_PrevStepY)/(m_StarDesc[i].m_nCurrentStep+1);
		

			m_StarDesc[i].m_CurrX = x_max;
			m_StarDesc[i].m_CurrY = y_max;
			m_StarDesc[i].m_nCurrentStep++;				
			m_nFramesCount++;


		   CPoint tmp( x_max, y_max, frame_index, frame_time, max_val, identType );
   		m_StarDesc[i].m_PosTab.Add( tmp );

			m_StarDesc[i].m_bNotEnd = FALSE;			
			if( x_max>10 && x_max<(m_SizeX-10) && y_max>10 && y_max<(m_SizeY-10) ){
				m_StarDesc[i].m_bNotEnd = TRUE;
				bRet = TRUE;
			}else{	
				_TRACE_PRINTF_3("INFO : Star : %d (%.2f,%.2f) already at edge : (%.2f,%.2f)\n",i,m_StarDesc[i].m_StartX,m_StarDesc[i].m_StartY,m_StarDesc[i].m_CurrX,m_StarDesc[i].m_CurrY);
				m_StarDesc[i].m_bNotEnd = FALSE;
			}
			m_StarDesc[i].m_MeanValue += max_val;
			m_StarDesc[i].UpdateRMS();
			_TRACE_PRINTF_5("New star position : (%.5f,%.5f)\n",m_StarDesc[i].m_CurrX,m_StarDesc[i].m_CurrY);
		}
/*		else{
			nNotInTolerance++;
			printf("Star %d not in tolerance range\n",i);
			printf("new value = %f , prev value = %f\n",max_val,m_StarDesc[i].m_PrevMAX);
			double r = fabs((max_val-m_StarDesc[i].m_PrevMAX)/m_StarDesc[i].m_PrevMAX);
			printf("ratio r = %f > %f\n",r,sizeTolerance );
		}*/
	}


	int nBad=0;
	for(int i=0;i<m_StarsCount;i++){
		if(!m_StarDesc[i].m_bIdentOK)
			nBad++;
	}
	if(nBad>=(m_StarsCount-1) && nBad>0){
		// almost all had Drastic change , re-intializing all !!! :
		ReInitMAXStarIfNeeded( p_data , frame_time, frame_index, TRUE );
		gCCDParams.m_nSkipNFramesAfterChange = gCCDParams.m_nMaxOfAverageOfPrevN;
		if( reInit ){
			(*reInit) = TRUE;
		}
	}

	return bRet;
}


BOOL_T CCDStarSpy::CalcRotation( double& angle_per_frame,  double& x_cross, double& y_cross,
											double& frameDX, double& frameDY, 
											double& frameDXPerSec, double& frameDYPerSec, 
											BOOL_T& bRotation, double& angle_per_sec, double endTime,
											int min_steps/* =5 */, int min_steps_for_single/* =20 */ )
{
	if(min_steps<1)
		min_steps = 1;
	
	x_cross = 0;
   y_cross = 0;
   angle_per_frame = 0;
	angle_per_sec = 0;
   frameDXPerSec = 0;
   frameDYPerSec = 0;


   m_EndTime = endTime;
	if(m_nFramesCount==0)
	   return FALSE;



	if(m_StarsCount==1){
		int half = (m_StarDesc->m_PosTab.size()/2);
   	int last = m_StarDesc->m_PosTab.size()-1;


		double totalShift = sqrt( ( m_StarDesc->m_PosTab[last].x - m_StarDesc->m_PosTab[0].x )*( m_StarDesc->m_PosTab[last].x - m_StarDesc->m_PosTab[0].x )
										  + ( m_StarDesc->m_PosTab[last].y - m_StarDesc->m_PosTab[0].y )*( m_StarDesc->m_PosTab[last].y - m_StarDesc->m_PosTab[0].y ) );

		frameDX = ( m_StarDesc->m_PosTab[last].x - m_StarDesc->m_PosTab[0].x )/m_StarDesc->m_nCurrentStep;
		frameDY = ( m_StarDesc->m_PosTab[last].y - m_StarDesc->m_PosTab[0].y )/m_StarDesc->m_nCurrentStep;


		double totalTime = 0;
		if(endTime>m_StartTime){
			totalTime = (endTime-m_StartTime);
			if(totalTime!=0){
				frameDXPerSec = ( m_StarDesc->m_PosTab[last].x - m_StarDesc->m_PosTab[0].x )/totalTime;
				frameDYPerSec = ( m_StarDesc->m_PosTab[last].y - m_StarDesc->m_PosTab[0].y )/totalTime;
			}
		}
	
		m_StarDesc->m_starDX = m_DX = frameDX;
      m_StarDesc->m_starDY = m_DY = frameDY;
   	m_StarDesc->m_starDXPerSec = m_DXPerSec = frameDXPerSec;
	   m_StarDesc->m_starDYPerSec = m_DYPerSec = frameDYPerSec;
					


		bRotation = FALSE;
		m_bRotation = FALSE;
	
		if(m_StarDesc->m_PosTab.size()>=min_steps_for_single ){	
			if( (fabs(frameDX)>gCCDParams.m_MinShiftToCalcRot || fabs(frameDY)>gCCDParams.m_MinShiftToCalcRot) ||
				 (totalShift > gCCDParams.m_MinTotalShiftToCalcRot) ){

				CPoint& beg1 = m_StarDesc->m_PosTab[0];
		      CPoint& end1 = m_StarDesc->m_PosTab[half];
	
				CPoint& beg2 = m_StarDesc->m_PosTab[half];
      		CPoint& end2 = m_StarDesc->m_PosTab[last];

				if( (fabs(beg1.x-end1.x)>1 && fabs(beg1.y-end1.y)>1) ||  (fabs(beg2.x-end2.x)>1 && fabs(beg2.y-end2.y)>1) ){
					// only in case >1 shift in each direction - otherwise translation assumed :
					bRotation = TRUE;
	
					double a,b,c,p,q,r;
		
		
					CMyCalcRot::CalcOrtoLine( beg1.x, beg1.y , end1.x, end1.y, a, b, c );
   	   		CMyCalcRot::CalcOrtoLine( beg2.x, beg2.y , end2.x, end2.y, p, q, r );

					if(b!=0)
			   	   _TRACE_PRINTF_5("line 1 is : y = %f*x+%f\n",-a/b,-c/b);
					else
						_TRACE_PRINTF_5("line 1 is : x = %f\n",(-c/a));

					if(q!=0)				
	   			   _TRACE_PRINTF_5("line 2 is : y = %f*x+%f\n",-p/q,-r/q);
					else
						_TRACE_PRINTF_5("line 1 is : x = %f\n",(-r/p));
				

	      		if(!CMyCalcRot::CalcCross( a, b, c, p, q, r, x_cross, y_cross )){
						return FALSE;
					}
			      printf(" Line 1 crosses with line 2 in point (%f,%f)\n",x_cross,y_cross);
					double angle = CMyCalcRot::CalcRotAngle( beg1.x, beg1.y, end2.x, end2.y, x_cross, y_cross );
		      	int frames_no = m_StarDesc->m_nCurrentStep;
					angle_per_frame = angle/frames_no;
					if( totalTime != 0 ){
						angle_per_sec = angle/totalTime;
					}
				
   		   	printf(" Rotation angle = %f ( %f deg ) = /%d => %.20f ( / frame ) \n",angle,angle*(180/PI_VALUE),frames_no,angle_per_frame);	

					m_bRotation = TRUE;
					m_RotCenterX = x_cross;
			      m_RotCenterY = y_cross;
	   		   m_dAlfa = angle_per_frame;
					m_dAlfaPerSec = angle_per_sec;
				}
				return TRUE;
			}else{
				printf("Path of star : %d (%d,%d) not long enough to calc rotation\n",0, (int)m_StarDesc->m_PosTab[0].x,(int)m_StarDesc->m_PosTab[0].y);
				printf("frameDX=%f, frameDY=%f, totalShift=%f\n",frameDX,frameDY,totalShift);
				return TRUE;
			}
		}
	}else{
		double totalDX = 0, totalDY = 0;
		double totalDXPerSec = 0, totalDYPerSec = 0;

		double totalTime = 0;
		if(endTime>m_StartTime){
			totalTime = (endTime-m_StartTime);
		}

		double totalShift=0.00;
		int goodCount=0;
		for(int i=0;i<m_StarsCount;i++){
			totalShift = totalShift + sqrt( ( m_StarDesc->m_PosTab.back().x - m_StarDesc->m_PosTab.front().x )*( m_StarDesc->m_PosTab.back().x - m_StarDesc->m_PosTab.front().x )
										  + ( m_StarDesc->m_PosTab.back().y - m_StarDesc->m_PosTab.front().y )*( m_StarDesc->m_PosTab.back().y - m_StarDesc->m_PosTab.front().y ) );


			if( m_StarDesc[i].m_PosTab.size()>min_steps ){
				m_StarDesc[i].m_starDX = (m_StarDesc[i].m_PosTab.back().x - m_StarDesc[i].m_PosTab.front().x )/(m_StarDesc[i].m_PosTab.size()-1);
				m_StarDesc[i].m_starDY = (m_StarDesc[i].m_PosTab.back().y - m_StarDesc[i].m_PosTab.front().y )/(m_StarDesc[i].m_PosTab.size()-1);
				totalDX += m_StarDesc[i].m_starDX;
				totalDY += m_StarDesc[i].m_starDY;

				m_StarDesc[i].m_starDXPerSec = 0;
				m_StarDesc[i].m_starDYPerSec = 0;
						
				if(endTime>m_StartTime){
					if(totalTime!=0){
						m_StarDesc[i].m_starDXPerSec = (m_StarDesc[i].m_PosTab.back().x - m_StarDesc[i].m_PosTab.front().x  )/totalTime;
						m_StarDesc[i].m_starDYPerSec = (m_StarDesc[i].m_PosTab.back().y - m_StarDesc[i].m_PosTab.front().y )/totalTime;
					}else{
						m_StarDesc[i].m_starDXPerSec = 0;
						m_StarDesc[i].m_starDYPerSec = 0;
					}
				}
				totalDXPerSec += m_StarDesc[i].m_starDXPerSec;
				totalDYPerSec += m_StarDesc[i].m_starDYPerSec;	
				goodCount++;
				m_StarDesc[i].m_bGood = TRUE;
			}else{
				m_StarDesc[i].m_bGood = FALSE;
			}
		}

		m_DX = frameDX = totalDX/goodCount;
		m_DY = frameDY = totalDY/goodCount;
		m_DXPerSec = frameDXPerSec = totalDXPerSec/goodCount;
		m_DYPerSec = frameDYPerSec = totalDYPerSec/goodCount;

		
		bRotation = FALSE;
		m_bRotation = FALSE;

		if(!CheckIfRotation() || 
			( fabs(frameDX)<gCCDParams.m_MinShiftToCalcRot && fabs(frameDY)<gCCDParams.m_MinShiftToCalcRot) ){
			return TRUE;
		}


		// calculate lines ortogonal to star sector :
		for(int i=0;i<m_StarsCount;i++){
			if(m_StarDesc[i].m_bGood){
				CPoint& beg1 = (m_StarDesc[i].m_PosTab.front());
				CPoint& end1 = (m_StarDesc[i].m_PosTab.back());

				CMyCalcRot::CalcOrtoLine( beg1.x, beg1.y , end1.x, end1.y, 
											     m_StarDesc[i].m_OrtoLineA, m_StarDesc[i].m_OrtoLineB, 
											     m_StarDesc[i].m_OrtoLineC );			
			}
		}



		// upper limit (we do not take those from same row and column )
		int nCrosses = m_nSize*(m_nSize+1)/2;

		CPoint* pCrossPoints = new CPoint[nCrosses];

		int nCrossLines=0;
		int halfSize=(m_nSize/2);
		for(int y=0;y<m_nSize;y++){
			for(int x=0;x<m_nSize;x++){
				CStarDesc& star1 = GetStarDesc( x, y );
				int pos = y*m_nSize+x;
				
				// now go through - avoiding repeating :
				for(int yy=y+1;yy<m_nSize;yy++){
					for(int xx=0;xx<m_nSize;xx++){
						int pos2 = yy*m_nSize+xx;							

						if(pos2<=pos)
							continue;
						if(y==yy || x==xx)
							continue;

						CStarDesc& star2 = GetStarDesc( xx, yy );

						if(!star1.m_bGood || !star2.m_bGood)
							continue;

						if(CMyCalcRot::CalcCross( star1.m_OrtoLineA, star1.m_OrtoLineB, 
							  							  star1.m_OrtoLineC, star2.m_OrtoLineA,
														  star2.m_OrtoLineB, star2.m_OrtoLineC,
														  pCrossPoints[nCrossLines].x,
														  pCrossPoints[nCrossLines].y )){
							nCrossLines++;
						}
					}
				}			
			}
		}
		
		if(nCrossLines==0){
			delete [] pCrossPoints;
			return TRUE;
		}

		double angle=0,angle_total=0;

		for(int i=0;i<nCrossLines;i++){
			x_cross += pCrossPoints[i].x;
			y_cross += pCrossPoints[i].y;
		}
		x_cross = x_cross/nCrossLines;
		y_cross = y_cross/nCrossLines;

		for(int i=0;i<m_StarsCount;i++){
			if(m_StarDesc[i].m_bGood){
				CPoint& beg1 = (m_StarDesc[i].m_PosTab.front());
   	      CPoint& end1 = (m_StarDesc[i].m_PosTab.back());

				angle += CMyCalcRot::CalcRotAngle( beg1.x, beg1.y, end1.x, end1.y, x_cross, y_cross )/(m_StarDesc[i].m_PosTab.size()-1);
				angle_total += CMyCalcRot::CalcRotAngle( beg1.x, beg1.y, end1.x, end1.y, x_cross, y_cross );
			}
		}
		angle_per_frame = angle/goodCount;		

		if(totalTime != 0){
			angle_per_sec = (angle_total/goodCount)/totalTime;
		}

		m_bRotation = bRotation = TRUE;
      m_RotCenterX = x_cross;
      m_RotCenterY = y_cross;
      m_dAlfa = angle_per_frame;
      m_dAlfaPerSec = angle_per_sec;

		delete [] pCrossPoints;	
		return TRUE;				
	}

	return FALSE;
}


BOOL_T CCDStarSpy::CheckIfRotation()
{
	BOOL_T bRet = TRUE;

	int nBad=0;

	for(int i=0;i<m_StarsCount;i++){
   	CPoint& beg = (m_StarDesc[i].m_PosTab.front());
   	CPoint& end = (m_StarDesc[i].m_PosTab.back());

		// first check if not horizontal line : -----
		long x0 = (int)beg.x;
		int nGood=0;
		for(int j=0;j<m_StarDesc[i].m_PosTab.size();j++){
			if(fabs(m_StarDesc[i].m_PosTab[j].x-x0)>2)
				nGood++;
		}
		if(nGood<1){
			nBad++;
			continue;
		}

		// then check if not line |
		long y0 = (int)beg.y;
		nGood=0;
		for(int j=0;j<m_StarDesc[i].m_PosTab.size();j++){
			if(fabs(m_StarDesc[i].m_PosTab[j].y-y0)>2)
				nGood++;
		}
		if(nGood<1){
			nBad++;
			continue;
		}

		// now check if not line /
		/*double a = (beg.y - end.y)/(beg.x - end.x);
		double b = (end.y*beg.x-beg.y*end.x)/(beg.x-end.x);
		
		int nAbove=0,nBelow=0;
		int before_last = (m_StarDesc[i].m_PosTab.size()-1);
		nGood = 0;
		double max_dist = -1000.00;
		int max_pos = -1;
		for(int j=1;j<before_last;j++){
			double xx = m_StarDesc[i].m_PosTab[j].x;
			double yy = m_StarDesc[i].m_PosTab[j].y;			
			
			double dist = sqrt( CMyFit::CalcDist2FromLine2Par( a, b, xx, yy ) );
			if(dist > 0.5)
				nGood++;
			if(dist > max_dist){
				max_dist = dist;
				max_pos = j;
			}

			// checking WYPUKLOSC :
			double y = a*xx+b;
			if( yy>=y )
				nAbove++;
			else
				nBelow++;
		}

		int twenty_percent_no = (m_StarDesc[i].m_PosTab.size()-2)*0.2;
		if(nGood<1 || (nBelow>twenty_percent_no && nAbove>twenty_percent_no)){
			nBad++;
			continue;
		}*/
	}

	return (nBad<1);
//	return (nBad<=2);
}

BOOL_T CCDStarSpy::SaveToFile( const char* szFileName, int CameraIndex, 
										 const char* load_file, CCDPipeline* pPipeline )
{
	MyOFile out( szFileName );
	out.Printf("#\tCamIdx\tRotCenterX\tRotCenterY\tdAlfaPerFrame\tFrameDX\tFrameDY\tRot\tdAlfaPerSec\n" );

	double x,y,dx,dy,dxsec,dysec,angle,angle_per_sec;

	x = m_RotCenterX;
   y = m_RotCenterY;
   angle = m_dAlfa;
	angle_per_sec = m_dAlfaPerSec;
   dx = m_DX;
   dy = m_DY;
   dxsec = m_DXPerSec;
   dysec = m_DYPerSec;

	BOOL_T bUseRot = FALSE;
   if( !gCCDParams.m_bDoNotUseRotation ){
   	bUseRot = m_bRotation;
   }

	out.Printf("%d\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%d\t%.10f\t%.10f\t%.10f\n",
					CameraIndex,x,y,angle,dx,dy,bUseRot,dxsec,dysec,angle_per_sec);		

	// saving stars used to calculate shifts :
	SaveStarsTracks( "starstrack.trace", CameraIndex, TRUE, pPipeline );

	if(load_file && load_file[0]){
		MyOFile out( load_file, "a" );

		out.Printf("# AUTO GENERATED BY startrace PROGRAM :\n");
		if(m_StarDesc->m_PosTab.size()<20){
			out.Printf("# WARNING : only %d points were used for calculation\n",m_StarDesc->m_PosTab.size());
		}

		out.Printf("CCD_USE_ROT_IN_AVER_OF_PREV=%d\n",bUseRot);
		out.Printf("CCD_ROT_CENTER_X=%.10f\n",x);
		out.Printf("CCD_ROT_CENTER_Y=%.10f\n",y);
		out.Printf("CCD_SINGLE_FRAME_D_ALFA=%.10f\n",angle);
		out.Printf("CCD_SINGLE_FRAME_D_ALFA_PER_SEC=%.10f\n",angle_per_sec);
		out.Printf("CCD_SINGLE_FRAME_DX=%.10f\n",dx);
		out.Printf("CCD_SINGLE_FRAME_DY=%.10f\n",dy);
		out.Printf("CCD_SINGLE_FRAME_DX_PER_SEC=%.10f\n",dxsec);
      out.Printf("CCD_SINGLE_FRAME_DY_PER_SEC=%.10f\n",dysec);			

	}

	return TRUE;
}


void CCDStarSpy::SaveStarsTracks( const char* fname, int cam_no, 
											 BOOL_T bSaveTrack, CCDPipeline* pPipeline )
{
	// saving stars used to calculate shifts :
	mystring szFileName = fname;
	if( pPipeline ){
		CCDPipeline::GetOutFileName( szFileName, "StarTrackLog", fname,
											  pPipeline, -1 , FALSE );
	}

	MyOFile out2( szFileName.c_str() );
	for(int i=0;i<m_StarsCount;i++){
		CPoint& beg1 = (m_StarDesc[i].m_PosTab.front());
      CPoint& end1 = (m_StarDesc[i].m_PosTab.back());

		 out2.Printf("Star#%d : (%d,%d)-(%d,%d)\n",i,(int)beg1.x,(int)beg1.y , (int)end1.x,(int)end1.y);


		 if( bSaveTrack ){
			  mystring szStarTraceNameTmp,szStarTraceName;
			  szStarTraceNameTmp << "cam" << cam_no << "_star" << i << "_track.trace";			
			  if( pPipeline ){
					CCDPipeline::GetOutFileName( szStarTraceName, "StarTrackLog",
											  szStarTraceNameTmp.c_str(),
											  pPipeline, -1 , FALSE );
			  }else{
				 szStarTraceName = szStarTraceNameTmp;
			  }
			  if( !MyOFile::DoesFileExist( szStarTraceName.c_str() ) ){
					MyOFile out3( szStarTraceName.c_str() );
			   	out3.Printf("Step# (x,y) Ident Frame#\n");		
			  }else{
					MyOFile out3( szStarTraceName.c_str() );
					out3.Printf("\n\n##################################\n\n");
			  }

		     MyOFile out3( szStarTraceName.c_str(), "a" );
			  for(int j=0;j<m_StarDesc[i].m_PosTab.size();j++){
			      out3.Printf("%d (%d,%d) %s %d\n",j,(int)m_StarDesc[i].m_PosTab[j].x,(int)m_StarDesc[i].m_PosTab[j].y,
						CPoint::get_ident_desc(m_StarDesc[i].m_PosTab[j].m_IdentType,TRUE),
						m_StarDesc[i].m_PosTab[j].m_FrameIndex);
			  }
		 }
	}
}

CCDStarSpyTab::CCDStarSpyTab( const char* szName )
: m_szName( szName ), m_pPipeline( NULL ), m_bFromFile(FALSE)
{}

BOOL_T CCDStarSpyTab::SaveToFile( const char* szFileName, const char* load_file,
											 CCDPipeline* pPipeline )
{	
	// MyOFile out( szFileName );
	// out.Printf("#\tCamIdx\tRotCenterX\tRotCenterY\tdAlfaPerFrame\tFrameDX\tFrameDY\tRot\tDXperSec\tDYperSec\n" );

	vector<CCDStarSpy>::iterator i;
	int j=0;

	for(i=begin();i!=end();j++,i++){
		int ccd_index=0;
		if( pPipeline ){
			ccd_index = pPipeline->GetPipelineIndex();
		}
		i->SaveToFile( szFileName, ccd_index, load_file, pPipeline );
		j++;
	}
	
	return TRUE;
}

void CCDStarSpyTab::ResetTracedStars()
{
	vector<CCDStarSpy>::iterator i;
	for(i=begin();i!=end();i++){
		i->ResetTracedStar();
	}
}


BOOL_T CCDStarSpyTab::ReadFromFile( const char* szFileName )
{
	MyIFile in;
	if( !in.Open(szFileName,"r",FALSE) )
		return FALSE;

	const char* pLine;
	
	while( pLine = in.GetLine() ){
		if( pLine[0]=='#' || pLine[0]=='\n' || strlen(pLine)==0)
			continue;
		CMyStrTable items;
		MyParser pars = pLine;
		pars.GetItems( items );
		if( items.size()>=7){
			double dxsec=0,dysec=0,angle_per_sec=0;

			if( items.size()>=9 ){
				dxsec = atof( items[7].c_str() );
				dysec = atof( items[8].c_str() );
			}
			if( items.size()>=10 ){
				angle_per_sec = atof( items[9].c_str() );
			}			

			CCDStarSpy tmp( gCCDParams.m_SizeX, gCCDParams.m_SizeY, atof(items[1]), atof(items[2]), 
								 atof(items[3]), atof(items[4]), atof(items[5]), atol(items[6]), 
								 dxsec, dysec, angle_per_sec, m_pPipeline  );
			push_back( tmp );			
			m_bFromFile = TRUE;
		}
	}

	return (size()>0);
}


BOOL_T CCDStarSpy::CheckIfHaveTwoWithFrame( int frame1, int frame2 )
{
	int nGood=0;
	for(int i=0;i<m_StarsCount;i++){
		if( m_StarDesc[i].m_PosTab.front().m_FrameIndex<=frame1 && m_StarDesc[i].m_PosTab.back().m_FrameIndex>=frame2 )
			nGood++;
	}
	return (nGood>=2);
}

BOOL_T CCDStarSpy::FindTransform( int frame1, int frame2, CCDConfig& params, int ccd_index, CCDPipeline* pPipeline )
{	
	int i=0;
	params.m_bTransformMatrixOK = FALSE;

	if(!CheckIfHaveTwoWithFrame( frame1, frame2 )){
		// in this case changing frame1 to best possible choice :
		int minFrame=10000;
		int minSecFrame=10000;
		int minPos=-1;
		int minSecPos=-1;

		// finding minimum :
		for(i=0;i<m_StarsCount;i++){
			if(m_StarDesc[i].m_PosTab.back().m_FrameIndex>=frame2){
				if(m_StarDesc[i].m_PosTab.front().m_FrameIndex<minFrame){
					minFrame = m_StarDesc[i].m_PosTab.front().m_FrameIndex;
					minPos = i;
				}
			}
		}		

		// finding sec - minimum :
		for(i=0;i<m_StarsCount;i++){
			if(m_StarDesc[i].m_PosTab.back().m_FrameIndex>=frame2 && i!=minPos){
				if(m_StarDesc[i].m_PosTab.front().m_FrameIndex<minSecFrame){
					minSecFrame = m_StarDesc[i].m_PosTab.front().m_FrameIndex;
					minSecPos = i;
				}
			}
		}		
		if(minPos>=0 && minSecPos>=0){
			frame1 = m_StarDesc[minSecPos].m_PosTab.front().m_FrameIndex;
		}
	}

	// finding first star existing on both frames (frame1 and frame2 ) :
	int fStar=-1;
	CPoint* pStar1Frame1 = NULL;
	CPoint* pStar1Frame2 = NULL;

	if(frame1!=frame2){
		for(i=0;i<m_StarsCount;i++){
			pStar1Frame1 = m_StarDesc[i].FindFrame( frame1 );
			pStar1Frame2 = m_StarDesc[i].FindFrame( frame2 );
			
			if(pStar1Frame1 && pStar1Frame2){
				fStar = i;
				break;
			}	
		}
	}


	if( fStar>=0 && fStar<(m_StarsCount-1)){
		// finding second star existing on both frames (frame1 and frame2 ) :
		CPoint* pStar2Frame1 = NULL;
		CPoint* pStar2Frame2 = NULL;
		for(i=(m_StarsCount-1);i>=(fStar+1);i--){
			pStar2Frame1 = m_StarDesc[i].FindFrame( frame1 );
      	pStar2Frame2 = m_StarDesc[i].FindFrame( frame2 );	

			if(pStar2Frame1 && pStar2Frame2){
				break;
			}
		}
	
		if( pStar1Frame1 && pStar1Frame2 && pStar2Frame1 && pStar2Frame2 ){
			// now finding transformation :
			CPoint Frame2_to_Frame1_Shift;
			Frame2_to_Frame1_Shift.Subtract( (*pStar1Frame1) , (*pStar1Frame2) );
		

			CPoint f2s2Shifted;
			f2s2Shifted.Add( (*pStar2Frame2) , Frame2_to_Frame1_Shift );


			// currently we have :
			// pStar2Frame1 - original position of star2 
			// Frame2Star2Shifted - star 2 on frame 2 after shift 
			// center of rotation in - star1 = pStar1Frame1
			// so now we find rotation angle :
			CPoint center = (*pStar1Frame1);


	/*		double r = sqrt( CMyMathFunc::mysqr( center.x - pStar2Frame1->x ) + CMyMathFunc::mysqr( center.y - pStar2Frame1->y ) );
			double d = CPoint::dist( (*pStar2Frame1) , f2s2Shifted );
			double sin_alfa = (0.5*d)/r;
			double alfa = 2.00*asin( sin_alfa );		*/
			double alfa = CMyCalcRot::CalcRotAngle( pStar2Frame1->x, pStar2Frame1->y, f2s2Shifted.x, f2s2Shifted.y, center.x, center.y );


			mystring szTrans;
			szTrans << "Frame# " << pPipeline->GetFrameIndex() << "\n";
			szTrans << "Transformation of : " << frame1 << " -> " << frame2 <<"\n";
			szTrans << "Translation = (" << Frame2_to_Frame1_Shift.x << "," << Frame2_to_Frame1_Shift.y << "\n";
			szTrans << "Rotation around star 1 , alfa=" << alfa << "\n";
			szTrans << "Transform matrix :\n";
			szTrans << cos(alfa) << "    " << -sin(alfa) << "    " << -Frame2_to_Frame1_Shift.x << "\n";
			szTrans << sin(alfa) << "    " << cos(alfa) << "    " << -Frame2_to_Frame1_Shift.y << "\n";
		
			double dt = ( pStar1Frame2->frame_time-pStar1Frame1->frame_time );
			double angle_per_sec = (alfa/dt);
			double dx_per_sec = -( Frame2_to_Frame1_Shift.x/dt );
			double dy_per_sec = -( Frame2_to_Frame1_Shift.y/dt );
			params.m_FrameTransformMatrix.Init( angle_per_sec, dx_per_sec, dy_per_sec, center.x, center.y );
			params.m_bTransformMatrixOK = TRUE;

			mystring szOutName,szInfo;
			szInfo << "With respect to star " << fStar << " at (" << (int)pStar1Frame1->x
					 << "," << (int)pStar1Frame1->y << ") and using star " << i
					 << " at (" <<  (int)pStar2Frame1->x << "," << (int)pStar2Frame2->x << ")";
			GetOutFileName( szOutName, "TransformMatrixPerTime", ccd_index, pPipeline );
			params.m_FrameTransformMatrix.Dump( szOutName.c_str(), pPipeline->GetFrameIndex(), szInfo.c_str() );
			params.m_nMatrixNotFoundCount = 0;
	

			for(int j=0;j<m_StarsCount;j++){
				CPoint* pFrame1 = m_StarDesc[j].FindFrame( frame1 );
				CPoint* pFrame2 = m_StarDesc[j].FindFrame( frame2 );
				if(pFrame1 && pFrame2){
					szTrans << "Star " << j << " (" << (int)pFrame1->x << "," << (int)pFrame1->y
							  << " -> (" << (int)pFrame2->x << "," << (int)pFrame2->y << ")\n";
				}else{
					if(pFrame1){
						szTrans << "Star " << j << " (" << (int)pFrame1->x << "," << (int)pFrame1->y
						<< " -> NOT ON FRAME =" << frame2 << "\n";
					}
				}
			}
			_TRACE_PRINTF_2("%s\n",szTrans.c_str());
			if(gCCDParams.m_bDumpTraceStarToFile){
				DumpToFile( szTrans.c_str(), "TransformMatrix.txt", ccd_index, pPipeline );
			}
			return TRUE;
		}
	}
	
	if(params.m_nMatrixNotFoundCount<gCCDParams.m_MaxFramesWithOldMatrix){
		params.m_bTransformMatrixOK = TRUE;
		params.m_nMatrixNotFoundCount++;
		mystring szOutName;
      GetOutFileName( szOutName, "TransformMatrixPerTime", ccd_index, pPipeline );
      params.m_FrameTransformMatrix.Dump( szOutName.c_str(), pPipeline->GetFrameIndex(), "RE-USING OLD MATRIX" );

	}

	return FALSE;
}

void CCDStarSpy::GetOutFileName( mystring& szOutName, const char* szBaseName, int ccd_index, CCDPipeline* pPipeline )
{
	CCDPipeline::GetOutFileName( szOutName, "Transform", szBaseName, pPipeline, ccd_index );
}

void CCDStarSpy::DumpToFile( const char* szReport, const char* szFileName, int ccd_index, CCDPipeline* pPipeline )
{
		mystring szFile;	
 		szFile << gCCDParams.GetOutputDir() << "/Transform/";
	   if( pPipeline->GetPipelineCount()>1 )
	      szFile << "Cam" << pPipeline->GetPipelineIndex() << "/";
		szFile << szFileName << "_ccd" << ccd_index << ".txt";

		MyOFile out( szFile.c_str(), "a+" );
		out.Printf("%s\n",szReport );		
}


void CCDStarSpy::LogCurrentPosition(  const char* szBaseName, int ccd_index, CCDPipeline* pPipeline )
{
	mystring szOutName;
	CCDPipeline::GetOutFileName( szOutName, "Transform", szBaseName, pPipeline, ccd_index );
	SaveCurrentPosition( szOutName );
}


BOOL_T CCDStarSpy::CalcPredicted( CCDPipeline* pPipeline, double curr_time, int ccd_index )
{
	for(int i=0;i<m_StarsCount;i++){
		CStarDesc& star = m_StarDesc[i];
		
		if( star.m_PosTab.size()>=2 ){
			CPoint& prevPos = star.m_PosTab[ (star.m_PosTab.size()-1) ];
			
			CCDDataResults::CalcStarPositionAuto( prevPos.x, prevPos.y, prevPos.frame_time,
															  curr_time,															  
															  (pPipeline->GetFrameIndex()-prevPos.m_FrameIndex),
															  star.m_PredictedX, star.m_PredictedY,
															  &(pPipeline->GetCCDInfoTab()[ccd_index]),
															  &(pPipeline->m_PipelineCfg) );
		}
	}
	return TRUE;
}

void CCDStarSpy::SaveCurrentPosition( const char* fname )
{
	if(!MyFile::DoesFileExist( fname )){
		MyOFile out(fname,"w");
		out.Printf("Star#\tFrame#\tCurrPos\tTime\tInFrame\tValue\tIdent\tPredictedAt\n");
	}

	MyOFile out(fname,"a");
	for(int i=0;i<m_StarsCount;i++){
		CStarDesc& star = m_StarDesc[i];

		char szDesc[1024];
		sprintf(szDesc,"%d %d %.2f %.2f %d %d %.2f %s %.2f %.2f\n",
					i, star.m_PosTab.back().m_FrameIndex, 
					star.m_PosTab.back().x, star.m_PosTab.back().y,
					star.m_PosTab.back().frame_time,
					star.m_bNotEnd, star.m_PosTab.back().m_Value, 
					CPoint::get_ident_desc( star.m_PosTab.back().m_IdentType, star.m_bNotEnd ),
					star.m_PredictedX, star.m_PredictedY );
		out.Printf("%s",szDesc);
	}		
}

void CCDStarSpyTab::DumpReportToFile()
{
	vector<CCDStarSpy>::iterator i;
   int j=0;

   for(i=begin();i!=end();j++,i++){
		mystring szFile;	
 		szFile << gCCDParams.GetOutputDir() << "/Transform/";
	   if( m_pPipeline->GetPipelineCount()>1 )
	      szFile << "Cam" << m_pPipeline->GetPipelineIndex() << "/";
		szFile << m_szName << "_ccd" << j << ".txt";

		MyOFile out( szFile.c_str(), "a+" );
		for(int k=0;k<i->m_StarsCount;k++){
			mystring szDesc;
			if( (i->m_StarDesc[k]).m_PosTab.size() ){
				szDesc << "Star:" << k << " Frame:" << (i->m_StarDesc[k]).m_PosTab.back().m_FrameIndex
					 << " Pos:(" << (int)(i->m_StarDesc[k]).m_PosTab.back().x << "," << (int)(i->m_StarDesc[k]).m_PosTab.back().y 
					 << ") Time:" << (i->m_StarDesc[k]).m_PosTab.back().frame_time
					 << "  InFrame:" << (i->m_StarDesc[k]).m_bNotEnd 
					 << "  Value:" << (i->m_StarDesc[k]).m_PosTab.back().m_Value
					 << "  Ident:" << CPoint::get_ident_desc( (i->m_StarDesc[k]).m_PosTab.back().m_IdentType, (i->m_StarDesc[k]).m_bNotEnd );
				out.Printf("%s\n",szDesc.c_str());
			}
		}
	}
}
