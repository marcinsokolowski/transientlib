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
// syscfgtxt.h: interface for the CSysCfgTxt class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_SYSCFGTXT_H__5212D431_3F44_11D6_B636_B42D07C10000__INCLUDED_)
#define AFX_SYSCFGTXT_H__5212D431_3F44_11D6_B636_B42D07C10000__INCLUDED_

#include <cfgfile.h>
#include "ccd_common_defines.h"
#include <mymatrix.h>

const char* GetGlobalParam( const char* cfgcode,BOOL_T bAllowNull=FALSE );
void SetGlobalParam( const char* cfgcode, const char* cfgval );
CCfgFile& GetGlobalParamFile();
CCfgFile* GetGlobalParamFilePtr();
void GetGlobalParams( mystring& szParams );
void DumpGlobalParams();
void RefreshGlobalParams();

// CCcdCfg* GetGlobalCccdCfg();

class CCDParams
{
public:
	CCDParams();

	void SetIgnoreEdges( double left, double right, double bottom, double up );

	char m_szBaseFileNameFITS[64];
	char m_szDarkFrameFile[256];
	char m_szFramesListFile[512];

	eCCDTYPE_T m_eCAMType;

	int m_SizeX;
	int m_SizeY;	

	int m_nAutoShiftsCalc;
	int m_nAutoShiftsCalcSav;

   double m_RotCenterX;
   double m_RotCenterY;
	double m_RotValueDAlfa;
	double m_RotValueDAlfaPerSec;
	BOOL_T m_bUseRotInAverageOfPrevN;
	double m_FrameDX;
	double m_FrameDY;
	double m_FrameDXPerSec;
	double m_FrameDYPerSec;		   	
   
                                       

	double m_TransformCCDFocus;
	double m_TransformCCDPixelSize;
	double m_TransformCCDOrientation;
	int m_CCDAxixXOrientation;
	int m_CCDAxixYOrientation;
	
	// transformation matrix - time dependent :
	CTransformMatrixInTime m_TransformMatrix;
	BOOL_T m_bTransformMatrixOK;
	int m_nMatrixNotFoundCount;
	
	
	// observation :
	double m_HorAzimuth;
	double m_HorAltitude;

	// corrections to position obtained from MOUNT to center of CCD :
	double m_HorAzimCorr;
	double m_HorAltCorr;
	
	double m_DecCorr;
	double m_RACorr;

	eObservationMode_T m_ObsMode;
	double m_DecObs;
	double m_RAObs;


	
	// communication :
	int m_PortNo;


	// left-bottom corner of second cam vs first cam :
	double m_CCDShiftX;
	double m_CCDShiftY;
	
	
	// threasholds are - different - because sigma is different :
	LONG_T m_nNewLaplace;
	LONG_T m_nMaxLaplaceOnOther;
	double m_nNewLaplaceInSigma;
	double m_nMaxLaplaceOnOtherInSigma;
	
	
	// ignore edge :
	LONG_T m_nIgnoreEdge;
	LONG_T m_nIgnoreEdgeRight;
	LONG_T m_nIgnoreEdgeLeft;
	LONG_T m_nIgnoreEdgeUp;
	LONG_T m_nIgnoreEdgeBottom;	            
};
         
         

class CCcdCfg
{
public:
	CCDParams m_CCDParams;

	CCcdCfg( const char* szFileName=DEFAULT_CFG_FILE );
	CCcdCfg( const CCcdCfg& right );
	CCcdCfg& operator=( const CCcdCfg& right );
	
	virtual ~CCcdCfg();
	BOOL_T Init( BOOL_T bAssert=TRUE, BOOL_T bUseAll=TRUE );
	void SetFileName ( const char* szFileName );

	// reading some known parameters to structure m_CCDParams
	void GetKnownParams();

	const char* GetParam(const char* cfgcode,BOOL_T bAllowNull=FALSE);
	const char* GetParamNoInit(const char* cfgcode,BOOL_T bAllowNull=FALSE);

	void SetParam( const char* cfgcode, const char* cfgval );
	void GetParams( mystring& szParams );
	void Dump();
	inline CCfgFile& GetParamFile() { return (*m_CfgFile); }
	inline CCfgFile* GetParamFilePtr() { return m_CfgFile; }

	mystring m_szFileName;
	BOOL_T m_bOverwriteMode;
	
	const char* GetName(){ return m_szFileName.c_str(); }
	BOOL_T IsInitialized(){ return m_bInitialized; }
protected :
	CCfgFile* m_CfgFile;
	BOOL_T m_bInitialized;
	int m_InitCount;
};

extern CCcdCfg* gCCDCfg;

#endif // !defined(AFX_SYSCFGTXT_H__5212D431_3F44_11D6_B636_B42D07C10000__INCLUDED_)
