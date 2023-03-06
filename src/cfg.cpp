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
#include "cfg.h"
#include <cexcp.h>
#include <myutil.h>
#include "ccd_piman_interf_defs.h"


// CCfgFile* CCcdCfg::m_CfgFile=NULL;

CCcdCfg* gCCDCfg = NULL;


CCDParams::CCDParams()
: m_RotCenterX(0), m_RotCenterY(0), m_bUseRotInAverageOfPrevN(FALSE), 
	m_TransformCCDFocus(0), m_TransformCCDPixelSize(0), m_TransformCCDOrientation(0),
	m_HorAzimuth(0), m_HorAltitude(0), m_RotValueDAlfa(0), m_SizeX(0), 
	m_SizeY(0), m_PortNo(DAQ_PORT), m_nAutoShiftsCalc(0), m_FrameDX(0), m_FrameDY(0),
	m_FrameDXPerSec(0), m_FrameDYPerSec(0), m_CCDShiftX(0), m_CCDShiftY(0), m_eCAMType( eFileSimulator ),
	m_nNewLaplace(0), m_nMaxLaplaceOnOther(0), m_RotValueDAlfaPerSec(0), m_nNewLaplaceInSigma(0),
	m_nMaxLaplaceOnOtherInSigma(0), m_nIgnoreEdge(0), m_nIgnoreEdgeRight(0), m_nIgnoreEdgeLeft(0),
	m_nIgnoreEdgeUp(0), m_nIgnoreEdgeBottom(0), m_bTransformMatrixOK(FALSE),m_nMatrixNotFoundCount(0),
	m_CCDAxixXOrientation(1), m_CCDAxixYOrientation(-1), m_ObsMode(eNoMovingMode),
   m_DecObs(0), m_RAObs(0), m_HorAzimCorr(0), m_HorAltCorr(0), m_DecCorr(0), m_RACorr(0)
{
	m_szDarkFrameFile[0] = '\0';
	m_szFramesListFile[0] = '\0';
	m_szBaseFileNameFITS[0] = '\0';
}

void CCDParams::SetIgnoreEdges( double left, double right, double bottom, double up )
{
   m_nIgnoreEdgeRight = my_round(right);
   m_nIgnoreEdgeLeft = my_round(left);
   m_nIgnoreEdgeUp = my_round(up);
   m_nIgnoreEdgeBottom = my_round(bottom);
}


void CCcdCfg::GetKnownParams()
{
	const char* szTmp = GetParam( "CCD_ROT_CENTER_X", TRUE );
	if(szTmp && szTmp[0])
		m_CCDParams.m_RotCenterX = atof( szTmp );

	szTmp = GetParam( "CCD_ROT_CENTER_Y", TRUE );
	if(szTmp && szTmp[0])
		m_CCDParams.m_RotCenterY = atof( szTmp );

	szTmp = GetParam( "CCD_SINGLE_FRAME_D_ALFA",TRUE);
	if(szTmp && szTmp[0])
		m_CCDParams.m_RotValueDAlfa = atof( szTmp );

	szTmp = GetParam("CCD_SINGLE_FRAME_D_ALFA_PER_SEC",TRUE);
   if(szTmp && szTmp[0])
   	m_CCDParams.m_RotValueDAlfaPerSec = atof( szTmp );


	szTmp = GetParam( "CCD_USE_ROT_IN_AVER_OF_PREV", TRUE );
	if(szTmp && szTmp[0])
		m_CCDParams.m_bUseRotInAverageOfPrevN = (atol(szTmp)>0);

	szTmp = GetParam("CCD_TRANSFORM_ORIENTATION",TRUE);
   if(szTmp && szTmp[0]){
   	m_CCDParams.m_TransformCCDOrientation = atof(szTmp);
   }

	szTmp = GetParam("CCD_AXIS_X_ORIENTATION",TRUE);
	if(szTmp && szTmp[0]){
		m_CCDParams.m_CCDAxixXOrientation = atol(szTmp);
	}

	szTmp = GetParam("CCD_AXIS_Y_ORIENTATION",TRUE);
	if(szTmp && szTmp[0]){
		m_CCDParams.m_CCDAxixYOrientation = atol(szTmp);
	}

	szTmp = GetParam("CCD_TRANSFORM_CCD_FOCUS",TRUE);
   if(szTmp && szTmp[0]){
   	m_CCDParams.m_TransformCCDFocus = atof( szTmp );
   }

	szTmp = GetParam("CCD_TRANSFORM_CCD_PIXEL_SIZE",TRUE);
   if(szTmp && szTmp[0]){
   	m_CCDParams.m_TransformCCDPixelSize = atof( szTmp );
   }


	// horizontal coordinates of observation :
   szTmp = GetParam("CCD_HOR_AZIMUTH",TRUE);
   if(szTmp && szTmp[0])
   	m_CCDParams.m_HorAzimuth = atof( szTmp );

   szTmp = GetParam("CCD_HOR_ALTITUDE",TRUE);
   if(szTmp && szTmp[0])
   	m_CCDParams.m_HorAltitude = atof( szTmp );

	szTmp = GetParam("CCD_HOR_ALT_CORR",TRUE);
	if(szTmp && szTmp[0])
		m_CCDParams.m_HorAltCorr = atof( szTmp );

	szTmp = GetParam("CCD_HOR_AZIM_CORR",TRUE);
	if(szTmp && szTmp[0])
		m_CCDParams.m_HorAzimCorr = atof( szTmp );


	szTmp = GetParam("CCD_RA_CORR",TRUE);
	if(szTmp && szTmp[0])
		m_CCDParams.m_RACorr = atof( szTmp );

	szTmp = GetParam("CCD_DEC_CORR",TRUE);
	if(szTmp && szTmp[0])
		m_CCDParams.m_DecCorr = atof( szTmp );

	szTmp = GetParam("CCD_OBS_MODE",TRUE);
	if(szTmp && szTmp[0])
		m_CCDParams.m_ObsMode = (eObservationMode_T)(atol(szTmp));

	szTmp = GetParam("CCD_DEC_OBS",TRUE);
	if(szTmp && szTmp[0])
		m_CCDParams.m_DecObs = atof(szTmp);

	szTmp = GetParam("CCD_RA_OBS",TRUE);
	if(szTmp && szTmp[0])
		m_CCDParams.m_RAObs = atof(szTmp);

	// sizes of pipeline:
	szTmp = GetParam("CCD_SIZE_X",TRUE);
	if(szTmp && szTmp[0]){
		m_CCDParams.m_SizeX = atol( szTmp );
	}
	szTmp = GetParam("CCD_SIZE_Y",TRUE);
	if(szTmp && szTmp[0]){
		m_CCDParams.m_SizeY = atol( szTmp );
	}

	// communication:
	szTmp = GetParam("CCD_PORT_NO",TRUE);
	if(szTmp && szTmp[0]){
		m_CCDParams.m_PortNo = atol( szTmp );
	}

	szTmp = GetParam("CCD_DARK_FRAME_FILE",TRUE);
	if(szTmp && szTmp[0]){
		strcpy(m_CCDParams.m_szDarkFrameFile,szTmp);
	}

	szTmp = GetParam("CCD_AUTO_SHIFTS_CALC",TRUE);
	if(szTmp && szTmp[0]){
		m_CCDParams.m_nAutoShiftsCalc = atol( szTmp );
	}		

	szTmp = GetParam("CCD_SHIFT_TO_OTHER_CCD_X",TRUE);
	if(szTmp && szTmp[0]){
		m_CCDParams.m_CCDShiftX = atof( szTmp );
	}
	szTmp = GetParam("CCD_SHIFT_TO_OTHER_CCD_Y",TRUE);
	if(szTmp && szTmp[0]){
		m_CCDParams.m_CCDShiftY = atof( szTmp );
	}


	szTmp = GetParam("CCD_CAM_TYPE",TRUE);
	if(szTmp && szTmp[0]){
		m_CCDParams.m_eCAMType = (eCCDTYPE_T)(atol(szTmp));
	}

	szTmp = GetParam("CCD_FULL_FRAMES_LIST",TRUE);
	if(szTmp && szTmp[0]){
		strcpy( m_CCDParams.m_szFramesListFile, szTmp);
	}

	szTmp = GetParam("CCD_SINGLE_FRAME_DX", TRUE );
   if(szTmp && szTmp[0])
   	m_CCDParams.m_FrameDX = atof( szTmp );

   szTmp = GetParam("CCD_SINGLE_FRAME_DY", TRUE );
   if(szTmp && szTmp[0])
   	m_CCDParams.m_FrameDY = atof( szTmp );

   szTmp = GetParam("CCD_SINGLE_FRAME_DX_PER_SEC",TRUE);
   if(szTmp && szTmp[0])
   	m_CCDParams.m_FrameDXPerSec = atof( szTmp );

   szTmp = GetParam("CCD_SINGLE_FRAME_DY_PER_SEC",TRUE);
   if(szTmp && szTmp[0])
   	m_CCDParams.m_FrameDYPerSec = atof( szTmp );

	szTmp = GetParam("CCD_BASE_FILE_NAME_FITS",TRUE);
	if(szTmp && szTmp[0]){
		strcpy( m_CCDParams.m_szBaseFileNameFITS, szTmp );
	}

	szTmp = GetParam("CCD_NEW_LAPLACE_TRESHOLD_IN_SIGMA",TRUE);
   if(szTmp && szTmp[0]){
   	m_CCDParams.m_nNewLaplaceInSigma = atof( szTmp );
   }

	szTmp = GetParam("CCD_MAX_LAPLACE_ON_OTHER_IN_SIGMA",TRUE);
   if(szTmp && szTmp[0]){
   	m_CCDParams.m_nMaxLaplaceOnOtherInSigma = atof( szTmp );
   }

	szTmp = GetParam("CCD_IGNORE_EDGE",TRUE);
	if(szTmp && szTmp[0]){
		m_CCDParams.m_nIgnoreEdge = atol( szTmp );
	}
	szTmp = GetParam("CCD_IGNORE_EDGE_RIGHT",TRUE);
   if(szTmp && szTmp[0]){
   	m_CCDParams.m_nIgnoreEdgeRight = atol( szTmp);
	}
   szTmp = GetParam("CCD_IGNORE_EDGE_LEFT",TRUE);
   if(szTmp && szTmp[0]){
   	m_CCDParams.m_nIgnoreEdgeLeft = atol( szTmp);
	}
   szTmp = GetParam("CCD_IGNORE_EDGE_UP",TRUE);
   if(szTmp && szTmp[0]){
   	m_CCDParams.m_nIgnoreEdgeUp = atol( szTmp);
	}
   szTmp = GetParam("CCD_IGNORE_EDGE_BOTTOM",TRUE);
   if(szTmp && szTmp[0]){
   	m_CCDParams.m_nIgnoreEdgeBottom = atol( szTmp);
	}
	


	if(m_InitCount==0){
		// first read saving value of m_nAutoShiftsCalc
		m_CCDParams.m_nAutoShiftsCalcSav = m_CCDParams.m_nAutoShiftsCalc;
	}
}

void InitGlobalParamObj()
{
	static BOOL_T bFirstCall=TRUE;
	if(bFirstCall){
		gCCDCfg = new CCcdCfg();
		bFirstCall = FALSE;
		// printf("FIRST CALL to InitGlobalParamObj - INITIALIZED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	}
}

class CGlobalCfgInitializer
{
public:
	CGlobalCfgInitializer();
	~CGlobalCfgInitializer();
};

CGlobalCfgInitializer::CGlobalCfgInitializer()
{
	InitGlobalParamObj();
}

CGlobalCfgInitializer::~CGlobalCfgInitializer()
{
	if(gCCDCfg)
		delete gCCDCfg;
}

CGlobalCfgInitializer gInitObj;

void RefreshGlobalParams()
{
	InitGlobalParamObj();
// 	gCCDCfg->
}

const char* GetGlobalParam( const char* cfgcode,BOOL_T bAllowNull ) 
{
	InitGlobalParamObj();
	return gCCDCfg->GetParam( cfgcode, bAllowNull );
}

void DumpGlobalParams()
{
	gCCDCfg->Dump();
}


void SetGlobalParam( const char* cfgcode, const char* cfgval )
{
	InitGlobalParamObj();
	gCCDCfg->SetParam( cfgcode, cfgval );
}

CCfgFile& GetGlobalParamFile()
{
	InitGlobalParamObj();
	return gCCDCfg->GetParamFile();
}

CCfgFile* GetGlobalParamFilePtr()
{
	InitGlobalParamObj();
	return gCCDCfg->GetParamFilePtr();
}


void GetGlobalParams( mystring& szParams )
{
	InitGlobalParamObj();
	gCCDCfg->GetParams( szParams );
}

CCcdCfg::CCcdCfg( const char* szFileName )
: m_CfgFile(NULL), m_szFileName(szFileName), m_bInitialized(FALSE), m_InitCount(0),
  m_bOverwriteMode(FALSE)
{}

CCcdCfg::CCcdCfg( const CCcdCfg& right )
: m_CfgFile(NULL), m_bInitialized(FALSE), m_InitCount(0), m_bOverwriteMode(FALSE)
{
	(*this) = right;
}

CCcdCfg& CCcdCfg::operator=( const CCcdCfg& right )
{
	memcpy( &m_CCDParams, &right.m_CCDParams, sizeof(CCDParams) );

	m_szFileName = right.m_szFileName;
	m_bOverwriteMode = right.m_bOverwriteMode;

	if( m_CfgFile ){
		delete m_CfgFile;
		m_CfgFile = NULL;
	}
	if( right.m_CfgFile ){
		m_CfgFile = new CCfgFile( (*right.m_CfgFile) );
		m_bInitialized = TRUE;
	}
	m_InitCount = right.m_InitCount;	

	return (*this);
}

CCcdCfg::~CCcdCfg()
{
	if(m_CfgFile){
		delete m_CfgFile;
	}
}

void CCcdCfg::SetFileName ( const char* szFileName )
{
	m_szFileName = szFileName;
}

BOOL_T CCcdCfg::Init( BOOL_T bAssert/*=TRUE*/, BOOL_T bUseAll/*=TRUE*/ )
{
	if(!m_bInitialized){
		BOOL_T bFound=FALSE;
		m_szFileName.env2str();
		// printf("1 Using file : %s\n",m_szFileName.c_str());
		if(!MyFile::DoesFileExist(m_szFileName.c_str())){			
			// printf("File %s not found !\n",m_szFileName.c_str());
			if(!m_bOverwriteMode){				
				if(bUseAll){
					if( strcmp(m_szFileName.c_str(),DEFAULT_CFG_FILE) ){
   		         if(MyFile::DoesFileExist(DEFAULT_CFG_FILE)){
							m_szFileName = DEFAULT_CFG_FILE;
							bFound = TRUE;
						}
					}
					if(!bFound){
						m_szFileName = "$(CFGFILE)";
						m_szFileName.env2str();
						if(MyFile::DoesFileExist( m_szFileName.c_str() )){
							bFound=TRUE;
						}
					}
				}
			}else{
				bAssert=FALSE;
			}
		}else{
			bFound = TRUE;
		}
		// printf("Using file : %s\n",m_szFileName.c_str());
		if(bAssert){
			Assert(bFound,"No configuration file can be found, %s - does not exist",m_szFileName.c_str());
		}

		if(bFound){
			m_CfgFile = new CCfgFile(m_szFileName.c_str());
		}
		m_bInitialized = TRUE;

		if(m_CfgFile)
			GetKnownParams();
		
		m_InitCount++;
	}
	return (m_CfgFile!=NULL);
}

const char* CCcdCfg::GetParam(const char* cfgcode,BOOL_T bAllowNull)
{
	Init();
	if(m_CfgFile)
		return m_CfgFile->GetParam(cfgcode,bAllowNull);
	else
		return NULL;
}

const char* CCcdCfg::GetParamNoInit(const char* cfgcode,BOOL_T bAllowNull)
{
	if(m_CfgFile)
		return m_CfgFile->GetParamNoInit(cfgcode,bAllowNull);
	else
		return NULL;
}


void CCcdCfg::SetParam( const char* cfgcode, const char* cfgval )
{
	Init();
	if(m_CfgFile)
		m_CfgFile->SetParam( cfgcode, cfgval );		
}


void CCcdCfg::GetParams( mystring& szParams )
{
	szParams = "";
	if(m_CfgFile)
		m_CfgFile->GetParams( szParams );
}


void CCcdCfg::Dump()
{
	if(m_CfgFile)
		m_CfgFile->Dump();
}

