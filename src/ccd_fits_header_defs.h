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
#ifndef _CCD_FITS_HEADER_DEFS_H__
#define _CCD_FITS_HEADER_DEFS_H__


#define DEFAULT_ORIGIN "PI OF THE SKY"
// #define DEFAULT_SITE   "LCO"
#define DEFAULT_SITE "BRWINOW"
#define DEFAULT_INSTRUME "PiOfTheSky 0"
#define DEFAULT_CAMOPTIC "50/1.4 Zeiss Planar"
#define DEFAULT_FILTER  "V"
#define DEFAULT_OBSERVER "PiMan"
#define DEFAULT_OBJECT "HETE"
#define CANON_85mm_FOV_SIZE 22.00

#define USHORT_BZERO   32768
#define USHORT_BSCALE  1

#define FH_SIMPLE "SIMPLE"
#define FH_BITPIX "BITPIX"
#define FH_NAXIS  "NAXIS"
#define FH_NAXIS1 "NAXIS1"
#define FH_NAXIS2 "NAXIS2"
#define FH_EXTEND "EXTEND"
#define FH_BZERO  "BZERO"
#define FH_BSCALE "BSCALE"
#define  FH_COMMENT "COMMENT"


// IMAGE INFO :
#define TAKENFROM "TAKENFROM"
#define PARTNAME "PARTNAME"
#define FLIP     "FLIP"
#define DRVFLIP  "DRVFLIP"
#define SECTION  "SECTIO"
#define COOLING  "COOLING"
#define SPEED    "SPEED"
#define REXPTIME "REXPTIME"
#define ROTIME   "RO-TIME"
#define CROTIME "CROTIME"
#define SHUTTER  "SHUTTER"
#define MPP_BC   "MPP_BC"
#define UT_START "UT-START"
#define DATE_OBS "DATE-OBS"
#define DATE_OBS2 "DATE_OBS"
#define LOCTIME  "LOCTIME"
#define LOCDATE  "LOCDATE"
#define UT_END   "UT-END"
#define DATE_END "DATE-END"
#define NIGHT    "NIGHT"
#define EXPOSURE "EXPOSURE"
#define MIPS_HI  "MIPS-HI"
#define MIPS_LO  "MIPS-LO"
#define CONTRAST "CONTRAST"
#define MEAN     "MEAN"
#define SIGMA    "SIGMA"
#define DATAMAX  "DATAMAX"
#define DATAMIN  "DATAMIN"
#define BGMEAN   "BGMEAN"  
#define BGSIGMA  "BGSIGMA"
#define CHIP_TEMP "CHIPTEMP"
#define CHIP_TEM  "CHIP_TEM"
#define CHISTEMP  "CHIPTSET"
#define CASE_TEMP "CASTEMP"
#define AMBT_TEMP "AMBTEMP"
#define ELECGAIN  "ELECGAIN"
#define ADCGAIN   "ADCGAIN"
#define ADCBIAS   "ADCBIAS"
#define ADCBSET	"ADCBSET"
#define ADCGSET	"ADCGSET"
#define ADCCLAMP  "ADCCLAMP"
#define RA_OBS   "RA"
#define DEC_OBS  "DEC"
#define ALT_OBS  "ALT"
#define AZIM_OBS "AZIM"
#define HA_OBS   "HA"
#define TIME_OBS "TIME_OBS"
#define TIME_UT  "TIME_UT"
#define TELLONG  "TELLONG"
#define TELALT   "TELALT"
#define FOCUS    "FOCUS"
#define FOCALLEN "FOCALLEN"
#define PIXSCALE "PIXSCALE"
#define PIXSIZE  "PIXSIZE"
#define NIMAGE   "NIMAGE"
#define DIMAGE   "DIMAGE"
#define JD       "JD"
#define HJD      "HJD"
#define ST       "ST"
#define AIRMASS  "AIRMASS"
#define EPOCH    "EPOCH"
#define ZENITH_D "ZENITH_D"
#define FILENAME "FILENAME"
#define INSTRUME "INSTRUME"
#define CAMOPTIC "CAMOPTIC"
#define FILTER   "FILTER"
#define SOFTVERSION "SOFTWARE"
#define SOFTBUILD "BUILD"
#define OBSERVER "OBSERVER"
#define ORIGIN "ORIGIN"
#define SITE "SITE"
#define SAVEAREA "SAVEAREA"
#define ROTATE "ROTATE"
#define CAMERA "CAMERA"
#define TELLAT "TELLAT"
#define EXPTIME "EXPTIME"
#define SHTIME  "SHTIME"
#define SPEEDSET "SPEEDST"
#define SPEEDMHZ "SPEEDMH"
#define CAMHUMID "CAMHUMID"
#define AMBHUMID "AMBHUMID"
#define INTRTEMP "INTRTEMP"
#define ABINN    "ABINN"
#define SBINN    "SBINN"
#define DRVTYPE  "DRVTYPE"
#define OBJECT   "OBJECT"
#define EPOCH    "EPOCH"
#define EQUINOX  "EQUINOX"
#define RNOISE   "RNOISE"
#define NAXIS3   "NAXIS3"
#define ADCRANGE "ADCRANGE"
#define USBMODE  "USBMODE"
#define FPGAVER  "FPGAVER"
#define CPRSVER  "CPRSVER"
#define VERDESC  "VERDESC"
#define HITLENS  "HITLENS"
#define CAMIIDX   "CAMIIDX"
#define CAMID     "CAMID"
#define LNAGAIN  "LNAGAIN"
#define RELNOISE "RELNOISE"
#define COMPRES "COMPRES"
#define BLOCKSZ "BLOCKSZ"
#define FH_DIVISOR "DIVISOR"
#define TAKEFROM "TAKEFROM"
#define ORGSIZEX "ORGSIZEX"
#define ORGSIZEY "ORGSIZEY"
#define MOUNTRA  "MOUNTRA"
#define MOUNTDEC "MOUNTDEC"
#define MOUNTAZIM "MOUNTAZ"
#define MOUNTALT "MOUNTALT"
#define MOUNTTRK "MOUNTTRK"
#define MOUNTHA  "MOUNTHA"
#define MOUNTTM  "MOUNTTM"
#define MOUNTDTM "MOUNTDTM"
#define DOME     "DOME"
#define PARTX    "PARTX"
#define PARTY    "PARTY"
#define INPUT_X0 "INPUT_X0"
#define INPUT_Y0 "INPUT_Y0"
#define PART_RA  "PART_RA"
#define PART_DEC "PART_DEC"
#define PARTMAG  "PARTMAG"
#define AVERAGEKEY "AVERAGE"
#define RMSKEY     "RMS" 
#define COLLMODE  "COLLMODE"

// photometry :
#define SIGMAL   "SIGMAL"
#define MEANL    "MEANL"
#define TRESHL   "TRESHL"
#define MAGLIMIT "MAGLIMIT"
#define FITA     "FITA"
#define FITB     "FITB"
#define FITCHI2  "FITCHI2"
#define FITSTARS "FITSTARS"
#define NSTARS   "NSTARS"

// astrometry :
#define POSANGLE "POSANGLE"
#define AST_ORD  "AST_ORD"
#define PARX     "PAR_X_"
#define PARY     "PAR_Y_"
#define PIXSCALE_AST "ASTPIXSC"
#define AST_UTTIME "ASTUTIME"
#define OBSMODE "OBSMODE"
#define ASTROOK "ASTROOK"
#define AST_ERR "AST_ERR"
#define AST_VER_KEYWORD "AST_VER"
#define N_APERT "N_APERT"
#define RA2000  "RA2000"
#define DEC2000  "DEC2000"

// frame quality :
#define QUALITY "QUALITY"
#define PIXEL_SHIFT_FIXED "REPAIRED_PIXEL_SHIFT"

// as input for function strptime
#define DATE_OBS_FMT  "'%Y-%m-%dT%H:%M:%S"
#define DATE_FMT "'%Y-%m-%d'"
#define TIME_FMT "'%H:%M:%S'"


// our private headers :
#define EVENT_X_ON_SMALL    "EVTXS"
#define EVENT_Y_ON_SMALL    "EVTYS"
#define EVENT_ORGINAL_X     "EVTOX"
#define EVENT_ORGINAL_Y     "EVTOY"
#define EVENT_CURRENT_FRAME "EVTCF"
#define EVENT_FRAME         "EVTFR"
#define TAKEN_AT_X          "EVTX0"
#define TAKEN_AT_Y          "EVTY0"
#define EVTX0               "EVTX0"
#define EVTY0               "EVTY0"
#define AVSTART             "AVSTART"
#define AVEND               "AVEND"
#define AVERF               "AVERF"
#define SHUTMODE            "SHUTMODE"

// extern const char* all_header_list[];

enum eFITSCompressionType { eFITSComprNone=0, eFITSComprASAS };

enum eHeaderKeyType { eSTRING=0, eFLOAT, eLONG, eEND, eUNKNOWN, eBOOL };

enum eHeaderSectionName { eSecUnknown=0, eSecStandard, eSecDataFormat, 
								  eSecObject, 
								  eSecObsSite, eSecInstrument, 
								  eSecExpID, eSecExpSettings, eSecExpEnv,
								  eSecExpDT, eSecAstrometry,
								  eSecPhotometry, eSecEvent,
								  eSecMount,
								  eSecOther };

extern eHeaderSectionName sectionList[];

struct sFITSHeaderDesc
{
	eHeaderKeyType keyType;
	const char* szKeyName;
	const char* szKeyComment;
	eHeaderSectionName eSection;
};

enum eFITSHeaderKey  { eUT_START=0, eDATE_OBS, eEXPTIME, eCCDPICNO,
				 			  eMPP_BC, eRO_TIME, eSHUTTER, eFOCUS, eELECGAIN,
							  eCHIP_TEMP, eCASE_TEMP, eAMBT_TEMP, eCOOLING,
							  eABINN, eSBINN, eELECOFFSET, eSPEED, eSAVEAREA, 
							  eORIGIN, eOBSERVER, eSITE, eTELLONG, eTELLAT,
							  eTELALT, eINSTRUME, eCAMERA, eCAMOPTIC, ePIXSCALE, 
		 					  eMIPS_HI, eMIPS_LO, eCONTRAST, eMEAN, eSIGMA, eDATAMAX, 
		 					  eDATAMIN, eBGMEAN, eBGSIGMA, eEVENT_FRAME, 
							  eEVENT_X_ON_SMALL, eEVENT_Y_ON_SMALL, eEVENT_ORGINAL_X, 
					  		  eEVENT_ORGINAL_Y, eEVENT_CURRENT_FRAME ,
					  		  eBZERO,eBSCALE,eBITPIX,
					  		  eNAXIS1, eNAXIS2, eTAKEN_AT_X, eTAKEN_AT_Y,
					  		  eRA_OBS,eDEC_OBS,eALT_OBS,eAZIM_OBS,eHA_OBS,
					  		  eTIME_OBS,eUT_TIME,
					  		  eNIMAGE,eJD,eHJD,eST,eAIRMASS,eEPOCH,
					  		  eZENITH_D,eFILENAME,eFILTER,eSOFTWARE,
					  		  eSOFTBUILD, eSHTIME, ePIXSIZE, eADCGAIN, eADCBIAS,
					  		  eSPEEDSET, eSPEEDMHZ, eCAMHUMID, eAMBHUMID,
					  		  eDRVTYPE, eROTATE, eOBJECT, eDATE_OBS2, eNAXIS3,
					  		  eEQUINOX, eRNOISE, eSIMPLE, eNAXIS, eEXTEND,
					  		  eCHISTEMP, eADCBSET, eADCGSET,	eFOCALLEN, eADCRANGE,
					  		  eADCCLAMPING, eUSBMODE, eFPGAVER, eCPRSVER,
					  		  eCAMIIDX, eCAMID, eINTRTEMP, eREXPTIME, 
					  		  eLNAGAIN, eCROTIME, eLOCTIME, eLOCDATE,
					  		  eVERDESC,eHITLENS,eRELNOISE,eDIMAGE,eFLIP,
					  		  ePOSANGLE, eASTROOK,
					  		  eAST_ORD, ePARX, ePARY, ePIXSCALE_AST,
					  		  eCCD_TRANSFORM_UT_TIME, eOBSMODE,
					  		  eCOMPRES, eBLOCKSZ, eDIVISOR,
							  eUT_END, eDATE_END, eCHIP_TEM, eORGSIZEX, eORGSIZEY,
							  eMOUNTRA, eMOUNTDEC, eMOUNTAZIM, eMOUNTALT, eMOUNTTRACK, eMOUNTHA,
							  eMOUNTTM, eDOME, eAVSTART, eAVEND,
							  ePARTNAME, ePARTX, ePARTY, eINPUT_X0, eINPUT_Y0,
							  eTAKENFROM, ePART_RA, ePART_DEC, ePART_MAG,
							  eSIGMAL, eMEANL, eTRESHL, eMAGLIMIT, eFITA, eFITB, eFITCHI2, eFITSTARS,
							  eNSTARS, eASTERR, eSHUTMODE, eDRVFLIP,
							  eAVERAGEKEY, eRMSKEY, eNIGHT,
					  		  eENDOfKEYLIST };

extern sFITSHeaderDesc fitsHeaderKeyDefTab[];
extern const char* szStandardHeaderKeyTab[];

// functions :
const char* GetKeyComment( const char* szKeyName );
sFITSHeaderDesc* GetKeyDesc( const char* szKeyName );
sFITSHeaderDesc* GetKeyDesc( const char* szKeyName  , int& start_pos );
eHeaderKeyType GetKeyType( const char* szKeyName );
int IsStandardKey( const char* key );
const char* GetSectionName( eHeaderSectionName eSection );

#endif
