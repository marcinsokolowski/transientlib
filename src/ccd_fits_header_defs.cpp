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
#include <stdlib.h>
#include <string.h>
#include "ccd_fits_header_defs.h"

//const char* all_header_list[] = { DATE_OBS, UT_START, EXPOSURE, MIPS_HI, MIPS_LO,
//		CONTRAST, MEAN, SIGMA, DATAMAX, DATAMIN, BGMEAN, BGSIGMA, EVENT_FRAME, 
//		EVENT_X_ON_SMALL, EVENT_Y_ON_SMALL, EVENT_ORGINAL_X, EVENT_ORGINAL_Y, EVENT_CURRENT_FRAME ,NULL };



sFITSHeaderDesc fitsHeaderKeyDefTab[] = 
 {
	{ eSTRING, "UT-START", "UT OF START FROM PC", eSecExpDT },
	{ eSTRING, "DATE-OBS", "UT DATE AT START OF FRAME", eSecExpDT },
	{ eFLOAT, "EXPTIME", "TRUE EXPOSURE TIME IN SECONDS", eSecExpSettings },
	{ eLONG , "CCDPICNO", "FRAME NUMBER OF IMAGE", eSecExpID },
	{ eSTRING, MPP_BC, "'MPP' OR 'BC' MODE", eSecExpSettings },
	{ eFLOAT, ROTIME, "USB transfer time [sec]", eSecExpSettings },
	{ eSTRING, "SHUTTER", "'OPEN' OR 'CLOSED' (FOR DARKS)", eSecExpSettings },
	{ eFLOAT, "FOCUS", "POSITION OF FOCUS MOTOR IN STEPS", eSecExpSettings },
	{ eFLOAT, "ELECGAIN", "ELECTRONS PER ADU", eSecExpSettings },
	{ eFLOAT, CHIP_TEMP, "CCD TEMPERATURE *C", eSecExpEnv },
	{ eFLOAT, CASE_TEMP, "CAMERA CASE TEMPERATURE *C", eSecExpEnv },
	{ eFLOAT, AMBT_TEMP, "AMBIENT TEMPERATURE *C", eSecExpEnv },
	{ eSTRING, COOLING, "'ON' OR 'OFF'", eSecExpSettings },
	{ eSTRING, "ABINN", "ANALOG BINNING", eSecExpSettings },
	{ eSTRING, "SBINN", "SOFTWARE BINNING", eSecExpSettings },
	{ eLONG,   "ELECOFFSET", "ELEC OFFSET", eSecExpSettings },
	{ eLONG,   SPEED     , "READOUT SPEED", eSecExpSettings  },
	{ eSTRING, "SAVEAREA" , "X0,Y0,X1,Y1", eSecExpSettings },
	{ eSTRING, "ORIGIN" , "ORGANIZATION/INSTITUTION WHICH CREATED FITS", eSecObsSite },
	{ eSTRING, "OBSERVER", "OBSERVER NAME/IDENTIFICATION", eSecExpID },
	{ eSTRING, "SITE", "OBSERVING SITE", eSecObsSite },
	{ eFLOAT, "TELLONG", "TELESCOPE LONGITUDE (DEG)", eSecObsSite },
	{ eFLOAT, "TELLAT", "TELESCOPE LATITUDE (DEG)", eSecObsSite },
	{ eFLOAT, "TELALT", "TELESCOPE ALTITUDE (M)", eSecObsSite },
	{ eSTRING, "INSTRUME", "INSTRUMENT USED TO ACQUIRE THE DATA", eSecInstrument },
	{ eSTRING, "CAMERA", "", eSecInstrument },
	{ eSTRING, "CAMOPTIC", "", eSecInstrument },
	{ eFLOAT, "PIXSCALE", "ANGULAR SIZE OF PIXEL ( arcsec )", eSecInstrument },
	{ eSTRING, MIPS_HI, MIPS_HI, eSecExpSettings },
	{ eSTRING, MIPS_LO, MIPS_LO, eSecExpSettings },
	{ eLONG, CONTRAST, CONTRAST, eSecExpSettings },
	{ eFLOAT, MEAN, MEAN, eSecOther },
	{ eFLOAT, SIGMA, SIGMA, eSecOther },
	{ eFLOAT, DATAMAX, DATAMAX, eSecOther },
	{ eFLOAT, DATAMIN, DATAMIN, eSecOther },
	{ eFLOAT, BGMEAN, BGMEAN, eSecOther },
	{ eFLOAT, BGSIGMA, BGSIGMA, eSecOther },
	{ eLONG, EVENT_FRAME, "Event frame number", eSecEvent },
	{ eLONG, EVENT_X_ON_SMALL, "X position on saved frame", eSecEvent },
	{ eLONG, EVENT_Y_ON_SMALL, "Y position on saved frame", eSecEvent },
	{ eLONG, EVENT_ORGINAL_X, "X position on original frame", eSecEvent },
	{ eLONG, EVENT_ORGINAL_Y, "Y position on original frame", eSecEvent },
	{ eLONG, EVENT_CURRENT_FRAME, EVENT_CURRENT_FRAME, eSecEvent },
	{ eFLOAT, FH_BZERO, "offset data range to that of unsigned short", eSecStandard },
	{ eFLOAT, FH_BSCALE, "default scaling factor", eSecStandard },
	{ eLONG, FH_BITPIX, "number of bits per data pixel", eSecStandard },
	{ eLONG, "NAXIS1", "length of data axis 1", eSecStandard },
	{ eLONG, "NAXIS2", "length of data axis 2", eSecStandard },
	{ eLONG, TAKEN_AT_X, "taken at X from big image", eSecEvent },
	{ eLONG, TAKEN_AT_Y, "taken at Y from big image", eSecEvent },
	{ eFLOAT, RA_OBS, "RightAscension - observed", eSecObject },
	{ eFLOAT, DEC_OBS, "Declination - observed", eSecObject },
	{ eFLOAT, ALT_OBS, "Altitude - observed", eSecObject },
	{ eFLOAT, AZIM_OBS, "Azimuth - observed", eSecObject },
	{ eFLOAT, HA_OBS,  "Hour angle", eSecObject },
	{ eSTRING, TIME_OBS, "GMT at end of exposoure", eSecExpDT },
	{ eLONG,   TIME_UT, "time_t - sec since 0 UTC 1/1/1970 at start exposure", eSecExpDT },
	{ eLONG,   NIMAGE,  "sequential frame index", eSecExpID },
	{ eFLOAT,  JD, "Julian Date ad mid_exposure", eSecExpDT },
	{ eFLOAT,  HJD, "Heliocentric JD at mid_exposure", eSecExpDT },
	{ eFLOAT,  ST, "siderial time", eSecExpDT },
	{ eFLOAT, AIRMASS, "air mass @ end exposure", eSecExpEnv },
	{ eFLOAT, EPOCH, "IRAF equinox", eSecExpDT },
	{ eFLOAT, ZENITH_D, "zenith distance @ end exposure", eSecObject },
	{ eSTRING, FILENAME, "FITS file name", eSecExpID },
	{ eSTRING, FILTER, "Filter", eSecInstrument },
	{ eSTRING, SOFTVERSION, "Software Version", eSecExpID },
	{ eSTRING, SOFTBUILD, "Software Build", eSecExpID },
	{ eFLOAT,  SHTIME , "Shutter Time Parameter", eSecExpSettings },
	{ eLONG, PIXSIZE, "pixel size in microns", eSecInstrument },
	{ eFLOAT, ADCGAIN, "ADC gain value", eSecExpSettings },
	{ eFLOAT, ADCBIAS, "ADC offset [mV]", eSecExpSettings },
	{ eSTRING, SPEEDSET, "READOUT SPEED setting Vertical/Horizontal", eSecExpSettings },
	{ eSTRING, SPEEDMHZ, "READOUT SPEED MHz Vertical/Horizontal", eSecExpSettings },
	{ eFLOAT, CAMHUMID, "CAMERA HUMIDITY %", eSecExpEnv },
	{ eFLOAT, AMBHUMID, "AMBIENT HUMIDITY %", eSecExpEnv },
	{ eSTRING, DRVTYPE, "CCD Driver type", eSecExpID },
	{ eLONG, ROTATE, "S is UP ( rotated FOV )", eSecObject },
	{ eSTRING, OBJECT, "Object", eSecObject },
	{ eSTRING, DATE_OBS2, "GMT", eSecExpDT },
	{ eSTRING, NAXIS3, "just for compatibility", eSecStandard },
	{ eSTRING, EQUINOX, "", eSecExpDT },
	{ eFLOAT, RNOISE, "Readout noise ( +gauss , -rms )", eSecExpSettings },
	{ eBOOL, FH_SIMPLE, "file does conform to FITS standard", eSecStandard },
	{ eLONG, FH_NAXIS, "number of data axes", eSecStandard },
	{ eBOOL, FH_EXTEND, "FITS dataset may contain extensions", eSecStandard },
	{ eFLOAT, CHISTEMP, "Chip temperature parameter", eSecExpEnv },
	{ eFLOAT, ADCBSET, "ADC offset setting", eSecExpSettings },
	{ eFLOAT, ADCGSET, "ADC gain setting", eSecExpSettings },
	{ eFLOAT, FOCALLEN, "focal length in mm", eSecInstrument },
	{ eFLOAT, ADCRANGE, "ADC range [V]", eSecExpSettings },
	{ eBOOL, ADCCLAMP, "ADC 4V clamping ON/OFF", eSecExpSettings },
	{ eSTRING, USBMODE, "USB mode", eSecExpSettings },
	{ eSTRING, FPGAVER, "FPGA soft version", eSecExpSettings },
	{ eSTRING, CPRSVER, "CPRS soft version", eSecExpSettings },
	{ eLONG,   CAMIIDX, "Camera internal index", eSecExpID },
	{ eSTRING, CAMID,   "Camera hardware ID", eSecExpID },
	{ eFLOAT,  INTRTEMP,"Internal temperature ( chamber )", eSecExpEnv },
	{ eFLOAT, REXPTIME, "Measured exposition time", eSecExpSettings },
	{ eSTRING, LNAGAIN, "LNA gain - preamplifier gain", eSecExpSettings},
	{ eFLOAT, CROTIME, "Chip read time [sec]", eSecExpSettings },
	{ eSTRING, LOCTIME, "Local time @ frame start", eSecExpDT },
	{ eSTRING, LOCDATE, "Local date @ frame start", eSecExpDT },
	{ eSTRING, VERDESC, "Camera software version desc", eSecExpSettings },
	{ eSTRING, HITLENS, "Lens hitting on/off", eSecExpSettings },
	{ eFLOAT, RELNOISE, "Readout noise in electrons", eSecExpSettings },
	{ eLONG, DIMAGE, "Current night image number", eSecExpID },
	{ eSTRING, FLIP, "Image flip FH-horiz, FV-vert" , eSecObject },
	{ eFLOAT, POSANGLE, "Astrometry FI angle", eSecAstrometry },
   { eLONG, ASTROOK,  "AstroOK 0-failed,1-ok", eSecAstrometry },
	{ eLONG, AST_ORD, "Astrometry transform order", eSecAstrometry },
	{ eFLOAT, PARX, "Astrometry param X", eSecAstrometry },
	{ eFLOAT, PARY, "Astrometry param Y", eSecAstrometry },
	{ eFLOAT, PIXSCALE_AST, "Astrometry pixscale", eSecAstrometry },
	{ eLONG,  AST_UTTIME, "Astrometry ut-time", eSecAstrometry },
	{ eLONG,  OBSMODE, "Obs mode  0-const, 1-tracking", eSecObject },
	{ eSTRING, COMPRES, "Compression type", eSecDataFormat },
	{ eLONG, BLOCKSZ, "Block size", eSecDataFormat },
	{ eLONG, FH_DIVISOR, "Divisor", eSecDataFormat },
	{ eSTRING, UT_END, "UT @ end of exposure", eSecExpDT },
	{ eSTRING, DATE_END, "UT date @ end of exposure", eSecExpDT },
	{ eFLOAT, CHIP_TEM, "Chip temperature - compatiblity", eSecExpEnv },
	{ eLONG, ORGSIZEX, "Original frame X size", eSecEvent },
	{ eLONG, ORGSIZEY, "Original frame Y size", eSecEvent },
	{ eFLOAT, MOUNTRA, "RA from mount", eSecMount },
	{ eFLOAT, MOUNTDEC, "DEC from mount", eSecMount },
	{ eFLOAT, MOUNTAZIM, "AZIM from mount", eSecMount },
   { eFLOAT, MOUNTDEC, "ALT from mount", eSecMount },
	{ eFLOAT, MOUNTTRK, "Mount tracking Yes/No", eSecMount },
	{ eFLOAT, MOUNTHA, "HA from mount", eSecMount },
	{ eSTRING, MOUNTTM, "Mount position time stamp" , eSecMount },
	{ eSTRING, DOME, "Dome status", eSecExpEnv },
	{ eLONG, AVSTART, "Average - start frame", eSecEvent },
	{ eLONG, AVEND, "Average - end   frame", eSecEvent },
	{ eSTRING, PARTNAME, "Part Name", eSecEvent },
	{ eLONG, PARTX, "Part X", eSecEvent },
	{ eLONG, PARTY, "Part Y", eSecEvent },
	{ eLONG, INPUT_X0, "x - low left", eSecEvent },
	{ eLONG, INPUT_Y0, "y - low left", eSecEvent },
	{ eSTRING, TAKENFROM, "from FITS", eSecEvent },
	{ eFLOAT, PART_RA, "Part RA [h]", eSecEvent },
	{ eFLOAT, PART_DEC, "Part RA [deg]", eSecEvent },
	{ eSTRING, PARTMAG, "Magnitudo", eSecEvent },
	{ eFLOAT, SIGMAL, "Sigma Laplace", eSecPhotometry },
	{ eFLOAT, MEANL , "Mean Laplace" , eSecPhotometry },
	{ eFLOAT, TRESHL, "Tresh Laplace", eSecPhotometry },
	{ eFLOAT, MAGLIMIT, "Mag Limit", eSecPhotometry },
	{ eFLOAT, FITA , "Fit a", eSecPhotometry },
	{ eFLOAT, FITB , "Fit b", eSecPhotometry },
	{ eFLOAT, FITCHI2 , "Fit chi2", eSecPhotometry },
	{ eLONG, FITSTARS , "Fit stars count", eSecPhotometry },
	{ eLONG, NSTARS   , "Frame stars", eSecPhotometry },
	{ eEND, NULL, NULL }
 };


sFITSHeaderDesc* GetKeyDesc( const char* szKeyName  , int& start_pos )
{
	// this are exceptions - not keeping whole PAR_X_1 etc ...
	if( strncmp( szKeyName, PARX, strlen( PARX ) )==0 ){
		return &(fitsHeaderKeyDefTab[ePARX]);
	}
	if( strncmp( szKeyName, PARY, strlen( PARY ) )==0 ){
		return &(fitsHeaderKeyDefTab[ePARY]);
	}

	int start_pos_sav=start_pos;
	if( start_pos<0 || start_pos>=eENDOfKEYLIST )
		start_pos = 0;
	while(  start_pos<eENDOfKEYLIST && fitsHeaderKeyDefTab[start_pos].keyType!=eEND ){
 		if(strcmp( fitsHeaderKeyDefTab[start_pos].szKeyName,szKeyName )==0){
			return &(fitsHeaderKeyDefTab[start_pos]);
		}
		start_pos++;
	}	
	if( start_pos_sav>0 ){
		start_pos = 0;
		return GetKeyDesc( szKeyName, start_pos );
	}

	return NULL;	
}

sFITSHeaderDesc* GetKeyDesc( const char* szKeyName )
{
	int p=0;
	return GetKeyDesc( szKeyName , p );		
}


const char* GetKeyComment( const char* szKeyName )
{
	int i=0;
	while( fitsHeaderKeyDefTab[i].keyType!=eEND ){
		if(strcmp(fitsHeaderKeyDefTab[i].szKeyName,szKeyName)==0)
			return fitsHeaderKeyDefTab[i].szKeyComment;
		i++;
	}	
	return "";
}

eHeaderKeyType GetKeyType( const char* szKeyName )
{
	int i=0;
   while( fitsHeaderKeyDefTab[i].keyType!=eEND ){
      if(strcmp(fitsHeaderKeyDefTab[i].szKeyName,szKeyName)==0)
         return fitsHeaderKeyDefTab[i].keyType;
		i++;
   }
   return eUNKNOWN;
}


const char* szStandardHeaderKeyTab[] = { FH_SIMPLE, FH_BITPIX, FH_NAXIS, FH_NAXIS1, FH_NAXIS2, FH_EXTEND, 
													  FH_BZERO, FH_BSCALE , NULL };

int IsStandardKey( const char* key )
{
	int i=0;
	while( szStandardHeaderKeyTab[i] ){
		if(strcmp( key, szStandardHeaderKeyTab[i])==0)
			return 1;
		i++;
	}
	return 0;
}


eHeaderSectionName sectionList[] = { 
								  eSecUnknown, 
								  eSecDataFormat , eSecObsSite, 
								  eSecInstrument, eSecObject,
                          eSecExpID, eSecExpDT, eSecExpSettings, eSecExpEnv, 
								  eSecStandard, eSecAstrometry,
								  eSecOther };


const char* GetSectionName( eHeaderSectionName eSection ){
	if(eSection==eSecDataFormat){
		return "Data Format";
	}
	if(eSection==eSecObsSite){
		return "Observing Site";
	}
	if(eSection==eSecInstrument){
		return "Instrument";
	}
	if(eSection==eSecObject){
		return "Object";
	}
	if(eSection==eSecExpID){
		return "Exposure id";
	}
	if(eSection==eSecExpDT){
		return "Exposure date/time";
	}
	if(eSection==eSecExpSettings){
		return "Exposure settings";
	}
	if(eSection==eSecExpEnv){
		return "Exposure environment";
	}
	if(eSection==eSecStandard){
		return "Standard Keys";
	}
	if(eSection==eSecAstrometry){
		return "Astrometry";
	}
	if(eSection==eSecPhotometry){
		return "Photometry";
	}
	if(eSection==eSecMount){
                return "Mount";
        }
	if( eSection==eSecEvent ){
		return "Event";
	}
	if(eSection==eSecOther){
		return "Others";
	}
	
	
	return "Others";
}
