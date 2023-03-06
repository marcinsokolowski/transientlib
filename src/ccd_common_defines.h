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
#ifndef _CCD_COMMON_DEFINES_H__
#define _CCD_COMMON_DEFINES_H__ 

#include<basedefines.h>

enum eActionOnTempWaitFail { eTempWaitFail_IGNORE=0, eTempWaitFail_SETFINAL };

enum eRunType { eRunOnlineCoic=0, eRunOfflineSum=1, eRunOfflineSN=2, 
				 	 eRunOfflineSingle=3, eRunOnlineSingle=4, eOfflineCoic=5,
					 eRunGetPart=10 };


enum eCheckIfMorePoint_T { eAfterCurrentFrameOnly=0, eAfterCoic, eBoth };

enum eCCDTYPE_T { e784x512_Audine=0, e2048x2062_Our, eFileSimulator, 
						e2048x2062_Version2, eDevice2K2KGigaEthType };

// oberving modes :
enum eObservationMode_T { eNoMovingMode=0,     // CONST MOUNT POSITION 
								  eEarthMovingMode };  // MOUNT ROTATES WITH EARTH 

extern const char* eObservationMode_DescTab[];
const char* GetObsModeDesc( eObservationMode_T obsMode );
eObservationMode_T GetObsMode( const char* szMode );

// #define ELEM_TYPE short
// #define BIG_ELEM_TYPE LONGLONG_T
// #define BIG_ELEM_TYPE long

#define ERROR_VALUE -1


// initialization files and datetime :
#define SAVE_DATA_DIR "$(HOME)/.pisoft/"
#define DEFAULT_INI_FILE "$(HOME)/.pisoft/ccdview.ini"
#define MOTOR_STEPS_SAVE_FILE "$(HOME)/.pisoft/motor%d.sav"
#define FRAMES_COUNTER_FILE "$(HOME)/.pisoft/fcounter%d.sav"
#define EVENT_COUNTER_FILE "$(HOME)/.pisoft/eventno.sav"

// event log types ( evt_type ) :
// defines moved from ccd_report.h :
enum eLogFileFormat { eVerif=0, eFinal, eUnknownLogFormat, eAllEventsLog, eGenEventsLog };

#endif
