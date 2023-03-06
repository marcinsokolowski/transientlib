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
#ifndef _CCD_PIMAN_INTERFACE_DEFS_H__
#define _CCD_PIMAN_INTERFACE_DEFS_H__

// CORBA DEFINES :
// in order to define it use ccd.cfg variables :
// CCD_CORBA_OPTIONS, CCD_CORBA_WITH_NAME_SERVICE, 
// and environement variable CORBA_NAME_SRVICE_PORT - will overwrite default CORBA_NAME_SERVER_PORTNO


// OLD :
// #define CORBA_NAME_SERVER_OPTIONS_SRV "corba_srv -ORBInitRef NameService=corbaloc::localhost.localdomain:%d/NameService"
// #define CORBA_NAME_SERVER_OPTIONS_CLN "corba_cln -ORBInitRef NameService=corbaloc::localhost.localdomain:%d/NameService"
// NEW :
#define CORBA_NAME_SERVER_OPTIONS_SRV "corba_srv -ORBInitRef NameService=corbaloc::%s:%d/NameService"
#define CORBA_NAME_SERVER_OPTIONS_CLN "corba_cln -ORBInitRef NameService=corbaloc::%s:%d/NameService"


#define CORBA_NORMAL_SRV "corba_srv -ORBIIOPAddr inet:localhost.localdomain:12123 -ORBNoCodeSets -ORBIIOPBlocking"
#define CORBA_NAME_SERVER_PORTNO 12456


// request :
// COMMAND.PIPELINE_INDEX.CCD_INDEX.PAR1.PAR2. ...
// COMMAND LIST

// type of request : 
// eCCDR_Global - can be executed globaly no special actions for 
//                single cameras is required
// eCCDR_Cam    - request for camera required 
enum eCCDRequest_T { eCCDR_Global=0, eCCDR_Cam };

// request type enumerator :
enum eCCDRequestType_T { eDAQReq_NotDef=0, 
								 eCCDReq_ParamChange,
                         eDAQReq_SatTriggerMode,
                         eDAQReq_MovingMove, 
                         eDAQReq_NormalMode,
                         eDAQReq_Reset, 
                         eDAQReq_EXIT,
                         eDAQReq_SaveCalibPicure, 
                         eDAQReq_SaveCurrentPicture, 
                         eCCDDRV_SetTemp, 
                         eCCDDRV_GetTemp, 
                         eCCDDRV_CoolingOnOff,
                         eCCDDRV_EXIT, 
                         eCCDDRV_GetFrame, 
                         eCCDReq_GetParam, 
                         eDAQReq_OnTriggerPosition,
                         eDAQReq_ChangeMountPosition, 
                         eDAQReq_StartAnalysis,
                         eDAQReq_StopAnalysis,
                         eDAQReq_DoDarks, 
                         eDAQReq_TakeNPictures,
                         eDAQReq_LoadParamFile,
                         eDAQReq_DoAstrometryNow,
                         eDAQReq_DoAstrometryMode,
                         eDAQReq_CorrectMountPosition,
                         eDAQReq_SatTriggerCancel
                      };

// common :
#define CCD_CHANGE_PARAM "CCD_CHANGE_PARAM"
#define CCD_GET_PARAM    "CCD_GET_PARAM"

#define DAQ_RESET               "DAQ_RESET"
#define DAQ_EXIT                "DAQ_EXIT"
#define DAQ_SAT_TRIGGER         "DAQ_SAT_TRIGGER"
#define DAQ_SEND_ALERT          "DAQ_SEND_ALERT"
#define DAQ_NORMAL_MODE         "DAQ_NORMAL_MODE"
#define DAQ_MOVING_MODE         "DAQ_MOVING_MODE"
#define DAQ_ON_TRIGGER_POSITION "DAQ_ON_TRIGGER_POSITION"
#define DAQ_SAVE_CALIB          "DAQ_SAVE_CALIB"
#define DAQ_SAVE_CURRENT        "DAQ_SAVE_CURRENT"
#define DAQ_GET_SHUTTER_TIME    "DAQ_GET_SHUTTER_TIME"
#define DAQ_SET_SHUTTER_TIME    "DAQ_SET_SHUTTER_TIME"
#define DAQ_GET_CURRENT_FRAME   "DAQ_GET_CURRENT_FRAME"
#define DAQ_GET_PICTURE_SIZE    "DAQ_GET_PICTURE_SIZE"
#define DAQ_SET_MOUNT_POSITION  "DAQ_SET_MOUNT_POSITION"
#define DAQ_START_ANALYSIS      "DAQ_START_ANALYSIS"
#define DAQ_STOP_ANALYSIS       "DAQ_STOP_ANALYSIS"
#define DAQ_DO_DARKS            "DAQ_DO_DARKS"
#define DAQ_TAKE_N_PICTURES     "DAQ_TAKE_N_PICTURES"
#define DAQ_GET_POSITION        "DAQ_GET_POSITION"
#define DAQ_LOAD_PARAM_FILE     "DAQ_LOAD_PARAM_FILE"
#define DAQ_DO_ASTROMETRY_NOW   "DAQ_DO_ASTROMETRY_NOW"
#define DAQ_DO_ASTROMETRY_MODE  "DAQ_DO_ASTROMETRY_MODE"
#define DAQ_START               "DAQ_START"
#define DAQ_SET_CUSTOM_KEY      "DAQ_SET_CUSTOM_KEY"
#define CCD_COMMAND_SEPARATOR " "



// requests to DRIVER :
#define CCDDRV_EXIT         "CCDDRV_EXIT"
#define CCDDRV_SET_TEMP     "CCDDRV_SET_TEMP"
#define CCDDRV_GET_TEMP     "CCDDRV_GET_TEMP"
#define CCDDRV_SET_COOLING  "CCDDRV_SET_COOLING"
#define CCDDRV_GET_COOLING  "CCDDRV_GET_COOLING"
#define CCDDRV_GET_NEXT_FRAME "CCDDRV_GET_NEXT_FRAME"


#define DAQ_PORT 51002
#define DAQ_PORT_1 51002
#define DAQ_PORT_2 51003

#define DAQ_PI_SYS_PORT_NO 51001

#define DRV_PI_SYS_PORT_NO 41001

#define CCDCAM_GET_NEXT_FRAME "CCDCAM_GET_NEXT_FRAME"
#define ANALYZE_PICTURE "AnalyzePicture"
#define LAST_FRAME      "LastFrame"

#define HEADER_SIZE 1000

#define READY_FLAG "READY"
#define CONNECT_FLAG "DRIVER_CONNECTED"




struct CPortAssigment {
	int ccd_index;
	int port_no;
};

extern CPortAssigment gCCDPortList[]; 

#endif
