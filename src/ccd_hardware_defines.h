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
#ifndef _CCD_HARDWARE_DEFINES_H__
#define _CCD_HARDWARE_DEFINES_H__

#include "mytypes.h"

// ENUMERATORS :
enum eADCRange { e2V=0, e4V };

const char* GetClampingDesc( BOOL_T bOnOff );

const char* GetADCRangeDesc( eADCRange adcRange );

enum eLNAGainValue_T { eLNA_Gain_x8=1, eLNA_Gain_x20 };

enum eMPP_T { eBC=0, eMPP };

enum eUSBMODE_T { eUSB_11=0, eUSB_20 };

// description functions
const char* GetUsbModeDesc( eUSBMODE_T  usbmode );


// ccd chip constatnts :


// default values :
// ADC :
#define DEFAULT_ADC_GAIN   63
#define DEFAULT_ADC_OFFSET 0
#define DEFAULT_ADC_CLAMPING TRUE
#define DEFAULT_ADC_RANGE e2V
#define DEFAULT_LNA_GAIN eLNA_Gain_x8
#define DEFAULT_ADC_CONF_REG0_SETTING 80
#define DEFAULT_ADC_CONF_REG1_SETTING 192

// READOUT SPEED :
#define DEFAULT_READOUT_SPEED_HORIZONTAL 1
#define DEFAULT_READOUT_SPEED_VERTICAL 1
#define DEFAULT_MPP_BC eBC

// cooling
#define DEFAULT_COOLING_ON_OFF 1
#define DEFAULT_TEMPERATURE 5

// TIMEOUT :
#define DEFAULT_DRIVER_ITER_TIMEOUT 5000
#define DEFAULT_RETRY_COUNT 5


#endif

