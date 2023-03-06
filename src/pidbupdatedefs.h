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
#ifndef _PIDB_UPDATE_DEFS_H__
#define _PIDB_UPDATE_DEFS_H__


enum eInstrumentsID {
	eUNKNOWN_INTRUMENT=0,
	eHETE_WXM=100,
	eHETE_FREGATE=101,
	eHETE_SXC=102,		
	eSPIACS_INSTRUMENT=201,
	eINTEGRAL_IBIS=202,
	eINTEGRAL_ISGRI=203,
	eKONUS_INSTRUMENT=400,
	eSWIFT_BAT_INSTRUMENT=501,
	eSWIFT_XRT=502,
	eSWIFT_UVOT=503
};


#define GCN_SOURCE_COUNT 11
enum eSourceDef { eUNKNOWN_SOURCE=0, eHETE=1, eINTEGRAL=2, eIPN=3, eKONUS=4, 
						eSWIFT=5, eSOURCE_EMPTY=6, eSAXWFC=7, eASM=8, eCGRO=9
					 };
extern const char* SourceIDDefTab[GCN_SOURCE_COUNT];

enum eAlertValidity { eValidityUnknown=0, 
							 ePossibleGRB=1, 
							 eDefiniteGRB=2, 
							 eShortSpike=3, 
							 eSolarActivity=4,
							 eValNoise=5 };
							 
							 
eAlertValidity GetValidityFlag( const char* szValidity );

enum eAlertStatus { eUNKNOWN_STATUS=0,
						  ePROBABLE=1,
						  eDEFINITE=2,
						  eTEST=3,
						  eNOT_GRB=4,
						  eFLARE=5,
						  eSOLAR_FLARE=6,
						  eSGR=7,
						  ePARTICLE=8,
						  eNOISE=9,
						  eCRAB=10,
						  eKNOWN_SOURCE=11,
						  eXRF=12 };
						  
eAlertStatus GetAlertStatus( const char* szStatus );						  
						  
int IsTestTrigger( int alert_type );

int GetInstrumentID( const char* szSourceID );

eSourceDef GetSourceID( const char* szSourceID );

enum eNoticeType { eNoticeUndefined=0,
                   eNoticeBatseOriginal=1,
                   eNoticeTest=2,
                   eNoticeImalive=3,
                   eNoticeKill=4 };
                                                                            

struct cNoticeTypeDef
{
	int type;
	const char* type_desc;
	const char* type_cmt;
	const char* full_desc;
};

int GetNoticeType( const char* szType );

extern cNoticeTypeDef NoticeDefTab[];

// convert functions :
int DegToDBDeg( double in_deg );
int MinToDBMin( double in_min );
int SecToDBSec( double in_sec );
double DBDegToDeg( int in_db_deg );
const char* IntToDBBool( int bVal );

// string DblToStr(double m_dbl);

// string BoolToStr(bool m_bool);

// SN table defines :
#define DISCOVERER_SIZE 64

#endif
