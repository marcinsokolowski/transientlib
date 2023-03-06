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
#ifndef _SAT_INTERFACE_H__
#define _SAT_INTERFACE_H__

#include <mystring.h>
#include <vector>

using namespace std;


enum eObjectType_T { eSatObj=0, eAstroObj };

struct CSatInfo
{
	mystring m_szInfoName;
	int m_TargetNumber;
	time_t m_StartTime;
	time_t m_EndTime;
	double m_RA_in_deg;
	double m_DEC_in_deg;
	mystring m_szComment;
};

class CSatInterface
{
public:
	static void (*m_PrintError)( const char* szError );
	static mystring m_szHttpAddress;
	static time_t   m_SatLastUpdateDTM;
	static int      m_IgnoreOlderThen;

	double m_FOV; // radius of object FOV in degrees , default = 1 arcmin 
	eObjectType_T m_eObjType;	
	mystring m_szSatInfoFile;
	mystring m_szSatInfoFileAlternative;
	mystring m_szSatName;
	int m_PosValidTime; // time of last position validity 
	int m_SkipOlder;    // skip older entries in pointing file

	// for star like objects with constant coordinates :
	double m_RA;  // right ascension [deg]
	double m_DEC; // declination [deg]
		
	vector<CSatInfo> m_SatInfoTab;
	CSatInfo* FindByTargetNumber( int traget_num );
	
	CSatInterface();
	CSatInterface( const char* szHttpAddress, const char* szName=NULL, 
						const char* szInfoFile=NULL, int pos_valid_time=1800,
						double fov=0.016, const char* szInfoFileAlternative=NULL,
						int skip_older=86400 );
	CSatInterface( const char* szName, double ra, double dec, eObjectType_T objtype,
						int pos_valid_time=1800, double fov=0.016, int skip_older=86400  );				

	static void SetHttpAddress( const char* szHttpAddress );
	void SetInfoFile( const char* szInfoFile );
	
	virtual BOOL_T GetInfo( time_t ut_time, double& ra, double& dec, 
									time_t& at_time );

	virtual BOOL_T GetRecentInfo( time_t ut_time, double& ra, double& dec, 
									time_t& at_time, int allow_before=0 );
									
	virtual CSatInfo* GetInfo( time_t ut_time, double& ra, double& dec,
									time_t& at_time, time_t& track_time,
									time_t& start_time, time_t& end_time );
									
	// by default skiping older then 1 day :
	virtual BOOL_T ParseInfoFile( int skip_older=86400, const char* info_file=NULL );
	virtual BOOL_T ParseInfoFileAll( int skip_older=86400 );	

	virtual BOOL_T Wget();
	
	virtual BOOL_T UpdateSatInfo( BOOL_T bWget );
};


#endif
