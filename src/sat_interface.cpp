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
#include "sat_interface.h"
#include <myfile.h>
#include <mystrtable.h>
#include <myparser.h>
#include <mydate.h>
#include <mymacros.h>

void (*CSatInterface::m_PrintError)( const char* szError ) = NULL;
mystring CSatInterface::m_szHttpAddress;
time_t   CSatInterface::m_SatLastUpdateDTM=0;
int      CSatInterface::m_IgnoreOlderThen=(24*3600*10); // 10 days older are ignored

CSatInterface::CSatInterface( const char* szName, double ra, double dec, eObjectType_T objtype,
										int pos_valid_time, double fov, int skip_older /*=86400*/  )
: m_szSatName( szName ) , m_RA(ra), m_DEC(dec), m_eObjType( objtype ), 
  m_PosValidTime(pos_valid_time),m_FOV(fov),m_SkipOlder(skip_older)
{

}

CSatInterface::CSatInterface()
: m_PosValidTime(1800), m_szSatName("UNKNOWN"), m_eObjType( eSatObj ), 
  m_RA(0), m_DEC(0), m_FOV(0.016),m_SkipOlder(24*3600)
{
}                                                                                
                                                                                
CSatInterface::CSatInterface( const char* szHttpAddress, 
										const char* szName,
										const char* szInfoFile,
										int pos_valid_time,
										double fov,
										const char* szInfoFileAlternative,
										int skip_older /*=86400*/ )
: m_PosValidTime(pos_valid_time), m_eObjType( eSatObj ), m_RA(0), m_DEC(0),
  m_FOV(fov),m_SkipOlder(skip_older)
{
	SetHttpAddress( szHttpAddress );
	if( szInfoFile && szInfoFile[0] ){
		SetInfoFile( szInfoFile );
	}
	if( szName && szName[0] ){
		m_szSatName = szName;
	}
	if( szInfoFileAlternative && szInfoFileAlternative[0] ){
		m_szSatInfoFileAlternative = szInfoFileAlternative;
	}
}

void CSatInterface::SetInfoFile( const char* szInfoFile )
{
	m_szSatInfoFile = szInfoFile;
}

void CSatInterface::SetHttpAddress( const char* szHttpAddress )
{
	m_szHttpAddress = szHttpAddress;
}

CSatInfo* CSatInterface::GetInfo( time_t ut_time, double& ra, double& dec,
	                      	       time_t& at_time, 
											 time_t& track_time, 
											 time_t& start_time, time_t& end_time )
{
	CSatInfo* pPrev=NULL;
	for(int i=0;i<m_SatInfoTab.size();i++){
		CSatInfo& info = m_SatInfoTab[i];
		if( ut_time>=info.m_StartTime && ut_time<=info.m_EndTime ){
			ra = info.m_RA_in_deg;
			dec = info.m_DEC_in_deg;
			at_time = ut_time;						
			track_time = ( info.m_EndTime - ut_time );
			start_time = info.m_StartTime;
			end_time = info.m_EndTime;
			return &(m_SatInfoTab[i]);
		}	
		if( pPrev ){
			if( ut_time>=pPrev->m_EndTime && ut_time<=info.m_StartTime ){
				ra = info.m_RA_in_deg;
	         dec = info.m_DEC_in_deg;
   	      at_time = ut_time;
				track_time = ( info.m_EndTime - ut_time );
				start_time = ut_time;
	         end_time = info.m_EndTime;
      	   return &(m_SatInfoTab[i]);
			}
		}
		pPrev = &info;
	}
	CSatInfo* pRet=NULL;
	if( m_SatInfoTab.size()>0 ){
		CSatInfo& info = m_SatInfoTab[0];
		if( ut_time<info.m_StartTime ){
			pRet = &info;
			ra = info.m_RA_in_deg;
			dec = info.m_DEC_in_deg;
			at_time = info.m_StartTime;
			track_time = ( info.m_EndTime - ut_time );
			start_time = ut_time;
         end_time = ut_time+3600*100;
		}

		info = m_SatInfoTab[ m_SatInfoTab.size()-1 ];
		if( ut_time>info.m_StartTime ){
			pRet = &info;
         ra = info.m_RA_in_deg;
         dec = info.m_DEC_in_deg;
         at_time = info.m_StartTime;
			start_time = info.m_StartTime;         	
			
// change on 20061017 - due to fact that last swift position 
// was used as good for long time and in fact it is not for sure !
//			track_time = 360000; // + inf
//         end_time = info.m_StartTime+3600*100;
			track_time = 7200; // 2 hours validity for last position is assumed
			end_time   = info.m_StartTime+7200; // 2 hours validity for last position is assumed			
		}
	}

	if( end_time > ut_time ){
		return pRet;
	}else{
		printf("Command should end at time=%d < curr_time=%d\n",end_time,ut_time);
		printf("Pointing info for time=%d is obsolate, not used\n",ut_time);
		return NULL;
	}	
}

BOOL_T CSatInterface::GetInfo( time_t ut_time, double& ra, double& dec, time_t& at_time )
{
	if( m_eObjType == eAstroObj ){
		// this handles constant astrophysical objects to be observed :
		at_time = ut_time;
		ra = m_RA;
		dec = m_DEC;
		return TRUE;
	}

	CSatInfo* pPrev=NULL;
	for(int i=0;i<m_SatInfoTab.size();i++){
		CSatInfo& info = m_SatInfoTab[i];
		if( ut_time>=info.m_StartTime && ut_time<=info.m_EndTime ){
			ra = info.m_RA_in_deg;
			dec = info.m_DEC_in_deg;
			at_time = ut_time;						
			return TRUE;
		}	
		if( pPrev ){
			if( ut_time>=pPrev->m_EndTime && ut_time<=info.m_StartTime ){
				ra = info.m_RA_in_deg;
	         dec = info.m_DEC_in_deg;
   	      at_time = ut_time;
      	   return TRUE;
			}
		}
		pPrev = &info;
	}
	BOOL_T bRet=FALSE;
	if( m_SatInfoTab.size()>0 ){
		CSatInfo& info = m_SatInfoTab[0];
		if( ut_time<info.m_StartTime ){
			bRet = TRUE;
			ra = info.m_RA_in_deg;
			dec = info.m_DEC_in_deg;
			at_time = info.m_StartTime;
		}

		info = m_SatInfoTab[ m_SatInfoTab.size()-1 ];
		if( ut_time>info.m_StartTime ){
			bRet = TRUE;
         ra = info.m_RA_in_deg;
         dec = info.m_DEC_in_deg;
         at_time = info.m_StartTime;
		}
	}

	return bRet;
}

BOOL_T CSatInterface::GetRecentInfo( time_t ut_time, double& ra, double& dec, 
												 time_t& at_time, int allow_before )
{
	if( m_eObjType == eAstroObj ){
		// this handles constant astrophysical objects to be observed :
		at_time = ut_time;
		ra = m_RA;
		dec = m_DEC;
		return TRUE;
	}

	CSatInfo* pPrev=NULL;
	for(int i=(m_SatInfoTab.size()-1);i>=0;i--){
		CSatInfo& info = m_SatInfoTab[i];
		if( ut_time>=(info.m_StartTime-allow_before) && ut_time<=info.m_EndTime ){
			ra = info.m_RA_in_deg;
			dec = info.m_DEC_in_deg;
			at_time = ut_time;						
			return TRUE;
		}	
		if( pPrev ){
			if( ut_time>=info.m_EndTime && ut_time<=pPrev->m_StartTime ){
				ra = pPrev->m_RA_in_deg;
	         dec = pPrev->m_DEC_in_deg;
   	      at_time = ut_time;
      	   return TRUE;
			}
		}
		pPrev = &info;
	}
	BOOL_T bRet=FALSE;
	if( m_SatInfoTab.size()>0 ){
		CSatInfo& info = m_SatInfoTab[0];
		if( ut_time<info.m_StartTime ){
			bRet = TRUE;
			ra = info.m_RA_in_deg;
			dec = info.m_DEC_in_deg;
			at_time = info.m_StartTime;
		}

		info = m_SatInfoTab[ m_SatInfoTab.size()-1 ];
		if( ut_time>info.m_StartTime ){
			bRet = TRUE;
         ra = info.m_RA_in_deg;
         dec = info.m_DEC_in_deg;
         at_time = info.m_StartTime;
		}
	}

	return bRet;
}


BOOL_T CSatInterface::ParseInfoFileAll( int skip_older )
{
	printf("DEBUG : skip_older = %d\n",skip_older);

	BOOL_T bRet=FALSE;
	if( !ParseInfoFile( skip_older ) ){
		printf("ERROR : could not parse pointing file : %s\n",m_szSatInfoFile.c_str());
		
		if( strlen( m_szSatInfoFileAlternative.c_str() )>0 ){
			printf("Retrying with file %s\n",m_szSatInfoFileAlternative.c_str());
			bRet = ParseInfoFile( skip_older , m_szSatInfoFileAlternative.c_str() );
			if( !bRet ){
				printf("ERROR : could not parse alternative pointing file %s\n",m_szSatInfoFileAlternative.c_str() );
			}			
		}else{
			printf("WARNING : no alternative poitning file provided for target %s\n",m_szSatName.c_str());
		}
	}else{
		bRet = TRUE;
	}
	
	return bRet;
}

BOOL_T CSatInterface::ParseInfoFile( int skip_older, const char* info_file )
{
	const char* pInfoFile = info_file;
	if( !pInfoFile ){
		pInfoFile = m_szSatInfoFile.c_str();
	}

	if( m_eObjType == eAstroObj ){
		return TRUE;
	}

	time_t curr_time=get_dttm();

	if( strlen( pInfoFile ) ){
		if( MyFile::DoesFileExist( pInfoFile ) ){
			MyIFile in( pInfoFile );
			const char* pLine;
			m_SatInfoTab.clear();
			int i=0;
			while( pLine = in.GetLine( TRUE )){	
				if(pLine && pLine[0]){
					i++;
					if( mystring::get_first_non_white( pLine )!='#' ){						
						CMyStrTable items;
						MyParser pars = pLine;
						pars.GetItems( items );
						if( items.size()>=5 ){
							CSatInfo info;
							info.m_szInfoName = "POSINFO";
							info.m_RA_in_deg = atof( items[3].c_str() );
							info.m_DEC_in_deg = atof( items[4].c_str() );

							mystring szStart = items[0].c_str();								
							info.m_StartTime = atol( szStart.c_str() );

							if( skip_older>0 ){
								if( (curr_time - info.m_StartTime) > skip_older ){
									continue;
								}
							}

							info.m_EndTime = 0;
							if( m_SatInfoTab.size() ){
								m_SatInfoTab.back().m_EndTime = info.m_StartTime;
							}
							CSatInterface::m_SatLastUpdateDTM = info.m_StartTime;							


							if( info.m_StartTime>0 ){
								m_SatInfoTab.push_back( info );
							}else{
								printf("ERROR in line %d : bad information in pointing file, unix_time=%d < 0 , (ra,dec)=(%.2f,%.2f), skiped\n",i,info.m_StartTime,info.m_RA_in_deg,info.m_DEC_in_deg);
							}
						}
					}
				}
			}
			
			if( m_SatInfoTab.size() && m_SatInfoTab.back().m_EndTime<10  ){
				m_SatInfoTab.back().m_EndTime = m_SatInfoTab.back().m_StartTime + (24*3600);
			}

			for(int i=0;i<m_SatInfoTab.size();i++){
				CSatInfo& info = m_SatInfoTab[i];

				_TRACE_PRINTF_4("%s %d %d %.5f %.5f\n",info.m_szInfoName.c_str(),
							info.m_StartTime,info.m_EndTime,
							info.m_RA_in_deg,
							info.m_DEC_in_deg);
			}

			printf("Found %d positions for target %s ( skip_older = %d sec)\n",m_SatInfoTab.size(),
					m_szSatName.c_str(),skip_older);fflush(stdout);
			return ( m_SatInfoTab.size() > 0 );
		}
	}else{
		printf("No satellite information file was specified, cannot continue\n");
		return FALSE;
	}		
	return FALSE;
}


BOOL_T CSatInterface::Wget()
{
	mystring szCmd;
	szCmd << "wget --cache=off -q -O " << m_szSatInfoFile.c_str()
			<< " " << m_szHttpAddress;
	printf("running command :\n");
	printf("%s\n",szCmd.c_str());
	int ret = system( szCmd.c_str() );

	if(ret<0){
		printf("wget - FAILED !\n");
	}

	return (ret>=0);
}

BOOL_T CSatInterface::UpdateSatInfo( BOOL_T bWget )
{
	if( bWget ){
		if( !Wget() ){
			printf("Wget failed ( file : %s )\n",m_szSatInfoFile.c_str());
			return FALSE;
		}
	}
	return ParseInfoFile();
}

CSatInfo* CSatInterface::FindByTargetNumber( int traget_num )
{
	int size_1 = m_SatInfoTab.size()-1;
	for(int i=size_1;i>=0;i--){
		if( m_SatInfoTab[i].m_TargetNumber == traget_num ){
			return &(m_SatInfoTab[i]);
		}
	}
	return NULL;
}
