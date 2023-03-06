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
#include "ccd_samples.h"
#include "ccd_globals.h"
#include "ccd_trace.h"
// #include "ccd_matrix.h"
#include <mystrtable.h>
#include <myparser.h>
#include <myutil.h>
#include <random.h>
#include <cexcp.h>
#include <fits_file.h>
#include "ccd_util.h"
#include <Astroangle.h>

// CDescFile CccdSamples::m_DescFile;

CMagSamples::CMagSamples()
: m_szMagnitude(""),m_pSamples(NULL),m_Count(0),m_bRead(FALSE)
{
	gCCDParams.InitParams();
	m_pSampleFilesTab = new CMyStrTable();
}

CMagSamples::~CMagSamples()
{
	if(m_pSampleFilesTab)
		delete m_pSampleFilesTab;
	if(m_pSamples)
		delete [] m_pSamples;
}



LONG_T CMagSamples::ReadSamples( const CEnvVar& sampleDesc, const char* sample_dir )
{
	LONG_T total_size=0;
	if(!m_bRead){
		m_bRead = TRUE;
		m_szMagnitude = sampleDesc.szName;
		m_SamplesList = sampleDesc;


		if(strstr(sampleDesc.szValue.c_str(),".fit")){
			printf("treating description as list of FITS files : %s\n",sampleDesc.szValue.c_str());
			MyParser pars = sampleDesc.szValue.c_str();
			pars.GetItems( *m_pSampleFilesTab );
		}else{
			printf("treating description as list file : %s\n",sampleDesc.szValue.c_str()); 
			mystring listfile;
			listfile << sample_dir << "/" << sampleDesc.szValue;
			CListFile magListFile( listfile.c_str() );
			(*m_pSampleFilesTab) = magListFile.GetListTable();
		}

		// allocating CCDMatrix table for samples :
//		m_pSamples = new Table2D<ELEM_SAMPLE_TYPE>[m_pSampleFilesTab->GetCount()](0,0,FALSE);

// version for gcc4.0 NEW :
     m_pSamples = new Table2D<ELEM_SAMPLE_TYPE>[m_pSampleFilesTab->GetCount()];
     for(int ii=0;ii<m_pSampleFilesTab->GetCount();ii++){
			m_pSamples[ii].InitConstructor( 0 , 0 , FALSE );
	  }

		CMyStrTable::iterator i;
		mystring szFileList;
		m_Count = 0;
		for(i=m_pSampleFilesTab->begin();i!=m_pSampleFilesTab->end();i++,m_Count++){
			mystring szSampleFName;
			szSampleFName << sample_dir << "/" << m_szMagnitude << "/" << i->c_str();
			szFileList << szSampleFName << ",";
			if(! CCDUtil::ReadFITSFile(m_pSamples[m_Count],szSampleFName.c_str() ) ){
				printf("Could not read sample file %s\n",szSampleFName.c_str());
				return 0;
			}
			m_pSamples[m_Count].CalcMaxAndPos();
			// printf("MAX : %d,%d,%d\n",m_pSamples[m_Count].m_MaxX,m_pSamples[m_Count].m_MaxY,m_pSamples[m_Count].m_MaxVal);
			total_size += (m_pSamples[m_Count].GetXSize()*m_pSamples[m_Count].GetYSize())*sizeof(ELEM_TYPE);
		}	
	}	
	mystring szSizeMsg;
	szSizeMsg << "Total size of samples of magnitude " << m_szMagnitude << " is " << (total_size/1000)
             << "." << (total_size % 1000) << "kB";
	MYTRACE3(gCCDTrace,szSizeMsg);
	printf("%s\n",szSizeMsg.c_str());
	return total_size;	
}

Table2D<ELEM_SAMPLE_TYPE>& CMagSamples::GetRandomSample( int* pIdx )
{
	Assert(m_Count>0,"Sample of magnitude %s not found",m_szMagnitude.c_str());
	LONG_T idx = CRandom::GetRandomInteger(0,m_Count);
	if( pIdx ){
		(*pIdx) = idx;
	}
	return m_pSamples[idx];
}

BOOL_T CMagSamples::CheckXY( int x, int y, int sample_x, int sample_y, 
									  int nFrame )
{
	int center_x = (gCCDParams.m_SizeX/2);
	int center_y = (gCCDParams.m_SizeY/2);

	double r = sqrt( (x-center_x)*(x-center_x) + (y-center_y)*(y-center_y) );
	double sample_r = sqrt( (sample_x-center_x)*(sample_x-center_x) + (sample_y-center_y)*(sample_y-center_y) );

	if ( fabs( r - sample_r ) < gCCDParams.m_fPutSampleCutOutInWidth || gCCDParams.m_fPutSampleCutOutInWidth<0 ){
		return TRUE;
	}
	return FALSE;
}

BOOL_T CMagSamples::GetRandomSample_ByName( int x, int y,
                                  Table2D<ELEM_SAMPLE_TYPE>& sample_in,
                                  Table2D<ELEM_SAMPLE_TYPE>& sample_out,
											 int nFrame )
{
	const char* szRA = sample_in.GetKeyTab().getKeyVal( PART_RA );
	const char* szDEC = sample_in.GetKeyTab().getKeyVal( PART_DEC );
	const char* szPartName = sample_in.GetKeyTab().getKeyVal( PARTNAME );
	const char* szPartMag = sample_in.GetKeyTab().getKeyVal( PARTMAG );
	const char* szPartNo = sample_in.GetKeyTab().getKeyVal( DIMAGE );
	
	if( szRA && szDEC && szPartName && szPartMag && 
		 szRA[0] && szDEC[0] && szPartName[0] && szPartMag[0] ){
		for(int i=0;i<m_Count;i++){
			const char* szCatRA = m_pSamples[i].GetKeyTab().getKeyVal( PART_RA );
			const char* szCatDEC = m_pSamples[i].GetKeyTab().getKeyVal( PART_DEC );
			const char* szCatName = m_pSamples[i].GetKeyTab().getKeyVal( PARTNAME );
			const char* szCatMag = m_pSamples[i].GetKeyTab().getKeyVal( PARTMAG );
			const char* szCatNo = sample_in.GetKeyTab().getKeyVal( DIMAGE );

			if( szCatRA && szCatDEC && szCatName && szCatMag && 
				 szCatRA[0] && szCatDEC[0] && szCatName[0] && szCatMag[0] ){
				if( strcmp( szPartMag , szCatMag )==0 ){
					if( strlen(szPartName)>10 && strlen(szCatName)>10 ){
						//if( strncmp(szPartName,szCatName,4)==0 && 
						//	 strcmp(szPartName+9,szCatName+9)==0 ){
						double ra_in_deg = atof(szRA)*15.00;
						double dec_in_deg = atof(szDEC);
						double cat_ra_in_deg = atof(szCatRA)*15.00;
						double cat_dec_in_deg = atof(szCatDEC);

						double diff_ra_in_arcsec = fabs(ra_in_deg-cat_ra_in_deg)*3600;
						double diff_dec_in_arcsec = fabs(dec_in_deg-cat_dec_in_deg)*3600;
						double coic_radius_in_sec = AstroAngle::rad2deg(gCCDParams.m_nCoicRadiusInRad)*3600;
	
						_TRACE_PRINTF_4("Sample with matching mag=%s found, diff_ra=%.2f arcsec, diff_dec=%.2f secsec\n",szCatMag,diff_ra_in_arcsec, diff_dec_in_arcsec);
						if( fabs(diff_ra_in_arcsec)<coic_radius_in_sec &&
                      fabs(diff_dec_in_arcsec)<coic_radius_in_sec ){
							const char* szPartX = m_pSamples[i].GetKeyTab().getKeyVal( PARTX );
							const char* szPartY = m_pSamples[i].GetKeyTab().getKeyVal( PARTY );
	
							if( szPartX && szPartY && szPartX[0] && szPartY[0] ){
								int part_x = atol( szPartX );
					         int part_y = atol( szPartY );
                                                                                
                                                                                
					         if( CheckXY( x , y, part_x, part_y ) ){
									BOOL_T bOK=TRUE;
			
									if( gCCDParams.m_nPutTakenFromFrameInRange>0 ){
										bOK = FALSE;
										if( szPartNo && szPartNo[0] && szCatNo && szCatNo[0] ){
											int part_no = atol( szPartNo );
											int cat_no = atol( szCatNo );
											if( part_no == cat_no ){
												bOK=TRUE;
											}
										}
									}
									if( bOK ){
										sample_out = m_pSamples[i];
										return TRUE;
									}
								}
							}
						}
					}
				}
			}
		}		
	}
	printf("could not find star matching mag=%s, ra=%s, dec=%s, name=%s, near (%d,%d)\n",szPartMag,szRA,szDEC,szPartName,x,y);
	return FALSE;
}

BOOL_T CMagSamples::GetRandomSample_CheckXY( int x, int y, 
													Table2D<ELEM_SAMPLE_TYPE>& sample_out,
													int nFrame )
{
	int nRejectedByFrame=0;
	int nClosest=1000000;
	int nClosestSample=-1;
	int nClosestFrame=-1;
	vector<int> matching_xy_table;

	Assert(m_Count>0,"Sample of magnitude %s not found",m_szMagnitude.c_str());
	for(int i=0;i<m_Count;i++){
		const char* szPartX = m_pSamples[i].GetKeyTab().getKeyVal( PARTX );
		const char* szPartY = m_pSamples[i].GetKeyTab().getKeyVal( PARTY );


		if( szPartX && szPartY && szPartX[0] && szPartY[0] ){
			int part_x = atol( szPartX );
			int part_y = atol( szPartY );			


			if( CheckXY( x , y, part_x, part_y ) ){
				BOOL_T bOK=TRUE;
				if( gCCDParams.m_nPutTakenFromFrameInRange>0 && nFrame>0 ){
					const char* szDIMAGE = m_pSamples[i].GetKeyTab().getKeyVal( DIMAGE );
					bOK=FALSE;
					if( szDIMAGE && szDIMAGE[0] ){
						int dimage_no = atol( szDIMAGE );
						if( abs( nFrame - dimage_no ) < gCCDParams.m_nPutTakenFromFrameInRange ){
							bOK = TRUE;							
						}
						if( abs( nFrame - dimage_no ) < nClosest ){
							nClosest = abs( nFrame - dimage_no );
							nClosestSample = i;
							nClosestFrame = dimage_no;
						}
					}
					if( !bOK ){
						nRejectedByFrame++;
					}
				}

				if( bOK ){
					matching_xy_table.push_back( i );
				}
			}
		}
	}

	if( matching_xy_table.size()==0 ){
		if( nClosestSample>0 && nRejectedByFrame>0 && nClosestFrame>0 ){
			// in case no sample matching frame is found closest one is used :
			for(int i=0;i<m_Count;i++){
				const char* szPartX = m_pSamples[i].GetKeyTab().getKeyVal( PARTX );
				const char* szPartY = m_pSamples[i].GetKeyTab().getKeyVal( PARTY );


				if( szPartX && szPartY && szPartX[0] && szPartY[0] ){
					int part_x = atol( szPartX );
					int part_y = atol( szPartY );			


					if( CheckXY( x , y, part_x, part_y ) ){
						BOOL_T bOK=TRUE;
						if( gCCDParams.m_nPutTakenFromFrameInRange>0 && nFrame>0 ){
							const char* szDIMAGE = m_pSamples[i].GetKeyTab().getKeyVal( DIMAGE );
							if( szDIMAGE && szDIMAGE[0] ){
		                  int dimage_no = atol( szDIMAGE );
								if( dimage_no == nClosestFrame ){
									matching_xy_table.push_back( i );
								}
							}
						}
					}
				}
			}
		}
	}


	if( matching_xy_table.size()>0 ){
		LONG_T idx = CRandom::GetRandomInteger(0,matching_xy_table.size());
		if( idx<matching_xy_table.size() ){
	   	sample_out = m_pSamples[ matching_xy_table[idx] ];	
			return TRUE;
		}else{
			return FALSE;
		}
	}
	return TRUE;
}

CccdSamples::CccdSamples() :m_bRead(FALSE),m_pMagSamples(NULL), m_pMagTab(NULL)
{
}

CccdSamples::~CccdSamples()
{
	if(m_pMagSamples)
		delete [] m_pMagSamples;
	if(m_pMagTab)
		delete [] m_pMagTab;
}

BOOL_T CccdSamples::ReadSamples( const char* samples_dir )
{
	if(m_bRead)
		return TRUE;
	LONG_T total_size=0;
	
	mystring szSamplesList,szSampleDir;
	if( samples_dir && samples_dir[0] ){
		szSampleDir = samples_dir;
	}else{
		szSamplesList = gCCDParams.m_szSampleDir;
	}
	
   szSamplesList << szSampleDir.c_str() << "/" << gCCDParams.m_szListName;	
	if( !MyFile::DoesFileExist( szSamplesList.c_str() ) ){
		printf("could not open list file : %s\n",szSamplesList.c_str());
		exit(0);
	}

	m_DescFile.Init( szSamplesList.c_str() );
	CSafeKeyTab& magList = m_DescFile.GetDescTab();

	if(m_pMagSamples)
		delete m_pMagSamples;
	m_pMagSamples = new CMagSamples[ magList.GetCount() ];
	m_pMagTab = new double[magList.GetCount()];

	for(int i=0;i<magList.GetCount();i++){		
		total_size += m_pMagSamples[i].ReadSamples( (m_DescFile.GetDescTab())[i], szSampleDir.c_str() );
		m_pMagTab[i] = atof( m_pMagSamples[i].m_szMagnitude.c_str() );
	}
	mystring szSizeMsg;
	szSizeMsg << "Total size of samples of all samples is " << (total_size/1000) 
             << "." << (total_size % 1000) << "kB";
	MYTRACE2(gCCDTrace,szSizeMsg);
	printf("%s\n",szSizeMsg.c_str());	
	if(total_size==0){
		mystring szErrMsg;
		szErrMsg << "No samples files available on path : " << gCCDParams.m_szSampleDir << ", exiting ...";
		MYTRACE2(gCCDTrace,szErrMsg);
		printf("%s\n",szErrMsg.c_str());		
		exit(-1);
	}

	my_sort_float( m_pMagTab, magList.GetCount() );

	printf("All available magnitudes :\n");
	for(int i=0;i<magList.GetCount();i++){
		printf("%i - %f\n",i,m_pMagTab[i]);
	}
	
	return TRUE;
}

CMagSamples* CccdSamples::GetMagSamples( int pos )
{
	CSafeKeyTab& magList = m_DescFile.GetDescTab();
	Assert( pos>=0 && pos<magList.GetCount(),"Sample range exceeded");
	return &(m_pMagSamples[pos]);
}

CMagSamples* CccdSamples::GetMagSamples( const char* szMag, int* pPos )
{
	CSafeKeyTab& magList = m_DescFile.GetDescTab();
	for(int i=0;i<magList.GetCount();i++){
      if(strcmp(m_pMagSamples[i].m_szMagnitude.c_str(),szMag)==0){
			if( pPos ){
				(*pPos) = i;
			}
			return &(m_pMagSamples[i]);
		}
   }
	return NULL;
}

Table2D<ELEM_SAMPLE_TYPE>& CccdSamples::GetRandomSample( const char* szMag, mystring* szSampleMag )
{
	CSafeKeyTab& magList = m_DescFile.GetDescTab();
	Assert(magList.GetCount()>0,"No samples available");
	if(!szMag || strlen(szMag)==0){
		LONG_T idx = CRandom::GetRandomInteger(0,magList.GetCount());
		if(szSampleMag)
			(*szSampleMag) = m_pMagSamples[idx].m_szMagnitude;
		return m_pMagSamples[idx].GetRandomSample();
	}else{
		int greaterCount = -1;
		double minMag = m_pMagTab[0];
		double mag = atof(szMag);
		Assert( mag>=minMag,"Cannot scale to requested magnitude - no samples lighter then %f available",minMag);
		for(int i=0;i<magList.GetCount();i++){
			if(m_pMagTab[i]>mag){
				break;
			}
			greaterCount = i;
		}
		Assert( greaterCount>=0, "No stars which can be scaled to %s",szMag);
		LONG_T idx = CRandom::GetRandomInteger(0,greaterCount);
		if(szSampleMag)
         (*szSampleMag) = m_pMagSamples[idx].m_szMagnitude;
		return m_pMagSamples[idx].GetRandomSample();
	}
}

