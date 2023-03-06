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
#include "ccd_util.h"
#include <myfile.h>
#include "ccd_matrix.h"
#include "ccd.h"
#include "ccd_globals.h"
#include <cfg.h>
#include <myutil.h>
#include <mystrtable.h>
#include <fits_file.h>
#include "ccd_trace.h"
#include <mykeytab.h>
#include <myprogress.h>
#include <cexcp.h>
#include <tab2Ddesc.h>
#include <mypixellist.h>
#include "ccd_photometry.h"

BOOL_T CCDUtil::bVerbose=FALSE;

CCDUtil::CCDUtil()
{}


CCDUtil::~CCDUtil()
{}



BOOL_T CCDUtil::BuildMedian( mystring szListFile, CCDMatrix& DarkFrame,
                			     mystring& szError, const char* szDark )
{
	if(!MyFile::DoesFileExist( szListFile.c_str() )){
      printf("list file %s, does not exist\n",szListFile.c_str());
     	return FALSE;
   }

	gCCDParams.InitParams();

   gCCDParams.SetParam("CCD_FULL_FRAMES_SAMPLE_DIR","");
   gCCDParams.SetParam( "CCD_FULL_FRAMES_LIST", szListFile.c_str() );
   gCCDParams.RefreshParams();

	CListFile frames_list(szListFile.c_str());
   CMyStrTable& list = frames_list.GetListTable();
	BOOL_T bRet = BuildMedian( list, DarkFrame, szError, szDark );
	return bRet;
}

BOOL_T CCDUtil::BuildMedian( CMyStrTable& list, CCDMatrix& DarkFrame,
                			     mystring& szError, const char* szDark )
{
   BOOL_T bOK = TRUE;

	CCDMatrix dark(0,0);
	if( szDark && strlen(szDark)){
		mystring szDarkFile = szDark;
		szDarkFile.env2str();
		if(!dark.ReadFITSFile( szDarkFile.c_str() )){
			printf("could not read frame : %s\n",szDarkFile.c_str());
			return FALSE;
		}				
	}

	CCDMatrix frame( 0, 0, FALSE );

	if(!frame.ReadFITSFile( list[0].c_str(), FALSE )){
      printf("could not read frame %s\n",list[0].c_str());
		return FALSE;
   }
	
	LONG_T xSize=frame.GetXSize();
	LONG_T ySize=frame.GetYSize();
	cCCD frames( xSize, ySize, list.size());
	

	int i=0;
	for(i=0;i<list.size() && bOK;i++){
      printf("reading frame# %d from  %s ...",i,list[i].c_str());fflush(0);
      if(frames[i].ReadFITSFile( list[i].c_str() )){
			printf("OK\n");fflush(0);
      }else{
         printf("count not read frame %s\n",list[i].c_str());fflush(0);
         bOK = FALSE;
      }
		if(szDark && strlen(szDark)){
			frames[i].Subtract( dark, frames[i], TRUE );
		}
   }

	printf("Frames initialized, calculating median frame ...\n");

/*	DarkFrame.Alloc(xSize,ySize);
	ELEM_TYPE** p_dark_data = DarkFrame.get_data_buffer_fast();
	LONG_T* all_values = new LONG_T[list.size()];
	LONG_T nCount = list.size();
	LONG_T median_pos=(nCount/2);
	
	CMyProgressBar bar( 0, ySize );	
	printf("median pos = %d\n",median_pos);
	for(register long y=0;y<ySize;y++){
      for(register long x=0;x<xSize;x++){
			for(i=0;i<nCount;i++){									
				all_values[i] = (frames[i].get_data_buffer_fast())[y][x];
			}
			my_qsort(all_values,nCount);
			p_dark_data[y][x] = all_values[ median_pos ];
		}
		bar.SetValue( y );
		bar.Update();
	}
	delete all_values;*/

	BOOL_T bRet = BuildMedian( frames.GetFrames(), frames.GetCount(), DarkFrame );

	DarkFrame.GetKeyTab() = frames[0].GetKeyTab();
	DarkFrame.GetKeyTab().Set( OBJECT, "DARKMEDIAN" );
	for(int i=0;i<list.GetCount();i++){
		mystring szFrame;
		szFrame << "FRAME" << i;
		DarkFrame.GetKeyTab().Add(szFrame.c_str(), list[i].c_str() );
	}			
	DarkFrame.GetKeyTab().Add( "MTIME" , (int)get_dttm() );
	mystring szDTM = get_date_time_string();
	DarkFrame.GetKeyTab().Add( "MDTM", szDTM.c_str() );

	return bRet;
}

BOOL_T CCDUtil::BuildMedian( CCDMatrix* pFrames , int nCount, CCDMatrix& DarkFrame )
{
	int xSize = pFrames[0].GetXSize();
	int ySize = pFrames[0].GetYSize();
	DarkFrame.Alloc(xSize,ySize);
	ELEM_TYPE** p_dark_data = DarkFrame.get_data_buffer_fast();
	LONG_T* all_values = new LONG_T[nCount];
	LONG_T median_pos=(nCount/2);

	CMyProgressBar bar( 0, ySize );	
	printf("median pos = %d\n",median_pos);
	for(register long y=0;y<ySize;y++){
      for(register long x=0;x<xSize;x++){
			for(register int i=0;i<nCount;i++){									
				all_values[i] = (pFrames[i].get_data_buffer_fast())[y][x];
			}
			my_qsort(all_values,nCount);
			p_dark_data[y][x] = all_values[ median_pos ];
		}
		bar.SetValue( y );
		bar.Update();
	}
	delete [] all_values;
	return TRUE;	
}

/*BOOL_T CCDUtil::BuildMedian( CCDMatrix* pFrames , int nCount, CCDMatrix& DarkFrame )
{
	int xSize = pFrames[0].GetXSize();
	int ySize = pFrames[0].GetYSize();
	DarkFrame.Alloc(xSize,ySize);
	ELEM_TYPE** p_dark_data = DarkFrame.get_data_buffer_fast();
	LONG_T* all_values = new LONG_T[nCount];
	LONG_T median_pos=(nCount/2);

	int max_elem_value=66000;
	int* p_values = new int[ max_elem_value ];
	CMyProgressBar bar( 0, ySize );	
	printf("median pos = %d\n",median_pos);
	for(register long y=0;y<ySize;y++){
      for(register long x=0;x<xSize;x++){
			for(register int i=0;i<max_elem_value;i++){
				p_values[i] = 0;
			}
			for(register int i=0;i<nCount;i++){									
				// all_values[i] = (pFrames[i].get_data_buffer_fast())[y][x];
				int val = (pFrames[i].get_data_buffer_fast())[y][x];	
				p_values[ val ]++;
			}
			// my_qsort(all_values,nCount);
			int smallerCount=0;
			for(register int i=0;i<max_elem_value;i++){
				smallerCount += p_values[i];
				if( smallerCount>=(nCount/2)){
					p_dark_data[y][x] = i;
					break;
				}
			}
			// p_dark_data[y][x] = all_values[ median_pos ];
		}
		bar.SetValue( y );
		bar.Update();
	}
	delete [] all_values;
	return TRUE;	
}*/


BOOL_T CCDUtil::BuildMedian( Table2D<double>* pCCDTab, int nCount,
									  Table2D<double>& DarkFrame )
{
	int xSize = pCCDTab[0].GetXSize();
	int ySize = pCCDTab[0].GetYSize();	

	DarkFrame.Alloc(xSize,ySize);
	double** p_dark_data = DarkFrame.get_data_buffer_fast();
	double* all_values = new double[ nCount ];
	LONG_T median_pos=(nCount/2);
	
	
	printf("median pos = %d\n",median_pos);
	for(register long y=0;y<ySize;y++){
      for(register long x=0;x<xSize;x++){
			for(register int i=0;i<nCount;i++){									
				all_values[i] = (pCCDTab[i].get_data_buffer_fast())[y][x];
			}
			my_sort_float(all_values,nCount);
			p_dark_data[y][x] = all_values[ median_pos ];
		}
	}
	delete all_values;
	return TRUE;

}

BOOL_T CCDUtil::BuildMedian( CMyStrTable& list, Table2D<double>& DarkFrame,
									  int xSize, int ySize,
                			     mystring& szError )
{
   BOOL_T bOK = TRUE;

	
	
//	Table2D<double>* pCCDTab = new Table2D<double>[list.size()]( xSize, ySize );
// version for gcc4.0 - NEW :
	Table2D<double>* pCCDTab = new Table2D<double>[list.size()];
	

	int i=0;
	for(i=0;i<list.size() && bOK;i++){
// version for gcc4.0 - NEW :
      pCCDTab[i].InitConstructor( xSize, ySize );

      printf("reading frame# %d from  %s ...",i,list[i].c_str());fflush(0);
		CFITSFile<double> in;
		mystring szError;
		if( in.ReadFITSFile( pCCDTab[i], szError, list[i].c_str() ) ){
			printf("OK\n");fflush(0);
      }else{
         printf("count not read frame %s\n",list[i].c_str());fflush(0);
			printf("error : %s\n",szError.c_str());
         bOK = FALSE;
      }
   }

	printf("Frames initialized, calculating median frame ...\n");
	
	BOOL_T bRet = BuildMedian( pCCDTab, list.size(), DarkFrame );
	return bRet;


/*	DarkFrame.Alloc(xSize,ySize);
	double** p_dark_data = DarkFrame.get_data_buffer_fast();
	double* all_values = new double[list.size()];
	LONG_T nCount = list.size();
	LONG_T median_pos=(nCount/2);
	
	
	printf("median pos = %d\n",median_pos);
	for(register long y=0;y<ySize;y++){
      for(register long x=0;x<xSize;x++){
			for(i=0;i<nCount;i++){									
				all_values[i] = (pCCDTab[i].get_data_buffer_fast())[y][x];
			}
			my_sort_float(all_values,nCount);
			p_dark_data[y][x] = all_values[ median_pos ];
		}
	}
	delete all_values;
	return TRUE;*/
}


void CCDUtil::CopyFromBigElemToElem(  Table2D<BIG_ELEM_TYPE>& bigElem, Table2D<ELEM_TYPE>& Elem )
{
	Elem.Alloc( bigElem.GetXSize(), bigElem.GetYSize() );
	ELEM_TYPE* p_elem = Elem.get_data_buffer();
	BIG_ELEM_TYPE* p_big_elem = bigElem.get_data_buffer();
	long size = bigElem.GetXSize()*bigElem.GetYSize();

	for(register int i=0;i<size;i++){
		if( p_big_elem[i]>=0 && p_big_elem[i]<=30000){
			p_elem[i] = (ELEM_TYPE)p_big_elem[i];
		}else{
			if( p_big_elem[i]<0 ){
				p_elem[i] = 0;
			}
			if( p_big_elem[i]>30000 ){
				p_elem[i] = 30000;
			}
		}
	}
	Elem.m_X_On_Big = bigElem.m_X_On_Big;
	Elem.m_Y_On_Big = bigElem.m_Y_On_Big;
}

void CCDUtil::CopyFromElemToBigElem(  const Table2D<ELEM_TYPE>& Elem, Table2D<BIG_ELEM_TYPE>& bigElem )
{
	bigElem.Alloc( Elem.GetXSize(), Elem.GetYSize() );
	ELEM_TYPE* p_elem = Elem.get_data_buffer();
	BIG_ELEM_TYPE* p_big_elem = bigElem.get_data_buffer();
	long size = bigElem.GetXSize()*bigElem.GetYSize();

	for(register int i=0;i<size;i++){
		p_big_elem[i] = (BIG_ELEM_TYPE)p_elem[i];
	}
	bigElem.m_X_On_Big = Elem.m_X_On_Big;
	bigElem.m_Y_On_Big = Elem.m_Y_On_Big;
}


BOOL_T CCDUtil::ReadFITSFile( Table2D<short>& sample, const char* fname)
{
	mystring szError,szFName = fname;
	szFName.env2str();

	
	CFITSFile<short> input;
	if(!input.ReadFITSFile( sample, szError, szFName.c_str() )){
		printf("could not read file : %s\n",szFName.c_str());
		return FALSE;
	}else{
		if(!szError.empty()){
			MYTRACE1(gCCDTrace,"Warning while reading FITS file : " << szFName
									  << ", error is : " << szError);
		}
	}
	return TRUE;
}

BOOL_T CCDUtil::ReadFITSFile( Table2D<float>& sample, const char* fname)
{
	mystring szError,szFName = fname;
	szFName.env2str();

	
	CFITSFile<float> input;
	if(!input.ReadFITSFile( sample, szError, szFName.c_str() )){
		MYTRACE1(gCCDTrace,"Error while reading FITS file : " << szFName
                         << ", error is : " << szError);
		return FALSE;
	}else{
		if(!szError.empty()){
			MYTRACE1(gCCDTrace,"Warning while reading FITS file : " << szFName
									  << ", error is : " << szError);
		}
	}
	return TRUE;
}



BOOL_T CCDUtil::ReadFITSFile( Table2D<BIG_ELEM_TYPE>& sample, const char* fname)
{
	mystring szError,szFName = fname;
	szFName.env2str();

	
	CFITSFile<BIG_ELEM_TYPE> input;
	if(!input.ReadFITSFile( sample, szError, szFName.c_str() )){
		MYTRACE1(gCCDTrace,"Error while reading FITS file : " << szFName
                         << ", error is : " << szError);
		return FALSE;
	}else{
		if(!szError.empty()){
			MYTRACE1(gCCDTrace,"Warning while reading FITS file : " << szFName
									  << ", error is : " << szError);
		}
	}
	return TRUE;

}

BOOL_T CCDUtil::Divide( CCDMatrix& frame,  Table2D<double>& flat, int border )
{
	ELEM_TYPE* pdata = frame.get_data_buffer();
	double* pflat = flat.get_data_buffer();
	
	int size1 = frame.GetXSize()*frame.GetYSize();
	int size2 = flat.GetXSize()*flat.GetYSize();
	Assert(size1==size2,"Sizes must be equal to divide frames");

	for(register int i=0;i<size1;i++){
		if( pflat[i]!=0 ){
			int tmp_data = (pdata[i]/pflat[i]);
			if( tmp_data > USHRT_MAX )
				tmp_data = USHRT_MAX;
			pdata[i] = (ELEM_TYPE)tmp_data;
		}
	}

	if(border>0){
		register int xSize = frame.GetXSize();
	   register int ySize = frame.GetYSize();
		ELEM_TYPE** pfast = frame.get_data_buffer_fast();
		for(register int y=0;y<ySize;y++){
			for(register int x=0;x<xSize;x++){
				if(x<border || x>(xSize-border) || y<border || y>(ySize-border)){
					 pfast[y][x] = 0;
				}
			}
		}
	}

	return (size1==size2);
}

BOOL_T CCDUtil::Divide( CCDMatrix& frame,  Table2D<float>& flat, int border )
{
	ELEM_TYPE* pdata = frame.get_data_buffer();
	float* pflat = flat.get_data_buffer();
	
	int size1 = frame.GetXSize()*frame.GetYSize();
	int size2 = flat.GetXSize()*flat.GetYSize();
	Assert(size1==size2,"Sizes must be equal to divide frames");


	for(register int i=0;i<size1;i++){
		if( pflat[i]!=0 ){
			int tmp_data = (pdata[i]/pflat[i]);
			if( tmp_data > USHRT_MAX )
				tmp_data = USHRT_MAX;
			pdata[i] = (ELEM_TYPE)tmp_data;
		}
	}

	if(border>0){
		register int xSize = frame.GetXSize();
	   register int ySize = frame.GetYSize();
		ELEM_TYPE** pfast = frame.get_data_buffer_fast();
		for(register int y=0;y<ySize;y++){
			for(register int x=0;x<xSize;x++){
				if(x<border || x>(xSize-border) || y<border || y>(ySize-border)){
					 pfast[y][x] = 0;
				}
			}
		}
	}

	return (size1==size2);
}

BOOL_T CCDUtil::WriteToFITSFile( Table2D<float>& float_image, const char* fname, CKeyTab* pHduTab)
{
	mystring szError,szFName = fname;
	szFName.env2str();	

	CFITSFile<float> out;
	
	if(!strstr(szFName.c_str(),".fit"))
		szFName << ".fit";

	CSafeKeyTab& keys = float_image.GetKeyTab();
	for( int i=0;i<keys.GetCount();i++ ){
		 out.AddHdu( keys[i].szName.c_str(), keys[i].szValue.c_str() );
	}

	if(pHduTab){
		CKeyTab::iterator i;
		for(i=pHduTab->begin();i!=pHduTab->end();i++){
			out.AddHdu(i->szName.c_str(),i->szValue.c_str());
		}
	}

	if (!out.WriteToFITSFile( float_image.get_data_buffer(),
									  float_image.GetXSize(), float_image.GetYSize(), 
                             szError, szFName.c_str() ) ){
		MYTRACE1(gCCDTrace,"Error while writing to FITS file : " << szFName
                         << ", error is : " << szError);
		return FALSE;
	}else{
		if(!szError.empty()){
			MYTRACE1(gCCDTrace,"Warning while writing to FITS file : " << szFName
									  << ", error is : " << szError);
		}

	}
	
	return TRUE;
}


BOOL_T CCDUtil::WriteToFITSFile( Table2D<short>& sample, const char* fname, CKeyTab* pHduTab)
{
	mystring szError,szFName = fname;
	szFName.env2str();	

	CFITSFile<ELEM_TYPE> out;
	
	if(!strstr(szFName.c_str(),".fit"))
		szFName << ".fit";
	if(pHduTab){
		CKeyTab::iterator i;
		for(i=pHduTab->begin();i!=pHduTab->end();i++){
			out.AddHdu(i->szName.c_str(),i->szValue.c_str());
		}
/*	}else{
		for(int i=0;i<m_KeyTab.GetCount();i++){
			if( strcmp( m_KeyTab[i].szName.c_str(), fitsHeaderKeyDefTab[eUT_START].szKeyName)==0 ||
				 strcmp( m_KeyTab[i].szName.c_str(), fitsHeaderKeyDefTab[eDATE_OBS].szKeyName)==0 )
				out.AddHdu( m_KeyTab[i].szName.c_str(), m_KeyTab[i].szValue.c_str() );
		}*/
	}

	if (!out.WriteToFITSFile( sample.get_data_buffer(),
									  sample.GetXSize(), sample.GetYSize(), 
                             szError, szFName.c_str() ) ){
		MYTRACE1(gCCDTrace,"Error while writing to FITS file : " << szFName
                         << ", error is : " << szError);
		return FALSE;
	}else{
		if(!szError.empty()){
			MYTRACE1(gCCDTrace,"Warning while writing to FITS file : " << szFName
									  << ", error is : " << szError);
		}

	}
	
	return TRUE;
}

BOOL_T CCDUtil::WriteToFITSFile( Table2D<BIG_ELEM_TYPE>& sample, const char* fname, CKeyTab* pHduTab)
{
	mystring szError,szFName = fname;
	szFName.env2str();	

	CFITSFile<BIG_ELEM_TYPE> out;
	
	if(!strstr(szFName.c_str(),".fit"))
		szFName << ".fit";
	if(pHduTab){
		CKeyTab::iterator i;
		for(i=pHduTab->begin();i!=pHduTab->end();i++){
			out.AddHdu(i->szName.c_str(),i->szValue.c_str());
		}
/*	}else{
		for(int i=0;i<m_KeyTab.GetCount();i++){
			if( strcmp( m_KeyTab[i].szName.c_str(), fitsHeaderKeyDefTab[eUT_START].szKeyName)==0 ||
				 strcmp( m_KeyTab[i].szName.c_str(), fitsHeaderKeyDefTab[eDATE_OBS].szKeyName)==0 )
				out.AddHdu( m_KeyTab[i].szName.c_str(), m_KeyTab[i].szValue.c_str() );
		}*/
	}

	if (!out.WriteToFITSFile( sample.get_data_buffer(),
									  sample.GetXSize(), sample.GetYSize(), 
                             szError, szFName.c_str() ) ){
		MYTRACE1(gCCDTrace,"Error while writing to FITS file : " << szFName
                         << ", error is : " << szError);
		return FALSE;
	}else{
		if(!szError.empty()){
			MYTRACE1(gCCDTrace,"Warning while writing to FITS file : " << szFName
									  << ", error is : " << szError);
		}

	}
	
	return TRUE;
}

/*BOOL_T CCDUtil::WriteToFITSFileInt( Table2D<int>& sample, const char* fname, CKeyTab* pHduTab)
{
	mystring szError,szFName = fname;
	szFName.env2str();	

	CFITSFile<int> out;
	
	if(!strstr(szFName.c_str(),".fit"))
		szFName << ".fit";
	if(pHduTab){
		CKeyTab::iterator i;
		for(i=pHduTab->begin();i!=pHduTab->end();i++){
			out.AddHdu(i->szName.c_str(),i->szValue.c_str());
		}
	}

	if (!out.WriteToFITSFile( sample.get_data_buffer(),
									  sample.GetXSize(), sample.GetYSize(), 
                             szError, szFName.c_str() ) ){
		MYTRACE1(gCCDTrace,"Error while writing to FITS file : " << szFName
                         << ", error is : " << szError);
		return FALSE;
	}else{
		if(!szError.empty()){
			MYTRACE1(gCCDTrace,"Warning while writing to FITS file : " << szFName
									  << ", error is : " << szError);
		}

	}
	
	return TRUE;
}*/

BOOL_T CCDUtil::SumFrames( CMyStrTable& list, CCDMatrix& out,
                          Table2D<BIG_ELEM_TYPE>& tmp, CCDMatrix* pDark )
{
	if( list.size()<=0 )
		return FALSE;
	CCDMatrix frame_tmp( 0, 0, FALSE  );

	// read frame to check size :
	if( !frame_tmp.ReadFITSFile( list[0].c_str() , FALSE ) ){
		printf("could not read init frame : %s\n",list[0].c_str() );
		return FALSE;
	}

	if( tmp.GetXSize()!=frame_tmp.GetXSize() ||
       tmp.GetYSize()!=frame_tmp.GetYSize() ){
		tmp.Alloc( frame_tmp.GetXSize(), frame_tmp.GetYSize() );		
	}
	int xFirstFrame = frame_tmp.GetXSize();
	int yFirstFrame = frame_tmp.GetYSize();

	BIG_ELEM_TYPE* tmp_data = tmp.get_data_buffer();	
	int size = tmp.GetXSize()*tmp.GetYSize();

	for(int i=0;i<list.size();i++){
		CCDMatrix frame( 0, 0, FALSE  );

		if( !frame.ReadFITSFile( list[i].c_str(), FALSE ) ){
			printf("could not read frame : %s\n",list[i].c_str() );
			return FALSE;
		}
		ELEM_TYPE* in_data = frame.get_data_buffer();

		if( frame.GetXSize()!=xFirstFrame || frame.GetYSize()!=yFirstFrame ){
			printf("cannot add frames of different sizes\n");
			return FALSE;
		}

		if( pDark ){
			if( frame.GetXSize()==pDark->GetXSize() && frame.GetYSize()==pDark->GetYSize() ){
				frame.Subtract( *pDark, frame, TRUE );
			}
		}		

		for(register int j=0;j<size;j++){
			tmp_data[j] += in_data[j];
		}

		if( i==0 ){ // NEW 20041025
			// copy FITS header from first averaged frame :
			out.GetKeyTab() = frame.GetKeyTab();
			const char* szDIMAGE = frame.getKeyValue( DIMAGE );
			if( szDIMAGE && szDIMAGE[0] ){
				out.GetKeyTab().Set( AVSTART, szDIMAGE );
			}
		}else{
			if( i==(list.size()-1) ){
				const char* szDIMAGE = frame.getKeyValue( DIMAGE );
				if( szDIMAGE && szDIMAGE[0] ){
					out.GetKeyTab().Set( AVEND, szDIMAGE );
				}
			}
		}		
	}	

	if( out.GetXSize()!=tmp.GetXSize() ||
       out.GetYSize()!=tmp.GetYSize() ){
      out.Alloc( xFirstFrame, yFirstFrame );
   }
	
	int nFrames = list.size();
	ELEM_TYPE* out_data = out.get_data_buffer();
	for(int j=0;j<size;j++){
		out_data[j] = ( tmp_data[j] / nFrames );
	}
	return TRUE;	
}

// static CMyMutex gSaveSumMutex;
BOOL_T CCDUtil::SaveSumFrame( CMyStrTable& list, const char* szOutName,
										CCDMatrix* pDark )
{
	// lock to protect optimization statics :
	// CMutexUnLock unlock( &gSaveSumMutex, TRUE );
	// gSaveSumMutex.Lock();

	CCDMatrix sumFrame( 0, 0, FALSE );
	Table2D<BIG_ELEM_TYPE> tmp(0,0);
	sumFrame.SetData( 0 );
	tmp.SetData( 0 );	

	BOOL_T bOK = SumFrames( list, sumFrame, tmp, pDark );
	if( bOK ){
		sumFrame.WriteToFITSFile( szOutName );
	}
	return bOK;
}


BOOL_T CCDUtil::SumFrames( Table2D<BIG_ELEM_TYPE>& sumFrame, CCDMatrix& newFrame )
{
	BIG_ELEM_TYPE* sum_data = sumFrame.get_data_buffer();	
	ELEM_TYPE* new_data = newFrame.get_data_buffer();

	if( newFrame.GetXSize()!=sumFrame.GetXSize() ||
       newFrame.GetYSize()!=sumFrame.GetYSize() ){
		return FALSE;
	}

	int size = newFrame.GetXSize()*newFrame.GetYSize();

	for(register int j=0;j<size;j++){
		/*int x = (j % newFrame.GetXSize());
      int y = (j / newFrame.GetXSize());
      if( x==1200 && y==972 )
         printf("odo");*/

		sum_data[j] += new_data[j];
	}	
	return TRUE;	
}

void CCDUtil::GetPartFromBigElemToElem( Table2D<BIG_ELEM_TYPE>& bigElem,
													  Table2D<ELEM_TYPE>& Elem,
													  int low_x, int low_y, int up_x, int up_y )
{
	int xSize = bigElem.GetXSize();
	int ySize = bigElem.GetYSize();
	if(low_x<0)
		low_x = 0;
	if(low_y<0)
		low_y = 0;
	if( up_x >= xSize )
		up_x = ( xSize-1 );
	if( up_x >= xSize )
		up_x = ( xSize-1 );
	
	
	int len_x = (up_x-low_x+1);
	int len_y = (up_y-low_y+1);

	// fist re-alloc :
	Elem.Alloc( len_x, len_y );

	BIG_ELEM_TYPE** in_data = bigElem.get_data_buffer_fast();
   ELEM_TYPE** out_data = Elem.get_data_buffer_fast();

	for(register int y=low_y;y<up_y;y++){
		for(register int x=low_x;x<up_x;x++){
			int x_small = ( x-low_x );
			int y_small = ( y-low_y );
			int val = in_data[y][x];	

			out_data[y_small][x_small] = val;
		}
	}
}

mystring CCDUtil::GetKeyValue( const char* fits_file, const char* key_name )
{
	mystring szRet;
	CSafeKeyTab keyTab;
	CFITSFile<ELEM_TYPE> in;
   if(!in.ReadFITSHeader( keyTab, fits_file )){
      printf("could not read FITS file : %s\n",fits_file );
		return szRet;
   }
	for(int i=0;i<keyTab.GetCount();i++){
		if( strcmp( keyTab[i].szName.c_str(), key_name )==0 ){
			szRet = keyTab[i].szValue;			
			break;
		}
	}
	return szRet;
}

BOOL_T CCDUtil::CalibrateFrame( CCDMatrix& image, Table2D<float>& calibrated )
{
	Table2D<float> a( image.GetXSize(), image.GetYSize() );
	Table2D<float> b( image.GetXSize(), image.GetYSize() );

	if( !ReadFITSFile( a, "a.fit" ) ){
		printf("could not read fits file a.fit\n");
		return FALSE;
	}	
	if( !ReadFITSFile( b, "b.fit" ) ){
		printf("could not read fits file b.fit\n");
		return FALSE;
	}	

	ELEM_TYPE* p_data = image.get_data_buffer();
	float* a_data = a.get_data_buffer();
	float* b_data = b.get_data_buffer();
	float* out_data = calibrated.get_data_buffer();

	register int size = image.GetXSize()*image.GetYSize();
	int xSize = image.GetXSize();
	int ySize = image.GetYSize();
	for( register int i=0;i<size;i++){
		int x = (i % xSize);
		int y = (i / xSize);


		// float tmp = ( ((float)p_data[i]) - b_data[i] )/( a_data[i] );
		float tmp = ( ((float)p_data[i])*( 1.00 - b_data[i] ) )/( a_data[i] );
		// float tmp = p_data[i];

		// OK :
		out_data[i] = (  ((float)p_data[i])/( a_data[i] )  )*10000;
		if( out_data[i]>60000 || out_data[i]<0 )
			out_data[i] = 0;

		if( ( x==1285 && y==1041 ) || ( x==1258 && y==1019 ) ){
         printf("odo value(%d,%d)=%.2f\n",x,y,out_data[i]);
      }


//		out_data[i] = ((float)(p_data[i]))/( a_data[i] );
		// out_data[i] =  0;
	}

	return TRUE;
}

BOOL_T CCDUtil::Normalize( CCDMatrix& image, Table2D<float>& normalized )
{
	double mean,rms;
	image.GetMeanAndRMS( mean, rms );
	ELEM_TYPE* data_in = image.get_data_buffer();
	float* data_out = normalized.get_data_buffer();

	if( image.GetXSize()!=normalized.GetXSize() || image.GetYSize()!=normalized.GetYSize() ){
		return FALSE;
	}

	register int size = image.GetXSize()*image.GetYSize();
	for(register int i=0;i<size;i++){
		data_out[i] = ((double)data_in[i])/mean;
	}

	return TRUE;
}

BOOL_T CCDUtil::AddFrame( Table2D<BIG_ELEM_TYPE>& out_frame, CCDMatrix& in_frame )
{
   BIG_ELEM_TYPE* p_data_int = out_frame.get_data_buffer();
   ELEM_TYPE* p_data = in_frame.get_data_buffer();

   int size1 = out_frame.GetXSize()*out_frame.GetYSize();
	int size2 = in_frame.GetXSize()*in_frame.GetYSize();
	int size = size1;
	if( size2 < size1 ){
		size = size2;
	}
	
	if( size1 != size2 ){
		printf("ERROR : CCDUtil::AddFrame adding frames of different sizes !!!\n");
		return FALSE;
	}

   for(register int i=0;i<size;i++){
      p_data_int[i] += p_data[i];
   }
	return TRUE;
}

void CCDUtil::NormalizeFrame( Table2D<BIG_ELEM_TYPE>& in_frame, CCDMatrix& out_frame, int cnt )
{
	if( CCDUtil::bVerbose ){
		printf_now2("Normalizing by factor %d ...",cnt);
	}

	BIG_ELEM_TYPE* p_data_int = in_frame.get_data_buffer();
	ELEM_TYPE* p_data_out = out_frame.get_data_buffer();
	int size = in_frame.GetXSize()*in_frame.GetYSize();

	for(register int i=0;i<size;i++){
		double tmp = ( double(p_data_int[i]) / double(cnt) );
		p_data_out[i] = my_round(tmp);
	}

	if( CCDUtil::bVerbose ){
		my_printf_now("OK\n");
	}
}


// 20061202 - changed frame with header - now all information is used from 
// first image and only UT-END and DATE-END are taken from last image 
// also not whole image is read but only header !
void CCDUtil::CalcAndSaveAver( CMyStrTable& file_tab, CMyStrTable& aver_file_list,
							 int first, int i,
							 Table2D<BIG_ELEM_TYPE>* aver_frame,
							 CCDMatrix& sav, float* flat_data, ELEM_TYPE* dark_data, 
							 time_t minTime, time_t maxTime,
							 mystring& szProcessedList, mystring& szOutDir,
							 int& gOutNo,BOOL_T bDARK,BOOL_T bFLAT,
							 BOOL_T bReduced, int x_size, int y_size,
							 BOOL_T gCutBorder, int* gBorderValues,
							 BOOL_T bFH, BOOL_T bFV,mystring& szOutList,
							 const char* szFramesDir, BOOL_T bSaveFlip,
							 CSafeKeyTab* pAddKeys )
{
				int last = ( first + aver_file_list.size() - 1 );
				if( last < 0 )	
					printf("ERROR IN CODE : CalcAndSaveAver !!!\n");
				
				printf("calculating average of frames from %s to %s...\n",
							file_tab[first].c_str(),file_tab[last].c_str());

/*				int middle = (i-1+first)/2;

				mystring szMiddleFrame=file_tab[middle].c_str();
				if( szFramesDir && szFramesDir[0] ){
					szMiddleFrame = szFramesDir;
					szMiddleFrame << "/" << file_tab[middle].c_str();
				}
				printf("reading middle frame (pos=%d) : %s\n",middle,szMiddleFrame.c_str() );
				if(!sav.ReadFITSFile( szMiddleFrame.c_str() )){
      		   printf("could not read frame : %s\n",szMiddleFrame.c_str());
		      	exit(0);
	      	}*/
				printf("reading first frame : %s\n",file_tab[0].c_str() );
//				CSafeKeyTab frame0_header;
//				CFITSFile<ELEM_TYPE> in;
//				if(!in.ReadFITSHeader( frame0_header, file_tab[first].c_str() )){

				CFITSFile<ELEM_TYPE> in;				
				// read first image to use its header :
				if(!sav.ReadFITSFile( file_tab[first].c_str() )){
// reading of header only was bad !
// wrong image size was allocated and program crashed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//				sav.GetKeyTab().Clear();
//				if(!in.ReadFITSHeader( sav.GetKeyTab(), file_tab[first].c_str() )){
					printf("could not read FITS file : %s\n",file_tab[first].c_str() );
					return;
				}

				CSafeKeyTab last_header;
				if(!in.ReadFITSHeader( last_header, file_tab[last].c_str() )){
					printf("could not read FITS file : %s\n",file_tab[last].c_str() );
					return;
				}

				for(int f=0;f<aver_file_list.size();f++){
					mystring szKey;
					szKey << "AVERF" << f;



					mystring szFName;
					szFName = getfname_with_ext( aver_file_list[f] );
					printf("aver file %d = %s / %s\n",f,aver_file_list[f].c_str(),szFName.c_str() );
					sav.GetKeyTab().Add( szKey.c_str(), szFName.c_str() );


					if( strlen( szProcessedList.c_str() ) ){
						MyOFile proc_list( szProcessedList.c_str() , "a+" );
						proc_list.Printf("%s\n",aver_file_list[f].c_str());
					}
				}
				// adding keywords :
				sav.GetKeyTab().Set( "DARK_SUB", "DARK" );
				sav.GetKeyTab().Set( "FLAT_DIV", "FLAT" );
				
				// min and max unix time :
				sav.GetKeyTab().Set( "MINTIME", (int)minTime );
				sav.GetKeyTab().Set( "MAXTIME", (int)maxTime );


				// add UT_END from last fram  :
				const char* szUTEND = last_header.getKeyVal( UT_END );
				if( szUTEND && szUTEND[0] ){
					sav.GetKeyTab().Set( UT_END , szUTEND );
				}
				const char* szDATEEND = last_header.getKeyVal( DATE_END );
				if( szDATEEND && szDATEEND[0] ){
					sav.GetKeyTab().Set( DATE_END, szDATEEND );
				}

//				double half_time = (double(minTime)+double(maxTime))/2.0;
//				time_t ut_time = (time_t)half_time;
				// to be consistant with driver convention that TIME_UT is 
				// start time of image taking :
				time_t ut_time = minTime;
				sav.GetKeyTab().Set( TIME_UT, (int)ut_time );

				sav.GetKeyTab().Set( DIMAGE, gOutNo );

				if( pAddKeys ){
					for(int n=0;n<pAddKeys->GetCount();n++){
						sav.GetKeyTab().Add( (*pAddKeys)[n].szName.c_str(), (*pAddKeys)[n].szValue.c_str() );
					}
				}

				mystring szBaseName = getfname( file_tab[first].c_str() );
				mystring szOutName;

				char szTmp[256];

				if( strlen( szOutDir.c_str() ) ){
   	       	szOutName << szOutDir.c_str() << "/";
      	   }
				mystring szOutFileName;
			
				sprintf(szTmp,"%.3d",gOutNo);
				gOutNo++;
				
				char szTmp2[128];
				strcpy(szTmp2,szBaseName.c_str());
				int len=strlen(szTmp2);
				szTmp2[len-5]='\0';
				szOutFileName << szTmp2 << szTmp << ".fit";		
				if( gCCDParams.m_eCompressFITS == eFITSComprASAS ){
					szOutFileName << "c";
				}

				szOutName << szOutFileName;


				if( bDARK || bFLAT ){
					int aver_count = aver_file_list.size();
					printf("subtracting dark and dividing by flat ...\n");
					if( bReduced ){
						printf("ERROR in code/call !!!!!\n");
						printf("Reduction already done - cannot repeat !\n");
					}
					BIG_ELEM_TYPE* data_aver = aver_frame->get_data_buffer();
					ELEM_TYPE* data = sav.get_data_buffer();
					int size=x_size*y_size;

					for(register int j=0;j<size;j++){
						double new_value = data_aver[j];
						new_value = new_value / aver_count;

						if( bDARK ){
							new_value = new_value - dark_data[j];
							if( new_value < 0 )
								new_value = 0;
						}
						if( bFLAT && flat_data[j]>0 ){
							new_value = new_value / double(flat_data[j]);
							if( new_value >= USHRT_MAX )
								new_value=USHRT_MAX;
						}
						int red_value = my_round( new_value );
						data[j] = (ELEM_TYPE)red_value;
					}										
				}else{
					CCDUtil::NormalizeFrame( *aver_frame, sav, aver_file_list.size() );
				}

				// cut border before flip , we know original borders !
				if( gCutBorder ){
					printf("removing border %d,%d,%d,%d ...\n",gBorderValues[0],
						gBorderValues[1],gBorderValues[2],gBorderValues[3]);

					CCDMatrix tmp;
					tmp.GetKeyTab() = sav.GetKeyTab();
					sav.CutBorder( gBorderValues[0], gBorderValues[1],gBorderValues[2],gBorderValues[3], tmp );
					sav = tmp;									
				}


				if( bFH ){
					printf("horizontal flip of averaged frame\n");fflush(0);
					sav.FH();
					sav.GetKeyTab().Set( "FLIP", "FH" );
				}
				if( bFV ){
					printf("vertical flip of averaged frame\n");fflush(0);
         		sav.FV();
					sav.GetKeyTab().Set( "FLIP", "FV" );
	   	   }
				if( (bFH || bFV) && bSaveFlip ){
					sav.GetKeyTab().Set( "FLIP", (int)eReverseImageNone );
				}

				double avg,rms;
				sav.GetStatistics( 30, 30, (sav.GetXSize()-30) ,(sav.GetYSize()-30), avg, rms );
				sav.setKeyValue( AVERAGEKEY, avg );
				sav.setKeyValue( RMSKEY , rms );

				if( strlen( szOutList.c_str() ) ){
					printf("writing file %s\n",szOutName.c_str() );
   	   	   sav.WriteToFITSFile( szOutName.c_str() );

					mystring szListFile;
         	   if( strlen( szOutDir.c_str() ) ){
            	   szListFile << szOutDir.c_str() << "/";
	            }
   	         szListFile << szOutList.c_str();
      	      MyOFile out( szListFile.c_str(), "a+" );
         	   out.Printf("%s\n",szOutFileName.c_str() );
				}

				// cleaning average frame :
				aver_frame->SetData(0);
}


int CCDUtil::CalcImageShift( const char* szImage1, const char* szImage2, CCDMatrix& dark,
										double& dx,double& dy,double& sigma_dx,double& sigma_dy,
										BOOL_T bVerb )
{
   CCDMatrix image1( 0 , 0 , FALSE, NOT_DEFINED, NULL, 0, FALSE, NULL, TRUE);
   CCDMatrix image2( 0 , 0 , FALSE, NOT_DEFINED, NULL, 0, FALSE, NULL, TRUE);


	if( !image1.ReadFITSFile( szImage1 , FALSE ) ){
		printf("could not read image : %s\n",szImage1);
		exit(-1);
	}	
	if( !image2.ReadFITSFile( szImage2 , FALSE ) ){
		printf("could not read image : %s\n",szImage2);
		exit(-1);
	}	

	image1.Subtract( dark );
	image2.Subtract( dark );

	image1.Laplace();
	image2.Laplace();

	Area2DInfo info;
	double lapm1,laps1;
	image1.GetVariableMeanAndSigma( gCCDParams.m_eLaplaceType,
											  lapm1, laps1, info, 20, 20,
											  (image1.GetXSize()-20), (image1.GetYSize()-20),
											  0, image1.get_laplace_data_fast() );
	printf("Laplace=%d mean=%.2f sigma=%.2f\n",gCCDParams.m_eLaplaceType,lapm1,laps1);

	double lapm2,laps2;
	image2.GetVariableMeanAndSigma( gCCDParams.m_eLaplaceType,
											  lapm2, laps2, info, 20, 20,
											  (image2.GetXSize()-20), (image2.GetYSize()-20), 
                                   0, image2.get_laplace_data_fast() );
	printf("Laplace=%d mean=%.2f sigma=%.2f\n",gCCDParams.m_eLaplaceType,lapm2, laps2);
											  

	time_t start=get_dttm();	
	CPixelAnalyseOut out;
	CPixelAnalyseIn in;
	CPixelList pixel_list(image1.GetXSize()*image1.GetYSize());
	CPixelList alreadyUsed(image1.GetXSize()*image1.GetYSize());
	in.pPixelList = &alreadyUsed;

	double n_sigma_min=10,n_sigma_max=15;
	vector<cStarCat> star_list1, star_list2;
	CPiPhotometry::GetFastPhotoList( image1, star_list1, out, in, pixel_list,
				  							   n_sigma_min,n_sigma_max,5,laps1, 20 );	
	pixel_list.Clear();
	alreadyUsed.Clear();
	CPiPhotometry::GetFastPhotoList( image2, star_list2, out, in, pixel_list,
				  							   n_sigma_min,n_sigma_max,5,laps2, 20 );

	printf("Number of found stars on image1 = %d\n",star_list1.size());
	printf("Number of found stars on image2 = %d\n",star_list2.size());
							
	int nStarCount=20;
   double mag_min=-1000;
   double mag_max=1000;
   int min_x=-1,min_y=-1,max_x=10000,max_y=10000;
   int min_dist=50;
   double delta_mag=0.5;
	int match_count = CPiPhotometry::CalcImageShift( star_list1, star_list2, 
											 dx, dy, sigma_dx, sigma_dy,
											 nStarCount, mag_min, mag_max, 
											 min_x, max_x, min_y, max_y, min_dist, 
											 delta_mag, bVerb, 2 );

	return 1;	
}