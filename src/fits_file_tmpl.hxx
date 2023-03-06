#include "fits_file_tmpl.h"
#include <ccd_fits_header_defs.h>
#include "fitslib_globals.h"
#include <mymacros.h>

#include <iostream.h>
#include <stdlib.h>

#include <mytypes.h>
#include <mystring.h>
#include <myparser.h>
#include <mystrtable.h>
#include <myfile.h>
#include <cexcp.h>
#include <basestructs.h>
#include <cfg.h>

// #include <fitsio.h>
// #include <config.h>
#include <cmath>
#include <vector>
#include <string>
// The library is enclosed in a namespace.

using namespace CCfits;
using namespace std;


template<class ARG_TYPE>
CFITSFileTemplate<ARG_TYPE>::CFITSFileTemplate()
{
}


//------------------------------------------------------------
//
//   ATTENTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
//  function from CCFitsio cannot save frames of ODD sizes -
//  they simply coredump ??????????!!!!!!!!!!!!!!!!!!
//
//
//
//------------------------------------------------------------
template<class ARG_TYPE>
BOOL_T CFITSFileTemplate<ARG_TYPE>::WriteToFITSFile( void* data, long SizeX, 
                                  long SizeY, mystring& szError, 
                                  const char* fname )
{
	gFITSLock.Lock();

	szError = "";
	szFName = "";
	// szFName << "!" << fname;
	szFName << fname;
	szFName.env2str();
	MyFile::CreateDir(szFName.c_str());
 	int nBits=0;
	int elem_size = sizeof(ARG_TYPE);
	int type = get_fits_format( (ARG_TYPE)0 );
	long naxis=2;
   long naxes[2];

	naxes[0] = SizeX;
	naxes[1] = SizeY;
	long nelements = SizeX*SizeY;	
	nBits = 8*elem_size;
	createFITS( szFName.c_str(), SizeX, SizeY, nBits );
	
	_TRACE_PRINTF_4("Copying header ...");
	for(int i=0;i<m_HduList.GetCount();i++){
		const char* szKey     = m_HduList[i].szName.c_str();
		const char* szKeyVal  = m_HduList[i].szValue.c_str();
		const char* szComment = m_HduList[i].szComment.c_str();
		if(!szComment || szComment[0]=='\0'){
			szComment = GetKeyComment( m_HduList[i].szName.c_str() );
		}

		string szKeyTmp = szKey; 
		string szCommentTmp = szComment;

		insertHeaderKeyword( i, (char*)m_HduList[i].szName.c_str() ,(char*)m_HduList[i].szValue.c_str() , "" );
	 	// pFits->pHDU().addKey( szKeyTmp , szKeyVal , szCommentTmp );
	}
	_TRACE_PRINTF_4("OK\n");
	

	_TRACE_PRINTF_4("writing FITS ...");fflush(0);
	saveFITS( data, SizeX, SizeY );
	close();
	_TRACE_PRINTF_4("OK\n");fflush(0);

	gFITSLock.UnLock();

	return TRUE;	
}

//------------------------------------------------------------
//
//   ATTENTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
//  function from CCFitsio cannot save frames of ODD sizes -
//  they simply coredump ??????????!!!!!!!!!!!!!!!!!!
//
//
//
//------------------------------------------------------------
template<class ARG_TYPE>
BOOL_T CFITSFileTemplate<ARG_TYPE>::WriteToFITSFile( void* data, 
										    long low_x, long low_y, long up_x, long up_y,
											 long FrameSizeX,long FrameSizeY,
                                  mystring& szError, 
                                  const char* fname )
{
	gFITSLock.Lock();
	
	szError = "";
	szFName = "";
	// szFName << "!" << fname;
	szFName << fname;
	szFName.env2str();
	MyFile::CreateDir(szFName.c_str());
 	int nBits=0;
	int elem_size = sizeof(ARG_TYPE);
	int type = get_fits_format( (ARG_TYPE)0 );
	long naxis    =   2;
   long naxes[2];
	nBits = 8*elem_size;

	long SizeX = (up_x - low_x + 1);
	long SizeY = (up_y - low_y + 1);

	if( (SizeX%2)!=0 ){
		printf("SizeX must be even to save - resizing window !!!\n");
		up_x--;
		Assert(up_x>low_x,"Cannot re-size window to save");
		SizeX = (up_x - low_x + 1);
	}
	if( (SizeY%2)!=0 ){
		printf("SizeY must be even to save - resizing window !!!\n");
		up_y--;
		Assert(up_y>low_y,"Cannot re-size window to save");
		SizeY = (up_y - low_y + 1);
	}

	naxes[0] = SizeX;
	naxes[1] = SizeY;	
	long nelements = SizeX*SizeY;

	

	createFITS( szFName.c_str(), SizeX, SizeY, nBits );
	

	
	// long fpixel(1);
	ARG_TYPE fpixel(1);

	for(int i=0;i<m_HduList.GetCount();i++){
		insertHeaderKeyword( i, (char*)m_HduList[i].szName.c_str() ,(char*)m_HduList[i].szValue.c_str() , "" );
	 	// pFits->pHDU().addKey(m_HduList[i].szName.c_str() ,m_HduList[i].szValue.c_str(),"");
	}	
	_TRACE_PRINTF_4("writing FITS ...");fflush(0);

	ARG_TYPE* pData = new ARG_TYPE[nelements];
	ARG_TYPE* start_pos = ((ARG_TYPE*)data)+low_y*FrameSizeX+low_x;
   register int line_size = SizeX*sizeof(ARG_TYPE);
   for(register int j=0;j<SizeY;j++){
      ARG_TYPE* pLineStart = start_pos+j*FrameSizeX;
      memcpy( &(pData[j*SizeX]), pLineStart, line_size );
   }


	saveFITS( pData, SizeX, SizeY );
	close();
	_TRACE_PRINTF_4("OK ...\n");fflush(0);	

	gFITSLock.UnLock();
	return TRUE;
}


template<class ARG_TYPE>
BOOL_T CFITSFileTemplate<ARG_TYPE>::ReadFITSFile( Table2D<ARG_TYPE>& matrix, mystring& szError,
					                           const char* fname )
{
	gFITSLock.Lock();

	Init();
	Clear();
	szFName = fname;
	szFName.env2str();
	// int type = get_fits_format( (ARG_TYPE)0 );
	int nBits = sizeof(ARG_TYPE)*8;

	if(!MyFile::IsReadable( szFName.c_str() )){
		printf("cannot access file : %s, exiting ...\n",szFName.c_str() );
		exit(-1);		
	}
	/*if( strstr( szFName.c_str(), ".fitc" )){
		szFName << "[compress]";
	}*/
	_TRACE_PRINTF_0("reading file %s ...\n",fname);fflush(0);

	
	BOOL_T bRet=FALSE;
	ARG_TYPE* data = matrix.get_data_buffer();
	if(open( szFName.c_str(), READONLY )){
		szError = "";
		if( load( data, matrix.GetXSize(), matrix.GetYSize(), nBits ) ){
			bRet = TRUE;
		}
	}

	gFITSLock.UnLock();
   return bRet;
}


template<class ARG_TYPE>
BOOL_T CFITSFileTemplate<ARG_TYPE>::ReadFITSHeader( CSafeKeyTab& keys, mystring& szError, const char* fname )
{
	gFITSLock.Lock();

	Init();
	Clear();
	szFName = fname;
	szFName.env2str();

	if(!MyFile::IsReadable( szFName.c_str() )){
		printf("cannot access file : %s, exiting ...\n",szFName.c_str() );
		exit(-1);		
	}
	_TRACE_PRINTF_0("reading file %s ...\n",fname);fflush(0);

	szError = "";
	std::auto_ptr<FITS> pInfile(new FITS(szFName.c_str(),Read,true));

	PHDU& image = pInfile->pHDU();
   std::valarray<ARG_TYPE>  contents;

   // read all user-specifed, coordinate, and checksum keys in the image
   image.readAllKeys();

	// std::cout << image << std::endl;		
	
	std::vector<string> keynames;
	std::vector<string> vals_string;

	fitsfile* fitsPointer = image.fitsPointer();

	keys.Clear();
	getHeader( keys, fitsPointer );
	
	gFITSLock.UnLock();
	return TRUE;	
}

