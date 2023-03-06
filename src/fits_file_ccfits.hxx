#include "fits_file_ccfits.h"
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
#include "fitsout_trace.h"

// #include <fitsio.h>
// #include <config.h>
#include <CCfits>
#include <cmath>
#include <vector>
#include <string>
// The library is enclosed in a namespace.

using namespace CCfits;
using namespace std;


template<class ARG_TYPE>
CFITSFileCCFITS<ARG_TYPE>::CFITSFileCCFITS()
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
BOOL_T CFITSFileCCFITS<ARG_TYPE>::WriteToFITSFile( void* data, long SizeX, 
                                  long SizeY, mystring& szError, 
                                  const char* fname, eFITSCompressionType compression )
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

	std::auto_ptr<FITS> pFits(0);
	// FITS* pFits = new FITS(0);
	// FITS pFits(0);

    try
    {
        // overwrite existing file if the file already exists.
        // Create a new FITS object, specifying the data type and axes for the
        // image. Simultaneously create the corresponding file.
        // this image is unsigned short data, demonstrating the cfitsio extensi
        // to the FITS standard.

        pFits.reset( new FITS(szFName.c_str() , type , naxis , naxes ) );
			// pFits.reset( new FITS(szFName.c_str()  , SHORT_IMG , naxis , naxes ) );
    }
    catch (FITS::CantCreate)
    {
          // ... or not, as the case may be.
		    gFITSLock.UnLock();
          return FALSE;
    }
	
	// long fpixel(1);
	ARG_TYPE fpixel(1);

	_TRACE_PRINTF_4("Copying header ...");
	SortHeaderBySections();
	for(int i=0;i<m_HduList.GetCount();i++){
		const char* szKey     = m_HduList[i].szName.c_str();
		const char* szKeyVal  = m_HduList[i].szValue.TrimApostrophs().c_str();
		const char* szComment = m_HduList[i].szComment.c_str();
		if(!szComment || szComment[0]=='\0'){
			szComment = GetKeyComment( m_HduList[i].szName.c_str() );
		}

		string szKeyTmp = szKey; 
		string szCommentTmp = szComment;

	 	pFits->pHDU().addKey( szKeyTmp , szKeyVal , szCommentTmp );
	}
	_TRACE_PRINTF_4("OK\n");
	
	_TRACE_PRINTF_4("Copying data ...");fflush(0);
	std::valarray<ARG_TYPE> array(nelements);
	ARG_TYPE* pData = (ARG_TYPE*)data;
	for(int j=0;j<SizeY;j++){		
		memcpy( &array[j*SizeX], pData+j*SizeX, ( SizeX*sizeof(ARG_TYPE)) );
		// printf("copying row y=%d, value=%d\n",j,(long)array[j*SizeX]);
	}
	_TRACE_PRINTF_4("OK\n");fflush(0);


	// pFits->pHDU().addKey("EXPOSURE", "xxx","Total Exposure Time");

	_TRACE_PRINTF_4("writing FITS ...");fflush(0);
	pFits->pHDU().write(fpixel,nelements,array);
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
BOOL_T CFITSFileCCFITS<ARG_TYPE>::WriteToFITSFile( void* data, 
										    long low_x, long low_y, long up_x, long up_y,
											 long FrameSizeX,long FrameSizeY,
                                  mystring& szError, 
                                  const char* fname, 
											 eFITSCompressionType compression )
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

	long SizeX = (up_x - low_x + 1);
	long SizeY = (up_y - low_y + 1);

	if( (SizeX%2)!=0 ){
		printf("SizeX must be even to save - resizing window !!!\n");fflush(0);
		up_x--;
		Assert(up_x>low_x,"Cannot re-size window to save");
		SizeX = (up_x - low_x + 1);
	}
	if( (SizeY%2)!=0 ){
		printf("SizeY must be even to save - resizing window !!!\n");fflush(0);
		up_y--;
		Assert(up_y>low_y,"Cannot re-size window to save");
		SizeY = (up_y - low_y + 1);
	}

	naxes[0] = SizeX;
	naxes[1] = SizeY;	
	long nelements = SizeX*SizeY;

	
	nBits = 8*elem_size;

	std::auto_ptr<FITS> pFits(0);
	// FITS* pFits = new FITS(0);
	// FITS pFits(0);

    try
    {
        // overwrite existing file if the file already exists.
        // Create a new FITS object, specifying the data type and axes for the
        // image. Simultaneously create the corresponding file.
        // this image is unsigned short data, demonstrating the cfitsio extensi
        // to the FITS standard.

        pFits.reset( new FITS(szFName.c_str() , type , naxis , naxes ) );
			// pFits.reset( new FITS(szFName.c_str()  , SHORT_IMG , naxis , naxes ) );
    }
    catch (FITS::CantCreate)
    {
          // ... or not, as the case may be.
			 gFITSLock.UnLock();
          return FALSE;
    }
	
	// long fpixel(1);
	ARG_TYPE fpixel(1);

	SortHeaderBySections();
	for(int i=0;i<m_HduList.GetCount();i++){
		const char* szKey     = m_HduList[i].szName.TrimApostrophs().c_str();
		const char* szKeyVal  = m_HduList[i].szValue.c_str();
		const char* szComment = m_HduList[i].szComment.c_str();
		if(!szComment || szComment[0]=='\0'){
			szComment = GetKeyComment( m_HduList[i].szName.c_str() );
		}

		string szKeyTmp = szKey; 
		string szCommentTmp = szComment;
		
		if(!IsStandardKey( szKeyTmp.c_str() ) ){
		 	pFits->pHDU().addKey( szKeyTmp , szKeyVal , szCommentTmp );
		}
	}
	
	std::valarray<ARG_TYPE> array(nelements);
	ARG_TYPE* pData = (ARG_TYPE*)data;
	ARG_TYPE* start_pos = pData+low_y*FrameSizeX+low_x;
	register int line_size = SizeX*sizeof(ARG_TYPE);
	for(int j=0;j<SizeY;j++){		
		ARG_TYPE* pLineStart = start_pos+j*FrameSizeX;
		memcpy( &array[j*SizeX], pLineStart, line_size );
	}

	// pFits->pHDU().addKey("EXPOSURE", "xxx","Total Exposure Time");
	pFits->pHDU().write(fpixel,nelements,array);

	gFITSLock.UnLock();
	return TRUE;
}


template<class ARG_TYPE>
BOOL_T CFITSFileCCFITS<ARG_TYPE>::ReadFITSFile( Table2D<ARG_TYPE>& matrix, mystring& szError,
					                           const char* fname )
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

	//keynames = m_HeaderKeys;
	// image.readKeys(keynames,vals_string);
	/*CSafeKeyTab& keyTab = matrix.GetKeyTab();
	keyTab.Clear();
	int i=0;
   for(i=0;i<vals_string.size();i++){
		keyTab.Add( m_HeaderKeys[i].data(),vals_string[i].data() );
   }*/
	fitsfile* fitsPointer = image.fitsPointer();
	CSafeKeyTab& keyTab = matrix.GetKeyTab();
	getHeader( keyTab, fitsPointer );

   image.read(contents);

	long ax1(image.axis(0));
   long ax2(image.axis(1));

	matrix.Alloc( ax1, ax2 );
	ARG_TYPE* pData = matrix.get_data_buffer();
	long line_size = ax1*sizeof(ARG_TYPE);
	long pos=0;
	for (register long j = 0; j < ax2; j++){
		memcpy( pData+pos, &contents[pos], line_size );
		pos += ax1;
   }


	gFITSLock.UnLock();
   return TRUE;
}


template<class ARG_TYPE>
BOOL_T CFITSFileCCFITS<ARG_TYPE>::ReadFITSHeader( CSafeKeyTab& keys, mystring& szError, const char* fname )
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

