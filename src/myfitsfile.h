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
#ifndef _MYFITSFILE_H__
#define _MYFITSFILE_H__


#include "fits_file_base.h"
#include "fits_file_standard.h"
#include "ccd_fits_header_defs.h"
#include "tab2D.h"

// enum eFITSCompressionType { eFITSComprNone=0, eFITSComprASAS };

template<class ARG_TYPE>
class CMyFITSFile : public CFITSFileBase
{
protected :
	int fd;
	int m_HeaderKeyCount;
	int m_HeaderSize;

	char* m_pHeader;
	int m_AllocatedSize;

	// compresion RICE :
	int m_ComprBlockSize;
	int m_ComprDivison;
	   
	
	BOOL_T Open( const char* fname, BOOL_T bWrite=FALSE );

	int get_default_bzero();
	float get_default_bscale();

	// calculation 
	BOOL_T AddBzero( void* data, int SizeX, int SizeY, int nBits, int bzero );

	// internal writing routines :
	const char* GetHeaderComment( const char* key );
	
	void ClearHeader();
	
	BOOL_T WriteHeaderLine( const char* szName, 
									const char* szValue, const char* szComment );
	
	BOOL_T WriteHeaderLine( const char* szName, int Value, const char* szComment );
	BOOL_T WriteHeaderLine( const char* szName, double Value, const char* szComment );

	BOOL_T WriteStandardHeader( int SizeX, int SizeY, int bzero, float bscale );
	
	BOOL_T AddHeaderKeys( CSafeKeyTab& keyTab, BOOL_T bWriteAll=FALSE );		
	
	BOOL_T AddAsasComprKey();
	BOOL_T AddAsasComprKeys();
	
	BOOL_T AddCommentLines( const char* comment="" );
	
	BOOL_T PadToBlockSize();
	
	BOOL_T WriteDataCompressed( void* data, int SizeX, int SizeY, int bz );		
	
	// internal read routines :
	BOOL_T ReadHeader( int& xSize, int& ySize,
	                   int& bitpix, float& bscale, float& bzero);

	BOOL_T ReadAsasCompressedData( void* data, int xSize, int ySize, float& bz, float& bs );
	BOOL_T GetCompressionKeys( int& blocksize, int& divisor, mystring& szComprType );

public :	
	int getFileDesc(){ return fd; }

	void Close();
	
	// function automaticaly checks if file is compressed :	                   
	BOOL_T ReadData( void* data, int xSize, int ySize, float& bzero, float& bscale );
	


	CMyFITSFile( const char* fname="" );
	~CMyFITSFile();


	// creating new FITS-like file :
	BOOL_T WriteFITSHeader( const char* fname,
									int SizeX, int SizeY, CSafeKeyTab& keyTab,
									eFITSCompressionType compr=eFITSComprNone );

	BOOL_T WriteToFile( void* data, int xSize, int ySize,
							  int size, CSafeKeyTab& keyTab, 
							  const char* fname, mystring& szError, 
							  eFITSCompressionType compr=eFITSComprNone );



	// creating new FITS file :	
	BOOL_T WriteToFITSFile( void* data, int SizeX, int SizeY, 
									mystring& szError, 
									eFITSCompressionType compr=eFITSComprNone );

	BOOL_T WriteToFITSFile( void* data, int SizeX, int SizeY, 
									mystring& szError, const char* fname,
									eFITSCompressionType compr=eFITSComprNone );

	BOOL_T WriteToFITSFile( void* data, int SizeX, int SizeY, 
									CSafeKeyTab& keyTab, mystring& szError,								
									const char* fname,											
									eFITSCompressionType compr=eFITSComprNone );
									
	 BOOL_T WriteToFITSFile( void* data,
	 								 long low_x, long low_y, long up_x, long up_y,
	 								 long FrameSizeX, long FrameSizeY, mystring& szError,
	 								 const char* fname,
	 								 eFITSCompressionType compr=eFITSComprNone );			
									
	// reading FITS file :
	BOOL_T ReadFITSFile( Table2D<ARG_TYPE>& matrix, mystring& szError,
								const char* fname );
								
	BOOL_T ReadFITSFile( Table2D<ARG_TYPE>& matrix, CSafeKeyTab& keyTab,
								mystring& szError, const char* fname );

	int ReadFITSFileForAsasPipeline( ARG_TYPE* data, int sizeX, int sizeY, 
													const char* fname );
	
	BOOL_T ReadFITSHeader( CSafeKeyTab& keyTab, const char* fname );

	mystring getKeyValue( const char* fname, const char* key );

	int ReadRaw( const char* fname, CSafeKeyTab& keyTab, char* data, int file_size );
	int WriteRaw( const char* fname, CSafeKeyTab& keyTab, char* data, int file_size );
};


#endif
