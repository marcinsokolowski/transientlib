#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "basestructs.h"
#include "ccd_fits_header_defs.h"
#include "myutil.h"
#include "mystrtable.h"
#include "myparser.h"
#include "myfile.h"
#include "fitslib_globals.h"
#include "mymacros.h"
#include "cexcp.h"
#include <errno.h>
extern "C" {
#include "asas_fitsio_exp.h"
}

/*extern "C" {
int write_fits_data( int fd,void* data,int nx,int ny,int bitpix,int flip );	
int close_fits(int fd,int nx,int ny,int bitpix);
void *convert_data(void* inp,int npixel,int bpi,float bscale,float bzero,int bpo);

int fits_rcomp(int a[],int nx,unsigned char *c,int clen,int nblock,int bsize);

};*/

#include "myfitsfile.h"

template<class ARG_TYPE>
int CMyFITSFile<ARG_TYPE>::get_default_bzero()
{
	int bpp=sizeof(ARG_TYPE);
	if(bpp==2)
		return 32768;
	return 32768;
}


template<class ARG_TYPE>
float CMyFITSFile<ARG_TYPE>::get_default_bscale()
{
	return 1;
}



template<class ARG_TYPE>
CMyFITSFile<ARG_TYPE>::CMyFITSFile( const char* fname )
: fd(-1), m_HeaderKeyCount(0), m_HeaderSize(0), m_ComprBlockSize(BLOCKSIZE),
	m_ComprDivison(DIVISOR)  
{
	szFName = fname;
	m_AllocatedSize = MAX_HEADER_KEY_COUNT*FITS_LINE_SIZE;
	m_pHeader = new char[ m_AllocatedSize ];	
	memset( m_pHeader, m_AllocatedSize , '\0' );
}

template<class ARG_TYPE>
void CMyFITSFile<ARG_TYPE>::ClearHeader()
{
	memset(m_pHeader,m_AllocatedSize , '\0' );
	m_HeaderKeyCount=0;
	m_HeaderSize=0;
}

template<class ARG_TYPE>
BOOL_T CMyFITSFile<ARG_TYPE>::Open( const char* fname, BOOL_T bWrite )
{
	if( bWrite ){
		if ( (fd = open(fname, O_WRONLY|O_TRUNC|O_CREAT, FITS_ACCESS_MODE)) <= 0){
			printf("could not open file : %s in WRITE mode !\n",fname);
			return FALSE;
		}
	}else{
		if ( (fd = open(fname, O_RDONLY, FITS_ACCESS_MODE)) <= 0){
			printf("could not open file : %s in READ mode !\n",fname);
			return FALSE;
		}
	}
	return TRUE;
}

template<class ARG_TYPE>
void CMyFITSFile<ARG_TYPE>::Close()
{
	if(fd > 0 ){
		close( fd );
		fd = -1;
	}
}

template<class ARG_TYPE>
CMyFITSFile<ARG_TYPE>::~CMyFITSFile()
{
	if(fd>0){
		Close();
	}
	if( m_pHeader ){
		delete [] m_pHeader;
	}
}

template<class ARG_TYPE>
BOOL_T CMyFITSFile<ARG_TYPE>::WriteHeaderLine( const char* szName, int Value, const char* szComment )
{
	char szTmp[64];
	sprintf(szTmp,"%d",Value);
	return WriteHeaderLine( szName, szTmp, szComment );
}

template<class ARG_TYPE>
BOOL_T CMyFITSFile<ARG_TYPE>::WriteHeaderLine( const char* szName, double Value, const char* szComment )
{
	char szTmp[64];
	sprintf(szTmp,"%.8f",Value);
	return WriteHeaderLine( szName, szTmp, szComment );
}

template<class ARG_TYPE>
BOOL_T CMyFITSFile<ARG_TYPE>::WriteHeaderLine( const char* szName, 
												 const char* szValue, const char* szComment )
{
	char line[FITS_LINE_SIZE+20];
	mystring szValOut;
	if(strlen(szValue) && !mystring::IsInApostrophs(szValue) 
		&& GetKeyType(szName)==eSTRING && !IsStandardKey(szName))
		szValOut << "'" << szValue << "'";
	else
		szValOut << szValue;
	int val_len=strlen(szValOut.c_str());
	int total_len=mystrlen(szName)+val_len+mystrlen(szComment)+4;
	if( strcmp( szName, FH_COMMENT )==0 ){
		total_len=mystrlen(szName)+val_len+3;
	}
	if(total_len>FITS_LINE_SIZE){
		printf("ERROR : WriteHeaderLine total_len>FITS_LINE_SIZE !\n");
		return FALSE;
	}

	strcpy(line,szName);
	if(strcmp(szName,"END") != 0){	
		fits_file_extend(line,' ',8,FITS_LINE_SIZE);	
		if(strcmp(szName,FH_COMMENT)){
			strcat(line,"= ");
		}else{
			strcat(line,"  ");
		}
		fits_file_extend(line,' ',30-val_len,FITS_LINE_SIZE);
		strcat(line,szValOut.c_str());
		// int total_len=strlen(line);
		fits_file_extend(line,' ',31,FITS_LINE_SIZE);

		if( strcmp(szName,FH_COMMENT)){
			strcat(line,"/ ");
			strcat(line,szComment);
		}
	}

	fits_file_extend(line,' ',80,FITS_LINE_SIZE);

	if(m_HeaderKeyCount<MAX_HEADER_KEY_COUNT){
		char* startPos = m_pHeader+m_HeaderKeyCount*FITS_LINE_SIZE;
		memcpy( startPos, line, FITS_LINE_SIZE );
		m_HeaderKeyCount++;
		m_HeaderSize += FITS_LINE_SIZE;
		return TRUE;		
	}

	printf("max number of keys (%d) exceeded : %d\n",MAX_HEADER_KEY_COUNT,m_HeaderKeyCount);
	return FALSE;	
}

template<class ARG_TYPE>
const char* CMyFITSFile<ARG_TYPE>::GetHeaderComment( const char* key )
{
	CEnvVar* pKey = m_HduList.Find( key );
	if( pKey && strlen(pKey->szComment.c_str()) ){
		return pKey->szComment.c_str();
	}
	const char* tmp = GetKeyComment( key );
	if(tmp && tmp[0])
		return tmp;
	return "";
}

template<class ARG_TYPE>
BOOL_T CMyFITSFile<ARG_TYPE>::WriteStandardHeader( int SizeX, int SizeY, 
																	int bzero, float bscale )
{
	ARG_TYPE argtype=0;
	int bitpix_sign = GetBitPixSign( argtype );
	int nBits=8*sizeof(ARG_TYPE)*bitpix_sign;
	int naxis=2;

	

	if(!WriteHeaderLine( FH_SIMPLE, "T", GetHeaderComment(FH_SIMPLE) ))
		return FALSE;
	if(!WriteHeaderLine( FH_BITPIX, nBits , GetHeaderComment(FH_BITPIX) ))
		return FALSE;
	if(!WriteHeaderLine( FH_NAXIS, naxis, GetHeaderComment(FH_NAXIS) ))
		return FALSE;
	if(!WriteHeaderLine( FH_NAXIS1, SizeX, GetHeaderComment(FH_NAXIS1) ))
		return FALSE;
	if(!WriteHeaderLine( FH_NAXIS2, SizeY, GetHeaderComment(FH_NAXIS2) ))
		return FALSE;
	if(!WriteHeaderLine( FH_EXTEND, "T" , GetHeaderComment(FH_EXTEND) ))
		return FALSE;
	if( bitpix_sign>= 0){
		if(!WriteHeaderLine( FH_BZERO, bzero, GetHeaderComment(FH_BZERO) ))
			return FALSE;
		if(!WriteHeaderLine( FH_BSCALE, bscale , GetHeaderComment(FH_BSCALE) ))
			return FALSE;
	}
	if(!WriteHeaderLine( FH_COMMENT, STD_COMMENT_LINE_1, NULL )){
		return FALSE;
	}
	if(!WriteHeaderLine( FH_COMMENT, STD_COMMENT_LINE_2, NULL )){
		return FALSE;
	}
	
	return TRUE;
}

template<class ARG_TYPE>
BOOL_T CMyFITSFile<ARG_TYPE>::AddAsasComprKey()
{
	if(!WriteHeaderLine( COMPRES, "RICE", "compresion type" ))
		return FALSE;
	if(!WriteHeaderLine( BLOCKSZ, m_ComprBlockSize, "" ))
		return FALSE;
	if(!WriteHeaderLine( FH_DIVISOR, m_ComprDivison, "" ))
		return FALSE;		
	
	return TRUE;
}

template<class ARG_TYPE>
BOOL_T CMyFITSFile<ARG_TYPE>::AddAsasComprKeys()
{
	m_HduList.Set( FITS_COMPRESS_KEY, "RICE", fitsHeaderKeyDefTab[eCOMPRES].szKeyComment );
	m_HduList.Set( BLOCKSZ, m_ComprBlockSize, fitsHeaderKeyDefTab[eBLOCKSZ].szKeyComment );
	m_HduList.Set( FH_DIVISOR, m_ComprDivison, fitsHeaderKeyDefTab[eDIVISOR].szKeyComment );
	
	return TRUE;
}


template<class ARG_TYPE>
BOOL_T CMyFITSFile<ARG_TYPE>::AddHeaderKeys( CSafeKeyTab& keyTab, BOOL_T bWriteAll ){
	for(int i=0;i<keyTab.GetCount();i++){		
		const char* szKey     = keyTab[i].szName.c_str();
		if(!IsStandardKey(szKey) || bWriteAll ){
	      const char* szKeyVal  = keyTab[i].szValue.c_str();
   	   const char* szComment = keyTab[i].szComment.c_str();
      	if(!szComment || szComment[0]=='\0'){
	         szComment = GetKeyComment( keyTab[i].szName.c_str() );
   	   }	
			if(!WriteHeaderLine( szKey, szKeyVal, szComment )){
				printf("error writing header key : %s=%s (%s)\n",szKey, szKeyVal, szComment);
				// exit(0);
				return FALSE;
			}
		}
	}
	return TRUE;
}

template<class ARG_TYPE>
BOOL_T CMyFITSFile<ARG_TYPE>::AddCommentLines( const char* comment )
{
	// we want to fill place to blocksize (2880 = 36*80) with comment lines :	
	int total_lines=m_HeaderKeyCount+1;
	if((total_lines%36)!=0){
		// unless by accidente we match block size :
		int block_count=(total_lines/36)+1;
		int lines_count=(block_count*36)-1; // minus 1 due to END line :
		for(int i=m_HeaderKeyCount;i<lines_count;i++){
			if(!WriteHeaderLine( "COMMENT", comment, "" ))
				return FALSE;
		}	
	}
	return TRUE;
}

template<class ARG_TYPE>
BOOL_T CMyFITSFile<ARG_TYPE>::PadToBlockSize()
{
	int block = (m_HeaderSize / FITS_HEADER_BLOCK_SIZE);
	int total_size=(block+1)*FITS_HEADER_BLOCK_SIZE;
	while(m_HeaderSize<total_size){
		if(m_HeaderSize>=m_AllocatedSize)
			return FALSE;
		m_pHeader[m_HeaderSize++]='\0';		
	}
	return TRUE;
}

template<class ARG_TYPE>
BOOL_T CMyFITSFile<ARG_TYPE>::ReadHeader( int& xSize, int& ySize,
														int& bitpix, float& bscale,
														float& bzero )
{
	BOOL_T bEnd=FALSE;
	char line[FITS_LINE_SIZE+20];
	char sName[FITS_LINE_SIZE];
	char sValue[FITS_LINE_SIZE];
	char sComment[FITS_LINE_SIZE];

	xSize=0;
	ySize=0;
	bitpix=16;
	bscale=1.00;
	bzero=0.00;

	int valLength=(30-(FIT_KEY_LEN+2));

	// clear key tab :
	m_HduList.Clear();
	char* ptr=m_pHeader;
	CMyStrTable items;
	int total_bytes_read=0;

	while( !bEnd ){
		if( ( m_AllocatedSize - total_bytes_read ) < FITS_HEADER_BLOCK_SIZE ){
			printf("Header to long already %d bytes read, allocated was %d !\n",total_bytes_read,m_AllocatedSize);
			printf("Cannot read next block of size %d bytes !\n",FITS_HEADER_BLOCK_SIZE);
			printf("Only %d bytes remains in buffer !\n",(m_AllocatedSize-total_bytes_read));
			return FALSE;
		}

		int read_bytes = read( fd , ptr, FITS_HEADER_BLOCK_SIZE );
		if(read_bytes!=FITS_HEADER_BLOCK_SIZE){
			printf("Incorrect block size of header !!!\n");
			return FALSE;
		}

		char* line_ptr=ptr;
		for(register int i=0;i<FITS_HEADER_LINES && !bEnd;i++){			
			strncpy( sName, line_ptr, FIT_KEY_LEN );
			sName[FIT_KEY_LEN]='\0';

			strncpy( sValue, line_ptr+10, (FITS_LINE_SIZE-10) );
			sValue[FITS_LINE_SIZE-10]='\0';

			/*strncpy( sComment, line_ptr+33, (FITS_LINE_SIZE-32) );
			sComment[FITS_LINE_SIZE-32]='\0';*/
			
			MyParser parsName=sName;
			mystring szName = parsName.GetNextItem();

			if(strcmp(szName.c_str(),FITS_COMMENT_KEY)){
				MyParser parsValue=sValue;
				MyParser szUpToComment=parsValue.GetItemUpTo("/");
				mystring szValueTmp = szUpToComment.TrimApostro();
				mystring szValue = szValueTmp.TrimApostro(' ');

				MyParser parsComment;
				char* pCommentStart=strstr(sValue,"/");
				if(pCommentStart){
					parsComment = pCommentStart+1;
				}				
				mystring szComment = parsComment.TrimApostro();
			

				if(strcmp(szName.c_str(),END_KEYWORD)==0){
					bEnd = TRUE;
				}else{
					m_HduList.Add( szName.c_str(), szValue.c_str(), szComment.c_str() );
				}
			}
			line_ptr += FITS_LINE_SIZE;						
		}
		ptr += FITS_HEADER_BLOCK_SIZE;
		total_bytes_read += read_bytes;
	}

	// now fill standard values :
	const char* szTMP=m_HduList.getKeyVal( FH_BITPIX );
	if(szTMP && szTMP[0])
		bitpix = atol( szTMP );
	szTMP = m_HduList.getKeyVal( FH_NAXIS1 );
	if( szTMP && szTMP[0])
		xSize = atol( szTMP );
	szTMP = m_HduList.getKeyVal( FH_NAXIS2 );
   if( szTMP && szTMP[0])
      ySize = atol( szTMP );
	szTMP = m_HduList.getKeyVal( FH_BSCALE );
	if( szTMP && szTMP[0])
		bscale = atof( szTMP );
	szTMP = m_HduList.getKeyVal( FH_BZERO );
	if( szTMP && szTMP[0])
		bzero = atof( szTMP );

	return TRUE;
}

template<class ARG_TYPE>
BOOL_T CMyFITSFile<ARG_TYPE>::GetCompressionKeys( int& blocksize, int& divisor, mystring& szComprType )
{
	const char* szTmp = 	m_HduList.getKeyVal( FITS_COMPRESS_KEY );
	if(szTmp && szTmp[0])
		szComprType = szTmp;
	else
		return FALSE;

	szTmp =  m_HduList.getKeyVal( FITS_COMPRESS_BLOCKSIZE );
	if(szTmp && szTmp[0])
		blocksize = atol( szTmp );
	else
		return FALSE;

	szTmp =  m_HduList.getKeyVal( FITS_COMPRESS_DIVISOR );
	if(szTmp && szTmp[0])
		divisor = atol( szTmp );
	else
		return FALSE;

	m_HduList.Delete( FITS_COMPRESS_KEY );
	m_HduList.Delete( FITS_COMPRESS_BLOCKSIZE );
	m_HduList.Delete( FITS_COMPRESS_DIVISOR );			
	return TRUE;
}

template<class ARG_TYPE>
BOOL_T CMyFITSFile<ARG_TYPE>::WriteFITSHeader( const char* fname,
															  int SizeX, int SizeY, CSafeKeyTab& keyTab,
															  eFITSCompressionType compr )
{
		if(fd<=0)
			return FALSE;

		if( compr==eFITSComprASAS ){
      	AddAsasComprKeys();
      }

		SortHeaderBySections();
		ClearHeader();

		const char* szBZERO = keyTab.getKeyVal( FH_BZERO );
		const char* szBSCALE = keyTab.getKeyVal( FH_BSCALE );
		int bzero=0;
		float bscale=1;
		if(szBZERO && szBZERO[0]){
			bzero = atol( szBZERO );
		}else{
			if( sizeof( ARG_TYPE )==2 ){
				bzero = USHORT_BZERO;
			}
		}
		if(szBSCALE && szBSCALE[0]){
			bscale = atof( szBSCALE );
		}else{
			if(sizeof( ARG_TYPE )==2){
				bscale = USHORT_BSCALE;
			}
		}

		if(!WriteStandardHeader( SizeX, SizeY, bzero, bscale )){
			Close();
			printf("Error writing in WriteStandardHeader\n");
			return FALSE;
		}
		if(!AddHeaderKeys( keyTab )){
			Close();
			printf("Error writing in AddHeaderKeys\n");
			return FALSE;
		}
		/*if( compr==eFITSComprASAS ){
			if(!AddAsasComprKey()){
				Close();
				return FALSE;
			}
		}*/
		if(!AddCommentLines()){
			Close();
			printf("Error writing in AddCommentLines\n");
			return FALSE;
		}
		if(!WriteHeaderLine( "END" ,"" ,"" )){
			Close();
			printf("Error writing END to FITS header\n");
			return FALSE;
		}

		if(fd>0){
			write(fd, m_pHeader,m_HeaderSize);
		}

		return TRUE;
}

template<class ARG_TYPE>
BOOL_T CMyFITSFile<ARG_TYPE>::WriteToFITSFile( void* data, int SizeX, int SizeY,
                    					    CSafeKeyTab& keyTab, mystring& szError,
												 const char* fname,
												 eFITSCompressionType compr )
{
	mystring szOutName = fname;
	szOutName.env2str();
	
	if(compr!=eFITSComprNone){
		if(!strstr( szOutName.c_str(), ".fitc" ))
			szOutName << "c";
	}

	MyFile::CreateDir( szOutName.c_str() );
	if( Open( szOutName.c_str(), TRUE ) ){
		if( compr==eFITSComprASAS ){
   	  	AddAsasComprKeys();
	   }
		
		SortHeaderBySections();
		ClearHeader();

		const char* szBZERO = keyTab.getKeyVal( FH_BZERO );
		const char* szBSCALE = keyTab.getKeyVal( FH_BSCALE );
		int bzero=0;
		float bscale=1;
		if(szBZERO && szBZERO[0]){
			bzero = atol( szBZERO );
		}else{
			if( sizeof( ARG_TYPE )==2 ){
				bzero = USHORT_BZERO;
			}
		}
		if(szBSCALE && szBSCALE[0]){
			bscale = atof( szBSCALE );
		}else{
			if(sizeof( ARG_TYPE )==2){
				bscale = USHORT_BSCALE;
			}
		}

		if(!WriteStandardHeader( SizeX, SizeY, bzero, bscale )){
			Close();
			return FALSE;
		}
		if(!AddHeaderKeys( keyTab )){
			Close();
			return FALSE;
		}
		/*if( compr==eFITSComprASAS ){
			if(!AddAsasComprKey()){
				Close();
				return FALSE;
			}
		}*/
		if(!AddCommentLines()){
			Close();
			return FALSE;
		}
		if(!WriteHeaderLine( "END" ,"" ,"" )){
			Close();
			return FALSE;
		}

		/*if(compr==eFITSComprNone){
			// pad to block size only for not compressed :
			if(!PadToBlockSize())
				return FALSE;
		}*/

		write(fd, m_pHeader,m_HeaderSize);
		
		int nBits = 8*sizeof(ARG_TYPE);
		/*int total_size = (SizeX*SizeY)*nBits;
		int bytes = write(fd, (char*)data, total_size );
		if( total_size!=bytes )
			return FALSE;*/

		if(compr==eFITSComprNone){
			if(bzero!=0){
				AddBzero( data, SizeX, SizeY, nBits, -bzero );
			}
			write_fits_data( fd , data, SizeX, SizeY, nBits, 0 );
			close_fits(fd , SizeX, SizeY, nBits );
			fd = -1;
			if(bzero!=0){
				// in order to keep original values now add bzero :
				AddBzero( data, SizeX, SizeY, nBits, bzero );
			}
		}else{
			if(compr==eFITSComprASAS){
				if(!WriteDataCompressed( data, SizeX, SizeY, bzero )){
					Close();
					return FALSE;
				}
			}
		}
		
		Close();
		return TRUE;
	}
	return FALSE;
}

template<class ARG_TYPE>
BOOL_T CMyFITSFile<ARG_TYPE>::WriteToFile( void* data, int xSize, int ySize,
							  int size, CSafeKeyTab& keyTab,
                       const char* fname, mystring& szError,
                       eFITSCompressionType compr )
{
	mystring szOutName = fname;
   szOutName.env2str();
                                                                                
   if( !Open( szOutName.c_str(), TRUE ) )
		return FALSE;

	if(!WriteFITSHeader( fname, xSize, ySize,  keyTab, compr ))
		return FALSE;
	
	if( size>0 && data ){
		int ret = write( fd, data, size );
		if( ret!=size )
			return FALSE;
	}

	return TRUE;		
}

template<class ARG_TYPE>
BOOL_T CMyFITSFile<ARG_TYPE>::WriteToFITSFile( void* data, 
										    long low_x, long low_y, long up_x, long up_y,
											 long FrameSizeX,long FrameSizeY,
                                  mystring& szError, 
                                  const char* fname,
											 eFITSCompressionType compr )
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

	
	ARG_TYPE* pDataPart = new ARG_TYPE[nelements];
	ARG_TYPE* pData = (ARG_TYPE*)data;
	ARG_TYPE* start_pos = pData+low_y*FrameSizeX+low_x;
	register int line_size = SizeX*sizeof(ARG_TYPE);
	for(register int j=0;j<SizeY;j++){		
		ARG_TYPE* pLineStart = start_pos+j*FrameSizeX;
		memcpy( &pDataPart[j*SizeX], pLineStart, line_size );
		/*int pos=j*SizeX;
		for(register int i=0;i<SizeX;i++){
			pDataPart[pos+i]=pLineStart[i];
		}*/
	}
	BOOL_T bRet = WriteToFITSFile( pDataPart, SizeX, SizeY, szError,
											 fname, compr );
	delete [] pDataPart;

	gFITSLock.UnLock();
	return bRet;
}


template<class ARG_TYPE>
BOOL_T CMyFITSFile<ARG_TYPE>::AddBzero( void* data, int SizeX, int SizeY, 
													 int nBits, int bzero )
{
	register int size=(SizeX*SizeY);
	if(nBits==16){
		for(register int i=0;i<size;i++){
			((unsigned short*)data)[i] += bzero;		
		}
	}
	return TRUE;
}

template<class ARG_TYPE>
BOOL_T CMyFITSFile<ARG_TYPE>::WriteDataCompressed( void* data, int SizeX, int SizeY, int bz )
{
	int bp = sizeof(ARG_TYPE)*8;
	int naxis=2, nx=SizeX, ny=SizeY, nz=1;
	char* obuf = (char *)malloc(4*nx);
   float bs = 1;
	bz=0;
   int nout = 0;

	ARG_TYPE* ptr=(ARG_TYPE*)data;
	for(register int j=0; j<ny; j++){
		// int* idata = (int*)convert_data(ptr,nx,bp,bs,bz,32);
		// int olen=fits_rcomp(idata,nx,(unsigned char*)obuf,4*nx,m_ComprBlockSize,2);
		int* idata=new int[nx];
		for(register int i=0;i<nx;i++){
			idata[i] = (int)(ptr[i]);
		}
		
		int olen=fits_rcomp(idata,nx,(unsigned char*)obuf,4*nx,m_ComprBlockSize,2);
		if(olen > 0){
			write(fd,obuf,olen);
			nout += olen;
		}else{
			printf("Compression error: %d\n",olen);
			return FALSE;
		}		
//		free( idata );
		delete [] idata;
		ptr += nx;
	}		
	_TRACE_PRINTF_5("compressed size = %d\n",nout);
	free( obuf );
	return TRUE;
}

template<class ARG_TYPE>
BOOL_T CMyFITSFile<ARG_TYPE>::ReadData( void* data, int xSize, int ySize, float& bzero, float& bscale )
{
	// find COMPRES key :
	const char* szCompres = m_HduList.getKeyVal( FITS_COMPRESS_KEY );
	eFITSCompressionType compresType = eFITSComprNone;
	if(szCompres && szCompres[0]){
		if(strcmp(szCompres,FITS_RICE_COMPRESS)==0){
			compresType=eFITSComprASAS;
		}
	}
	
	if( compresType==eFITSComprNone ){
		// no compression simply read data :
		int nBits=8*sizeof(ARG_TYPE);
		int size = (xSize*ySize*sizeof(ARG_TYPE));
		// int read_bytes = read( fd, (char*)data, size );
		int read_bytes = read_fits_data_new( fd, xSize, ySize, nBits, data );
		if(read_bytes!=size){
			printf("Error while reading FITS pixel data : %s\n",strerror(errno));
			return FALSE;
		}
		AddBzero( data, xSize, ySize, nBits, bzero );
		return TRUE;		
	}else{
		if( compresType==eFITSComprASAS ){
			if(!ReadAsasCompressedData( data, xSize, ySize, bzero, bscale ))
				return FALSE;
			AddBzero( data, xSize, ySize, sizeof(ARG_TYPE)*8, (int)(-bzero) );
			return TRUE;
		}else{
			printf("ERROR : unknown compression type : %d\n",compresType);
			return FALSE;
		}
	}
	return FALSE;
}

template<class ARG_TYPE>
BOOL_T CMyFITSFile<ARG_TYPE>::ReadAsasCompressedData( void* data, int xSize, int ySize, float& bz, float& bs )
{
	int blocksize,divisor;
	mystring szComprType;
	if(!GetCompressionKeys( blocksize, divisor, szComprType ))
		return FALSE;
	bs *= divisor;
   bz *= divisor;
	
	int n=0;
	int left = 0;
   int more = 1;
	int nx = xSize;
	int ny = ySize;
	char* obuf=(char *)malloc(4*nx);
	char* obuf_tmp=(char *)malloc(4*nx);
	short* odata = (short *)malloc(nx*sizeof(short));

	ARG_TYPE* ptr = (ARG_TYPE*)data;
	for(register int j=0; j<ny; j++){
		if(more){
			n = read (fd, obuf+left,4*nx-left);
			if(n <=0 ) {
        		more = 0;
	   		n = left;
	      }else{
 				n += left;
			}
		}else{
			n = left;
		}
		int* idata = (int *)malloc(nx*sizeof(int));
		left=fits_rdecomp((unsigned char*)obuf,n,(unsigned int*)idata,nx,blocksize,2);
		if(left < 0){
     		free(obuf);
     		free(obuf_tmp);
	      free(odata);
	      free(idata);
			return FALSE;
	   }
	   if(left){
			// NEW correction - due to overlap found be mpatrol :
			memcpy(obuf_tmp, obuf+n-left, left);
			memcpy(obuf, obuf_tmp, left);
		}
	   for(register int i=0; i<nx; i++){
      	if(idata[i] >= 65535){
	      	odata[i] = (u_short)65535-bz;
	      }else{
      		if(idata[i] > 65534 ){
					odata[i] =(short)(65534-bz);
				}else{
			 		if(idata[i] > 0 ){
						odata[i] =(short)(idata[i]+0.5-bz);
		         }else{
						odata[i] = 0-bz;
					}
				}
      	}
	   }
		memcpy( ptr, odata, nx*sizeof(short));
		ptr += nx;
   	free(idata);		
	}
	free(odata);
   free(obuf);	
	free(obuf_tmp);
	
	return TRUE;
}

template<class ARG_TYPE>
BOOL_T CMyFITSFile<ARG_TYPE>::WriteToFITSFile( void* data, int SizeX, int SizeY,
												 mystring& szError, 
												 eFITSCompressionType compr )
{
	return WriteToFITSFile( data, SizeX, SizeY, m_HduList, szError, szFName.c_str(), compr );
}

template<class ARG_TYPE>
BOOL_T CMyFITSFile<ARG_TYPE>::WriteToFITSFile( void* data, int SizeX, int SizeY,
												 mystring& szError, const char* fname,
												 eFITSCompressionType compr )
{
	return WriteToFITSFile( data, SizeX, SizeY, m_HduList, szError, fname, compr );
}

template<class ARG_TYPE>
BOOL_T CMyFITSFile<ARG_TYPE>::ReadFITSFile( Table2D<ARG_TYPE>& matrix, mystring& szError,
														  const char* fname )
{
	return ReadFITSFile( matrix, matrix.GetKeyTab(), szError, fname );
}

template<class ARG_TYPE>
BOOL_T CMyFITSFile<ARG_TYPE>::ReadFITSHeader( CSafeKeyTab& keyTab, const char* fname )
{
	mystring szInName=fname;
	szInName.env2str();

	int xSize=0,ySize=0,bpp=16;
	float bscale=1.00,bzero=0.00;


	if(Open( szInName.c_str() )){
		// first read header :
		if(!ReadHeader(xSize,ySize,bpp,bscale,bzero)){
			Close();
			return FALSE;
		}
		keyTab = m_HduList;
		return TRUE;		
	}
	return FALSE;
}


template<class ARG_TYPE>
mystring CMyFITSFile<ARG_TYPE>::getKeyValue( const char* fname, const char* key )
{
	mystring szRet;
	CSafeKeyTab keyTab;
	if( ReadFITSHeader( keyTab, fname ) ){
		const char* szVal = keyTab.getKeyVal( key );
		if( szVal && szVal[0] )
			szRet = szVal;
	}
	return szRet;
}

template<class ARG_TYPE>								
BOOL_T CMyFITSFile<ARG_TYPE>::ReadFITSFile( Table2D<ARG_TYPE>& matrix, CSafeKeyTab& keyTab,
														  mystring& szError, const char* fname )
{
	mystring szInName=fname;
	szInName.env2str();

	int xSize=0,ySize=0,bpp=16;
	float bscale=1.00,bzero=0.00;

	int x_matrix_size = matrix.GetXSize();
	int y_matrix_size = matrix.GetYSize();


	if(Open( szInName.c_str() )){
		// first read header :
		if(!ReadHeader(xSize,ySize,bpp,bscale,bzero)){
			Close();
			return FALSE;
		}
		if( (!x_matrix_size && !y_matrix_size) || gDoAutoAllocData ){
			// alloc buffer here :
			matrix.Alloc( xSize, ySize );			
		}else{
			Assert( xSize==x_matrix_size && ySize==y_matrix_size,"CCDMatrix have incorrect size (%d,%d) , FITS file has (%d,%d)",
					  x_matrix_size,y_matrix_size,xSize,ySize);						
		}
		
		// now read data :
		if(!ReadData( matrix.get_data_buffer(), xSize, ySize, bzero, bscale )){
			Close();
			return FALSE;
		}

		matrix.GetKeyTab() = m_HduList;
		Close();
		return TRUE;		
	}
	return FALSE;
}

template<class ARG_TYPE>
int CMyFITSFile<ARG_TYPE>::ReadFITSFileForAsasPipeline( ARG_TYPE* data, int sizeX, int sizeY,
                        							              const char* fname )
{
	mystring szInName=fname;
	szInName.env2str();

	int xSize=0,ySize=0,bpp=16;
	float bscale=1.00,bzero=0.00;

	int x_matrix_size = sizeX;
	int y_matrix_size = sizeY;


	if(Open( szInName.c_str() )){
		// first read header :
		if(!ReadHeader(xSize,ySize,bpp,bscale,bzero)){
			Close();
			return FALSE;
		}
		Assert( xSize==x_matrix_size && ySize==y_matrix_size,"CCDMatrix have incorrect size (%d,%d) , FITS file has (%d,%d)",
					x_matrix_size,y_matrix_size,xSize,ySize);						
		
		// now read data :
		if(!ReadData( data, xSize, ySize, bzero, bscale )){
			Close();
			return FALSE;
		}

		Close();
		return TRUE;		
	}
	return FALSE;
	
}


template<class ARG_TYPE>
int CMyFITSFile<ARG_TYPE>::ReadRaw( const char* fname, CSafeKeyTab& keyTab, char* data, int file_size )
{
	mystring szInName=fname;
   szInName.env2str();


	int xSize=0,ySize=0,bpp=16;
	float bscale=1.00,bzero=0.00;


	if(Open( szInName.c_str() )){
		// first read header :
		if(!ReadHeader(xSize,ySize,bpp,bscale,bzero)){
			Close();
			return -1;
		}
		keyTab = m_HduList;
		
		int ret = read( fd, data, file_size );
		close(fd);
		return ret;
	}
	return -1;	
}

template<class ARG_TYPE>
int CMyFITSFile<ARG_TYPE>::WriteRaw( const char* fname, CSafeKeyTab& keyTab, char* data, int file_size )
{
	if( file_size>0 && data ){
		mystring szOutName = fname;
   	szOutName.env2str();
                                                                                
	   if( !Open( szOutName.c_str(), TRUE ) )
			return -1;

		if(!AddHeaderKeys( keyTab, TRUE )){
			Close();
			return -1;
		}

		if(!AddCommentLines()){
			Close();
			printf("Error writing in AddCommentLines\n");
			return -1;
		}
		if(!WriteHeaderLine( "END" ,"" ,"" )){
			Close();
			printf("Error writing END to FITS header\n");
			return -1;
		}

		int ret=0;
		if(fd>0){
			write(fd, m_pHeader,m_HeaderSize);
	
			ret = write( fd, data, file_size );
			close(fd);
		}
	
		return ret;	
	}

	return -1;	
}
