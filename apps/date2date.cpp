#include <stdio.h>
#include "mydate.h"
#include "mystring.h"
#include "myutil.h"
#include "AstroCCD.h"
#include "mystrtable.h"
#include "myparser.h"

double gRA=0.00;
double gDEC=0.00;

void print_usage()
{
	printf("date2date -night2unix=20050419 -jd2ux=2 -jd2ut -ux2ut -nonewline -localnight2unix=20050419 -lco2ux=20060831_211223 -frac2hjd=20070411.324 -ra=12.432 -dec=-12.00 -ut2ux=20070830_130220 -jd2hjd=1234.56789 -ut2doy=20070830_130220 -local2ux=20141208_131232 -ux2declochour=14343255325\n");
}

int main( int argc , char* argv[] )
{
	if(argc<2 || strncmp(argv[1],"-h",2)==0){
		print_usage();	
		exit(0);
	}

	mystring szNight2Unix,szUnix2UT,szJD,szJDUT,szUX2UT,szLocalNight2Unix,szLco2Ux,szJD2HJD,szFRACT,szUt2Ux,szUt2DOY,szLocal2Ux,szUs2DecimalLocalHour;
	BOOL_T bNoNewLine=FALSE;

	for(int i=1;i<argc;i++){
		int tmp_int;
		mystring szTmp;


		if( check_param( argv[i], "-localnight2unix", szLocalNight2Unix )){
			continue;
		}
		if( check_param( argv[i], "-night2unix", szNight2Unix )){
			continue;
		}
		if( check_param( argv[i], "-unix2ut", szUnix2UT )){
			continue;
		}
		if( check_param( argv[i], "-jd2ux", szJD ) ){
			continue;
		}
		if( check_param( argv[i], "-jd2ut", szJDUT ) ){
			continue;
		}
		if( check_param( argv[i], "-ux2ut", szUX2UT ) ){
			continue;
		}
		if( check_param( argv[i], "-nonewline", szTmp ) ){
			bNoNewLine = TRUE;
			continue;
		}
		if( check_param( argv[i], "-lco2ux", szLco2Ux ) ){
			continue;
		}
		if( check_param( argv[i], "-local2ux", szLocal2Ux ) ){
			continue;
		}
		if( check_param( argv[i], "-jd2hjd", szJD2HJD ) ){
			continue;
		}
		if( check_param( argv[i], "-frac2hjd",szFRACT ) ){
			continue;
		}
		if( check_param( argv[i], "-ra", gRA ) ){
			continue;
		}
		if( check_param( argv[i], "-dec", gDEC ) ){
			continue;
		}
		if( check_param( argv[i], "-ut2ux",szUt2Ux ) ){
			continue;
		}
		if( check_param( argv[i], "-ut2doy",szUt2DOY ) ){
			continue;
		}
		if( check_param( argv[i], "-ux2declochour", szUs2DecimalLocalHour ) ){
			continue;
		}
	}	

	if( strlen( szLocalNight2Unix.c_str() ) ){
		mystring szDTM;
		szDTM << szLocalNight2Unix.c_str() << "_000000";		
//		printf("szDTM = %s\n",szDTM.c_str());
		time_t dtm = get_unixtime_from_local_string( szDTM.c_str() );

		if( bNoNewLine ){
			printf("%d",dtm);
		}else{
			printf("%d\n",dtm);
		}
		exit(0);
	}


	if( strlen( szNight2Unix.c_str() ) ){
		mystring szDTM;
		szDTM << szNight2Unix.c_str() << "_000000";		
//		printf("szDTM = %s\n",szDTM.c_str());
		time_t dtm = get_gmtime_from_string( szDTM.c_str() );

		if( bNoNewLine ){
			printf("%d",dtm);
		}else{
			printf("%d\n",dtm);
		}
		exit(0);
	}

	if( strlen( szJD.c_str() ) ){
		time_t ux = AstroCCD::jd2ux( atof( szJD.c_str() ) );
		mystring szUT = get_gmtime_string( ux );
		if( bNoNewLine ){
			printf("%d",ux);
		}else{
			printf("%d\n",ux);
		}
		exit(0);
	}

	if( strlen( szJDUT.c_str() ) ){
		time_t ux = AstroCCD::jd2ux( atof( szJDUT.c_str() ) );
		mystring szUT = get_gmtime_string( ux );
		if( bNoNewLine ){
			printf("%s",szUT.c_str());
		}else{
			printf("%s\n",szUT.c_str());
		}
		exit(0);
	}

	if( strlen( szUX2UT.c_str() ) ){
		mystring szUT = get_gmtime_string( (int)atof( szUX2UT.c_str() ) );
		if( bNoNewLine ){
			printf("%s",szUT.c_str());
		}else{
			printf("%s\n",szUT.c_str());
		}
		exit(0);
	}
	if( strlen( szLco2Ux.c_str() ) ){
		CMyDate::SetCLT();
		time_t ux = get_unixtime_from_local_string( szLco2Ux.c_str() );
		printf("%d\n",ux);
	}
	if( strlen( szLocal2Ux.c_str() ) ){
		time_t ux = get_unixtime_from_local_string( szLocal2Ux.c_str() );
		printf("%d\n",ux);
	}
		
	if( strlen( szJD2HJD.c_str() ) ){
		CMyStrTable items;
		MyParser pars = szJD2HJD.c_str();
		pars.GetItems( items, "," );
		if( items.size()>=3 ){
			double jd = atof( items[0].c_str() );
			double ra = atof( items[1].c_str() );
			double dec = atof( items[2].c_str() );

			double hjd = AstroCCD::jd2hjd( jd, ra, dec );
			printf("hjd( jd=%.8f ) = %.8f\n",jd,hjd);
		}
	}
	if( strlen( szFRACT.c_str() ) ){
		double dday=atof( szFRACT.c_str() );
		int day=(int)atol( szFRACT.c_str() );
		double fract = dday-day;
		printf("fract = %.2f\n",fract);
		char tmp[128];
		sprintf(tmp,"%d_000000",day);
		time_t t = get_gmtime_from_string( tmp );
		printf("unix_time(%s) = %d\n",tmp,t);
		t = t + (int)(fract*(3600.00*24.00));
		printf("unix_time(%s) = %d\n",tmp,t);
		double hjd = AstroCCD::getHJD( t, gRA, gDEC );
		printf("hjd = %.8f\n",hjd);
	}

	if( strlen( szUt2Ux.c_str() ) ){
		time_t ux = get_gmtime_from_string( szUt2Ux.c_str() );
		printf("unixtime = %d\n",ux);
	}
	
	if( strlen( szUt2DOY.c_str() ) ){
		time_t ux = get_gmtime_from_string( szUt2DOY.c_str() );
		double  doy = get_doyd( ux );
		printf("doy = %.8f\n",doy);
	}
	
	if( strlen( szUs2DecimalLocalHour.c_str() ) ){
	   
	}
}
