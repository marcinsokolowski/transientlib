#include <stdio.h>
#include <stdlib.h>
#include <mydate.h>
#include <sat_info.h>
#include <Astroangle.h>
#include <AstroCCD.h>
#include <myutil.h>
#include <myfile.h>
#include <myprogress.h>

/*#include <ccd_globals.h>
#include <ccd_matrix.h>
#include <ccd_util.h>*/

extern "C" {
#include <predict.h>
}

#include <vector>

using namespace std;

mystring szTleFile="satelitesdb.tle";
mystring szQthFile="satelitesdb.qth";
int time_step=10;
int time_interval=0;
mystring gOutFile;
BOOL_T gDumpOnlyNonEmptyFiles=TRUE;
BOOL_T gPrintHeader=FALSE;
BOOL_T gPrintFullName=FALSE;
BOOL_T gFormatSatName=TRUE;
double gMinElevation=-1000;
double gMaxAltitude=1e20;

void print_usage()
{
   printf("usage:\n");
   printf("sattest ut_time=143544354 -all -tle=%s -qth=%s -ra=60[deg] -dec=0[deg] -ra_in_h=4 -radius=10 -interval=%d -time_interval=%d -time_step=%d -break_on_first -format_ra_in_deg -silent -show_closest -satname=PHOBOS -outfile=satlist.txt -shift_azim=0 -above_hor_only -print_header -print_full_name -show_geostat_only -fits=FITS.fits -outregfile=out.reg -mwa -min_elevation=%.4f -max_alt=%e\n",
               szTleFile.c_str(),szQthFile.c_str(),
               time_interval,time_interval,
               time_step,gMinElevation,gMaxAltitude);
   printf("INFO : options -interval and -time_range are the same !\n");
   printf("Format of the output file :\n");
   // szDTM_UT.c_str(),sat_ra_deg,sat_dec_deg,final_azim, i->sat_ele,i->sat_alt,ux_time,szSatName.c_str()
   printf("UT(Date-Time) RA[deg] DEC[deg] AZIM[deg] ELEVETION[deg] ALT[km] UNIX_TIME SAT_NAME X Y\n");
   printf("X Y - only if -fits option is provided\n");
   printf("-no_spaces : do not format satname with padding space to make it 30 characters wide\n");
}

BOOL_T gFitsRead=FALSE;
/*CCDMatrix gFitsFrame(0,0,FALSE,NOT_DEFINED, NULL, 0, FALSE, NULL, TRUE );

int ad2xy( const char* fitsfile, double ra_deg, double dec_deg, double& x, double& y )
{
   if( !gFitsRead ){
      gCCDParams.InitParams();
   
      printf("Reading fits file %s ...\n",fitsfile);
      if(!gFitsFrame.ReadFITSHeader( fitsfile )){
         printf("could not read header from file  : %s\n",fitsfile);
         exit(-1);
      }
      gFitsRead=TRUE;
      printf("Fits file %s read OK\n",fitsfile);fflush(stdout);
   }
      
   double ra_h = ra_deg/15.00;
   CCDUtil::ad2xy( gFitsFrame, ra_h, dec_deg, x, y, TRUE );                     
   
   return 1;
}*/

int main(int argc,char* argv[] )
{
   if(argc<2 || (argc>=2 && strncmp(argv[1],"-h",2)==0) )
   {
      print_usage();
      exit(0);
   }

   int idx=0;
   int bVisibleOnly=1;
   time_t ut_time=get_dttm();
   if(argc>=2)
      ut_time=atol(argv[1]);
   double ra=-10000;
   double dec=-10000;
   double radius=10.00;   
   int bBreakOnFirst=0;
   int bRA_In_DEG=0;
   int bVerb=1;
   int bShowClosest=0;
   double shift_azim=180.00;
   mystring gSatName,szTmp;
   BOOL_T bAboveHorOnly=FALSE;
   BOOL_T bGeoStatOnly=FALSE;
   mystring szFitsFile,szOutRegFile="sattest.reg";
   BOOL_T gAstroAzim=FALSE;

   for(int i=1;i<argc;i++){
      int tmp_int;
      double tmp_double;
      mystring szTmp;

      if( check_param( argv[i],"-min_elevation", gMinElevation ) ){
         continue;
      }
      if( check_param( argv[i],"-mwa", szTmp ) ){
         gAstroAzim=TRUE;
         continue;
      }
      if( check_param( argv[i],"-mwa", szFitsFile ) ){
         continue;
      }
      if( check_param( argv[i],"-outregfile", szOutRegFile ) ){
         continue;
      }
      if( check_param( argv[i],"-show_geostat_only", szTmp ) ){
         bGeoStatOnly=TRUE;
         continue;
      }
      if( check_param( argv[i],"-tle", szTleFile ) ){
         continue;
      }
      if( check_param( argv[i],"-qth", szQthFile ) ){
         continue;
      }
      if( check_param( argv[i],"-all", tmp_int ) ){
         bVisibleOnly=0;
         continue;
      }
      if( check_param( argv[i],"-radius", radius ) ){
         continue;
      }
      if( check_param( argv[i],"-ra_in_h", tmp_double ) ){
         ra = AstroAngle::hours2deg( tmp_double );
         continue;
      }
      if( check_param( argv[i],"-ra", ra ) ){
         continue;
      }
      if( check_param( argv[i],"-dec", dec ) ){
         continue;
      }
      if( check_param( argv[i],"-interval", time_interval ) ){
         continue;
      }
      if( check_param( argv[i],"-time_step", time_step ) ){
         continue;
      }
      if( check_param( argv[i],"-break_on_first",tmp_int ) ){
         bBreakOnFirst=1;
         continue;
      }
      if( check_param( argv[i],"-format_ra_in_deg",tmp_int ) ){
         bRA_In_DEG = TRUE;
         continue;
      }
      if( check_param( argv[i],"-silent",tmp_int ) ){
         bVerb = 0;
         continue;
      }
      if( check_param( argv[i],"-show_closest",tmp_int ) ){
         bShowClosest = 1;
         continue;
      }
      if( check_param( argv[i],"-satname", gSatName ) ){
         continue;
      }
      if( check_param( argv[i],"-outfile", gOutFile ) ){
         continue;
      }
      if( check_param( argv[i],"-shift_azim", shift_azim ) ){
         continue;
      }
      
      if( check_param( argv[i],"-above_hor_only", szTmp ) ){
         bAboveHorOnly=TRUE;
         continue;
      }                              

      if( check_param( argv[i],"-print_full_name", szTmp ) ){
         gPrintFullName = TRUE;
         continue;
      }                              

      if( check_param( argv[i],"-no_spaces", szTmp ) ){
         gFormatSatName = FALSE;
         continue;
      }                              

      if( check_param( argv[i],"-max_alt", gMaxAltitude ) ){
         continue;
      }                              

      if( check_param( argv[i],"-print_header", szTmp ) ){
                  gPrintHeader=TRUE;
             continue;
             }                              
   }

   CSatInfo::InitSateliteLibrary( szQthFile.c_str(), szTleFile.c_str() );      

   if( bVerb ){      
      printf("checking satelites at : %d\n",(int)ut_time);
      printf("using tle file        : %s\n",szTleFile.c_str() );
      printf("using qth file        : %s\n",szQthFile.c_str() );
      printf("only visible          : %d\n",bVisibleOnly );
      printf("(ra,dec)              : (%.2f,%.2f)\n",ra,dec);
      printf("radius                : %.2f [deg]\n",radius);
      printf("time_interval         : %d [sec]\n",time_interval);
      printf("time_step             : %d [sec]\n",time_step);
      printf("bBreakOnFirst         : %d\n",bBreakOnFirst);
      printf("bRA_In_DEG            : %d\n",bRA_In_DEG);
      printf("shift_azim            : %.2f [deg]\n",shift_azim);
      if( strlen(gSatName.c_str()) > 0 ){
         printf("satname               : %s\n",gSatName.c_str());
      }
      if( strlen(gOutFile.c_str()) > 0 ){
         printf("Out File Name         : %s\n",gOutFile.c_str());
      }
      printf("Astro azim            : %d\n",gAstroAzim);
      printf("Min elevation         : %.4f [deg]\n",gMinElevation);
      printf("Padding spaces        : %d\n",gFormatSatName);
      printf("Max altitude          : %.3f [m]",gMaxAltitude);
   }

   int ret=1;
//   int ret = CalcCurrentPosition( idx );
//   int ret = CalcCurrentPositionOfSingle( idx, ut_time );
//   CheckAll( ut_time );

   if( bVerb ){
      printf("NAME\tVISIBLILITY\tAZIM\tELEVATION\tALT\tDEC\tRA\tGEO\n");
      printf("-----------------------------------------------------\n");
   }
   vector<satInfo> satList;

   int found_count=0;
   int shown_count=0;


   time_t ux_time = ut_time;

   if( bVerb ){
      printf("SatName Visible Azim[deg] Ele[deg] Altitude[km] Dec RA GeoStat Dist[deg] SatCode HA[deg]\n");
      printf("-----------------------------------------------------\n");
   }

   double min_dist=10000000.00;
   satInfo* pClosest=NULL;

        if( gPrintHeader ){
           if( strlen(gOutFile.c_str()) ){
              MyOFile OutFile( gOutFile.c_str(), "w" );                           
              OutFile.Printf("# UT-DATE-TIME SAT_RA[deg] SAT_DEC[deg] AZIM[deg] ELE[deg] ALT[km] UXTIME SATTNAME HA[deg] GeoStat[flag] X Y\n");
           }           
        }
   
   CMyProgressBar* pProgress=NULL;
   if( time_interval > 3600 ){
      pProgress = new CMyProgressBar(0,time_interval);
   }

   while( ux_time <= (ut_time+time_interval) ){
      if(CSatInfo::GetSatInfo( satList, ux_time, bVisibleOnly )>0){
         vector<satInfo>::iterator i;
         for(i=satList.begin();i!=satList.end();i++){
            if(!bVisibleOnly || i->visibility==SAT_VISIBLE){
               mystring szSatName=i->sat_name;
               szSatName.replace_char(' ','_');
               if( gPrintFullName ){
                  szSatName += "_";
                  szSatName += i->sat_code;
               }
               
               if( i->sat_alt > gMaxAltitude ){
                  if( bVerb ){
                     printf("DEBUG : skipping %s due to altitude %.3f km > limit = %.3f [km]\n",szSatName.c_str(),i->sat_alt,gMaxAltitude);
                  }
                  continue;
               }


               double sat_ra_deg = AstroAngle::rad2deg( i->sat_ra );
               double sat_dec_deg = AstroAngle::rad2deg( i->sat_dec );
/*               if( ( ra<-1000 || fabs(ra-sat_ra_deg)<radius ) && ( dec<-1000 || fabs(dec-sat_dec_deg)<radius ) ){
                  double distance = 0;
                  if( ra>-1000 && dec>-1000 ){
                     distance = AstroCCD::CalcDistInDeg( ra, dec, sat_ra_deg, sat_dec_deg );
                  }*/
               double distance = 0;   
               if( ra>-1000 && dec>-1000 ){   
                  distance = AstroCCD::CalcDistInDeg( ra, dec, sat_ra_deg, sat_dec_deg );   
               }else{ 
                  if( ra >= -1000 ){ 
                     distance = fabs(ra-sat_ra_deg);
                  }
                  if( dec >= -1000 ){ 
                     double distance_dec = fabs(dec-sat_dec_deg);
                     if( distance_dec > distance ){
                        distance = distance_dec;
                     }
                  }
               }

               if( distance < radius ){
                  char szRA[64];
                  sprintf(szRA,"%s",AstroAngle::toString(i->sat_ra, ANGLE_RA_TYPE).c_str());
                  if( bRA_In_DEG ){
                      sprintf(szRA,"%.4f", AstroAngle::rad2deg( i->sat_ra ));
                  }
         
                  found_count++;

                  double final_azim= i->sat_azim;
                  final_azim = final_azim + shift_azim;
                  while( final_azim > 360 ){
                     final_azim = final_azim - 360;
                  }
                  while( final_azim < 0 ){
                     final_azim = final_azim + 360;
                  }
                  
                  if( gAstroAzim ){
                     double out_az = final_azim + 180.00;
                     if ( out_az > 360 ){
                         out_az = out_az - 360.00;
                     }
                     
                     final_azim = out_az;
                  }
                  
                  double ha_rad,ha_deg;
                  double geo_long_rad = obs_geodetic.lon;
                  AstroCCD::calculateHourAngle( i->sat_ra, i->sat_dec, ux_time, geo_long_rad , ha_rad);
                  ha_deg = AstroAngle::rad2deg(ha_rad);

                  if( !bShowClosest && (!bAboveHorOnly || i->sat_ele>0) && (!bGeoStatOnly || i->sat_geo_stat>0) && (i->sat_ele > gMinElevation) ){
                     if( strlen(gSatName.c_str())==0 || strstr(szSatName.c_str(),gSatName.c_str()) ){
                        printf("%s\t%c\t%.2f\t%.2f\t%.2f\t%.2f\t%s\t%d\t%.2f\t%s\t%.2f\n",szSatName.c_str(),
                           i->visibility, final_azim, i->sat_ele,
                           i->sat_alt, 
                           AstroAngle::rad2deg( i->sat_dec ),
//                           AstroAngle::toString(i->sat_dec, ANGLE_DEC_TYPE).c_str(),
                           szRA,
                           i->sat_geo_stat, distance, i->sat_code, ha_deg );
                           
                        double x=-1000,y=-1000;
//                        if( strlen(szFitsFile.c_str()) > 0 ){
//                           ad2xy( szFitsFile.c_str(), sat_ra_deg, sat_dec_deg, x, y );
//                        }
                           
                        shown_count++;      
                        if( strlen(gOutFile.c_str()) ){
                           MyOFile* pOutFile = new MyOFile( gOutFile.c_str(), "a" );                           
                           mystring szDTM_UT = get_gmtime_string( ux_time );
                           
                           char szOutName[128];
                           memset(szOutName,'\0',128);
                           sprintf(szOutName,"%s_%s",szSatName.c_str(),i->sat_code);

                           if( gFormatSatName ){
                              if( strlen(szOutName) < 30 ){
                                 int add_spaces = (30 - strlen(szOutName));
                                 mystring szSpaces;
                                 for(int s=0;s<add_spaces;s++){
                                    szSpaces += " ";
                                 }
                                 sprintf(szOutName,"%s%s",szSatName.c_str(),szSpaces.c_str());
                              }
                           }
                           pOutFile->Printf("%s %012.8f %012.8f %012.8f %012.8f %012.3f %d %s %012.8f %d %d %d\n",szDTM_UT.c_str(),sat_ra_deg,sat_dec_deg,final_azim, i->sat_ele,i->sat_alt,ux_time,szOutName,ha_deg,i->sat_geo_stat,(int)x,(int)y);
                           pOutFile->Close();
                           delete pOutFile;
                        }   
                     }
                  }

                  if( bBreakOnFirst ){
                     ux_time = (ut_time+time_interval)+1000;
                     break;
                  }
      
                  if( distance <  min_dist ){
                     min_dist = distance;
                     pClosest = &(*i);
                  }
               }
            }
         }
      }
      
      if( pProgress ){
         pProgress->Update( (int)(ux_time - ut_time) );
      }

      ux_time = ux_time + time_step;
   }

   if( bShowClosest && pClosest ){
      mystring szSatName=pClosest->sat_name;
      szSatName.replace_char(' ','_');
      if( gPrintFullName ){
         szSatName += "_";
         szSatName += pClosest->sat_code;
      }


      double sat_ra_deg = AstroAngle::rad2deg( pClosest->sat_ra );
      double sat_dec_deg = AstroAngle::rad2deg( pClosest->sat_dec );

      char szRA[64];
      sprintf(szRA,"%s",AstroAngle::toString(pClosest->sat_ra, ANGLE_RA_TYPE).c_str());
      if( bRA_In_DEG ){
          sprintf(szRA,"%.4f", AstroAngle::rad2deg( pClosest->sat_ra ));
      }

      double final_azim= pClosest->sat_azim;
      final_azim = final_azim + shift_azim;
      while( final_azim > 360 ){
         final_azim = final_azim - 360;
      }
      while( final_azim < 0 ){
         final_azim = final_azim + 360;
      }


      double x=-1000,y=-1000;
//      if( strlen(szFitsFile.c_str()) > 0 ){
//         ad2xy( szFitsFile.c_str(), sat_ra_deg, sat_dec_deg, x, y );
//      }

      printf("%s\t%c\t%.2f\t%.2f\t%.2f\t%.2f\t%s\t%d\t%.2f\t%s\t%d\t%d\n",szSatName.c_str(),
              pClosest->visibility, final_azim, pClosest->sat_ele,
              pClosest->sat_alt, 
            AstroAngle::rad2deg( pClosest->sat_dec ),
//            AstroAngle::toString(pClosest->sat_dec, ANGLE_DEC_TYPE).c_str(),
            szRA,
            pClosest->sat_geo_stat, min_dist, pClosest->sat_code, (int)x, (int)y );
   }

   if( bVerb ){
      printf("-----------------------------------------------------\n");
      printf("FOUND_SATELLITES_NUMBER = %d\n",found_count);
      printf("SHOWED = %d ( at UX-TIME = %d )\n",shown_count,(int)ux_time);
      printf("-----------------------------------------------------\n");
   }
      
      
   return ret;
}

