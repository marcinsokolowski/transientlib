#ifndef _MYMACROS_H__
#define _MYMACROS_H__

#include "mydate.h"
#include "mycmnglobals.h"

#define tmp_printf printf
#define printf_now( msg ) printf( msg );fflush(stdout);
#define printf_now2( fmt, msg ) printf( fmt, msg );fflush(stdout);
#define printf_now3( fmt, msg1, msg2 ) printf( fmt, msg1, msg2 );fflush(stdout);
#define printf_now4( fmt, msg1, msg2, mag3 ) printf( fmt, msg1, msg2, mag3 );fflush(stdout);

#define _TRACE_PRINTF(level) if(gPrintfLevel>level)printf


#define _STDOUT_TRACE_1( msg ) { printf("[ %s ]\t",CMyDate::getdate());printf( msg );fflush(stdout); }
#define _STDOUT_TRACE_ _STDOUT_TRACE_1
#define _STDOUT_TRACE_2( msg1, msg2 ) { printf("[ %s ]\t",CMyDate::getdate());printf( msg1, msg2 );fflush(stdout); }
#define _STDOUT_TRACE_3( msg1, msg2, msg3 ) { printf("[ %s ]\t",CMyDate::getdate());printf( msg1, msg2, msg3 );fflush(stdout); }

#define _STDOUT_TRACE_NL_1( msg ) { printf("\n[ %s ]\t",CMyDate::getdate());printf( msg );fflush(stdout); }
#define _STDOUT_TRACE_NL _STDOUT_TRACE_1
#define _STDOUT_TRACE_NL_2( msg1, msg2 ) { printf("\n[ %s ]\t",CMyDate::getdate());printf( msg1, msg2 );fflush(stdout); }
#define _STDOUT_TRACE_NL_3( msg1, msg2, msg3 ) { printf("\n[ %s ]\t",CMyDate::getdate());printf( msg1, msg2, msg3 );fflush(stdout); }


#define _STDOUT_TRACE_1_LEVEL_0( msg ) if(gPrintfLevel>=0){ printf("[ %s ]\t",CMyDate::getdate());printf( msg );fflush(stdout); }
#define _STDOUT_TRACE_LEVEL_0 _STDOUT_TRACE_1_LEVEL_0
#define _STDOUT_TRACE_2_LEVEL_0( msg1, msg2 ) if(gPrintfLevel>=0){ printf("[ %s ]\t",CMyDate::getdate());printf( msg1, msg2 );fflush(stdout); }
#define _STDOUT_TRACE_3_LEVEL_0( msg1, msg2, msg3 ) if(gPrintfLevel>=0){ printf("[ %s ]\t",CMyDate::getdate());printf( msg1, msg2, msg3 );fflush(stdout); }

#define _STDOUT_TRACE_1_LEVEL_1( msg ) if(gPrintfLevel>=1){ printf("[ %s ]\t",CMyDate::getdate());printf( msg );fflush(stdout); }
#define _STDOUT_TRACE_LEVEL_1 _STDOUT_TRACE_1_LEVEL_1
#define _STDOUT_TRACE_2_LEVEL_1( msg1, msg2 ) if(gPrintfLevel>=1){ printf("[ %s ]\t",CMyDate::getdate());printf( msg1, msg2 );fflush(stdout); }
#define _STDOUT_TRACE_3_LEVEL_1( msg1, msg2, msg3 ) if(gPrintfLevel>=1){ printf("[ %s ]\t",CMyDate::getdate());printf( msg1, msg2, msg3 );fflush(stdout); }

#define _STDOUT_TRACE_1_LEVEL_2( msg ) if(gPrintfLevel>=2){ printf("[ %s ]\t",CMyDate::getdate());printf( msg );fflush(stdout); }
#define _STDOUT_TRACE_LEVEL_2 _STDOUT_TRACE_1_LEVEL_2
#define _STDOUT_TRACE_2_LEVEL_2( msg1, msg2 ) if(gPrintfLevel>=2){ printf("[ %s ]\t",CMyDate::getdate());printf( msg1, msg2 );fflush(stdout); }
#define _STDOUT_TRACE_3_LEVEL_2( msg1, msg2, msg3 ) if(gPrintfLevel>=2){ printf("[ %s ]\t",CMyDate::getdate());printf( msg1, msg2, msg3 );fflush(stdout); }

#define _TRACE_PRINTF_  if(gPrintfLevel>=-1)printf
#define _TRACE_PRINTF_0 if(gPrintfLevel>=0)printf
#define _TRACE_PRINTF_1 if(gPrintfLevel>=1)printf
#define _TRACE_PRINTF_2 if(gPrintfLevel>=2)printf
#define _TRACE_PRINTF_3 if(gPrintfLevel>=3)printf
#define _TRACE_PRINTF_4 if(gPrintfLevel>=4)printf
#define _TRACE_PRINTF_5 if(gPrintfLevel>=5)printf
#define _TRACE_PRINTF_6 if(gPrintfLevel>=6)printf
#define _TRACE_PRINTF_7 if(gPrintfLevel>=7)printf
#define _TRACE_PRINTF_8 if(gPrintfLevel>=8)printf


#ifdef _ENABLE_PROFILER_
#define PROFILER_START clock_t t1=clock();
#define PROFILER_RESTART t1=clock();
#define PROFILER_STOP clock_t t2=clock();
#define PROFILER_STOP2 t2=clock();
#define PROFILER_END(desc) clock_t t2=clock(); \
					      mystring msg=get_clock_in_sec_string( t2-t1 );\
							MYTRACE2(gCCDTrace,"PROFILER :" << desc << msg );\
							_TRACE_PRINTF_0("PROFILER %s %s\n",desc,msg.c_str());
#define PROFILER_END_CMN(desc) clock_t t2=clock(); \
					      mystring msg=get_clock_in_sec_string( t2-t1 );\
							_TRACE_PRINTF_0("PROFILER %s %s\n",desc,msg.c_str());
#define PROFILER_REEND(desc) t2=clock(); \
										msg=get_clock_in_sec_string( t2-t1 );\
										MYTRACE2(gCCDTrace,desc << msg );\
										_TRACE_PRINTF_3("PROFILER %s %s\n",desc,msg.c_str());													
#define PROFILER_PUT(desc,total) mystring _msg=get_clock_in_sec_string( total );\
											MYTRACE2(gCCDTrace,"PROFILER :" << desc << _msg );\
											_TRACE_PRINTF_3("PROFILER %s %s\n",desc,_msg.c_str());


#define START_TIMER 

											
#else
#define PROFILER_START
#define PROFILER_RESTART
#define PROFILER_END(desc)
#define PROFILER_REEND(desc)
#define PROFILER_PUT
#endif


#endif
