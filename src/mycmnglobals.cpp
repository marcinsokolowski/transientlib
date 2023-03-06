#include "mycmnglobals.h"
#include "laplace_info.h"

int gDebugTracks=1;
int gPrintfLevel=4;
int gGlobalDebugLevel=0;
BOOL_T gDoDumpBadFit=FALSE;
BOOL_T gDoDumpAllHisto=FALSE;
BOOL_T gSafeLoadWarning=TRUE;
BOOL_T gChangeSeed=FALSE;
int gExternalSeed=-1;
BOOL_T gShowChi2Points3=FALSE;
BOOL_T gShowChi2_NPoints=FALSE;
BOOL_T gShowChi2OfAdded=FALSE;

class CCmnGlobalsCreator
{
public :
	CCmnGlobalsCreator();
	~CCmnGlobalsCreator(){};
};

CCmnGlobalsCreator::CCmnGlobalsCreator()
{
	InitGlobalLaplaceInfoTable();
}


CCmnGlobalsCreator gCmnGlobalsCreator;


int is_number(const char* string)
{
	if( string && string[0] ){
		int i=0;
		while( string[i] ){
			if( string[i]!='0' && string[i]!='1' && string[i]!='2' && 
			    string[i]!='3' && string[i]!='4' && string[i]!='5' &&
				 string[i]!='6' && string[i]!='7' && string[i]!='8' &&
			    string[i]!='9' ){
				return 0;
			}
			i++;
		}
		return 1;
	}
	return 0;
}

int is_digit( char znak )
{
	if( znak=='0' || znak=='1' || znak=='2' || znak=='3' || znak=='4' || 
		 znak=='5' || znak=='6' || znak=='7' || znak=='8' || znak=='9' ){
		return 1;
	}
	return 0;
}
