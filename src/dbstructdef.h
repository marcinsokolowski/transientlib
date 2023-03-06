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
#ifndef _DBSTRUCT_DEF_H__
#define _DBSTRUCT_DEF_H__

#include <string>
#include <vector>
#include "pidbupdatedefs.h"
#include <ccd_common_defines.h>
#include <string.h>

using namespace std;

#define ASTRO_PARAMS 21
#define ASTRO_PARAMS_ACTIVE 15

struct cGRBLC
{
	cGRBLC() : glc_grb_id(0),glc_type(0),glc_time(0),glc_magnitude(0),glc_error(0),glc_limit(0),
				glc_ra(0),glc_dec(0),glc_trgnum(0),glc_obstype(0),lumdist(0),luminosity(0), flux(0)
		{
			glc_name[0] = '\0';
			glc_filter[0] = '\0';
			glc_sat[0] = '\0';
			strcpy( glc_system , "JOHNSON" );
		};
		
	char glc_name[14];			
	int glc_grb_id;
	int glc_type;
	double glc_time;
	double glc_magnitude;
	double glc_error;
	char glc_filter[20];
	char glc_system[64];
	int glc_limit;
	double glc_ra;
	double glc_dec;
	char glc_sat[20];
	int glc_trgnum;
	int glc_obstype;


	// for calculations 
	double lumdist; // luminocity distance 
	double luminosity; // luminosity
	double flux; //  
};

class CGRBInfo {
public:
	CGRBInfo();
	CGRBInfo( const CGRBInfo& right );		  	
	CGRBInfo& operator=( const CGRBInfo& right );

	double ra;
	double dec;
	time_t grbTime;
	int source_id;
	int gcn_id; // nie wiem jeszcze czym sie roznia wiec moze
	int grb_id; // na razie oba
	char foundgrbname[20];	                                                
	char z[20];	
	bool has_ot;

	// galactic extinction correction :
	double galext_u;
	double galext_b;
	double galext_v;
	double galext_r;
	double galext_i;
	double galext_j;
	double galext_h;
	double galext_k;
	double galext_lprim;

	
	// ot light curve :
	vector<cGRBLC> m_OtLc;
};

#define EVT_SLT_REASON_SIZE 127
struct CEvent {
	CEvent():
		frame_no(0),iFrameSecondNo(0),iDayNight(0),iCamId(0),event_no(0),
		x(0),y(0),lap_new(0),lap_prev(0),lap_max_new(0),s_raw(0),
		sphericity(0),siginificance(0),cluster_size(0),coic_radius_pix(0),
		aver_in_pixel(0),black_ratio(0),
		coic_radius_sec(0),ra(0),dec(0),evt_time(0),mag(0),evt_slt_ok(false),
		evt_runtype(0),evt_online(1),lap_next(0),evt_type(eFinal),
		evt_min_dist_cone(0)
	{}

	int 	frame_no;	//Frame no [uzywane do odnalezienia if_frm]
	int 	iFrameSecondNo;	//Second Frame number [uzywane do odnalezienia if_frm, nie zapisywane]
	int 	iDayNight;		//RUN id [uzywane do odnalezienia if_frm, nie zapisywane]
	int	iCamId;		//Camera id	[uzywane do odnalezienia if_frm, nie zapisywane]
	int 	event_no;	//Event no on given frame
	int 	x;		//x position on chip
	int 	y;		//y position on chip
	double	lap_new;
	double   lap_next;
	double	lap_prev;
	double	lap_max_new;
	double	s_raw;
	double	sphericity;
	double	siginificance;
	double	cluster_size;
	double	coic_radius_pix;
	double	aver_in_pixel;
	double	black_ratio;
	string	event_type;
	double	coic_radius_sec;
	double	ra;
	double	dec;
	time_t	evt_time;
	string	evtPath;
	string	evtList;
	double	mag;
	string evtid;
	bool     evt_slt_ok;
	int      evt_online;
	int      evt_runtype;
	int      evt_type;
	string   evt_slt_reason;
	double   evt_min_dist_cone;
};

struct CFrameDetails {
	
	CFrameDetails():
		iCAMIIDX(0),fEXPTIME(0),fREXPTIME(0),bSHUTTER(0),
		fADCGAIN(0),fADCBIAS(0),fADCGSET(0),fLNAGAIN(0),fADCBSET(0),
		fADCRANGE(0),fADCCLAMP(0),fELECGAIN(0),bCOOLING(0),fABINN(0),
		fSBINN(0),fSPEED(0),bMPP_BC(0),fRO_TIME(0),fCROTIME(0),iFOCUS(0),bHITLENS(0),
		fRNOISE(0),fRELNOISE(0),
		fCHIPTSET(0),fCHIPTEMP(0),fCASTEMP(0),fAMBTEMP(0),fCAMHUMID(0),fAMBHUMID(0),
		fINTRTEMP(0),fAIRMASS(0)
		{}
				

	// ---------- Exposure id ---------- //
	string sOBSERVER;//eg.'PiMan' / OBSERVER NAME/IDENTIFICATION
	string sSOFTWARE;//eg.'PISOFT v 1.00 beta' / Software Version
	string sBUILD;//eg.'Sun Jun 6 22:36:21 CEST 2004 '/ Software Build
	string sDRVTYPE;//eg.'2K2K Camera Driver' / CCD Driver type
	int iCAMIIDX;//eg.1 / Camera internal index
	// ---------- Exposure settings ---------- //
	double fEXPTIME;//eg.5.000000 / TRUE EXPOSURE TIME IN SECONDS
	double fREXPTIME;//eg.5.209000 / Measured exposition time
	bool   bSHUTTER;//true='OPEN'  false='CLOSED' (FOR DARKS)
	double fADCGAIN;//eg. 1.000000 / ADC gain value
	double fADCBIAS;//eg. 0.000000 / ADC offset [mV]
	double fADCGSET;//eg. 0.000000 / ADC gain setting
	double fLNAGAIN;//eg. 8.000000 / LNA gain - preamplifier gain
	double fADCBSET; //eg. 0.000000 / ADC offset setting
	double fADCRANGE;//eg.4 / ADC range [V]
	double fADCCLAMP;//eg.4V / ADC 4V clamping ON/OFF
	double fELECGAIN;//eg.4.359654 / ELECTRONS PER ADU
	bool bCOOLING;//'ON'=true OR 'OFF'=falas
	double fABINN;//eg.'1' / ANALOG BINNING
	double fSBINN;//eg.'1' / SOFTWARE BINNING
	double fSPEED; //eg.0 / READOUT SPEED
	string sSPEEDMH;//eg.'500kHz / 2 MHz' / READOUT SPEED MHz Vertical/Horizontal
	bool bMPP_BC; //'MPP'=ture OR 'BC'=false
	double fRO_TIME;//eg. 0.802000 / USB transfer time [sec]
	double fCROTIME;//eg.2.570000 / Chip read time [sec]
	int iFOCUS;//eg.100 / POSITION OF FOCUS MOTOR IN STEPS
	bool bHITLENS;//eg.'OFF' / Lens hitting on/off
	string sSAVEAREA;//eg.'0 0 2061 2047' / X0,Y0,X1,Y1
	string sUSBMODE;//eg.'2.0' / USB mode
	string sFPGAVER;//eg.'2004-05-04 / 4b' / FPGA soft version
	string sCPRSVER;//eg.'2004-04-1a / 01' / CPRS soft version
	string sVERDESC;//eg.'Kasprowicz Stankiewicz'/ Camera software version desc
	double fRNOISE;//eg.52.431057 / Readout noise ( +gauss , -rms )
	double fRELNOISE;//eg.228.581268 / Readout noise in electrons
	//---------- Exposure environment ----------//
	double fCHIPTSET;//eg.-10.000000 / Chip temperature parameter
	double fCHIPTEMP;//eg.-10.000000 / CCD TEMPERATURE *C
	double fCASTEMP;//eg.27.000000 / CAMERA CASE TEMPERATURE *C
	double fAMBTEMP;//eg.71.000000 / AMBIENT TEMPERATURE *C
	double fCAMHUMID;//eg.0.000000 / CAMERA HUMIDITY %
	double fAMBHUMID;//eg.0.000000 / AMBIENT HUMIDITY %
	double fINTRTEMP;//eg.0.000000 / Internal temperature ( chamber )
	double fAIRMASS;//eg.1.00000000 / air mass @ end exposure
};


struct CFrame {
	CFrame():
		iFrameNo(0),frmDet(NULL),istarcount(0),ifitssize(0),
		posangle(0),ast_ord(0),asttime(0),ast_err(0),
		ra2000(0),dec2000(0),flip(0),matchedstars(0),errcode(0),
		shutter_mode(2),fwhm(0),avg(0),rms(0),star_fraction(0),cat_stars(0),
		astrook(0),avg_shape(0),pixscale(0)
	{
		for(int i=0;i<ASTRO_PARAMS;i++){
			px[i] = 0;
			py[i] = 0;
		}
	}			
	
	int iFrameNo;		//Frame number in given RUN
	CFrameDetails *frmDet;	//details given Frame
	int iFrameSecondNo;	//Second Frame number [not used]
	int iDayNight;		//RUN id
	string sPathToFile;	//Path to File
	int iCamId;		//Camera id
	string sCamFilter;	//Camera filter
	int iNaxis1;		// length of data axis 1
	int iNaxis2;		// length of data axis 2
	//---------- Object ----------//
	string sOBJECT;		//eg.'SKY' / Object
	double fROTATE;		//eg.1 / S is UP ( rotated FOV )
	double fRA;		//eg.20.08735108 / RightAscension - observed
	double fHA;		//eg.23.76398872 / Hour angle
	double fDEC;		//eg.43.48200639 / Declination - observed
	double fALT;		//eg.81.01740673 / Altitude - observed
	double fZENITH_D;	//eg.8.98259250 / zenith distance @ end exposure
	double fAZIM;		//eg.343.32409065 / Azimuth - observed
	//---------- Exposure date/time ----------//
	string sDATE_OBS;	//eg.'07/06/04' / GMT
	time_t tTIME_UT;	//eg.1086571503 / time_t - sec since 0 UTC 1/1/1970 at start expo
	string sLOCTIME;	//eg.'03:25:03.199' / Local time @ frame start
	string sLOCDATE;	//eg.'2004-06-07' / Local date @ frame start
	double fJD;		//eg.2453163.55906250 / Julian Date ad mid_exposure
	double fHJD;		//eg.2453163.55906250 / Heliocentric JD at mid_exposure
	
	int istarcount;
	int ifitssize;
	int matchedstars;
	int errcode;
	
	double posangle;
	int ast_ord;
	time_t asttime;
	string ast_ver;
	double ast_err;
	double ra2000;
	double dec2000;
	double px[21];
	double py[21];	  
	int flip;
	int shutter_mode;
	double fwhm;
	double avg;
	double rms;
	double star_fraction;
	int cat_stars;
	int astrook;
	double avg_shape;	
	double pixscale;
	
	vector<string> aver_tab;
};


struct cSN
{
	string sn_name;
	string sn_host_glx;
	string sn_date;
	double sn_ra;
	double sn_dec;
	string sn_offset;
	double sn_mag;
	string sn_type;
	string sn_discoverer;
};


struct cAlert
{
	cAlert()
		: alert_id(0),source_id(0),alert_type(0),trg_num(0),seq_num(0),instrument_id(0),
		  alert_status(0),datetjd(0),time(0),flax(0),duration(0),gamma_rate(0),
		  coordra(0),coorddec(0),coord1(0),coord2(0),grb_id(0),validity(eValidityUnknown),
		  sigma(0),unix_time(0), sod(0), coorderr(0), gal_long(0), gal_lat(0),
		  ecl_long(0),ecl_lat(0),flux1(0),flux2(0),flux3(0),flux4(0),flux(0)
	{}		  	

	 int alert_id;
    int source_id;
    int alert_type;
    int trg_num;
    int seq_num;
    int instrument_id;
    int alert_status;
    int datetjd;
    string date;
    double time;
    string rays;
    string band;
    int flax;
    int duration;
    int gamma_rate;
    double coordra;
    double coorddec;
    double coorderr;
    double coord1;
    double coord2;
    string fname;
    int grb_id;
	 eAlertValidity validity;	
	 double sigma;
	 double unix_time;

	// additional - NOT DB fields:
	double sod; //seconds of day 
	
	double gal_long;
	double gal_lat;
	
	double ecl_long;
	double ecl_lat;

	double flux1;
	double flux2;
	double flux3;
	double flux4;
	double flux;
	
	string comment;
};


struct cDBStarCat
{
	cDBStarCat()
	: id_star(0),ra(0),dec(0),mag(0),mag_cat(0),asas_star_id(0),camid(0)
	{
		name[0]='\0';
	}

	int id_star;
	double ra;
	double dec;
	double mag;
	double mag_cat;
	char name[24];
	int asas_star_id;
	int camid;
};

struct cDBStarObs
{
	cDBStarObs()
	: id_star(0),id_frm(0),ra(0),dec(0),mag(0),mag_piphoto(0),time_hjd(0),
	  x(0),y(0)
	{}

	int id_star;
	int id_frm;
	double ra;
	double dec;
	double mag;	
	double mag_piphoto;
	double time_hjd;
	double x;
	double y;
};

struct cDBTrack
{
	cDBTrack()
		: tr_id(0),tr_night(0),tr_a(0),tr_b(0),tr_frame_start(0),
		  tr_frame_end(0),tr_vx(0),tr_vy(0),tr_type(0)
		{}
	int tr_id;
	int tr_night;
	double tr_a;
	double tr_b;
	int tr_frame_start;
	int tr_frame_end;
	double tr_vx;
	double tr_vy;
	int tr_type;			
	double tr_radec_a;
	double tr_radec_b;
	double tr_v_ra;
	double tr_v_dec;
};

struct cDBEventOnTrack
{
	cDBEventOnTrack()
		: eotr_tr_id(0),eotr_tr_night(0),eotr_frame(0),eotr_x(0),eotr_y(0)
	{}		

	int eotr_tr_id;
	int eotr_tr_night;
	int eotr_frame;
	double eotr_x;
	double eotr_y;
	double eotr_ra;
	double eotr_dec;		
};

struct cDBPointingInfo 
{
	cDBPointingInfo() 
		: pt_sat(0),pt_start_time(0),pt_end_time(0),pt_target_no(0),pt_ra(0),
		  pt_dec(0)
	{
		pt_object[0] = '\0';
	}
	
	int pt_sat;
	int pt_start_time;
	int pt_end_time;
	int pt_target_no;
	char pt_object[24];
	double pt_ra;
	double pt_dec;	
};

struct cTargetToObs
{
	cTargetToObs()
		: to_ra(0),to_dec(0),to_obshour(0),to_obstime(0),to_mountid(0),to_active(0),to_creatdtm(0)
	{
		to_name[0] = '\0';
		to_obssite[0] = '\0';
		to_object[0] = '\0';
		to_desc[0] = '\0';
	}
		
	char to_name[128];
	double to_ra;
	double to_dec;
	int to_obshour;
	char to_obssite[16];
	int to_obstime;
	char to_object[16];
	char to_desc[128];
	int to_mountid;
	int to_active;
	int to_creatdtm;
};

struct cMountDef 
{
	cMountDef() : 
		md_mount_id(-1),md_shift_ra_eq(0),md_shift_dec_eq(0),md_shift_radius_deg(0)
	{
		md_comment[0] = '\0';
	}

	int md_mount_id; // MOUNT ID
   int md_shift_ra_eq; // RA shift from main target at equator in degrees
   int md_shift_dec_eq; // DEC shift from main target at equator in degrees
   int md_shift_radius_deg; // radius of shift from main target in degrees
   char md_comment[1024]; // comment 
};



#endif
