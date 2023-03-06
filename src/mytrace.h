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
// mytrace.h: interface for the CMyTrace class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MYTRACE_H__5212D424_3F44_11D6_B636_B42D07C10000__INCLUDED_)
#define AFX_MYTRACE_H__5212D424_3F44_11D6_B636_B42D07C10000__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "myfile.h"
#include "mydate.h"

class BASELIB_EI CMyTrace : public MyOFile  
{
protected:
	long m_TrLevel;
public:
	CMyTrace(const char* filename,long TrLevel);
	virtual ~CMyTrace();
	long GetTrLevel(){ return m_TrLevel;}
	void SetTrLevel( int TrLevel ){ m_TrLevel = TrLevel; }
};

#ifdef _ENABLE_TRACING_

#define MYTRACE(TrObj,Cmd)  { TrObj << "[" << CMyDate::getdate() << "] " << Cmd <<"\n"; TrObj.Flush(); TrObj.Close();}
#define MYTRACE0(TrObj,Cmd) if(TrObj.GetTrLevel()>=0){ TrObj << "[" << CMyDate::getdate() << "] " << Cmd <<"\n"; TrObj.Flush(); TrObj.Close();}
#define MYTRACE1(TrObj,Cmd) if(TrObj.GetTrLevel()>=1){ TrObj << "[" << CMyDate::getdate() << "] " << Cmd <<"\n"; TrObj.Flush(); TrObj.Close();}
#define MYTRACE2(TrObj,Cmd) if(TrObj.GetTrLevel()>=2){ TrObj << "[" << CMyDate::getdate() << "] " << Cmd <<"\n"; TrObj.Flush(); TrObj.Close();}
#define MYTRACE3(TrObj,Cmd) if(TrObj.GetTrLevel()>=3){ TrObj << "[" << CMyDate::getdate() << "] " << Cmd <<"\n"; TrObj.Flush(); TrObj.Close();}
#define MYTRACE4(TrObj,Cmd) if(TrObj.GetTrLevel()>=4){ TrObj << "[" << CMyDate::getdate() << "] " << Cmd <<"\n"; TrObj.Flush(); TrObj.Close();}
#define MYTRACE5(TrObj,Cmd) if(TrObj.GetTrLevel()>=5){ TrObj << "[" << CMyDate::getdate() << "] " << Cmd <<"\n"; TrObj.Flush(); TrObj.Close();}
#else
#define MYTRACE0(TrObj,Cmd)
#define MYTRACE1(TrObj,Cmd)
#define MYTRACE2(TrObj,Cmd)
#define MYTRACE3(TrObj,Cmd)
#define MYTRACE4(TrObj,Cmd)
#define MYTRACE5(TrObj,Cmd)
#endif
#define ALWAYS_TRACE(TrObj,Cmd) { TrObj << "[" << CMyDate::getdate() << "] " << Cmd <<"\n"; TrObj.Flush(); TrObj.Close();}


class CCmnTrace : public CMyTrace
{
public :
   CCmnTrace();
};
      
extern CCmnTrace gCmnTrace;
   

#endif // !defined(AFX_MYTRACE_H__5212D424_3F44_11D6_B636_B42D07C10000__INCLUDED_)
