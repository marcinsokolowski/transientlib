#include "mystring.h"
#include "myprogress.h"


		 
CMyProgressBar::CMyProgressBar(double min,double max)
: m_min(min),m_max(max),m_curr(m_min),m_prev(m_min)
{		
}

void CMyProgressBar::SetValue( double new_val )
{
	m_prev = m_curr;
	m_curr = new_val;
}
   
void CMyProgressBar::Update( double new_val )
{
	SetValue( new_val );
	Update();
}

void CMyProgressBar::Update()
{
	long prev_percent = (long)(((m_prev-m_min)/(m_max-m_min))*100);
	long curr_percent = (long)(((m_curr-m_min)/(m_max-m_min))*100);

	if(curr_percent!=prev_percent){
		mystring szProg;
		for(register int i=prev_percent+1;i<=curr_percent;i++){
			szProg << "/";
		}
		printf("%s",szProg.c_str());
		fflush(stdout);
	}
}
