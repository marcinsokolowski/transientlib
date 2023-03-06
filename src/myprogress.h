#ifndef _MY_PROGRESS_H__
#define _MY_PROGRESS_H__



class CMyProgressBar
{
public :
	double m_min;
	double m_max;
	double m_curr;
	double m_prev;

	CMyProgressBar(double min,double max);
	void SetValue( double new_val );
	void Update();
	void Update( double new_val );
	
};


#endif
