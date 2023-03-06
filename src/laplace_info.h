#ifndef _LAPLACE_INFO_H__
#define _LAPLACE_INFO_H__

#include "mytypes.h"
#include "mypoints.h"
#include "basedefines.h"


int InitGlobalLaplaceInfoTable();

class CLaplaceInfo
{
public:
	CLaplaceInfo();

	CLongPoint m_PlusList[MAX_PIXELS_IN_LAPLACE];
	CLongPoint m_MinusList[MAX_PIXELS_IN_LAPLACE];		
	LONG_T m_PlusCount;
	LONG_T m_MinusCount;
	double m_MinusFactor;
	double m_TypicalSigma;
		
};

extern CLaplaceInfo gLaplaceInfo[MAX_LAPLACE_DEFINED];
extern double gSigmaGCoef[MAX_LAPLACE_DEFINED];
extern double gSigmaBacgroundG[MAX_LAPLACE_DEFINED];
extern double gAverageBacgroundG[MAX_LAPLACE_DEFINED];

// sigma dependence on number of averaged frames :
struct sigmaDep
{
	double mean;
	double sigma;
	double mean_aver;
	double a;
	double b;
	double c;
};

extern sigmaDep gSigmaAverOfPrevNOMAX[MAX_LAPLACE_DEFINED];
extern sigmaDep gSigmaAverOfPrevMAX[MAX_LAPLACE_DEFINED];

double get_sigma( eLaplaceType_T lap, int nPrev, int bDoMax3x3, double sigmaNew );

double get_sigma_and_mean( eLaplaceType_T lap, int nPrev, int bDoMax3x3, 
						double sigmaNew, double& mean );

double get_gt_gauss_prob( double mean, double sigma, double val );
double get_lt_gauss_prob( double mean, double sigma, double val );
double get_significance( double lapNew, double lapPrev, double sigmaNew, double meanNew,
                         eLaplaceType_T lap, int nPrev, int bDoMax3x3, int bCoic );         

#endif
