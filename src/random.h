#ifndef _RANDOM_H__
#define _RANDOM_H__

#include "mytypes.h"

#define DEFAULT_SEED 477737
 
class CRandom
{
	private:
		static BOOL_T m_bInitialized;
		static unsigned long m_seed;
	public:
		static void Initialize(unsigned long seed=DEFAULT_SEED);
		CRandom();
		~CRandom();
		
		// returns random integer from range lowerband - upperband
		static ULONG_T GetRandomInteger(ULONG_T lowerbound,ULONG_T upperbound);
		static double GetRandom();
		
		static double GetGauss( double sigma, double mean );
		static double GetFastGauss( double sigma, double mean );			
};

#endif
