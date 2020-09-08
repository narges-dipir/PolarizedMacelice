#include <stdio.h>
#include <stdlib.h>	//srand
#include <time.h>	//time
#include <string.h>
#include "opMath.h"
#include "mat.h"
#include "rm.h"
#include "rcsv.h"
#include "msg.h"



opMath::opMath() {
}

opMath::opMath(const opMath& orig) {
}

opMath::~opMath() {
}

unsigned int combination(unsigned int n, const unsigned int& k)
{
	if (k > n) return 0;

	unsigned int result = 1;
	for (unsigned int i = 1; i <= k; i++)
	{
		result *= n--;
		result /= i;
	}

	return result;
}

unsigned int power(const unsigned int& exponent)
{
	return (1 << exponent);
}
