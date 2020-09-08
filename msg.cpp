#include <stdio.h>
#include <stdlib.h>	//srand
#include <time.h>	//time
#include <string.h>
#include "msg.h"
#include "mat.h"
#include "rm.h"
#include "rcsv.h"
#include "opMath.h"


void msgEncrypt(const unsigned int& k, unsigned int** msg, unsigned int** publicKey, const unsigned int& n, unsigned int** error, unsigned int** prime, unsigned int** c)
{
	mtrixMultiplication(1, k, msg, k, n, publicKey, prime);
	errorAdd(n, prime[0], error[0], c[0]);
}
void msgDecrypt(const unsigned int& n, unsigned int** c, unsigned int** PInverse, const unsigned int& k, unsigned int** SInverse, unsigned int** msg, unsigned int** GRM, const unsigned int& r, const unsigned int& m)
{
	unsigned int** cPrime = (unsigned int**) malloc(sizeof(unsigned int*) * 1);
	cPrime[0] = (unsigned int*) malloc(sizeof(unsigned int) * n);

	mtrixMultiplication(1, n, c, n, n, PInverse, cPrime);
	unsigned int** mPrime = (unsigned int**)malloc(sizeof(unsigned int*) * 1);
	mPrime[0] = (unsigned int*)malloc(sizeof(unsigned int) * k);

	mtrixMultiplication(1, k, mPrime, k, k, SInverse, msg);
	free(cPrime[0]);
	free(cPrime);

	free(mPrime[0]);
	free(mPrime);
}
