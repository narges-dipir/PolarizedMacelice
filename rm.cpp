#include <stdio.h>
#include <stdlib.h>	//srand
#include <time.h>	//time
#include <string.h>
#include "rm.h"
#include "mat.h"
#include "rcsv.h"
#include "opMath.h"
#include "msg.h"


void rmRUpdate(unsigned int& start, unsigned int& end, unsigned int **recovery, unsigned int **G, unsigned int **codeword, const unsigned int& n)
{
	for (int i = 0; i < n; i++)
	{
		unsigned int sum = 0;
		for (int j = start; j < end; j++)
		{
			sum += (recovery[0][j] * G[j][i]);
		}

		codeword[0][i] = (codeword[0][i] + sum) % 2;
	}
}

unsigned int rmDecision(unsigned int *array, unsigned int arrayLength)
{
	short zeroCount = 0;
	short oneCount = 0;

	for (int i = 0; i < arrayLength; i++)
	{
		if (0 == array[i])
		{
			zeroCount++;
		}
		else
		{
			oneCount++;
		}
	}

	if (zeroCount > oneCount)
		return 0;
	else if (zeroCount < oneCount)
		return 1;
	else
		return rand() % 2;
}


void rmIdxCalculator(unsigned int *array, const int& arrayLen, unsigned int& arrayPointer, unsigned int *result, unsigned int& resultIdx, unsigned int& resultValue)
{
	if (arrayLen == arrayPointer)
	{
		arrayPointer--;

		result[resultIdx] = resultValue;
		resultIdx++;

		return;
	}

	for (int i = 0; i < 2; i++)
	{
		unsigned short temp = power(array[arrayPointer]) * i;
		resultValue += temp;

		arrayPointer++;
		rmIdxCalculator(array, arrayLen, arrayPointer, result, resultIdx, resultValue);

		resultValue -= temp;
	}
	arrayPointer--;
}

void rmSCalculate(unsigned int *row, const int& rowLength, unsigned int *S)
{
	unsigned int *rowCopy = (unsigned int *)malloc(sizeof(unsigned int) * rowLength);

	for (int i = 0; i < rowLength; i++)
	{
		rowCopy[i] = row[i] - 1;
	}

	unsigned int SValue = 0;
	unsigned int SIdx = 0;

	unsigned int rowPointer = 0;

	rmIdxCalculator(rowCopy, rowLength, rowPointer, S, SIdx, SValue);

	free(rowCopy);

}

void rmECalculate(unsigned int *row, const int& rowLength, const int& m, unsigned int *E)
{

	int k = 0;
	for (int i = 0; i < m; i++)
	{
		int temp = 0;
		for (int j = 0; j < rowLength; j++)
		{
			if (i == (row[j] - 1))
			{
				temp = 1;
				break;
			}
		}

		if (temp == 0)
		{
			E[k] = i;
			k++;
		}
	}
}


void rmSCCalculate(unsigned int *row, const int& rowLength, unsigned int *SC)
{
	unsigned int SCValue = 0;
	unsigned int SCIdx = 0;

	unsigned int rowPointer = 0;

	rmIdxCalculator(row, rowLength, rowPointer, SC, SCIdx, SCValue);
}

unsigned int rmACalculate(unsigned int& SLength, unsigned int *S, const unsigned int& SCLength, unsigned int *SC, unsigned int **codeword)
{
	unsigned int *A = (unsigned int *)malloc(sizeof(unsigned int) * SCLength);

	unsigned int val;
	for (int i = 0; i < SCLength; i++)
	{
		unsigned int sum = 0;

		for (int j = 0; j < SLength; j++)
		{
			sum += codeword[0][S[j] + SC[i]];
		}
		A[i] = sum % 2;
	}

	val = rmDecision(A, SCLength);

	free(A);

	return val;
}



void rmDecoding(const unsigned int& n, const unsigned int& m, const unsigned int& rowLength, unsigned int* row, unsigned int rowIdx, unsigned int& rowValue, unsigned int **codeword, unsigned int **recovery, unsigned int& recoveryIdx)
{
	if (rowLength <= rowIdx)
	{
		/**/
		unsigned int SLength = power(rowLength);
		unsigned int *S = (unsigned int *)malloc(sizeof(unsigned int) * SLength);
		rmSCalculate(row, rowLength, S);
		unsigned int ELength = m - (rowLength);
		unsigned int *E = (unsigned int *)malloc(sizeof(unsigned int) * ELength);
		rmECalculate(row, rowLength, m, E);

		unsigned int SCLength = power(ELength);
		unsigned int *SC = (unsigned int *)malloc(sizeof(unsigned int) * SCLength);
		rmSCCalculate(E, ELength, SC);

		recovery[0][recoveryIdx] = rmACalculate(SLength, S, SCLength, SC, codeword);

		free(SC);
		free(E);
		free(S);

		recoveryIdx++;
		return;
	}

	for (unsigned int i = rowValue - 1; i >= 1; i--)
	{
		row[rowIdx] = i;
		rmDecoding(n, m, rowLength, row, rowIdx + 1, i, codeword, recovery, recoveryIdx);
	}

}


void rmDecoder(const unsigned int& r, const unsigned int& m, unsigned int **GRM, const unsigned int& n, unsigned int **codeword, unsigned int **recovery)
{
	unsigned int* multiplRows = (unsigned int*)malloc(sizeof(unsigned int) * r);

	for (unsigned int i = r; i >= 1; i--)
	{
		static unsigned int recoveryIdx = m + 1;
		static unsigned int upIdx = recoveryIdx;

		for (unsigned int j = m; j >= 1; j--)
		{
			int rowIdx = 0;
			if (i == 1)
			{
				recoveryIdx = j;
			}

			multiplRows[rowIdx] = j;
			rmDecoding(n, m, i, multiplRows, rowIdx + 1, j, codeword, recovery, recoveryIdx);
		}


		if (i == 1)
		{
			upIdx = 1;
			recoveryIdx = m + 1;
		}
		rmRUpdate(upIdx, recoveryIdx, recovery, GRM, codeword, n);
		upIdx = recoveryIdx;
	}

	recovery[0][0] = rmDecision(codeword[0], n);


	free(multiplRows);
}
