#include "mat.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>	//time
#include <string.h>
#include "rm.h"
#include "rcsv.h"
#include "opMath.h"
#include "msg.h"



mat::mat() {
}

mat::mat(const mat& orig) {
}

mat::~mat() {
}

int inverseMatrix(unsigned int** matrix, unsigned int** inverse, const unsigned int& size)
{
	unsigned int** clone = (unsigned int**)malloc(sizeof(unsigned int*) * size);

	for (int i = 0; i < size; i++)
	{
		clone[i] = (unsigned int*)malloc(sizeof(unsigned int) * size);
		for (int j = 0; j < size; j++)
		{
			clone[i][j] = matrix[i][j];
		}
	}

	for (int i = 0; i < size; i++)
	{
		for (int j = i + 1; j < size; j++)
		{
			if (clone[i][i] == 0 && clone[j][i] != 0)
			{
				unsigned int* swap = clone[i];
				clone[i] = clone[j];
				clone[j] = swap;

				swap = inverse[i];
				inverse[i] = inverse[j];
				inverse[j] = swap;
			}
			else if (clone[j][i] != 0)
			{
				for (int k = 0; k < size; k++)
				{
					clone[j][k] = clone[j][k] ^ clone[i][k];
					inverse[j][k] = inverse[j][k] ^ inverse[i][k];
				}
			}
		}

		if (clone[i][i] == 0)
		{
			return 0;
		}

		for (int j = 0; j < i; j++)
		{
			if (clone[j][i] != 0)
			{
				for (int k = 0; k<size; k++)
				{
					clone[j][k] = clone[j][k] ^ clone[i][k];
					inverse[j][k] = inverse[j][k] ^ inverse[i][k];
				}
			}
		}
	}


	/*printf("inverse matrix\n");
	for (int i = 0; i < size; i++)
	{
	for (int j = 0; j < size; j++)
	{
	printf("%d ", inverse[i][j]);
	}
	printf("\n");
	}
	printf("\n"); printf("\n");

	printf("verify inverse matrix\n");
	for (int i = 0; i < size; i++)
	{
	for (int j = 0; j < size; j++)
	{
	int sum; sum = 0;
	for (int k = 0; k < size; k++)
	{
	sum ^= matrix[i][k] * inverse[k][j];
	}
	printf("%d ", (sum));

	}
	printf("\n");
	}
	printf("\n"); printf("\n");*/


	for (int i = 0; i < size; i++)
	{
		free(clone[i]);
	}
	free(clone);

	return 1;
}

void nonSingularMatrix(unsigned int** matrix, unsigned int** identity, const unsigned int& size)
{
	srand((unsigned)time(NULL));

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			matrix[i][j] = rand() % 2;

			identity[i][j] = 0;
			if (i == j)
			{
				identity[i][j] = 1;
			}
		}
	}
}

void permutationMatrix(unsigned int** matrix, unsigned int** identity, const unsigned int& size)
{
	srand((unsigned)time(NULL));

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			matrix[i][j] = 0;
			identity[i][j] = 0;
			if (i == j)
			{
				matrix[i][j] = 1;
				identity[i][j] = 1;
			}
		}

		int idx = rand() % (i + 1);
		if (i != idx)
		{
			unsigned int* swap = matrix[i];
			matrix[i] = matrix[idx];
			matrix[idx] = swap;
		}
	}
}

void mtrixMultiplication(const unsigned int& pRow, const unsigned int& pCol, unsigned int** prev, const unsigned int& nRow, const unsigned int& nCol, unsigned int** next, unsigned int** result)
{
    
	for (int i = 0; i < pRow; i++)
	{
		for (int j = 0; j < nCol; j++)
		{
            int sum=0;
			for (int z = 0; z < pCol; z++)
			{
                sum = (int)(prev[i][z] * next[z][j]);
			}
			result[i][j] = sum;
			//printf("%d ", result[i][j]);
		}
		//printf("\n");
	}

}

void errorAdd(const unsigned int& size, unsigned int* original, unsigned int* error, unsigned int* result)
{
	for (int i = 0; i < size; i++)
	{
		result[i] = original[i];

		if (error[i] == 1)
		{
			if (original[i] == 0)
			{
				result[i] = 1;
			}
			else
			{
				result[i] = 0;
			}
			//result[i] = (~original[i])%2;
		}
	}

}

void errorGenerator(const unsigned int& weight, const unsigned int& size, unsigned int* error)
{
	srand((unsigned)time(NULL));

	int num = 0;

	while (1)
	{
		unsigned int idx = rand() % size;

		if (error[idx] == 0)
		{
			error[idx] = 1;
			num++;
		}

		if (num >= weight)
		{
			break;
		}
	}
}

void generatePublicKey(const unsigned int& k, const unsigned int& n, unsigned int** S,unsigned int** G, unsigned int** P, unsigned int** key,const unsigned int& t, unsigned int* z)
{
	unsigned int** pub = (unsigned int**)malloc(sizeof(unsigned int*) * k);

	for (int i = 0; i < k; i++)
	{
		pub[i] = (unsigned int*)malloc(sizeof(unsigned int) * n);
	}

	mtrixMultiplication(k, k, S, k, n, G, pub);
	mtrixMultiplication(k, n, pub, n, n, P, key);

	errorGenerator(t, n, z);

	for (int i = 0; i < k; i++)
	{
		free(pub[i]);
	}
	free(pub);
}

void generateMatrixMultiple(const int& n, const int& m, const int& lim, unsigned int** matrix, unsigned int& idx, int* row, int rowIdx, int& rowValue)
{
	if (lim <= rowIdx)
	{
		int i = 0;
		while (i != rowIdx)
		{
			for (int j = 0; j < n; j++)
			{
				matrix[idx][j] *= matrix[row[i]][j];
			}
			i++;
		}
		idx++;
		return;
	}

	for (int i = rowValue - 1; i >= 1; i--)
	{
		row[rowIdx] = i;
		generateMatrixMultiple(n, m, lim, matrix, idx, row, rowIdx + 1, i);
	}

}
