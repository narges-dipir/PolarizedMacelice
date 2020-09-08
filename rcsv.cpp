#include <stdio.h>
#include <stdlib.h>	//srand
#include <time.h>	//time
#include <string.h>
#include "rcsv.h"
#include "mat.h"
#include "rm.h"
#include "opMath.h"
#include "msg.h"

rcsv::rcsv() {
}

rcsv::rcsv(const rcsv& orig) {
}

rcsv::~rcsv() {
}


void rcsvDecisionZero(double *orgMsg, int orgMsgLength, double *result, int resultLength)
{
	double sum = 0;
	for (int i = 0; i<resultLength; i++)
		sum += orgMsg[i];

	double everage = sum / orgMsgLength;
	for (int i = 0; i<resultLength; i++)
	{
		if (everage > 0)
			result[i] = 1;
		else if (everage < 0)
			result[i] = -1;
		else
			result[i] = (rand() % 2) == 0 ? 1 : -1;
	}
}

void rcsvDecisionEqual(double *orgMsg, int orgMsgLength, double *result, int resultLength)
{
	for (int i = 0; i<resultLength; i++)
	{
		if (orgMsg[i] > 0)
			result[i] = 1;
		else if (orgMsg[i] < 0)
			result[i] = -1;
		else
			result[i] = (rand() % 2) == 0 ? 1 : -1;
	}
}
double *rcsvLeft(int r, int m, double *orgMsg, int orgMsgLength, double *v, int vLength, double *result, int resultLength)
{
	double *yu;
	int yuLength = resultLength;

	yu = (double *)malloc(sizeof(double) * yuLength);

	for (int i = 0; i<yuLength; i++)
	{
		yu[i] = (orgMsg[i] + (orgMsg[i + yuLength] * v[i])) / 2;
	}

	result = (double *)malloc(sizeof(double) * resultLength);
	if (r == 0)
	{
		rcsvDecisionZero(yu, yuLength, result, resultLength);
	}
	else if (r == m)
	{
		rcsvDecisionEqual(yu, yuLength, result, resultLength);
	}
	else
	{
		double *v = NULL;
		int vLength = yuLength / 2;

		double *u = NULL;
		int uLength = yuLength / 2;

		v = rcsvRight(r - 1, m - 1, yu, yuLength, v, vLength);
		u = rcsvLeft(r, m - 1, yu, yuLength, v, vLength, u, uLength);

		for (int i = 0; i<uLength; i++)
		{
			result[i + uLength] = u[i] * v[i];
			result[i] = u[i];
		}

		free(v);
		free(u);
	}

	free(yu);

	return result;
}

double *rcsvRight(int r, int m, double *orgMsg, int orgMsgLength, double *result, int resultLength)
{
	double *yv;
	int yvLength = resultLength;

	yv = (double *)malloc(sizeof(double) * yvLength);
	for (int i = 0; i<yvLength; i++)
	{
		yv[i] = orgMsg[i] * orgMsg[i + yvLength];
	}

	result = (double *)malloc(sizeof(double) * resultLength);
	if (r == 0)
	{
		rcsvDecisionZero(yv, yvLength, result, resultLength);
	}
	else if (r == m)
	{
		rcsvDecisionEqual(yv, yvLength, result, resultLength);
	}
	else
	{
		double *v = NULL;
		int vLength = yvLength / 2;

		double *u = NULL;
		int uLength = yvLength / 2;

		v = rcsvRight(r - 1, m - 1, yv, yvLength, v, vLength);
		u = rcsvLeft(r, m - 1, yv, yvLength, v, vLength, u, uLength);

		for (int i = 0; i<uLength; i++)
		{
			result[i + (uLength)] = u[i] * v[i];
			result[i] = u[i];
		}

		free(v);
		free(u);
	}

	free(yv);

	return result;
}




void rcsvDecodingHard(int r, int m, unsigned int *recvCodeword, int recvCodewordLength)
{
	double *codeword;
	int codewordLength = power(m);
	codeword = (double *)malloc(sizeof(double) * codewordLength);

	double *result;
	result = (double *)malloc(sizeof(double) * power(m));

	for (int i = 0; i<recvCodewordLength; i++)
	{
		if (recvCodeword[i] == 1)
			codeword[i] = -1;
		else
			codeword[i] = 1;
	}

	double *v = NULL;
	int vLength = codewordLength / 2;
	double *u = NULL;
	int uLength = codewordLength / 2;

	v = rcsvRight(r - 1, m - 1, codeword, codewordLength, v, vLength);
	u = rcsvLeft(r, m - 1, codeword, codewordLength, v, vLength, u, uLength);

	for (int i = 0; i<uLength; i++)
	{
		result[i + uLength] = u[i] * v[i];
		result[i] = u[i];
	}

	for (int i = 0; i<recvCodewordLength; i++)
	{
		if (result[i] == -1)
		{
			recvCodeword[i] = 1;
		}
		else
		{
			recvCodeword[i] = 0;
		}
	}

	free(v);
	free(u);

	free(result);
	free(codeword);
}