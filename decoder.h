#pragma once
#include <stack>
#include <iostream>
#include <vector>
#include <set>
using namespace std;
class decoder{
private:
	int L;
	int m;
	long double dispersion;
	int n;/*size of y = 2^m*/
	set<int> u; 
	vector <int> inactivePathIndices;
	vector<bool>activePath;
	long double ****arrayPointer_P;
	int ****arrayPointer_C;
	int **pathIndexToArrayIndex;
	vector<vector<int> > inactiveArrayIndices;
	int **arrayReferenceCount;


	vector<int>row;
	long double **probForks;
	bool **contForks;
	vector<long double> temp_cont;
	int* c;
	int lambda;
	int s,s2;
	int l, l2;
	int beta;
	int u2;
	int rho;
	int h;
	int j;
	long double def;
	int phi, psi;
	int p;

public:
	decoder(int L, int n, set<int> u, long double sigma);

	~decoder();

	void initialization();

	int assignInitialPath();

	long double** getArrayPointer_P(int lambda, int l);

	int** getArrayPointer_C(int lambda, int l);

	void recursivelyCalcP(int lambda, int phi);

	void continuePaths_UnfrozenBit(int phi);

	void recursivelyUpdateC(int lambda, int phi);

	int Partition(vector<long double>& vector, int start, int end);

	void Sort(vector<long double>& vector, int start, int end);

	void killPath(int l);

	int clonePath(int l);

	int* decode(long double* y);

	void check();
};