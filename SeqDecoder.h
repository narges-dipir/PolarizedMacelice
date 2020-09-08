#pragma once
#include <stack>
#include <iostream>
#include <vector>
#include <set>
#include <queue> 
#include <list>
using namespace std;
typedef pair<double, int> Pairs;
struct Order
{
	bool operator()(Pairs const& a, Pairs const& b) const
	{
		return b.first<a.first;
	}
};
class SeqDecoder {
private:
	int L;
	int m;
	double dispersion;
	int n;/*size of y = 2^m*/
	set<int> u;
	//int* inactivePathIndices;
	vector <int> inactivePathIndices;
	vector<bool>activePath;
	int ****arrayPointer_C;
	int **pathIndexToArrayIndex;
	//int* inactiveArrayIndices;
	vector<vector<int> > inactiveArrayIndices;
	int **arrayReferenceCount;
	double *S;
	double *D;
	multiset<Pairs, Order> prQueue;
	int* q;			// how many paths are on each phase
	Pairs temppair;
	int *phi_l;  //phases of different active paths
	int P;
	int* out;

	vector<int>row;	// for what do I have it?
	
	int u2;
	int rho;
	int h;
	int j;
	int phi, psi;
	
	double* P_j;
	int theta;
	int** C_m;
	double *S_m;
public:
	

	SeqDecoder(int _L, int _n, set<int> _u, double _dispersion, double* _P_j, int _theta);

	~SeqDecoder();

	void initialization();

	int assignInitialPath();

	double* getArray_S(int lambda, int l);
	int** getArrayPointer_C(int lambda, int l);

	//void continuePaths_UnfrozenBit(int phi);

	void recursivelyUpdateC(int lambda, int phi, int l);

	void killPath(int l);

	int clonePath(int l);

	void recursivelyCalcS(int l, int lambda, int phi_l);

	int* decode(long double* y);

	double omega(int i);

	void check();

	inline double& getS(int l, int lambda, int beta);

	inline double* getSPointer(int l, int lambda);
	inline double& getD(int l, int m);

	//void copyArrays(int lambda, int s, int s2);
};