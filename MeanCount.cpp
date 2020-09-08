#include <iostream>
#include <algorithm>

#include <set>
#include <vector>
#define alpha -0.4527
#define beta 0.0218
#define gama 0.86

using namespace std;
#include <math.h>
struct Order2
{
	bool operator()(pair<int, double> const& a, pair<int, double> const& b) const
	{
		return b.second<a.second;
	}
} my_str;

double phi(double x) {
	return exp(alpha*pow(x, gama) + beta);
}

double phi_1(double x) {
	return pow((log(x) - beta) / alpha, 1/gama);
}



/*void meanCount(int n, double dispersion, double *P) {
	//pair<int, double>* Lp = new pair<int, double>[n];
	double * L = new double[n];
	L[0] = 2 / dispersion;
	//get<1>(Lp[0]) = 2 / dispersion;
	double T;
	for (int i = 0; i < log2(n); i++) {
		int u = pow(2, i);
		for (int t = 0; t < n; t = t + n / u) {
			T = L[t];
			L[t] = phi_1((1 - pow(1 - phi(T), 2)));
			L[t + n / (2 * u)] = 2 * T;
			//T = get<1>(Lp[t]);
			//get<1>(Lp[t]) = phi_1((1 - pow(1 - phi(T), 2)));
			//get<1>(Lp[t + n / (2 * u)]) = 2 * T;
		}
	}
	//double *P = new double[n];				//��������� �� ���, � ������ �� ������, ������� ������������?
	for (int i = 0; i < n; i++) {
		P[i] = 0.5*erfc(sqrt(L[i]) / 2);
	}
}
*/

void meanCount(int n, double dispersion, double *P, int n_frozen, set<int>* frozen) {
	vector<pair<int, double> > sortarr;
	sortarr.resize(n);
	double * L = new double[n];
	L[0] =  2 / dispersion;
	double T;
	for (int i = 0; i < log2(n); i++) {
		int u = pow(2, i);
		for (int t = 0; t < n; t = t + n / u) {
			T = L[t];
			L[t] = phi_1((1 - pow(1 - phi(T), 2)));
			L[t + n / (2 * u)] = 2 * T;
		}
	}
										
	for (int i = 0; i < n; i++) {
		P[i] = 0.5*erfc(sqrt(L[i])/2);
		sortarr[i] = make_pair(i, P[i]);
	}
	sort(sortarr.begin(), sortarr.end(), my_str);

	for (int i = 0; i < n_frozen; i++) {
		frozen->insert(sortarr[i].first);
	}
	sortarr.clear();
	sortarr.shrink_to_fit();
	delete[] L;
}