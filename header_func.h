#pragma once
#include <vector>
using namespace std;
void encode(vector<int> &u);
void cTransform(long double p, int n_inf, int n_frozen, set<int>* frozen);
void sortTest();
void noise(long double mean, long double sigma, int n, long double *a);
void meanCount(int n, double dispersion, double *P, int n_frozen, set<int>* frozen);
