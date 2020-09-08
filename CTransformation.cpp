#include <map>
#include <iostream>
#include <set>
#include <vector>
#include <array>
#include <math.h>
using namespace std; 
int Partition(vector <long double>& vect, int start, int end, int* helper) {

	int i = start - 1;
	long double elem = vect[end];
	for (int j = start; j <= end - 1; j++) {
		if (vect[j] <= elem) {
			i++;
			long double temp = vect.at(j);
			vect.at(j) = vect.at(i);
			vect.at(i) = temp;
			swap(helper[j], helper[i]);
		}
	}
	long double temp = vect.at(i+1);
	vect.at(i+1) = vect.at(end);
	vect.at(end) = temp;
	swap(helper[i+1], helper[end]);
	return i + 1;
}

void Sort(vector <long double>& vect, int start, int end, int* helper) {
	if (start < end) {
		int q = Partition(vect, start, end, helper);
		Sort(vect, start, q - 1, helper);
		Sort(vect, q + 1, end, helper);
	}
}
void cTransform(long double p, int n, int n_frozen, set<int>* frozen) {
	
	vector<long double> Z(n);
	Z[0] = p;
	long double T;
	int *helper = new int[n];
	for (int i = 0; i < n; i++) {
		helper[i] = i;
	}
	for (int i = 0; i < log2(n); i++) {
		int u = pow(2, i);
		for (int t = 0; t < n; t = t+n/u) {
			T = Z[t];
			Z[t] = 2 * T - pow(T, 2);
			Z[t+n/(2*u)] = pow(T, 2);
		}
	}
	Sort(Z, 0, n - 1,helper);
	for (int i = 0; i < n_frozen; i++) {
		frozen->insert(helper[n - 1 - i]);
	}
}



