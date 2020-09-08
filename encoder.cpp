#include <vector>
#include <iostream>
#include <math.h>
using namespace std;

void  encode(vector<int> &u) {
	int* temp = new int[u.size()/2];
	int* temp2 = new int[u.size() / 2];
	int m = log2(u.size());
	int N;
	int b;
	int position;
	for (int i = 0; i < m; i++) {
		N = pow(2, m - i);
		b = pow(2, i);
		for (int j = 0; j < b; j++) {
			position = j*N;
			for (int l = 0; l < N / 2; l++) {
				u[2*l+position] = u[2*l+position] ^ u[2*l+1+position];
				//u[position + l] = u[position + l] ^ u[position + N / 2 + l];
				temp[l] = u[2 * l + 1 + position];
				temp2[l] = u[2 * l + position];
			}
			for (int l = 0; l < N / 2; l++) {
				u[position + l] = temp2[l];
				u[position + N / 2 + l] = temp[l];
			}
		}
	}
	delete[] temp;
	delete[] temp2;
}


