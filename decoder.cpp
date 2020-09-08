#include "decoder.h"
#define _USE_MATH_DEFINES
#include <math.h>
using namespace std;
decoder::decoder(int _L, int _n, set<int> _u, long double _dispersion) {
	L = _L;
	n = _n;
	m = log2(n);
	u = _u;
	dispersion = _dispersion;
	inactivePathIndices.reserve(L);
	activePath.resize(L);

	arrayPointer_P = new long double***[m + 1];
	for (s = 0; s < m + 1; s++)
		arrayPointer_P[s] = new long double**[L];

	arrayPointer_C = new int***[m + 1];
	for (s = 0; s < m + 1; s++)
		arrayPointer_C[s] = new int**[L];

	pathIndexToArrayIndex = new int*[m + 1];
	for (s = 0; s < m + 1; s++)
		pathIndexToArrayIndex[s] = new int[L];

	row.reserve(L);
	for (s = 0; s < m + 1; s++) {
		inactiveArrayIndices.push_back(row);
	}

	arrayReferenceCount = new int*[m + 1];
	for (s = 0; s< m + 1; s++)
		arrayReferenceCount[s] = new int[L];

	for (lambda = 0; lambda <= m; lambda++) {
		for (s = 0; s < L; s++) {
			arrayPointer_P[lambda][s] = new long double*[(int)pow(2, m - lambda)];
			//cout << "(int)pow(2, m - lambda) = " << (int)pow(2, m - lambda) << endl;
			arrayPointer_C[lambda][s] = new int*[(int)pow(2, m - lambda)];
			for (int i = 0; i < (int)pow(2, m - lambda); i++) {
				arrayPointer_P[lambda][s][i] = new long double[2];
				arrayPointer_C[lambda][s][i] = new int[2];
			}
			//arrayReferenceCount[lambda][s] = 0;
			//inactiveArrayIndices[lambda].push_back(s);
		}
	}
	probForks = new long double*[L];
	for (l = 0; l< L; l++)
		probForks[l] = new long double[2];

	contForks = new bool*[L];
	for (l = 0;l < L; l++)
		contForks[l] = new bool[2];

	temp_cont.resize(2 * L);

	c = new int[n];
}
decoder::~decoder() {
	for (s = 0; s < inactiveArrayIndices.size(); s++){
		inactiveArrayIndices.at(s).clear();
		inactiveArrayIndices.at(s).shrink_to_fit();
	}
	inactiveArrayIndices.clear();
	inactiveArrayIndices.shrink_to_fit();
	activePath.clear();
	activePath.shrink_to_fit();

	for (s = 0; s < m + 1; s++) {
		delete[] pathIndexToArrayIndex[s];
	}
	delete[] pathIndexToArrayIndex;

	for (lambda = 0; lambda < m + 1; lambda++) {
		for (s = 0; s < L; s++) {
			for (int i = 0; i < (int)pow(2, m - lambda); i++) {
				delete[] arrayPointer_P[lambda][s][i];
				delete[] arrayPointer_C[lambda][s][i];
			}
			delete[] arrayPointer_P[lambda][s];
			delete[] arrayPointer_C[lambda][s];
			
		}
		delete[] arrayPointer_P[lambda];
		delete[] arrayPointer_C[lambda];
		delete[] arrayReferenceCount[lambda];
	}
	delete[] arrayReferenceCount;

	delete[] arrayPointer_P;
	
	delete[] arrayPointer_C;
	inactivePathIndices.clear();
	inactivePathIndices.shrink_to_fit();
	delete[] c;
	temp_cont.clear();
	temp_cont.shrink_to_fit();
	delete[] probForks;
	delete[] contForks;

	row.clear();
	row.shrink_to_fit();
}

void decoder::initialization() {

	for (lambda = 0; lambda <= m; lambda++) {
		for (s = 0; s < L; s++) {
			arrayReferenceCount[lambda][s] = 0;
			inactiveArrayIndices[lambda].push_back(s);
		}
	}
	for (l = 0; l < L; l++) {
		activePath[l] = false;
		inactivePathIndices.push_back(l);
	}

}
int decoder::assignInitialPath() {
	l = inactivePathIndices[inactivePathIndices.size() - 1];
	inactivePathIndices.pop_back();
	activePath[l] = true;
	for (lambda = 0; lambda <= m; lambda++) {
		s = inactiveArrayIndices[lambda][inactiveArrayIndices[lambda].size() - 1];
		inactiveArrayIndices[lambda].pop_back();
		pathIndexToArrayIndex[lambda][l] = s;
		
		arrayReferenceCount[lambda][s] = 1;
	}
	return l;
}

long double** decoder::getArrayPointer_P(int lambda, int l) {
	s = pathIndexToArrayIndex[lambda][l];
	if (arrayReferenceCount[lambda][s] == 1) {
		s2 = s;
	}
	else {
		s2 = inactiveArrayIndices[lambda][inactiveArrayIndices[lambda].size() - 1];
		inactiveArrayIndices[lambda].pop_back();
		for (int i = 0; i < (int)pow(2, m - lambda); i++) { //(int)pow(2, m - lambda(+1 or not));
			arrayPointer_P[lambda][s2][i][0] = arrayPointer_P[lambda][s][i][0];
			arrayPointer_P[lambda][s2][i][1] = arrayPointer_P[lambda][s][i][1];

			arrayPointer_C[lambda][s2][i][0] = arrayPointer_C[lambda][s][i][0];
			arrayPointer_C[lambda][s2][i][1] = arrayPointer_C[lambda][s][i][1];
		}

		arrayReferenceCount[lambda][s]--;
		arrayReferenceCount[lambda][s2] = 1;
		pathIndexToArrayIndex[lambda][l] = s2;
	}
	return arrayPointer_P[lambda][s2];
}
int** decoder::getArrayPointer_C(int lambda, int l) {
	s = pathIndexToArrayIndex[lambda][l];
	if (arrayReferenceCount[lambda][s] == 1) {
		s2 = s;
	}
	else {
		s2 = inactiveArrayIndices[lambda][inactiveArrayIndices[lambda].size() - 1];
		inactiveArrayIndices[lambda].pop_back();
		for (int i = 0; i < (int)pow(2, m - lambda); i++) {
			arrayPointer_P[lambda][s2][i][0] = arrayPointer_P[lambda][s][i][0];
			arrayPointer_P[lambda][s2][i][1] = arrayPointer_P[lambda][s][i][1];

			arrayPointer_C[lambda][s2][i][0] = arrayPointer_C[lambda][s][i][0];
			arrayPointer_C[lambda][s2][i][1] = arrayPointer_C[lambda][s][i][1];
		}
		arrayReferenceCount[lambda][s]--;
		arrayReferenceCount[lambda][s2] = 1;
		pathIndexToArrayIndex[lambda][l] = s2;
	}
	return arrayPointer_C[lambda][s2];
}
void decoder::recursivelyCalcP(int lambda, int phi) {
	long double **P_lambda, **P_prelambda;
	int **C_lambda;
	if (lambda == 0) return;
	psi = phi / 2;
	if (phi % 2 == 0)
		recursivelyCalcP(lambda - 1, psi);
	long double sigma2 = 0;
	for (l = 0; l < L; l++) {
		if (activePath[l] == false)
			continue;
		P_lambda = getArrayPointer_P(lambda, l);
		P_prelambda = getArrayPointer_P(lambda - 1, l);
		C_lambda = getArrayPointer_C(lambda, l);
		int tml = pow(2, m - lambda);
		for (beta = 0; beta < tml; beta++) {
			if (phi % 2 == 0) {
				//P_lambda[beta][0] = (P_prelambda[2*beta][0] * P_prelambda[2*beta + 1][0] + P_prelambda[2*beta][1] * P_prelambda[2*beta + 1][1]) / 2;
				P_lambda[beta][0] = (P_prelambda[beta][0] * P_prelambda[beta + tml][0] + P_prelambda[beta][1] * P_prelambda[beta + tml][1]) / 2;
				sigma2 = sigma2 > P_lambda[beta][0] ? sigma2 : P_lambda[beta][0];

				//P_lambda[beta][1] = (P_prelambda[2*beta][1] * P_prelambda[2*beta + 1][0] + P_prelambda[2*beta][0] * P_prelambda[2*beta + 1][1]) / 2;
				P_lambda[beta][1] = (P_prelambda[beta][1] * P_prelambda[beta + tml][0] + P_prelambda[beta][0] * P_prelambda[beta + tml][1]) / 2;
				sigma2 = sigma2 > P_lambda[beta][1] ? sigma2 : P_lambda[beta][1];
			}
			else {
				u2 = C_lambda[beta][0];
				//P_lambda[beta][0] = P_prelambda[2 * beta][u] * P_prelambda[2 * beta + 1][0] / 2;
				P_lambda[beta][0] = P_prelambda[beta][u2] * P_prelambda[beta + tml][0] / 2;
				sigma2 = sigma2 > P_lambda[beta][0] ? sigma2 : P_lambda[beta][0];

				//P_lambda[beta][1] = P_prelambda[2 * beta][u2 ^ 1] * P_prelambda[2 * beta + 1][1] / 2;
				P_lambda[beta][1] = P_prelambda[beta][u2 ^ 1] * P_prelambda[beta + tml][1] / 2;
				sigma2 = sigma2 > P_lambda[beta][1] ? sigma2 : P_lambda[beta][1];
			}
		}
	}
	//sigma2 = 1;
	for (l = 0; l < L; l++) {
		if (activePath[l] == false)
			continue;
		P_lambda = getArrayPointer_P(lambda, l);
		for (beta = 0; beta < pow(2, m - lambda); beta++) {
			P_lambda[beta][0] = P_lambda[beta][0] / sigma2;
			P_lambda[beta][1] = P_lambda[beta][1] / sigma2;
		}
	}

}

void decoder::continuePaths_UnfrozenBit(int phi) {
	//long double **probForks = new long double*[L];
	//for (int count = 0; count < L; count++)
	//	probForks[count] = new long double[2];

	int i = 0;
	//
	for (l = 0; l < L; l++) {
		if (activePath[l] == true) {
			long double **P_m = getArrayPointer_P(m, l);
			probForks[l][0] = P_m[0][0];
			probForks[l][1] = P_m[0][1];
			i++;
		}
		else {
			probForks[l][0] = -1;
			probForks[l][1] = -1;
		}
	}
	//bool **contForks = new bool*[L];
//	for (int count = 0; count < L; count++)
	//	contForks[count] = new bool[2];
	rho = (2 * i <= L) ? 2 * i : L;// min(2 * i, L);
	i = 0;
	//vector<long double> temp_cont;
	//temp_cont.resize(2 * L);
	//long double* temp = new long double[2 * L];
	for (h = 0; h < L; h++) {
		for (j = 0; j < 2; j++) {
			temp_cont[i] = probForks[h][j];
		
			i++;
		}
	}
/*
	for (int j = 0; j < 2; j++) {
	for (int h = 0; h < L; h++) {
		
			cout << probForks[h][j] << " ";
		}
		cout << endl;
	}*/
	Sort(temp_cont, 0, i - 1);
	def = temp_cont[2 * L - rho];
	i = 0;
	for (h = 0; h < L; h++) {
		for (j = 0; j < 2; j++) {
			if ((probForks[h][j] > def) && (i < rho)) {  //add stoppin' condition
				contForks[h][j] = true;
				i++;
			}
			else
				contForks[h][j] = false;
		}
	}
	//in still < rho, look for some equal
	for (h = 0; h < L; h++) {
		for (j = 0; j < 2; j++) {
			if (i < rho) {
				if (probForks[h][j] == def) {
					contForks[h][j] = true;
					i++;
				}
			}
		}
	}
	/*for (int j = 0; j < 2; j++) {
	for (int h = 0; h < L; h++) {
		
			cout << contForks[h][j] << " ";
		}
		cout << endl;
	}
	*/
	for (l = 0; l < L; l++) {
		if (activePath[l] == false)
			continue;
		if ((contForks[l][0] == false) && (contForks[l][1] == false))
			killPath(l);
	}
	for (l = 0; l < L; l++) {
		if ((contForks[l][0] == false) && (contForks[l][1] == false))
			continue;
		int** C_m = getArrayPointer_C(m, l);
		if ((contForks[l][0] == true) && (contForks[l][1] == true)) {
			C_m[0][phi % 2] = 0;
			l2 = clonePath(l);
			C_m = getArrayPointer_C(m, l2);
			C_m[0][phi % 2] = 1;
		}
		//exactly one fork is good
		else {
			if (contForks[l][0] == true) {
				C_m[0][phi % 2] = 0;
			}
			else
				C_m[0][phi % 2] = 1;
		}
	}
}

int decoder::clonePath(int l) {
	l2 = inactivePathIndices[inactivePathIndices.size() - 1];
	inactivePathIndices.pop_back();
	activePath[l2] = true;

	for (lambda = 0; lambda <= m; lambda++) {
		s = pathIndexToArrayIndex[lambda][l];
		pathIndexToArrayIndex[lambda][l2] = s;
		arrayReferenceCount[lambda][s]++;
	}
	return l2;
}
void decoder::killPath(int l) {

	activePath[l] = false;
	inactivePathIndices.push_back(l);
	for (lambda = 0; lambda <= m; lambda++) {
		s = pathIndexToArrayIndex[lambda][l];
		arrayReferenceCount[lambda][s]--;
		if (arrayReferenceCount[lambda][s] == 0)
			inactiveArrayIndices[lambda].push_back(s);
	}
}
int decoder::Partition(vector<long double>& vector, int start, int end) {

	s = start - 1;
	long double elem = vector[end];
	//exchange [start + index] with [end]
	for (j = start; j < end; j++) {
		if (vector[j] <= elem) {
			//exchange [j] and [i]
			s++;
			swap(vector[j], vector[s]);
		}
	}
	//exchange [i + 1] and [end] and exchange [ i + 1 + .size()/2] and [end + .size()/2]
	swap(vector[s + 1], vector[end]);
	return s + 1;
}

void decoder::Sort(vector<long double>& vector, int start, int end) {
	if (start < end) {
		h = Partition(vector, start, end);
		Sort(vector, start, h - 1);
		Sort(vector, h + 1, end);
	}
}

void decoder::recursivelyUpdateC(int lambda, int phi) {
	psi = phi / 2;
	for (l = 0; l < L; l++) {
		if (activePath[l] == false)
			continue;
		int** C_lambda = getArrayPointer_C(lambda, l);
		int** C_prelambda = getArrayPointer_C(lambda - 1, l);
		h = pow(2, m - lambda);
		for (beta = 0; beta < h; beta++) {
			C_prelambda[beta][psi % 2] = C_lambda[beta][0] ^ C_lambda[beta][1];
			C_prelambda[beta + h][psi % 2] = C_lambda[beta][1];
		}
	}
	if (psi % 2 == 1)
		recursivelyUpdateC(lambda - 1, psi);

}
int* decoder::decode(long double* y) {
	initialization();

	int** C_m;
	long double** P_m;
	l = assignInitialPath();
	long double** P_0 = getArrayPointer_P(0, l);

	for (beta = 0; beta < n; beta++) {
		//double L = 2 * y[beta] / dispersion;
		
		//P_0[beta][1] = 1. / (exp(L) + 1);
		//P_0[beta][0] = 1. - P_0[beta][1];
		
		P_0[beta][0] = exp(-pow((y[beta] - 1), 2) / (2 * dispersion)) / sqrt(2 * M_PI * dispersion);
		
		P_0[beta][1] = exp(-pow((y[beta] + 1), 2) / (2 * dispersion)) / sqrt(2 * M_PI * dispersion);
	}
	
	for (phi = 0; phi < n; phi++) {
		recursivelyCalcP(m, phi);
		
		if (u.find(phi) != u.end()/*u[phi] == false*/) {
			for (l = 0; l < L; l++) {
				if (activePath[l] == false)
					continue;
				C_m = getArrayPointer_C(m, l);
				C_m[0][phi % 2] = 0; // zero or not zero? or what?!
			}
			
		}
		else {
			continuePaths_UnfrozenBit(phi);
		}
		
		if (phi % 2 == 1)
			recursivelyUpdateC(m, phi);
	
	}
	
	l2 = 0;
	p = 0;

	for (l = 0; l < L; l++) {
		if (activePath[l] == false)
			continue;
		C_m = getArrayPointer_C(m, l);
		P_m = getArrayPointer_P(m, l);

		if (p < P_m[0][C_m[0][1]]) {
			l2 = l;
			p = P_m[0][C_m[0][1]];
		}
	}

	int** C_0 = getArrayPointer_C(0, l2);
	//int* c = new int[n];
	
	for (s = 0; s < n; s++)
		c[s] = C_0[s][0];
	return c;
}

void decoder::check() {
	
	/*cout << "arryReferenceCount content: " << endl;
	for (int lambda = 0; lambda < m + 1; lambda++) {
		for (int s = 0; s < L; s++) {
			cout << "arrayReferenceCount[" << lambda << "][" << s << "]: " << arrayReferenceCount[lambda][s] << endl;
		}
		cout << "--------" << endl;
	}*/
	/*
	cout << "inactivePathIndices content (stack):" << endl;
	for (int i = 0; i < inactivePathIndices.size(); i++)
		cout << inactivePathIndices[i] << " ";
	cout << endl;


	cout << "pathIndexToArrayIndex content:" << endl;
	for (int lambda = 0; lambda < m + 1; lambda++) {
		for (int s = 0; s < L; s++)
			cout << "pathIndexToArrayIndex[" << lambda << "][" << s << " ] is " << pathIndexToArrayIndex[lambda][s] << endl;
	}
	cout << "------" << endl;

	cout << "activePath" << endl;
	for (int a = 0; a < L; a++) {
		cout << activePath[a] << " ";
	}
	cout << endl;
	*/
	cout << "Array Pointer P content: " << endl;
	for (int lambda = 0; lambda < m + 1; lambda++) {
		for (int s = 0; s < L; s++) {
			cout << "arrayPointer_P[ " << lambda << " ][ " << s << " ]" << endl;
			for (int i = 0; i < (int)pow(2.0, m - lambda); i++) {
				cout << arrayPointer_P[lambda][s][i][0] << " ";
				//cout << endl;
				//for (int i = 0; i < (int)pow(2.0, m - lambda); i++)
				cout << arrayPointer_P[lambda][s][i][1] << " ";
				cout << endl;
			}
			cout << endl;
			cout << "--------" << endl;
		}
	}

	
	cout << "Array Pointer C content: " << endl;
	for (int lambda = 0; lambda < m + 1; lambda++) {
		for (int s = 0; s < L; s++) {
			cout << "arrayPointer_C[ " << lambda << " ][ " << s << " ]" << endl;
			int** te = arrayPointer_C[lambda][s];
			for (int i = 0; i < (int)pow(2.0, m - lambda); i++) {
				cout << te[i][0] << " ";
				//cout << arrayPointer_C[lambda][s][i][0] << " ";
				//cout << endl;
				//for (int i = 0; i < (int)pow(2.0, m - lambda); i++)
				cout << te[i][1] << " ";
				//cout << arrayPointer_C[lambda][s][i][1] << " ";
				cout << endl;
			}
			cout << endl;
			cout << "--------" << endl;
		}
	}
}
