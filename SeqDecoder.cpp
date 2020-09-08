#include "SeqDecoder.h"
#define _USE_MATH_DEFINES
#include <math.h>
using namespace std;

double SeqDecoder::omega(int i) {
	double result = 1;
	for (int k = i+1; k < n; k++) {
		if (u.find(k) != u.end()/*u[phi] == false*/) {
			if (P_j[k] > 0.31)
				P_j[k] = 0.31;
			result = result*(1 - P_j[k]);
		}
	}
	return result;
}
inline double& SeqDecoder::getS(int l, int lambda, int beta) {
	return S[(lambda*theta + l)*n + beta];
}

inline double* SeqDecoder::getSPointer(int l, int lambda) {
	return &S[(lambda*theta + l)*n];
}

inline double& SeqDecoder::getD(int l, int m) {
	int s = pathIndexToArrayIndex[m][l];
	return D[m*theta + s];
}
SeqDecoder::SeqDecoder(int _L, int _n, set<int> _u, double _dispersion, double* _P_j, int _theta) {
	L = _L;
	n = _n;
	m = log2(n);
	u = _u;
	dispersion = _dispersion;
	theta = _theta;
	inactivePathIndices.reserve(theta);
	activePath.resize(theta);
	int s;
		arrayPointer_C = new int***[m + 1];
	for (s = 0; s < m + 1; s++)
		arrayPointer_C[s] = new int**[theta];

	pathIndexToArrayIndex = new int*[m + 1];
	for (s = 0; s < m + 1; s++)
		pathIndexToArrayIndex[s] = new int[theta];

	//inactiveArrayIndices = new int[theta*(m + 1)];

	row.reserve(theta);
	for (s = 0; s < m + 1; s++) {
		inactiveArrayIndices.push_back(row);
	}

	arrayReferenceCount = new int*[m + 1];
	for (s = 0; s< m + 1; s++)
		arrayReferenceCount[s] = new int[theta];

	for (int lambda = 0; lambda <= m; lambda++) {
		for (s = 0; s < theta; s++) {
			arrayPointer_C[lambda][s] = new int*[(int)pow(2, m - lambda)];
			for (int i = 0; i < (int)pow(2, m - lambda); i++) {
				arrayPointer_C[lambda][s][i] = new int[2];
			}
		}
	}

	out = new int[n];

	S = new double[theta*(m+1)*n]; //S_l_lambda[beta] = S[(lambda*theta+l)*n+beta]

	D = new double[theta*(m+1)]; //D_l_m = D[m*theta+l]
	q = new int[n+1];	

	P_j = _P_j;
	
	phi_l = new int[theta];
}

SeqDecoder::~SeqDecoder() {
	//cout << "Destructor!" << endl;
	delete[] q;
	delete[] S;
	delete[] D;
	
	for (int s = 0; s < inactiveArrayIndices.size(); s++) {
		inactiveArrayIndices.at(s).clear();
		inactiveArrayIndices.at(s).shrink_to_fit();
	}
	inactiveArrayIndices.clear();
	inactiveArrayIndices.shrink_to_fit();
	activePath.clear();
	activePath.shrink_to_fit();


delete[] phi_l;
for (int s = 0; s < m + 1; s++) {
		delete[] pathIndexToArrayIndex[s];
	}
	delete[] pathIndexToArrayIndex;
	
	for (int lambda = 0; lambda < m + 1; lambda++) {
		for (int s = 0; s < theta; s++) {
			for (int i = 0; i < (int)pow(2, m - lambda); i++) {
				delete[] arrayPointer_C[lambda][s][i];
			}	
			delete[] arrayPointer_C[lambda][s];
		}
		delete[] arrayPointer_C[lambda];
		delete[] arrayReferenceCount[lambda];
	}
	delete[] arrayReferenceCount;
	delete[] arrayPointer_C;
	inactivePathIndices.clear();
	inactivePathIndices.shrink_to_fit();
	delete[] out;
	row.clear();
	row.shrink_to_fit();
}

void SeqDecoder::initialization() {
	int l;
	int s;
	for (int lambda = 0; lambda <= m; lambda++) {
		for (s = 0; s < theta; s++) {
			arrayReferenceCount[lambda][theta -s-1] = 0;
			inactiveArrayIndices[lambda].push_back(theta -s-1);
		}
	}
	for (l = 0; l < theta; l++) {
		activePath[theta - l - 1] = false;
		inactivePathIndices.push_back(theta -l-1);
		if (l > theta)
			cout << "ALERT!" << endl;
		phi_l[l] = 0;
	}
	for (l = 0; l < n; l++) {
		q[l] = 0;
	}
}

int SeqDecoder::assignInitialPath() {
	int l;
	int s;
	l = inactivePathIndices[inactivePathIndices.size() - 1];
	inactivePathIndices.pop_back();
	activePath[l] = true;
	for (int lambda = 0; lambda <= m; lambda++) {
		s = inactiveArrayIndices[lambda][inactiveArrayIndices[lambda].size() - 1];
		inactiveArrayIndices[lambda].pop_back();
		if ((lambda>m) || (l>=theta))
			cout << "ALERT2" << endl;
		pathIndexToArrayIndex[lambda][l] = s;

		arrayReferenceCount[lambda][s] = 1;
	}
	return l;
}

double* SeqDecoder::getArray_S(int lambda, int l) {
	int s, s2;
	if ((lambda>m) || (l>=theta))
		cout << "ALERT2" << endl;
	s = pathIndexToArrayIndex[lambda][l];
	if (arrayReferenceCount[lambda][s] == 1) {
		s2 = s;
	}
	else {
		s2 = inactiveArrayIndices[lambda][inactiveArrayIndices[lambda].size() - 1];
		inactiveArrayIndices[lambda].pop_back();
		for (int i = 0; i < (int)pow(2, m - lambda); i++) { 
			arrayPointer_C[lambda][s2][i][0] = arrayPointer_C[lambda][s][i][0];			//make a good function for that stuff finally!
			arrayPointer_C[lambda][s2][i][1] = arrayPointer_C[lambda][s][i][1];
			S[(lambda*theta + s2)*n + i] = S[(lambda*theta + s)*n + i];
			D[lambda*theta + s2] = D[lambda*theta + s];
		}

		arrayReferenceCount[lambda][s]--;
		arrayReferenceCount[lambda][s2] = 1;
		if ((lambda>m) || (l>=theta))
			cout << "ALERT2" << endl;
		pathIndexToArrayIndex[lambda][l] = s2;
	}
	return getSPointer(s2, lambda);
}

/*void SeqDecoder::copyArrays(int lambda, int s, int s2) {
	for (int i = 0; i < (int)pow(2, m - lambda); i++) {
		memcpy(*(arrayPointer_C[lambda][s2]), *(arrayPointer_C[lambda][s]), pow(2,(i+1))*sizeof(int));
	//	arrayPointer_C[lambda][s2][i][0] = arrayPointer_C[lambda][s][i][0];		// awful thing like "NO MEMCPY! NONONO"
	//	arrayPointer_C[lambda][s2][i][1] = arrayPointer_C[lambda][s][i][1];
		//memcpy(&S[(lambda*theta + s2)*n], &S[(lambda*theta + s)*n], (i + 1)*sizeof(double));
		S[(lambda*theta + s2)*n + i] = S[(lambda*theta + s)*n + i];
		D[lambda*theta + s2] = D[lambda*theta + s];
	}
}*/

int** SeqDecoder::getArrayPointer_C(int lambda, int l) {
	if ((lambda>m) || (l>=theta))
		cout << "ALERT2" << endl;
	int s = pathIndexToArrayIndex[lambda][l];		// lambda becomes -1 sometimes. WHY. well, cause you're dummy
	int s2;
	if (arrayReferenceCount[lambda][s] == 1) {
		s2 = s;
	}
	else {
		s2 = inactiveArrayIndices[lambda][inactiveArrayIndices[lambda].size() - 1];
		inactiveArrayIndices[lambda].pop_back();

		//copyArrays(lambda, s, s2);
		for (int i = 0; i < (int)pow(2, m - lambda); i++) {
			arrayPointer_C[lambda][s2][i][0] = arrayPointer_C[lambda][s][i][0];		// awful thing like "NO MEMCPY! NONONO"
			arrayPointer_C[lambda][s2][i][1] = arrayPointer_C[lambda][s][i][1];
			S[(lambda*theta + s2)*n + i] = S[(lambda*theta + s)*n + i];
			D[lambda*theta + s2] = D[lambda*theta + s];
		}
		arrayReferenceCount[lambda][s]--;
		arrayReferenceCount[lambda][s2] = 1;
		if ((lambda>m) || (l>=theta))
			cout << "ALERT2" << endl;
		pathIndexToArrayIndex[lambda][l] = s2;
	}
	return arrayPointer_C[lambda][s2];
}

int SeqDecoder::clonePath(int l) {
	int l2;
	int s;
	l2 = inactivePathIndices[inactivePathIndices.size() - 1];
	inactivePathIndices.pop_back();
	activePath[l2] = true;
	
	for (int lambda = 0; lambda <= m; lambda++) {
		s = pathIndexToArrayIndex[lambda][l];
		if ((lambda>m) || (l2>=theta))
			cout << "ALERT2" << endl;
		pathIndexToArrayIndex[lambda][l2] = s;
		arrayReferenceCount[lambda][s]++;
	}
	return l2;
}
void SeqDecoder::killPath(int l) {
	int s;
	activePath[l] = false;
	inactivePathIndices.push_back(l);
	for (int lambda = 0; lambda <= m; lambda++) {
		s = pathIndexToArrayIndex[lambda][l];
		if ((lambda>m) || (s>=theta))
			cout << "ALERT2" << endl;
		arrayReferenceCount[lambda][s]--;
		if (arrayReferenceCount[lambda][s] == 0)
			inactiveArrayIndices[lambda].push_back(s);
	}
}
void SeqDecoder::recursivelyUpdateC(int lambda, int phi, int l) {
	psi = phi / 2;
	if (activePath[l] == false)		return;

		int** C_lambda = getArrayPointer_C(lambda, l);
		int** C_prelambda = getArrayPointer_C(lambda - 1, l);
		h = pow(2, m - lambda);
		for (int beta = 0; beta < h; beta++) {
			C_prelambda[2*beta][psi % 2] = C_lambda[beta][0] ^ C_lambda[beta][1];
			C_prelambda[2*beta + 1][psi % 2] = C_lambda[beta][1];
		}
	if (psi % 2 == 1)
		recursivelyUpdateC(lambda - 1, psi,l);

}



double Q(double a, double b) {
	if ((a < 0 && b < 0) || (a>0 && b>0))
		return (abs(a) < abs(b)) ? abs(a) : abs(b);
	return (abs(a) < abs(b)) ? abs(a)*(-1) : abs(b)*(-1);
}

void SeqDecoder::recursivelyCalcS(int l, int lambda, int phi) {
	if (lambda == 0) return;
	if (phi % 2 == 0) recursivelyCalcS(l, lambda - 1, phi / 2); 
	double* S_lambda = getArray_S(lambda, l);
	double* S_prelambda = getArray_S(lambda-1, l); // this is for copying
	getD(l, lambda) = getD(l, lambda - 1);
	for (int beta = 0; beta < pow(2, m - lambda); beta++) {
		if (phi % 2 == 0) {
			*(S_lambda + beta) = Q(*(S_prelambda + 2 * beta), *(S_prelambda + 2 * beta + 1));
			getD(l, lambda) = getD(l, lambda) + max(*(S_prelambda + 2 * beta), *(S_prelambda + 2 * beta + 1));
		}
		else{
		if (arrayPointer_C[lambda][pathIndexToArrayIndex[lambda][l]][beta][0] == 0) {
			*(S_lambda + beta) = *(S_prelambda + 2 * beta) + *(S_prelambda + 2 * beta + 1);
		}
		else {

			*(S_lambda + beta) = *(S_prelambda + 2 * beta + 1) - *(S_prelambda + 2 * beta);
			getD(l, lambda) = getD(l, lambda) + *(S_prelambda + 2 * beta);
		
		}
	}
	}
}

int* SeqDecoder::decode(long double* y) {
	int s;
	initialization();
	//int** C_m;
//	int _l;
	//long double** P_m;
	int l, l2, l_temp;
	l = assignInitialPath();
	//double *S_m;


	for (int yt = 0; yt < theta; yt++) {
		if (yt > theta)
			cout << "ALERT!" << endl; 
		phi_l[yt] = 0; }

	for (int beta = 0; beta < n; beta++) {		
		getS(l,0,beta) = 2 * y[beta] / dispersion;
		q[beta] = 0;		
	}
	getD(l, 0) = 0;
	pair<double, int> pair(0, l);	
	prQueue.emplace(pair);
	P = 1;
	while (P > 0) {
		temppair = *prQueue.begin();
		l = temppair.second;
		prQueue.erase(prQueue.begin());
		P = P - 1;
		if (phi_l[l] > n) {
			cout << "ALEEEEEEEERT!" << endl;
		}
		q[phi_l[l]] ++;
		if (phi_l[l] == n) {
			for (s = 0; s < n; s++)
				out[s] = arrayPointer_C[0][pathIndexToArrayIndex[0][l]][s][0];



			prQueue.clear();															//Strange stuff
			for (int s = 0; s < inactiveArrayIndices.size(); s++) {
				inactiveArrayIndices.at(s).clear();
				inactiveArrayIndices.at(s).shrink_to_fit();
			}
			inactivePathIndices.clear();
			inactivePathIndices.shrink_to_fit();



			return out;
		}
		recursivelyCalcS(l, m, phi_l[l]);
		if (u.find(phi_l[l]) != u.end()) {
			getArrayPointer_C(m, l)[0][phi_l[l] % 2] = 0;
			
			prQueue.emplace(Pairs(getS(pathIndexToArrayIndex[m][l], m, 0) + getD(l, m) + log(omega(phi_l[l])), l));
			P++;
			if (phi_l[l] % 2 == 1) {
				recursivelyUpdateC(m, phi_l[l], l);			
			}
			if (l >= theta)
				cout << "ALERT!" << endl;
			phi_l[l]++;
			
		}
		else {
			while (P >= theta - 1) {
				std::set<Pairs>::iterator it = prQueue.end();
				it--;
				temppair = *it;
				prQueue.erase(it); //PopMin
				l_temp = temppair.second;
				killPath(l_temp);
				P = P - 1;
			}
			getArrayPointer_C(m, l)[0][phi_l[l] % 2] = 0;
			//arrayPointer_C[m][l][0][phi_l % 2] = 0;
			l2 = clonePath(l);
			getArrayPointer_C(m, l2)[0][phi_l[l] % 2] = 1;
			//arrayPointer_C[m][l2][0][phi_l % 2] = 1;
		
			S_m = getArray_S(m, l);
		
			prQueue.emplace(Pairs(*S_m + getD(l, m) + log(omega(phi_l[l])), l));
			prQueue.emplace(Pairs(getD(l, m) + log(omega(phi_l[l])), l2));
			P = P + 2;

			if (phi_l[l] % 2 == 1) {
				recursivelyUpdateC(m, phi_l[l], l);
				recursivelyUpdateC(m, phi_l[l], l2);
			}
			if (l >= theta)
				cout << "ALERT!" << endl;
				phi_l[l]++;			
				if (l2 >= theta)
					cout << "ALERT!" << endl;
				phi_l[l2] = phi_l[l];
			
		}
		if (q[phi_l[l]] >= L) {
			std::set<Pairs>::iterator el = prQueue.begin();
			while(el != prQueue.end()){
				l2 = el->second;
				if (phi_l[l2] < phi_l[l]) {
					prQueue.erase(el++);
					killPath(l2);
					P = P - 1;
					continue;
				}
				++el;
			}
		}
	}
		return out; //not sure
}

void SeqDecoder::check() {
	/*
	cout << "arryReferenceCount content: " << endl;
	for (int lambda = 0; lambda < m + 1; lambda++) {
	for (int s = 0; s < L; s++) {
	cout << "arrayReferenceCount[" << lambda << "][" << s << "]: " << arrayReferenceCount[lambda][s] << endl;
	}
	cout << "--------" << endl;
	}

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

	*/
	cout << "Array Pointer C content: " << endl;
	int lambda2;

	for (lambda2 = 0; lambda2 < m + 1; lambda2++) {
		for (auto&&el : prQueue) {
			cout << "array C[" << el.second << "]" << endl;
			for (int i = 0; i < (int)pow(2.0, m - lambda2); i++) {
				cout << arrayPointer_C[lambda2][pathIndexToArrayIndex[lambda2][el.second]][i][0] << " ";
				cout << arrayPointer_C[lambda2][pathIndexToArrayIndex[lambda2][el.second]][i][1] << " ";
				cout << endl;
			}
			cout << endl;
			cout << "--------" << endl;
		}
	}
}
	//S = new double[L*(m+1)*n]; //S_l_lambda[beta] = S[(lambda*L+l)*n+beta]
	/*
	cout << "Array Pointer S content: " << endl;
	

	for (lambda2 = 0; lambda2 < m + 1; lambda2++) {
		for (auto&&el : prQueue) {
			cout << "array S[" << el.second << "]" << endl;
			for (int i = 0; i < (int)pow(2.0, m - lambda2); i++) {
				cout << S[(lambda2*theta + pathIndexToArrayIndex[lambda2][el.second])*n + i];
				//cout << S[lambda2][pathIndexToArrayIndex[lambda2][el.second]][i] << " ";
				
				cout << endl;
			}
			cout << endl;
			cout << "--------" << endl;
		}
	}


	cout << "Array D content: " << endl;

	for (lambda2 = 0; lambda2 < m + 1; lambda2++) {
		for (auto&&el : prQueue) {
			cout << "array D[" << el.second << "]" << endl;
			cout << D[lambda2*theta + pathIndexToArrayIndex[lambda2][el.second]] << " ";
				//cout << S[lambda2][pathIndexToArrayIndex[lambda2][el.second]][i] << " ";

				cout << endl;
			}
			cout << endl;
			cout << "--------" << endl;
		}
	}

/*
	cout << "Array D content: " << endl;
	for (lambda2 = 0; lambda2 < m + 1; lambda2++) {
		for (int t = 0; t < L; t++) {
			cout << "array D[ " << lambda2 << " ][ " << t << " ]" << endl;
				cout << D[lambda2*L + t] << " ";
				cout << endl;
			cout << endl;
			cout << "--------" << endl;
		}
	}
	
}
*/
