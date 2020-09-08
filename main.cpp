#include <vector>
#include <iostream>
#include <set>
#include <fstream>
#include <functional>
#include "SeqDecoder.h"
#include "decoder.h"
#include "header_func.h"
#include <chrono>
#include <random>
#include <math.h>
#include <queue>
#include "mat.h"
#include "msg.h"

using namespace std;
typedef pair<int, int> P;
int main(int argc, char** argv) {
	
	unsigned int n = 2048;
        unsigned int k = 1024;
	int L = 32;
	int iteration_count = 1;
	double SNR =2;
	long double dispersion;
	long double sigma;
	double E_s = 1;
	set<int> frozen;
//	int* input;
	double h = 1;
	int* random_gen = new int[k];					// random vector
	int error_count = 0;
	vector<int> cw_frozen(n);					//word with frozen symbols
	//long double *ccw_frozen = new long double[n];
	//long double *x = new long double[n];					// codeword
	auto *x = new int[n];
	auto *x2 = new long double[n];	// noisy codeword;
	int flag = 0;							// for error count
	int index = 0;							// for making symbols frozen
	int* out;								//word from decoder
	int theta = 1024;
	int g;
	auto * P_j = new double[n];
    unsigned int** publicKeyG = (unsigned int**)malloc(sizeof(unsigned int*) * k);
    unsigned int** z1 = (unsigned int**)malloc(sizeof(unsigned int*) * 1);
    const unsigned int publicKeyt = 38; //31
    unsigned int** S = (unsigned int**)malloc(sizeof(unsigned int*) * k);
    unsigned int** SInverse = (unsigned int**)malloc(sizeof(unsigned int*) * k);
    unsigned int** P = (unsigned int**)malloc(sizeof(unsigned int*) * n);
    unsigned int** PInverse = (unsigned int**)malloc(sizeof(unsigned int*) * n);
    unsigned int** msg = (unsigned int**)malloc(sizeof(unsigned int*) * 1);
    unsigned int** cPrm = (unsigned int**)malloc(sizeof(unsigned int*) * 1);
    unsigned int** c = (unsigned int**)malloc(sizeof(unsigned int*) * 1);
    unsigned int** G = (unsigned int**)malloc(sizeof(unsigned int*) * k);
    unsigned int** decipher = (unsigned int**)malloc(sizeof(unsigned int*) * 1);
    decipher[0] = (unsigned int*)malloc(sizeof(unsigned int) * k);

	while (SNR <=3) {
		dispersion = E_s * n / (2 * k* pow(10, (SNR / 10)));
		sigma = sqrt(dispersion);
		SNR = SNR + h;
		//cTransform(epsilon, n, n - k, &frozen);		//find frozen symbols
		error_count = 0;
		std::mt19937 generator;
		std::uniform_int_distribution<int> distribution(0, 1);
		meanCount(n, static_cast<double>(dispersion), P_j, n - k, &frozen);
		//meanCount(n, dispersion, P_j);
		SeqDecoder decoder(L, n, frozen, static_cast<double>(dispersion), P_j, theta);
		//decoder dec(L, n, frozen, dispersion);
        for (int i = 0; i < k; i++)
        {
            publicKeyG[i] = (unsigned int*)malloc(sizeof(unsigned int) * n);
            G[i] = (unsigned int*)malloc(sizeof(unsigned int) * n);



            S[i] = (unsigned int*)malloc(sizeof(unsigned int) * k);
            SInverse[i] = (unsigned int*)malloc(sizeof(unsigned int) * k);
        }
        while(1)
        {
            nonSingularMatrix(S, SInverse, k); //generate S

            if (inverseMatrix(S, SInverse, k) != 0) //calcurate S inverse
                break;
        }
		for (g = 0; g < iteration_count; g++) {
			
			flag = 0;
			if (error_count == 20) {
				cout << SNR - h<<" "<< (double)20/g;
				cout << endl;
				break;
			}
			for (int i = 0; i < k; i++) {
				random_gen[i] = distribution(generator);		//random generator
			}
            msg[0] = (unsigned int*)malloc(sizeof(unsigned int) * k);
            for (int l = 0; l <k ; ++l)
                        msg[0][l] = random_gen[l];


			index = 0;
			for (int i = 0; i < n; i++) {
				if (frozen.find(i) != frozen.end())
					cw_frozen[i] = 0;
				else {
					cw_frozen[i] = random_gen[index];
					index++;
				}
			}

			encode(cw_frozen);
            for (int i = 0; i < n; i++)
                x[i] = -2 * cw_frozen[i] + 1;


            for(int j = 0; j < k; j++)
            {
                for(int i = 0; i < n; i++)
                {
                    G[j][i] = static_cast<unsigned int>(cw_frozen[j * n + i]);
                }
            }

            for (int i = 0; i < n; i++)
            {
                P[i] = (unsigned int*)malloc(sizeof(unsigned int) * n);
                PInverse[i] = (unsigned int*)malloc(sizeof(unsigned int) * n);
            }
            permutationMatrix(P, PInverse, n); //generate P
            inverseMatrix(P, PInverse, n); //calcurate P inverse




            z1[0] = (unsigned int*)malloc(sizeof(unsigned int) * n);
            for (int i = 0; i < n; i++)
            {
                z1[0][i] = 0;
            }

            generatePublicKey(k, n, SInverse, G, P, publicKeyG, publicKeyt, z1[0]);
            cPrm[0] = (unsigned int*)malloc(sizeof(unsigned int) * n);
            c[0] = (unsigned int*)malloc(sizeof(unsigned int) * n);
            msgEncrypt(k, msg, publicKeyG, n, z1, cPrm, c); //encrypt meassage
            msgDecrypt(n, c, PInverse, k, SInverse, decipher, G, 5, 12); //decrypt recived meassage

			for (int j = 0; j <n ; ++j) {
				x2[j] =(long double)decipher[0][j];
			}
			out = decoder.decode(x2);	// decoded word;
			//out = dec.decode(x2);
			for (int i = 0; i < n; i++) {
				if (out[i] != cw_frozen[i])
					flag = 1;
			}
			if (flag == 1) {
				error_count++;
			}


		}
        cout << "errors number : "<< error_count;
		if (g == iteration_count) {
			cout << SNR - h;
			cout << endl;
		}
		frozen.clear();
	}
    
        free(decipher[0]);
        free(decipher);
        free(cPrm[0]);
        free(cPrm);
        free(c[0]);
        free(c);
        free(msg[0]);
        free(msg);
	return 0;
}

