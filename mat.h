#ifndef MAT_H
#define MAT_H

int inverseMatrix(unsigned int **matrix, unsigned int **inverse, const unsigned int& size);
void nonSingularMatrix(unsigned int **matrix, unsigned int **identity, const unsigned int& size);
void permutationMatrix(unsigned int **matrix, unsigned int **identity, const unsigned int& size);
void mtrixMultiplication(const unsigned int& pRow, const unsigned int& pCol, unsigned int **prev, const unsigned int& nRow, const unsigned int& nCol, unsigned int **next, unsigned int **result);
void errorAdd(const unsigned int& size, unsigned int *original, unsigned int *error, unsigned int *result);
void errorGenerator(const unsigned int& weight, const unsigned int& size, unsigned int *error);
void generatePublicKey(const unsigned int& k, const unsigned int& n, unsigned int **S, unsigned int**G, unsigned int **P, unsigned int **key,const unsigned int& t, unsigned int *z);
void generateMatrixMultiple(const int& n, const int& m, const int& lim, unsigned int **matrix, unsigned int& idx, int *row, int rowIdx, int& rowValue);
void generateMatrix(const int& r, const int& m, const int& n, unsigned int**matrix);


class mat {
public:
    mat();
    mat(const mat& orig);
    virtual ~mat();
private:

};

#endif /* MAT_H */

