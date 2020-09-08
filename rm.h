
#ifndef RM_H
#define RM_H
void rmRUpdate(unsigned int& start, unsigned int& end, unsigned int **recovery, unsigned int **G, unsigned int **codeword, const unsigned int& n);
unsigned int rmDecision(unsigned int *array, unsigned int arrayLength);
void rmIdxCalculator(unsigned int *array, const int& arrayLen, unsigned int& arrayPointer, unsigned int *result, unsigned int& resultIdx, unsigned int& resultValue);
void rmSCalculate(unsigned int *row, const int& rowLength, unsigned int *S);
void rmECalculate(unsigned int *row, const int& rowLength, const int& m, unsigned int *E);
void rmSCCalculate(unsigned int *row, const int& rowLength, unsigned int *SC);
unsigned int rmACalculate(unsigned int& SLength, unsigned int *S, const unsigned int& SCLength, unsigned int *SC, unsigned int **codeword);
void rmDecoding(const unsigned int& n, const unsigned int& m, const unsigned int& rowLength, unsigned int* row, unsigned int rowIdx, unsigned int& rowValue, unsigned int **codeword, unsigned int **recovery, unsigned int& recoveryIdx);
void rmDecoder(const unsigned int& r, const unsigned int& m, unsigned int **GRM, const unsigned int& n, unsigned int **codeword, unsigned int **recovery);



class rm {
public:
    rm();
    rm(const rm& orig);
    virtual ~rm();
private:

};

#endif /* RM_H */

