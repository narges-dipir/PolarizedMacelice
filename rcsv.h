#ifndef RCSV_H
#define RCSV_H

double* rcsvRight(int  r, int  m, double *orgMsg, int orgMsgLength, double *result, int resultLength);
double* rcsvLeft(int  r, int  m, double *orgMsg, int orgMsgLength, double *v, int vLength, double *result, int resultLength);
void rcsvDecisionZero(double *orgMsg, int orgMsgLength, double *result, int resultLength);
void rcsvDecisionEqual(double *orgMsg, int orgMsgLength, double *result, int resultLength);
void rcsvDecodingHard(int  r, int m, unsigned int *recvCodeword, int recvCodewordLength);

class rcsv {
public:
    rcsv();
    rcsv(const rcsv& orig);
    virtual ~rcsv();
private:

};

#endif /* RCSV_H */

