#ifndef MSG_H
#define MSG_H

void msgEncrypt(const unsigned int& k, unsigned int **msg, unsigned int **publicKey, const unsigned int& n, unsigned int **error, unsigned int **prime, unsigned int **c);
void msgDecrypt(const unsigned int& n, unsigned int **c, unsigned int **PInverse, const unsigned int& k, unsigned int **SInverse, unsigned int **msg, unsigned int **GRM, const unsigned int& r, const unsigned int& m);


class msg {
public:
    msg();
    msg(const msg& orig);
    virtual ~msg();
private:

};

#endif /* MSG_H */

