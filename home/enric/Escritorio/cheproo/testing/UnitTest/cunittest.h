#ifndef CUNITTEST_H
#define CUNITTEST_H

#include <QString>

#include "base.h"


extern void gInit_cheproo(void);
extern void gEnd(void);

#define MSG_UNSUCCESFULL_TEST " has NOT passed test  "
#define MSG_UNSUCCESFULL_TEST_STEP " step "

class CUnitTest 
{
public:
    static void RunBasicUT(const char* aTestName, int aUnitTest(void), int testNumber);

    static void RunExceptionUT(const char* aTestName, 
                               int aUnitTest(void), 
                               int testNumber, 
                               QString aExceptionMessage);

};

#endif
