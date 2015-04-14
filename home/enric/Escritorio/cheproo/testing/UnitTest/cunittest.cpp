#include "cunittest.h"

#include "objectmanager.h"
#include <iostream>
using namespace std;

void 
CUnitTest::RunBasicUT(const char* aTestName, 
                           int aUnitTest(void), 
                           int testNumber)
{
    gInit_cheproo();
	gTimeStamp = 0;
        
    int returnValue = -1;

    try
	{
        returnValue = aUnitTest();
    }
    catch (GeneralException& ex)
	{
		cout << ex.MessageString() << endl;
        returnValue = 100;
	}
	catch(...)
	{
        returnValue = 200;
	}
	gEnd();

    // if an error was generated, then the object manager may not have been cleaned properly. 
    ObjectManager* om = ObjectManager::instance();
    om->EraseAllLists();

    cout <<            aTestName
         << endl;


    if(returnValue != 0)
    {
        cout << aTestName 
             << MSG_UNSUCCESFULL_TEST 
             << testNumber 
             << MSG_UNSUCCESFULL_TEST_STEP 
             << returnValue
             << endl;
    }
}

void 
CUnitTest::RunExceptionUT(const char* aTestName, 
                               int aUnitTest(void), 
                               int testNumber,
                               QString aExceptionMessage)
{
    gInit_cheproo();
	gTimeStamp = 0;
        
    int returnValue = -1;

    try
	{
        aUnitTest();
    }
    catch (GeneralException& ex)
	{
        QString thrownMessage = ex.Message();
        if(thrownMessage.endsWith("."))
        {
            thrownMessage.chop(1);
            aExceptionMessage.chop(1);
        }
        thrownMessage = thrownMessage.section('.', -1);
        if(thrownMessage == aExceptionMessage)
        {
            returnValue = 0;
            ObjectManager* om = ObjectManager::instance();

        }
        else
        {
            cout << ex.MessageString() << endl;
            returnValue = 100;
        }
    }	
	catch(...)
	{
        returnValue = 200;
	}
	gEnd();

    if(returnValue != 0)
    {
        cout << aTestName 
             << MSG_UNSUCCESFULL_TEST 
             << testNumber 
             << MSG_UNSUCCESFULL_TEST_STEP 
             << returnValue
             << endl;
    }
}
