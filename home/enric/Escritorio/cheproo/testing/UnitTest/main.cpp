#include <iostream>
#include <vector>
using namespace std;

#include <QCoreApplication>
#include <QString>
#include <QStringList>
#include <QPair>
#include <QDir>

#include "generalexception.h"

using namespace std;

extern void UTCAqueous(void);
extern void UTCAqueousDebyeHuckel(void);
extern void UTCAqueousIdeal(void);
extern void UTCCheprooPlusPlus(void);
extern void UTCLocalChemicalSystem(void);
extern void UTCLocalChemicalSystemSaaltink(void);

extern int cheprooMain(QStringList& cheprooParams);

//extern void getSystestInputs(vector<STestSpec>&);


//int ExecCheprooMain(const QString &path, 
//                   const QString &inputFilename, 
//                   vector<CFileComp > aComparefiles, 
//                   const QString &execPath)
//				   
//{
//    QString fullFileName = path + inputFilename;
//    fullFileName = QDir(fullFileName).absolutePath();
//
//
//	cout << "Running: '" << qPrintable(inputFilename) << "'" << endl;
//
//	int i = main(argc, );
//
//    if ( i == 0 ) 
//	{
//       bool allpassed = true;
//       for (int j = 0; j < (int) aComparefiles.size(); j++)
//       {
//            allpassed = allpassed && aComparefiles[j].Compare();
//       }
//       if (allpassed)
//       {
//            cout << "Finished" << endl << endl; ;
//       }
//		else
//		{
//			cout << "Failed: '" << qPrintable(fullFileName) << "'" << " with error:" << i << endl << endl;
//		}
//    }
//
//
//	return i;
//}


int main ( int argc, char *argv[] )
{
    //QCoreApplication app(argc, argv);

    //-----------------------------------------------------------
	cout << "Starting Cheproo Unit Tests...\n"

         << "Running:\n";



    try
    {
        // CCheprooPlusPlus Tests
        //-----------------------------------------
        UTCCheprooPlusPlus();


        // CLocalChemicalSystem Tests
        //-----------------------------------------
        UTCLocalChemicalSystem();
        UTCLocalChemicalSystemSaaltink();

        // CPhase Tests
        //-----------------------------------------
        UTCAqueous();
	    UTCAqueousDebyeHuckel();
        UTCAqueousIdeal();

        //-----------------------------------------
    }
    catch (GeneralException& ex)
	{
		cout << ex.MessageString() << endl;
		return 1;
	}
	
    cout << "...Ended Cheproo Unit Tests\n";
    //-----------------------------------------------------------
    //sys_tests:
//    //-----------------------------------------------------------
//    cout << "\nStarting Cheproo System Tests...\n";
//    vector<STestSpec> tests;
//	getSystestInputs(tests);
//    foreach(STestSpec test, tests)
//    {
//        ExecProostMain( test.mPath, test.mInputfile, test.mComparefiles, test.mExecdir  );
//    }
//
//	cout << "...Ended Proost System Test\n";
    //-----------------------------------------------------------

	cin.get();
	return 0;
}


