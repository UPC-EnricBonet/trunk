#include <iostream>
#include <QFile>

#include "ccheprooplusplus.h"
#include "cheproomain.h"

#include "xmlconstants.h"
#include "xmlqt.h"
#include "generalexception.h"
#include "objectmanager.h"
#include "classnameconstants.h"
//#include "proostglobals.h"


#define CHEPROO_VERSION "0.0.1"

extern void gInit_cheproo(void);
extern void gEnd(void);


static const char* START_MESSAGE = QT_TRANSLATE_NOOP("main", "Started Cheproo version ");
static const char* END_MESSAGE = QT_TRANSLATE_NOOP("main", "Ended Cheproo");
static const char* ERROR_FILE_DOES_NOT_EXIST = QT_TRANSLATE_NOOP("main", "Input file does not exist: ");
static const char* INPUT_FILENAME_MSG = QT_TRANSLATE_NOOP("main", "Input filename: ");
static const char* PROB_TAG_MUST_EXIST = QT_TRANSLATE_NOOP("main", "CheprooPlusPlus input file must have a <problem> tag" );
static const char* CANNOT_SOLVE_SIM = QT_TRANSLATE_NOOP("main", "Cannot solve simulation");

using namespace std;

/*!
\brief This function starts the problem \n
The reading of the problem from XML is initiated in this function.
If everything ran succesfully, 0 is returned.
*/

int cheprooMain ( QStringList& cheprooParams )
{
	try{
		int retVal = 0; // Return 0 if program ends succesfully
		bool showOutput = true;

		// The number of parameters is 2
		// 1. the executable name
		// 2. (optional) the filename of the xml file that contains the problem to be solved
		// If no xml filename is given, the default "cheprooplusplus.xml" is used
		int numOfParams = cheprooParams.size();
		QString inputFile;
		if( numOfParams < 2 )
		{
			inputFile = "cheprooplusplus.xml";
		}
		else if( numOfParams == 2)
		{
			inputFile = cheprooParams[1];
		}
		else
		{
			inputFile = cheprooParams[1];
			if(cheprooParams[2] == "-silent")
			{
				showOutput = false;
			}
		}

		if(showOutput)
		{
			cout << PrTr("main", START_MESSAGE) << CHEPROO_VERSION << endl;
			cout << "--------------------------------" << endl;
		}
	    try
	    {		
		    // Check if the xml file exists, if not show an error and return 1
		    QFile tstFile;

		    tstFile.setFileName(inputFile);

		    if( ! tstFile.exists() )
		    {
			    cout << PrTr("main", ERROR_FILE_DOES_NOT_EXIST) <<  qPrintable(inputFile) << endl;
			    return 1;
		    }
            else
		    {
			    cout << PrTr("main", INPUT_FILENAME_MSG) <<  qPrintable(inputFile) << endl;
		    }

		    // Open inputFile and parse the xml. The result is stored in doc
		    QDomDocument doc("mydocument");
		    XmlQt::xmlParseFile(inputFile, &doc);

		    // Instantiate the object manager and register all the lists and classes
		    gInit_cheproo();

            // Get a pointer to the object manager instance
		    ObjectManager* om = ObjectManager::instance();

		    // Fill all the lists with empty pointers. 
		    // These pointers will be allocated by the specific Read() functions of each class
		    bool mustDelete = true;
            bool typeAttributeMustExist = false;
            
		    om->FillAllListsFromXML(doc.documentElement(), typeAttributeMustExist, mustDelete);

    	    // Retrieve the CheprooPlusPlus element in the xml and the name of the problem
		    QDomElement nodeCheprooPlusPlus = XmlQt::getFirstElemByTag(doc.documentElement(),XML_TAG_CHEPROOPLUSPLUS);
            QString aProblemName;
            XmlQt::getXMLAttribute(nodeCheprooPlusPlus, XML_ATTR_NAME, aProblemName); 

		    bool CheprooPlusPlusExists = ! nodeCheprooPlusPlus.isNull();
		    if( CheprooPlusPlusExists ) 
		    {
			    // Create a Problem object and read its properties from the xml
                CCheprooPlusPlus* cheprooPlusPlusPointer = dynamic_cast<CCheprooPlusPlus*>(om->newInstanceOfClassName(CLASS_NAME_CHEPROOPLUSPLUS, aProblemName));

                // Read and initialize the problem
                cheprooPlusPlusPointer->ReadAndInitialize(nodeCheprooPlusPlus);

                // Retrieve the name of the function
                string aFunctionName = cheprooPlusPlusPointer->mFunctionName;


                bool realizedSimulation = false;

                if(aFunctionName == "speciateWaters")
                {
                    bool realizedSimulation = cheprooPlusPlusPointer->SpeciateWaters();
                }
     
                if(aFunctionName == "mixWaters")
                {
                    bool realizedSimulation = cheprooPlusPlusPointer->MixWaters();
                }
                else if(aFunctionName == "evaporateWaters")
                {
                    bool realizedSimulation = cheprooPlusPlusPointer->EvaporateWaters();
                }
                else
                {
                    // do nothing
                }

			    //if( ! realizedSimulation )
			    //{
				   // cout << PrTr("main", CANNOT_SOLVE_SIM) << endl;
				   // retVal = 1;
			    //}
		    }

		    else
		    {
			    // If no problem tag exist, return an error
			    cout << PrTr("main", PROB_TAG_MUST_EXIST) << endl;
			    retVal = 1;
    		
            }
        }

   		    catch (GeneralException& ex)
		    {
			    cout << qPrintable(ex.Message()) << endl;

                // Delete all the lists
                ObjectManager* om = ObjectManager::instance();
		        om->EraseAllLists();
                return 1;
            }
        }
	    catch (GeneralException& ex)
	    {
		    cout << qPrintable(ex.Message()) << endl;

            // Delete all the lists
            ObjectManager* om = ObjectManager::instance();
	        om->EraseAllLists();
            return 2;
        }

		// Unregister all the classes
		gEnd();
        
		cout << "Ended CheprooPlusPlus" << endl;

        return(0);
        }



