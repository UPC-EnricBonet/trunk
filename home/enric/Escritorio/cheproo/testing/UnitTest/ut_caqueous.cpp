// Unit Testing CAqueous

#include "cunittest.h"
#include "xmlconstants.h"
#include "base.h"
#include "xmlqt.h"
#include "objectlist.h"
#include "objectmanager.h"
#include "valarraytools.h"
#include "classnameconstants.h"
#include "ccheprooplusplus.h"
#include "cglobalchemicalsystem.h"

#include "caqueous.h"
#include "cphase.h"

#include <vector>
#include <QFile>

static QString baseDir = "UnitTest/";

// Method to evaluate the ionic strength
int UTCAqueous_1()
{
    QString inputFile = baseDir + "in/UTCheprooInput.xml";
    QString refFile = baseDir + "ref/UTCAqueous1.xml";

    QDomDocument doc("mydocument");
	XmlQt::xmlParseFile(inputFile, &doc);


    // Get a pointer to the object manager instance
    ObjectManager* om = ObjectManager::instance();

    // Fill all the lists with empty pointers. 
    bool mustDelete = true;
    bool typeAttributeMustExist = false;

    om->FillAllListsFromXML(doc.documentElement(), typeAttributeMustExist, mustDelete);

	// Retrieve the CheprooPlusPlus element in the xml and the name of the problem
    QDomElement nodeCheprooPlusPlus = XmlQt::getFirstElemByTag(doc.documentElement(),XML_TAG_CHEPROOPLUSPLUS);
    QString aProblemName;
    XmlQt::getXMLAttribute(nodeCheprooPlusPlus, XML_ATTR_NAME, aProblemName); 

    bool CheprooPlusPlusExists = ! nodeCheprooPlusPlus.isNull();

    // Declare vector to comprare the results
    valarray<valarray<double> > ionicStrArray(1);
    ionicStrArray[0].resize(1,0);

    if( CheprooPlusPlusExists ) 
    {
	    // Create a Problem object and read its properties from the xml
        CCheprooPlusPlus* cheprooPlusPlusPointer = dynamic_cast<CCheprooPlusPlus*>(om->newInstanceOfClassName(CLASS_NAME_CHEPROOPLUSPLUS, aProblemName));

        // Read and initialize the problem
        cheprooPlusPlusPointer->ReadAndInitialize(nodeCheprooPlusPlus, false);
        
        // Pointer to LocalChemicalSystem
        CLocalChemicalSystem* myLocalPointer = cheprooPlusPlusPointer->mGlobalChemSysPointer->mLocChemSysVector[0];

        // Call the function to speciate initial waters
        myLocalPointer->SpeciateInitialWater(cheprooPlusPlusPointer->mGlobalChemSysPointer->mLocChemSysVector[0]->mRelChemicalCompositions[0]);

        // Declare variable to compare: Ionic Strength
        double ionicStrength = 0.0;

        // Check the method
        CPhase* myAqPhase = cheprooPlusPlusPointer->mGlobalChemSysPointer->mAllPhases["aqueous1"];
        myAqPhase->ComputeIonicStrength(myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration, ionicStrength,
                                        myLocalPointer->mSpecies[eAqueous1nc], myLocalPointer->mSpecies[eAqueous2],
                                        myLocalPointer->mSpeciesIndices);

        ionicStrArray[0][0] = ionicStrength;
        
    }

    //Obtains from sRefBal the content of reference file to compare for Balance
	QFile RefFileBal(refFile);
	RefFileBal.open(QIODevice::ReadOnly | QIODevice::Text);
	QString sRefBal = RefFileBal.readAll();
    RefFileBal.close();
    valarray< valarray< double > > ref;
    ValarrayTools::fromString(sRefBal, ref);

    
	if ( !ValarrayTools::IsEqual(ionicStrArray,ref, 1e-12)) return 1;
    
	return 0;
}

// Method to evaluate the charge balance
int UTCAqueous_2()
{
    QString inputFile = baseDir + "in/UTCheprooInput.xml";
    QString refFile = baseDir + "ref/UTCAqueous2.xml";

    QDomDocument doc("mydocument");
	XmlQt::xmlParseFile(inputFile, &doc);


    // Get a pointer to the object manager instance
    ObjectManager* om = ObjectManager::instance();

    // Fill all the lists with empty pointers. 
    bool mustDelete = true;
    bool typeAttributeMustExist = false;

    om->FillAllListsFromXML(doc.documentElement(), typeAttributeMustExist, mustDelete);

	// Retrieve the CheprooPlusPlus element in the xml and the name of the problem
    QDomElement nodeCheprooPlusPlus = XmlQt::getFirstElemByTag(doc.documentElement(),XML_TAG_CHEPROOPLUSPLUS);
    QString aProblemName;
    XmlQt::getXMLAttribute(nodeCheprooPlusPlus, XML_ATTR_NAME, aProblemName); 

    bool CheprooPlusPlusExists = ! nodeCheprooPlusPlus.isNull();

    // Declare vector to comprare the results
    valarray<valarray<double> > chargeBalArray(1);
    chargeBalArray[0].resize(1,0);

    if( CheprooPlusPlusExists ) 
    {
	    // Create a Problem object and read its properties from the xml
        CCheprooPlusPlus* cheprooPlusPlusPointer = dynamic_cast<CCheprooPlusPlus*>(om->newInstanceOfClassName(CLASS_NAME_CHEPROOPLUSPLUS, aProblemName));

        // Read and initialize the problem
        cheprooPlusPlusPointer->ReadAndInitialize(nodeCheprooPlusPlus, false);

        // Pointer to GlobalChemicalSystem
        CLocalChemicalSystem* myLocalPointer = cheprooPlusPlusPointer->mGlobalChemSysPointer->mLocChemSysVector[0];

        // Call the function to speciate initial waters
        myLocalPointer->SpeciateInitialWater(cheprooPlusPlusPointer->mGlobalChemSysPointer->mLocChemSysVector[0]->mRelChemicalCompositions[0]);

        // Declare variable to compare: Ionic Strength
        double chargeBalance = 0.0;

        // Check the method
        CPhase* myAqPhase = cheprooPlusPlusPointer->mGlobalChemSysPointer->mAllPhases["aqueous1"];
        myAqPhase->ComputeChargeBalance(myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration, chargeBalance,
                                        myLocalPointer->mSpecies[eAqueous1nc], myLocalPointer->mSpecies[eAqueous2],
                                        myLocalPointer->mSpeciesIndices);

        chargeBalArray[0][0] = chargeBalance;
        
    }

    //Obtains from sRefBal the content of reference file to compare for Balance
	QFile RefFileBal(refFile);
	RefFileBal.open(QIODevice::ReadOnly | QIODevice::Text);
	QString sRefBal = RefFileBal.readAll();
    RefFileBal.close();
    valarray< valarray< double > > ref;
    ValarrayTools::fromString(sRefBal, ref);

    
	if ( !ValarrayTools::IsEqual(chargeBalArray,ref, 1e-12)) return 1;
    
	return 0;
}


void UTCAqueous()
{
    const char* testName = "UTCAqueous";

    CUnitTest::RunBasicUT(testName, UTCAqueous_1, 1);
    CUnitTest::RunBasicUT(testName, UTCAqueous_2, 2);

    //CUnitTest::RunExceptionUT(testName, UTBase_1, 1, ERROR_UNINITIALIZED_XML_NODE);
}
