// Unit Testing CLocalChemicalSystem

#include "cunittest.h"
#include "xmlconstants.h"
#include "base.h"
#include "xmlqt.h"
#include "objectlist.h"

#include "objectmanager.h"
#include "valarraytools.h"
#include "classnameconstants.h"
#include "ccheprooplusplus.h"

#include <vector>
#include <QFile>
#include <math.h>

static QString baseDir = "UnitTest/";


// Method to verify the sequence of calls to the speciation method (the results have already been verified by means of other system tests)
int UTCCheprooPlusPlus_1()
{
    QString inputFile = baseDir + "in/UTCheprooInput1.xml";
    QString refFile = baseDir + "ref/UTCLocalChemicalSystem1.xml";

    QDomDocument doc("mydocument");
	XmlQt::xmlParseFile(inputFile, &doc);


    // Get a pointer to the object manager instance
    ObjectManager* om = ObjectManager::instance();

    // Create a CheprooPlusPlus object
    CCheprooPlusPlus* cheprooPointer = dynamic_cast<CCheprooPlusPlus*>(om->newInstanceOfClassName(CLASS_NAME_CHEPROOPLUSPLUS, "CheprooPlusPlus"));

    // Create water objects
    CChemicalComposition* chemComp1 = dynamic_cast<CChemicalComposition*> (om->newInstanceOfClassName(CLASS_NAME_CHEMICALCOMPOSITION, "initial_water_1"));
    CChemicalComposition* chemComp2 = dynamic_cast<CChemicalComposition*> (om->newInstanceOfClassName(CLASS_NAME_CHEMICALCOMPOSITION, "boundary_water_1"));

    vector<CChemicalComposition*> inWaters;
    vector<CChemicalComposition*> boundWaters;
    inWaters.push_back(chemComp1);
    boundWaters.push_back(chemComp2);

    // Try the function
    //cheprooPointer->ReadAndInitialize("cheprooplusplusForProost.xml", inWaters, boundWaters);

    // If everything was good then return 0
	return 0;
}

// Method to evaluate the mixWaters method with CLocalChemicalSystemPrimSpecies
int UTCCheprooPlusPlus_2()
{
    QString inputFile = baseDir + "in/UTCLocalChemicalSystem/UTCLocalChemicalSystemMixWaters.xml";
    QString refFile = baseDir + "ref/UTCLocalChemicalSystem/UTCLocalChemicalSystemMixWaters.xml";

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
    valarray<valarray<double> > concVector;

    if( CheprooPlusPlusExists ) 
    {
	    // Create a Problem object and read its properties from the xml
        CCheprooPlusPlus* cheprooPlusPlusPointer = dynamic_cast<CCheprooPlusPlus*>(om->newInstanceOfClassName(CLASS_NAME_CHEPROOPLUSPLUS, aProblemName));

        // Read and initialize the problem
        cheprooPlusPlusPointer->ReadAndInitialize(nodeCheprooPlusPlus, false);

        // Pointer to LocalChemicalSystem
        CLocalChemicalSystem* myLocalPointer = cheprooPlusPlusPointer->mGlobalChemSysPointer->mLocChemSysVector[0];

        // Speciate initial waters
        myLocalPointer->SpeciateInitialWater(cheprooPlusPlusPointer->mGlobalChemSysPointer->mLocChemSysVector[0]->mRelChemicalCompositions[0]);
        myLocalPointer->SpeciateInitialWater(cheprooPlusPlusPointer->mGlobalChemSysPointer->mLocChemSysVector[0]->mRelChemicalCompositions[1]);

        // Declare variable to compare: Concentrations
        vector<double> concentrations;

        // Check the method
        cheprooPlusPlusPointer->MixWaters();

        // Verify concentrations
        concVector.resize(myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration.size());
        for(int i=0; i!=concVector.size(); i++)
        {
            concVector[i].resize(1);
            concVector[i][0] = myLocalPointer->mRelChemicalCompositions.back()->mConcentration[i];
        }
        
    }

    //Obtains from sRefBal the content of reference file to compare for Balance
	QFile RefFileBal(refFile);
	RefFileBal.open(QIODevice::ReadOnly | QIODevice::Text);
	QString sRefBal = RefFileBal.readAll();
    RefFileBal.close();
    valarray< valarray< double > > ref;
    ValarrayTools::fromString(sRefBal, ref);

    
	if ( !ValarrayTools::IsEqual(concVector,ref, 1e-12)) return 1;
    
	return 0;
}

// Method to evaluate the mixChemicalComposition method with CLocalChemicalSystemPrimSpecies (gypsum example, it has to give the same results as MixWaters)
int UTCCheprooPlusPlus_3()
{
    QString inputFile = baseDir + "in/UTCLocalChemicalSystem/UTCLocalChemicalSystemMixWaters.xml";
    QString refFile = baseDir + "ref/UTCLocalChemicalSystem/UTCLocalChemicalSystemMixWaters.xml";

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
    valarray<valarray<double> > concVector;

    if( CheprooPlusPlusExists ) 
    {
	    // Create a Problem object and read its properties from the xml
        CCheprooPlusPlus* cheprooPlusPlusPointer = dynamic_cast<CCheprooPlusPlus*>(om->newInstanceOfClassName(CLASS_NAME_CHEPROOPLUSPLUS, aProblemName));

        // Read and initialize the problem
        cheprooPlusPlusPointer->ReadAndInitialize(nodeCheprooPlusPlus, false);

        // Pointer to LocalChemicalSystem
        CLocalChemicalSystem* myLocalPointer = cheprooPlusPlusPointer->mGlobalChemSysPointer->mLocChemSysVector[0];

        // Speciate initial waters
        myLocalPointer->SpeciateInitialWater(cheprooPlusPlusPointer->mGlobalChemSysPointer->mLocChemSysVector[0]->mRelChemicalCompositions[0]);
        myLocalPointer->SpeciateInitialWater(cheprooPlusPlusPointer->mGlobalChemSysPointer->mLocChemSysVector[0]->mRelChemicalCompositions[1]);

        // Declare variable to compare: Concentrations
        vector<double> concentrations;

        // Check the method
        cheprooPlusPlusPointer->MixChemicalCompositions();

        // Verify concentrations
        concVector.resize(myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration.size());
        for(int i=0; i!=concVector.size(); i++)
        {
            concVector[i].resize(1);
            concVector[i][0] = myLocalPointer->mRelChemicalCompositions.back()->mConcentration[i];
        }
        
    }

    //Obtains from sRefBal the content of reference file to compare for Balance
	QFile RefFileBal(refFile);
	RefFileBal.open(QIODevice::ReadOnly | QIODevice::Text);
	QString sRefBal = RefFileBal.readAll();
    RefFileBal.close();
    valarray< valarray< double > > ref;
    ValarrayTools::fromString(sRefBal, ref);

    
	if ( !ValarrayTools::IsEqual(concVector,ref, 1e-12)) return 1;
    
	return 0;
}


void UTCCheprooPlusPlus()
{
    const char* testName = "UTCCheprooPlusPlus";

    //CUnitTest::RunBasicUT(testName, UTCCheprooPlusPlus_1, 1);
    CUnitTest::RunBasicUT(testName, UTCCheprooPlusPlus_2, 2);
    CUnitTest::RunBasicUT(testName, UTCCheprooPlusPlus_3, 3);
}
