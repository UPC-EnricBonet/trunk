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
#include "cglobalchemicalsystem.h"

#include "clocalchemicalsystem.h"

#include <vector>
#include <QFile>
#include <math.h>

static QString baseDir = "UnitTest/";

// Method to evaluate the speciation method to initialize with CLocalChemicalSystemConstPrimary
int UTCLocalChemicalSystem_1()
{
    QString inputFile = baseDir + "in/UTCLocalChemicalSystem/UTCLocalChemicalSystem1.xml";
    QString refFile = baseDir + "ref/UTCLocalChemicalSystem/UTCLocalChemicalSystem1.xml";

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

        // Call the function
        myLocalPointer->SpeciateInitialWater(cheprooPlusPlusPointer->mGlobalChemSysPointer->mLocChemSysVector[0]->mRelChemicalCompositions[0]);

        // Declare variable to compare: concentrations
        vector<double> concentrations;

        concVector.resize(myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration.size());
        for(int i=0; i!=concVector.size(); i++)
        {
            concVector[i].resize(1);
            concVector[i][0] = myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration[i];
        }
        
    }

    //Obtains from sRefBal the content of reference file to compare for Balance
	QFile RefFileBal(refFile);
	RefFileBal.open(QIODevice::ReadOnly | QIODevice::Text);
	QString sRefBal = RefFileBal.readAll();
    RefFileBal.close();
    valarray< valarray< double > > ref;
    ValarrayTools::fromString(sRefBal, ref);

    
	if ( !ValarrayTools::IsEqual(concVector,ref, 1e-9)) return 1;
    
	return 0;
}

// Method to evaluate the speciation method with charge balance costrain
int UTCLocalChemicalSystem_2()
{
    QString inputFile = baseDir + "in/UTCLocalChemicalSystem/UTCLocalChemicalSystem2.xml";
    QString refFile = baseDir + "ref/UTCLocalChemicalSystem/UTCLocalChemicalSystem2.xml";

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

        // Call the function
        myLocalPointer->SpeciateInitialWater(cheprooPlusPlusPointer->mGlobalChemSysPointer->mLocChemSysVector[0]->mRelChemicalCompositions[0]);

        // Declare variable to compare: Activity coefficients
        vector<double> concentrations;

        concVector.resize(myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration.size());
        for(int i=0; i!=concVector.size(); i++)
        {
            concVector[i].resize(1);
            concVector[i][0] = myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration[i];
        }
        
    }

    //Obtains from sRefBal the content of reference file to compare for Balance
	QFile RefFileBal(refFile);
	RefFileBal.open(QIODevice::ReadOnly | QIODevice::Text);
	QString sRefBal = RefFileBal.readAll();
    RefFileBal.close();
    valarray< valarray< double > > ref;
    ValarrayTools::fromString(sRefBal, ref);

    
	if ( !ValarrayTools::IsEqual(concVector,ref, 1e-10)) return 1;
    
	return 0;
}

// Method to evaluate the speciation method with CLocalChemicalSystemSaaltink. All homogeneous aqueous reactions
int UTCLocalChemicalSystem_3()
{
    QString inputFile = baseDir + "in/UTCLocalChemicalSystem/UTCLocalChemicalSystemSaaltink2.xml";
    QString refFile = baseDir + "ref/UTCLocalChemicalSystem/UTCLocalChemicalSystemSaaltink2.xml";

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

        // Call the function
        myLocalPointer->SpeciateInitialWater(cheprooPlusPlusPointer->mGlobalChemSysPointer->mLocChemSysVector[0]->mRelChemicalCompositions[0]);

        // Declare variable to compare: Concentrations
        vector<double> concentrations;

        concVector.resize(myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration.size());
        for(int i=0; i!=concVector.size(); i++)
        {
            concVector[i].resize(1);
            concVector[i][0] = myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration[i];
        }
        
    }

    //Obtains from sRefBal the content of reference file to compare for Balance
	QFile RefFileBal(refFile);
	RefFileBal.open(QIODevice::ReadOnly | QIODevice::Text);
	QString sRefBal = RefFileBal.readAll();
    RefFileBal.close();
    valarray< valarray< double > > ref;
    ValarrayTools::fromString(sRefBal, ref);

    
	if ( !ValarrayTools::IsEqual(concVector,ref, 1e-10)) return 1;
    
	return 0;
}

// Method to evaluate the speciation method with CLocalChemicalSystemSaaltink. All homogeneous aqueous reactions
int UTCLocalChemicalSystem_4()
{
    QString inputFile = baseDir + "in/UTCLocalChemicalSystem/UTCLocalChemicalSystemSaaltink3.xml";
    QString refFile = baseDir + "ref/UTCLocalChemicalSystem/UTCLocalChemicalSystemSaaltink3.xml";

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

        // Call the function
        myLocalPointer->SpeciateInitialWater(cheprooPlusPlusPointer->mGlobalChemSysPointer->mLocChemSysVector[0]->mRelChemicalCompositions[0]);

        // Declare variable to compare: Concentrations
        vector<double> concentrations;

        concVector.resize(myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration.size());
        for(int i=0; i!=concVector.size(); i++)
        {
            concVector[i].resize(1);
            concVector[i][0] = myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration[i];
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

// Method to evaluate the RISA method with CLocalChemicalSystemConstPrimary. Gypsum example
int UTCLocalChemicalSystem_5()
{
    QString inputFile = baseDir + "in/RISA/gypsum.xml";
    QString refFile = baseDir + "ref/RISA/gypsum.xml";

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

        // Call the function
        myLocalPointer->RISA(cheprooPlusPlusPointer->mGlobalChemSysPointer->mLocChemSysVector[0]->mRelChemicalCompositions[0]);

        // Declare variable to compare: Concentrations
        vector<double> concentrations;

        concVector.resize(myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration.size());
        for(int i=0; i!=concVector.size(); i++)
        {
            concVector[i].resize(1);
            concVector[i][0] = myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration[i];
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

// Method to evaluate the RISA method with CLocalChemicalSystemConstPrimary. Rezaei example
int UTCLocalChemicalSystem_6()
{
    QString inputFile = baseDir + "in/RISA/gypsum.xml";
    QString refFile = baseDir + "ref/RISA/gypsum.xml";

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

        // Call the function
        myLocalPointer->RISA(cheprooPlusPlusPointer->mGlobalChemSysPointer->mLocChemSysVector[0]->mRelChemicalCompositions[0]);

        // Declare variable to compare: Concentrations
        vector<double> concentrations;

        concVector.resize(myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration.size());
        for(int i=0; i!=concVector.size(); i++)
        {
            concVector[i].resize(1);
            concVector[i][0] = myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration[i];
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


void UTCLocalChemicalSystem()
{
    const char* testName = "UTCLocalChemicalSystem";

    CUnitTest::RunBasicUT(testName, UTCLocalChemicalSystem_1, 1);
    CUnitTest::RunBasicUT(testName, UTCLocalChemicalSystem_2, 2);
    CUnitTest::RunBasicUT(testName, UTCLocalChemicalSystem_3, 3);
	//CUnitTest::RunBasicUT(testName, UTCLocalChemicalSystem_4, 4); 
    CUnitTest::RunBasicUT(testName, UTCLocalChemicalSystem_5, 5); 
    CUnitTest::RunBasicUT(testName, UTCLocalChemicalSystem_6, 6); 

    //CUnitTest::RunExceptionUT(testName, UTBase_1, 1, ERROR_UNINITIALIZED_XML_NODE);
}
