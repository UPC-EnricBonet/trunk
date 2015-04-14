// Unit Testing CAqueousIdeal

#include "cunittest.h"
#include "xmlconstants.h"
#include "xmlqt.h"
#include "objectlist.h"
#include "commonglobals.h"
#include "objectmanager.h"
#include "valarraytools.h"
#include "classnameconstants.h"
#include "ccheprooplusplus.h"
#include "cglobalchemicalsystem.h"

#include "caqueous.h"
#include "cphase.h"
#include "caqueousideal.h"

#include <vector>
#include <QFile>

static QString baseDir = "UnitTest/";

// Method to evaluate the activity coefficients
int UTCAqueousIdeal_1()
{
    QString inputFile = baseDir + "in/UTCheprooInput2.xml";
    QString refFile = baseDir + "ref/UTCAqueousIdeal1.xml";

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
    valarray<valarray<double> > actCoeffs;

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
        vector<double> activityCoeffs;
        activityCoeffs.resize(myLocalPointer->mPrimSpeciesNum + myLocalPointer->mEqReactions.size());


        // Check the method
        CPhase* myAqPhase = cheprooPlusPlusPointer->mGlobalChemSysPointer->mAllPhases["aqueous1"];
        myAqPhase->ComputeActivityCoeff(myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration, activityCoeffs, 
                                        myLocalPointer->mSpecies[eAqueous1nc], myLocalPointer->mSpecies[eAqueous2],
                                        myLocalPointer->mSpeciesIndices);

        actCoeffs.resize(activityCoeffs.size());
        for(int i=0; i!=activityCoeffs.size(); i++)
        {
            actCoeffs[i].resize(1);
            actCoeffs[i][0] = activityCoeffs[i];
        }

    }

    //Obtains from sRefBal the content of reference file to compare for Balance
	QFile RefFileBal(refFile);
	RefFileBal.open(QIODevice::ReadOnly | QIODevice::Text);
	QString sRefBal = RefFileBal.readAll();
    RefFileBal.close();
    valarray< valarray< double > > ref;
    ValarrayTools::fromString(sRefBal, ref);

    
	if ( !ValarrayTools::IsEqual(actCoeffs,ref, 1e-12)) return 1;
    
	return 0;
}

void UTCAqueousIdeal()
{
    const char* testName = "UTCAqueousIdeal";

    CUnitTest::RunBasicUT(testName, UTCAqueousIdeal_1, 1);

    //CUnitTest::RunExceptionUT(testName, UTBase_1, 1, ERROR_UNINITIALIZED_XML_NODE);
}
