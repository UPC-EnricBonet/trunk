// Unit Testing CLocalChemicalSystemSaaltink

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
#include "clocalchemicalsystemsaaltink.h"

#include <vector>
#include <QFile>
#include <math.h>

static QString baseDir = "UnitTest/";

// We check the method to build the component matrix (EU) and Uc
int UTCLocalChemicalSystemSaaltink_1()
{
    QString inputFile = baseDir + "in/UTCLocalChemicalSystem/UTCLocalChemicalSystemSaaltink1.xml";
    QString refFile = baseDir + "ref/UTCLocalChemicalSystem/UTCLocalChemicalSystemSaaltink1.xml";

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
    valarray<valarray<double> > comp_matrix;

    if( CheprooPlusPlusExists ) 
    {
	    // Create a Problem object and read its properties from the xml
        CCheprooPlusPlus* cheprooPlusPlusPointer = dynamic_cast<CCheprooPlusPlus*>(om->newInstanceOfClassName(CLASS_NAME_CHEPROOPLUSPLUS, aProblemName));

        // Read and initialize the problem
        cheprooPlusPlusPointer->ReadAndInitialize(nodeCheprooPlusPlus, false);

        // Pointer to LocalChemicalSystem
        CLocalChemicalSystem* myLocalPointer = cheprooPlusPlusPointer->mGlobalChemSysPointer->mLocChemSysVector[0];

        MatrixXd ComponentMatrix = myLocalPointer->GetComponentMatrix();

        // Declare variable to compare: Activity coefficients
        vector<double> concentrations;

        comp_matrix.resize(myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration.size());
        for(int i=0; i!=comp_matrix.size(); i++)
        {
            comp_matrix[i].resize(1);
            comp_matrix[i][0] = ComponentMatrix(0,i);
        }
        
    }

    //Obtains from sRefBal the content of reference file to compare for Balance
	QFile RefFileBal(refFile);
	RefFileBal.open(QIODevice::ReadOnly | QIODevice::Text);
	QString sRefBal = RefFileBal.readAll();
    RefFileBal.close();
    valarray< valarray< double > > ref;
    ValarrayTools::fromString(sRefBal, ref);

    
	if ( !ValarrayTools::IsEqual(comp_matrix,ref, 1e-9)) return 1;
    
	return 0;
}


void UTCLocalChemicalSystemSaaltink()
{
    const char* testName = "UTCLocalChemicalSystemSaaltink";

    CUnitTest::RunBasicUT(testName, UTCLocalChemicalSystemSaaltink_1, 1);

    //CUnitTest::RunExceptionUT(testName, UTBase_1, 1, ERROR_UNINITIALIZED_XML_NODE);
}
