#include "ccheprooplusplus.h"


CCheprooPlusPlus::CCheprooPlusPlus()
{
    this->mClassName="CCheprooPlusPlus";

}

CCheprooPlusPlus::CCheprooPlusPlus(QDomElement aNode)
: CCheprooBase(aNode)
{
    this->mClassName="CCheprooPlusPlus";

}

CCheprooPlusPlus::~CCheprooPlusPlus(void)
{
}

void
CCheprooPlusPlus::ReadAndInitialize(const QDomElement aNode, bool isNotForUnitTest)
{
    // Find the Function Tag and read the name of the method is going to be used
    const QDomElement nodeFunction =XmlQt::getFirstElemByTag(aNode, XML_TAG_FUNCTION); 

    XmlQt::getXMLAttribute(nodeFunction, XML_ATTR_NAME, mFunctionName); 

    // Read list of referenced waters
    vector<QString> listOfIDs;
    vector<QDomElement> waterTags;
    this->ReadReferences(nodeFunction, 
                         XML_TAG_CHEMCOMP, 
                         waterTags, 
                         &listOfIDs);

	//Check that water refs are not void
	if(waterTags.size() == 0)
	{
		this->PrThrow(ERROR_WATERS_MISSING);
	}

    // Put the names of the referenced waters in a vector
    vector<QString> aWaterNames;
    double attributeVal;

    // Add the waters referenced to a list that will be passed to LocalChemicalSystem
    for (unsigned i = 0 ; i < waterTags.size(); i++ )
	{
		// Retrieve the water and add it to a vector that will be passed to the LocalChemicalSystem
		const QString curNodeText = waterTags[i].text();
 
        // Add the name of the water to the vector of all the names of the waters
        aWaterNames.push_back(curNodeText);

        if(this->mFunctionName == "mixWaters")
        {
            // Get the mixing ratio value
            XmlQt::getXMLAttribute(waterTags[i], XML_ATTR_MIXINGRATIO, attributeVal);

            // Save them in the map
            this->mWatersToMixRatios.insert(pair<string,double>(curNodeText.toStdString(), attributeVal));
        }

 
	}
    waterTags.clear();

    // Create a GlobalChemicalSystem
    this->CreateGlobal();

    // Read and initialize the problem
    this->mGlobalChemSysPointer->ReadAndInitialize(aNode, isNotForUnitTest);


}

void 
CCheprooPlusPlus::CreateGlobal()
{
    ObjectManager* om = ObjectManager::instance();

    this->mGlobalChemSysPointer = dynamic_cast<CGlobalChemicalSystem*>(om->newInstanceOfClassName(CLASS_NAME_GLOBALCHEMICALSYSTEM, "GlobalChemicalSystem"));
}

bool
CCheprooPlusPlus::SpeciateWaters()
{
    // Call the speciation function
    this->mGlobalChemSysPointer->Speciate();

    return true;
}

bool
CCheprooPlusPlus::MixWaters()
{
    ObjectManager* om = ObjectManager::instance();

    // Create a chemical composition for the mixed water
    mAuxChemComp = dynamic_cast<CChemicalComposition*> (om->newInstanceOfClassName(CLASS_NAME_CHEMICALCOMPOSITION, "mixedWater"));

    this->mGlobalChemSysPointer->MixWaters(mWatersToMixRatios, mAuxChemComp);
 
    return true;

}

bool
CCheprooPlusPlus::MixChemicalCompositions()
{
    ObjectManager* om = ObjectManager::instance();

    // Create a chemical composition for the mixed water
    mAuxChemComp = dynamic_cast<CChemicalComposition*> (om->newInstanceOfClassName(CLASS_NAME_CHEMICALCOMPOSITION, "mixedWater"));

    this->mGlobalChemSysPointer->MixChemicalCompositions(mWatersToMixRatios, mAuxChemComp);
 
    return true;

}

bool
CCheprooPlusPlus::EvaporateWaters()
{
    return true;
}

void 
CCheprooPlusPlus::ReadAndInitialize(QString aFileName)
{
    // Check if the xml file exists, if not show an error and return 1
	QFile tstFile;

	tstFile.setFileName(aFileName);

	if( ! tstFile.exists() )
	{
		cout << PrTr("main", ERROR_INPUT_FILE_DOES_NOT_EXIST) <<  qPrintable(aFileName) << endl;
        return;
	}

	// Open inputFile and parse the xml. The result is stored in doc
	QDomDocument doc("mydocument");
	XmlQt::xmlParseFile(aFileName, &doc);

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
	QDomElement nodeGlobal = XmlQt::getFirstElemByTag(doc.documentElement(),XML_TAG_GLOBALCHEMICALSYSTEM);

	if( !nodeGlobal.isNull() ) 
	{
        // Create a GlobalChemicalSystem
        this->CreateGlobal();

        // Read and initialize the problem
        this->mGlobalChemSysPointer->ReadAndInitialize(nodeGlobal);
    }
}

CChemicalComposition* 
CCheprooPlusPlus::GetChemicalComposition(string aName)
{
    CChemicalComposition* aChemComp = 0;

    aChemComp = mGlobalChemSysPointer->GetChemicalComposition(aName);

    return aChemComp;
}
