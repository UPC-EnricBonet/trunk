#include "cglobalchemicalsystem.h"
#include "ccheprooplusplus.h"

vector<CChemicalComposition*> generatedChemComp; // for RISA simulations - delete later
CLocalChemicalSystem* localPointer; // for RISA simulations - delete later

CGlobalChemicalSystem::CGlobalChemicalSystem(QDomElement aNode)
: CCheprooBase(aNode)
{
    this->mClassName = "CGlobalChemicalSystem";
    this->mAllPhases.clear();
    this->mAllReactions.clear();
    this->mAllSpecies.clear();
}

CGlobalChemicalSystem::CGlobalChemicalSystem()
{
    this->mClassName = "CGlobalChemicalSystem";
    this->mAllPhases.clear();
    this->mAllReactions.clear();
    this->mAllSpecies.clear();

}

CGlobalChemicalSystem::~CGlobalChemicalSystem()
{
    this->mLocalChemicalSystems.clear();
    this->mAllPhases.clear();
    this->mAllReactions.clear();
    this->mAllSpecies.clear();
}

void 
CGlobalChemicalSystem::ReadAndInitialize(const QDomElement aNode, bool isNotForUnitTest)
{
    // Get a pointer to the object manager instance
	ObjectManager* om = ObjectManager::instance();

    // Returns the document to which the node aNode belongs
    QDomDocument doc = aNode.ownerDocument();

    // Find the tag of GlobalChemicalSystem
	QDomElement nodeGlobal = XmlQt::getFirstElemByTag(doc.documentElement(),XML_TAG_GLOBALCHEMICALSYSTEM);

    // Read the option to build the GlobalChemicalSystem and the thermodynamic database name
    XmlQt::getXMLAttribute(nodeGlobal, XML_ATTR_BUILDOPTION, mBuildGCSOption); 
    XmlQt::getXMLAttribute(nodeGlobal, XML_ATTR_DATABASE, this->mDatabaseName); 

    if(this->mBuildGCSOption=="readFromMasterTemp")
    {
        // Read all the phases and related species - create phases, species and reactions objects with their attributes
        // Begin to fill the list of phases, species and reactions

		this->CreatePhases(nodeGlobal);

    }
    else cout << "error: build option for GlobalChemicalSystem not allowed";

    // Retrieve the LocalChemicalSystems element
	QDomElement kinReactions = XmlQt::getFirstElemByTag(doc.documentElement(),XML_TAG_KINETICREACTIONS);

    if(!kinReactions.isNull())
    {
        this->CreateRateLaw(kinReactions);
    }

    // Retrieve tag <proostList type = "chemicalComposition"> in the input file
    vector<QDomElement> proostListElem =  XmlQt::GetElemsByNameAndAttrib(doc.documentElement(), "proostList", 
                                                                            XML_ATTR_TYPE, CLASS_NAME_CHEMICALCOMPOSITION); 
    // Retrieve the current ChemicalComposition element in the xml 
    vector<QDomElement> currChemCompElems =  XmlQt::getElementsByTagName(proostListElem.at(0), XML_TAG_CHEMCOMP);

    vector<string> primSpecies;
    set<string> cas;  // constant activity species
    QString name;

	//*****WARNING: MODIFICATIONS FOR RISA SIMULATIONS
	// Read file with perturbed values

    vector<double> perturbed_values1;
	vector<double> perturbed_values2;
	vector<double> perturbed_values3;
	vector<double> perturbed_values4;


    for(int i=0; i!=currChemCompElems.size(); i++) // Loop over initial chemical compositions
    {
        // Retrieve chemical composition name
        XmlQt::getXMLAttribute(currChemCompElems.at(i), XML_ATTR_NAME, name);

        // Create chemical composition
		CChemicalComposition* chemComp = dynamic_cast<CChemicalComposition*> (om->newInstanceOfClassName(CLASS_NAME_CHEMICALCOMPOSITION, name));

        // Read chemical composition
        chemComp->Read(currChemCompElems.at(i), primSpecies, cas);

        QDomElement locChemSysElem;
        CLocalChemicalSystem* currLocal = 0;

        if(i==0)
        {
            // Retrieve the LocalChemicalSystems element
		    locChemSysElem = XmlQt::getFirstElemByTag(doc.documentElement(),XML_TAG_LOCALCHEMICALSYSTEMS);

            // Create local
            this->CreateLocal(locChemSysElem, primSpecies, cas, chemComp, i, true);
        }
        else
        {
            // First check if this water belongs to some locals previously created
            bool itBelongs = false;
            itBelongs = this->CheckChemComp(cas, currLocal);

            // If not, create new local and set its attributes
            if(!itBelongs) 
            {
                this->CreateLocal(locChemSysElem, primSpecies, cas, chemComp, i);
            }
            else // If yes, build the association
            {
                currLocal->mRelChemicalCompositions.push_back(chemComp);
                chemComp->mLocChemSysPointer = currLocal;
                this->mChemCompToLocal.insert(pair<CChemicalComposition*,CLocalChemicalSystem*>(chemComp,currLocal));
            }
        }

        // Set the chemical composition attributes
        chemComp->SetChemComp();

        //*****WARNING: MODIFICATIONS FOR RISA SIMULATIONS
        
        //for(double p=0.05; p!=2.0500000000000007; p+=0.05)
        //{
        //    CChemicalComposition* newChemComp = new CChemicalComposition(chemComp, p);
        //    generatedChemComp.push_back(newChemComp);
        //    //this->mChemCompToLocal.insert(pair<CChemicalComposition*,CLocalChemicalSystem*>(newChemComp,newChemComp->mLocChemSysPointer));
        //}
     //       QString fileName = "perturbed_values_"  + chemComp->name() + ".out";
	    //    //QString fileName = "perturbed_values.out";
     //       QFile inputFile(fileName);

     //       // Open the file, if is not open
     //       if(!inputFile.isOpen())
     //       {
     //           if ( !inputFile.open(QFile::ReadOnly | QIODevice::Text))
     //           {
     //               this->PrThrow(ERROR_FILEOPEN);
     //           }
     //       }

     //       QTextStream input_stream(&inputFile);

     //       double value;

     //       while(!input_stream.atEnd())
     //       {
     //           QString line = input_stream.readLine();
		   //     QStringList list = line.split(QRegExp("\\s+"));

		   //     perturbed_values1.push_back(list.at(0).toDouble());
     //           //perturbed_values2.push_back(list.at(1).toDouble());

     //           if(chemComp->name() == "water2")
     //           {
     //               perturbed_values2.push_back(list.at(1).toDouble());

     //           }
     //           else if(chemComp->name() == "water3")
     //           {
		   //         perturbed_values2.push_back(list.at(1).toDouble());
     //               perturbed_values3.push_back(list.at(2).toDouble());
     //           }
		   //     else if(chemComp->name() == "water4")
     //           {
     //               perturbed_values2.push_back(list.at(1).toDouble());
     //               perturbed_values3.push_back(list.at(2).toDouble());
     //               perturbed_values4.push_back(list.at(3).toDouble());
     //           }

     //       }
	    //
     //   // Create chemical compositions from the perturbed values that have just been read
     //   for(int k=0; k!= perturbed_values1.size();k++)
     //   {
     //       double p = 0.0;
     //       double q = 0.0;
     //       double r = 0.0;
     //       double s = 0.0;
     //       if(chemComp->name() == "water1")
     //       {
     //           CChemicalComposition* newChemComp = new CChemicalComposition(chemComp, p, perturbed_values1[k], q, r, s, true);
     //           //CChemicalComposition* newChemComp = new CChemicalComposition(chemComp, perturbed_values1[k], perturbed_values2[k], q, r, s, true);
     //           generatedChemComp.push_back(newChemComp);
     //       }
     //       else if(chemComp->name() == "water2")
     //       {
     //           CChemicalComposition* newChemComp = new CChemicalComposition(chemComp, p, perturbed_values1[k], perturbed_values2[k], q, r, true);
     //           //CChemicalComposition* newChemComp = new CChemicalComposition(chemComp, perturbed_values1[k], perturbed_values2[k], perturbed_values3[k], q, r, true);
     //           generatedChemComp.push_back(newChemComp);
     //       }
     //       else if(chemComp->name() == "water3")
     //       {
     //           CChemicalComposition* newChemComp = new CChemicalComposition(chemComp, p, perturbed_values1[k], perturbed_values2[k], perturbed_values3[k], r, true);
     //           //CChemicalComposition* newChemComp = new CChemicalComposition(chemComp, perturbed_values1[k], perturbed_values2[k], perturbed_values3[k], perturbed_values4[k], r, true);
     //           generatedChemComp.push_back(newChemComp);
     //       }
     //       else if(chemComp->name() == "water4")
     //       {
     //           CChemicalComposition* newChemComp = new CChemicalComposition(chemComp, s, perturbed_values1[k], perturbed_values2[k], perturbed_values3[k], perturbed_values4[k], true);
     //           //CChemicalComposition* newChemComp = new CChemicalComposition(chemComp, s, perturbed_values1[k], perturbed_values2[k], perturbed_values3[k], perturbed_values4[k], true);
     //           generatedChemComp.push_back(newChemComp);
     //       }
     //   }

     //   cout << "finished creating waters from:" << chemComp->name().toStdString() << endl;
     //       
	    //perturbed_values1.resize(0);
	    //perturbed_values2.resize(0);
	    //perturbed_values3.resize(0);
	    //perturbed_values4.resize(0);
	    //*****WARNING: END OF MODIFICATIONS FOR RISA SIMULATIONS

    } 



    // Speciate initial waters
    if(isNotForUnitTest)
    {
        // Delete output files that eventually have been generated before
        this->DeleteOutputFiles();

	    map<CChemicalComposition*, CLocalChemicalSystem*>::iterator it;
	    for(it=this->mChemCompToLocal.begin(); it!=this->mChemCompToLocal.end(); it++)
	    {
           it->second->SpeciateInitialWater(it->first);
           //it->second->RISA(it->first); // To check RISA algorithm
           //localPointer = it->second; // For RISA simulations only - delete later
        }

        //*****WARNING: MODIFICATIONS FOR RISA SIMULATIONS
        //for(int k=0; k!=generatedChemComp.size(); k++)
        //{
        //    localPointer->RISA(generatedChemComp.at(k), false);

        //    cout<< k << "\n";
        //}
        //*****WARNING: END OF MODIFICATIONS FOR RISA SIMULATIONS

		
	}
    else
    {
        // Don't do anything, Unit Test is calling the specific function
    }
}

void 
CGlobalChemicalSystem::CreateLocal(const QDomElement aNode, vector<string>& primSpecies, set<string>& cas, CChemicalComposition* chemComp, int i, bool isFirstLocal)
{
    // Get a pointer to the object manager instance
	ObjectManager* om = ObjectManager::instance();

    if(isFirstLocal) // Then I read its specialization
    {
        // Retrieve the LocalChemicalSystems defined in the xml (only one)
        vector<QDomElement> localChemSysElems =  XmlQt::getElementsByTagName(aNode,XML_TAG_LOCALCHEMICALSYSTEM);
 
        // Create the LocalChemicalSystem and put it in the vector
        this->mLocChemSysVector.push_back(dynamic_cast<CLocalChemicalSystem*>(om->newInstanceFromXML(localChemSysElems.at(0))));

		// Assign pointer to current GlobalChemicalSystem
		this->mLocChemSysVector.back()->mGlobChemSysPointer = this;

        // Call the Read function of CLocalChemicalSystem
        this->mLocChemSysVector.back()->Read(localChemSysElems.at(0), primSpecies, cas);
    }
    else
    {
        //Generate the name for the local
        stringstream sstm;
        string className = "LocalChemicalSystem_";
        sstm << className << i;
        QString name = QString::fromStdString(sstm.str());
        
        // Create the LocalChemicalSystem NOT from xml (by default, it creates a Local Chemical System Saaltink)
        this->mLocChemSysVector.push_back(dynamic_cast<CLocalChemicalSystem*>(om->newInstanceOfClassName(CLASS_NAME_LOCALCHEMICALSYSTEMSAALTINK, name)));

		// Assign pointer to current GlobalChemicalSystem
		this->mLocChemSysVector.back()->mGlobChemSysPointer = this;

        // Copy convergence parameters from a previously defined local (by default, the first that has been defined and read from xml)
        this->mLocChemSysVector.back()->SetConvParams(this->mLocChemSysVector[0]);

        // Call the Set function of CLocalChemicalSystem
        this->mLocChemSysVector.back()->SetLocal(primSpecies, cas);
    }

    // Assign correspondence between this local and the chemical composition
    this->mLocChemSysVector.back()->mRelChemicalCompositions.push_back(chemComp);
    chemComp->mLocChemSysPointer = this->mLocChemSysVector.back();
    this->mChemCompToLocal.insert(pair<CChemicalComposition*,CLocalChemicalSystem*>(chemComp,this->mLocChemSysVector[0]));



}

void 
CGlobalChemicalSystem::CreatePhases(const QDomElement aNode)
{
	// Get a pointer to the object manager instance
	ObjectManager* om = ObjectManager::instance();

	// Retrieve all the nodes with XML_TAG_PHASE
    QDomNodeList xmlNodeList = aNode.elementsByTagName(XML_TAG_PHASE);
    if(!xmlNodeList.isEmpty())
    {
        const int nodeListSize = xmlNodeList.size();

        for(unsigned i=0;i!=nodeListSize;i++)
        {
            QDomElement aPhaseNode = xmlNodeList.at(i).toElement();

            // Create the correct specialization of the phase and associate the GlobalChemicalSystem
            CPhase* aPhase = dynamic_cast<CPhase*>(om->newInstanceFromXML(aPhaseNode));
            aPhase->mGlobChemSysPointer = this;

            // Retrieve the phase model
            QString aPhaseModel;
            XmlQt::getXMLAttribute(aPhaseNode, XML_ATTR_MODEL, aPhaseModel); 
            aPhase->mPhaseKind = CPhase::str2PhasesKind(aPhaseModel);


            // Read the phase (note that this will trigger the reading of the species and the insertion of the species in the mAllSpecies list)
            aPhase->Read(aPhaseNode);

			this->mAllPhases.insert(pair<string, CPhase*>(aPhase->name().toStdString(),aPhase));

        }

    }
}

void 
CGlobalChemicalSystem::CreateReaction(QString aRelSpeciesName)
{
    ObjectManager* om = ObjectManager::instance();

    // It creates the reaction, whose name is the name of the secondary species associated
    CReaction* aReaction = dynamic_cast<CReaction*>(om->newInstanceOfClassName(CLASS_NAME_REACTION, aRelSpeciesName));

	// Convert QString to std string
    string aReactionName = aRelSpeciesName.toStdString();

    // Initialize the pointer of the species to the GlobalChemicalSystem (CHECK THIS!!!)
    aReaction->mGlobChemSysPointer = this;

    // Add the species to the mAllSpecies list
    this->mAllReactions.insert(pair<string, CReaction*>(aReactionName,aReaction));

    // Check if the xml database exists, if not show an error and return 1

    QString inputFile = this->mDatabaseName;

	QFile tstFile;

	tstFile.setFileName(mDatabaseName);

	if( ! tstFile.exists() )
	{
		cout << PrTr("CGlobalChemicalSystem", "Thermodynamic database file does not exist: ") <<  qPrintable(inputFile) << endl;
	}

	// Open inputFile and parse the xml. The result is stored in doc
	QDomDocument doc("mydocument");
	XmlQt::xmlParseFile(inputFile, &doc);

    // Now read the species attributes from the thermodynamic database
    aReaction->Read(doc.documentElement());    

}

void
CGlobalChemicalSystem::CreateSpecies(QString aSpeciesName, QString aSpeciesModel)
{
	CSpecies* spc = new CSpecies(aSpeciesModel);

    // Assign the correct name to the species
    spc->SetName(aSpeciesName);

	// Convert QString to std string
    string aSpcName = aSpeciesName.toStdString();

    // Initialize the pointer of the species to the GlobalChemicalSystem (CHECK THIS!!!)
    spc->mGlobChemSysPointer = this;

    // Add the species to the mAllSpecies list
    this->mAllSpecies.insert(pair<string, CSpecies*>(aSpcName,spc));

    // Check if the xml database exists, if not show an error and return 1

    QString inputFile = this->mDatabaseName;

	QFile tstFile;

	tstFile.setFileName(mDatabaseName);

	if( ! tstFile.exists() )
	{
		cout << PrTr("CGlobalChemicalSystem", "Thermodynamic database file does not exist: ") <<  qPrintable(inputFile) << endl;
	}

	// Open inputFile and parse the xml. The result is stored in doc
	QDomDocument doc("mydocument");
	XmlQt::xmlParseFile(inputFile, &doc);

    // Now read the species attributes from the thermodynamic database
    spc->Read(doc.documentElement());    

}

void
CGlobalChemicalSystem::Speciate()
{
    for(int i=0; i!=mLocChemSysVector.size(); i++)
    {
        //Speciate LocalChemicalSystems 
        this->mLocChemSysVector.at(i)->Speciate();

    }

}

void
CGlobalChemicalSystem::CopyGlobalChemSysIntoLocal(CLocalChemicalSystem* aLocalChemSystem)
{

}

void 
CGlobalChemicalSystem::ComputeSI(CChemicalComposition* aCurrChemComp, map<CSpecies*, int> speciesIndices, set<string> &cas)
{
    // Check over CAS
    set<string>::iterator it;
    for(it=cas.begin(); it!=cas.end(); it++)
    {
        // Skip water
        if(*it == "h2o") continue;

        // Retrieve corresponding mineral reaction
		CReaction* reaction = this->mAllReactions[*it];


    }

	// Declare and initialize variables
	double value = 0;
	double min_conc = 0;
	int key = 0;

	// Loop over all possible minerals
	for(int i=0; i!=this->mMinNames.size(); i++)
	{
		// Retrieve corresponding mineral reaction
		CReaction* reaction = this->mAllReactions[this->mMinNames[i]];

		// Evaluate the saturation index of the current mineral
		// value = ...

		// Get the mineral concentration from aCurrChemComp
		// min_conc = ...

		// Update key of the LocalChemicalSystem this ChemicalComposition should belong
		if(value > 0 || min_conc > 0) key += pow(2.0,i);
	}

	// Check if aCurrChemComp needs to change LocalChemicalSystem
	if (this->mLocalChemicalSystems[aCurrChemComp->mLocChemSysPointer] != key)
	{
		// Change LocalChemicalSystem
		this->ChangeLocal(aCurrChemComp, this->mLocalChemicalSystems[aCurrChemComp->mLocChemSysPointer], key);
	}

}

void 
CGlobalChemicalSystem::ChangeLocal(CChemicalComposition* aChemComp, int key_old, int key_new)
{
	// Delete aChemComp from the previous LocalChemicalSystem and add it to the new LocalChemicalSystem
	map<CLocalChemicalSystem*, int>::iterator it;
	for(it=this->mLocalChemicalSystems.begin(); it!=this->mLocalChemicalSystems.end(); it++)
	{
		if(it->second == key_old) // Delete it from the old 
		{
			it->first->mChemicalCompositions.erase(aChemComp->name().toStdString());
		}

		if(it->second == key_new) // Add it to the new
		{
			it->first->mChemicalCompositions.insert(pair<string,CChemicalComposition*>(aChemComp->name().toStdString(), aChemComp));
			aChemComp->mLocChemSysPointer = it->first;
		}
	}
}

void
CGlobalChemicalSystem::MixWaters(map<string,double> &aWatersToMixRatios, CChemicalComposition* aMixedWater)
{
    map<CChemicalComposition*, double> watersMap;

    // Look for waters in local chemical systems
    map<string,double>::iterator it;
    for(it=aWatersToMixRatios.begin(); it!=aWatersToMixRatios.end(); it++)
    {
        for(int i=0; i!=this->mLocChemSysVector.size(); i++)
        {
            for(int j=0; j!= this->mLocChemSysVector[i]->mRelChemicalCompositions.size(); j++)
            {
                if(it->first == this->mLocChemSysVector[i]->mRelChemicalCompositions[j]->name().toStdString())
                {
                    watersMap.insert(pair<CChemicalComposition*, double>(this->mLocChemSysVector[i]->mRelChemicalCompositions[j], it->second));
                }
            }
        }
    }

    // Add this water in the mChemCompToLocal
    this->mChemCompToLocal.insert(pair<CChemicalComposition*,CLocalChemicalSystem*>(aMixedWater,this->mLocChemSysVector[0]));

    this->DeleteOutputFiles();

    this->mLocChemSysVector[0]->MixWaters(watersMap, aMixedWater, false);
}

void
CGlobalChemicalSystem::MixChemicalCompositions(map<string,double> &aWatersToMixRatios, CChemicalComposition* aMixedWater)
{
    map<CChemicalComposition*, double> watersMap;

    // Look for waters in local chemical systems
    map<string,double>::iterator it;
    for(it=aWatersToMixRatios.begin(); it!=aWatersToMixRatios.end(); it++)
    {
        for(int i=0; i!=this->mLocChemSysVector.size(); i++)
        {
            for(int j=0; j!= this->mLocChemSysVector[i]->mRelChemicalCompositions.size(); j++)
            {
                if(it->first == this->mLocChemSysVector[i]->mRelChemicalCompositions[j]->name().toStdString())
                {
                    watersMap.insert(pair<CChemicalComposition*, double>(this->mLocChemSysVector[i]->mRelChemicalCompositions[j], it->second));
                }
            }
        }
    }

    this->mLocChemSysVector[0]->MixChemicalCompositions(watersMap, aMixedWater, false);
}

bool
CGlobalChemicalSystem::CheckChemComp(set<string>& cas, CLocalChemicalSystem* &local)
{
    bool returnValue = false;

    for(int i=0; i!=this->mLocChemSysVector.size(); i++)
    {
        if(cas == this->mLocChemSysVector.at(i)->mCAS) 
        {
            local = this->mLocChemSysVector.at(i);
            returnValue = true;
        }
    }

    return returnValue;
}



void
CGlobalChemicalSystem::DeleteOutputFiles()
{
    QString fileName = "output.out";
    fileName = QDir(fileName).absolutePath();

    QFile* outputFile = new QFile(fileName);

    outputFile->remove();

	map<CChemicalComposition*, CLocalChemicalSystem*>::iterator it;
    QString fileName1;
    QFile* outputFile1;
    QString fileName2;
    QFile* outputFile2;
    QString fileName3;
    QFile* outputFile3;
    QString fileName4;
    QFile* outputFile4;
	for(it=this->mChemCompToLocal.begin(); it!=this->mChemCompToLocal.end(); it++)
	{
        fileName1 = "output_RISA_" + it->first->name() + ".out";
        fileName2 = "output_" + it->first->name() + ".out";
        fileName3 = "output_" + it->first->name() + "_conc.out";
        fileName4 = "output_" + it->first->name() + "_act_coeffs.out";

        fileName1 = QDir(fileName1).absolutePath();
        fileName2 = QDir(fileName2).absolutePath();
        fileName3 = QDir(fileName3).absolutePath();
        fileName4 = QDir(fileName4).absolutePath();

        outputFile1 = new QFile(fileName1);
        outputFile2 = new QFile(fileName2);
        outputFile3 = new QFile(fileName3);
        outputFile4 = new QFile(fileName4);

        outputFile1->remove();
        outputFile2->remove();
        outputFile3->remove();
        outputFile4->remove();


        delete outputFile1;
        delete outputFile2;
        delete outputFile3;
        delete outputFile4;
    }

    delete outputFile;
}

void
CGlobalChemicalSystem::CreateRateLaw(const QDomElement aNode)
{
 //   // Get a pointer to the object manager instance
	//ObjectManager* om = ObjectManager::instance();

 //   // Retrieve the kinetic reactions defined in the xml
 //   vector<QDomElement> kinReacElems =  XmlQt::getElementsByTagName(aNode,XML_TAG_KINETICREACTION);

 //   for(int i=0; i!=kinReacElems.size(); i++)
 //   {
 //       // Create the reaction
 //       CReaction* aReaction = dynamic_cast<CReaction*>(om->newInstanceFromXML(kinReacElems[i]));

 //       // Read the reaction
 //       aReaction->Read(kinReacElems[i], true);

 //       // Put it in the map of reactions
 //       this->mAllReactions.insert(pair<string,CReaction*>(aReaction->name().toStdString(),aReaction));
 //   }
      

}

CChemicalComposition* 
CGlobalChemicalSystem::GetChemicalComposition(string aName)
{
    ObjectManager* om = ObjectManager::instance();

    CChemicalComposition* aChemComp = 0;

    // Look up the chemicalcomposition
    map<CChemicalComposition*, CLocalChemicalSystem*>::iterator it;
    for(it=this->mChemCompToLocal.begin(); it!=this->mChemCompToLocal.end(); it++)
    {
        if(it->first->name().toStdString() == aName)
        {
            // Create new chemical composition
            aChemComp = dynamic_cast<CChemicalComposition*> (om->newInstanceOfClassName(CLASS_NAME_CHEMICALCOMPOSITION, it->first->name()));

            // Copy it->first in the new chemical composition
            aChemComp->mLocChemSysPointer = it->first->mLocChemSysPointer;
            aChemComp->CopyForRT(it->first);
        }
    }

    return aChemComp;
}
