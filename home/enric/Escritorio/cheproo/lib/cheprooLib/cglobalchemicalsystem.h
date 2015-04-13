#ifndef CGLOBALCHEMICALSYSTEM_H
#define CGLOBALCHEMICALSYSTEM_H

#include "ccheproobase.h"

#include "clocalchemicalsystem.h"
#include "cphase.h"
#include "cspecies.h"
#include "creaction.h"

#include <QString>
#include <QFile>
#include <string>
#include <vector>
#include <map>
using namespace std;

class CCheprooPlusPlus;

class CGlobalChemicalSystem : public CCheprooBase
{
public:

    map<string,CReaction*> mAllReactions;  // Map of all reaction in the system
                                            // The key is the name of the secondary species corresponding to the reaction

    map<string,CSpecies*> mAllSpecies;  // Map of all species in the system
                                         // The key is the name of the species

    map<string,CPhase*> mAllPhases;  // Map of all phases in the system   
                                      // The key is the name of the phase (aqueous, mineral, gas or surface)

    vector<CLocalChemicalSystem*> mLocChemSysVector;  // Vector of LocalChemicalSystems

	vector<string> mMinNames; // Vector with names of possible minerals ***to be filled, just an idea***

	map<CLocalChemicalSystem*, int> mLocalChemicalSystems; // Map with LocalChemicalSystems. The Key is a number corresponding to the minerals that are present in the system

    map<CChemicalComposition*, CLocalChemicalSystem*> mChemCompToLocal; // Map with the correspondence between chemical compositions and locals

    /// \brief Constructor from a XML file.
	CGlobalChemicalSystem(
		   QDomElement aNode ///< XML node with attributes
		   );


	CGlobalChemicalSystem();    
	~CGlobalChemicalSystem();

private:

    bool mLocChemSysHasBeenCreated;  // bool that indicates if the LocalChemicalSystem has been created or not

    QString mDatabaseName;

    QString mBuildGCSOption;


public:

    /// brief Initialize the problem for Cheproo++ standalone version
    void ReadAndInitialize(const QDomElement aNode, bool isNotForUnitTest = true);

    /// brief Creates and assign the local chemical systems for the cheproo standalone version
    void CreateLocal(const QDomElement aNode,
                     vector<string>& primSpecies, 
                     set<string>& cas,
                     CChemicalComposition* chemComp,
                     int i,
                     bool isFirstLocal = false);

	/// brief Creates phases
	void CreatePhases(const QDomElement aNode);

    /// brief Creates reactions
    void CreateReaction(QString aRelSpeciesName);

	/// brief Creates phases
	void CreateRateLaw(const QDomElement aNode);

    /// brief Creates reactions
    void CreateSpecies(QString aSpeciesName,
                       QString aSpeciesModel);


    /// brief Speciates local chemical systems
    void Speciate();

    /// brief Mix waters given mixing ratios values
    void MixWaters(map<string,double> &aWatersToMixRatios,
                   CChemicalComposition* aMixedWater);

    /// brief Mix waters given mixing ratios values
    void MixChemicalCompositions(map<string,double> &aWatersToMixRatios,
                                 CChemicalComposition* aMixedWater);

	/// brief Evaluate saturation index for a given ChemicalComposition
	void ComputeSI(CChemicalComposition* aCurrChemComp,   ///< ChemicalComposition pointer
				   map<CSpecies*, int> speciesIndices,    ///< Species indices;
                   set<string> &cas);          

    /// brief Updates local chemical system for a given ChemicalComposition
    void ChangeLocal(CChemicalComposition* aChemComp,     ///< ChemicalComp
		             int key_old,                         ///< Key of the old LocalChemicalSystem
					 int key_new);                        ///< Key of the new LocalChemicalSystem

    /// \brief Returns a pointer to a chemical composition with a specific name
    CChemicalComposition* GetChemicalComposition(string aName);


private:

    /// brief Function that checks if a chemical composition belongs to a local chemical system already present
    bool CheckChemComp(set<string>& cas, 
                       CLocalChemicalSystem* &local); // pointer to chemical composition necessary for checking after calculations
    
    /// brief Function that copies all the reactions, species and phases in one CLocalChemicalSystem
    void CopyGlobalChemSysIntoLocal(CLocalChemicalSystem* aLocalChemSystem);

    /// brief Delete output files  
    void DeleteOutputFiles();

};

#endif
