#ifndef CCHEPROOPLUSPLUS_H
#define CCHEPROOPLUSPLUS_H

#include "ccheproobase.h"

#include "cglobalchemicalsystem.h"

#include <string>
#include <vector>
using namespace std;

extern void gInit_cheproo(void);

static const char* ERROR_FUNCTIONATTRIBUTE_MISSING = QT_TRANSLATE_NOOP("CCheprooPlusPlus", "Error: Type of function not allowed");
static const char* ERROR_WATERS_MISSING = QT_TRANSLATE_NOOP("CCheprooPlusPlus", "Error: Waters not found");
static const char* ERROR_NUMBER_OF_MIXING_STEPS_MISSING = QT_TRANSLATE_NOOP("CCheprooPlusPlus", "Error: Number of mixing steps not found");

static const char* ERROR_INPUT_FILE_DOES_NOT_EXIST = QT_TRANSLATE_NOOP("main", "Input file does not exist: ");
static const char* INPUT_FILENAME_NAME = QT_TRANSLATE_NOOP("main", "Input filename: ");


class CCheprooPlusPlus : public CCheprooBase
{
public:

    CGlobalChemicalSystem* mGlobalChemSysPointer;  // Pointer to the CGlobalChemicalSystem

    string mFunctionName;

    map<string, double> mWatersToMixRatios;


private:

    int stepsNumber; 
    CChemicalComposition* mAuxChemComp;

public:

    /// \brief Constructor from a XML file.
	CCheprooPlusPlus(
		   QDomElement aNode ///< XML node with attributes
           );

     /// brief Speciate the waters defined in the input file
    bool SpeciateWaters();

	/// brief Speciate the waters defined in the input file from c1
    bool Speciate_from_c1(vector<double> &c1);

    /// brief Mixes more waters and puts them in a water with name "aMixedWater"
    bool MixWaters();

    /// brief Mixes more waters and puts them in a water with name "aMixedWater"
    bool MixChemicalCompositions();

    /// brief Evaporate waters 
    bool EvaporateWaters();

    /// brief Calls the Read function of GlobalChemicalSystem that reads everything and creates LocalChemicalSystems
    void ReadAndInitialize(QString aFileName);

    /// brief Read and initialize the problem for Cheproo++ standalone version
    void ReadAndInitialize(const QDomElement aNode, 
                           bool isNotForUnitTest = true);

    /// \brief Returns a pointer to a chemical composition with a specific name
    CChemicalComposition* GetChemicalComposition(string aName);

    CCheprooPlusPlus();
	~CCheprooPlusPlus();


private:

    void CreateGlobal();

};

#endif
