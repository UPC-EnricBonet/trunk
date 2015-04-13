#ifndef CSPECIES_H
#define CSPECIES_H

#include "ccheproobase.h"
#include "classnameconstants.h"
#include "objectmanager.h"
#include "xmlconstants.h"
#include "xmlqt.h"
//#include "proostglobals.h"

#include <iostream>
#include <QString>
#include <QFile>
#include <string>
#include <vector>
using namespace std;


static const char* ERROR_SPECIES_NOT_FOUND = QT_TRANSLATE_NOOP("CSpecies", "Error! Species not found: ");
static const char* ERROR_SPECIES_KIND_NAME = QT_TRANSLATE_NOOP("CSpecies", "Species kind name not allowed: ");
static const char* ERROR_SPECIES_KIND_NUMBER = QT_TRANSLATE_NOOP("CSpecies", "Species kind number not allowed: ");


class CGlobalChemicalSystem;

/// \brief Enumeration containing a list of species kind
enum eSpeciesKind
{
    eAqueous = 0,                ///< 0: Aqueous
    eMineral = 1,                ///< 1: Mineral
    eGas = 2,                    ///< 2: Gas
    eSurface = 3,                ///< 3: Surface
    eLastSpeciesKind = 4                ///< 4: Exception
};


class CSpecies : public CCheprooBase
{
public:

    eSpeciesKind mSpeciesKind; // Type of the species
  
    double mIonicCharge;

    double mIonicRadius;

    double mMolWeight;

    double mLimMolConduct;

    double mMolVolume;

    double mCEC;

    double mElectricCharge;

    bool mIsConstActivity;          // True if the species has constant activity

    CGlobalChemicalSystem* mGlobChemSysPointer;

private:

    QString mDatabaseName;


public:

    void Read(		
        const QDomElement aNode ///< XML node with attributes
                      );

    /// \brief Constructor.
	CSpecies(QString aSpeciesType);

    /// \brief Constructor from a XML file.
	CSpecies(
		   QDomElement aNode ///< XML node with attributes
		   );

    /// \brief Converts a string into a species kind
	static eSpeciesKind str2SpeciesKind(QString s   ///< String that has to be converted into a time target
                                        );
    /// \brief Converts a species kind into a string
	static QString SpeciesKind2str(eSpeciesKind e ///< Time target that has to be converted into a string
                                   );

    /// \brief Set the thermodynamic database name
    void SetDatabase(QString aDatabaseName);

    /// \brief Destructor.
    ~CSpecies();

private:



};

#endif
