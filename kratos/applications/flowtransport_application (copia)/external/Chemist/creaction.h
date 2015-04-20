#ifndef CREACTION_H
#define CREACTION_H

#include "ccheproobase.h"

#include "cspecies.h"
#include "creactionratelaw.h"

#include <QString>
#include <string>
#include <vector>
#include <map>
#include <math.h>
using namespace std;

static const char* ERROR_PRIMARY_SPECIES_NUM = QT_TRANSLATE_NOOP("CReaction", "Error: Dimension of primary species vector not correct");

class CGlobalChemicalSystem;

class CReaction : public CCheprooBase
{
public:

    CGlobalChemicalSystem* mGlobChemSysPointer; // Pointer to GlobalChemicalSystem

    CReactionRateLaw* mReacRateLawPointer;  // Pointer to reaction rate law class

    map<CSpecies*, double> mRelSpeciesVsStoichCoeffs;  // Map of species involved in the reaction vs corresponding stoichiometric coefficients

	map<double, double> mTempVsLogK; // Map containing the values of the equilibrium constant at several temperatures

    VectorXd mCoeffsLogK;  // Vector containing the coefficients for the calculation of logK at different temp values

    double mLogK;

    bool mIsKinetic;

private:

    QString mDatabaseName;



public:

    void Read(		
        const QDomElement aNode, ///< XML node with attributes
        bool isKinetic = false
                      );

    /// brief Evaluates concentration of a secondary species solving the Mass Action Law (without adsorbed species)
    void SolveMALExplicit(double &aConcValue,                               // Value of concentration that returns this method
                          int &aIndex,                                      // Index of the species we want to evaluate the concentration
					      vector<double> &aConcVector,                       // Vector of concentrations of species
						  vector<double> &aActCoeffsVector,                  // Vector of activity coefficients of primary species
                          double aTemp,                                      // Temperature value
						  map<CSpecies*, int> mSpeciesIndices    // Map containing the indices of the species that belong to the LocalChemicalSystem
	                      );

    /// \brief Set the thermodynamic database name
    void SetDatabase(QString aDatabaseName);

    /// \brief Set the value of logK
    void SetLogK(double &temperature);

    /// \brief Constructor.
	CReaction();

    /// \brief Constructor from a XML file.
	CReaction(
		   QDomElement aNode ///< XML node with attributes
		   );

    /// \brief Destructor.
    ~CReaction();
};

#endif
