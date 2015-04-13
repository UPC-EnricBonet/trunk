#ifndef CPHASE_H
#define CPHASE_H

#include "ccheproobase.h"
#include "classnameconstants.h"
#include "objectmanager.h"
#include "xmlconstants.h"
#include "xmlqt.h"
//#include "proostglobals.h"

#include "cspecies.h"
#include "cchemicalcomposition.h"

#include <QString>
#include <string>
#include <vector>
#include <map>
#include <math.h>
using namespace std;

#include <Eigen/Dense>
using namespace Eigen;

class CGlobalChemicalSystem;
class CLocalChemicalSystem;
class CChemicalComposition;

static const char* ERROR_PHASE_KIND_NAME = QT_TRANSLATE_NOOP("CPhase", "Phase kind name not allowed");
static const char* ERROR_PHASE_KIND_NUMBER = QT_TRANSLATE_NOOP("CPhase", "Phase kind number not allowed");

/// \brief Enumeration containing a list of phases kind 
enum ePhaseKind
{
    eAqueousPhase = 0,                ///< 0: Aqueous 
    eMineralPhase = 1,                ///< 1: Mineral 
    eGasPhase = 2,                    ///< 2: Gas 
    eSurfacePhase = 3,                ///< 3: Surface
    eLastPhaseKind = 4           ///< 4: Exception
};


class CPhase : public CCheprooBase
{
public:

    ePhaseKind mPhaseKind;  // Type of the phase (aqueous, mineral, gas or surface)

    CGlobalChemicalSystem* mGlobChemSysPointer; // Pointer to the GlobalChemicalSystem

	vector<string> mSpeciesNames;

public:
    
    virtual void Read(		
        const QDomElement aNode ///< XML node with attributes
                      );  // Virtual function because the attributes of each phase are different and must be read in a different way

    void ComputeDensity(double &density, CChemicalComposition* aChemComp, bool mustCalcDerivs);

    void ComputeViscosity();

    void ComputeDensityDerivs(){};

    void ComputeViscosityDerivs(){};

    void ComputeIonicStrength(vector<double> &aConcVector, // Vector of molalities of all species (in all phases)
                              double &mIonicStrength,       // Ionic strength
                              vector<CSpecies*> &aSV1,
                              vector<CSpecies*> &aSV2,
                              map<CSpecies*, int> &mSpeciesIndices
                              );

    void ComputeIonicStrengthDerivsWrtConc(vector<double> &aIonicStrengthDerivs,  // Vector of derivatives of ionic strength wrt concentrations
                                           vector<CSpecies*> &aSV1,
                                           vector<CSpecies*> &aSV2,
                                           map<CSpecies*, int> &mSpeciesIndices
                                           );

    virtual void ComputeActivityCoeff(vector<double> &aConcVector, 
                                      vector<double> &aActCoeffsVector,
                                      vector<CSpecies*> &aSV1,
                                      vector<CSpecies*> &aSV2,
                                      map<CSpecies*, int> &mSpeciesIndices,
                                      double mIonicStrength = 0.0){}; // N.B. the size of the actCoeffVector is the number of species in the specific phase!!! 
                                                                          // (That means that we can ask the evaluate them FOR phase, separately)
    virtual void ComputeActivityCoeffDerivs(vector<double> &aSVVector, 
                                            MatrixXd &MatrixToContrib,
                                            MatrixXd &aS1Mod,
                                            vector<CSpecies*> &aSV1,
                                            vector<CSpecies*> &aSV2,
                                            map<CSpecies*, int> &mSpeciesIndices,
                                            bool isContribToJacobian,
                                            bool isForRISA = false){}; 

    void ComputeChargeBalance(vector<double> &aConcVector,
                              double &aChargeBalance,
                              vector<CSpecies*> &aSV1,
                              vector<CSpecies*> &aSV2,
                              map<CSpecies*, int> &mSpeciesIndices);

    void ComputeChargeBalanceDerivs(MatrixXd &dc2_dc1,
                                    VectorXd &aChargeBalanceDerivs,
                                    vector<CSpecies*> &aSV1,
                                    vector<CSpecies*> &aSV2,
                                    map<CSpecies*, int> &mSpeciesIndices);

    virtual void ComputeMolarity(vector<double> &c,               // Vector of concentration in all phases
                                 vector<CSpecies*> &aSV1,
                                 vector<CSpecies*> &aSV2,
                                 map<CSpecies*, int> &mSpeciesIndices,
                                 CChemicalComposition* aChemComp){};

    virtual void ComputeConcFromMolarity(vector<double> &c,               // Vector of molarity in all phases
                                 vector<CSpecies*> &aSV1,
                                 vector<CSpecies*> &aSV2,
                                 map<CSpecies*, int> &mSpeciesIndices,
                                 CChemicalComposition* aChemComp){};

    /// \brief Set the attributes of the phase, such as the coefficients to evaluate the activity coefficients
    virtual void Set(CChemicalComposition* aChemComp){};
        
    /// \brief Converts a string into a phase kind
	static ePhaseKind str2PhasesKind(QString s   ///< String that has to be converted into a time target
                                        );
    /// \brief Converts a phase kind into a string
	static QString phasesKind2str(ePhaseKind e ///< Time target that has to be converted into a string
                                   );


    /// \brief Constructor.
	CPhase();

    /// \brief Constructor from a XML file.
	CPhase(
		   QDomElement aNode ///< XML node with attributes
		   );

    /// \brief Destructor.
	~CPhase();

protected:

    virtual void pComputeDensity(
            double &density,
            CChemicalComposition* aChemComp,
            bool mMustCalcDerivs  ///<If TRUE the derivatives of the properties wrt to temperature, pressure and composition must be evaluated
            ){};

    virtual void pComputeViscosity(
            bool mustCalcDerivs  ///<If TRUE the derivatives of the properties wrt to temperature, pressure and composition must be evaluated
            ){};

    virtual void pComputeDensityDerivs(){};

    virtual void pComputeViscosityDerivs(){};

    virtual void pComputeIonicStrength(vector<double> &aConcVector, double &mIonicStrength, vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2,
                                       map<CSpecies*, int> &mSpeciesIndices){};


    virtual void pComputeIonicStrengthDerivsWrtConc(vector<double> &aIonicStrengthDerivs,
                                                    vector<CSpecies*> &aSV1,
                                                    vector<CSpecies*> &aSV2,
                                                    map<CSpecies*, int> &mSpeciesIndices
                                                    ){};

    virtual void pComputeChargeBalance(vector<double> &aConcVector,
                                       double &aChargeBalance,
                                       vector<CSpecies*> &aSV1,
                                       vector<CSpecies*> &aSV2,
                                       map<CSpecies*, int> &mSpeciesIndices){};

    virtual void pComputeChargeBalanceDerivsWrtConc(MatrixXd &dc2_dc1,
                                                    VectorXd &aChargeBalanceDerivs,
                                                    vector<CSpecies*> &aSV1,
                                                    vector<CSpecies*> &aSV2,
                                                    map<CSpecies*, int> &mSpeciesIndices){};


};

#endif
