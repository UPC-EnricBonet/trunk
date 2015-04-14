#ifndef CGAS_H
#define CGAS_H

#include "cphase.h"
#include "cglobalchemicalsystem.h"

#include <string>
#include <vector>
using namespace std;

#define Rgas 8.314472  // gas constant (Pa * m3)/(K * mol)

class CGas : public CPhase
{
public:

    CPhase* mPhasesFromGasPointer;

    double mPressure;

    double mVolume;

    double mDensity;

    double mViscosity;

    double mDensityDerivWrtTemp;

    double mDensityDerivWrtPress;

    double mDensityDerivWrtComp;

    double mViscosityDerivWrtTemp;

    double mViscosityDerivWrtPress;

    double mViscosityDerivWrtComp;


public:

    void ComputeActivityCoeff(vector<double> &aConcVector, 
                              vector<double> &aActCoeffsVector,
                              vector<CSpecies*> &aSV1,
                              vector<CSpecies*> &aSV2,
                              map<CSpecies*, int> &mSpeciesIndices,
                              double mIonicStrength=0.0);

    void ComputeActivityCoeffDerivs(vector<double> &aSVVector, 
                                    MatrixXd &MatrixToContrib,
                                    MatrixXd &aS1Mod,
                                    vector<CSpecies*> &aSV1,
                                    vector<CSpecies*> &aSV2,
                                    map<CSpecies*, int> &mSpeciesIndices,
                                    bool isContribToJacobian,
                                    bool isForRISA = false);

    void ComputeMolarity(vector<double> &c,
                         vector<CSpecies*> &aSV1,
                         vector<CSpecies*> &aSV2,
                         map<CSpecies*, int> &mSpeciesIndices,
                         CChemicalComposition* aChemComp);

    void ComputeConcFromMolarity(vector<double> &c,               // Vector of molarity in all phases
                         vector<CSpecies*> &aSV1,
                         vector<CSpecies*> &aSV2,
                         map<CSpecies*, int> &mSpeciesIndices,
                         CChemicalComposition* aChemComp);


    void Set(CChemicalComposition* aChemComp);

    void ComputeDensityDerivs(){};

    void ComputeViscosityDerivs(){};

    /// \brief Constructor.
	CGas();

    /// \brief Constructor from a XML file.
	CGas(
		   QDomElement aNode ///< XML node with attributes
		   );

    /// \brief Destructor.
	~CGas();

protected:

    virtual void pComputeDensity(double &density, CChemicalComposition* aChemComp, bool mustCalcDerivs){};
    virtual void pComputeViscosity(bool mustCalcDerivs){};
    virtual void pComputeDensityDerivs(){};
    virtual void pComputeViscosityDerivs(){};
    
    virtual void Read(
        const QDomElement aNode ///< XML node with attributes
              ); 

    virtual void pComputeFugacityCoeff(vector<double> &aConcVector, 
                                       vector<double> &aActCoeffsVector,
                                       vector<CSpecies*> &aSV1,
                                       vector<CSpecies*> &aSV2,
                                       map<CSpecies*, int> &mSpeciesIndices,
                                       double &mIonicStrength){};  // Specific function that evaluates the activity coefficients for each type of aqueous phase. 
                                                 // N.B. the size of the actCoeffVector is the number of species in the specific phase!!! 
                                                 // (That means that we can ask the evaluate them FOR phase, separately)

    virtual void pComputeFugacityCoeffDerivsWrtSV(vector<double> &aSVVector, 
                                                  MatrixXd &MatrixToContrib,
                                                  MatrixXd &aS1Mod,
                                                  vector<CSpecies*> &aSV1,
                                                  vector<CSpecies*> &aSV2,
                                                  map<CSpecies*, int> &mSpeciesIndices,
                                                  bool isContribToJacobian,
                                                  bool isForRISA = false){};

    virtual void pSet(CChemicalComposition* aChemComp){}; // Function that set the phase attributes

    virtual void pComputeMolarity(vector<double> &c,
                                  vector<CSpecies*> &aSV1,
                                  vector<CSpecies*> &aSV2,
                                  map<CSpecies*, int> &mSpeciesIndices,
                                  CChemicalComposition* aChemComp){};

    virtual void pComputeConcFromMolarity(vector<double> &c,
                                  vector<CSpecies*> &aSV1,
                                  vector<CSpecies*> &aSV2,
                                  map<CSpecies*, int> &mSpeciesIndices,
                                  CChemicalComposition* aChemComp){};
};

#endif
