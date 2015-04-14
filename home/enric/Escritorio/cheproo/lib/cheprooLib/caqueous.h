#ifndef CAQUEOUS_H
#define CAQUEOUS_H

#include "cphase.h"
#include "cglobalchemicalsystem.h"

#include <vector>

/// \brief Enumeration containing a list of density models
enum eDensityKind
{
    eDensityCB1 = 0,             ///< 0: Code Bright model 1
    eDensityCB2 = 1,             ///< 1: Code Bright model 2
    eDensityDefault = 2           ///< 2: Exception
};

class CAqueous : public CPhase
{
public:

    eDensityKind mDensityKind;  // Model for density

private:

    double m_rho_ref;
    double mCompressibility;
    double mVolumThermCoeff;
    double mSoluteVar;
    double mRefPressure;

public:

    void ComputeActivityCoeff(vector<double> &aConcVector, 
                              vector<double> &aActCoeffsVector,
                              vector<CSpecies*> &aSV1,
                              vector<CSpecies*> &aSV2,
                              map<CSpecies*, int> &mSpeciesIndices,
                              double mIonicStrength);

    void ComputeActivityCoeffDerivs(vector<double> &aSVVector, 
                                    MatrixXd &MatrixToContrib,
                                    MatrixXd &aS1Mod,
                                    vector<CSpecies*> &aSV1,
                                    vector<CSpecies*> &aSV2,
                                    map<CSpecies*, int> &mSpeciesIndices,
                                    bool isContribToJacobian,
                                    bool isForRISA = false);

    void ComputeActivityCoeffDerivs(CChemicalComposition* aChemComp,
                                    vector<CSpecies*> &aSV1,
                                    vector<CSpecies*> &aSV2,
                                    map<CSpecies*, int> &mSpeciesIndices, 
                                    VectorXd aDerivatives,
                                    bool isWrtTemp = true);

    void Set(CChemicalComposition* aChemComp);

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

    void ComputeViscosity(){};

    /// \brief Compute water content in liquid 
    double ComputeWH2OLiq(vector<double> &c,
                         vector<CSpecies*> &aSV1,
                         vector<CSpecies*> &aSV2,
                         map<CSpecies*, int> &mSpeciesIndices);

    /// \brief Constructor.
	CAqueous();

    /// \brief Constructor from a XML file.
	CAqueous(
		   QDomElement aNode ///< XML node with attributes
		   );

    /// \brief Destructor.
    ~CAqueous();


protected: 


    void pComputeIonicStrength(vector<double> &aConcVector, // Vector of molalities of all species (in all phases)
                               double &mIonicStrength,       // Ionic strength
                               vector<CSpecies*> &aSV1,
                               vector<CSpecies*> &aSV2,
                               map<CSpecies*, int> &mSpeciesIndices
                               );
 
    void pComputeIonicStrengthDerivsWrtConc(vector<double> &aIonicStrengthDerivs,  // Vector of derivatives of ionic strength wrt concentrations
                                            vector<CSpecies*> &aSV1,
                                            vector<CSpecies*> &aSV2,
                                            map<CSpecies*, int> &mSpeciesIndices
                                            );


    virtual void pComputeActivityCoeff(vector<double> &aConcVector, 
                                       vector<double> &aActCoeffsVector,
                                       vector<CSpecies*> &aSV1,
                                       vector<CSpecies*> &aSV2,
                                       map<CSpecies*, int> &mSpeciesIndices,
                                       double &mIonicStrength){};  // Specific function that evaluates the activity coefficients for each type of aqueous phase. 
                                                 // N.B. the size of the actCoeffVector is the number of species in the specific phase!!! 
                                                 // (That means that we can ask the evaluate them FOR phase, separately)

    virtual void pComputeActivityCoeffDerivsWrtSV(vector<double> &aSVVector, 
                                                  MatrixXd &MatrixToContrib,
                                                  MatrixXd &aS1Mod,
                                                  vector<CSpecies*> &aSV1,
                                                  vector<CSpecies*> &aSV2,
                                                  map<CSpecies*, int> &mSpeciesIndices,
                                                  bool isContribToJacobian,
                                                  bool isForRISA = false){}; 

    virtual void pComputeActivityCoeffDerivs(CChemicalComposition* aChemComp,
                                             vector<CSpecies*> &aSV1,
                                             vector<CSpecies*> &aSV2,
                                             map<CSpecies*, int> &mSpeciesIndices, 
                                             VectorXd aDerivatives,
                                             bool isWrtTemp = true){};  // Method that evaluates dLnGamma/dT

    void pComputeChargeBalance(vector<double> &aConcVector,
                               double &aChargeBalance,
                               vector<CSpecies*> &aSV1,
                               vector<CSpecies*> &aSV2,
                               map<CSpecies*, int> &mSpeciesIndices);

    void pComputeChargeBalanceDerivsWrtConc(MatrixXd &dc2_dc1,
                                            VectorXd &aChargeBalanceDerivs,
                                            vector<CSpecies*> &aSV1,
                                            vector<CSpecies*> &aSV2,
                                            map<CSpecies*, int> &mSpeciesIndices); // Compute the derivative of the charge balance wrt primary species, 
                                                                                   // which contributes to the Jacobian matrix for the evaluation of c1

    virtual void pSet(CChemicalComposition* aChemComp){}; // Function that set the phase attributes

    virtual void Read(
        const QDomElement aNode ///< XML node with attributes
              ); 

    void pComputeDensity(double &density, CChemicalComposition* aChemComp, bool mustCalcDerivs);

};

#endif
