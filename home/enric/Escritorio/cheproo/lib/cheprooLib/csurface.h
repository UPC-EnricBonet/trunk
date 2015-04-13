#ifndef CSURFACE_H
#define CSURFACE_H

#include "cphase.h"
#include "cglobalchemicalsystem.h"

#include <vector>

class CSurface : public CPhase
{
public:


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

    void Set(CChemicalComposition* aChemComp);


    void ComputeMolarity(vector<double> &c,
                         vector<CSpecies*> &aSV1,
                         vector<CSpecies*> &aSV2,
                         map<CSpecies*, int> &mSpeciesIndices,
                         CChemicalComposition* aChemComp);

    void ComputeConcFromMolarity(vector<double> &c,
                         vector<CSpecies*> &aSV1,
                         vector<CSpecies*> &aSV2,
                         map<CSpecies*, int> &mSpeciesIndices,
                         CChemicalComposition* aChemComp);

    void ComputeDensity(){};

    void ComputeViscosity(){};

    /// \brief Constructor.
	CSurface();

    /// \brief Constructor from a XML file.
	CSurface(
		   QDomElement aNode ///< XML node with attributes
		   );

    /// \brief Destructor.
    ~CSurface();


protected: 

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

    virtual void pSet(CChemicalComposition* aChemComp){}; // Function that set the phase attributes


    virtual void Read(
        const QDomElement aNode ///< XML node with attributes
              ); 

    virtual void pComputeDensity(double &density, CChemicalComposition* aChemComp, bool mustCalcDerivs){};

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
