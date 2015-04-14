#ifndef CMINERALPURE_H
#define CMINERALPURE_H

#include "cmineral.h"
#include "cglobalchemicalsystem.h"

#include <vector>

class CMineralPure : public CMineral
{

public:

    /// \brief Constructor.
	CMineralPure();

    /// \brief Constructor from a XML file.
	CMineralPure(
		   QDomElement aNode ///< XML node with attributes
		   );

    /// \brief Destructor.
    ~CMineralPure();


protected: 

    void pComputeActivityCoeff(vector<double> &aConcVector, 
                               vector<double> &aActCoeffsVector,
                               vector<CSpecies*> &aSV1,
                               vector<CSpecies*> &aSV2,
                               map<CSpecies*, int> &mSpeciesIndices,
                               double &mIonicStrength
                               );

    void pComputeActivityCoeffDerivsWrtSV(vector<double> &aSVVector, 
                                          MatrixXd &MatrixToContrib,
                                          MatrixXd &aS1Mod,
                                          vector<CSpecies*> &aSV1,
                                          vector<CSpecies*> &aSV2,
                                          map<CSpecies*, int> &mSpeciesIndices,
                                          bool isContribToJacobian,
                                          bool isForRISA = false); 

    void pComputeMolarity(vector<double> &c,
                          vector<CSpecies*> &aSV1,
                          vector<CSpecies*> &aSV2,
                          map<CSpecies*, int> &mSpeciesIndices,
                          CChemicalComposition* aChemComp);

    void pComputeConcFromMolarity(vector<double> &c,
                          vector<CSpecies*> &aSV1,
                          vector<CSpecies*> &aSV2,
                          map<CSpecies*, int> &mSpeciesIndices,
                          CChemicalComposition* aChemComp);

    void Read(
        const QDomElement aNode ///< XML node with attributes
              ); 


};

#endif
