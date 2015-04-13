#ifndef CSURFACECEXCH_H
#define CSURFACECEXCH_H

#include "csurface.h"
#include "cglobalchemicalsystem.h"

#include <vector>

class CSurfaceCExchange : public CSurface
{

public:

    /// \brief Constructor.
	CSurfaceCExchange();

    /// \brief Constructor from a XML file.
	CSurfaceCExchange(
		   QDomElement aNode ///< XML node with attributes
		   );

    /// \brief Destructor.
    ~CSurfaceCExchange();


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
