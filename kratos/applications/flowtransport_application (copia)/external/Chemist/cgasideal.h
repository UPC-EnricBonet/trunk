#ifndef CGASIDEAL_H
#define CGASIDEAL_H

#include "cgas.h"
#include "cglobalchemicalsystem.h"

#include <string>
#include <vector>
using namespace std;


/*!
\brief This class represents an ideal gas behavior

*/

class CGasIdeal : public CGas
{
public:


public:

    CGasIdeal();
	~CGasIdeal();

    // \brief Constructor from a XML file.
	CGasIdeal(
		QDomElement aNode ///< XML node with attributes
        ){}


protected:

    void pComputeFugacityCoeff(vector<double> &aConcVector, 
                               vector<double> &aActCoeffsVector,
                               vector<CSpecies*> &aSV1,
                               vector<CSpecies*> &aSV2,
                               map<CSpecies*, int> &mSpeciesIndices,
                               double &mIonicStrength);  // Specific function that evaluates the activity coefficients for each type of aqueous phase. 
                                                 // N.B. the size of the actCoeffVector is the number of species in the specific phase!!! 
                                                 // (That means that we can ask the evaluate them FOR phase, separately)

    void pComputeFugacityCoeffDerivsWrtSV(vector<double> &aSVVector, 
                                          MatrixXd &MatrixToContrib,
                                          MatrixXd &aS1Mod,
                                          vector<CSpecies*> &aSV1,
                                          vector<CSpecies*> &aSV2,
                                          map<CSpecies*, int> &mSpeciesIndices,
                                          bool isContribToJacobian,
                                          bool isForRISA = false); 

    void pComputeDensity(double &density, CChemicalComposition* aChemComp, bool mustCalcDerivs);

    void pComputeViscosity(bool mustCalcDerivs);

    void ComputeDensityDerivs(){};

    void ComputeViscosityDerivs(){};
    
    void Read(        
        const QDomElement aNode ///< XML node with attributes
              ); 

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

private:

    double mMolVol;  // Molar volume of the phase

};

#endif
