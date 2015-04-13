#ifndef CAQUEOUSIDEAL_H
#define CAQUEOUSIDEAL_H

#include "caqueous.h"
#include "cglobalchemicalsystem.h"

#include <vector>
#include <math.h>

class CAqueousIdeal : public CAqueous
{
public:

    /// \brief Constructor.
	CAqueousIdeal();

    /// \brief Constructor from a XML file.
	CAqueousIdeal(
		   QDomElement aNode ///< XML node with attributes
		   );

    /// \brief Destructor.
	~CAqueousIdeal();


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

    void pComputeActivityCoeffDerivs(CChemicalComposition* aChemComp,
                                    vector<CSpecies*> &aSV1,
                                    vector<CSpecies*> &aSV2,
                                    map<CSpecies*, int> &mSpeciesIndices, 
                                    VectorXd aDerivatives,
                                    bool isWrtTemp = true);


    void pRead(        
        const QDomElement aNode ///< XML node with attributes
              ); 


};

#endif
