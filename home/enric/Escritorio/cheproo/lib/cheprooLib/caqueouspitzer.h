#ifndef CAQUEOUSPITZER_H
#define CAQUEOUSPITZER_H

#include "caqueous.h"
#include "cglobalchemicalsystem.h"

#include <vector>
#include <math.h>


class CAqueousPitzer : public CAqueous
{
public:


    /// \brief Constructor.
	CAqueousPitzer();

    /// \brief Constructor from a XML file.
	CAqueousPitzer(
		   QDomElement aNode ///< XML node with attributes
		   );

    /// \brief Destructor.
	~CAqueousPitzer();


protected: 

    void pComputeActivityCoeff(vector<double> &aConcVector, 
                               vector<double> &aActCoeffsVector,
                               vector<CSpecies*> &aSV1,
                               vector<CSpecies*> &aSV2,
                               map<CSpecies*, int> &mSpeciesIndices,
                               double &mIonicStrength
                               );  // Specific function that evaluates the activity coefficients for each type of aqueous phase. 
                                   // The dimension of the vector is the number of species that belong to one ChemicalComposition

    void pComputeActivityCoeffDerivsWrtSV(vector<double> &aSVVector, 
                                          MatrixXd &MatrixToContrib,
                                          MatrixXd &aS1Mod,
                                          vector<CSpecies*> &aSV1,
                                          vector<CSpecies*> &aSV2,
                                          map<CSpecies*, int> &mSpeciesIndices,
                                          bool isContribToJacobian,
                                          bool isForRISA = false); 

    void pRead(        
        const QDomElement aNode ///< XML node with attributes
              );

    void ComputeZ(vector<double> &aConcVector, // Vector of molalities of all species (in all phases)
                   double &aZ,                  // Z value
                   vector<CSpecies*> &aSV1,
                   vector<CSpecies*> &aSV2,
                   map<CSpecies*, int> &mSpeciesIndices
                   );

    void ComputeM(vector<double> &aConcVector, // Vector of molalities of all species (in all phases)
                   double &aM) ;                 // M value
 


};

#endif
