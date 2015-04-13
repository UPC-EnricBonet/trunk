#ifndef CGASBINARYMIXTURE_H
#define CGASBINARYMIXTURE_H

#include "cgas.h"
#include "cglobalchemicalsystem.h"

#include <string>
#include <vector>
using namespace std;


/*!
\brief This class represents a mixture of CO2 and water. 
The termodynamic properties are evaluated as in Spicher et al.(2003),
by means of a third order equation of state.

*/

class CGasBinaryMixture : public CGas
{
public:

    CGas* mGasFromBinMixPointer;

public:

    CGasBinaryMixture();
	~CGasBinaryMixture();

    // \brief Constructor from a XML file.
	CGasBinaryMixture(
		QDomElement aNode ///< XML node with attributes
        ){}


protected:

    void pComputeFugacityCoeff(vector<double> &aConcVector, 
                               vector<double> &aActCoeffsVector, 
                               double &mIonicStrength){};  // Specific function that evaluates the activity coefficients for each type of aqueous phase. 
                                                 // N.B. the size of the actCoeffVector is the number of species in the specific phase!!! 
                                                 // (That means that we can ask the evaluate them FOR phase, separately)

    void pComputeFugacityCoeffDerivsWrtSV(vector<double> &aSVVector, 
                                            vector<vector<double> > &aActCoeffsDerivWrtSV,
                                            MatrixXd &aS1Mod,
                                            vector<int> &aSV1GlobalIndices,
                                            vector<int> &aSV2GlobalIndices){}; 



    void pComputeDensity(double &density, CChemicalComposition* aChemComp, bool mustCalcDerivs);

    void pComputeViscosity(bool mustCalcDerivs);

    void ComputeDensityDerivs(){};

    void ComputeViscosityDerivs(){};
    
    void pRead(        
        const QDomElement aNode ///< XML node with attributes
              ); 



private:

    double mMolVol;  // Molar volume of the phase

    double mMixtureMolWeight;  // Molar weight of the mixture

    double Amix; // Parameter representing molecular attraction

    double Bmix; // Parameter representing molecular repulsion


};

#endif
