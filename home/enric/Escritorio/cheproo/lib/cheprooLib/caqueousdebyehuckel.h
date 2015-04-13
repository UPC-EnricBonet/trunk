#ifndef CAQUEOUSDEBYEHUCKEL_H
#define CAQUEOUSDEBYEHUCKEL_H

#include "caqueous.h"
#include "cglobalchemicalsystem.h"

#include <vector>
#include <math.h>

#define ALT1 4.910300000e-01
#define ALT2 5.808571429e-04
#define ALT3 5.527142857e-06
#define ALT4 -4.857142857e-09
#define ALT5 0.000000000e+00
#define AHT1 6.440000000e-01
#define AHT2 -3.436166667e-03
#define AHT3 4.408833333e-05
#define AHT4 -1.691333333e-07
#define AHT5 2.766666667e-10
#define BLT1 3.247000000e-01
#define BLT2 1.309285714e-04
#define BLT3 5.502380952e-07
#define BLT4 -1.095238095e-09
#define BLT5 0.000000000e+00
#define BHT1 3.302000000e-01
#define BHT2 -1.650000000e-05
#define BHT3 1.991666667e-06
#define BHT4 -7.400000000e-09
#define BHT5 1.133333333e-11
#define DLT1 1.740000000e-02
#define DLT2 1.509047619e-03
#define DLT3 -2.605904762e-05
#define DLT4 1.382857143e-07
#define DLT5 0.000000000e+00
#define DHT1 1.090000000e-01
#define DHT2 -1.483333333e-03
#define DHT3 1.173333333e-05
#define DHT4 -3.466666667e-08
#define DHT5 2.666666667e-11

class CAqueousDebyeHuckel : public CAqueous
{
public:

	double A;           // Coefficient A temperature depending

    double B;           // Coefficient B temperature depending

    double Bdot;        // Coefficient Bdot temperature depending

    /// \brief Constructor.
	CAqueousDebyeHuckel();

    /// \brief Constructor from a XML file.
	CAqueousDebyeHuckel(
		   QDomElement aNode ///< XML node with attributes
		   );

    /// \brief Destructor.
	~CAqueousDebyeHuckel();


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

    void pComputeActivityCoeffDerivs(CChemicalComposition* aChemComp,
                                    vector<CSpecies*> &aSV1,
                                    vector<CSpecies*> &aSV2,
                                    map<CSpecies*, int> &mSpeciesIndices, 
                                    VectorXd aDerivatives,
                                    bool isWrtTemp = true);

    void pSet(CChemicalComposition* aChemComp);            // Function that set the phase attributes

    void pComputeCoeffAB(CChemicalComposition* aChemComp); // Function that evaluates the coefficients A and B

    void pComputeCoeffABDerivsWrtT(CChemicalComposition* aChemComp, // Function that evaluates coefficients A and B derivatives wrt temperature
                                   double &aA_deriv,
                                   double &aB_deriv,
                                   double &aBdot_deriv); 

    void pRead(        
        const QDomElement aNode ///< XML node with attributes
              ); 


};

#endif
