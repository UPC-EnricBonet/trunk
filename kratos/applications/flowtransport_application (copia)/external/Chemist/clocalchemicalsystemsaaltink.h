#ifndef CLOCALCHEMICALSYSTEMSAALTINK_H
#define CLOCALCHEMICALSYSTEMSAALTINK_H

#include "ccheproobase.h"

#include "clocalchemicalsystem.h"

#include <QString>
#include <string>
#include <vector>
using namespace std;

/*! \brief This class represent a LocalChemicalSystem that is characterized by evaluating the elimination matrix to get rid of constant activity species
  
  The component matrix that is evaluated
  For further details: Saaltink et al. "A mathematical formulation for reactive transport", WRR Vol.34, 7, 1649-1656, 1998.

*/

class CLocalChemicalSystemSaaltink : public CLocalChemicalSystem
{
public:

    MatrixXd mEliminationMatrix;  // Elimination matrix E

public:

    void SetLocal(vector<string>& primSpecies, set<string>& cas);

    void ComputeS1Mod();

    void ComputeComponentMatrix();   // The method computes the component matrix as EU
    
    /// \brief Constructor from a XML file.
	CLocalChemicalSystemSaaltink(
		   QDomElement aNode ///< XML node with attributes
           );


    CLocalChemicalSystemSaaltink();
	~CLocalChemicalSystemSaaltink();


private:

    void ComputeEliminationMatrix();

    void UpdateLocal();

};

#endif
