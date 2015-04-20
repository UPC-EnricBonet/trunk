#ifndef CLOCALCHEMICALSYSTEMCONSTPRIMARY_H
#define CLOCALCHEMICALSYSTEMCONSTPRIMARY_H

#include "ccheproobase.h"

#include "clocalchemicalsystem.h"

#include <QString>
#include <string>
#include <vector>
using namespace std;

/*! \brief This class represent a chemical system where minerals are treated as primary species and eliminated
is applied.
*/

class CLocalChemicalSystemConstPrimary : public CLocalChemicalSystem
{
private:



public:

    void SetLocal(vector<string>& primSpecies, set<string>& cas);

    /// \brief Constructor from a XML file.
	CLocalChemicalSystemConstPrimary(
		   QDomElement aNode ///< XML node with attributes
           );


    void ComputeS1Mod();

    void ComputeComponentMatrix();

    MatrixXd Get_dc_dc1_ByPhase(QString aVarType);

    CLocalChemicalSystemConstPrimary();
	~CLocalChemicalSystemConstPrimary();




private:

    MatrixXd ComputeComponentMatrixInit();

	MatrixXd GetComponentInit(map<CSpecies*,int>& tempIndices);


};

#endif
