#include "caqueousideal.h"


CAqueousIdeal::CAqueousIdeal()
{
    this->mClassName = "CAqueousIdeal";
}

CAqueousIdeal::CAqueousIdeal(QDomElement aNode)
: CAqueous(aNode)
{
    this->mClassName = "CAqueousIdeal";
}


CAqueousIdeal::~CAqueousIdeal()
{
}

void
CAqueousIdeal::pComputeActivityCoeff(vector<double> &aConcVector, vector<double> &aActCoeffsVector, vector<CSpecies*> &aSV1,
                              vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, double &mIonicStrength)
{
    for(int i=0; i!=aSV1.size(); i++)
    {
        aActCoeffsVector.at(mSpeciesIndices[aSV1[i]]) = 1.0;
    }

    for(int i=0; i!=aSV2.size(); i++)
    {
        aActCoeffsVector.at(mSpeciesIndices[aSV2[i]]) = 1.0;
    }
};


void 
CAqueousIdeal::pComputeActivityCoeffDerivsWrtSV(vector<double> &aSVVector, MatrixXd &MatrixToContrib, MatrixXd &aS1Mod,
                                          vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices,
                                          bool isContribToJacobian, bool isForRISA)
{
}


void 
CAqueousIdeal::pRead(const QDomElement aNode)
{
    CAqueous::Read(aNode);
}

void
CAqueousIdeal::pComputeActivityCoeffDerivs(CChemicalComposition* aChemComp, vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, 
                                           VectorXd aDerivatives, bool isWrtTemp)
{
   // Do nothing, derivatives are zero
}
