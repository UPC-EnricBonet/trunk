#include "cmineralpure.h"

CMineralPure::CMineralPure()
{
    this->mClassName = "CMineralPure";

}

CMineralPure::CMineralPure(QDomElement aNode)
: CMineral(aNode)
{
    this->mClassName = "CMineralPure";

}

CMineralPure::~CMineralPure()
{
}

void
CMineralPure::pComputeActivityCoeff(vector<double> &aConcVector, vector<double> &aActCoeffsVector, vector<CSpecies*> &aSV1,
                              vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, double &mIonicStrength)
{
    if(aSV1.size()!=0)
    {
        for(int i=0; i!=aSV1.size(); i++)
        {
            aActCoeffsVector.at(mSpeciesIndices[aSV1[i]]) = 1.0;
        }
    }

    if(aSV2.size()!=0)
    {
        for(int i=0; i!=aSV2.size(); i++)
        {
            aActCoeffsVector.at(mSpeciesIndices[aSV2[i]]) = 1.0;
        }
    }
}

void 
CMineralPure::pComputeActivityCoeffDerivsWrtSV(vector<double> &aSVVector, MatrixXd &MatrixToContrib, MatrixXd &aS1Mod,
                                                vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices,
                                                bool isContribToJacobian, bool isForRISA)
{

}

void 
CMineralPure::Read(const QDomElement aNode)
{
    CMineral::Read(aNode);
}

void 
CMineralPure::pComputeMolarity(vector<double> &c, vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, CChemicalComposition* aChemComp)
{
    // First loop over primary species
    for(int i=0; i!=aSV1.size(); i++)        
    {
        c[mSpeciesIndices[aSV1[i]]] = 1 / aSV1[i]->mMolVolume / 1e-6;
    }

    // Second loop over secondary species
    for(int i=0; i!=aSV2.size(); i++)        
    {
        c[mSpeciesIndices[aSV2[i]]] = 1 / aSV2[i]->mMolVolume / 1e-6;
    }
}

void 
CMineralPure::pComputeConcFromMolarity(vector<double> &c, vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, CChemicalComposition* aChemComp)
{
    // First loop over primary species
    for(int i=0; i!=aSV1.size(); i++)        
    {
        c[mSpeciesIndices[aSV1[i]]] = aSV1[i]->mMolVolume * 1e-6;
    }

    // Second loop over secondary species
    for(int i=0; i!=aSV2.size(); i++)        
    {
        c[mSpeciesIndices[aSV2[i]]] = aSV2[i]->mMolVolume * 1e-6;
    }
}