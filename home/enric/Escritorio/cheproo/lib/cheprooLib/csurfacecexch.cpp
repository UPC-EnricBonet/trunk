#include "csurfacecexch.h"

CSurfaceCExchange::CSurfaceCExchange()
{
    this->mClassName = "CMineralPure";

}

CSurfaceCExchange::CSurfaceCExchange(QDomElement aNode)
: CSurface(aNode)
{
    this->mClassName = "CSurfaceCExchange";

}

CSurfaceCExchange::~CSurfaceCExchange()
{
}

void
CSurfaceCExchange::pComputeActivityCoeff(vector<double> &aConcVector, vector<double> &aActCoeffsVector, vector<CSpecies*> &aSV1,
                              vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, double &mIonicStrength)
{


}

void 
CSurfaceCExchange::pComputeActivityCoeffDerivsWrtSV(vector<double> &aSVVector, MatrixXd &MatrixToContrib, MatrixXd &aS1Mod,
                                                vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices,
                                                bool isContribToJacobian, bool isForRISA)
{

}

void 
CSurfaceCExchange::Read(const QDomElement aNode)
{
    CSurface::Read(aNode);
}

void 
CSurfaceCExchange::pComputeMolarity(vector<double> &c, vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, CChemicalComposition* aChemComp)
{
    // First loop over primary species
    for(int i=0; i!=aSV1.size(); i++)        
    {
        c[mSpeciesIndices[aSV1[i]]] = c[mSpeciesIndices[aSV1[i]]] * aSV1[i]->mCEC;
    }

    // Second loop over secondary species
    for(int i=0; i!=aSV2.size(); i++)        
    {
        c[mSpeciesIndices[aSV2[i]]] = c[mSpeciesIndices[aSV1[i]]] * aSV1[i]->mCEC;
    }
}

void 
CSurfaceCExchange::pComputeConcFromMolarity(vector<double> &c, vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, CChemicalComposition* aChemComp)
{
    // First loop over primary species
    for(int i=0; i!=aSV1.size(); i++)        
    {
        c[mSpeciesIndices[aSV1[i]]] = c[mSpeciesIndices[aSV1[i]]] / aSV1[i]->mCEC;
    }

    // Second loop over secondary species
    for(int i=0; i!=aSV2.size(); i++)        
    {
        c[mSpeciesIndices[aSV2[i]]] = c[mSpeciesIndices[aSV1[i]]] / aSV1[i]->mCEC;
    }
}