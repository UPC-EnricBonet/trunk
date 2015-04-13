#include "cgasideal.h"


CGasIdeal::CGasIdeal(void)
{
    this->mClassName = "CGasIdeal";

}

CGasIdeal::~CGasIdeal(void)
{
    this->mClassName = "CGasIdeal";

}

void
CGasIdeal::pComputeFugacityCoeff(vector<double> &aConcVector, vector<double> &aActCoeffsVector, vector<CSpecies*> &aSV1,
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
CGasIdeal::pComputeFugacityCoeffDerivsWrtSV(vector<double> &aSVVector, MatrixXd &MatrixToContrib, MatrixXd &aS1Mod,
                                                vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices,
                                                bool isContribToJacobian, bool isForRISA)
{
    //if(aS1Mod.size()!=0 && aSV2.size()!=0)
    //{
    //    if (isContribToJacobian) // I'm evaluating contribution to Jacobian
    //    {
    //        // First loop over secondary species
    //        for(int i=0; i!=aSV2.size(); i++)
    //        {
    //            // Second loop over secondary species
    //            for(int j=0; j!=aSV2.size(); j++)
    //            {
    //                MatrixToContrib(mSpeciesIndices[aSV2[i]], mSpeciesIndices[aSV2[j]]) = 0.0;
    //            }
    //        }
    //    }
    // }
    //else  // I'm evaluating contribution to residual
    //{
    //    // First loop over secondary species
    //    for(int i=0; i!=aSV2.size(); i++)
    //    {
    //        // Second loop over primary species
    //        for(int j=0; j!=aSV1.size(); j++)
    //        {
    //            MatrixToContrib(mSpeciesIndices[aSV2[i]], mSpeciesIndices[aSV1[j]]) = 0.0;
    //        }

    //    }
    //}

}

void
CGasIdeal::pComputeDensity(double &density, CChemicalComposition* aChemComp, bool mustCalcDerivs)
{
}

void
CGasIdeal::pComputeViscosity(bool mustCalcDerivs)
{
}

void 
CGasIdeal::Read(const QDomElement aNode)
{
    CGas::Read(aNode);
}

void 
CGasIdeal::pComputeMolarity(vector<double> &c, vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, CChemicalComposition* aChemComp)
{
    // First loop over primary species
    for(int i=0; i!=aSV1.size(); i++)        
    {
        c[mSpeciesIndices[aSV1[i]]] = c[mSpeciesIndices[aSV1[i]]] / Rgas / aChemComp->mTemp;
    }

    // Second loop over secondary species
    for(int i=0; i!=aSV2.size(); i++)        
    {
        c[mSpeciesIndices[aSV2[i]]] = c[mSpeciesIndices[aSV2[i]]] / Rgas / aChemComp->mTemp;
    }
}

void 
CGasIdeal::pComputeConcFromMolarity(vector<double> &c, vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, CChemicalComposition* aChemComp)
{
    // First loop over primary species
    for(int i=0; i!=aSV1.size(); i++)        
    {
        c[mSpeciesIndices[aSV1[i]]] = Rgas * aChemComp->mTemp * c[mSpeciesIndices[aSV1[i]]];
    }

    // Second loop over secondary species
    for(int i=0; i!=aSV2.size(); i++)        
    {
        c[mSpeciesIndices[aSV2[i]]] = Rgas * aChemComp->mTemp * c[mSpeciesIndices[aSV2[i]]];
    }
}
