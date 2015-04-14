#include "caqueouspitzer.h"


CAqueousPitzer::CAqueousPitzer()
{
    this->mClassName = "CAqueousPitzer";
}

CAqueousPitzer::CAqueousPitzer(QDomElement aNode)
: CAqueous(aNode)
{
    this->mClassName = "CAqueousPitzer";
}

CAqueousPitzer::~CAqueousPitzer()
{
}

void
CAqueousPitzer::pComputeActivityCoeff(vector<double> &aConcVector, vector<double> &aActCoeffsVector, vector<CSpecies*> &aSV1,
                                           vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, double &mIonicStrength)
{
    // Evaluate ionic strength
    this->ComputeIonicStrength(aConcVector, mIonicStrength, aSV1, aSV2, mSpeciesIndices);

    // Evaluate Z
    double z = 0.0;
    this->ComputeZ(aConcVector, z, aSV1, aSV2, mSpeciesIndices);

    // Evaluate M
    double m = 0.0;
    this->ComputeM(aConcVector, m);



};

void
CAqueousPitzer::pComputeActivityCoeffDerivsWrtSV(vector<double> &aSVVector, MatrixXd &MatrixToContrib, MatrixXd &aS1Mod,
                                                      vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices,
                                                      bool isContribToJacobian, bool isForRISA)
{
    double mIonicStrength = 0.0;

    vector<double> aIonicStrengthDerivs;
    aIonicStrengthDerivs.resize(mSpeciesIndices.size());

    //vector<double> aActCoeffsVector;
    //aActCoeffsVector.resize(mSpeciesIndices.size());

    //// Compute activity coefficients
    //this->pComputeActivityCoeff(aSVVector, aActCoeffsVector, aSV1, aSV2, mSpeciesIndices, mIonicStrength);

    //// Compute derivative of ionic strength wrt concentrations
    //this->ComputeIonicStrengthDerivsWrtConc(aIonicStrengthDerivs, aSV1, aSV2, mSpeciesIndices);

    //// Matrix where I put dGamma/dc1 and dGamma/dc2 for all AQUEOUS species
    //vector<vector<double> > actCoeffsDerivAux;

    //// If I'm evaluating the contribution to jacobian and residual (i.e., I'm evaluating dGamma/dc1 and dGamma/dc2)
    //if(aS1Mod.size()!=0 && aSV2.size()!=0)
    //{
    //    if (isContribToJacobian) // I'm evaluating contribution to Jacobian
    //    {
    //        // Resize the aux matrix (Ns x N2)
    //        actCoeffsDerivAux.resize(mSpeciesIndices.size()); 
    //        for(int i=0; i!=actCoeffsDerivAux.size(); i++)
    //        {
    //            actCoeffsDerivAux[i].resize(aS1Mod.rows());
    //        }

    //        // First I evaluate dGamma1/dc2

    //        // First loop over primary species
    //        for(int i=0; i!=aSV1.size(); i++)
    //        {
    //            int SV1Index = mSpeciesIndices[aSV1[i]];
    //            double aIonSize = aSV1[i]->mIonicRadius;
    //            double aCharge = aSV1[i]->mIonicCharge;
    //            
    //            double dGammadI = 2*sqrt(mIonicStrength)*(1.0 + aIonSize*B*sqrt(mIonicStrength))*(1.0 + aIonSize*B*sqrt(mIonicStrength));
    //            dGammadI = (-(A*aCharge*aCharge/dGammadI) + Bdot) * log(10.0);

    //            // Second loop over secondary species
    //            for(int j=0; j!=aSV2.size(); j++)
    //            {
    //                int SV2Index = mSpeciesIndices[aSV2[j]] - aS1Mod.cols();

    //                actCoeffsDerivAux[SV1Index][SV2Index] = dGammadI * aIonicStrengthDerivs[SV2Index] * aSVVector[SV2Index];
    //            }
    //        }

    //        // Then I evaluate dGamma2/dc2

    //        // First loop over secondary species
    //        for(int i=0; i!=aSV2.size(); i++)
    //        {
    //            int SV2Index1 = mSpeciesIndices[aSV2[i]];
    //            double aIonSize = aSV2[i]->mIonicRadius;
    //            double aCharge = aSV2[i]->mIonicCharge;
    //            
    //            double dGammadI = 2*sqrt(mIonicStrength)*(1.0 + aIonSize*B*sqrt(mIonicStrength))*(1.0 + aIonSize*B*sqrt(mIonicStrength));
    //            dGammadI = (-(A*aCharge*aCharge/dGammadI) + Bdot) * log(10.0);

    //            // Second loop over secondary species
    //            for(int j=0; j!=aSV2.size(); j++)
    //            {
    //                int SV2Index2 = mSpeciesIndices[aSV2[j]] - aS1Mod.cols();

    //                actCoeffsDerivAux[SV2Index1][SV2Index2] = dGammadI * aIonicStrengthDerivs[SV2Index2] * aSVVector[SV2Index2];
    //            }
    //        }

    //        // Add the contribution to the Jacobian
    //        for(int i=0; i!=aSV2.size(); i++)
    //        {
    //            int secGlobIndex1 = mSpeciesIndices[aSV2[i]] - aS1Mod.cols();
    //            for(int j=0; j!=aSV2.size(); j++)
    //            {
    //                double sum = 0.0;
    //                int secGlobIndex2 = mSpeciesIndices[aSV2[j]] - aS1Mod.cols();
    //                for(int k=0; k!=aSV1.size(); k++)
    //                {
    //                    sum = sum + aS1Mod(secGlobIndex1,k) * actCoeffsDerivAux[k][secGlobIndex2];
    //                }

    //                MatrixToContrib(secGlobIndex1,secGlobIndex2) = MatrixToContrib(secGlobIndex1,secGlobIndex2) + actCoeffsDerivAux[mSpeciesIndices[aSV2[i]]][secGlobIndex2] - sum;
    //            }
    //            
    //        }
    //    }

    //    else // I'm evaluating the contribution to the residual
    //    {
    //        // Resize the aux matrix (Ns x N1)
    //        actCoeffsDerivAux.resize(mSpeciesIndices.size()); 
    //        for(int i=0; i!=actCoeffsDerivAux.size(); i++)
    //        {
    //            actCoeffsDerivAux[i].resize(aS1Mod.cols());
    //        }

    //        // First I evaluate dGamma1/dc1

    //        // First loop over primary species
    //        for(int i=0; i!=aSV1.size(); i++)
    //        {
    //            int SV1Index1 = mSpeciesIndices[aSV1[i]];
    //            double aIonSize = aSV1[i]->mIonicRadius;
    //            double aCharge = aSV1[i]->mIonicCharge;
    //            
    //            double dGammadI = 2*sqrt(mIonicStrength)*(1.0 + aIonSize*B*sqrt(mIonicStrength))*(1.0 + aIonSize*B*sqrt(mIonicStrength));
    //            dGammadI = (-(A*aCharge*aCharge/dGammadI) + Bdot) * log(10.0);

    //            // Second loop over primary species
    //            for(int j=0; j!=aSV1.size(); j++)
    //            {
    //                int SV1Index2 = mSpeciesIndices[aSV1[j]];

    //                actCoeffsDerivAux[SV1Index1][SV1Index2] = dGammadI * aIonicStrengthDerivs[SV1Index2] * aSVVector[SV1Index2];
    //            }
    //        }

    //        // Then I evaluate dGamma2/dc1

    //        // First loop over secondary species
    //        for(int i=0; i!=aSV2.size(); i++)
    //        {
    //            int SV2Index = mSpeciesIndices[aSV2[i]];
    //            double aIonSize = aSV2[i]->mIonicRadius;
    //            double aCharge = aSV2[i]->mIonicCharge;

    //            double dGammadI = 2*sqrt(mIonicStrength)*(1.0 + aIonSize*B*sqrt(mIonicStrength))*(1.0 + aIonSize*B*sqrt(mIonicStrength));
    //            dGammadI = (-(A*aCharge*aCharge/dGammadI) + Bdot) * log(10.0);

    //            // Second loop over primary species
    //            for(int j=0; j!=aSV1.size(); j++)
    //            {
    //                int SV1Index = mSpeciesIndices[aSV1[j]];

    //                actCoeffsDerivAux[SV2Index][SV1Index] = dGammadI * aIonicStrengthDerivs[SV1Index] * aSVVector[SV1Index];
    //            }
    //        } 

    //        // Add the contribution to the residual
    //        for(int i=0; i!=aSV2.size(); i++)
    //        {
    //            int secGlobIndex = mSpeciesIndices[aSV2[i]] - aS1Mod.cols();
    //            for(int j=0; j!=aSV1.size(); j++)
    //            {
    //                double sum = 0.0;
    //                int primGlobIndex = mSpeciesIndices[aSV1[j]];
    //                for(int k=0; k!=aSV1.size(); k++)
    //                {
    //                    sum = sum + aS1Mod(secGlobIndex,k) * actCoeffsDerivAux[k][primGlobIndex];
    //                }

    //                MatrixToContrib(secGlobIndex,primGlobIndex) = MatrixToContrib(secGlobIndex,primGlobIndex) + sum - actCoeffsDerivAux[mSpeciesIndices[aSV2[i]]][primGlobIndex];
    //            }
    //            
    //        }
    //    }
    // }

    //actCoeffsDerivAux.clear();
};

void 
CAqueousPitzer::pRead(const QDomElement aNode)
{
    CAqueous::Read(aNode);


}


void
CAqueousPitzer::ComputeZ(vector<double> &aConcVector, double &aZ,vector<CSpecies*> &aSV1,
                          vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices)
{
    double sum = 0.0;

    // Loop over primary species
    for(int i=0; i!=aSV1.size(); i++)
    {
        // Get the ionic charge of the species
        double aCharge = aSV1[i]->mIonicCharge;

        // Get global index of the species
        int species1Index = mSpeciesIndices[aSV1[i]];

        sum = sum + aConcVector[species1Index] * abs(aCharge);
    }

    // Loop over secondary species
    for(int i=0; i!=aSV2.size(); i++)
    {
        // Get the ionic charge of the species
        double aCharge = aSV2[i]->mIonicCharge;

        // Get global index of the species
        int species2Index = mSpeciesIndices[aSV2[i]];

        sum = sum + aConcVector[species2Index] * abs(aCharge);
    }

    // Evaluate Z
    aZ = sum;

}

void
CAqueousPitzer::ComputeM(vector<double> &aConcVector, double &aM)
{
    double sum = 0.0;
    for(int i=0; i!= aConcVector.size(); i++)
    {
        sum = sum +  aConcVector[i];
    }
    aM = sum;

}