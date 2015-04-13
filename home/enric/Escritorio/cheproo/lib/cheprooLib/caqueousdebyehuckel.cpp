#include "caqueousdebyehuckel.h"


CAqueousDebyeHuckel::CAqueousDebyeHuckel()
{
    this->mClassName = "CAqueousDebyeHuckel";
    this->A = 0.0;
    this->B = 0.0;
    this->Bdot = 0.0;
}

CAqueousDebyeHuckel::CAqueousDebyeHuckel(QDomElement aNode)
: CAqueous(aNode)
{
    this->mClassName = "CAqueousDebyeHuckel";
    this->A = 0.0;
    this->B = 0.0;
    this->Bdot = 0.0;
}

CAqueousDebyeHuckel::~CAqueousDebyeHuckel()
{
}

void
CAqueousDebyeHuckel::pComputeActivityCoeff(vector<double> &aConcVector, vector<double> &aActCoeffsVector, vector<CSpecies*> &aSV1,
                                           vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, double &mIonicStrength)
{
    this->ComputeIonicStrength(aConcVector, mIonicStrength, aSV1, aSV2, mSpeciesIndices);

    // First loop over primary species
    for(int i=0; i!=aSV1.size(); i++)        
    {
       double aIonSize = aSV1[i]->mIonicRadius;
       double aCharge = aSV1[i]->mIonicCharge;
       
       double aDenomin = 1.0 + aIonSize*B*sqrt(mIonicStrength);
       double aNumerator = -A*aCharge*aCharge*sqrt(mIonicStrength);
       double aLogGamma = aNumerator / aDenomin + Bdot*mIonicStrength;

       aActCoeffsVector.at(mSpeciesIndices[aSV1[i]])=pow(10.0,aLogGamma);

    }

    // Second loop over secondary species
    for(int i=0; i!=aSV2.size(); i++)        
    {
       double aIonSize = aSV2[i]->mIonicRadius;
       double aCharge = aSV2[i]->mIonicCharge;

       double aDenomin = 1.0 + aIonSize*B*sqrt(mIonicStrength);
       double aNumerator = -A*aCharge*aCharge*sqrt(mIonicStrength);
       double aLogGamma = aNumerator / aDenomin + Bdot*mIonicStrength;

       aActCoeffsVector.at(mSpeciesIndices[aSV2[i]])=pow(10.0,aLogGamma);

    }

};

void
CAqueousDebyeHuckel::pComputeActivityCoeffDerivsWrtSV(vector<double> &aSVVector, MatrixXd &MatrixToContrib, MatrixXd &aS1Mod,
                                                      vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices,
                                                      bool isContribToJacobian, bool isForRISA)
{
    double mIonicStrength = 0.0;

    vector<double> aIonicStrengthDerivs;
    aIonicStrengthDerivs.resize(mSpeciesIndices.size());

    vector<double> aActCoeffsVector;
    aActCoeffsVector.resize(mSpeciesIndices.size());

    // Compute activity coefficients
    this->pComputeActivityCoeff(aSVVector, aActCoeffsVector, aSV1, aSV2, mSpeciesIndices, mIonicStrength);

    // Compute derivative of ionic strength wrt concentrations
    this->ComputeIonicStrengthDerivsWrtConc(aIonicStrengthDerivs, aSV1, aSV2, mSpeciesIndices);

    // Matrix where I put dlnGamma/dlnc1 and dlnGamma/dlnc2 for all AQUEOUS species
    vector<vector<double> > actCoeffsDerivAux;

    // If I'm evaluating the contribution to jacobian and residual or RISA (i.e., I'm evaluating dlnGamma/dlnc1 and dlnGamma/dlnc2)
    if(aS1Mod.size()!=0 && aSV2.size()!=0)
    {
        if (isForRISA)
        {
            // Resize the aux matrix (Ns x N1)
            actCoeffsDerivAux.resize(mSpeciesIndices.size()); 
            for(int i=0; i!=actCoeffsDerivAux.size(); i++)
            {
                actCoeffsDerivAux[i].resize(aSV1.size());
            }

            // First I evaluate dlnGamma1/dlnc1

            // First loop over primary species
            for(int i=0; i!=aSV1.size(); i++)
            {
                int SV1Index1 = mSpeciesIndices[aSV1[i]];
                double aIonSize = aSV1[i]->mIonicRadius;
                double aCharge = aSV1[i]->mIonicCharge;
                
                double dGammadI = 2*sqrt(mIonicStrength)*(1.0 + aIonSize*B*sqrt(mIonicStrength))*(1.0 + aIonSize*B*sqrt(mIonicStrength));
                dGammadI = (-(A*aCharge*aCharge/dGammadI) + Bdot) * log(10.0);

                // Second loop over primary species
                for(int j=0; j!=aSV1.size(); j++)
                {
                    int SV1Index2 = mSpeciesIndices[aSV1[j]];

                    actCoeffsDerivAux[SV1Index1][SV1Index2] = dGammadI * aIonicStrengthDerivs[SV1Index2] * aSVVector[SV1Index2];
                }
            }

            // Then I evaluate dlnGamma2/dlnc1

            // First loop over secondary species
            for(int i=0; i!=aSV2.size(); i++)
            {
                int SV2Index = mSpeciesIndices[aSV2[i]];
                double aIonSize = aSV2[i]->mIonicRadius;
                double aCharge = aSV2[i]->mIonicCharge;

                double dGammadI = 2*sqrt(mIonicStrength)*(1.0 + aIonSize*B*sqrt(mIonicStrength))*(1.0 + aIonSize*B*sqrt(mIonicStrength));
                dGammadI = (-(A*aCharge*aCharge/dGammadI) + Bdot) * log(10.0);

                // Second loop over primary species
                for(int j=0; j!=aSV1.size(); j++)
                {
                    int SV1Index = mSpeciesIndices[aSV1[j]];

                    actCoeffsDerivAux[SV2Index][SV1Index] = dGammadI * aIonicStrengthDerivs[SV1Index] * aSVVector[SV1Index];
                }
            } 

            // Add the contribution to dfd2/dnlc1
            for(int i=0; i!=MatrixToContrib.rows(); i++)
            {
                int secGlobIndex = mSpeciesIndices[aSV2[i]] - aS1Mod.cols();
                for(int j=0; j!=aSV1.size(); j++)
                {
                    double sum1 = 0.0;
                    int primGlobIndex = mSpeciesIndices[aSV1[j]];
                    for(int k=0; k!=aSV1.size(); k++)
                    {
                        sum1 += aS1Mod(i,k) * actCoeffsDerivAux[k][primGlobIndex];
                    }
					sum1 += aS1Mod(i,j);

                    double sum2 = 0.0;
                    for(int k=0; k!=aSV2.size(); k++)
                    {
                        sum2 += aS1Mod(i,k) * actCoeffsDerivAux[mSpeciesIndices[aSV2[k]]][primGlobIndex];
                    }

                    MatrixToContrib(i,primGlobIndex) = MatrixToContrib(i,primGlobIndex) + sum1 + sum2;
                }
                
            }

        }
        else
        {
            if (isContribToJacobian) // I'm evaluating contribution to Jacobian
            {
                // Resize the aux matrix (Ns x N2)
                actCoeffsDerivAux.resize(mSpeciesIndices.size()); 
                for(int i=0; i!=actCoeffsDerivAux.size(); i++)
                {
                    actCoeffsDerivAux[i].resize(aS1Mod.rows());
                }

                // First I evaluate dlnGamma1/dlnc2

                // First loop over primary species
                for(int i=0; i!=aSV1.size(); i++)
                {
                    int SV1Index = mSpeciesIndices[aSV1[i]];
                    double aIonSize = aSV1[i]->mIonicRadius;
                    double aCharge = aSV1[i]->mIonicCharge;
                
                    double dGammadI = 2*sqrt(mIonicStrength)*(1.0 + aIonSize*B*sqrt(mIonicStrength))*(1.0 + aIonSize*B*sqrt(mIonicStrength));
                    dGammadI = (-(A*aCharge*aCharge/dGammadI) + Bdot) * log(10.0);

                    // Second loop over secondary species
                    for(int j=0; j!=aSV2.size(); j++)
                    {
                        int SV2Index = mSpeciesIndices[aSV2[j]] - aS1Mod.cols();

                        actCoeffsDerivAux[SV1Index][SV2Index] = dGammadI * aIonicStrengthDerivs[SV2Index] * aSVVector[SV2Index];
                    }
                }

                // Then I evaluate dlnGamma2/dlnc2

                // First loop over secondary species
                for(int i=0; i!=aSV2.size(); i++)
                {
                    int SV2Index1 = mSpeciesIndices[aSV2[i]];
                    double aIonSize = aSV2[i]->mIonicRadius;
                    double aCharge = aSV2[i]->mIonicCharge;
                
                    double dGammadI = 2*sqrt(mIonicStrength)*(1.0 + aIonSize*B*sqrt(mIonicStrength))*(1.0 + aIonSize*B*sqrt(mIonicStrength));
                    dGammadI = (-(A*aCharge*aCharge/dGammadI) + Bdot) * log(10.0);

                    // Second loop over secondary species
                    for(int j=0; j!=aSV2.size(); j++)
                    {
                        int SV2Index2 = mSpeciesIndices[aSV2[j]] - aS1Mod.cols();

                        actCoeffsDerivAux[SV2Index1][SV2Index2] = dGammadI * aIonicStrengthDerivs[SV2Index2] * aSVVector[SV2Index2];
                    }
                }

                // Add the contribution to the Jacobian
                for(int i=0; i!=aSV2.size(); i++)
                {
                    int secGlobIndex1 = mSpeciesIndices[aSV2[i]] - aS1Mod.cols();
                    for(int j=0; j!=aSV2.size(); j++)
                    {
                        double sum = 0.0;
                        int secGlobIndex2 = mSpeciesIndices[aSV2[j]] - aS1Mod.cols();
                        for(int k=0; k!=aSV1.size(); k++)
                        {
                            sum = sum + aS1Mod(secGlobIndex1,k) * actCoeffsDerivAux[k][secGlobIndex2];
                        }

                        MatrixToContrib(secGlobIndex1,secGlobIndex2) = MatrixToContrib(secGlobIndex1,secGlobIndex2) + actCoeffsDerivAux[mSpeciesIndices[aSV2[i]]][secGlobIndex2] - sum;
                    }
                
                }
            }

            else // I'm evaluating the contribution to the residual
            {
                // Resize the aux matrix (Ns x N1)
                actCoeffsDerivAux.resize(mSpeciesIndices.size()); 
                for(int i=0; i!=actCoeffsDerivAux.size(); i++)
                {
                    actCoeffsDerivAux[i].resize(aS1Mod.cols());
                }

                // First I evaluate dGamma1/dc1

                // First loop over primary species
                for(int i=0; i!=aSV1.size(); i++)
                {
                    int SV1Index1 = mSpeciesIndices[aSV1[i]];
                    double aIonSize = aSV1[i]->mIonicRadius;
                    double aCharge = aSV1[i]->mIonicCharge;
                
                    double dGammadI = 2*sqrt(mIonicStrength)*(1.0 + aIonSize*B*sqrt(mIonicStrength))*(1.0 + aIonSize*B*sqrt(mIonicStrength));
                    dGammadI = (-(A*aCharge*aCharge/dGammadI) + Bdot) * log(10.0);

                    // Second loop over primary species
                    for(int j=0; j!=aSV1.size(); j++)
                    {
                        int SV1Index2 = mSpeciesIndices[aSV1[j]];

                        actCoeffsDerivAux[SV1Index1][SV1Index2] = dGammadI * aIonicStrengthDerivs[SV1Index2] * aSVVector[SV1Index2];
                    }
                }

                // Then I evaluate dlnGamma2/dlnc1

                // First loop over secondary species
                for(int i=0; i!=aSV2.size(); i++)
                {
                    int SV2Index = mSpeciesIndices[aSV2[i]];
                    double aIonSize = aSV2[i]->mIonicRadius;
                    double aCharge = aSV2[i]->mIonicCharge;

                    double dGammadI = 2*sqrt(mIonicStrength)*(1.0 + aIonSize*B*sqrt(mIonicStrength))*(1.0 + aIonSize*B*sqrt(mIonicStrength));
                    dGammadI = (-(A*aCharge*aCharge/dGammadI) + Bdot) * log(10.0);

                    // Second loop over primary species
                    for(int j=0; j!=aSV1.size(); j++)
                    {
                        int SV1Index = mSpeciesIndices[aSV1[j]];

                        actCoeffsDerivAux[SV2Index][SV1Index] = dGammadI * aIonicStrengthDerivs[SV1Index] * aSVVector[SV1Index];
                    }
                } 

                // Add the contribution to the residual
                for(int i=0; i!=aSV2.size(); i++)
                {
                    int secGlobIndex = mSpeciesIndices[aSV2[i]] - aS1Mod.cols();
                    for(int j=0; j!=aSV1.size(); j++)
                    {
                        double sum = 0.0;
                        int primGlobIndex = mSpeciesIndices[aSV1[j]];
                        for(int k=0; k!=aSV1.size(); k++)
                        {
                            sum = sum + aS1Mod(secGlobIndex,k) * actCoeffsDerivAux[k][primGlobIndex];
                        }

                        MatrixToContrib(secGlobIndex,primGlobIndex) = MatrixToContrib(secGlobIndex,primGlobIndex) + sum - actCoeffsDerivAux[mSpeciesIndices[aSV2[i]]][primGlobIndex];
                    }
                
                }
            }
        }
     }

    actCoeffsDerivAux.clear();
};

void 
CAqueousDebyeHuckel::pRead(const QDomElement aNode)
{
    CAqueous::Read(aNode);
}

void 
CAqueousDebyeHuckel::pSet(CChemicalComposition* aChemComp)
{
    this->pComputeCoeffAB(aChemComp);
}

void 
CAqueousDebyeHuckel::pComputeCoeffAB(CChemicalComposition* aChemComp)
{
    if(aChemComp->mTemp < 100)
    {
        this->A = ALT1+(ALT2+(ALT3+(ALT4+ALT5*aChemComp->mTemp)*aChemComp->mTemp)*aChemComp->mTemp)*aChemComp->mTemp;
 
        this->B = BLT1+(BLT2+(BLT3+(BLT4+BLT5*aChemComp->mTemp)*aChemComp->mTemp)*aChemComp->mTemp)*aChemComp->mTemp;
 
        this->Bdot = DLT1+(DLT2+(DLT3+(DLT4+DLT5*aChemComp->mTemp)*aChemComp->mTemp)*aChemComp->mTemp)*aChemComp->mTemp;
    }

    else
    {
        this->A = AHT1+(AHT2+(AHT3+(AHT4+AHT5*aChemComp->mTemp)*aChemComp->mTemp)*aChemComp->mTemp)*aChemComp->mTemp;
 
        this->B = BHT1+(BHT2+(BHT3+(BHT4+BHT5*aChemComp->mTemp)*aChemComp->mTemp)*aChemComp->mTemp)*aChemComp->mTemp;
 
        this->Bdot = DHT1+(DHT2+(DHT3+(DHT4+DHT5*aChemComp->mTemp)*aChemComp->mTemp)*aChemComp->mTemp)*aChemComp->mTemp;
    }

}

void
CAqueousDebyeHuckel::pComputeActivityCoeffDerivs(CChemicalComposition* aChemComp, vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, 
                                                 VectorXd aDerivatives, bool isWrtTemp)
{
    if(isWrtTemp) // Derivative wrt temperature
    {
        double coeffADeriv = 0;
        double coeffBDeriv = 0;
        double coeffBdotDeriv = 0;
        double ionStrength = 0;

        this->pComputeCoeffABDerivsWrtT(aChemComp, coeffADeriv, coeffBDeriv, coeffBdotDeriv);

        this->ComputeIonicStrength(aChemComp->mConcentration, ionStrength, aSV1, aSV2, mSpeciesIndices);

        // First loop over primary species
        for(int i=0; i!=aSV1.size(); i++)
        {
            int SV1Index = mSpeciesIndices[aSV1[i]];
            double aIonSize = aSV1[i]->mIonicRadius;
            double aCharge = aSV1[i]->mIonicCharge;

            double aDenom = 1.0 + aIonSize* this->B * sqrt(ionStrength);
            double aValue = (this->A * aCharge * aCharge * aIonSize * ionStrength * coeffBDeriv - aCharge * aCharge * sqrt(ionStrength) * coeffADeriv * aDenom) / aDenom;
            aValue = aValue + ionStrength * coeffBdotDeriv * log(10.0);

            aDerivatives(SV1Index) = aValue;
        }

        for(int i=0; i!=aSV2.size(); i++)
        {
            int SV2Index = mSpeciesIndices[aSV2[i]];
            double aIonSize = aSV2[i]->mIonicRadius;
            double aCharge = aSV2[i]->mIonicCharge;

            double aDenom = 1.0 + aIonSize* this->B * sqrt(ionStrength);
            double aValue = (this->A * aCharge * aCharge * aIonSize * ionStrength * coeffBDeriv - aCharge * aCharge * sqrt(ionStrength) * coeffADeriv * aDenom) / aDenom;
            aValue = aValue + ionStrength * coeffBdotDeriv * log(10.0);

            aDerivatives(SV2Index) = aValue;
        }

         


    }
    else    // Derivative wrt pressure
    {
        cout << "Aqueous activity coefficients don't depend on pressure!";
    }
}

void 
CAqueousDebyeHuckel::pComputeCoeffABDerivsWrtT(CChemicalComposition* aChemComp, double &aA_deriv, double &aB_deriv, double &aBdot_deriv)
{
    if(aChemComp->mTemp < 100)
    {
        aA_deriv = ALT2+(2*ALT3+(3*ALT4+4*ALT5*aChemComp->mTemp)*aChemComp->mTemp)*aChemComp->mTemp;

        aB_deriv = BLT2+(2*BLT3+(3*BLT4+4*BLT5*aChemComp->mTemp)*aChemComp->mTemp)*aChemComp->mTemp;

        aBdot_deriv = DLT2+(2*DLT3+(3*DLT4+4*DLT5*aChemComp->mTemp)*aChemComp->mTemp)*aChemComp->mTemp;
    }

    else
    {
        aA_deriv = AHT2+(2*AHT3+(3*AHT4+4*AHT5*aChemComp->mTemp)*aChemComp->mTemp)*aChemComp->mTemp;

        aB_deriv = BHT2+(2*BHT3+(3*BHT4+4*BHT5*aChemComp->mTemp)*aChemComp->mTemp)*aChemComp->mTemp; 

        aBdot_deriv = DHT2+(2*DHT3+(3*DHT4+4*DHT5*aChemComp->mTemp)*aChemComp->mTemp)*aChemComp->mTemp;
    }
}