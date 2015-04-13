#include "caqueous.h"

CAqueous::CAqueous()
{
    this->mClassName = "CAqueousPhase";
    this->mDensityKind = eDensityDefault;
}

CAqueous::CAqueous(QDomElement aNode)
: CPhase(aNode)
{
    this->mClassName = "CAqueousPhase";
    this->mDensityKind = eDensityDefault;
}

CAqueous::~CAqueous()
{
}

void
CAqueous::pComputeIonicStrength(vector<double> &aConcVector, double &mIonicStrength, vector<CSpecies*> &aSV1,
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

        sum = sum + aConcVector[species1Index] * aCharge * aCharge;
    }

    // Loop over secondary species
    for(int i=0; i!=aSV2.size(); i++)
    {
        // Get the ionic charge of the species
        double aCharge = aSV2[i]->mIonicCharge;

        // Get global index of the species
        int species2Index = mSpeciesIndices[aSV2[i]];

        sum = sum + aConcVector[species2Index] * aCharge * aCharge;
    }

    // Evaluate ionic strength
    mIonicStrength = 0.5 * sum;

};

void
CAqueous::pComputeIonicStrengthDerivsWrtConc(vector<double> &aIonicStrengthDerivs, vector<CSpecies*> &aSV1,
                                             vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices)
{
    // Loop over primary species
    for(int i=0; i!=aSV1.size(); i++)
    {
        // Get the ionic charge of the species
        double aCharge = aSV1[i]->mIonicCharge;

        // Get global index of the species
        int species1Index = mSpeciesIndices[aSV1[i]];

            aIonicStrengthDerivs[species1Index] = 0.5 * aCharge * aCharge;
    }

    // Loop over secondary species
    for(int i=0; i!=aSV2.size(); i++)
    {
        // Get the ionic charge of the species
        double aCharge = aSV2[i]->mIonicCharge;

        // Get global index of the species
        int species2Index = mSpeciesIndices[aSV2[i]];

            aIonicStrengthDerivs[species2Index] = 0.5 * aCharge * aCharge;
    }


};

void
CAqueous::pComputeChargeBalance(vector<double> &aConcVector, double &aChargeBalance, vector<CSpecies*> &aSV1,
                                vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices)
{
    double aCharge = 0.0;
    int aSpcIndex = 0;

    // Loop over primary species
    for(int i=0; i!=aSV1.size(); i++)
    {
        // Get charge of the species
        aCharge = aSV1[i]->mIonicCharge;

        // Get global index of the species
        aSpcIndex = mSpeciesIndices[aSV1[i]];

        // Evaluate charge balance
        aChargeBalance = aChargeBalance + aConcVector[aSpcIndex]*aCharge;
    }
        
    // Loop over secondary species
    for(int i=0; i!=aSV2.size(); i++)
    {
        // Get charge of the species
        aCharge = aSV2[i]->mIonicCharge;

        // Get global index of the species
        aSpcIndex = mSpeciesIndices[aSV2[i]];

        // Evaluate charge balance
        aChargeBalance = aChargeBalance + aConcVector[aSpcIndex]*aCharge;
    }
}

void 
CAqueous::pComputeChargeBalanceDerivsWrtConc(MatrixXd &dc2_dc1, VectorXd &aChargeBalanceDerivs, vector<CSpecies*> &aSV1,
                                        vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices)
{
    double aCharge = 0.0;
    int aSpc1Index = 0;
    int aSpc2Index = 0;

    // Loop over primary species
    for(int i=0; i!=aSV1.size(); i++)
    {
        // Get charge of the species
        aCharge = aSV1[i]->mIonicCharge;

        // Get global index of the species
        aSpc1Index = mSpeciesIndices[aSV1[i]];

        // Evaluate the charge balance derivative
        aChargeBalanceDerivs[aSpc1Index] = aCharge;
    }


    // Loop over primary species
    for(int i=0; i!=aSV1.size(); i++)
    {
        // Get global index of the species
        aSpc1Index = mSpeciesIndices[aSV1[i]];

        double sum = 0.0;

        // Loop over secondary species
        for(int j=0; j!=aSV2.size(); j++)
        {
            // Get charge of the species
            aCharge = aSV2[j]->mIonicCharge;

            // Get global index of the species
            aSpc2Index = mSpeciesIndices[aSV2[j]] - aSV1.size();
            
            // Evaluate charge balance
            sum = sum + aCharge * dc2_dc1(aSpc2Index,aSpc1Index);

        }
        aChargeBalanceDerivs[aSpc1Index] = aChargeBalanceDerivs[aSpc1Index] + sum;
    }
}

void
CAqueous::ComputeActivityCoeff(vector<double> &aConcVector, vector<double> &aActCoeffsVector, vector<CSpecies*> &aSV1,
                              vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, double mIonicStrength)
{
    this->pComputeActivityCoeff(aConcVector, aActCoeffsVector, aSV1, aSV2, mSpeciesIndices, mIonicStrength);
};

void
CAqueous::ComputeActivityCoeffDerivs(vector<double> &aSVVector, MatrixXd &MatrixToContrib, MatrixXd &aS1Mod,
                                     vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices,
                                     bool isContribToJacobian, bool isForRISA)
{
    this->pComputeActivityCoeffDerivsWrtSV(aSVVector, MatrixToContrib, aS1Mod, aSV1, aSV2, mSpeciesIndices, isContribToJacobian, isForRISA);
};

void 
CAqueous::ComputeActivityCoeffDerivs(CChemicalComposition* aChemComp, vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, 
                                VectorXd aDerivatives, bool isWrtTemp)
{
    this->pComputeActivityCoeffDerivs(aChemComp, aSV1, aSV2, mSpeciesIndices, aDerivatives, isWrtTemp);
}

void 
CAqueous::Read(const QDomElement aNode)
{
    //// Read parameters relative to density model
    //string name;
    //XmlQt::getXMLAttribute(aNode, XML_ATTR_MODEL_NAME, name);
    //if(name == "")this->mDensityKind = eDensityCB1;
    //else if(name == "")this->mDensityKind = eDensityCB2;

    //XmlQt::getXMLAttribute(aNode, XML_ATTR_MODEL, m_rho_ref); 
    //XmlQt::getXMLAttribute(aNode, XML_ATTR_MODEL, mCompressibility); 
    //XmlQt::getXMLAttribute(aNode, XML_ATTR_MODEL, mVolumThermCoeff); 
    //XmlQt::getXMLAttribute(aNode, XML_ATTR_MODEL, mSoluteVar); 
    //XmlQt::getXMLAttribute(aNode, XML_ATTR_MODEL, mRefPressure); 
     

    CPhase::Read(aNode);
}

void 
CAqueous::Set(CChemicalComposition* aChemComp)
{
    this->pSet(aChemComp);
}

void 
CAqueous::pComputeDensity(double &density, CChemicalComposition* aChemComp, bool mustCalcDerivs)
{
    double salinity = 0; // Salinity must be evaluated yet!

    if(mDensityKind = eDensityCB1)
    {  
        density = this->m_rho_ref * exp(this->mCompressibility * (aChemComp->mLiquidPressure - this->mRefPressure) + 
                  this->mVolumThermCoeff * aChemComp->mTemp + this->mSoluteVar * salinity);  
    }
    else if(mDensityKind = eDensityCB2)
    {
        density = this->m_rho_ref * (1 + this->mCompressibility * (aChemComp->mLiquidPressure - this->mRefPressure) + 
                  this->mVolumThermCoeff * aChemComp->mTemp + this->mSoluteVar * salinity);  
    }
    else if(mDensityKind = eDensityDefault) density = 1000;
}

void 
CAqueous::ComputeMolarity(vector<double> &c, vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, CChemicalComposition* aChemComp)
{
    double sum = 0;

    // Compute liquid water content
    double wH2OLiq = this->ComputeWH2OLiq(c, aSV1, aSV2, mSpeciesIndices);

    // Get density value
    double density = aChemComp->mLiquidDensity;
    if(density == 0)
    {
        //Option to evaluate it
        this->ComputeDensity(density, aChemComp, false);
    }

    // First loop over primary species
    for(int i=0; i!=aSV1.size(); i++)        
    {
        c[mSpeciesIndices[aSV1[i]]] = density * wH2OLiq * c[mSpeciesIndices[aSV1[i]]];
    }

    // Second loop over secondary species
    for(int i=0; i!=aSV2.size(); i++)        
    {
        c[mSpeciesIndices[aSV2[i]]] = density * wH2OLiq * c[mSpeciesIndices[aSV2[i]]];
    }


}

void 
CAqueous::ComputeConcFromMolarity(vector<double> &c, vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, CChemicalComposition* aChemComp)
{
    double sum = 0;

    // Compute liquid water content
    double wH2OLiq = 0;
    wH2OLiq = this->ComputeWH2OLiq(c, aSV1, aSV2, mSpeciesIndices);

    // Get density value
    double density = aChemComp->mLiquidDensity;
    if(density == 0)
    {
        //Option to evaluate it
        this->ComputeDensity(density, aChemComp, false);
    }

    // First loop over primary species
    for(int i=0; i!=aSV1.size(); i++)        
    {
        c[mSpeciesIndices[aSV1[i]]] = c[mSpeciesIndices[aSV1[i]]] / density / wH2OLiq;
    }

    // Second loop over secondary species
    for(int i=0; i!=aSV2.size(); i++)        
    {
        c[mSpeciesIndices[aSV2[i]]] = c[mSpeciesIndices[aSV2[i]]] / density / wH2OLiq;
    }


}

double 
CAqueous::ComputeWH2OLiq(vector<double> &c, vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices)
{
    double sum = 0;
    double w_H2O_liq = 0;

    // First loop over primary species
    for(int i=0; i!=aSV1.size(); i++)        
    {
        sum += aSV1[i]->mMolWeight * c[mSpeciesIndices[aSV1[i]]];
    }

    // Second loop over secondary species
    for(int i=0; i!=aSV2.size(); i++)        
    {
        sum += aSV2[i]->mMolWeight * c[mSpeciesIndices[aSV2[i]]];
    }

    w_H2O_liq = pow(1+sum, -1);

    return w_H2O_liq;
}
