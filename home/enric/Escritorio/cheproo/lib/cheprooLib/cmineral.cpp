#include "cmineral.h"

CMineral::CMineral()
{
    this->mClassName = "CMineralPhase";
}

CMineral::CMineral(QDomElement aNode)
: CPhase(aNode)
{
    this->mClassName = "CMineralPhase";
}

CMineral::~CMineral()
{
}

void
CMineral::ComputeActivityCoeff(vector<double> &aConcVector, vector<double> &aActCoeffsVector, vector<CSpecies*> &aSV1,
                              vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, double mIonicStrength)
{
    this->pComputeActivityCoeff(aConcVector, aActCoeffsVector, aSV1, aSV2, mSpeciesIndices, mIonicStrength);
};

void
CMineral::ComputeActivityCoeffDerivs(vector<double> &aSVVector, MatrixXd &MatrixToContrib, MatrixXd &aS1Mod,
                                     vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices,
                                     bool isContribToJacobian, bool isForRISA)
{
    this->pComputeActivityCoeffDerivsWrtSV(aSVVector, MatrixToContrib, aS1Mod, aSV1, aSV2, mSpeciesIndices, isContribToJacobian, isForRISA);
};

void 
CMineral::Read(const QDomElement aNode)
{
    CPhase::Read(aNode);
}

void 
CMineral::Set(CChemicalComposition* aChemComp)
{
    this->pSet(aChemComp);
}

void 
CMineral::ComputeMolarity(vector<double> &c, vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, CChemicalComposition* aChemComp)
{
    this->pComputeMolarity(c, aSV1, aSV2, mSpeciesIndices, aChemComp);
}

void 
CMineral::ComputeConcFromMolarity(vector<double> &c, vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, CChemicalComposition* aChemComp)
{
    this->pComputeConcFromMolarity(c, aSV1, aSV2, mSpeciesIndices, aChemComp);
}

