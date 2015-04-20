#include "cgas.h"

CGas::CGas()
{
    this->mClassName = "CGasPhase";
}

CGas::CGas(QDomElement aNode)
: CPhase(aNode)
{
    this->mClassName = "CGasPhase";
}
CGas::~CGas()
{
}

void 
CGas::Read(const QDomElement aNode)
{
    CPhase::Read(aNode);
}

void
CGas::ComputeActivityCoeff(vector<double> &aConcVector, vector<double> &aActCoeffsVector, vector<CSpecies*> &aSV1,
                              vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, double mIonicStrength)
{
    this->pComputeFugacityCoeff(aConcVector, aActCoeffsVector, aSV1, aSV2, mSpeciesIndices, mIonicStrength);
};

void
CGas::ComputeActivityCoeffDerivs(vector<double> &aSVVector, MatrixXd &MatrixToContrib, MatrixXd &aS1Mod,
                                     vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices,
                                     bool isContribToJacobian, bool isForRISA)
{
    this->pComputeFugacityCoeffDerivsWrtSV(aSVVector, MatrixToContrib, aS1Mod, aSV1, aSV2, mSpeciesIndices, isContribToJacobian, isForRISA);
};

void 
CGas::Set(CChemicalComposition* aChemComp)
{
    this->pSet(aChemComp);
}

void 
CGas::ComputeMolarity(vector<double> &c, vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, CChemicalComposition* aChemComp)
{
    this->pComputeMolarity(c, aSV1, aSV2, mSpeciesIndices, aChemComp);
}

void 
CGas::ComputeConcFromMolarity(vector<double> &c, vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, CChemicalComposition* aChemComp)
{
    this->pComputeConcFromMolarity(c, aSV1, aSV2, mSpeciesIndices, aChemComp);
}
