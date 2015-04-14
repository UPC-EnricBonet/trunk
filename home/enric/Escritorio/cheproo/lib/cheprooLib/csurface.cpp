#include "csurface.h"

CSurface::CSurface()
{
    this->mClassName = "CSurfacePhase";
}

CSurface::CSurface(QDomElement aNode)
: CPhase(aNode)
{
    this->mClassName = "CSurfacePhase";
}

CSurface::~CSurface()
{
}

void
CSurface::ComputeActivityCoeff(vector<double> &aConcVector, vector<double> &aActCoeffsVector, vector<CSpecies*> &aSV1,
                              vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, double mIonicStrength)
{
    this->pComputeActivityCoeff(aConcVector, aActCoeffsVector, aSV1, aSV2, mSpeciesIndices, mIonicStrength);
};

void
CSurface::ComputeActivityCoeffDerivs(vector<double> &aSVVector, MatrixXd &MatrixToContrib, MatrixXd &aS1Mod,
                                     vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices,
                                     bool isContribToJacobian, bool isForRISA)
{
    this->pComputeActivityCoeffDerivsWrtSV(aSVVector, MatrixToContrib, aS1Mod, aSV1, aSV2, mSpeciesIndices, isContribToJacobian, isForRISA);
};

void 
CSurface::Read(const QDomElement aNode)
{
    CPhase::Read(aNode);
}

void 
CSurface::Set(CChemicalComposition* aChemComp)
{
    this->pSet(aChemComp);
}

void 
CSurface::ComputeMolarity(vector<double> &c, vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, CChemicalComposition* aChemComp)
{
    this->pComputeMolarity(c, aSV1, aSV2, mSpeciesIndices, aChemComp);
}

void 
CSurface::ComputeConcFromMolarity(vector<double> &c, vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2, map<CSpecies*, int> &mSpeciesIndices, CChemicalComposition* aChemComp)
{
    this->pComputeConcFromMolarity(c, aSV1, aSV2, mSpeciesIndices, aChemComp);
}

