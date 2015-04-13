#include "cgasbinarymixture.h"


CGasBinaryMixture::CGasBinaryMixture(void)
{
    this->mClassName = "CGasBinaryMixture";

}

CGasBinaryMixture::~CGasBinaryMixture(void)
{
    this->mClassName = "CGasBinaryMixture";

}

void
CGasBinaryMixture::pComputeDensity(double &density, CChemicalComposition* aChemComp, bool mustCalcDerivs)
{
}

void
CGasBinaryMixture::pComputeViscosity(bool mustCalcDerivs)
{
}

void 
CGasBinaryMixture::pRead(const QDomElement aNode)
{
    CGas::Read(aNode);
}
