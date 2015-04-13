#include "creactionratelaw.h"


CReactionRateLaw::CReactionRateLaw()
{
    this->mClassName = "CReactionRateLaw";
}

CReactionRateLaw::CReactionRateLaw(QDomElement aNode)
: CCheprooBase(aNode)
{
    this->mClassName = "CReactionRateLaw";
}

CReactionRateLaw::~CReactionRateLaw()
{
}

void
CReactionRateLaw::Read(QDomElement aNode)
{
    // Assign the name of the LocalChemicalSystem
    QString aName;
    XmlQt::getXMLAttribute(aNode, XML_ATTR_NAME, aName); 
    this->SetName(aName);
    
    // Get the value of the activation energy
    XmlQt::getXMLAttribute(aNode, XML_ATTR_EA, this->mEa); // If it's not specified then is zero
    
}

void 
CReactionRateLaw::ComputeOmegaTerm(double &omega, double &eta, double &theta, double &omegaTerm, double &omegaTermDeriv)
{
    if(omega<1)
    {
        double parenthesis = 0;
        parenthesis = 1.0 - pow(omega,theta);
        if(parenthesis > 1e-5)
        {
            omegaTerm = - pow(parenthesis, eta);
            omegaTermDeriv = eta * theta * (parenthesis / omega) * (pow(parenthesis,eta-1));
        }
        else
        {
            omegaTerm = 0.0;
            omegaTermDeriv = 1.0;
        }
    }
    else
    {
        double parenthesis = 0;
        parenthesis = pow(omega,theta);
        omegaTerm = pow(parenthesis - 1, eta);
        omegaTermDeriv = eta * theta * (parenthesis / omega) * (omegaTerm / (parenthesis - 1));
    }
}

void 
CReactionRateLaw::ComputeArrheniusTerm(double &temp, double &ea, double &arrheniusTerm)
{
    arrheniusTerm = exp(-(1/(temp+273.15) - 1/298.15) / Rgas * ea);
}

