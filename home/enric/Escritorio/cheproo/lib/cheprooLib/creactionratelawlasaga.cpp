#include "creactionratelawlasaga.h"


CReactionRateLawLasaga::CReactionRateLawLasaga()
{
    this->mClassName = "CReactionRateLawLasaga";
    this->mArea = 0;
    this->omegaThreshold = 0;
    this->mCatalystsNames.clear();
    this->m_k.clear();
    this->m_n.clear();
    this->m_p.clear();
    this->m_q.clear();
}

CReactionRateLawLasaga::CReactionRateLawLasaga(QDomElement aNode)
: CReactionRateLaw(aNode)
{
    this->mClassName = "CReactionRateLawLasaga";
    this->mArea = 0;
    this->omegaThreshold = 0;
    this->mCatalystsNames.clear();
    this->m_k.clear();
    this->m_n.clear();
    this->m_p.clear();
    this->m_q.clear();

}

CReactionRateLawLasaga::~CReactionRateLawLasaga()
{
}

void
CReactionRateLawLasaga::Read(QDomElement aNode)
{
    CReactionRateLaw::Read(aNode);

    // Get tags with experimental terms (sum)
    vector<QDomElement> experimentalTerms = XmlQt::getElementsByTagName(aNode, XML_TAG_TERM);

    double k, n, p, q = 0;

    for(int i=0; i!=experimentalTerms.size(); i++)
    {
        XmlQt::getXMLAttribute(experimentalTerms.at(i), "k", k);
        XmlQt::getXMLAttribute(experimentalTerms.at(i), "k", p);
        XmlQt::getXMLAttribute(experimentalTerms.at(i), "k", q);
        
        this->m_k.push_back(k);
        this->m_p.push_back(p);
        this->m_q.push_back(q);


        // Get references to catalytic terms
        vector<QString> listOfIDs;
        vector<QDomElement> speciesTags;
        this->ReadReferences(experimentalTerms.at(i), 
                             XML_TAG_SPECIES, 
                             speciesTags, 
                             &listOfIDs);


        if(speciesTags.size()!=0)
        {
            for(int j=0; j!=speciesTags.size(); j++)
            {

                // For each j-th catalytic term put its name and n(j) in the relative attribute
                this->mCatalystsNames.push_back(speciesTags[j].text().toStdString());

                XmlQt::getXMLAttribute(experimentalTerms.at(i), "n", n);
                if(n=0) n=1.0; // default value
                this->m_n.push_back(n);
                n=0;

            }
        }
    }

    // Get tags with parameters (area)
    QDomElement params = XmlQt::getFirstElemByTag(aNode, XML_TAG_PARAM);
    XmlQt::getXMLAttribute(params, "area", mArea);
     
}

void 
CReactionRateLawLasaga::ComputeReactionRate(CChemicalComposition* aChemComp, map<CSpecies*, int> mSpeciesIndices, double &omega, double &reactionRate, double area)
{
    double omegaTerm, omegaTermDeriv, p, q, expTerm = 0.0;

    reactionRate = 0;

    VectorXd catalyTerm;

    if(area==0)
    {
        this->ComputeArrheniusTerm(aChemComp->mTemp, this->mEa, expTerm);
    }
    else
    {
        this->ComputeArrheniusTerm(aChemComp->mTemp, area, expTerm);
    }

    // Loop over experimental terms
    for(int i = 0; i != m_k.size(); i++)
    {
        p = this->m_p[i];
        q = this->m_q[i];

        // If there are catalytic terms
        if(this->mCatalystsNames.size() != 0)
        {
            catalyTerm.resize(mCatalystsNames.size());
            map<CSpecies*,int>::iterator it;

            // Get concentration and activity coefficients of the catalytic species
            for(int j=0; j!=mCatalystsNames.size(); j++)
            {
                // Get species indices
                for(it=mSpeciesIndices.begin(); it!=mSpeciesIndices.end(); it++)
                {
                    if(it->first->name().toStdString() == mCatalystsNames[j])
                    {
                        // Save in the vector their activities
                        catalyTerm(j) = aChemComp->mConcentration.at(it->second) * aChemComp->mActivityCoeff.at(it->second);
                    }
                }
            }
        }
        else
        {
            catalyTerm.resize(1);
            catalyTerm(0) = 1.0;
        }

        this->ComputeOmegaTerm(omega, p, q, omegaTerm, omegaTermDeriv);

        for(int j=0; j!= catalyTerm.size(); j++)
        {
            reactionRate *= pow(catalyTerm(j), m_n[i+j]);
        }
        reactionRate += this->m_k[i] * omegaTerm;
    }

    if(area!=0) reactionRate *= area * expTerm;
    else reactionRate *= expTerm;
}

void
CReactionRateLawLasaga::ComputeReactionRateDerivs(CChemicalComposition* aChemComp, map<CSpecies*, int> mSpeciesIndices, vector<CSpecies*> &aSV1, 
                                                  vector<CSpecies*> &aSV2, MatrixXd dc2_dc1, double &omega, VectorXd &rateDeriv, double area)
{
}