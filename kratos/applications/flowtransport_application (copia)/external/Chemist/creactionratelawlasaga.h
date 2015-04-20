#ifndef CREACTIONRATELAWLASAGA_H
#define CREACTIONRATELAWLASAGA_H

#include "creactionratelaw.h"
#include "cchemicalcomposition.h"


class CReactionRateLawLasaga : public CReactionRateLaw
{
public:

    double mArea;                    ///< Specific area

    double omegaThreshold;           ///< Vector containing supersaturation threshold for the precipitation to start

    vector<double> m_k;              ///< Vector containing the k values (dimensions: number of exponential terms Nt)

    vector<double> m_p;              ///< Vector containing the p values (dimensions: number of exponential terms Nt)
            
    vector<double> m_q;              ///< Vector containing the q values (dimensions: number of exponential terms Nt)

    vector<string> mCatalystsNames;  ///< Vector containing the names of the catalysts (dimensions: Ns)

    vector<double> m_n;              ///< Vector containing n values (dimensions: Nt x Ns)


public:

    void ComputeReactionRate(CChemicalComposition* aChemComp,
                             map<CSpecies*, int> mSpeciesIndices,
                             double &omega,
                             double &reactionRate,
                             double area = 0);

    void ComputeReactionRateDerivs(CChemicalComposition* aChemComp,
                                           map<CSpecies*, int> mSpeciesIndices,
                                           vector<CSpecies*> &aSV1,
                                           vector<CSpecies*> &aSV2,
                                           MatrixXd dc2_dc1,
                                           double &omega,
                                           VectorXd &rateDeriv,
                                           double area = 0);



    void Read(
        QDomElement aNode ///< XML node with attributes
        );

    // \brief Constructor from a XML file.
	CReactionRateLawLasaga(
		QDomElement aNode ///< XML node with attributes
        );

	CReactionRateLawLasaga();
	~CReactionRateLawLasaga();
};

#endif
