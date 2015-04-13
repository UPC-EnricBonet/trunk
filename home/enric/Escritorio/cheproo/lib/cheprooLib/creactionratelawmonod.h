#ifndef CREACTIONRATELAWMONOD_H
#define CREACTIONRATELAWMONOD_H

#include "creactionratelaw.h"

class CReactionRateLawMonod : public CReactionRateLaw
{
public:


public:

    void ComputeReactionRate(CChemicalComposition* aChemComp,
                             map<CSpecies*, int> mSpeciesIndices,
                             double &omega,
                             double &reactionRate){};

    void ComputeReactionRateDerivs(CChemicalComposition* aChemComp,
                                           map<CSpecies*, int> mSpeciesIndices,
                                           vector<CSpecies*> &aSV1,
                                           vector<CSpecies*> &aSV2,
                                           MatrixXd dc2_dc1,
                                           double &omega,
                                           VectorXd &rateDeriv,
                                           double area = 0){};

    void Read(
        QDomElement aNode ///< XML node with attributes
        );

    // \brief Constructor from a XML file.
	CReactionRateLawMonod(
		QDomElement aNode ///< XML node with attributes
        );

	CReactionRateLawMonod();
	~CReactionRateLawMonod();
};

#endif
