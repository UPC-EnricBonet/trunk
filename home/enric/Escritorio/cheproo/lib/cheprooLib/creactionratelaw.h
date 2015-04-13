#ifndef CREACTIONRATELAW_H
#define CREACTIONRATELAW_H

#include "ccheproobase.h"
#include "cspecies.h"

#include <QString>
#include <string>
#include <vector>
#include <map>
#include <math.h>

#include <Eigen/Dense>
using namespace Eigen;

using namespace std;

class CReaction;
class CChemicalComposition;

#define Rgas 8.314472

class CReactionRateLaw : public CCheprooBase
{

public:

    double mEa; ///< Activation energy


public:

    virtual void ComputeReactionRate(CChemicalComposition* aChemComp,
                                     map<CSpecies*, int> mSpeciesIndices,
                                     double &omega,
                                     double &reactionRate,
                                     double area = 0){};

    virtual void ComputeReactionRateDerivs(CChemicalComposition* aChemComp,
                                           map<CSpecies*, int> mSpeciesIndices,
                                           vector<CSpecies*> &aSV1,
                                           vector<CSpecies*> &aSV2,
                                           MatrixXd dc2_dc1,
                                           double &omega,
                                           VectorXd &rateDeriv,
                                           double area = 0){};

    void ComputeOmegaTerm(double &omega,
                          double &eta,
                          double &theta,
                          double &omegaTerm,         ///< Omega term of the reaction law
                          double &omegaTermDeriv);   ///< Derivative of the omega term wrt omega

    void ComputeArrheniusTerm(double &temp,
                              double &ea,
                              double &arrheniusTerm);

    void Read(
        QDomElement aNode ///< XML node with attributes
        );

        // \brief Constructor from a XML file.
	CReactionRateLaw(
		QDomElement aNode ///< XML node with attributes
        );

	CReactionRateLaw();
	~CReactionRateLaw();
};

#endif
