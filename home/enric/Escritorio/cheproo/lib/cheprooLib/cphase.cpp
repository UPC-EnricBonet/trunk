#include "cphase.h"
#include "cglobalchemicalsystem.h"
#include "clocalchemicalsystem.h"

CPhase::CPhase()
{
    this->mClassName = "CPhase";
	this->mSpeciesNames.clear();
    mGlobChemSysPointer = 0;
}

CPhase::CPhase(QDomElement aNode)
: CCheprooBase(aNode)
{
    this->mClassName = "CPhase";
	this->mSpeciesNames.clear();
    mGlobChemSysPointer = 0;
}

CPhase::~CPhase()
{
}

ePhaseKind CPhase::str2PhasesKind(QString s)
{
	if ( s == "aqueous" ) return eAqueousPhase;
	if ( s == "mineral" ) return eMineralPhase;
	if ( s == "gas" ) return eGasPhase;
    if ( s == "surface" ) return eSurfacePhase;
	throw GeneralException(ERROR_PHASE_KIND_NAME);
}

QString CPhase::phasesKind2str(ePhaseKind e)
{
	if ( e == eAqueousPhase ) return "aqueous";
	if ( e == eMineralPhase) return "mineral";
	if ( e == eGasPhase ) return "gas";
    if ( e == eSurfacePhase ) return "surface";
	throw GeneralException(ERROR_PHASE_KIND_NUMBER);
}

void 
CPhase::Read(const QDomElement aNode)
{
	ObjectManager* om = ObjectManager::instance();

    QString aSpeciesModel;
    XmlQt::getXMLAttribute(aNode, XML_ATTR_MODEL, aSpeciesModel); 

    // Load all the species that are contained in this phase.
    vector<QString> listOfIDs;
    vector<QDomElement> speciesTags;
    this->ReadReferences(aNode,
                         XML_TAG_SPECIES, 
                         speciesTags, 
                         &listOfIDs);
	for (unsigned i = 0 ; i < speciesTags.size(); i++ )
	{
		// Retrieve the species name
		const QString aSpeciesName = speciesTags[i].text();

		// Add to the relative species names
		this->mSpeciesNames.push_back(aSpeciesName.toStdString());

        // Call CGlobalChemicalSystem to create the species and read their properties 
        // (this is something that is always done, regardless of the number of LocalChemicalSystem)
        this->mGlobChemSysPointer->CreateSpecies(aSpeciesName, aSpeciesModel);

    }

    speciesTags.clear();
    listOfIDs.clear();

}
void
CPhase::ComputeDensity(double &density, CChemicalComposition* aChemComp, bool mustCalcDerivs)
{
    this->pComputeDensity(density, aChemComp, false);
}

void
CPhase::ComputeViscosity()
{
    this->pComputeViscosity(false);
}

void 
CPhase::ComputeIonicStrength(vector<double> &aConcVector, double &mIonicStrength,  vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2,
                              map<CSpecies*, int> &mSpeciesIndices)
{
      this->pComputeIonicStrength(aConcVector, mIonicStrength, aSV1, aSV2, mSpeciesIndices);
}

void
CPhase::ComputeIonicStrengthDerivsWrtConc(vector<double> &aIonicStrengthDerivs, vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2,
                                            map<CSpecies*, int> &mSpeciesIndices)

{
    this->pComputeIonicStrengthDerivsWrtConc(aIonicStrengthDerivs,aSV1,aSV2,mSpeciesIndices);
}

void 
CPhase::ComputeChargeBalance(vector<double> &aConcVector, double &aChargeBalance, vector<CSpecies*> &aSV1, vector<CSpecies*> &aSV2,
                              map<CSpecies*, int> &mSpeciesIndices)
{
    this->pComputeChargeBalance(aConcVector, aChargeBalance, aSV1, aSV2, mSpeciesIndices);
}

void 
CPhase::ComputeChargeBalanceDerivs(MatrixXd &dc2_dc1,
                                    VectorXd &aChargeBalanceDerivs,
                                    vector<CSpecies*> &aSV1,
                                    vector<CSpecies*> &aSV2,
                                    map<CSpecies*, int> &mSpeciesIndices)
{
    this->pComputeChargeBalanceDerivsWrtConc(dc2_dc1, aChargeBalanceDerivs, aSV1, aSV2, mSpeciesIndices);
}