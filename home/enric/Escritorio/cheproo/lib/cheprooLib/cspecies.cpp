#include "cspecies.h"
#include "cglobalchemicalsystem.h"
#include "generalexception.h"


CSpecies::CSpecies(QString aSpeciesType)
{
    this->mClassName = "CSpecies";
	this->mIsConstActivity = false;

    if(aSpeciesType=="aqueous")
    {
        this->mSpeciesKind = eAqueous;
        this->mIonicCharge = 0.0;
        this->mIonicRadius = 0.0;
        this->mElectricCharge = 0.0;
        this->mLimMolConduct = 0.0;
        this->mMolWeight = 0.0;
    }
    else if(aSpeciesType=="mineral")
    {
        this->mSpeciesKind = eMineral;
        this->mMolVolume = 0.0;
    }
    else if(aSpeciesType=="gas")
    {
        this->mSpeciesKind = eGas;
        this->mMolVolume = 0.0;
    }
    else if(aSpeciesType=="surface")
    {
        this->mSpeciesKind = eSurface;
        this->mElectricCharge = 0.0;
        this->mCEC = 0.0;
        this->mMolWeight = 0.0;
    }
    else throw GeneralException(ERROR_SPECIES_KIND_NAME);

}

CSpecies::CSpecies(QDomElement aNode)
: CCheprooBase(aNode)
{
    this->mClassName = "CSpecies";
    this->mSpeciesKind = eLastSpeciesKind;
    this->mCEC = 0.0;
    this->mElectricCharge = 0.0;
    this->mLimMolConduct = 0.0;
    this->mIonicCharge = 0.0;
    this->mIonicRadius = 0.0;
    this->mMolVolume = 0.0;
    this->mMolWeight = 0.0;
}

CSpecies::~CSpecies()
{
}

eSpeciesKind CSpecies::str2SpeciesKind(QString s)
{
	if ( s == "aqueous" ) return eAqueous;
	if ( s == "mineral" ) return eMineral;
	if ( s == "gas" ) return eGas;
    if ( s == "surface" ) return eSurface;
	throw GeneralException(ERROR_SPECIES_KIND_NAME);
}

QString CSpecies::SpeciesKind2str(eSpeciesKind e)
{
	if ( e == eAqueous ) return "aqueous";
	if ( e == eMineral) return "mineral";
	if ( e == eGas ) return "gas";
    if ( e == eSurface ) return "surface";
	throw GeneralException(ERROR_SPECIES_KIND_NUMBER);
}

void 
CSpecies::SetDatabase(QString aDatabaseName)
{
    this->mDatabaseName = aDatabaseName;
}

void 
CSpecies::Read(const QDomElement aNode)
{
    // Retrieve the phases element
	QDomElement phasesElem = XmlQt::getFirstElemByTag(aNode,XML_TAG_PHASES);

    // Aqueous species
    if(this->mSpeciesKind == eAqueous)
    {
        // Retrieve the aqueous phase element in the xml 
        vector<QDomElement> aqueousPhaseElem =  XmlQt::GetElemsByNameAndAttrib(phasesElem, XML_TAG_PHASE, XML_ATTR_TYPE, "aqueous"); 

        if(!aNode.isNull())
        {

            // Retrieve the species element in the xml 
            vector<QDomElement> speciesRef =  XmlQt::GetElemsByNameAndAttrib(aqueousPhaseElem.at(0), XML_TAG_SPECIES, XML_ATTR_NAME, this->mName.trimmed()); 

            // Find the species tag with the specified name and read the properties
            if(speciesRef.size()!=1) 
            {
                cout << PrTr("CSpecies", ERROR_SPECIES_NOT_FOUND) <<  qPrintable(this->mName) << endl;            
            }

            // Get ionic charge
            XmlQt::getXMLAttribute(speciesRef.at(0), XML_ATTR_IONICCHARGE, this->mIonicCharge);

            // Get ionic radius
            XmlQt::getXMLAttribute(speciesRef.at(0), XML_ATTR_IONICRADIUS, this->mIonicRadius); 

            // Get molecular weight
            XmlQt::getXMLAttribute(speciesRef.at(0), XML_ATTR_MOLWEIGHT, this->mMolWeight); 

            // Get molecular weight
            XmlQt::getXMLAttribute(speciesRef.at(0), XML_ATTR_LIMMOLCONDUCTIVITY, this->mLimMolConduct); 

            // Default value for non tabulated species..it can be changed!
            if(this->mLimMolConduct==0.0) this->mLimMolConduct = 80.0;

            // If it's water, its activity is constant
            if(this->mName == "h2o") this->mIsConstActivity = true;

        }

    }
    else if (this->mSpeciesKind == eMineral)
    {
        // Retrieve the mineral phase element in the xml 
        vector<QDomElement> mineralPhaseElem =  XmlQt::GetElemsByNameAndAttrib(phasesElem, XML_TAG_PHASE, XML_ATTR_TYPE, "mineral"); 

        if(!aNode.isNull())
        {
                // Retrieve the species element in the xml 
                vector<QDomElement> speciesRef =  XmlQt::GetElemsByNameAndAttrib(mineralPhaseElem.at(0), XML_TAG_SPECIES, XML_ATTR_NAME, this->mName.trimmed()); 

				// Find the species tag with the specified name and read the properties
				if(speciesRef.size()!=1) 
				{
					cout << PrTr("CSpecies", ERROR_SPECIES_NOT_FOUND) <<  qPrintable(this->mName) << endl;            
				}

                // If it's a pure phase, mineral is constant activity species (think about non-pure minerals)
                this->mIsConstActivity = false;

                // Get molar volume
                XmlQt::getXMLAttribute(speciesRef.at(0), XML_ATTR_MOLARVOLUME, this->mMolVolume); 

        }
    }
    //else if (this->mSpeciesKind = eSurface)
    //{
    //    // Read surface species
    //}
    else if (this->mSpeciesKind = eGas)
    {
        // Retrieve the aqueous phase element in the xml 
        vector<QDomElement> gasPhaseElem =  XmlQt::GetElemsByNameAndAttrib(phasesElem, XML_TAG_PHASE, XML_ATTR_TYPE, "gas"); 

        if(!aNode.isNull())
        {
            // Loop over gas phases
            for(int i = 0; i!= gasPhaseElem.size(); i++)
            {
                // Retrieve the species element in the xml 
                vector<QDomElement> speciesRef =  XmlQt::GetElemsByNameAndAttrib(gasPhaseElem.at(0), XML_TAG_SPECIES, XML_ATTR_NAME, this->mName.trimmed()); 

                // Find the species tag with the specified name and read the properties
                if(speciesRef.size()!=1)
                {
                    // Then the phase is a mixture: read properties of the gases costituting the mixture
                }

                // Get molar volume
                XmlQt::getXMLAttribute(speciesRef.at(0), XML_ATTR_MOLARVOLUME, this->mMolVolume); 

            }
        }
    }
}

