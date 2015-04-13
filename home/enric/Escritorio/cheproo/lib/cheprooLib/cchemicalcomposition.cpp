#include "cchemicalcomposition.h"
#include "clocalchemicalsystem.h"
#include "cglobalchemicalsystem.h"
#include "cspecies.h"

CChemicalComposition::CChemicalComposition()
{
    this->mClassName="CChemicalComposition";
    this->mConcentration.clear(); 
    this->mActivityCoeff.clear();
    this->mNodesVolume = 0.0;
    this->mTemp = 0.0;
    this->mLiquidDensity = 1000;

}

CChemicalComposition::CChemicalComposition(QDomElement aNode)
: CCheprooBase(aNode)
{
    this->mClassName="CChemicalComposition";
    this->mConcentration.clear();
    this->mActivityCoeff.clear();
    this->mNodesVolume = 0.0;
    this->mTemp = 0.0;
    this->mLiquidDensity = 1000;
}

CChemicalComposition::CChemicalComposition(CChemicalComposition* aChemComp, double& p, double& q, double& r, double& s, double& t, bool isPerturbedValue)
: CCheprooBase()
{
    this->mClassName="CChemicalComposition";

    this->mName = aChemComp->name();

    this->mLocChemSysPointer = aChemComp->mLocChemSysPointer;

    this->mConcentration.resize(aChemComp->mConcentration.size());

    this->mActivityCoeff.resize(aChemComp->mActivityCoeff.size());

    map<string, double>::iterator it;
    for(it=aChemComp->mVolumetricFractions.begin(); it!=aChemComp->mVolumetricFractions.end(); it++)
    {
        this->mVolumetricFractions.insert(pair<string, double>(it->first, it->second));
    }

    this->mTemp = aChemComp->mTemp;

    this->mLiquidDensity = aChemComp->mLiquidDensity;

   if(aChemComp->iCon.size()!=0)
    {
        this->iCon.resize(aChemComp->iCon.size());
        this->constrValue.resize(aChemComp->constrValue.size());
        this->constraints.resize(aChemComp->constraints.size());
        this->refSpecies.resize(aChemComp->refSpecies.size());

        int spGlobIndex = 0;

        for(int i=0; i!=this->iCon.size(); i++)
        {
            this->iCon.at(i) = aChemComp->iCon.at(i);
            this->constrValue(i) = aChemComp->constrValue(i);
            this->constraints.at(i) = aChemComp->constraints.at(i);
            this->refSpecies.at(i) = aChemComp->refSpecies.at(i);

            string aSpeciesName;

            if(iCon.at(i)=="eqgas") // Assign concentration of gases and minerals in equilibrium
            {
                // Retrieve the species name
		        aSpeciesName = constraints.at(i);
        
                // Retrieve its index in the local chemical system
                CSpecies* currSpecies = mLocChemSysPointer->mGlobChemSysPointer->mAllSpecies[aSpeciesName];
                spGlobIndex = mLocChemSysPointer->mSpeciesIndices[currSpecies];

                // Put the concentration value in the vector
                mConcentration.at(spGlobIndex) = constrValue(i);
                mActivityCoeff.at(spGlobIndex) = 1.0;
            }
            else if(iCon.at(i)=="eqmin")
            {
                // Retrieve the species name
		        aSpeciesName = constraints.at(i);
        
                // Retrieve its index in the local chemical system
                CSpecies* currSpecies = mLocChemSysPointer->mGlobChemSysPointer->mAllSpecies[aSpeciesName];
                spGlobIndex = mLocChemSysPointer->mSpeciesIndices[currSpecies];

                // Put the concentration value in the vector
                mConcentration.at(spGlobIndex) = 1.0; 
                mActivityCoeff.at(spGlobIndex) = 1.0; 
            }

        }

        // Assign c and gamma for h2o = 1 (by now here)
        // Retrieve its index in the local chemical system
        CSpecies* currSpecies = mLocChemSysPointer->mGlobChemSysPointer->mAllSpecies["h2o"];
        if(currSpecies != 0)
        {
            spGlobIndex = mLocChemSysPointer->mSpeciesIndices[currSpecies];
            mConcentration.at(spGlobIndex) = 1.0; 
            mActivityCoeff.at(spGlobIndex) = 1.0; 
        }
    }

    if(aChemComp->iCon_uncertain.size()!=0)
    {
        if(isPerturbedValue) // If I am modifying the constrValue_uncertain with the perturbed value
        {
            this->errors.resize(aChemComp->errors.size());
            this->iCon_uncertain.resize(aChemComp->errors.size());
            this->constrValue_uncertain.resize(aChemComp->errors.size());
            this->constraints_uncertain.resize(aChemComp->errors.size());
            this->refSpecies_uncertain.resize(aChemComp->errors.size());
            this->isLogData.resize(aChemComp->errors.size());
            this->isRelError.resize(aChemComp->errors.size());

            for(int i=0; i!=this->errors.size(); i++)
            {
                this->errors.at(i) = aChemComp->errors.at(i);
                this->iCon_uncertain.at(i) = aChemComp->iCon_uncertain.at(i);
                this->constraints_uncertain.at(i) = aChemComp->constraints_uncertain.at(i);
                this->refSpecies_uncertain.at(i) = aChemComp->refSpecies_uncertain.at(i);
                this->isLogData.at(i) = aChemComp->isLogData.at(i);
                this->isRelError.at(i) = aChemComp->isRelError.at(i);

                if(iCon_uncertain.at(i)=="alk") this->measured_value = aChemComp->constrValue_uncertain(i);  
            }
            if(this->name()=="water1")
			{
				this->constrValue_uncertain(this->errors.size()-2) = p;
                this->constrValue_uncertain(this->errors.size()-1) = q;
			}
			else if (this->name()=="water2")
			{
				this->constrValue_uncertain(this->errors.size()-3) = p;
				this->constrValue_uncertain(this->errors.size()-2) = q;
                this->constrValue_uncertain(this->errors.size()-1) = r;
			}
			else if (this->name()=="water3")
			{
				this->constrValue_uncertain(this->errors.size()-4) = p;
				this->constrValue_uncertain(this->errors.size()-3) = q;
				this->constrValue_uncertain(this->errors.size()-2) = r;
                this->constrValue_uncertain(this->errors.size()-1) = s;
			}
			else if (this->name()=="water4")
			{
				this->constrValue_uncertain(this->errors.size()-5) = p;
				this->constrValue_uncertain(this->errors.size()-4) = q;
				this->constrValue_uncertain(this->errors.size()-3) = r;
				this->constrValue_uncertain(this->errors.size()-2) = s;
                this->constrValue_uncertain(this->errors.size()-1) = t;

			}
        }
        else // If I am modifying the error
        {
            this->errors.resize(aChemComp->errors.size());
            this->iCon_uncertain.resize(aChemComp->errors.size());
            this->constrValue_uncertain.resize(aChemComp->errors.size());
            this->constraints_uncertain.resize(aChemComp->errors.size());
            this->refSpecies_uncertain.resize(aChemComp->errors.size());
            this->isLogData.resize(aChemComp->errors.size());
            this->isRelError.resize(aChemComp->errors.size());


            for(int i=0; i!=this->errors.size(); i++)
            {
                this->errors.at(i) = aChemComp->errors.at(i);
                this->iCon_uncertain.at(i) = aChemComp->iCon_uncertain.at(i);
                this->constrValue_uncertain(i) = aChemComp->constrValue_uncertain(i);
                this->constraints_uncertain.at(i) = aChemComp->constraints_uncertain.at(i);
                this->refSpecies_uncertain.at(i) = aChemComp->refSpecies_uncertain.at(i);
                this->isLogData.at(i) = aChemComp->isLogData.at(i);
                this->isRelError.at(i) = aChemComp->isRelError.at(i);
            }

            // Modify the last error value
            this->errors.back() = p;
        }
    }

    // Fill concentration vector with first guess
    string aSpeciesName;
    int spGlobIndex;
    map<string,double>::iterator it1;
    for(it1=aChemComp->cGuess.begin(); it1!=aChemComp->cGuess.end(); it1++)
    {
        // Retrieve the primary species name
		aSpeciesName = it1->first;

        // Retrieve its index in the local chemical system
        CSpecies* currSpecies = mLocChemSysPointer->mGlobChemSysPointer->mAllSpecies[aSpeciesName];
        spGlobIndex = mLocChemSysPointer->mSpeciesIndices[currSpecies];

        // Put the concentration value in the vector
        mConcentration.at(spGlobIndex) = aChemComp->cGuess[aSpeciesName];
    }
}
void 
CChemicalComposition::Read(const QDomElement aNode, vector<string> &primSpecies, set<string> &cas)
{
    // Read the temperature 
    XmlQt::getXMLAttribute(aNode, XML_ATTR_TEMP, this->mTemp); 

    // Read referenced first guess
    vector<QString> listOfIDs;
    vector<QDomElement> firstGuessTags;
    this->ReadReferences(aNode,
                         XML_TAG_FIRSTGUESS, 
                         firstGuessTags, 
                         &listOfIDs);

    // Get cGuess for relative primary species
    double guess = 0;
    string primSpName;
    for(unsigned i=0; i!= firstGuessTags.size(); i++)
    {
        XmlQt::getXMLAttribute(firstGuessTags.at(i), XML_ATTR_CGUESS, guess);
        XmlQt::getXMLAttribute(firstGuessTags.at(i), XML_ATTR_PRIMARYSPECIES, primSpName);

        this->cGuess.insert(pair<string, double>(primSpName,guess));
        primSpecies.push_back(primSpName);
    }
    listOfIDs.clear();
    firstGuessTags.clear();

    // Read referenced species for the water
    vector<QDomElement> constraintTags;
    this->ReadReferences(aNode,
                         XML_TAG_CONSTRAINT, 
                         constraintTags, 
                         &listOfIDs);

    this->constrValue.resize(constraintTags.size());
    this->iCon.resize(constraintTags.size());
    this->constraints.resize(constraintTags.size());
	this->refSpecies.resize(constraintTags.size());

    this->constrValue_uncertain.resize(constraintTags.size());
    this->iCon_uncertain.resize(constraintTags.size());
    this->constraints_uncertain.resize(constraintTags.size());
	this->refSpecies_uncertain.resize(constraintTags.size());
    this->errors.resize(constraintTags.size());
    this->isLogData.resize(constraintTags.size());
    this->isRelError.resize(constraintTags.size());

    double error = 0;
    for(unsigned i=0; i!= constraintTags.size(); i++)
    {
        // Get the initial conditions
        XmlQt::getXMLAttribute(constraintTags.at(i), XML_ATTR_STD, error);
                
        // If an error on the condition is given, then put data for RISA in the corresponding vectors
        if(error!=0) 
        {
            string isLogValue;
            string isRel;
            this->errors.at(i) = error;
            XmlQt::getXMLAttribute(constraintTags.at(i), XML_ATTR_ICON, iCon_uncertain.at(i));

            XmlQt::getXMLAttribute(constraintTags.at(i), XML_ATTR_ISLOGVALUE, isLogValue);
            XmlQt::getXMLAttribute(constraintTags.at(i), XML_ATTR_ISRELERROR, isRel);
            if(isLogValue!="")
            {
                if(isLogValue=="true")
                { 
                    isLogData.at(i)=true;
                    if(iCon_uncertain.at(i)=="chgbal") throw GeneralException(ERROR_LOGSCALE_NOT_ALLOWED); // Check if logValue is given for cghbal constraint
                }
                else
                {
                    isLogData.at(i)=false;
                }
            }
            else
            {
                if(iCon_uncertain.at(i)=="chgbal" || iCon_uncertain.at(i)=="cTot" || iCon_uncertain.at(i)=="ec" || iCon_uncertain.at(i)=="alk")
                {
                    isLogData.at(i)=false; // aritmetic default for these conditions
                }
                else
                {
                    isLogData.at(i)=true; // logaritmic default for these conditions
                }
            }
            if(isRel!="")
            {
                if(isRel=="true")
                { 
                    isRelError.at(i)=true;
                }
                else
                {
                    isRelError.at(i)=false;
                }
            }
            else
            {
                if(isLogData.at(i)==true)
                {
                    isRelError.at(i)=false; // 
                }
                else
                {
                    isRelError.at(i)=true; // 
                }
            }


            XmlQt::getXMLAttribute(constraintTags.at(i), XML_ATTR_VALUE, constrValue_uncertain(i));
            XmlQt::getXMLAttribute(constraintTags.at(i), XML_ATTR_NAME, constraints_uncertain.at(i));
			XmlQt::getXMLAttribute(constraintTags.at(i), XML_ATTR_REFSPECIES, refSpecies_uncertain.at(i));

            //***OJO: delete later after RISA simulations!
            if(iCon_uncertain.at(i)=="alk") this->measured_value = constrValue_uncertain(i);  
        }
        else
        {
            XmlQt::getXMLAttribute(constraintTags.at(i), XML_ATTR_ICON, iCon.at(i));
            XmlQt::getXMLAttribute(constraintTags.at(i), XML_ATTR_VALUE, constrValue(i));
            XmlQt::getXMLAttribute(constraintTags.at(i), XML_ATTR_NAME, constraints.at(i));
			XmlQt::getXMLAttribute(constraintTags.at(i), XML_ATTR_REFSPECIES, refSpecies.at(i));

            // If system is in equilibrium with a gas (pgas is fixed) then this gas is CAS as well
            if(iCon.at(i) == "eqgas" || iCon.at(i) == "activity") cas.insert(constraints.at(i));

        }
        error = 0;      
    }

    listOfIDs.clear();
    constraintTags.clear();

    // Read referenced phase contents
    vector<QDomElement> phaseContentTags;
    this->ReadReferences(aNode,
                         XML_TAG_PHASECONTENT, 
                         phaseContentTags, 
                         &listOfIDs);

    string phaseName;
    double thetaValue;
    string model;
    for(unsigned i=0; i!= phaseContentTags.size(); i++)
    {
        XmlQt::getXMLAttribute(phaseContentTags.at(i), XML_ATTR_NAME, phaseName);
        XmlQt::getXMLAttribute(phaseContentTags.at(i), XML_ATTR_MODEL, model);
        XmlQt::getXMLAttribute(phaseContentTags.at(i), XML_ATTR_VALUE, thetaValue);

        // Fill the mVolumetricFractions attribute and cas vector
        this->mVolumetricFractions.insert(pair<string,double>(phaseName,thetaValue));
        if(model=="mineral" && thetaValue>0) cas.insert(phaseName);
    }

}

CChemicalComposition::~CChemicalComposition(void)
{
}

void 
CChemicalComposition::GetSetOfConcVector(VectorXd &concVector, vector<int> &speciesIndices)
{
    for(int i=0; i!=speciesIndices.size(); i++)
    {
        concVector(i) = this->mConcentration.at(speciesIndices[i]); 
    }

}

vector<double> 
CChemicalComposition::GetVector(eCheprooVariableType aVar, bool useClassicComponents)
{
	vector<double> values;

	if(aVar == eSecondaryConc)
	{
		// Get secondary species indices
		vector<int> &spcIndices = mLocChemSysPointer->GetSecondarySpcIndices();

		// Get secondary concentration values
		values.resize(spcIndices.size());
		for(int i=0; i!=spcIndices.size(); i++)
		{
			values[i] = this->mConcentration.at(spcIndices[i]); 
		}
	}

	else if(aVar == eAqComponent)
	{
        if(useClassicComponents) values = mLocChemSysPointer->Get_u_aq_classic(this);
		else values = mLocChemSysPointer->Get_u_aq(this);
	}
    
    else if(aVar == eGasComponent)
	{
		values = mLocChemSysPointer->Get_u_gas(this);

	}

    else if(aVar == eMinComponent)
	{
		values = mLocChemSysPointer->Get_u_min(this);

	}
    else if(aVar == eSurfComponent)
	{
		//values = mLocChemSysPointer->Get_u_surf(this->mConcentration);

	}
    else if(aVar == eKineticSinkSource)
    {
        values = mLocChemSysPointer->GetUSkrk(this);
    }
	else if(aVar == ePrimaryConc)
	{
		// Get secondary species indices
		vector<int> &spcIndices = mLocChemSysPointer->GetPrimarySpcIndices();

		// Get secondary concentration values
		values.resize(spcIndices.size());
		for(int i=0; i!=spcIndices.size(); i++)
		{
			values[i] = this->mConcentration.at(spcIndices[i]); 
		}
	}

	return values;
}

double 
CChemicalComposition::GetScalar(eCheprooVariableType aVar)
{
    double result = 0;

    if(aVar == eLiqDensity)
    {
        result = mLiquidDensity;
    }
    else if(aVar == eGasDensity)
    {
        result = mGasDensity;
    }
    if(aVar == eLiqViscosity)
    {
        result = mLiquidViscosity;
    }
    else if(aVar == eGasViscosity)
    {
        result = mGasViscosity;
    }

    return result;
}

void 
CChemicalComposition::GetSetOfLnConcVector(VectorXd &concVector, vector<int> &speciesIndices)
{
    for(int i=0; i!=speciesIndices.size(); i++)
    {
        if(concVector(i) != 0)
        {
            concVector(i) = log(this->mConcentration.at(speciesIndices[i])); 
        }
    }

}


void 
CChemicalComposition::GetSetOfActivityCoeffVector(VectorXd &actCoeffVector, vector<int> &speciesIndices)
{
    for(int i=0; i!=speciesIndices.size(); i++)
    {
        actCoeffVector(i) = this->mActivityCoeff.at(speciesIndices[i]); 
    }

}

void 
CChemicalComposition::UpdateConcentrations(VectorXd &concVector, vector<int> &speciesIndices)
{
    for(int i=0; i!=speciesIndices.size(); i++)
    {
       this->mConcentration.at(speciesIndices[i]) = concVector(i);
    }


}

void 
CChemicalComposition::UpdateConcentration(double &concValue, int &speciesIndices)
{
    this->mConcentration.at(speciesIndices) = concValue;
}

void 
CChemicalComposition::UpdateConcentrationsFromLn(VectorXd &concVector, vector<int> &speciesIndices)
{
    for(int i=0; i!=speciesIndices.size(); i++)
    {
       this->mConcentration.at(speciesIndices[i]) = exp(concVector[i]);
    }


}

void 
CChemicalComposition::UpdateActivityCoeffs(vector<double> &actCoeffsVector, vector<int> &speciesIndices)
{
    for(int i=0; i!=speciesIndices.size(); i++)
    {
       this->mActivityCoeff.at(speciesIndices[i]) = actCoeffsVector[i];
    }

}

void 
CChemicalComposition::Set_c1(vector<double> &c1)
{
	// First set the value of c1:
	for(int i=0; i!=c1.size(); i++)
	{
		this->mConcentration.at(i) = c1[i]; 
	}

	// Then comute c2:
	mLocChemSysPointer->ComputeSecondaries(this, c1);
}

MatrixXd 
CChemicalComposition::Get_dc_dc1_byPhase(QString aVarType)
{
    MatrixXd result = this->mLocChemSysPointer->Compute_dc_dc1_ByPhase(aVarType, this);

    return result;
}

MatrixXd 
CChemicalComposition::GetUByPhase(QString aVarType)
{
    MatrixXd result = this->mLocChemSysPointer->GetComponentMatrixByPhase(aVarType);

    return result;
}

void 
CChemicalComposition::UpdateVolumFraction(double value, string phaseName)
{
    this->mVolumetricFractions[phaseName] = value;
}

void
CChemicalComposition::Copy(CChemicalComposition* aChemComp)
{
    this->mConcentration.resize(aChemComp->mConcentration.size());
    this->mConcentration = aChemComp->mConcentration;

    this->mActivityCoeff.resize(aChemComp->mActivityCoeff.size());
    this->mActivityCoeff = aChemComp->mActivityCoeff;

    map<string, double>::iterator it;
    for(it=aChemComp->mVolumetricFractions.begin(); it!=aChemComp->mVolumetricFractions.end(); it++)
    {
        this->mVolumetricFractions.insert(pair<string, double>(it->first, it->second));
    }

    this->mTemp = aChemComp->mTemp;
    this->mLiquidDensity = aChemComp->mLiquidDensity;

}

void
CChemicalComposition::CopyForRT(CChemicalComposition* aChemComp)
{
    this->mConcentration.resize(aChemComp->mConcentration.size());
    this->mConcentration = aChemComp->mConcentration;

    this->mActivityCoeff.resize(aChemComp->mActivityCoeff.size());
    this->mActivityCoeff = aChemComp->mActivityCoeff;

    map<string, double>::iterator it;
    for(it=aChemComp->mVolumetricFractions.begin(); it!=aChemComp->mVolumetricFractions.end(); it++)
    {
        this->mVolumetricFractions.insert(pair<string, double>(it->first, it->second));
    }

    this->mTemp = aChemComp->mTemp;
    this->mLiquidDensity = aChemComp->mLiquidDensity;
    this->mGasDensity = aChemComp->mGasDensity;
    this->mGasViscosity = aChemComp->mGasViscosity;
    this->mLiquidViscosity = aChemComp->mLiquidViscosity;
    this->mLiquidPressure = aChemComp->mLiquidPressure;

    // Think about those attributes, especially constrValue, if they need to be public or if they can be deleted at some point
    
    this->iCon.resize(aChemComp->iCon.size());
    this->iCon = aChemComp->iCon;
    this->constrValue.resize(aChemComp->constrValue.size());
    this->constrValue = aChemComp->constrValue;
    this->constraints.resize(aChemComp->constraints.size());
    this->constraints = aChemComp->constraints;
    this->refSpecies.resize(aChemComp->refSpecies.size());
    this->refSpecies = aChemComp->refSpecies;


}

void
CChemicalComposition::SetChemComp()
{
    set<string>::iterator it;

    vector<string> iConTemp;
    VectorXd constrValueTemp;
    vector<string> constraintsTemp;
	vector<string> refSpeciesTemp;

    // Resize useful temporary vectors
    int aSpeciesNum = this->mLocChemSysPointer->mSpeciesIndices.size();
    this->mConcentration.resize(aSpeciesNum);
    this->mActivityCoeff.resize(aSpeciesNum);

    iConTemp.resize(aSpeciesNum);
    constrValueTemp = VectorXd::Zero(aSpeciesNum);
    constraintsTemp.resize(aSpeciesNum);
	refSpeciesTemp.resize(aSpeciesNum);

    string aSpeciesName;
    int spGlobIndex;
    map<string,double>::iterator it1;

    // Fill concentration vector with first guess
    for(it1=cGuess.begin(); it1!=cGuess.end(); it1++)
    {
        // Retrieve the primary species name
		aSpeciesName = it1->first;

        // Retrieve its index in the local chemical system
        CSpecies* currSpecies = mLocChemSysPointer->mGlobChemSysPointer->mAllSpecies[aSpeciesName];
        spGlobIndex = mLocChemSysPointer->mSpeciesIndices[currSpecies];

        // Put the concentration value in the vector
        mConcentration.at(spGlobIndex) = this->cGuess[aSpeciesName];
    }

    // Fill constraints for the "certain" system
    for(unsigned i=0; i!=iCon.size(); i++)
    {
        if(iCon.at(i)=="") continue;

        //if(iCon.at(i)=="cTot" || iCon.at(i)=="activity")
        //{
            // Retrieve the species name
		    aSpeciesName = constraints.at(i);
        
            // Retrieve its index in the local chemical system
            CSpecies* currSpecies = mLocChemSysPointer->mGlobChemSysPointer->mAllSpecies[aSpeciesName];
            spGlobIndex = mLocChemSysPointer->mSpeciesIndices[currSpecies];

            // Fill vectors
            iConTemp.at(spGlobIndex) = iCon.at(i);
            constrValueTemp(spGlobIndex) = constrValue(i);
            constraintsTemp.at(spGlobIndex) = constraints.at(i);
		    refSpeciesTemp.at(spGlobIndex) = refSpecies.at(i);
     /*   }*/
        if(iCon.at(i)=="eqgas") // Assign concentration of gases and minerals in equilibrium
        {
            // Retrieve the species name
		    aSpeciesName = constraints.at(i);
        
            // Retrieve its index in the local chemical system
            CSpecies* currSpecies = mLocChemSysPointer->mGlobChemSysPointer->mAllSpecies[aSpeciesName];
            spGlobIndex = mLocChemSysPointer->mSpeciesIndices[currSpecies];

            // Put the concentration value in the vector
            mConcentration.at(spGlobIndex) = constrValue(i);
            mActivityCoeff.at(spGlobIndex) = 1.0;
        }
        else if(iCon.at(i)=="eqmin")
        {
            // Retrieve the species name
		    aSpeciesName = constraints.at(i);
        
            // Retrieve its index in the local chemical system
            CSpecies* currSpecies = mLocChemSysPointer->mGlobChemSysPointer->mAllSpecies[aSpeciesName];
            spGlobIndex = mLocChemSysPointer->mSpeciesIndices[currSpecies];

            // Put the concentration value in the vector
            mConcentration.at(spGlobIndex) = 1.0; 
            mActivityCoeff.at(spGlobIndex) = 1.0; 
        }
	}

    // Clear and Resize attributes
    iCon.clear();
    constrValue.resize(0);
    constraints.clear();
	refSpecies.clear();

    iCon.resize(aSpeciesNum);
    constrValue = VectorXd::Zero(aSpeciesNum);
    constraints.resize(aSpeciesNum);
	refSpecies.resize(aSpeciesNum);

    // Copy vectors into attributes
    for(int i=0; i!=aSpeciesNum; i++)
    {
        iCon.at(i) = iConTemp.at(i);
        constrValue(i) = constrValueTemp(i);
        constraints.at(i) = constraintsTemp.at(i);
		refSpecies.at(i) = refSpeciesTemp.at(i);

        // Assign c and gamma for h2o = 1 (by now here)
        // Retrieve its index in the local chemical system
        CSpecies* currSpecies = mLocChemSysPointer->mGlobChemSysPointer->mAllSpecies["h2o"];
        if(currSpecies != 0)
        {
            spGlobIndex = mLocChemSysPointer->mSpeciesIndices[currSpecies];
            mConcentration.at(spGlobIndex) = 1.0; 
            mActivityCoeff.at(spGlobIndex) = 1.0; 
        }
    }

    iConTemp.resize(0);
    constrValueTemp.resize(0);
    constraintsTemp.resize(0);
	refSpeciesTemp.clear();
}


void
CChemicalComposition::ChemicalStepSIA(vector<double> &u_aq, double &f_aq, double &f_s, double &f_tot, bool useClassicComponents)
{
	string phaseName ="";
    multimap<ePhaseKind,CPhase*>::iterator it;
    pair<multimap<ePhaseKind,CPhase*>::iterator,multimap<ePhaseKind,CPhase*>::iterator> it1;

    // Put theta where they have to go

    // - Access the aqueous phase
    it1 = this->mLocChemSysPointer->mPhases.equal_range(eAqueousPhase); 

    // - Get the name of the phase
    phaseName = it1.first->second->name().toStdString();

	this->mVolumetricFractions[phaseName] = f_aq / f_tot;
	phaseName ="";

	// - Get mineral phase name
    it1 = this->mLocChemSysPointer->mPhases.equal_range(eMineralPhase); 

    if(it1.first!=this->mLocChemSysPointer->mPhases.end()) 
    {
        for (it=it1.first; it!=it1.second; it++) // OSS: it won't work for more than one mineral phase...same for ads and gases
        {
            this->mVolumetricFractions[phaseName] = f_s / f_tot;
        }
    }
    phaseName ="";


	// - Get surface phase name
    it1 = this->mLocChemSysPointer->mPhases.equal_range(eSurfacePhase); 

    if(it1.first!=this->mLocChemSysPointer->mPhases.end()) 
    {
        for (it=it1.first; it!=it1.second; it++) // OSS: it won't work for more than one mineral phase...same for ads and gases
        {
            this->mVolumetricFractions[phaseName] = f_s / f_tot;
        }
    }
    phaseName ="";

	// Calculate u_tot_classic for SIA
    VectorXd u_tot_classic = this->mLocChemSysPointer->Calculate_utot_classic_SIA(this, u_aq);

	// Speciate from u_tot
	if(useClassicComponents)
	{
        // Put the new reduced primary species in the vector
        for(int i=0; i!=u_tot_classic.size(); i++)
        {
            this->constrValue(i) = u_tot_classic(i);
        }

        // Speciate
		this->mLocChemSysPointer->Speciate_utot_classic(this);
	}
	else
	{
		// First calculate u_tot_red from u_tot_classic
        VectorXd u_tot_red;
		u_tot_red = this->mLocChemSysPointer->Calculate_ured_from_utot_classic(this, u_tot_classic);

        // Put the new reduced primary species in the vector
        for(int i=0; i!=u_tot_red.size(); i++)
        {
            this->constrValue(i) = u_tot_red(i);
        }

        // Speciate
		this->mLocChemSysPointer->Speciate_utot(this);
	}


}	