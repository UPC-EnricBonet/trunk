#include "creaction.h"
#include "cglobalchemicalsystem.h"

CReaction::CReaction()
{
    this->mClassName = "CReaction";
    this->mCoeffsLogK.resize(0);
    this->mTempVsLogK.clear();
    this->mIsKinetic = false;

}

CReaction::CReaction(QDomElement aNode)
: CCheprooBase(aNode)
{
    this->mClassName = "CReaction";
    this->mCoeffsLogK.resize(0);
    this->mTempVsLogK.clear();
    this->mIsKinetic = false;
}


CReaction::~CReaction()
{
}

void 
CReaction::SetDatabase(QString aDatabaseName)
{
    this->mDatabaseName = aDatabaseName;
}


void 
CReaction::Read(const QDomElement aNode, bool isKinetic)
{
    QDomElement reactionsElem = XmlQt::getFirstElemByTag(aNode,XML_TAG_REACTIONS);
    if (!reactionsElem.isNull())
    {
        // Retrieve the reaction element in the xml 
        vector<QDomElement> reactionElem =  XmlQt::GetElemsByNameAndAttrib(reactionsElem, XML_TAG_REACTION, XML_ATTR_NAME, this->mName); 

        // Read referenced species for the reaction, where the IDs are the stoichiometric coefficients of the reaction
        vector<QString> listOfIDs;
        vector<QDomElement> speciesTags;
        this->ReadReferences(reactionElem.at(0),
                             XML_TAG_SPECIES, 
                             speciesTags, 
                             &listOfIDs);

        vector<string> aRelSpeciesNames;
        int definedSpeciesNum = 0;
		QString name;
        for(unsigned i=0; i!= speciesTags.size(); i++)
        {
            //check if all the species referenced have been previously defined and created
			XmlQt::getXMLAttribute(speciesTags.at(i), XML_ATTR_VALUE, name);
			//XmlQt::getXMLAttribute(speciesTags.at(i), XML_ATTR_NAME, name);
            aRelSpeciesNames.push_back(name.trimmed().toStdString());
            
            if ( this->mGlobChemSysPointer->mAllSpecies[aRelSpeciesNames[i]] != 0 )
            {
                definedSpeciesNum += 1;
            }
        }

        // If all the species referenced have been defined in the GlobalChemicalSystem, then the reaction is possible.
        // Read the reaction properties
        if ( definedSpeciesNum == speciesTags.size() )
        {
            int nrRelspecies = speciesTags.size();
            vector<double> aStoichCoeffs;
            aStoichCoeffs.resize(nrRelspecies);

            for(int i = 0 ; i != nrRelspecies; i++ )
            {
                double aStoichCoeff = listOfIDs.at(i).toDouble();

                // Get the species from the list of all the species
                CSpecies* aRelSpecies = this->mGlobChemSysPointer->mAllSpecies[aRelSpeciesNames[i]];
             
                // Insert elements in the mRelSpeciesVsStoichCoeffs map
                this->mRelSpeciesVsStoichCoeffs.insert(pair<CSpecies*, double>(aRelSpecies,aStoichCoeff));

            }
                
            // Assigns the equilibrium constant to the current reaction 
            QDomElement logKElements = XmlQt::getFirstElemByTag(reactionElem.at(0),XML_TAG_LOGK);
                        
            if(!logKElements.isNull())
            {
                vector<QDomElement> tempTags;
                listOfIDs.clear();
                this->ReadReferences(logKElements,
                                     XML_TAG_TEMP, 
                                     tempTags, 
                                     &listOfIDs);
                
                double logKValue = 0.0;
                double tempValue = 0.0;

                vector<double> tempValues;
                VectorXd logKValues = VectorXd::Zero(tempTags.size());

                for(int itemp = 0; itemp!=tempTags.size(); itemp++)
                {
                    XmlQt::getXMLAttribute(tempTags[itemp],XML_ATTR_VALUE,logKValue);
                    tempValues.push_back(listOfIDs.at(itemp).toDouble());
                    logKValues(itemp) = logKValue;
                    
                    // This will have to be deleted when I delete mTempVsLogK
                    this->mTempVsLogK.insert(pair<double,double>(listOfIDs.at(itemp).toDouble(),logKValue));
                }

                if(tempValues.size()>1)
                {
                    this->mCoeffsLogK = VectorXd::Zero(tempValues.size());
                    MatrixXd coeffs = MatrixXd::Zero(5,tempValues.size());
                    for(int i=0; i!= tempValues.size(); i++)
                    {
                        coeffs(0,i)=log(tempValues.at(i) + 273.15);
                        coeffs(1,i)=1.0;
                        coeffs(2,i)=tempValues.at(i) + 273.15;
                        coeffs(3,i)=1/(tempValues.at(i) + 273.15);
                        coeffs(4,i)=1/((tempValues.at(i) + 273.15)*(tempValues.at(i) + 273.15));
                    }
                    //this->mCoeffsLogK = coeffs.jacobiSvd(ComputeThinU | ComputeThinV).solve(logKValues); // CHECK THIS WHEN I HAVE TIME!
                }
                else
                {
                    this->mLogK = logKValue;
                }
			    
            }
            
        }

    }


}


void
CReaction::SolveMALExplicit(double &aConcValue, int &aIndex, vector<double> &aConcVector, vector<double> &aActCoeffsVector, double aTemp, map<CSpecies*, int> mSpeciesIndices)
{

   map<CSpecies*, double>::iterator it;
   double product = 1.0;
   int index = 0;
   CSpecies* mySpecies; 
   for(it = this->mRelSpeciesVsStoichCoeffs.begin(); it!=this->mRelSpeciesVsStoichCoeffs.end(); it++)
   {
       // Get index of the related species
       index = mSpeciesIndices[it->first];

       if(index == aIndex ) 
       {
           mySpecies = it->first;
           continue;
       }
       if(it->first->mIsConstActivity) continue;

       product *= pow(aConcVector[index]*aActCoeffsVector[index],it->second);
   }

   // Update the concentration value
   double k = pow(10.0, this->mTempVsLogK[aTemp]);
   aConcValue = k / (product * pow(aActCoeffsVector[aIndex],mSpeciesIndices[mySpecies]));

   aConcValue = pow(aConcValue,1/mSpeciesIndices[mySpecies]);

}

void 
CReaction::SetLogK(double &temperature)
{
    if(mCoeffsLogK.size()!=0)
    {
        double temp = temperature + 273.15;

        this->mLogK = this->mCoeffsLogK[0] * log10(temp) + this->mCoeffsLogK[1] + this->mCoeffsLogK[2] * temp + this->mCoeffsLogK[3] / temp + this->mCoeffsLogK[4] / (temp*temp);
    }
}