// Unit Testing CAqueousDebyeHuckel

#include "cunittest.h"
#include "xmlconstants.h"
#include "xmlqt.h"
#include "objectlist.h"
#include "commonglobals.h"
#include "objectmanager.h"
#include "valarraytools.h"
#include "classnameconstants.h"
#include "ccheprooplusplus.h"
#include "cglobalchemicalsystem.h"

#include "caqueous.h"
#include "cphase.h"
#include "caqueousdebyehuckel.h"

#include <vector>
#include <QFile>
#include <Eigen/Dense>
#include <math.h>

using namespace Eigen;

static QString baseDir = "UnitTest/";

// Method to evaluate the activity coefficients
int UTCAqueousDebyeHuckel_1()
{
    QString inputFile = baseDir + "in/UTCheprooInput.xml";
    QString refFile = baseDir + "ref/UTCAqueousDebyeHuckel1.xml";

    QDomDocument doc("mydocument");
	XmlQt::xmlParseFile(inputFile, &doc);


    // Get a pointer to the object manager instance
    ObjectManager* om = ObjectManager::instance();

    // Fill all the lists with empty pointers. 
    bool mustDelete = true;
    bool typeAttributeMustExist = false;

    om->FillAllListsFromXML(doc.documentElement(), typeAttributeMustExist, mustDelete);

	// Retrieve the CheprooPlusPlus element in the xml and the name of the problem
    QDomElement nodeCheprooPlusPlus = XmlQt::getFirstElemByTag(doc.documentElement(),XML_TAG_CHEPROOPLUSPLUS);
    QString aProblemName;
    XmlQt::getXMLAttribute(nodeCheprooPlusPlus, XML_ATTR_NAME, aProblemName); 

    bool CheprooPlusPlusExists = ! nodeCheprooPlusPlus.isNull();

    // Declare vector to comprare the results
    valarray<valarray<double> > actCoeffs;

    if( CheprooPlusPlusExists ) 
    {
	    // Create a Problem object and read its properties from the xml
        CCheprooPlusPlus* cheprooPlusPlusPointer = dynamic_cast<CCheprooPlusPlus*>(om->newInstanceOfClassName(CLASS_NAME_CHEPROOPLUSPLUS, aProblemName));

        // Read and initialize the problem
        cheprooPlusPlusPointer->ReadAndInitialize(nodeCheprooPlusPlus, false);

        // Pointer to GlobalChemicalSystem
        CLocalChemicalSystem* myLocalPointer = cheprooPlusPlusPointer->mGlobalChemSysPointer->mLocChemSysVector[0];

        // Call the function to speciate initial waters
        myLocalPointer->SpeciateInitialWater(cheprooPlusPlusPointer->mGlobalChemSysPointer->mLocChemSysVector[0]->mRelChemicalCompositions[0]);

        // Declare variable to compare: Activity coefficients
        vector<double> activityCoeffs;
        activityCoeffs.resize(myLocalPointer->mPrimSpeciesNum + myLocalPointer->mEqReactions.size());

        // Obtain a pointer to the current phase
        CPhase* myAqPhase = cheprooPlusPlusPointer->mGlobalChemSysPointer->mAllPhases["aqueous1"];

        // Set the aqueous phase
        myAqPhase->Set(myLocalPointer->mRelChemicalCompositions.at(0));
        // Check the method
        myAqPhase->ComputeActivityCoeff(myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration, activityCoeffs,
                                        myLocalPointer->mSpecies[eAqueous1nc], myLocalPointer->mSpecies[eAqueous2],
                                        myLocalPointer->mSpeciesIndices);

        actCoeffs.resize(activityCoeffs.size());
        for(int i=0; i!=activityCoeffs.size(); i++)
        {
            actCoeffs[i].resize(1);
            actCoeffs[i][0] = activityCoeffs[i];
        }
        
    }

    //Obtains from sRefBal the content of reference file to compare for Balance
	QFile RefFileBal(refFile);
	RefFileBal.open(QIODevice::ReadOnly | QIODevice::Text);
	QString sRefBal = RefFileBal.readAll();
    RefFileBal.close();
    valarray< valarray< double > > ref;
    ValarrayTools::fromString(sRefBal, ref);
    
	if ( !ValarrayTools::IsEqual(actCoeffs,ref, 1e-12)) return 1;
    
	return 0;
}

// Method to evaluate the activity coefficients derivatives dGamma/dc1
int UTCAqueousDebyeHuckel_2()
{
    QString inputFile = baseDir + "in/UTCheprooInput.xml";
    QString refFile = baseDir + "ref/UTCAqueousDebyeHuckel2.xml";

    QDomDocument doc("mydocument");
	XmlQt::xmlParseFile(inputFile, &doc);


    // Get a pointer to the object manager instance
    ObjectManager* om = ObjectManager::instance();

    // Fill all the lists with empty pointers. 
    bool mustDelete = true;
    bool typeAttributeMustExist = false;

    om->FillAllListsFromXML(doc.documentElement(), typeAttributeMustExist, mustDelete);

	// Retrieve the CheprooPlusPlus element in the xml and the name of the problem
    QDomElement nodeCheprooPlusPlus = XmlQt::getFirstElemByTag(doc.documentElement(),XML_TAG_CHEPROOPLUSPLUS);
    QString aProblemName;
    XmlQt::getXMLAttribute(nodeCheprooPlusPlus, XML_ATTR_NAME, aProblemName); 

    bool CheprooPlusPlusExists = ! nodeCheprooPlusPlus.isNull();

    // Declare vector to comprare the results

    valarray<double> actCoeffDerivs;
    vector<double> upperturbationActCoeff;
    vector<double> downperturbationActCoeff;
    vector<double> upperturbationConc;
    vector<double> downperturbationConc;
    
    //now we perturbate the prescribed value field and compute derivatives by finite differences
    int werebad = 0;


    if( CheprooPlusPlusExists ) 
    {
	    // Create a Problem object and read its properties from the xml
        CCheprooPlusPlus* cheprooPlusPlusPointer = dynamic_cast<CCheprooPlusPlus*>(om->newInstanceOfClassName(CLASS_NAME_CHEPROOPLUSPLUS, aProblemName));

        // Read and initialize the problem
        cheprooPlusPlusPointer->ReadAndInitialize(nodeCheprooPlusPlus, false);

        // Pointer to GlobalChemicalSystem
        CLocalChemicalSystem* myLocalPointer = cheprooPlusPlusPointer->mGlobalChemSysPointer->mLocChemSysVector[0];

        // Call the function to speciate initial waters
        myLocalPointer->SpeciateInitialWater(cheprooPlusPlusPointer->mGlobalChemSysPointer->mLocChemSysVector[0]->mRelChemicalCompositions[0]);

        // Obtain a pointer to the aqueous phase and all other quantities I need to check
        CPhase* myAqPhase = cheprooPlusPlusPointer->mGlobalChemSysPointer->mAllPhases["aqueous1"];
        vector<CSpecies*> &aSV1 = myLocalPointer->mSpecies[eAqueous1nc];
        vector<CSpecies*> &aSV2 = myLocalPointer->mSpecies[eAqueous2];
        map<CSpecies*, int> &mSpeciesIndices = myLocalPointer->mSpeciesIndices;

        // Set the aqueous phase
        myAqPhase->Set(myLocalPointer->mRelChemicalCompositions.at(0));
        // Get secondary species global indices
        vector<int>& secondaryGlobIndices = myLocalPointer->GetSecondarySpcIndices();

        // Begin the verification of the function
        double mIonicStrength = 0.0;

        vector<double> aIonicStrengthDerivs;
        aIonicStrengthDerivs.resize(myLocalPointer->mSpeciesIndices.size());

        vector<double> aActCoeffsVector;
        aActCoeffsVector.resize(myLocalPointer->mSpeciesIndices.size());
        upperturbationActCoeff.resize(aActCoeffsVector.size());
        downperturbationActCoeff.resize(aActCoeffsVector.size());
        upperturbationConc.resize(aActCoeffsVector.size());
        downperturbationConc.resize(aActCoeffsVector.size());

        // Compute activity coefficients
        myAqPhase->ComputeActivityCoeff(myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration, 
                                        aActCoeffsVector, aSV1, aSV2, mSpeciesIndices, mIonicStrength);

        myAqPhase->ComputeIonicStrength(myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration, mIonicStrength,
                                        myLocalPointer->mSpecies[eAqueous1nc], myLocalPointer->mSpecies[eAqueous2],
                                        myLocalPointer->mSpeciesIndices);

        // Compute derivative of ionic strength wrt concentrations
        myAqPhase->ComputeIonicStrengthDerivsWrtConc(aIonicStrengthDerivs, aSV1, aSV2, mSpeciesIndices);

        // Matrix where I put dGamma/dc1, which is the variable to compare
        vector<vector<double> > dGammadc1;

        // First I evaluate dGamma1/dc1

        // First loop over primary species
        dGammadc1.resize(myLocalPointer->mSpeciesIndices.size());
        for(int i=0; i!=dGammadc1.size(); i++)
        {
            dGammadc1[i].resize(myLocalPointer->mSpecies[eAqueous1nc].size());
        }

        for(int i=0; i!=aSV1.size(); i++)
        {
            
            int SV1Index1 = mSpeciesIndices[aSV1[i]];
            double aIonSize = aSV1[i]->mIonicRadius;
            double aCharge = aSV1[i]->mIonicCharge;
                
            double dGammadI = 2*sqrt(mIonicStrength)*(1.0 + aIonSize*0.328299999999266*sqrt(mIonicStrength))*(1.0 + aIonSize*0.328299999999266*sqrt(mIonicStrength));
            dGammadI = (-(0.508930000000984*aCharge*aCharge/dGammadI) + 0.0409999999984375) * aActCoeffsVector[SV1Index1] * log(10.0);

            // Second loop over primary species
            for(int j=0; j!=aSV1.size(); j++)
            {
                int SV1Index2 = mSpeciesIndices[aSV1[j]];

                dGammadc1[SV1Index1][SV1Index2] = dGammadI * aIonicStrengthDerivs[SV1Index2];
            }
        }

        // Then I evaluate dGamma2/dc1

        // First loop over secondary species
        for(int i=0; i!=aSV2.size(); i++)
        {
            int SV2Index = mSpeciesIndices[aSV2[i]];
            double aIonSize = aSV2[i]->mIonicRadius;
            double aCharge = aSV2[i]->mIonicCharge;

            double dGammadI = 2*sqrt(mIonicStrength)*(1.0 + aIonSize*0.328299999999266*sqrt(mIonicStrength))*(1.0 + aIonSize*0.328299999999266*sqrt(mIonicStrength));
            dGammadI = (-(0.508930000000984*aCharge*aCharge/dGammadI) + 0.0409999999984375) * aActCoeffsVector[SV2Index] * log(10.0);

            // Second loop over primary species
            for(int j=0; j!=aSV1.size(); j++)
            {
                int SV1Index = mSpeciesIndices[aSV1[j]];

                dGammadc1[SV2Index][SV1Index] = dGammadI * aIonicStrengthDerivs[SV1Index];
            }
        } 

        vector<double> delta;
        delta.resize(myLocalPointer->mSpecies[eAqueous1nc].size());

        // Copy concentration vector in upperturbationConc and downperturbationConc
        for(unsigned int i=0; i< upperturbationConc.size(); i++)
        {
            upperturbationConc[i] = myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration[i];
            downperturbationConc[i] = myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration[i];
        }

        for (unsigned primSp = 0; primSp < myLocalPointer->mSpecies[eAqueous1nc].size(); primSp++)
        {
            // define perturbation size
            delta[primSp] = 0.001;
            if (  myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration[primSp] != 0)
            {
                delta[primSp] = abs(myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration[primSp]) *0.001;
            }
            
            // perturbate positive and get back the original value of the previous specues
            upperturbationConc[primSp] = myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration[primSp] + delta[primSp];
            if(primSp!=0) upperturbationConc[primSp-1] = myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration[primSp-1];

            // perturbate negative
            downperturbationConc[primSp] = myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration[primSp] - delta[primSp];            
            if(primSp!=0) downperturbationConc[primSp-1] = myLocalPointer->mRelChemicalCompositions.at(0)->mConcentration[primSp-1];

            // Compute activity coefficients with perturbations
            myAqPhase->ComputeActivityCoeff(upperturbationConc, upperturbationActCoeff, aSV1, aSV2, mSpeciesIndices, mIonicStrength);
            myAqPhase->ComputeActivityCoeff(downperturbationConc, downperturbationActCoeff, aSV1, aSV2, mSpeciesIndices, mIonicStrength);
       
            for (unsigned ithsp = 0; ithsp < aActCoeffsVector.size(); ithsp++)
            {
                double analyDeriv = dGammadc1[ithsp][primSp];
                double numderiv = (upperturbationActCoeff[ithsp] - aActCoeffsVector[ithsp])/delta[primSp];
                double HOT = (0.5 / delta[primSp]) *( upperturbationActCoeff[ithsp]  - 2*aActCoeffsVector[ithsp] + downperturbationActCoeff[ithsp]);
                double maxdif = abs( 1.1* HOT)+ 1e-7;
                double diff =  abs( analyDeriv - numderiv);

                if ( diff > maxdif)
                {
                    werebad++;
                }
            }
        }
    }

    return werebad;
}

void UTCAqueousDebyeHuckel()
{
    const char* testName = "UTCAqueousDebyeHuckel";

    CUnitTest::RunBasicUT(testName, UTCAqueousDebyeHuckel_1, 1);
    CUnitTest::RunBasicUT(testName, UTCAqueousDebyeHuckel_2, 2);

    //CUnitTest::RunExceptionUT(testName, UTCProost_1, 1, ERROR_UNINITIALIZED_XML_NODE);
}
