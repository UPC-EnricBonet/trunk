#include "clocalchemicalsystem.h"
#include "cglobalchemicalsystem.h"

// Declare useful iterators over Phases, used almost in every method of the class
multimap<ePhaseKind,CPhase*>::iterator it;
pair<multimap<ePhaseKind,CPhase*>::iterator,multimap<ePhaseKind,CPhase*>::iterator> it1;


CLocalChemicalSystem::CLocalChemicalSystem(QDomElement aNode)
: CCheprooBase(aNode)
{
    this->mClassName = "CLocalChemicalSystem";
    mGlobChemSysPointer = 0;
    this->mEqReactions.clear();
    this->mPrimRedSpecies.clear();
    this->mPhases.clear();
    this->mRelPrimarySpcIndices.clear();
    this->mRelPrimSpcIndicesWithoutCA.clear();
    this->mRelSecondarySpcIndices.clear();
    this->mRelConstActSpcIndices.clear();
    this->mRelChemicalCompositions.clear();
    this->mSpeciesIndices.clear();
    this->mSpecies.clear();
    this->mMaxIterNum = 0;
	this->mConstActSpeciesNum = 0;
    this->mSpeciateWithNR = false;


}

CLocalChemicalSystem::CLocalChemicalSystem()
{
    this->mClassName = "CLocalChemicalSystem";
    mGlobChemSysPointer = 0;
    this->mEqReactions.clear();
    this->mPrimRedSpecies.clear();
    this->mPhases.clear();
    this->mRelPrimarySpcIndices.clear();
    this->mRelPrimSpcIndicesWithoutCA.clear();
    this->mRelSecondarySpcIndices.clear();
    this->mRelConstActSpcIndices.clear();
    this->mRelChemicalCompositions.clear();
    this->mSpeciesIndices.clear();
    this->mSpecies.clear();
    this->mMaxIterNum = 0;
	this->mConstActSpeciesNum = 0;
    this->mSpeciateWithNR = false;

}

CLocalChemicalSystem::~CLocalChemicalSystem()
{
    this->mChemicalCompositions.clear();
}

eSpcToPhaseKind CLocalChemicalSystem::str2PhaseKind(QString s)
{
	if ( s == "aqueousPrimary" ) return eAqueous1;
	if ( s == "aqueousSecondary" ) return eAqueous2;	
    if ( s == "mineralPrimary" ) return eMineral1;
    if ( s == "mineralSecondary" ) return eMineral2;
	if ( s == "gasPrimary" ) return eGas1;
	if ( s == "gasSecondary" ) return eGas2;
    if ( s == "surfacePrimary" ) return eSurface1;
    if ( s == "surfaceSecondary" ) return eSurface2;
	throw GeneralException(ERROR_SPCTOPHASE_KIND_NAME);
}

QString CLocalChemicalSystem::phaseKind2str(eSpcToPhaseKind e)
{
	if ( e == eAqueous1 ) return "aqueousPrimary";
	if ( e == eAqueous2) return "aqueousSecondary";
	if ( e == eMineral1 ) return "mineralPrimary";
    if ( e == eMineral2 ) return "mineralSecondary";
	if ( e == eGas1 ) return "gasPrimary";
	if ( e == eGas2) return "gasSecondary";
	if ( e == eSurface1 ) return "surfacePrimary";
    if ( e == eSurface2 ) return "surfaceSecondary";
	throw GeneralException(ERROR_SPCTOPHASE_KIND_NUMBER);
}

void 
CLocalChemicalSystem::Read(const QDomElement aNode, vector<string>& primSpecies, set<string>& cas)
{
    //// Read references to kinetic reactions, if there are any
    //vector<QString> listOfIDs;
    //vector<QDomElement> reactionTags;
    //this->ReadReferences(aNode,
    //                     XML_TAG_KINETICREACTION, 
    //                     reactionTags, 
    //                     &listOfIDs);

    // Assign the name of the LocalChemicalSystem
    QString aLocChemSysName;
    XmlQt::getXMLAttribute(aNode, XML_ATTR_NAME, aLocChemSysName); 
    this->SetName(aLocChemSysName);

    // Read convergence parameters

    string solverType;
	XmlQt::getXMLAttribute(aNode, XML_ATTR_TYPESOLVER, solverType);
    if(solverType == "N-R") mSpeciateWithNR = true; // Picard default

    XmlQt::getXMLAttribute(aNode, XML_ATTR_MAXRELATIVEERROR, mMaxRelativeError); 
    XmlQt::getXMLAttribute(aNode, XML_ATTR_MAXRESIDUAL, mMaxResidual); 
    XmlQt::getXMLAttribute(aNode, XML_ATTR_MAXITERNUM, mMaxIterNum); 

    // Set local attributes (primary, secondary, indices, matrices...)

    this->SetLocal(primSpecies, cas);

}



void 
CLocalChemicalSystem::ComputeLogkMod()
{
    this->mLogkMod = VectorXd::Zero(mEqReactions.size());

    VectorXd aLogK = VectorXd::Zero(mEqReactions.size());

	// If S2 has been built, then mLogKMod = (1/S2) * logK
	if(this->mS2Inverse.size() != 0)
	{
		for ( int i = 0 ; i != mRelSecondarySpcIndices.size() ; i++ ) 
		{
			// Obtain pointer to current reaction
			CReaction* aCurrReaction = mEqReactions.at(i);

            // Build the vector of logK
			double temp = 25;
			aLogK(i) = aCurrReaction->mTempVsLogK[temp];
		}

        // Evaluate mLogKMod
        this->mLogkMod = this->mS2Inverse * aLogK;
	}
	else
	{
        this->mS2Inverse = MatrixXd::Zero(mEqReactions.size(), mEqReactions.size());

        // Loop over secondary species indices (the species must be in order of index in S2)
        for(int i=0; i!=mRelSecondarySpcIndices.size(); i++)
        {
            // Get the species correspondent to the index
            for(map<CSpecies*, int>::iterator it= this->mSpeciesIndices.begin() ; it != this->mSpeciesIndices.end() ; it++)
            {
                if(mSpeciesIndices[it->first] == mRelSecondarySpcIndices[i])
                {
		            for ( int j = 0 ; j != mEqReactions.size() ; j++ ) 
		            {
			            // Obtain pointer to current reaction
			            CReaction* aCurrReaction = mEqReactions.at(j);

                        // Obtain stoichiometric coefficient for the current species in the current reaction
			            double stoichCoeff = aCurrReaction->mRelSpeciesVsStoichCoeffs[it->first];

			            // Put it in the matrix
			            this->mS2Inverse(j,i) = stoichCoeff;

                        // Evaluate mLogKMod
						double temp = 25;
			            aLogK(j) = aCurrReaction->mTempVsLogK[temp];
		            }
                }
            }
	    }

        this->mS2Inverse = this->mS2Inverse.inverse();

        // Evaluate mLogKMod
        this->mLogkMod = this->mS2Inverse * aLogK;
    }

   aLogK.resize(0);
}


void
CLocalChemicalSystem::Speciate()
{
	// Declare vector of c1 and deltac1
    VectorXd deltac1 = VectorXd::Zero(this->mPrimSpRedNum);
    VectorXd c1 = VectorXd::Zero(this->mPrimSpRedNum);

    // Loop over all ChemicalCompositions
    for(int aChemCompNum=0; aChemCompNum!=mRelChemicalCompositions.size(); aChemCompNum++)
    {
        CChemicalComposition* aCurrChemComp = this->mRelChemicalCompositions[aChemCompNum];

        // Main N-R loop over global iterations
        for(int iter = 0, itHasConvergedToc1 = false; iter < mMaxIterNum && !itHasConvergedToc1; iter++ )
        {
            // Get the initial guess of primary species (without constant activity species) concentrations from the current ChemicalComposition
            aCurrChemComp->GetSetOfConcVector(c1, this->mRelPrimSpcIndicesWithoutCA);
            if(IsZeroVector(c1)) cout<<"Primary concentrations can't be zero";
                                        
            // Evaluate the concentrations of aqueous secondary species at iter+1 by means of N-R
            this->ComputeSecondaries(aCurrChemComp, c1, true, iter);

            // Compute the residual -fu(iter) = -U1*c1-U2*c2+u
            this->ComputeResidualOfComponents(aCurrChemComp, c1, this->mComponentMatrixInit);

            // Evaluate the Jacobian
            this->ComputeJacobian(this->mComponentMatrixInit);

            // Solve the system (evaluate deltac1)
            deltac1 = jacobian.fullPivLu().solve(residual);

            // Update the solution
			c1 = c1 + deltac1;

            // Check convergence
            itHasConvergedToc1 = this->CheckConvergence(c1, deltac1, residual);  

            // Update the concentrations of the current chemical composition if it has converged
            aCurrChemComp->UpdateConcentrations(c1, this->mRelPrimSpcIndicesWithoutCA);

            // Assign water molality to all chemical composition of the CLocalChemicalSystem (by now here, change it later)
            if(iter == 0 && this->mSpecies[eAqueous1c].size() != 0) 
            {
                double waterMolality = 1/H2O_mol_mass;
                aCurrChemComp->UpdateConcentration(waterMolality, mSpeciesIndices[this->mSpecies[eAqueous1c].at(0)]);
            }

            //cout << iter; 
        }

        // Assign water molality to all chemical composition of the CLocalChemicalSystem (by now here, change it later)


		// Calculate the concentration of non-eliminated species of other phases

		// - Select the gas phases
		it1 = this->mPhases.equal_range(eGasPhase);  

		// Compute c2 gas
		for (it=it1.first; it!=it1.second; it++)
		{
			//this->ComputeSecondaries(aCurrChemComp, primConcVector, false, iter, it->second);
		}

		// - Select all the surface phases
		it1 = this->mPhases.equal_range(eSurfacePhase);   

		// Compute c2 surface
		for (it=it1.first; it!=it1.second; it++)
		{
			//this->ComputeSecondaries(aCurrChemComp, primConcVector, false, iter, it->second);
		}

    }
    c1.resize(0);
    deltac1.resize(0);
}


void
CLocalChemicalSystem::ComputeSecondaries(CChemicalComposition* aCurrentChemicalComposition, VectorXd &c1, bool mustCalcDerivs, int iterc1)
{
    if(this->mSecSpeciesNum!=0) // Only if there are secondary species (if there are reactions)
    {
        // Declare and resize the vector of secondary species
        VectorXd c2 = VectorXd::Zero(this->mSecSpeciesNum);

        // Declare the iterator used to access the multimap containing the phases
        multimap<ePhaseKind,CPhase*>::iterator it;
        pair<multimap<ePhaseKind,CPhase*>::iterator,multimap<ePhaseKind,CPhase*>::iterator> it1;

        // Initialize mSystemMatrix, mSystemRHS and dc2_dc1
        this->dc2_dc1 = MatrixXd::Zero(this->mSecSpeciesNum,this->mPrimSpRedNum);
        this->mSystemRHS = MatrixXd::Zero(this->mSecSpeciesNum, this->mPrimSpRedNum);
        this->mSystemMatrix = MatrixXd::Identity(this->mSecSpeciesNum, this->mSecSpeciesNum);

        if (this->mSpeciateWithNR) // Solve with Newton-Raphson
        {
            // Select the aqueous phase
            it1 = this->mPhases.equal_range(eAqueousPhase);               
            CPhase* acurrPhase = it1.first->second;

            // Set the aqueous phase (i.e. dependent-temperature coefficients)
            acurrPhase->Set(aCurrentChemicalComposition);

            // Compute first estimation of activity coefficients
            acurrPhase->ComputeActivityCoeff(aCurrentChemicalComposition->mConcentration, 
                                             aCurrentChemicalComposition->mActivityCoeff,
                                             mSpecies[eAqueous1nc], 
                                             mSpecies[eAqueous2],
                                             mSpeciesIndices
                                             );

            // Select the mineral phase
            it1 = this->mPhases.equal_range(eMineralPhase);  
            if(it1.first!=this->mPhases.end())
            {
                acurrPhase = it1.first->second;

                acurrPhase->ComputeActivityCoeff(aCurrentChemicalComposition->mConcentration, 
                                                 aCurrentChemicalComposition->mActivityCoeff,
                                                 mSpecies[eMineral1], 
                                                 mSpecies[eMineral2],
                                                 mSpeciesIndices
                                                 );
            }

            // Select the gas phase
            it1 = this->mPhases.equal_range(eGasPhase); 
            if(it1.first!=this->mPhases.end())
            {
                acurrPhase = it1.first->second;
                acurrPhase->ComputeActivityCoeff(aCurrentChemicalComposition->mConcentration, 
                                                 aCurrentChemicalComposition->mActivityCoeff,
                                                 mSpecies[eGas1], 
                                                 mSpecies[eGas2],
                                                 mSpeciesIndices
                                                 );
            }

            // Select the surface phase
            it1 = this->mPhases.equal_range(eSurfacePhase); 
            if(it1.first!=this->mPhases.end())
            {
                acurrPhase = it1.first->second;
                acurrPhase->ComputeActivityCoeff(aCurrentChemicalComposition->mConcentration, 
                                                 aCurrentChemicalComposition->mActivityCoeff,
                                                 mSpecies[eSurface1], 
                                                 mSpecies[eSurface2],
                                                 mSpeciesIndices
                                                 );
            }

            // Check if some secondary species have some constraint (because they were defined as primary at first for example. Only for eqmin or activity though)
            this->ApplyConstraints(aCurrentChemicalComposition, c1, false);

            // Transform c1 in ln(c1)
            VectorXd lnc1 = VectorXd::Zero(this->mPrimSpRedNum);
            for(int i=0; i!=this->mPrimSpRedNum; i++)
            {
                lnc1(i) = log(c1(i));
            }

            // Get set of ln of secondary concentrations of the current chemical composition
            aCurrentChemicalComposition->GetSetOfLnConcVector(c2, mRelSecondarySpcIndices);

            // At first iteration, initialize secondaries
            this->InitializeSecondaries(c2, lnc1, aCurrentChemicalComposition);

            for(int iterNR = 0, itHasConverged = false; iterNR < mMaxIterNum && !itHasConverged; iterNR++ )
            {
 
                // Initialize mSystemMatrix
                mSystemMatrix = MatrixXd::Identity(this->mSecSpeciesNum, this->mSecSpeciesNum);

                // Compute derivatives of activity coefficients wrt concentrations for the aqueous phase and add their contribution to the jacobian (mSystemMatrix)
                acurrPhase->ComputeActivityCoeffDerivs(aCurrentChemicalComposition->mConcentration, mSystemMatrix, this->mS1Mod, 
                                                       mSpecies[eAqueous1nc], mSpecies[eAqueous2], 
                                                       mSpeciesIndices, true);

                // Select the mineral phase
                it1 = this->mPhases.equal_range(eMineralPhase);  
                if(it1.first!=this->mPhases.end())
                {
                    acurrPhase = it1.first->second;

                    // Compute derivatives of activity coefficients wrt concentrations for the mineral phase and add their contribution to the jacobian (mSystemMatrix)
                    acurrPhase->ComputeActivityCoeffDerivs(aCurrentChemicalComposition->mConcentration, mSystemMatrix, this->mS1Mod, 
                                                           mSpecies[eMineral1], mSpecies[eMineral2], 
                                                           mSpeciesIndices, true);
                }

                // Select the gas phase
                it1 = this->mPhases.equal_range(eGasPhase);  
                if(it1.first!=this->mPhases.end())
                {
                    acurrPhase = it1.first->second;

                    // Compute derivatives of activity coefficients wrt concentrations for the gas phase and add their contribution to the jacobian (mSystemMatrix)
                    acurrPhase->ComputeActivityCoeffDerivs(aCurrentChemicalComposition->mConcentration, mSystemMatrix, this->mS1Mod, 
                                                           mSpecies[eGas1], mSpecies[eGas2], 
                                                           mSpeciesIndices, true);
                }

                // Select the surface phase
                it1 = this->mPhases.equal_range(eSurfacePhase);  
                if(it1.first!=this->mPhases.end())
                {
                    acurrPhase = it1.first->second;
                    // Compute derivatives of activity coefficients wrt concentrations for the gas phase and add their contribution to the jacobian (mSystemMatrix)
                    acurrPhase->ComputeActivityCoeffDerivs(aCurrentChemicalComposition->mConcentration, mSystemMatrix, this->mS1Mod, 
                                                           mSpecies[eSurface1], mSpecies[eSurface2], 
                                                           mSpeciesIndices, true);
                }

                // Build the residual (-fMAL)
                VectorXd residual = VectorXd::Zero(this->mSecSpeciesNum);
                this->ComputeResidualOfMAL(residual, lnc1, c2, aCurrentChemicalComposition);

                // Solve the system
                VectorXd deltac2 = VectorXd::Zero(this->mSecSpeciesNum);
                deltac2 = mSystemMatrix.fullPivLu().solve(residual);

                // Update the solution and check convergence
                c2 = c2 + deltac2;

                // Check convergence
                itHasConverged = this->CheckConvergence(c2, deltac2, residual);  

                // Update the concentrations of the current chemical composition if it has converged
                aCurrentChemicalComposition->UpdateConcentrationsFromLn(c2, mRelSecondarySpcIndices);

                // Evaluate activity coefficients
			    acurrPhase->ComputeActivityCoeff(aCurrentChemicalComposition->mConcentration, 
											     aCurrentChemicalComposition->mActivityCoeff,
											     mSpecies[eAqueous1nc], 
											     mSpecies[eAqueous2],
											     mSpeciesIndices
											     );
                // Select the mineral phase
                it1 = this->mPhases.equal_range(eMineralPhase);  
                if(it1.first!=this->mPhases.end())
                {
                    acurrPhase = it1.first->second;

                    acurrPhase->ComputeActivityCoeff(aCurrentChemicalComposition->mConcentration, 
                                                 aCurrentChemicalComposition->mActivityCoeff,
                                                 mSpecies[eMineral1], 
                                                 mSpecies[eMineral2],
                                                 mSpeciesIndices
                                                 );
                }

                // Select the gas phase
                it1 = this->mPhases.equal_range(eGasPhase);  
                if(it1.first!=this->mPhases.end())
                {
                    acurrPhase->ComputeActivityCoeff(aCurrentChemicalComposition->mConcentration, 
                                                 aCurrentChemicalComposition->mActivityCoeff,
                                                 mSpecies[eGas1], 
                                                 mSpecies[eGas2],
                                                 mSpeciesIndices
                                                 );
                }

                // Select the surface phase
                it1 = this->mPhases.equal_range(eSurfacePhase);  
                if(it1.first!=this->mPhases.end())
                {
                    acurrPhase->ComputeActivityCoeff(aCurrentChemicalComposition->mConcentration, 
                                                 aCurrentChemicalComposition->mActivityCoeff,
                                                 mSpecies[eSurface1], 
                                                 mSpecies[eSurface2],
                                                 mSpeciesIndices
                                                 );
                }

            }

        }
        else // Solve with Picard
        {
            // Declare and resize the vector of secondary species
            VectorXd secondaryConcentrations = VectorXd::Zero(this->mSecSpeciesNum);

            // Declare deltac2
            VectorXd deltac2 = VectorXd::Zero(this->mSecSpeciesNum);

            // Select the aqueous phase
            it1 = this->mPhases.equal_range(eAqueousPhase);               
            CPhase* aAqPhase = it1.first->second;

            // Set the aqueous phase (i.e. dependent-temperature coefficients)
            aAqPhase->Set(aCurrentChemicalComposition);

            // Compute first estimation of activity coefficients
            aAqPhase->ComputeActivityCoeff(aCurrentChemicalComposition->mConcentration, 
                                             aCurrentChemicalComposition->mActivityCoeff,
                                             mSpecies[eAqueous1nc], 
                                             mSpecies[eAqueous2],
                                             mSpeciesIndices
                                             );

            // Select the mineral phase
            it1 = this->mPhases.equal_range(eMineralPhase);  
            if(it1.first!=this->mPhases.end())
            {
                CPhase* acurrPhase = it1.first->second;

                acurrPhase->ComputeActivityCoeff(aCurrentChemicalComposition->mConcentration, 
                                                 aCurrentChemicalComposition->mActivityCoeff,
                                                 mSpecies[eMineral1], 
                                                 mSpecies[eMineral2],
                                                 mSpeciesIndices
                                                 );
            }

            // Select the gas phase
            it1 = this->mPhases.equal_range(eGasPhase);  
            if(it1.first!=this->mPhases.end())
            {
                CPhase* gasPhase = it1.first->second;

                gasPhase->ComputeActivityCoeff(aCurrentChemicalComposition->mConcentration, 
                                               aCurrentChemicalComposition->mActivityCoeff,
                                               mSpecies[eGas1], 
                                               mSpecies[eGas2],
                                               mSpeciesIndices
                                               );
            }

            // Select the surface phase
            it1 = this->mPhases.equal_range(eSurfacePhase);  
            if(it1.first!=this->mPhases.end())
            {
                CPhase* gasPhase = it1.first->second;

                gasPhase->ComputeActivityCoeff(aCurrentChemicalComposition->mConcentration, 
                                               aCurrentChemicalComposition->mActivityCoeff,
                                               mSpecies[eSurface1], 
                                               mSpecies[eSurface2],
                                               mSpeciesIndices
                                               );
            }
            
            // Check if some secondary species have some constraint (because they were defined as primary at first for example. Only for eqmin or activity though)
            this->ApplyConstraints(aCurrentChemicalComposition, c1, false);

            // Get vector of activity coefficients of primary species for the current ChemicalComposition
            VectorXd primaryActCoeffs = VectorXd::Zero(this->mPrimSpRedNum);
            aCurrentChemicalComposition->GetSetOfActivityCoeffVector(primaryActCoeffs, mRelPrimSpcIndicesWithoutCA);

            // Get vector of concentration of primary species for the current ChemicalComposition
            VectorXd c1_c = VectorXd::Zero(this->mConstActSpeciesNum);
            aCurrentChemicalComposition->GetSetOfConcVector(c1_c, mRelConstActSpcIndices);

            // Get vector of activity coefficients of primary species for the current ChemicalComposition
            VectorXd primaryActCoeffs_c = VectorXd::Zero(this->mConstActSpeciesNum);
            aCurrentChemicalComposition->GetSetOfActivityCoeffVector(primaryActCoeffs_c, mRelConstActSpcIndices);

		    // Declare a map with species of concentrations defined by the user (that the user defined as primary but became secondary): bool is true and indicates that c has been defined
		    map<int,bool> definedConc;

            // Declare and resize vector logc1g1
            VectorXd logc1g1 = VectorXd::Zero(primaryActCoeffs.size());
            VectorXd logc1cg1c = VectorXd::Zero(c1_c.size());


            for(int k=0; k!=c1_c.size(); k++)
		    {
			    logc1cg1c(k) = log10(c1_c(k)*(primaryActCoeffs_c(k)));
		    }

            for(int iterPicard = 0, itHasConverged = 0; iterPicard < mMaxIterNum && !itHasConverged; iterPicard++ )
            {

                for(int j=0; j!=c1.size(); j++)
			    {
				    logc1g1(j) = log10(c1(j)*(primaryActCoeffs(j)));
			    }
                for(int i=0; i!=this->mSecSpeciesNum; i++)
                {
                    // Index of the secondary variable
                    int secSpeciesGloblIndex = mRelSecondarySpcIndices[i];

                    // Get old value of the concentration
                    double c2old = aCurrentChemicalComposition->mConcentration.at(secSpeciesGloblIndex);

				    if(iterPicard == 0 && c2old != 0 && iterc1 == 0) // If the concentration has been previously defined by the user
				    {
                        if(aCurrentChemicalComposition->iCon.size()!=0)
                        {
    					    // Put the concentration defined by the user in the vector of secondary concentrations
    					    if(aCurrentChemicalComposition->constrValue(secSpeciesGloblIndex)!=0) 
                            {
                                secondaryConcentrations(i) = aCurrentChemicalComposition->constrValue(secSpeciesGloblIndex);
                                definedConc.insert(pair<int, bool>(secSpeciesGloblIndex,true));
                            }
                            else
                            {
    					        // Put the concentration defined by the user in the vector of secondary concentrations
    					        secondaryConcentrations(i) = aCurrentChemicalComposition->mConcentration.at(secSpeciesGloblIndex);
                            }
                            
                        }
                        else
                        {
                            secondaryConcentrations(i) = aCurrentChemicalComposition->mConcentration.at(secSpeciesGloblIndex);
                        }
				    }
				    else
				    {
                        if(!definedConc[secSpeciesGloblIndex])
					    {
                            double c2 = 0.0;

						    // Evaluate new value of c2
                            if(this->mSecMod.size()!=0)
                            {
						        c2 = this->mLogkMod(i) + this->mS1Mod.row(i).dot(logc1g1) + this->mSecMod.row(i).dot(logc1cg1c) - log10(aCurrentChemicalComposition->mActivityCoeff.at(secSpeciesGloblIndex));
                            }
                            else
                            {
                                c2 = this->mLogkMod(i) + this->mS1Mod.row(i).dot(logc1g1) - log10(aCurrentChemicalComposition->mActivityCoeff.at(secSpeciesGloblIndex));
                            }
						    c2 = pow(10.0, c2);

						    // Update vector of secondary concentrations and deltac2
						    secondaryConcentrations(i) = c2;
						    deltac2(i) = c2old - c2;
					    }

				    }


                }

			    // Update the concentrations of the current chemical composition if it has converged
			    aCurrentChemicalComposition->UpdateConcentrations(secondaryConcentrations, mRelSecondarySpcIndices);

			    // Evaluate activity coefficients
			    aAqPhase->ComputeActivityCoeff(aCurrentChemicalComposition->mConcentration, 
											     aCurrentChemicalComposition->mActivityCoeff,
											     mSpecies[eAqueous1nc], 
											     mSpecies[eAqueous2],
											     mSpeciesIndices
											     );

                // Select the mineral phase
                it1 = this->mPhases.equal_range(eMineralPhase);  
                if(it1.first!=this->mPhases.end())
                {
                    CPhase* acurrPhase = it1.first->second;

                    acurrPhase->ComputeActivityCoeff(aCurrentChemicalComposition->mConcentration, 
                                                     aCurrentChemicalComposition->mActivityCoeff,
                                                     mSpecies[eMineral1], 
                                                     mSpecies[eMineral2],
                                                     mSpeciesIndices
                                                     );
                }

                // Select the gas phase
                it1 = this->mPhases.equal_range(eGasPhase);  
                if(it1.first!=this->mPhases.end())
                {
                    CPhase* gasPhase = it1.first->second;

                    gasPhase->ComputeActivityCoeff(aCurrentChemicalComposition->mConcentration, 
                                                   aCurrentChemicalComposition->mActivityCoeff,
                                                   mSpecies[eGas1], 
                                                   mSpecies[eGas2],
                                                   mSpeciesIndices
                                                   );
                }

                // Select the surface phase
                it1 = this->mPhases.equal_range(eSurfacePhase);  
                if(it1.first!=this->mPhases.end())
                {
                    CPhase* gasPhase = it1.first->second;

                    gasPhase->ComputeActivityCoeff(aCurrentChemicalComposition->mConcentration, 
                                                   aCurrentChemicalComposition->mActivityCoeff,
                                                   mSpecies[eSurface1], 
                                                   mSpecies[eSurface2],
                                                   mSpeciesIndices
                                                   );
                }

                // Check convergence
                VectorXd residual = VectorXd::Zero(this->mSecSpeciesNum);
                itHasConverged = this->CheckConvergence(secondaryConcentrations, deltac2, residual);


            }
            logc1g1.resize(0);
        }


        // Compute also derivatives if necessary

        if(mustCalcDerivs) this->Compute_dc2_dc1(aCurrentChemicalComposition, c1);

        c2.resize(0);
    }
}
void
CLocalChemicalSystem::Compute_dc2_dc1(CChemicalComposition* aCurrentChemicalComposition, VectorXd& c1)
{
    if(this->dc2_dc1.size() != 0) dc2_dc1.resize(0,0);

    // Declare and resize the vector of secondary species
    VectorXd c2 = VectorXd::Zero(this->mSecSpeciesNum);

    // Initialize the right hand side matrix (equal to S1*)
    this->mSystemRHS = this->mS1Mod;

    aCurrentChemicalComposition->GetSetOfConcVector(c2,mRelSecondarySpcIndices);

    // Select the aqueous phase
    it1 = this->mPhases.equal_range(eAqueousPhase);               
    CPhase* acurrPhase = it1.first->second;

    if(!this->mSpeciateWithNR)
    {
        // Initialize mSystemMatrix
        this->mSystemMatrix = MatrixXd::Identity(this->mSecSpeciesNum, this->mSecSpeciesNum);

        // Compute dlnGamma/dlnc2 and add their contribution to the jacobian (mSystemMatrix)
	    acurrPhase->ComputeActivityCoeffDerivs(aCurrentChemicalComposition->mConcentration, mSystemMatrix, this->mS1Mod, 
										    mSpecies[eAqueous1nc], mSpecies[eAqueous2], 
										    mSpeciesIndices, true);

        // Select the gas phase
        it1 = this->mPhases.equal_range(eGasPhase); 
        if(it1.first!=this->mPhases.end())
        {
            CPhase* gasPhase = it1.first->second;

            // Compute dlnGamma/dlnc2 and add their contribution to the jacobian (mSystemMatrix)
	        gasPhase->ComputeActivityCoeffDerivs(aCurrentChemicalComposition->mConcentration, mSystemMatrix, this->mS1Mod, 
										        mSpecies[eAqueous1nc], mSpecies[eAqueous2], 
										        mSpeciesIndices, true);
        }
    }
  
    //Compute dlnGamma/dlnc1 and add their contribution to the right hand side (mSystemRHS)
    acurrPhase->ComputeActivityCoeffDerivs(aCurrentChemicalComposition->mConcentration, mSystemRHS, this->mS1Mod, 
                                            mSpecies[eAqueous1nc], mSpecies[eAqueous2], 
                                            mSpeciesIndices, false);


    // Select the gas phase
    it1 = this->mPhases.equal_range(eGasPhase); 
    if(it1.first!=this->mPhases.end())
    {
        CPhase* gasPhase = it1.first->second;

        // Compute dlnGamma/dlnc2 and add their contribution to the jacobian (mSystemMatrix)
	    gasPhase->ComputeActivityCoeffDerivs(aCurrentChemicalComposition->mConcentration, mSystemMatrix, this->mS1Mod, 
										    mSpecies[eAqueous1nc], mSpecies[eAqueous2], 
										    mSpeciesIndices, true);
    }

	// Solve the system 
	dc2_dc1 = mSystemMatrix.fullPivLu().solve(mSystemRHS);	

	// Transform ln(c) to c
	for(int i=0; i!=this->mSecSpeciesNum; i++)
	{
		for(int j=0; j!=this->mPrimSpRedNum; j++)
		{
			dc2_dc1(i,j) = dc2_dc1(i,j) * c2(i) / c1(j);
		}
	}

    c2.resize(0);
}
void
CLocalChemicalSystem::Speciate_uaq(CChemicalComposition* aCurrChemComp)
{
	// Convert molality (mol/kgw) to molarity (mol/m3l) for aqueous phase
	it1 = this->mPhases.equal_range(eAqueousPhase);          
    CPhase* aAqPhase = it1.first->second;

	//aAqPhase->ComputeMolarity(aCurrChemComp->mConcentration, mSpecies[eAqueous1nc], mSpecies[eAqueous2], 
    //                          mSpeciesIndices, aCurrChemComp);

	// Declare vector of c1 and deltac1
    VectorXd c1 = VectorXd::Zero(this->mPrimSpRedNum);
    VectorXd deltac1 = VectorXd::Zero(this->mPrimSpRedNum);        
                    
    // Get the initial guess of primary species (without constant activity species) concentrations from aCurrChemComp
	this->CalculateInitialGuess(aCurrChemComp, c1);

    // Main N-R loop over global iterations
    for(int iter = 0, itHasConvergedToc1 = false; iter < mMaxIterNum && !itHasConvergedToc1; iter++ )
    {
        // Evaluate the concentrations of aqueous secondary species at iter+1 by means of N-R
        this->ComputeSecondaries(aCurrChemComp, c1, true, iter);

        // Compute the residual -fu(iter) = -U1*c1-U2*c2+u
        this->ComputeResidualOfComponents(aCurrChemComp, c1, this->mComponentMatrix);

        // Evaluate the Jacobian
        this->ComputeJacobian(this->mComponentMatrix);

        // Solve the system (evaluate deltac1)
        deltac1 = this->jacobian.fullPivLu().solve(this->residual);

        // Update the solution
		c1 = c1 + deltac1;

        // Update the concentrations of the current chemical composition
        aCurrChemComp->UpdateConcentrations(c1, this->mRelPrimSpcIndicesWithoutCA);

        // Check convergence
        itHasConvergedToc1 = this->CheckConvergence(c1, deltac1, residual);  

        //cout << iter; 
    }

	// Convert molarity (mol/m3l) to molality (mol/kgw) for aqueous phase
	//aAqPhase->ComputeConcFromMolarity(aCurrChemComp->mConcentration, mSpecies[eAqueous1nc], mSpecies[eAqueous2], 
 //                             mSpeciesIndices, aCurrChemComp);

    // Assign water molality to all chemical composition of the CLocalChemicalSystem (by now here, change it later)
    if(this->mSpecies[eAqueous1c].size() != 0) 
    {
        double waterMolality = 1/H2O_mol_mass;
        aCurrChemComp->UpdateConcentration(waterMolality, mSpeciesIndices[this->mSpecies[eAqueous1c].at(0)]);
    }

    // Print chemical informations on an output file
    this->WriteChemCompInfos(aCurrChemComp, true);

    c1.resize(0);
    deltac1.resize(0);
}


void
CLocalChemicalSystem::Speciate_utot(CChemicalComposition* aCurrChemComp)
{
    // Convert concentration (mol/kgw) to molarity (mol/m3_phase) for all phases

    // - Aqueous
	it1 = this->mPhases.equal_range(eAqueousPhase);
    if(it1.first!=this->mPhases.end())
    {
        for (it=it1.first; it!=it1.second; it++)
        {
            CPhase* aAqPhase = it1.first->second;
		    //aAqPhase->ComputeMolarity(aCurrChemComp->mConcentration,mSpecies[eAqueous1nc], mSpecies[eAqueous2], 
      //                                mSpeciesIndices, aCurrChemComp);
        }
    }
    // - Mineral
	it1 = this->mPhases.equal_range(eMineralPhase);
    if(it1.first!=this->mPhases.end())
    {
        for (it=it1.first; it!=it1.second; it++)
        {
            CPhase* aMinPhase = it->second;
		    //aMinPhase->ComputeMolarity(aCurrChemComp->mConcentration,mSpecies[eMineral1], mSpecies[eMineral2], 
            //                             mSpeciesIndices, aCurrChemComp);
        }
    }

    // - Gas
	it1 = this->mPhases.equal_range(eGasPhase);
    for (it=it1.first; it!=it1.second; it++)
    {
        CPhase* aGasPhase = it1.first->second;
		//aGasPhase->ComputeMolarity(aCurrChemComp->mConcentration,mSpecies[eGas1], mSpecies[eGas2], 
        //                             mSpeciesIndices, aCurrChemComp);
    }
    // - Surface
	it1 = this->mPhases.equal_range(eSurfacePhase);
    if(it1.first!=this->mPhases.end())
    {
        for (it=it1.first; it!=it1.second; it++)
        {
            CPhase* aSurfPhase = it1.first->second;
		    //aSurfPhase->ComputeMolarity(aCurrChemComp->mConcentration,mSpecies[eSurface1], mSpecies[eSurface2], 
            //                             mSpeciesIndices, aCurrChemComp);
        }
    }

	// Declare vector of c1 and deltac1
    VectorXd c1 = VectorXd::Zero(this->mPrimSpRedNum);
    VectorXd deltac1 = VectorXd::Zero(this->mPrimSpRedNum);
                    
    // Get the initial guess of primary species (without constant activity species) concentrations from the current ChemicalComposition
    this->CalculateInitialGuess(aCurrChemComp, c1);

    // Main N-R loop over global iterations
    for(int iter = 0, itHasConvergedToc1 = false; iter < mMaxIterNum && !itHasConvergedToc1; iter++ )
    {
        // Evaluate the concentrations of aqueous secondary species at iter+1 by means of N-R
        this->ComputeSecondaries(aCurrChemComp, c1, true, iter);

        // Compute the residual -fu(iter) = -U1*c1-U2*c2+u
        this->ComputeResidualOfComponents_allPhases(aCurrChemComp, c1);

        // Evaluate the Jacobian
        this->ComputeJacobian_allPhases(aCurrChemComp);

        // Solve the system (evaluate deltac1)
        deltac1 = this->jacobian.fullPivLu().solve(this->residual);

        // Update the solution
		c1 = c1 + deltac1;

        // Check convergence
        itHasConvergedToc1 = this->CheckConvergence(c1, deltac1, residual);  

        // Update the concentrations of the current chemical composition (here because it's useful for the get_u's methods)
        aCurrChemComp->UpdateConcentrations(c1, this->mRelPrimSpcIndicesWithoutCA);

        //cout << iter; 
    }


    // Update the volumetric fraction of the aqueous phase
    //string nameAqPhase = it1.first->second->name().toStdString();
    //aCurrChemComp->UpdateVolumFraction(c1(this->mPrimSpRedNum), nameAqPhase);

	// Evaluate constant acivity species concentrations
	// this->CalculateCASConcentration(aCurrChemComp);

    // Convert molarity to concentration
    // - Aqueous
	it1 = this->mPhases.equal_range(eAqueousPhase);
    if(it1.first!=this->mPhases.end())
    {
        for (it=it1.first; it!=it1.second; it++)
        {
            CPhase* aAqPhase = it1.first->second;
 /*           aAqPhase->ComputeConcFromMolarity(aCurrChemComp->mConcentration,mSpecies[eAqueous1nc], mSpecies[eAqueous2], 
                                              mSpeciesIndices, aCurrChemComp);*/
        }
    }

    // Assign water molality to all chemical composition of the CLocalChemicalSystem (by now here, change it later)
    if(this->mSpecies[eAqueous1c].size() != 0) 
    {
        double waterMolality = 1/H2O_mol_mass;
        aCurrChemComp->UpdateConcentration(waterMolality, mSpeciesIndices[this->mSpecies[eAqueous1c].at(0)]);
    }

    c1.resize(0);
    deltac1.resize(0);
}

void
CLocalChemicalSystem::SpeciateInitialWater(CChemicalComposition* aCurrChemComp)
{
	// Declare vector of c1 and deltac1
    VectorXd c1 = VectorXd::Zero(this->mPrimSpRedNum);
    VectorXd deltac1 = VectorXd::Zero(this->mPrimSpRedNum);        
                    
    // Get the initial guess of primary species (without constant activity species) concentrations from aCurrChemComp
	this->CalculateInitialGuess(aCurrChemComp, c1);

	// Convert molality (mol/kgw) to molarity (mol/m3l) for aqueous phase
	it1 = this->mPhases.equal_range(eAqueousPhase);          
    CPhase* aAqPhase = it1.first->second;
	//aAqPhase->ComputeMolarity(aCurrChemComp->mConcentration,mSpecies[eAqueous1nc], mSpecies[eAqueous2], 
 //                             mSpeciesIndices, aCurrChemComp);

    // Get part of the "classic" component matrix relative to the reduced primary species
    MatrixXd compInitPrim = MatrixXd::Zero(this->mPrimSpRedNum,this->mPrimSpRedNum+this->mEqReactions.size());
    VectorXd row = VectorXd::Zero(this->mEqReactions.size());

    int index_classic = 0;
    for(int i=0; i!=this->mPrimSpRedNum; i++)
    {
       index_classic = this->mClassicComponents[this->mPrimRedSpecies.at(i)->name().toStdString()];

       compInitPrim(i,i) = 1.0;
       row = this->mComponentMatrixInit.row(index_classic).tail(this->mEqReactions.size());
       compInitPrim.row(i).tail(this->mEqReactions.size()) = row;
       
    }

    // Main N-R loop over global iterations
    for(int iter = 0, itHasConvergedToc1 = false; iter < mMaxIterNum && !itHasConvergedToc1; iter++ )
    {
        // Evaluate the concentrations of aqueous secondary species at iter+1 by means of N-R
        this->ComputeSecondaries(aCurrChemComp, c1, true, iter);

        // Compute the residual -fu(iter) = -U1*c1-U2*c2+u
        this->ComputeResidualOfComponents(aCurrChemComp, c1, compInitPrim, true);

        // Evaluate the Jacobian
        this->ComputeJacobian(compInitPrim);

        // Apply constraints
        this->ApplyConstraints(aCurrChemComp, c1, true);

        // Solve the system (evaluate deltac1)
        deltac1 = this->jacobian.fullPivLu().solve(this->residual);

        // Update the solution
		c1 = c1 + deltac1;

        // Update the concentrations of the current chemical composition
        aCurrChemComp->UpdateConcentrations(c1, this->mRelPrimSpcIndicesWithoutCA);

        // Check convergence
        itHasConvergedToc1 = this->CheckConvergence(c1, deltac1, residual);  

        //cout << iter; 
    }

    // Assign water molality to all chemical composition of the CLocalChemicalSystem (by now here, change it later)
    if(this->mSpecies[eAqueous1c].size() != 0) 
    {
        double waterMolality = 1/H2O_mol_mass;
        aCurrChemComp->UpdateConcentration(waterMolality, mSpeciesIndices[this->mSpecies[eAqueous1c].at(0)]);
    }

    c1.resize(0);
    deltac1.resize(0);
}


void 
CLocalChemicalSystem::ComputeResidualOfComponents(CChemicalComposition* aCurrChemComp, VectorXd &c1, MatrixXd &compMatrixPrim, bool isInizialization)
{    
    // Get secondary concentrations of the current chemical composition
    VectorXd c2 = VectorXd::Zero(this->mSecSpeciesNum);
    aCurrChemComp->GetSetOfConcVector(c2,mRelSecondarySpcIndices);

    // Set the residual to zero
    this->residual = VectorXd::Zero(this->mPrimSpRedNum);

    // Obtain constrValue vector for reduced primary species
    VectorXd cTotRed = VectorXd::Zero(this->mPrimSpRedNum);
    for(int i=0; i!=this->mPrimSpRedNum; i++)
    {
        cTotRed(i) = aCurrChemComp->constrValue(this->mRelPrimSpcIndicesWithoutCA[i]);
    }

    // Evaluate the residual
    if(isInizialization && this->mComponentMatrixInit.size()!=0)
    {
        this->residual = cTotRed - compMatrixPrim.block(0,0,this->mPrimSpRedNum,this->mPrimSpRedNum) * c1 - 
            compMatrixPrim.block(0,this->mPrimSpRedNum,this->mPrimSpRedNum,this->mSecSpeciesNum) * c2;
    }
    else
    {
        this->residual = cTotRed - this->mComponentMatrix.block(0,0,this->mPrimSpRedNum,this->mPrimSpRedNum) * c1 - 
                   this->mComponentMatrix.block(0,this->mPrimSpRedNum,this->mPrimSpRedNum,this->mSecSpeciesNum) * c2;
    }

    c2.resize(0);
    cTotRed.resize(0);
}

void 
CLocalChemicalSystem::ComputeResidualOfComponents_allPhases(CChemicalComposition* aCurrChemComp, VectorXd &c1)
{    
    // Get secondary concentrations of the current chemical composition
    VectorXd c2 = VectorXd::Zero(this->mSecSpeciesNum);
    aCurrChemComp->GetSetOfConcVector(c2,mRelSecondarySpcIndices);

    // Set the residual to zero
    this->residual = VectorXd::Zero(this->mPrimSpRedNum);

    // Obtain constrValue vector for reduced primary species
    VectorXd cTotRed = VectorXd::Zero(this->mPrimSpRedNum);
    for(int i=0; i!=this->mPrimSpRedNum; i++)
    {
        cTotRed(i) = aCurrChemComp->constrValue(this->mRelPrimSpcIndicesWithoutCA[i]);
    }

    // - First declare iterators and useful vector
    map<string, double>::iterator it3;
    vector<double> u;

        // - Access the aqueous phase
    it1 = this->mPhases.equal_range(eAqueousPhase); 

    // - Get the name of the phase
    string phaseName = it1.first->second->name().toStdString();

    // - Retrieve the volumetric fraction from CChemicalComposition
    it3 = aCurrChemComp->mVolumetricFractions.find(phaseName);
    double theta = it3->second;

    // - Get u_aq 
    u = this->Get_u_aq(aCurrChemComp);

    // - Add contribution to residual vector
    for(int i=0; i!=this->mPrimSpRedNum; i++)
    {
        this->residual(i) = u[i] * theta;
    }
    u.resize(0);
    theta = 0;

    // Add contribution of other phases

    // - Get u_min 
    u = this->Get_u_min(aCurrChemComp);

    if(u.size()!=0)
    {
        // - Access the mineral phase
        it1 = this->mPhases.equal_range(eMineralPhase); 

        if(it1.first!=this->mPhases.end()) 
        {
            for (it=it1.first; it!=it1.second; it++) // OSS: it won't work for more than one mineral phase...same for ads and gases
            {
                // - Get the name of the phase
                phaseName = it1.first->second->name().toStdString();

                // - Retrieve the volumetric fraction from CChemicalComposition
                it3 = aCurrChemComp->mVolumetricFractions.find(phaseName);
                theta = it3->second;
            }
        }
        // - Add contribution to residual vector
        for(int i=0; i!=this->mPrimSpRedNum; i++)
        {
            this->residual(i) += u[i] * theta;
        }
        u.resize(0);
    }

    // - Get u_gas
    u = this->Get_u_gas(aCurrChemComp);

    // - Add contribution to residual vector
    if(u.size()!=0)
    {
        // - Access the gas phase
        it1 = this->mPhases.equal_range(eGasPhase); 
        double theta = 0;

        if(it1.first!=this->mPhases.end()) 
        {
            for (it=it1.first; it!=it1.second; it++) // OSS: it won't work for more than one mineral phase...same for ads and gases
            {
                // - Get the name of the phase
                phaseName = it1.first->second->name().toStdString();

                // - Retrieve the volumetric fraction from CChemicalComposition
                it3 = aCurrChemComp->mVolumetricFractions.find(phaseName);
                theta = it3->second;
            }
        }

        // - Add contribution to residual vector
        for(int i=0; i!=this->mPrimSpRedNum; i++)
        {
            this->residual(i) += u[i];
        }
        u.resize(0);
    }

    // - Get u_surf 
    u = this->Get_u_surf(aCurrChemComp);

    // - Add contribution to residual vector
    if(u.size()!=0)
    {
        // - Access the surface phase
        it1 = this->mPhases.equal_range(eSurfacePhase); 
        double theta = 0;

        if(it1.first!=this->mPhases.end()) 
        {
            for (it=it1.first; it!=it1.second; it++) // OSS: it won't work for more than one mineral phase...same for ads and gases
            {
                // - Get the name of the phase
                phaseName = it1.first->second->name().toStdString();

                // - Retrieve the volumetric fraction from CChemicalComposition
                it3 = aCurrChemComp->mVolumetricFractions.find(phaseName);
                theta = it3->second;
            }
        }
        // - Add contribution to residual vector
        for(int i=0; i!=this->mPrimSpRedNum; i++)
        {
            this->residual(i) += u[i];
        }
        u.resize(0);
    }
 
    // Evaluate the residual
    for(int i=0; i!=this->mPrimSpRedNum; i++)
    {
        this->residual(i) = cTotRed(i) - this->residual(i);
    }

    c2.resize(0);
    cTotRed.resize(0);
}

void 
CLocalChemicalSystem::ComputeResidualOfMAL(VectorXd &residual, VectorXd &c1, VectorXd &c2, CChemicalComposition* aCurrChemComp)
{
    // Vector of activity coefficients of primary species for the current ChemicalComposition
    VectorXd primaryActCoeffs = VectorXd::Zero(this->mPrimSpRedNum);
    aCurrChemComp->GetSetOfActivityCoeffVector(primaryActCoeffs, this->mRelPrimSpcIndicesWithoutCA);

    // Vector of activity coefficients of secondary species for the current ChemicalComposition
    VectorXd secActCoeffs = VectorXd::Zero(this->mSecSpeciesNum);
    aCurrChemComp->GetSetOfActivityCoeffVector(secActCoeffs, mRelSecondarySpcIndices);

    // Get c1,c
    VectorXd lnc1_c = VectorXd::Zero(this->mConstActSpeciesNum);
    aCurrChemComp->GetSetOfLnConcVector(lnc1_c, mRelConstActSpcIndices);

    // Vector of activity coefficients of c1,c for the current ChemicalComposition
    VectorXd primaryActCoeffs_c = VectorXd::Zero(this->mConstActSpeciesNum);
    aCurrChemComp->GetSetOfActivityCoeffVector(primaryActCoeffs_c, this->mRelConstActSpcIndices);

    // Transform gamma in ln(gamma)
    for(int i=0; i!=this->mPrimSpRedNum; i++)
    {
        primaryActCoeffs(i) = log(primaryActCoeffs(i));
    }
    for(int i=0; i!=this->mConstActSpeciesNum; i++)
    {
        primaryActCoeffs_c(i) = log(primaryActCoeffs_c(i));
    }
    for(int i=0; i!=this->mSecSpeciesNum; i++)
    {
        secActCoeffs(i) = log(secActCoeffs(i));
    }

    // Evaluate residual
    if(this->mSecMod.size()!=0)
    {
        residual = this->mS1Mod * primaryActCoeffs + this->mS1Mod * c1 + this->mSecMod * primaryActCoeffs_c + this->mSecMod * lnc1_c - c2 - secActCoeffs + 
                   this->mLogkMod * 2.302585093;
    }
    else
    {
        residual = this->mS1Mod * primaryActCoeffs + this->mS1Mod * c1 - c2 - secActCoeffs + 
                   this->mLogkMod * 2.302585093;

    }

    primaryActCoeffs.resize(0);
    secActCoeffs.resize(0);

}

int 
CLocalChemicalSystem::CheckConvergence(VectorXd &cnew, VectorXd &deltac, VectorXd &residual)
{
    // Check max relative error in the solution
    double maxRelError = abs((deltac(0))/(cnew(0)));

    for(int i=0; i!=cnew.size(); i++)
    {
        double relError = abs((deltac(i))/(cnew(i)));
        if( relError > maxRelError) 
        {
            maxRelError = relError;
        }
        
    }

    // Check max residual
    double maxResidual = 1.0e+10;
    if (residual.size()!=0)
    {
        maxResidual = abs(residual(0));
        for(int j=1; j!=residual.size(); j++)
        {
            double jthresidual = abs(residual(j));
            if( jthresidual > maxResidual) maxResidual = jthresidual;
        }
    }

    if(maxRelError < mMaxRelativeError && maxResidual < mMaxResidual)
    {
        return 1;
    }

    return 0;
}

int 
CLocalChemicalSystem::CheckConvergence(VectorXd &cnew, VectorXd &deltac, VectorXd &residual, double& relErrMax, int index)
{
    VectorXd lnc = VectorXd::Zero(cnew.size());

    // Transform c in lnc
    for(int i=0; i!=cnew.size(); i++)
    {
        lnc(i) = log(cnew(i));
    }

    // Check max relative error in the solution
    double maxRelError = abs((deltac(0))/(lnc(0)));
    relErrMax = maxRelError;

    for(int i=1; i!=cnew.size(); i++)
    {
        double relError = abs((deltac(i))/(lnc(i)));
        if( relError > maxRelError) 
        {
            maxRelError = relError;
            relErrMax = relError;
        }
        
    }

    // Check max residual
    //double maxResidual = 1.0e+10;
    //if (residual.size()!=0)
    //{
    //    maxResidual = abs(residual(0));
    //    for(int j=1; j!=residual.size(); j++)
    //    {
    //        double jthresidual = abs(residual(j));
    //        if( jthresidual > maxResidual) maxResidual = jthresidual;
    //    }
    //}

    //if(maxRelError < mMaxRelativeError && maxResidual < mMaxResidual)
    //{
    //    return 1;
    //}

    if(maxRelError < mMaxRelativeError) return 1;

    return 0;
}

void 
CLocalChemicalSystem::ComputeJacobian(MatrixXd &compMatrixPrim)
{
    // Set the jacobian to zero
    this->jacobian = MatrixXd::Zero(this->mPrimSpRedNum, this->mPrimSpRedNum);

    // Evaluate the jacobian
    this->jacobian = compMatrixPrim.block(0,0,this->mPrimSpRedNum,this->mPrimSpRedNum) + 
                     compMatrixPrim.block(0,this->mPrimSpRedNum,this->mPrimSpRedNum,this->mSecSpeciesNum) * dc2_dc1;
 
}

void 
CLocalChemicalSystem::ComputeJacobian_allPhases(CChemicalComposition* aCurrChemComp)
{
    // Retrieve the volumetric fraction corresponding to the aqueous phase:

    // - First declare iterators and useful vector
    map<string, double>::iterator it3;
    map<CSpecies*, int>::iterator it4;
    VectorXd thetaVector = VectorXd::Zero(this->mPrimSpRedNum+this->mEqReactions.size());

    // - Access to the aqueous phase
    it1 = this->mPhases.equal_range(eAqueousPhase); 

    // - Get the name of the phase
    string phaseName = it1.first->second->name().toStdString();

    // - Retrieve the volumetric fraction from CChemicalComposition
    it3 = aCurrChemComp->mVolumetricFractions.find(phaseName);
    double theta_liq = it3->second;

    //- Loop over species belonging to this phase
    for(int j=0; j!=it1.first->second->mSpeciesNames.size(); j++)
    {
        //- Loop over species index
        map<CSpecies*, int>::iterator it4;
        for(it4=this->mSpeciesIndices.begin(); it4!=this->mSpeciesIndices.end(); it4++)
        {
            if(it4->first->name().toStdString()==it1.first->second->mSpeciesNames[j])
            {
                thetaVector(it4->second) = theta_liq;
            }
        }
    }

    // - Access the mineral phase
    it1 = this->mPhases.equal_range(eMineralPhase); 
    double theta_min = 0;


    if(it1.first!=this->mPhases.end() && this->mCAS.find(it1.first->second->name().toStdString())==this->mCAS.end()) 
    {
        for (it=it1.first; it!=it1.second; it++) // OSS: it won't work for more than one mineral phase...same for ads and gases
        {
            // - Get the name of the phase
            phaseName = it1.first->second->name().toStdString();

            // - Retrieve the volumetric fraction from CChemicalComposition
            it3 = aCurrChemComp->mVolumetricFractions.find(phaseName);
            theta_min = it3->second;
             
            //- Loop over species belonging to this phase
            for(int j=0; j!=it1.first->second->mSpeciesNames.size(); j++)
            {
                //- Loop over species index
                map<CSpecies*, int>::iterator it4;
                for(it4=this->mSpeciesIndices.begin(); it4!=this->mSpeciesIndices.end(); it4++)
                {
                    if(it4->first->name().toStdString()==it1.first->second->mSpeciesNames[j])
                    {
                        int index = it4->second;
                        thetaVector(index) = theta_min;
                    }
                }
            }
            theta_min = 0;
        }
    }

    // - Access the gas phase
    it1 = this->mPhases.equal_range(eGasPhase); 
    double theta_gas = 0;
    if(it1.first!=this->mPhases.end())
    {
        for (it=it1.first; it!=it1.second; it++)
        {
            // - Get the name of the phase
            phaseName = it1.first->second->name().toStdString();

            // - Retrieve the volumetric fraction from CChemicalComposition
            it3 = aCurrChemComp->mVolumetricFractions.find(phaseName);
            theta_gas = it3->second;

            //- Loop over species belonging to this phase
            for(int j=0; j!=it1.first->second->mSpeciesNames.size(); j++)
            {
                //- Loop over species index
                map<CSpecies*, int>::iterator it4;
                for(it4=this->mSpeciesIndices.begin(); it4!=this->mSpeciesIndices.end(); it4++)
                {
                    if(it4->first->name().toStdString()==it1.first->second->mSpeciesNames[j])
                    {
                        thetaVector(it4->second) = theta_min;
                    }
                }
            }
            theta_gas = 0;
        }
    }

    // - Access to the surface phase
    it1 = this->mPhases.equal_range(eSurfacePhase); 
    double theta_surf = 0;
    if(it1.first!=this->mPhases.end())
    {
        for (it=it1.first; it!=it1.second; it++)
        {
            // - Get the name of the phase
            phaseName = it1.first->second->name().toStdString();

            // - Retrieve the volumetric fraction from CChemicalComposition
            it3 = aCurrChemComp->mVolumetricFractions.find(phaseName);
            theta_surf = it3->second;

            //- Loop over species belonging to this phase
            for(int j=0; j!=it1.first->second->mSpeciesNames.size(); j++)
            {
                //- Loop over species index
                map<CSpecies*, int>::iterator it4;
                for(it4=this->mSpeciesIndices.begin(); it4!=this->mSpeciesIndices.end(); it4++)
                {
                    if(it4->first->name().toStdString()==it1.first->second->mSpeciesNames[j])
                    {
                        thetaVector(it4->second) = theta_min;
                    }
                }
            }
            theta_surf = 0;
        }
    }

    // Set the jacobian to the identity matrix
    this->jacobian = MatrixXd::Zero(this->mPrimSpRedNum, this->mPrimSpRedNum);

    // Evaluate the jacobian
	if(dc2_dc1.size()!=0)
	{
		this->jacobian = this->mComponentMatrix.block(0,0,this->mPrimSpRedNum,this->mPrimSpRedNum) * thetaVector.head(this->mPrimSpRedNum) + 
				   this->mComponentMatrix.block(0,this->mPrimSpRedNum,this->mPrimSpRedNum,this->mSecSpeciesNum) * dc2_dc1 * thetaVector.tail(this->mSecSpeciesNum) ;
	}
	else
	{
		this->jacobian = this->mComponentMatrix.block(0,0,this->mPrimSpRedNum,this->mPrimSpRedNum) * thetaVector.head(this->mPrimSpRedNum);
	}
 
}

vector<int>& 
CLocalChemicalSystem::GetPrimarySpcIndices()
{
	return (this->mRelPrimSpcIndicesWithoutCA);
}

vector<int>& 
CLocalChemicalSystem::GetSecondarySpcIndices()
{
    return (this->mRelSecondarySpcIndices);
}

MatrixXd& 
CLocalChemicalSystem::GetS1Mod()
{
    return (this->mS1Mod);
}

bool 
CLocalChemicalSystem::IsZeroVector(VectorXd &vector)
{
    for(int i=0; i!=vector.size(); i++)
    {
        if(vector(i)!=0.0) return false;
    }
    return true;
}

void 
CLocalChemicalSystem::InitializeSecondaries(VectorXd &c2, VectorXd &c1, CChemicalComposition* currChemComp)
{
    // Evaluate ln(c2) = S1*ln(c1) + log(10) logK*
    c2 = this->mS1Mod * c1 + this->mLogkMod * 2.302585093;
 
}

void 
CLocalChemicalSystem::ApplyConstraints(CChemicalComposition* aCurrChemComp, VectorXd &c1, bool onPrimary, bool ispHAssigned, bool isPGasAssigned)
{
    if(onPrimary && aCurrChemComp->iCon.size()!=0)
    {
        // Loop over reduced primary species
        for(int i=0; i!=this->mPrimSpRedNum; i++)
        {
            // Evaluate concentration from total concentration values
            if(aCurrChemComp->iCon[i] == "cTot")
            {
                this->jacobian.row(i) = VectorXd::Zero(this->mPrimSpRedNum);

                vector<CSpecies*> spc1, spc2;
                map<CSpecies*,int>::iterator it;
                map<string,int>::iterator it1;
                map<string,CSpecies*>::iterator it2;
                int index1, index2, coeff, indexClassic = 0;

                // Get component name and index in the local
                string compName = aCurrChemComp->constraints[i];
                it2 = this->mGlobChemSysPointer->mAllSpecies.find(compName);
                it = this->mSpeciesIndices.find(it2->second);
                index1 = it->second;                   
                it1 = this->mClassicComponents.find(compName);
                if(it1 == this->mClassicComponents.end()) cout << PrTr("CLocalChemicalSystem", ERROR_CLASSICCOMPONENT_NOT_FOUND) <<  qPrintable(spc1[i]->name()) << endl; 
                indexClassic = it1->second;

                // Loop over primary species
                spc1 = this->mSpecies[eAqueous1nc];
                for(int k=0; k!=spc1.size(); k++)
                {
                    // Get coefficient of the species in classic component matrix
                    it1 = this->mClassicComponents.find(spc1[i]->name().toStdString());
                    coeff = this->mComponentMatrixInit(indexClassic,it1->second);

                    // Put the jacobian value
                    this->jacobian(i,index1) += coeff;
                    
                }
                
                // Loop over primary species
                spc2 = this->mSpecies[eAqueous2];

                double sum = 0.0;
                // Loop over secondary species
                for(int j=0; j!=spc2.size(); j++)
                {
                    // Get coefficient of the species in classic component matrix
                    it1 = this->mClassicComponents.find(spc2[j]->name().toStdString());
                    coeff = this->mComponentMatrixInit(indexClassic,it1->second);

                    index2 = this->mSpeciesIndices[spc2[j]] - this->mPrimSpRedNum;
 
                    // Evaluate charge balance
                    sum += coeff * dc2_dc1(index2,index1);

                }
                // Put the jacobian value
                this->jacobian(i,index1) += sum;
                this->residual(i) = aCurrChemComp->constrValue(i) - c1(i);
            }

            // Evaluate concentration from total concentration values
            else if(aCurrChemComp->iCon[i] == "cpri")
            {
                this->residual(i) = 0.0;
                this->jacobian(i,i) = 1.0;
            }
            // Evaluate concentration from charge balance
            else if(aCurrChemComp->iCon[i] == "chgbal")
            {
                // Select the aqueous phase
                it1 = this->mPhases.equal_range(eAqueousPhase);               
                CPhase* acurrPhase = it1.first->second;
                    
                // Set the aqueous phase (i.e. dependent-temperature coefficients)
                acurrPhase->Set(aCurrChemComp);

                // Declare and initialize the charge balance and its derivative
                double chargeBalance = 0.0;
                VectorXd chargeBalanceDeriv = VectorXd::Zero(this->mPrimSpRedNum);

                // Evaluate the charge balance
                acurrPhase->ComputeChargeBalance(aCurrChemComp->mConcentration, chargeBalance, this->mSpecies[eAqueous1nc], 
                                                 this->mSpecies[eAqueous2], this->mSpeciesIndices);

                // Evaluate the derivative of the charge balance
                acurrPhase->ComputeChargeBalanceDerivs(this->dc2_dc1, chargeBalanceDeriv, this->mSpecies[eAqueous1nc], 
                                                       this->mSpecies[eAqueous2], this->mSpeciesIndices);

                // Substitute the charge balance in the residual
                this->residual(i) = - chargeBalance;

                // Substitute the charge balance derivative in the jacobian
                this->jacobian.row(i) = chargeBalanceDeriv;

                chargeBalanceDeriv.resize(0);
            }

            // Evaluate concentration from equilibrium with a mineral
            else if(aCurrChemComp->iCon[i] == "eqmin")
            {
                // Get the name of the mineral
                string minName = aCurrChemComp->constraints[i];

                for(int j=0; j!=this->mEqReactions.size(); j++)
                {
                    if(this->mEqReactions[j]->name().toStdString() == minName)
                    {
                        // Get the value of the mineral SI
                        // double si = aCurrChemComp->mSatIndices[minName];

                        // Substitute the SI in the residual (by now assuming that the mineral is in equilibrium)
                        //this->residual(i) = 1.0 - si(j);
                        this->residual(i) = 1.0;

                        // Substitute the charge balance derivative in the jacobian
                        //this->jacobian.row(i) = dsi(j);
                        this->jacobian.row(i).setZero();
                    }
                }
 
            }
            //// Evaluate concentration from activity values
            //else if(aCurrChemComp->iCon[i] == "activity" && !ispHAssigned)
            //{
            //    // Assign pH value
            //    aCurrChemComp->mConcentration[i] = -log(aCurrChemComp->constrValue(i));

            //    // Eliminate the row relative to H+
            //    MatrixXd auxMatrix = MatrixXd::Zero(this->mComponentMatrixInit.rows()-1,this->mComponentMatrixInit.cols()-1);
            //    for(int j=0; j!=this->mComponentMatrixInit.rows(); j++)
            //    {
            //        for(int k=0; k!=this->mComponentMatrixInit.cols(); k++)
            //        {
            //            if(j==i || k==i) continue;

            //            auxMatrix(j,k) = this->mComponentMatrixInit(j,k);
            //        }
            //    }
            //    this->mComponentMatrixInit.conservativeResize(auxMatrix.rows(),auxMatrix.cols());
            //    this->mComponentMatrixInit = auxMatrix;
            //    auxMatrix.resize(0,0);
            //    ispHAssigned = true;
            //}
            //else if(aCurrChemComp->iCon[i] == "eqgas" && !isPGasAssigned)
            //{
            //    // Assign pressure value
            //    aCurrChemComp->mConcentration[i] = aCurrChemComp->constrValue(i);

            //    // Eliminate the row relative to gas species
            //    MatrixXd auxMatrix = MatrixXd::Zero(this->mComponentMatrixInit.rows()-1,this->mComponentMatrixInit.cols()-1);
            //    for(int j=0; j!=this->mComponentMatrixInit.rows(); j++)
            //    {
            //        for(int k=0; k!=this->mComponentMatrixInit.cols(); k++)
            //        {
            //            if(j==i || k==i) continue;

            //            auxMatrix(j,k) = this->mComponentMatrixInit(j,k);
            //        }
            //    }
            //    this->mComponentMatrixInit.conservativeResize(auxMatrix.rows(),auxMatrix.cols());
            //    this->mComponentMatrixInit = auxMatrix;
            //    auxMatrix.resize(0,0);
            //    isPGasAssigned = true;
            //}
        }
    }
    else   // For primary species that are treated as secondary because the constant activity species eliminated
    {
        if(aCurrChemComp->iCon.size()!=0) 
        {
            // Loop over secondary species that were primary in the input
            for(int i=0; i!=this->mConstActSpeciesNum; i++)
            {
                if(aCurrChemComp->iCon[this->mRelSecondarySpcIndices[i]] == "eqmin" || aCurrChemComp->iCon[this->mRelSecondarySpcIndices[i]] == "eqgas")
                {
                    // Get the name of the mineral/gas
                    string spcName = aCurrChemComp->constraints[this->mRelSecondarySpcIndices[i]];

                    for(int j=0; j!=this->mEqReactions.size(); j++)
                    {
                        if(this->mEqReactions[j]->name().toStdString() == spcName) // Select the  mineral reaction
                        {
                            double concValue = 0;
                            this->mEqReactions[j]->SolveMALExplicit(concValue, this->mRelSecondarySpcIndices[i], aCurrChemComp->mConcentration, aCurrChemComp->mActivityCoeff, aCurrChemComp->mTemp, this->mSpeciesIndices);

                            // Substitute the concentration value in aCurrChemComp, in the concentration vector and in constrValue
                            aCurrChemComp->UpdateConcentration(concValue, this->mRelSecondarySpcIndices[i]);
                            aCurrChemComp->constrValue(this->mRelSecondarySpcIndices[i]) = concValue;
                        }
                    }
                }
                else if(aCurrChemComp->iCon[i] == "activity")
                {
                }
                else if(aCurrChemComp->iCon[i] == "ctot")
                {
                    // do nothing
                }
                else
                {
                    // do nothing
                }
            }
        }

    }
}

MatrixXd 
CLocalChemicalSystem::GetComponentMatrix()
{
    return (this->mComponentMatrix);
}

void 
CLocalChemicalSystem::ComputeSecondaries(CChemicalComposition* aCurrChemComp,vector<double> &c1)
{
	VectorXd cPrimary = VectorXd::Zero(c1.size());
	for(int i=0; i!=c1.size(); i++)
	{
		cPrimary(i) = c1[i];
	}

	this->ComputeSecondaries(aCurrChemComp, cPrimary, true, 0);
}

vector<double> 
CLocalChemicalSystem::Get_u_aq(CChemicalComposition* aChemComp)
{
    vector<double> u_aq;

    vector<CSpecies*> aqueous1 = this->mSpecies[eAqueous1nc];
    vector<CSpecies*> aqueous2 = this->mSpecies[eAqueous2];

    vector<int> indices;

    // First evaluate molarity from concentration
	it1 = this->mPhases.equal_range(eAqueousPhase);
    if(it1.first!=this->mPhases.end())
    {
        for (it=it1.first; it!=it1.second; it++)
        {
            CPhase* aAqPhase = it1.first->second;
		    //aAqPhase->ComputeMolarity(aChemComp->mConcentration,mSpecies[eAqueous1nc], mSpecies[eAqueous2], 
      //                                mSpeciesIndices, aChemComp);
        }
    }

    // Fill the vector with the species indices
    map<CSpecies*, int>::iterator it;
    for(int i=0; i!= aqueous1.size(); i++)
    {
        it = this->mSpeciesIndices.find(aqueous1.at(i));
        indices.push_back(it->second);
    }
    for(int i=0; i!= aqueous2.size(); i++)
    {
        it = this->mSpeciesIndices.find(aqueous2.at(i));
        indices.push_back(it->second);
    }

    // Evaluate u_aq
    for(int i=0; i!=this->mPrimSpRedNum; i++)
    {
        double u_val = 0;
        for(int j=0; j!=indices.size(); j++)
        {
            u_val = u_val + this->mComponentMatrix(i,indices.at(j)) * aChemComp->mConcentration[indices.at(j)];
        }
        u_aq.push_back(u_val);
    }


    return u_aq;
}

vector<double> 
CLocalChemicalSystem::Get_u_aq_classic(CChemicalComposition* aChemComp)
{
    vector<double> u_aq;

    vector<CSpecies*> aqueous1 = this->mSpecies[eAqueous1nc];
    vector<CSpecies*> aqueous2 = this->mSpecies[eAqueous2];

    vector<int> indices;

    // First evaluate molarity from concentration
	it1 = this->mPhases.equal_range(eAqueousPhase);
    if(it1.first!=this->mPhases.end())
    {
        for (it=it1.first; it!=it1.second; it++)
        {
            CPhase* aAqPhase = it1.first->second;
		    //aAqPhase->ComputeMolarity(aChemComp->mConcentration,mSpecies[eAqueous1nc], mSpecies[eAqueous2], 
      //                                mSpeciesIndices, aChemComp);
        }
    }

    // Fill the vector with the species indices
    map<string, int>::iterator it;
    for(int i=0; i!= aqueous1.size(); i++)
    {
        it = this->mClassicComponents.find(aqueous1.at(i)->name().toStdString());
        indices.push_back(it->second);
    }
    for(int i=0; i!= aqueous2.size(); i++)
    {
        it = this->mClassicComponents.find(aqueous2.at(i)->name().toStdString());
        indices.push_back(it->second);
    }

    // Evaluate u_aq
    for(int i=0; i!=this->mComponentMatrixInit.rows(); i++)
    {
        double u_val = 0;
        for(int j=0; j!=indices.size(); j++)
        {
            u_val = u_val + this->mComponentMatrixInit(i,indices.at(j)) * aChemComp->mConcentration[indices.at(j)];
        }
        u_aq.push_back(u_val);
    }

    // If there is water in the system, delete relative component which is always the first (by now we don't calculate water concentration)
    it = this->mClassicComponents.find("h2o");
    if(it!=this->mClassicComponents.end())
    {
        for(int i=0; i!=u_aq.size(); i++)
        {
            u_aq[i] = u_aq[i+1]; // Shift all the components of one
        }
        u_aq.pop_back();  // Delete the last element
    }


    return u_aq;
}


vector<double> 
CLocalChemicalSystem::Get_u_gas(CChemicalComposition* aChemComp)
{
    vector<double> u_gas;

    vector<CSpecies*> gas1 = this->mSpecies[eGas1];
    vector<CSpecies*> gas2 = this->mSpecies[eGas2];

    vector<int> indices;

   // First evaluate molarity from concentration
	it1 = this->mPhases.equal_range(eGasPhase);
    if(it1.first!=this->mPhases.end() &&  this->mCAS.find(it1.first->second->name().toStdString())==this->mCAS.end())
    {
        for (it=it1.first; it!=it1.second; it++)
        {
            CPhase* aGasPhase = it1.first->second;
		    //aGasPhase->ComputeMolarity(aChemComp->mConcentration,mSpecies[eAqueous1nc], mSpecies[eAqueous2], 
      //                                mSpeciesIndices, aChemComp);
        }

        // Fill the vector with the species indices
        map<CSpecies*, int>::iterator it;
        for(int i=0; i!= gas1.size(); i++)
        {
            it = this->mSpeciesIndices.find(gas1.at(i));
            indices.push_back(it->second);
        }
        for(int i=0; i!= gas2.size(); i++)
        {
            it = this->mSpeciesIndices.find(gas2.at(i));
            indices.push_back(it->second);
        }

        // Evaluate u_gas
        for(int i=0; i!=this->mPrimSpRedNum; i++)
        {
            double u_val = 0;
            for(int j=0; j!=indices.size(); j++)
            {
                u_val = u_val + this->mComponentMatrix(i,indices.at(j)) * aChemComp->mConcentration[indices.at(j)];
            }
            u_gas.push_back(u_val);
        }
    }

    return u_gas;

}

vector<double> 
CLocalChemicalSystem::Get_u_gas_classic(CChemicalComposition* aChemComp)
{
    vector<double> u_gas;

    vector<CSpecies*> gas1 = this->mSpecies[eGas1];
    vector<CSpecies*> gas2 = this->mSpecies[eGas2];

    vector<int> indices;

   // First evaluate molarity from concentration
	it1 = this->mPhases.equal_range(eGasPhase);
    if(it1.first!=this->mPhases.end() &&  this->mCAS.find(it1.first->second->name().toStdString())==this->mCAS.end())
    {
        for (it=it1.first; it!=it1.second; it++)
        {
            CPhase* aGasPhase = it1.first->second;
		    //aGasPhase->ComputeMolarity(aChemComp->mConcentration,mSpecies[eAqueous1nc], mSpecies[eAqueous2], 
      //                                mSpeciesIndices, aChemComp);
        }

        // Fill the vector with the species indices
        map<string, int>::iterator it;
        for(int i=0; i!= gas1.size(); i++)
        {
            it = this->mClassicComponents.find(gas1.at(i)->name().toStdString());
            indices.push_back(it->second);
        }
        for(int i=0; i!= gas2.size(); i++)
        {
            it = this->mClassicComponents.find(gas2.at(i)->name().toStdString());
            indices.push_back(it->second);
        }

        // Evaluate u_gas
        for(int i=0; i!=this->mComponentMatrixInit.rows(); i++)
        {
            double u_val = 0;
            for(int j=0; j!=indices.size(); j++)
            {
                u_val = u_val + this->mComponentMatrixInit(i,indices.at(j)) * aChemComp->mConcentration[indices.at(j)];
            }
            u_gas.push_back(u_val);
        }
    }

    return u_gas;

}

vector<double> 
CLocalChemicalSystem::Get_u_min(CChemicalComposition* aChemComp)
{
    vector<double> u_min;

    vector<CSpecies*> solid1 = this->mSpecies[eMineral1];
    vector<CSpecies*> solid2 = this->mSpecies[eMineral2];

    vector<int> indices;

   // First evaluate molarity from concentration
	it1 = this->mPhases.equal_range(eMineralPhase);
    if(it1.first!=this->mPhases.end() &&  this->mCAS.find(it1.first->second->name().toStdString())==this->mCAS.end())
    {
        for (it=it1.first; it!=it1.second; it++)
        {
            CPhase* aMinPhase = it1.first->second;
		    //aMinPhase->ComputeMolarity(aChemComp->mConcentration,mSpecies[eAqueous1nc], mSpecies[eAqueous2], 
      //                                mSpeciesIndices, aChemComp);
        }
        // Fill the vector with the species indices
        map<CSpecies*, int>::iterator it;
        for(int i=0; i!= solid1.size(); i++)
        {
            it = this->mSpeciesIndices.find(solid1.at(i));
            indices.push_back(it->second);
        }
        for(int i=0; i!= solid2.size(); i++)
        {
            it = this->mSpeciesIndices.find(solid2.at(i));
            indices.push_back(it->second);
        }

        // Evaluate u_min
        for(int i=0; i!=this->mPrimSpRedNum; i++)
        {
            double u_val = 0;
            for(int j=0; j!=indices.size(); j++)
            {
                u_val = u_val + this->mComponentMatrix(i,indices.at(j)) * aChemComp->mConcentration[indices.at(j)];
            }
            u_min.push_back(u_val);
        }
    }

    return u_min;

}

vector<double> 
CLocalChemicalSystem::Get_u_min_classic(CChemicalComposition* aChemComp)
{
    vector<double> u_min;

    vector<CSpecies*> solid1 = this->mSpecies[eMineral1];
    vector<CSpecies*> solid2 = this->mSpecies[eMineral2];

    vector<int> indices;

   // First evaluate molarity from concentration
	it1 = this->mPhases.equal_range(eMineralPhase);
    if(it1.first!=this->mPhases.end() &&  this->mCAS.find(it1.first->second->name().toStdString())==this->mCAS.end())
    {
        for (it=it1.first; it!=it1.second; it++)
        {
            CPhase* aMinPhase = it1.first->second;
		    //aMinPhase->ComputeMolarity(aChemComp->mConcentration,mSpecies[eAqueous1nc], mSpecies[eAqueous2], 
      //                                mSpeciesIndices, aChemComp);
        }
        // Fill the vector with the species indices
        map<string, int>::iterator it;
        for(int i=0; i!= solid1.size(); i++)
        {
            it = this->mClassicComponents.find(solid1.at(i)->name().toStdString());
            indices.push_back(it->second);
        }
        for(int i=0; i!= solid2.size(); i++)
        {
            it = this->mClassicComponents.find(solid2.at(i)->name().toStdString());
            indices.push_back(it->second);
        }

        // Evaluate u_min
        for(int i=0; i!=this->mComponentMatrixInit.rows(); i++)
        {
            double u_val = 0;
            for(int j=0; j!=indices.size(); j++)
            {
                u_val = u_val + this->mComponentMatrixInit(i,indices.at(j)) * aChemComp->mConcentration[indices.at(j)];
            }
            u_min.push_back(u_val);
        }
    }

    return u_min;

}

vector<double> 
CLocalChemicalSystem::Get_u_surf(CChemicalComposition* aChemComp)
{
    vector<double> u_surf;

    vector<CSpecies*> solid1 = this->mSpecies[eSurface1];
    vector<CSpecies*> solid2 = this->mSpecies[eSurface2];

    vector<int> indices;

   // First evaluate molarity from concentration
	it1 = this->mPhases.equal_range(eSurfacePhase);
    if(it1.first!=this->mPhases.end())
    {
        for (it=it1.first; it!=it1.second; it++)
        {
            CPhase* aSurfPhase = it1.first->second;
		    //aSurfPhase->ComputeMolarity(aChemComp->mConcentration,mSpecies[eAqueous1nc], mSpecies[eAqueous2], 
      //                                mSpeciesIndices, aChemComp);
        }

        // Fill the vector with the species indices
        map<string, int>::iterator it;
        for(int i=0; i!= solid1.size(); i++)
        {
            it = this->mClassicComponents.find(solid1.at(i)->name().toStdString());
            indices.push_back(it->second);
        }
        for(int i=0; i!= solid2.size(); i++)
        {
            it = this->mClassicComponents.find(solid2.at(i)->name().toStdString());
            indices.push_back(it->second);
        }

        // Evaluate u_surf
        for(int i=0; i!=this->mComponentMatrixInit.rows(); i++)
        {
            double u_val = 0;
            for(int j=0; j!=indices.size(); j++)
            {
                u_val = u_val + this->mComponentMatrixInit(i,indices.at(j)) * aChemComp->mConcentration[indices.at(j)];
            }
            u_surf.push_back(u_val);
        }
    }

    return u_surf;

}

vector<double> 
CLocalChemicalSystem::Get_u_surf_classic(CChemicalComposition* aChemComp)
{
    vector<double> u_surf;

    vector<CSpecies*> solid1 = this->mSpecies[eSurface1];
    vector<CSpecies*> solid2 = this->mSpecies[eSurface2];

    vector<int> indices;

   // First evaluate molarity from concentration
	it1 = this->mPhases.equal_range(eSurfacePhase);
    if(it1.first!=this->mPhases.end())
    {
        for (it=it1.first; it!=it1.second; it++)
        {
            CPhase* aSurfPhase = it1.first->second;
		    //aSurfPhase->ComputeMolarity(aChemComp->mConcentration,mSpecies[eAqueous1nc], mSpecies[eAqueous2], 
      //                                mSpeciesIndices, aChemComp);
        }

        // Fill the vector with the species indices
        map<CSpecies*, int>::iterator it;
        for(int i=0; i!= solid1.size(); i++)
        {
            it = this->mSpeciesIndices.find(solid1.at(i));
            indices.push_back(it->second);
        }
        for(int i=0; i!= solid2.size(); i++)
        {
            it = this->mSpeciesIndices.find(solid2.at(i));
            indices.push_back(it->second);
        }

        // Evaluate u_surf
        for(int i=0; i!=this->mPrimSpRedNum; i++)
        {
            double u_val = 0;
            for(int j=0; j!=indices.size(); j++)
            {
                u_val = u_val + this->mComponentMatrix(i,indices.at(j)) * aChemComp->mConcentration[indices.at(j)];
            }
            u_surf.push_back(u_val);
        }
    }

    return u_surf;

}

MatrixXd 
CLocalChemicalSystem::GetComponentMatrixByPhase(QString aVarType)
{
    MatrixXd result;
       
    vector<CSpecies*> species1;
    vector<CSpecies*> species2;

    vector<int> indices;

    if(aVarType == "component_aq")
	{
		species1 = this->mSpecies[eAqueous1];
        species2 = this->mSpecies[eAqueous2];
	}
    
    else if(aVarType == "component_gas")
	{
		species1 = this->mSpecies[eGas1];
        species2 = this->mSpecies[eGas2];
	}

    else if(aVarType == "component_solid")
	{
        species1 = this->mSpecies[eMineral1];
        species2 = this->mSpecies[eMineral2];
	}

    // Fill the vector with the species indices
    map<CSpecies*, int>::iterator it;
    for(int i=0; i!= species1.size(); i++)
    {
        it = this->mSpeciesIndices.find(species1.at(i));
        indices.push_back(it->second);
    }
    for(int i=0; i!= species2.size(); i++)
    {
        it = this->mSpeciesIndices.find(species2.at(i));
        indices.push_back(it->second);
    }


    result = MatrixXd::Zero(this->mComponentMatrix.rows(), indices.size()); 
                
    for(int i=0; i!=result.rows(); i++)
    {
        for(int j=0; j!=indices.size(); j++)
        {
            result(i,j) = this->mComponentMatrix(i,indices.at(j));
        }
    }
 
    return result;
}

MatrixXd 
CLocalChemicalSystem::Compute_dc_dc1_ByPhase(QString aVarType, CChemicalComposition* aCurrChemComp)
{
    MatrixXd result = MatrixXd::Zero(0,0);

    // Get the vector of primary concentrations
    VectorXd c1 = VectorXd::Zero(this->mPrimSpRedNum);
    aCurrChemComp->GetSetOfConcVector(c1, this->mRelPrimSpcIndicesWithoutCA);

    // First compute dc2/dc1 total, if it's not been computed already
    //if(this->dc2_dc1.size() == 0)
    //{
        this->Compute_dc2_dc1(aCurrChemComp, c1);
    //}

    // Then get dc2/dc1 by phases
    result = this->Get_dc_dc1_ByPhase(aVarType);

    return result;
}

MatrixXd 
CLocalChemicalSystem::Get_dc_dc1_ByPhase(QString aVarType)
{
    MatrixXd result = MatrixXd::Zero(0,0);
    return result;
}

void 
CLocalChemicalSystem::CalculateInitialGuess(CChemicalComposition* aCurrChemComp, VectorXd &c)
{
	// If the user has given a first guess, it's already in aCurrChemComp (it has been read)
	aCurrChemComp->GetSetOfConcVector(c, this->mRelPrimSpcIndicesWithoutCA);

	// If not, I evaluate it
    if(IsZeroVector(c)) 
	{
        cout << "first guess of c1 is zero!";
	}
}

void 
CLocalChemicalSystem::MixChemicalCompositions(map<CChemicalComposition*,double> aInitialWaters, CChemicalComposition* aMixedWater, bool checkSI)
{
    // By now this method is equal to the method of mixWaters

    // Add the mixed water to the chemical composition associated to this local chemical system
    this->mRelChemicalCompositions.push_back(aMixedWater);

	// Declare useful iterators
	map<CChemicalComposition*, double>::iterator itw;
	map<string, CChemicalComposition*>::iterator itw1;   

    double prevMixRatio = 0;
    double denom = 0;
    vector<double> u_aq;

    // First evaluate the sum of the volumes of the initial waters
	for(itw = aInitialWaters.begin(); itw != aInitialWaters.end(); itw++)
	{
        denom += itw->second;
    }

	// Loop over initial waters
	for(itw = aInitialWaters.begin(); itw != aInitialWaters.end(); itw++)
	{
        // Copy the water with smallest mixing ratio in the mixedWater, for initial guess values
        if(itw == aInitialWaters.begin())
        {
            aMixedWater->constrValue = VectorXd::Zero(itw->first->constrValue.size());
            aMixedWater->mLocChemSysPointer = this;
            aMixedWater->Copy(itw->first);

            prevMixRatio = itw->second;
        }

        else if(itw != aInitialWaters.begin() && itw->second < prevMixRatio)
        {
            prevMixRatio = itw->second;
            aMixedWater->Copy(itw->first);
        }

		// Retrieve name of the water
		string waterName = itw->first->name().toStdString();

		// Check if this water already belongs to this LocalChemicalSystem
		itw1 = this->mChemicalCompositions.find(waterName);

        // Evaluate u*_aq

        // - Get the aqueous components vector
        u_aq = this->Get_u_aq(itw->first);

        // - Loop over components
	    for(int i=0; i!= this->mPrimSpRedNum; i++)
	    {
		    // - Loop over aqueous concentrations (change loop over constrValue.size())
		    for(int j=0; j!=u_aq.size(); j++)
		    {
			    aMixedWater->constrValue(i) += itw->second * u_aq[j] / denom;
		    }
	    }
        u_aq.resize(0);
    }

	// Speciate 
    this->Speciate_utot(aMixedWater);

    // Evaluate SI of minerals to check if some of them precipitated
	if(checkSI) this->mGlobChemSysPointer->ComputeSI(aMixedWater, this->mSpeciesIndices, this->mCAS);
}

void 
CLocalChemicalSystem::MixWaters(map<CChemicalComposition*,double> &aInitialWaters, CChemicalComposition* aMixedWater, bool checkSI)
{
    // Add the mixed water to the chemical composition associated to this local chemical system
    this->mRelChemicalCompositions.push_back(aMixedWater);

	// Declare useful iterators
	map<CChemicalComposition*, double>::iterator itw;
	map<string, CChemicalComposition*>::iterator itw1;   

    double prevMixRatio = 0;
    double denom = 0;
    vector<double> u_aq;

    // First evaluate the sum of the volumes of the initial waters
	for(itw = aInitialWaters.begin(); itw != aInitialWaters.end(); itw++)
	{
        denom += itw->second;
    }

	// Loop over initial waters
	for(itw = aInitialWaters.begin(); itw != aInitialWaters.end(); itw++)
	{
        // Copy the water with smallest mixing ratio in the mixedWater, for initial guess values
        if(itw == aInitialWaters.begin())
        {
            aMixedWater->constrValue = VectorXd::Zero(itw->first->constrValue.size());
            aMixedWater->mLocChemSysPointer = this;
            aMixedWater->Copy(itw->first);

            prevMixRatio = itw->second;
        }

        else if(itw != aInitialWaters.begin() && itw->second < prevMixRatio)
        {
            prevMixRatio = itw->second;
            aMixedWater->Copy(itw->first);
        }

		// Retrieve name of the water
		string waterName = itw->first->name().toStdString();

		// Check if this water already belongs to this LocalChemicalSystem
		itw1 = this->mChemicalCompositions.find(waterName);

        // Evaluate u*_aq

        // - Get the aqueous components vector
        u_aq = this->Get_u_aq(itw->first);

        // - Loop over components
	    for(int i=0; i!= this->mPrimSpRedNum; i++)
	    {
		    // - Loop over aqueous concentrations (change loop over constrValue.size())
		    for(int j=0; j!=u_aq.size(); j++)
		    {
			    aMixedWater->constrValue(i) += itw->second * u_aq[j] / denom;
		    }
	    }
        u_aq.resize(0);
    }

	// Speciate 
  	this->Speciate_uaq(aMixedWater);

    // Evaluate SI of minerals to check if some of them precipitated
	if(checkSI) this->mGlobChemSysPointer->ComputeSI(aMixedWater, this->mSpeciesIndices, this->mCAS);
}

vector<double>
CLocalChemicalSystem::GetUSkrk(CChemicalComposition* aChemComp)
{
    vector<double> result_vector;

    VectorXd result = VectorXd::Zero(this->mPrimSpRedNum);
    if (this->mKinReactions.size())
    {

    VectorXd rk = VectorXd::Zero(this->mKinReactions.size());

    double omega; // Think about where to evaluate it!
    double area = 0;  // Think about where to evaluate it!
    double rate;


    for(int i=0; i!=this->mKinReactions.size(); i++)
    {
        this->mKinReactions[i]->mReacRateLawPointer->ComputeReactionRate(aChemComp, this->mSpeciesIndices, omega, rate, area);
        rk(i) = rate;
    }

    result = this->mComponentMatrix * this->mSk.transpose() * rk;
    }

    result_vector.resize(result.size());

    for(int j=0; j!=result.size(); j++)
    {
        result_vector[j] = result(j);
    }

    return result_vector;
}

MatrixXd 
CLocalChemicalSystem::GetdUSkrkdt(const double dt, CChemicalComposition* aChemComp)
{
    MatrixXd result = MatrixXd::Zero(this->mPrimSpRedNum, this->mPrimSpRedNum);

    return result;

}

void 
CLocalChemicalSystem::RISA(CChemicalComposition* aCurrChemComp, bool writeHeader)
{
    // Declare matrices Xi and Lambda and fill them
    MatrixXd aXi;
    MatrixXd aLambda;    
    VectorXd aX;
    VectorXd aVinv_diag;

    this->SetRISAInputs(aCurrChemComp, aXi, aLambda, aX, aVinv_diag);

	// Declare vectors and matrices
    VectorXd c1 = VectorXd::Zero(this->mPrimSpRedNum);
    VectorXd deltalnc1 = VectorXd::Zero(this->mPrimSpRedNum);  
    VectorXd fd = VectorXd::Zero(aXi.rows()+aLambda.rows());
    MatrixXd jacobian = MatrixXd::Zero(aXi.rows()+aLambda.rows(),this->mPrimSpRedNum);

    // Vectors to be printed
    double objfun,  norm, relErr = 0;
    VectorXd grad = VectorXd::Zero(aXi.rows()+aLambda.rows());
    QString speciesMaxRelErr;

    // Write header on output files
    if(writeHeader)
    {
        this->WriteChemCompInfos_conc(aCurrChemComp, writeHeader, true);
        this->WriteChemCompInfos_gamma(aCurrChemComp, writeHeader, true);
    }
                    
    // Get the initial guess of primary species (without constant activity species) concentrations from aCurrChemComp
	this->CalculateInitialGuess(aCurrChemComp, c1);

	// Convert molality (mol/kgw) to molarity (mol/m3l) for aqueous phase
	it1 = this->mPhases.equal_range(eAqueousPhase);          
    CPhase* aAqPhase = it1.first->second;
	//aAqPhase->ComputeMolarity(aCurrChemComp->mConcentration,mSpecies[eAqueous1nc], mSpecies[eAqueous2], 
 //                             mSpeciesIndices, aCurrChemComp);

    // Main N-R loop over global iterations
    for(int iter = 0, itHasConvergedToc1 = false; iter < mMaxIterNum && !itHasConvergedToc1; iter++ )
    {
        // Evaluate the concentrations of aqueous secondary species at iter+1 by means of N-R
        this->ComputeSecondaries(aCurrChemComp, c1, true, iter);

        // Compute fd
        this->Computefd(aCurrChemComp, c1, aX, aXi, aLambda, fd);

        // Evaluate the Jacobian
        this->ComputeJacobian(aCurrChemComp, c1, aXi, aLambda, jacobian);

        // Evaluate objective function
        objfun = fd.transpose()*aVinv_diag.asDiagonal()*fd;

        // Compute RHS and LHS
        fd = -jacobian.transpose()*aVinv_diag.asDiagonal()*fd;
        jacobian = jacobian.transpose()*aVinv_diag.asDiagonal()*jacobian;

        // Evaluate gradient and norm of gradient
        grad = -2*fd;
        norm = grad.transpose()*grad;

        // Solve the system (evaluate deltac1)
        deltalnc1 = jacobian.fullPivLu().solve(fd);

        // Update the solution
        for(int i=0; i!= c1.size(); i++)
        {
		    c1(i) = exp(log(c1(i)) + deltalnc1(i));
        }

        // Check convergence
        int index = 0;
        itHasConvergedToc1 = this->CheckConvergence(c1, deltalnc1, fd, relErr, index);  

        // Find species with maximum relative error
        map<CSpecies*, int>::iterator it2;
        for(it2 = this->mSpeciesIndices.begin(); it2!=this->mSpeciesIndices.end(); it2++)
        {
            if(it2->second == index) speciesMaxRelErr = it2->first->name();
        }

        // Write RISA infos in the output file
        this->WriteRISAInfos(aCurrChemComp, iter, objfun, grad, norm, relErr, speciesMaxRelErr);
    }

    // Update the concentrations of the current chemical composition if it has converged
    aCurrChemComp->UpdateConcentrations(c1, this->mRelPrimSpcIndicesWithoutCA);

    // Assign water molality to all chemical composition of the CLocalChemicalSystem (by now here, change it later)
    if(this->mSpecies[eAqueous1c].size() != 0) 
    {
        double waterMolality = 1/H2O_mol_mass;
        aCurrChemComp->UpdateConcentration(waterMolality, mSpeciesIndices[this->mSpecies[eAqueous1c].at(0)]);
    }

    // Print chemical informations on an output file
    this->WriteChemCompInfos(aCurrChemComp, true);

    // Write concentrations on output files
    this->WriteChemCompInfos_conc(aCurrChemComp, false, true);

    // Write concentrations on output files
    this->WriteChemCompInfos_gamma(aCurrChemComp, false, true);

    c1.resize(0);
    deltalnc1.resize(0);
    fd.resize(0);
    jacobian.resize(0,0);
    grad.resize(0);
}

void 
CLocalChemicalSystem::SetRISAInputs(CChemicalComposition* aChemComp, MatrixXd& aXi, MatrixXd& aLambda, VectorXd& aX, VectorXd& aV_diag)
{
    int aIndex = 0;
    vector<CSpecies*> speciesVec;
    map<CSpecies*, int>::iterator it;
	map<string,int>::iterator it1;
	map<string,CSpecies*>::iterator it2;
    vector<double> aV;
    vector<double> aData;
    vector<bool> tempLogXi;
    vector<bool> tempLogLambda;
    int aSpcNum = this->mPrimSpRedNum + this->mSecSpeciesNum;

    // Get primary aq species
    int k=0;
    int z=0;
    double std2;

    // Loop over iCon_uncertain
    for(int j=0; j!=aChemComp->iCon_uncertain.size(); j++)
    {
        // Flag: if iCon == “cTot”
        if(aChemComp->iCon_uncertain.at(j)=="cTot")
        {
            k++;

            // Resize the matrix aXi and initialize to zero relative row
	        aXi.conservativeResize(k,aSpcNum);
            aXi.row(k-1) = VectorXd::Zero(aSpcNum);

            tempLogXi.resize(k);
            tempLogXi.at(k-1) = aChemComp->isLogData.at(j);

            // Get species relative to the condition and its index in the classic component definition
            string spcName = aChemComp->constraints_uncertain.at(j);
			it1 = this->mClassicComponents.find(spcName);

			if(it1 == this->mClassicComponents.end()) cout << PrTr("CLocalChemicalSystem", ERROR_CLASSICCOMPONENT_NOT_FOUND) <<  qPrintable(QString::fromStdString(spcName)) << endl; 

			int indexClassic = it1->second;

			double coeff = 0;

			string nameClassic;
			for(it1=this->mClassicComponents.begin(); it1!=this->mClassicComponents.end(); it1++)
			{
				// Get spc coeff in the classic definition
				nameClassic = it1->first;
				if(nameClassic == "h2o") continue; // Skip water
				it2 = this->mGlobChemSysPointer->mAllSpecies.find(nameClassic);
				it = this->mSpeciesIndices.find(it2->second);
				
				if(it!=this->mSpeciesIndices.end()) // If this species is defined in this local chemical system
				{
					// Put the row j of the component matrix in the row k of aXi	
					aXi(k-1,it->second) = this->mComponentMatrixInit(indexClassic, it1->second);
				}

			}

            // Update aX (data vector) and aV_diag
            if(aChemComp->isLogData.at(j) == false)
            {
                aData.push_back(aChemComp->constrValue_uncertain(j)); // Store the value as it is 

                if(aChemComp->isRelError.at(j) == true)
                {
                    std2 = 1/((aChemComp->errors.at(j) * aChemComp->constrValue_uncertain(j))*(aChemComp->errors.at(j) * aChemComp->constrValue_uncertain(j)));
                }
                else
                {
                    std2 = 1/((aChemComp->errors.at(j) *aChemComp->errors.at(j)));
                }
            }
            else
            {
                aData.push_back(aChemComp->constrValue_uncertain(j)/0.432969291); // If is logaritmic, then transform log10 to ln

                if(aChemComp->isRelError.at(j) == true)
                {
                    std2 = 1/((aChemComp->errors.at(j) * aChemComp->constrValue_uncertain(j)/0.432969291)*(aChemComp->errors.at(j) * aChemComp->constrValue_uncertain(j)/0.432969291));
                }
                else
                {
                    std2 = 1/((aChemComp->errors.at(j) *aChemComp->errors.at(j)));
                }
            }
            aV.push_back(std2);
            std2 = 0;
        }

        // Flag: if iCon == "alk" 
        else if(aChemComp->iCon_uncertain.at(j)=="alk")
        {
            k++;

            // Resize the matrix aXi and initialize to zero relative row
	        aXi.conservativeResize(k,aSpcNum);
            aXi.row(k-1) = VectorXd::Zero(aSpcNum);

            tempLogXi.resize(k);
            tempLogXi.at(k-1) = aChemComp->isLogData.at(j);

            for(it=this->mSpeciesIndices.begin(); it!=this->mSpeciesIndices.end(); it++)
            {
                if(it->first->name()=="hco3-") aXi(k-1,it->second) = 1.0;
                if(it->first->name()=="co3-2") aXi(k-1,it->second) = 2.0;
                if(it->first->name()=="oh-") aXi(k-1,it->second) = 1.0;
                if(it->first->name()=="h+") aXi(k-1,it->second) = -1.0;
            }

            // Update aX (data vector) and aV_diag
            if(aChemComp->isLogData.at(j) == false)
            {
                aData.push_back(aChemComp->constrValue_uncertain(j)); // Store the value as it is

                if(aChemComp->isRelError.at(j) == true)
                {
                    // Restore when RISA simulations are finished!!!
                    /*std2 = 1/((aChemComp->errors.at(j) * aChemComp->constrValue_uncertain(j))*(aChemComp->errors.at(j) * aChemComp->constrValue_uncertain(j)));*/

                    // Delete when RISA simulations are finished!!!
                    std2 = 1/((aChemComp->errors.at(j) * aChemComp->measured_value)*(aChemComp->errors.at(j) * aChemComp->measured_value));
                }
                else
                {
                    std2 = 1/((aChemComp->errors.at(j) *aChemComp->errors.at(j)));
                }
            }
            else
            {
                aData.push_back(aChemComp->constrValue_uncertain(j)/0.4343); // If is logaritmic, then transform log10 to ln

                if(aChemComp->isRelError.at(j) == true)
                {
                    std2 = 1/((aChemComp->errors.at(j) * aChemComp->constrValue_uncertain(j)/0.4343)*(aChemComp->errors.at(j) * aChemComp->constrValue_uncertain(j)/0.4343));
                }
                else
                {
                    std2 = 1/((aChemComp->errors.at(j) *aChemComp->errors.at(j)));
                }
            }
            aV.push_back(std2);
            std2 = 0;
        }

        // Flag: if iCon == “chbal” o “ec”
        else if(aChemComp->iCon_uncertain.at(j)=="chgbal" || aChemComp->iCon_uncertain.at(j)=="ec")
        {
            k++;

            // Resize the matrix aXi and initialize to zero relative row
	        aXi.conservativeResize(k,aSpcNum);
            aXi.row(k-1) = VectorXd::Zero(aSpcNum);

            tempLogXi.resize(k);
            tempLogXi.at(k-1) = aChemComp->isLogData.at(j);

            vector<CSpecies*> aSpcVector;
            double value = 0;
            double sum = 0;

            // Loop over primary aqueous species
            aSpcVector = this->mSpecies[eAqueous1nc];
            for(int m=0; m!=aSpcVector.size(); m++)
            {
                // Get species index
                aIndex = this->mSpeciesIndices[aSpcVector.at(m)];

                // Get species attribute: charge for "chbal" and limiting conductivity for "ec"
                if(aChemComp->iCon_uncertain.at(j)=="chgbal")
                {
                    value = aSpcVector.at(m)->mIonicCharge;
                    sum += value * aChemComp->mConcentration.at(aIndex) * value * aChemComp->mConcentration.at(aIndex);
                }

                // Put the value of the attribute in the right position
                if(aChemComp->iCon_uncertain.at(j)=="ec")
                {
                    value = aSpcVector.at(m)->mLimMolConduct;
                }

                // Put the value in the right position of aXi
                aXi(k-1,aIndex) = value;

            }
            aSpcVector.clear();

            // Loop over secondary aqueous species
            aSpcVector = this->mSpecies[eAqueous2];
            for(int n=0; n!=aSpcVector.size(); n++)
            {
                // Get species index
                aIndex = this->mSpeciesIndices[aSpcVector.at(n)];

                // Get species attribute: charge for "chbal" and limiting conductivity for "ec"
                if(aChemComp->iCon_uncertain.at(j)=="chgbal")
                {
                    value = aSpcVector.at(n)->mIonicCharge;
                    sum += value * aChemComp->mConcentration.at(aIndex) * value * aChemComp->mConcentration.at(aIndex);
                }

                // Put the value of the attribute in the right position
                if(aChemComp->iCon_uncertain.at(j)=="ec")
                {
                    value = aSpcVector.at(n)->mLimMolConduct;
                }

                // Put the value in the right position of aXi
                aXi(k-1,aIndex) = value;

            }
            aSpcVector.clear();
            sum = sqrt(sum);

            // Update aX (data vector) and aV_diag
            aData.push_back(0.0);
            if(aChemComp->iCon_uncertain.at(j)=="chgbal")
            {
                if(aChemComp->isRelError.at(j) == true)
                {
                    std2 = 1/(aChemComp->errors.at(j)* sum * aChemComp->errors.at(j)* sum);
                }
                else
                {
                    std2 = 1/((aChemComp->errors.at(j) *aChemComp->errors.at(j)));
                }
            }
            else // for EC values
            {
                // Update aX (data vector) and aV_diag
                if(aChemComp->isLogData.at(j) == false)
                {
                    aData.push_back(aChemComp->constrValue_uncertain(j)); // Store the value as it is
                    
                    if(aChemComp->isRelError.at(j) == true)
                    {
                        std2 = 1/((aChemComp->errors.at(j) * aChemComp->constrValue_uncertain(j))*(aChemComp->errors.at(j) * aChemComp->constrValue_uncertain(j)));
                    }
                    else
                    {
                        std2 = 1/((aChemComp->errors.at(j) *aChemComp->errors.at(j)));
                    }
                }
                else
                {
                    aData.push_back(aChemComp->constrValue_uncertain(j)/0.4343); // If is logaritmic, then transform log10 to ln
                    
                    if(aChemComp->isRelError.at(j) == true)
                    {
                        std2 = 1/((aChemComp->errors.at(j) * aChemComp->constrValue_uncertain(j)/0.4343)*(aChemComp->errors.at(j) * aChemComp->constrValue_uncertain(j)/0.4343));
                    }
                    else
                    {
                        std2 = 1/((aChemComp->errors.at(j) *aChemComp->errors.at(j)));
                    }
                }
            }
            aV.push_back(std2);
            std2 = 0;

         }
        // Flag: if iCon == “activity”
        else if(aChemComp->iCon_uncertain.at(j)=="activity")
        {
            z++;
	        aLambda.conservativeResize(z,aSpcNum);
            aLambda.row(z-1) = VectorXd::Zero(aSpcNum);

            tempLogLambda.resize(z);
            tempLogLambda.at(z-1) = aChemComp->isLogData.at(j);

			// Get species indices
			string spcName = aChemComp->constraints_uncertain.at(j);
			CSpecies* spcPointer = this->mGlobChemSysPointer->mAllSpecies[spcName];
			int index = this->mSpeciesIndices[spcPointer];

            aLambda(z-1,index)=1.0;

            // Update aX (data vector) and aV_diag
			aData.push_back(aChemComp->constrValue_uncertain(j));

            std2 = 1/((aChemComp->errors.at(j))*(aChemComp->errors.at(j)));
      
            aV.push_back(std2);
            std2 = 0;
        }
        // Flag: if iCon == “eqmin”
        else if(aChemComp->iCon_uncertain.at(j)=="eqmin" || aChemComp->iCon_uncertain.at(j)=="eqgas")
        {
            z++;
	        aLambda.conservativeResize(z,aSpcNum);
            aLambda.row(z-1) = VectorXd::Zero(aSpcNum);

            tempLogLambda.resize(z);
            tempLogLambda.at(z-1) = aChemComp->isLogData.at(j);

 			// Get species indices
			string spcName = aChemComp->constraints_uncertain.at(j);
			CSpecies* spcPointer = this->mGlobChemSysPointer->mAllSpecies[spcName];
			int index = this->mSpeciesIndices[spcPointer];

            aLambda(z-1,index)=1.0;

            if(aChemComp->isLogData.at(j) == false)
            {
                if(aChemComp->isRelError.at(j) == true)
                {
                    std2 = 1/((aChemComp->errors.at(j) * aChemComp->constrValue_uncertain(j))*(aChemComp->errors.at(j) * aChemComp->constrValue_uncertain(j)));
                }
                else
                {
                    std2 = 1/((aChemComp->errors.at(j) *aChemComp->errors.at(j)));
                }
            }
            else
            {
                if(aChemComp->isRelError.at(j) == true)
                {
                    std2 = 1/((aChemComp->errors.at(j) * aChemComp->constrValue_uncertain(j)/0.4343)*(aChemComp->errors.at(j) * aChemComp->constrValue_uncertain(j)/0.4343));
                }
                else
                {
                    std2 = 1/((aChemComp->errors.at(j) *aChemComp->errors.at(j)));
                }
            }

            // Update aX (data vector) and aV_diag
            //aData.push_back(data);
            aData.push_back(aChemComp->constrValue_uncertain(j));

            aV.push_back(std2);

            std2 = 0;
        }
        // Flag: if iCon == “userComp”
        


    }
    aV_diag.resize(aV.size());
    aX.resize(aData.size());
    for(int i=0; i!=aV.size(); i++)
    {
        aX(i) = aData.at(i);
        aV_diag(i) = aV.at(i);
    }

    if(aXi.size()==0 && aLambda.size()==0) this->PrThrow(ERROR_RISA_INITIALIZATION);

    // Update isLogData vector
    aChemComp->isLogData.resize(0);
    aChemComp->isLogData.resize(aXi.rows()+aLambda.rows());
    for(int k=0; k!=aXi.rows(); k++)
    {
        aChemComp->isLogData.at(k) = tempLogXi.at(k);
    }
    for(int j=0; j!=aLambda.rows(); j++)
    {
        aChemComp->isLogData.at(j+aXi.rows()) = tempLogLambda.at(j);
    }
    tempLogXi.resize(0);
    tempLogLambda.resize(0);
}

void 
CLocalChemicalSystem::Computefd(CChemicalComposition* aCurrChemComp, VectorXd &c1, VectorXd &x, MatrixXd& aXi, MatrixXd& aLambda, VectorXd &fd)
{    
    // Get secondary concentrations of the current chemical composition
    VectorXd c2 = VectorXd::Zero(this->mSecSpeciesNum);
    aCurrChemComp->GetSetOfConcVector(c2,mRelSecondarySpcIndices);

	// Get vector of activity coefficients of primary species for the current ChemicalComposition
    VectorXd primaryActCoeffs = VectorXd::Zero(this->mPrimSpRedNum);
    aCurrChemComp->GetSetOfActivityCoeffVector(primaryActCoeffs, this->mRelPrimSpcIndicesWithoutCA);

    // Get vector of activity coefficients of secondary species for the current ChemicalComposition
    VectorXd secActCoeffs = VectorXd::Zero(this->mSecSpeciesNum);
    aCurrChemComp->GetSetOfActivityCoeffVector(secActCoeffs, mRelSecondarySpcIndices);

	fd = VectorXd::Zero(aXi.rows()+aLambda.rows());

    // Fill the first part relative to Xi

    if(aXi.size()!=0) 
    {
		for(int i=0; i!=aXi.rows(); i++)
		{
            if(aCurrChemComp->isLogData.at(i) == true) // is cTot and logaritmic data is given
            {
                /*if(aCurrChemComp->iCon_uncertain.at(i) == "chgbal")
                {
                    throw GeneralException(ERROR_LOGDATA_NOT_ALLOWED);
                }
			    else
                {*/    
                    fd(i) = aXi.row(i).head(this->mPrimSpRedNum)* c1;
			        fd(i) += aXi.row(i).tail(this->mSecSpeciesNum)* c2; 
                    fd(i) = log(fd(i));
			        fd(i) -= x(i);	
                //}
            }
            else // not logaritmic data
            {
			    fd(i) = aXi.row(i).head(this->mPrimSpRedNum) * c1;
			    fd(i) += aXi.row(i).tail(this->mSecSpeciesNum) * c2; 
			    fd(i) -= x(i);				 
            }
		}
    }

    // Fill the second part relative to Lambda
    if(aLambda.size()!=0) 
    {
		VectorXd lnc1 = VectorXd::Zero(this->mPrimSpRedNum);
		// Transform gamma in ln(gamma)
		for(int i=0; i!=this->mPrimSpRedNum; i++)
		{
			lnc1(i) = log(c1(i));
			primaryActCoeffs(i) = log(primaryActCoeffs(i));
		}
		for(int i=0; i!=this->mSecSpeciesNum; i++)
		{
			c2(i) = log(c2(i));
			secActCoeffs(i) = log(secActCoeffs(i));
		}

		for(int k=0; k!=aLambda.rows(); k++)
		{
			fd(k+aXi.rows()) = aLambda.row(k).head(this->mPrimSpRedNum) * (lnc1 + primaryActCoeffs);
			fd(k+aXi.rows()) += aLambda.row(k).tail(this->mSecSpeciesNum) * (c2 + secActCoeffs);
			fd(k+aXi.rows()) -=	log(x(k+aXi.rows()));
		}
		lnc1.resize(0);
    }


    c2.resize(0);
	

}

void
CLocalChemicalSystem::ComputeJacobian(CChemicalComposition* aCurrChemComp, VectorXd &c1, MatrixXd& aXi, MatrixXd& aLambda, MatrixXd& jacobian)
{
    // Get secondary concentrations of the current chemical composition
    VectorXd c2 = VectorXd::Zero(this->mSecSpeciesNum);
    aCurrChemComp->GetSetOfConcVector(c2,mRelSecondarySpcIndices);

	jacobian = MatrixXd::Zero(aXi.rows()+aLambda.rows(),this->mPrimSpRedNum);

    // Fill jacobian part from Xi
    if(aXi.size()!=0)
    {
        jacobian.block(0,0,aXi.rows(),this->mPrimSpRedNum) = MatrixXd::Zero(aXi.rows(),this->mPrimSpRedNum);

        // First evaluates the jacobian as there was no log data (to take advantage of the block operation)
        jacobian.block(0,0,aXi.rows(),this->mPrimSpRedNum) = aXi.block(0,0,aXi.rows(),this->mPrimSpRedNum) + aXi.block(0,this->mPrimSpRedNum,aXi.rows(),this->mSecSpeciesNum) * dc2_dc1;

        // Multiply times c1,j
        for(int i=0; i!=aXi.rows(); i++)
        {
            for(int j=0; j!= jacobian.cols(); j++)
            {
                jacobian(i,j) = jacobian(i,j) * c1(j);
            }
        }

        // Check if some of the data is given in log scale (BY NOW: ONLY ONE CONDITION ALLOWED...TO TRY CTOT!)
        bool someIsLog = false;
        int logRow = 0;
        for(int k=0; k!=aXi.rows(); k++)
        {
            if(aCurrChemComp->isLogData.at(k)==true)
            {
                logRow = k;
                double sum = aXi.row(logRow).head(this->mPrimSpRedNum) * c1;
			    sum += aXi.row(logRow).tail(this->mSecSpeciesNum) * c2; 
                jacobian.row(logRow) = jacobian.row(logRow) / sum;
            }
        }
    }

    // Fill jacobian part from Lambda
    if(aLambda.size()!=0)
    {
        // - Calculate dlnGamma/dlnc1 and put it in jacobian_block
        MatrixXd jacobian_block = MatrixXd::Zero(aLambda.rows(),this->mPrimSpRedNum);

        // - Select the aqueous phase
        it1 = this->mPhases.equal_range(eAqueousPhase);               
        CPhase* acurrPhase = it1.first->second;

        // - Compute dfd2/dlnc1 and add their contribution to the jacobian block
        acurrPhase->ComputeActivityCoeffDerivs(aCurrChemComp->mConcentration, jacobian_block, aLambda, 
                                               mSpecies[eAqueous1nc], mSpecies[eAqueous2], 
                                               mSpeciesIndices, false, true);

        // - Transform dc2/dc1 in dlnc2/dlnc1
        for(int i=0; i!=this->mSecSpeciesNum; i++)
        {
            for(int j=0; j!=this->mPrimSpRedNum; j++)
            {
                dc2_dc1(i,j) = dc2_dc1(i,j) * c1(j) / c2(i);
            }
        }

        // - Add contributions of jacobian_block and derivatives dlnc2/dlnc1
        jacobian.block(aXi.rows(),0,aLambda.rows(),this->mPrimSpRedNum) = jacobian_block + 
                                               aLambda.block(0,this->mPrimSpRedNum,aLambda.rows(),this->mSecSpeciesNum) * dc2_dc1;
    }


}

void
CLocalChemicalSystem::ComputeS2(MatrixXd s2)
{
    // Loop over secondary species indices (the species must be in order of index in S2)
    for(int i=0; i!=mRelSecondarySpcIndices.size(); i++)
    {
        // Get the species correspondent to the index
        for(map<CSpecies*, int>::iterator it= this->mSpeciesIndices.begin() ; it != this->mSpeciesIndices.end() ; it++)
        {
            if(mSpeciesIndices[it->first] == mRelSecondarySpcIndices[i])
            {
                for ( int j = 0 ; j!=mEqReactions.size() ; j++ ) 
	            {
                    // Obtain pointer to current reaction
                    CReaction* aCurrReaction = mEqReactions.at(j);
			        
                    // Obtain stoichiometric coefficient for the current species in the current reaction
			        double stoichCoeff = aCurrReaction->mRelSpeciesVsStoichCoeffs[it->first];

			        // Put it in the matrix
			        this->mS2Inverse(j,i) = stoichCoeff;

                }
             }
         }
     } 
}

void 
CLocalChemicalSystem::SetConvParams(CLocalChemicalSystem* aLocal)
{
    this->mSpeciateWithNR = aLocal->mSpeciateWithNR;
    this->mMaxRelativeError = aLocal->mMaxRelativeError;
    this->mMaxResidual = aLocal->mMaxResidual;
    this->mMaxIterNum = aLocal->mMaxIterNum;

}

void 
CLocalChemicalSystem::WriteChemCompInfos(CChemicalComposition* aCurrChemComp, bool chemCompOnSeparateFiles)
{
    QString fileName; 

    if(chemCompOnSeparateFiles)
    {
        fileName = "output_" + aCurrChemComp->name() + ".out";
    }
    else
    {
        fileName = "output.out";
    }

    fileName = QDir(fileName).absolutePath();

    QFile* outputFile = new QFile(fileName);

    // Open the file, if is not open
    if(!outputFile->isOpen())
    {
        if ( !outputFile->open(QFile::WriteOnly | QIODevice::Append))
        {
            this->PrThrow(ERROR_FILEOPEN);
        }
    }

    QTextStream* output_stream = new QTextStream(outputFile);

    // Write name of chemical composition
    *output_stream << endl << "-----------------------------------------------------------------------------------------------------" ;
    *output_stream << endl << "CHEMICAL COMPOSITION:" << "\t" << aCurrChemComp->name() << endl;
    *output_stream << "-----------------------------------------------------------------------------------------------------" ;
    
    *output_stream << endl << "Certain conditions:" << endl;
    for(int i=0; i!=aCurrChemComp->iCon.size(); i++)
    {
        if(aCurrChemComp->iCon.at(i)!="") *output_stream << QString::fromStdString(aCurrChemComp->iCon.at(i)) << "\t" << QString::fromStdString(aCurrChemComp->constraints.at(i)) << endl;
    }

    *output_stream << endl << "Uncertain conditions:" << endl;
    for(int i=0; i!=aCurrChemComp->iCon_uncertain.size(); i++)
    {
        if(aCurrChemComp->iCon_uncertain.at(i)!="") *output_stream << QString::fromStdString(aCurrChemComp->iCon_uncertain.at(i)) << ": " << QString::fromStdString(aCurrChemComp->constraints_uncertain.at(i)) << "\t" << "value: " << aCurrChemComp->constrValue_uncertain[i] << "\t" << "std: " << aCurrChemComp->errors.at(i) <<endl;
    }
    *output_stream << "-----------------------------------------------------------------------------------------------------" ;

    // Write species concentrations and activities
    *output_stream << endl << qSetFieldWidth(22) << left << "name";
    *output_stream << qSetFieldWidth(22) << left << "concentration";
    *output_stream << qSetFieldWidth(22) << left << "gamma";
    *output_stream << qSetFieldWidth(22) << left << "activity";
    *output_stream << qSetFieldWidth(22) << left << "log(activity)";
    *output_stream << qSetFieldWidth(0) << endl;
    *output_stream << "-----------------------------------------------------------------------------------------------------"  << endl;

    map<CSpecies*, int>::iterator it;
    for(int k=0; k!=mSpeciesIndices.size(); k++)
    {
        for(it = this->mSpeciesIndices.begin(); it != this->mSpeciesIndices.end(); it++)
        {
            if(k==it->second)
            {
                *output_stream << qSetFieldWidth(22) << left << it->first->name();
                *output_stream << qSetFieldWidth(22) << qSetRealNumberPrecision(15) << left << aCurrChemComp->mConcentration[it->second];
                *output_stream << qSetFieldWidth(22) << qSetRealNumberPrecision(15) << left << aCurrChemComp->mActivityCoeff[it->second];
                if(it->first->name()=="h2o")
                {
                    *output_stream << qSetFieldWidth(22) << qSetRealNumberPrecision(15) << left << 1.0;
                    *output_stream << qSetFieldWidth(22) << qSetRealNumberPrecision(15) << left << 0;
                }
                else 
                {
                    *output_stream << qSetFieldWidth(22) << qSetRealNumberPrecision(15) << left << aCurrChemComp->mConcentration[it->second] * aCurrChemComp->mActivityCoeff[it->second];
                    *output_stream << qSetFieldWidth(22) << qSetRealNumberPrecision(15) << left << log10(aCurrChemComp->mConcentration[it->second] * aCurrChemComp->mActivityCoeff[it->second]);
                }
        
            *output_stream << qSetFieldWidth(0) << endl;
            }
        }
    }
    *output_stream << "-----------------------------------------------------------------------------------------------------" ;

    outputFile->close();

}

void 
CLocalChemicalSystem::WriteRISAInfos(CChemicalComposition* aChemComp, int& iter, double& objFun, VectorXd& grad, double& norm, double& maxRelError, QString name)
{
    QString fileName = "output_RISA_" + aChemComp->name() + ".out";
    fileName = QDir(fileName).absolutePath();

    QFile* outputFile = new QFile(fileName);

    // Open the file, if is not open
    if(!outputFile->isOpen())
    {
        if ( !outputFile->open(QFile::WriteOnly | QIODevice::Append))
        {
            this->PrThrow(ERROR_FILEOPEN);
        }
    }

    QTextStream* output_stream = new QTextStream(outputFile);

    if(iter == 0)
    {
        QTextStream* output_stream = new QTextStream(outputFile);

        // Write name of chemical composition
        *output_stream << "-----------------------------------------------------------------------------------------------------" ;
        *output_stream << endl << "-----------------------------------------------------------------------------------------------------" ;
        *output_stream << endl << "RISA INFOS" << endl;
        *output_stream << "-----------------------------------------------------------------------------------------------------" ;
        
        *output_stream << endl << "CHEMICAL COMPOSITION:" << "\t" << aChemComp->name() << endl;
        *output_stream << "Certain conditions:" << endl;
        for(int i=0; i!=aChemComp->iCon.size(); i++)
        {
            if(aChemComp->iCon.at(i)!="") *output_stream << QString::fromStdString(aChemComp->iCon.at(i)) << "\t" << QString::fromStdString(aChemComp->constraints.at(i)) << endl;
        }

        *output_stream << endl << "Uncertain conditions:" << endl;
        for(int i=0; i!=aChemComp->iCon_uncertain.size(); i++)
        {
            if(aChemComp->iCon_uncertain.at(i)!="") *output_stream << QString::fromStdString(aChemComp->iCon_uncertain.at(i)) << ": " << QString::fromStdString(aChemComp->constraints_uncertain.at(i)) << "\t" << "value: " << aChemComp->constrValue_uncertain[i] << "\t" << "std: " << aChemComp->errors.at(i) <<endl;
        }
        *output_stream << "-----------------------------------------------------------------------------------------------------" ;

        // Write species concentrations and activities
        *output_stream << endl << qSetFieldWidth(15) << left << "iter";
        *output_stream << qSetFieldWidth(15) << left << "Obj. Fun.";
        *output_stream << qSetFieldWidth(15) << left << "Grad.";
        *output_stream << qSetFieldWidth(15) << left << "Norm";
        *output_stream << qSetFieldWidth(15) << left << "Max Rel. Err.";
        *output_stream << qSetFieldWidth(15) << left << "Rel. Species";
        *output_stream << qSetFieldWidth(0) << endl;
        *output_stream << "-----------------------------------------------------------------------------------------------------"  << endl;
    }

    *output_stream << qSetFieldWidth(15) << left << iter;
    *output_stream << qSetFieldWidth(15) << left << objFun;
    *output_stream << qSetFieldWidth(15) << left << grad(0); // Think about when gradient.size() > 1!
    *output_stream << qSetFieldWidth(15) << left << norm;
    *output_stream << qSetFieldWidth(15) << left << maxRelError;
    *output_stream << qSetFieldWidth(15) << left << name;
    *output_stream << qSetFieldWidth(0) << endl;

    outputFile->close();

}

void 
CLocalChemicalSystem::WriteChemCompInfos_conc(CChemicalComposition* aCurrChemComp, bool isHeader, bool chemCompOnSeparateFiles)
{
    QString fileName; 

    if(chemCompOnSeparateFiles)
    {
        fileName = "output_" + aCurrChemComp->name() + "_conc.out";
    }
    else
    {
        fileName = "output_conc.out";
    }

    fileName = QDir(fileName).absolutePath();

    QFile* outputFile = new QFile(fileName);

    // Open the file, if is not open
    if(!outputFile->isOpen())
    {
        if ( !outputFile->open(QFile::WriteOnly | QIODevice::Append))
        {
            this->PrThrow(ERROR_FILEOPEN);
        }
    }

    QTextStream* output_stream = new QTextStream(outputFile);

    if(isHeader)
    {
        // Write name of chemical composition
        *output_stream << endl << "-----------------------------------------------------------------------------------------------------" ;
        *output_stream << endl << "CONCENTRATIONS OF CHEMICAL COMPOSITION:" << "\t" << aCurrChemComp->name() << endl;
        *output_stream << "-----------------------------------------------------------------------------------------------------" ;
        *output_stream << "-----------------------------------------------------------------------------------------------------" << endl;
        map<CSpecies*, int>::iterator it;
        for(int k=0; k!=mSpeciesIndices.size(); k++)
        {
            for(it = this->mSpeciesIndices.begin(); it != this->mSpeciesIndices.end(); it++)
            {
                if(k==it->second) *output_stream << qSetFieldWidth(22) << left << it->first->name();
            }
        }
        *output_stream << qSetFieldWidth(0) << endl;
    }
    else
    {
        map<CSpecies*, int>::iterator it;
        for(int k=0; k!=mSpeciesIndices.size(); k++)
        {
            for(it = this->mSpeciesIndices.begin(); it != this->mSpeciesIndices.end(); it++)
            {
                if(k==it->second) *output_stream << qSetFieldWidth(22) << qSetRealNumberPrecision(15) << left << aCurrChemComp->mConcentration[it->second];
            }
        }
        *output_stream << qSetFieldWidth(0) << endl;
        
    }

    outputFile->close();

}

void 
CLocalChemicalSystem::WriteChemCompInfos_gamma(CChemicalComposition* aCurrChemComp, bool isHeader, bool chemCompOnSeparateFiles)
{
    QString fileName; 

    if(chemCompOnSeparateFiles)
    {
        fileName = "output_" + aCurrChemComp->name() + "_act_coeffs.out";
    }
    else
    {
        fileName = "output_act_coeffs.out";
    }

    fileName = QDir(fileName).absolutePath();

    QFile* outputFile = new QFile(fileName);

    // Open the file, if is not open
    if(!outputFile->isOpen())
    {
        if ( !outputFile->open(QFile::WriteOnly | QIODevice::Append))
        {
            this->PrThrow(ERROR_FILEOPEN);
        }
    }

    QTextStream* output_stream = new QTextStream(outputFile);

    if(isHeader)
    {
        // Write name of chemical composition
        *output_stream << endl << "-----------------------------------------------------------------------------------------------------" ;
        *output_stream << endl << "ACTIVITY COEFFICIENTS OF CHEMICAL COMPOSITION:" << "\t" << aCurrChemComp->name() << endl;
        *output_stream << "-----------------------------------------------------------------------------------------------------" ;
        *output_stream << "-----------------------------------------------------------------------------------------------------" << endl;
        map<CSpecies*, int>::iterator it;
        for(int k=0; k!=mSpeciesIndices.size(); k++)
        {
            for(it = this->mSpeciesIndices.begin(); it != this->mSpeciesIndices.end(); it++)
            {
                if(k==it->second) *output_stream << qSetFieldWidth(22) << left << it->first->name();
            }
        }
        *output_stream << qSetFieldWidth(0) << endl;
    }
    else
    {
        map<CSpecies*, int>::iterator it;
        for(int k=0; k!=mSpeciesIndices.size(); k++)
        {
            for(it = this->mSpeciesIndices.begin(); it != this->mSpeciesIndices.end(); it++)
            {
                if(k==it->second) *output_stream << qSetFieldWidth(22) << qSetRealNumberPrecision(15) << left << aCurrChemComp->mActivityCoeff[it->second];
            }
        }
        *output_stream << qSetFieldWidth(0) << endl;
        
    }

    outputFile->close();

}

VectorXd 
CLocalChemicalSystem::Calculate_utot_classic_SIA(CChemicalComposition* aCurrChemComp, vector<double> &u_aq)
{
    VectorXd utot = VectorXd::Zero(u_aq.size());

    // - Access the aqueous phase
    it1 = this->mPhases.equal_range(eAqueousPhase); 

    // - Get the name of the phase
    string phaseName = it1.first->second->name().toStdString();

    // - Retrieve the volumetric fraction from CChemicalComposition
    map<string,double>::iterator it3;
    it3 = aCurrChemComp->mVolumetricFractions.find(phaseName);
    double theta_liq = it3->second;

    for(int i=0; i!=utot.size(); i++)
    {
        utot(i) = theta_liq * u_aq[i];
    }

    // Add contribution of other phases

    // - Get u_min 
    vector<double> u;
    u = this->Get_u_min_classic(aCurrChemComp);

    if(u.size()!=0)
    {
        // - Access the mineral phase
        it1 = this->mPhases.equal_range(eMineralPhase); 
        double theta = 0;

        if(it1.first!=this->mPhases.end()) 
        {
            for (it=it1.first; it!=it1.second; it++) // OSS: it won't work for more than one mineral phase...same for ads and gases
            {
                // - Get the name of the phase
                phaseName = it1.first->second->name().toStdString();

                // - Retrieve the volumetric fraction from CChemicalComposition
                it3 = aCurrChemComp->mVolumetricFractions.find(phaseName);
                theta = it3->second;
            }
        }

        // - Add contribution to residual vector
        for(int i=0; i!=utot.size(); i++)
        {
            utot(i) += u[i]*theta;
        }
        u.resize(0);
        theta = 0;
    }

    // - Get u_gas
    u = this->Get_u_gas_classic(aCurrChemComp);

    // - Add contribution to residual vector
    if(u.size()!=0)
    {
        // - Access the gas phase
        it1 = this->mPhases.equal_range(eGasPhase); 
        double theta = 0;

        if(it1.first!=this->mPhases.end()) 
        {
            for (it=it1.first; it!=it1.second; it++) // OSS: it won't work for more than one mineral phase...same for ads and gases
            {
                // - Get the name of the phase
                phaseName = it1.first->second->name().toStdString();

                // - Retrieve the volumetric fraction from CChemicalComposition
                it3 = aCurrChemComp->mVolumetricFractions.find(phaseName);
                theta = it3->second;
            }
        }

        // - Add contribution to residual vector
        for(int i=0; i!=utot.size(); i++)
        {
            utot(i) += u[i]*theta;
        }
        u.resize(0);
        theta = 0;
    }

    // - Get u_surf 
    u = this->Get_u_surf_classic(aCurrChemComp);

    // - Add contribution to residual vector
    if(u.size()!=0)
    {
        // - Access the surface phase
        it1 = this->mPhases.equal_range(eSurfacePhase); 
        double theta = 0;

        if(it1.first!=this->mPhases.end()) 
        {
            for (it=it1.first; it!=it1.second; it++) // OSS: it won't work for more than one mineral phase...same for ads and gases
            {
                // - Get the name of the phase
                phaseName = it1.first->second->name().toStdString();

                // - Retrieve the volumetric fraction from CChemicalComposition
                it3 = aCurrChemComp->mVolumetricFractions.find(phaseName);
                theta = it3->second;
            }
        }

        // - Add contribution to residual vector
        for(int i=0; i!=utot.size(); i++)
        {
            utot(i) += u[i]*theta;
        }
        u.resize(0);
        theta = 0;
    }

    // Multiply times f_tot
    utot = utot;

    return utot;

}

void
CLocalChemicalSystem::Speciate_utot_classic(CChemicalComposition* aCurrChemComp)
{
}

void
CLocalChemicalSystem::CalculateCASConcentration(CChemicalComposition* aCurrChemComp)
{
    // Get vector of secondary species
    VectorXd c2 = VectorXd::Zero(this->mSecSpeciesNum);
    aCurrChemComp->GetSetOfConcVector(c2,mRelSecondarySpcIndices);
        
    VectorXd theta_c = VectorXd::Zero(this->mCAS.size());
    theta_c = - this->mSecMod.transpose() * c2;

    // Get water index to skip it


    for(int i=0; i!=this->mRelConstActSpcIndices.size(); i++)
    {
        // OJO: mConcentration has to contain concentration or molarity...I still have to deduce the theta...think about how to do it.
        aCurrChemComp->mConcentration[i+this->mRelPrimSpcIndicesWithoutCA.size()+this->mRelSecondarySpcIndices.size()] = theta_c(i);
    }
}

VectorXd
CLocalChemicalSystem::Calculate_ured_from_utot_classic(CChemicalComposition* aCurrChemComp, VectorXd u_tot_classic)
{
    VectorXd result; // vector containing the "reduced" components

    result.resize(this->mPrimSpRedNum);

    // Declare elimination matrix
    MatrixXd e_matrix;

    map<string,int>::iterator it5;
    it5 = this->mClassicComponents.find("h2o") ;
    if(it5 != this->mClassicComponents.end())
    {
        // Get rid of the first and last column of the "classic" component matrix (which is the water component) and put the new matrix in "comp_matrix_without_water"
        e_matrix = MatrixXd::Zero(this->mPrimSpRedNum, this->mSpeciesIndices.size() - this->mSecSpeciesNum - 1);

        MatrixXd comp_matrix_without_water = MatrixXd::Zero(this->mComponentMatrixInit.rows()-1, this->mComponentMatrixInit.cols()-1);
         
        for(int i=0; i!=comp_matrix_without_water.rows(); i++)
        {
            for(int j=0; j!=comp_matrix_without_water.cols(); j++)
            {
                comp_matrix_without_water(i,j) = this->mComponentMatrixInit(i+1, j+1);
            }
        }

        // Calculate elimination  matrix
        e_matrix = this->mComponentMatrix * this->mComponentMatrixInit.transpose() * (this->mComponentMatrixInit*this->mComponentMatrixInit.transpose()).inverse();

        comp_matrix_without_water.resize(0,0);
    }
    else
    {
        e_matrix = MatrixXd::Zero(this->mPrimSpRedNum, this->mSpeciesIndices.size() - this->mSecSpeciesNum);

        e_matrix = this->mComponentMatrix * this->mComponentMatrixInit.transpose() * (this->mComponentMatrixInit*this->mComponentMatrixInit.transpose()).inverse();
    }

    result = e_matrix * u_tot_classic;

    return result;
}