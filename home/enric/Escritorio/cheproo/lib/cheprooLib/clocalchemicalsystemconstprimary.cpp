#include "clocalchemicalsystemconstprimary.h"
#include "cglobalchemicalsystem.h"


CLocalChemicalSystemConstPrimary::CLocalChemicalSystemConstPrimary(void)
{
    this->mClassName = "CLocalChemicalSystemConstPrimary";
}


CLocalChemicalSystemConstPrimary::CLocalChemicalSystemConstPrimary(QDomElement aNode)
: CLocalChemicalSystem(aNode)
{
    this->mClassName = "CLocalChemicalSystemConstPrimary";
}


CLocalChemicalSystemConstPrimary::~CLocalChemicalSystemConstPrimary(void)
{
}

void
CLocalChemicalSystemConstPrimary::SetLocal(vector<string>& primSpecies, set<string>& cas)
{
    // Add species to the map of maps mSpecies
    vector<CSpecies*> auxVector1;  // Aqueous primary + cas
    vector<CSpecies*> auxVector2;  // Surface primary
    vector<CSpecies*> auxVector3;  // Mineral primary	
    vector<CSpecies*> auxVector4;  // Gas primary
    vector<CSpecies*> auxVector9;  // Aqueous primary without constant activity species
    vector<CSpecies*> auxVector10; // Aqueous primary constant activity species (typically, water or h+)
    vector<CSpecies*> constActSpcToLocChemSys;

    set<string>::iterator myit, myitFind;
	set<string> primSpeciesAux;
    int casPrimary = 0;
    bool isProtonCAS = false;

    for (int j=0; j!=primSpecies.size(); j++)
	{
		// Retrieve the species name
		const string aSpeciesName = primSpecies.at(j);

        // Check if h+ is a CAS
        if(aSpeciesName == "h+")
        {
            myitFind = cas.find(aSpeciesName);
            if(myitFind != cas.end()) isProtonCAS = true;
        }

        // Check if it's water then is a CAS
        if(aSpeciesName == "h2o" || isProtonCAS == true) // Careful if I give h+ a cTot condition then is not true!
        {
            this->mConstActSpeciesNum++;

            // Add it to the constant activity species pointers vector and to the auxVector10 (vector of primary aq constant activity species)
            constActSpcToLocChemSys.push_back(mGlobChemSysPointer->mAllSpecies[aSpeciesName]);
            mGlobChemSysPointer->mAllSpecies[aSpeciesName]->mIsConstActivity = true;

            auxVector10.push_back(mGlobChemSysPointer->mAllSpecies[aSpeciesName]);
            auxVector1.push_back(mGlobChemSysPointer->mAllSpecies[aSpeciesName]);
            primSpeciesAux.insert(aSpeciesName);
            casPrimary++;

            // Update indices
            this->mRelPrimSpcIndicesWithoutCA.push_back(j);
            // Fill vectors with primary species indices, with and without constant activity species
            this->mRelPrimarySpcIndices.push_back(j);

        }
        else
        {
            // substitute mPrim...size() with mPrimSpeciesNum, everywhere!
            mPrimRedSpecies.push_back(mGlobChemSysPointer->mAllSpecies[aSpeciesName]);
		    primSpeciesAux.insert(aSpeciesName);

            // Retrieve this species from the mAllSpecies list in GlobalChemicalSystem and put it in the auxVectors
            if(mGlobChemSysPointer->mAllSpecies[aSpeciesName]->mSpeciesKind == eAqueous)
            {
                auxVector1.push_back(mGlobChemSysPointer->mAllSpecies[aSpeciesName]);

                // If it's a constant activity species
                myitFind = cas.find(aSpeciesName);
                if(myitFind!=cas.end()) 
                {
                    constActSpcToLocChemSys.push_back(mGlobChemSysPointer->mAllSpecies[aSpeciesName]);
                    auxVector10.push_back(mGlobChemSysPointer->mAllSpecies[aSpeciesName]);
                }
                else
                {
                    if(!mGlobChemSysPointer->mAllSpecies[aSpeciesName]->mIsConstActivity)
                    {
                        auxVector9.push_back(mGlobChemSysPointer->mAllSpecies[aSpeciesName]);
                    }

                    // Update indices
                    this->mRelPrimSpcIndicesWithoutCA.push_back(j);
                }

            }
            else if(mGlobChemSysPointer->mAllSpecies[aSpeciesName]->mSpeciesKind == eSurface)
            {
                auxVector2.push_back(mGlobChemSysPointer->mAllSpecies[aSpeciesName]);

                // Update indices
                this->mRelPrimSpcIndicesWithoutCA.push_back(j);
            }

            // Start filling the map of species indices
            this->mSpeciesIndices.insert(pair<CSpecies*, int>(mGlobChemSysPointer->mAllSpecies[aSpeciesName],j-casPrimary));

            // Update number of primary species
            this->mPrimSpeciesNum = j+1;

            // Fill vectors with primary species indices, with and without constant activity species
            this->mRelPrimarySpcIndices.push_back(j);

        }
    }

    if(auxVector2.size() != 0)
    {
        this->mSpecies.insert(pair<eSpcToPhaseKind, vector<CSpecies*> >(eSurface1,auxVector2));
    }


    // Keep filling the map of maps mSpecies with the remaining species
    vector<CSpecies*> auxVector5; // Aq secondary
    vector<CSpecies*> auxVector6; // Surf secondary
    vector<CSpecies*> auxVector7; // Min secondary
    vector<CSpecies*> auxVector8; // Gas secondary

    int i = this->mPrimSpeciesNum - casPrimary;

    map<string,CSpecies*>::iterator it;
    for(it = this->mGlobChemSysPointer->mAllSpecies.begin() ; it != this->mGlobChemSysPointer->mAllSpecies.end() ; it++)
    {   
        // Check if the current species is already present in a reaction
        bool isNotPresent = true;
        bool isNotWaterThere = true;
        map<string,CReaction*>::iterator it1;
        map<CSpecies*,int>::iterator it2;
        for(it1 = this->mGlobChemSysPointer->mAllReactions.begin(); it1!=this->mGlobChemSysPointer->mAllReactions.end(); it1++)
        {
            int m=0;
            CReaction* currReaction = it1->second;
            m = currReaction->mRelSpeciesVsStoichCoeffs.count(it->second);
            if(m!=0 && it1->first==it->first) isNotPresent = false;

            // Check if water reaction is already present
            if(it->first == "h2o" && it1->first == "oh-") isNotWaterThere = false; 
            if(it->first == "oh-" && it1->first == "h2o") isNotWaterThere = false;

        }
        // Check if I can find this species in the vector of primary species names
        myitFind = primSpeciesAux.find(it->first); 
        if(myitFind==primSpeciesAux.end() && isNotPresent && isNotWaterThere) // Read the reaction if is not present
        {
            QString aReactionName = QString::fromStdString(it->first);
            this->mGlobChemSysPointer->CreateReaction(aReactionName);
            this->mEqReactions.push_back(this->mGlobChemSysPointer->mAllReactions[it->first]);

        }

        // Add the reaction to the corresponding vector
        //if(this->mGlobChemSysPointer->mAllReactions[it->first]->mIsKinetic)
        //{
        //    this->mEqReactions.push_back(this->mGlobChemSysPointer->mAllReactions[it->first]);
        //}
        //else
        //{
  /*      map<string, CReaction*>::iterator itreac = this->mGlobChemSysPointer->mAllReactions.find(it->first);
        if(itreac!=this->mGlobChemSysPointer->mAllReactions.end())*/
            //this->mEqReactions.push_back(this->mGlobChemSysPointer->mAllReactions[it->first]);
        //}

        if(myitFind==primSpeciesAux.end())
        {
            if(it->second->mSpeciesKind == eAqueous)
            {
                // If it's a constant activity species, update the attribute mConstActSpeciesNum
                if(it->second->mIsConstActivity) 
                {
                    this->mConstActSpeciesNum++;

                    // Add it to the constant activity species pointers vector and to the auxVector10 (vector of primary aq constant activity species)
                    constActSpcToLocChemSys.push_back(it->second);

                    auxVector10.push_back(it->second);
                    auxVector1.push_back(it->second);

                    // Delete one of the primary species and add it to the secondary vectors, depending on their mSpeciesKind
                    if(mPrimRedSpecies.back()->mSpeciesKind == eAqueous)
                    {
                        auxVector5.push_back(mPrimRedSpecies.back());
                    }
                    else if(mPrimRedSpecies.back()->mSpeciesKind == eSurface)
                    {
                        auxVector6.push_back(mPrimRedSpecies.back());
                    }

                    // Erase the last element that I added to secondary species vector
                    mPrimRedSpecies.pop_back();

                }
                else
                {
					// Continue filling the map of species indices
					this->mSpeciesIndices.insert(pair<CSpecies*, int>(it->second,i));

                    this->mRelSecondarySpcIndices.push_back(i);

					i++;

                    auxVector5.push_back(it->second);
                }
                    

            }
            else if(it->second->mSpeciesKind == eSurface)
            {
                auxVector6.push_back(it->second);
            }
            else if(it->second->mSpeciesKind == eMineral)
            {
                // Check if the mineral is CAS or not (THIS WON'T WORK FOR NON PURE MINERAL PHASES, OJO!)
				// Find the phase with gypsum
				map<string, CPhase*>::iterator itPhase;
				for(myit=cas.begin(); myit!=cas.end(); myit++)  // Loop over phases defined in cas
				{
					itPhase = this->mGlobChemSysPointer->mAllPhases.find(*myit); // Find the cas phase in global phases
					if(itPhase!=this->mGlobChemSysPointer->mAllPhases.end())
					{
						if(itPhase->second->mSpeciesNames.at(0) == it->second->name().toStdString()) it->second->mIsConstActivity = true;
						
					}
				}

				// If it's a constant activity species, update the attribute mConstActSpeciesNum
                if(it->second->mIsConstActivity) 
                {
                    this->mConstActSpeciesNum++;

					// Add it to the constant activity species pointers vector and to the auxVector3 (vector of primary species)
                    constActSpcToLocChemSys.push_back(it->second);

                    // Delete one of the primary species and add it to the secondary vectors, depending on their mSpeciesKind
                    if(mPrimRedSpecies.back()->mSpeciesKind == eAqueous)
                    {
                        auxVector5.push_back(mPrimRedSpecies.back());
                    }
                    else if(mPrimRedSpecies.back()->mSpeciesKind == eSurface)
                    {
                        auxVector6.push_back(mPrimRedSpecies.back());
                    }

                    // Erase the last element that I added to secondary species vector
                    mPrimRedSpecies.pop_back();

                    auxVector3.push_back(it->second);
                }
                else
                {
                    auxVector7.push_back(it->second);
                }

                

                }

                else if(it->second->mSpeciesKind == eGas)
                {
                    // Check if the mineral is CAS or not
                    myitFind = cas.find(it->second->name().toStdString());         // If it's present in the minerals of the chemical composition
                    if(myitFind!=cas.end()) it->second->mIsConstActivity = true;   // then is a CAS

                    // If it's a constant activity species, update the attribute mConstActSpeciesNum
                    if(it->second->mIsConstActivity) 
                    {
                        this->mConstActSpeciesNum++;

					    // Add it to the constant activity species pointers vector and to the auxVector4 (vector of primary species)
                        constActSpcToLocChemSys.push_back(it->second);

                        // Delete one of the primary species and add it to the secondary vectors, depending on their mSpeciesKind
                        if(mPrimRedSpecies.back()->mSpeciesKind == eAqueous)
                        {
                            auxVector5.push_back(mPrimRedSpecies.back());
                        }
                        else if(mPrimRedSpecies.back()->mSpeciesKind == eSurface)
                        {
                            auxVector6.push_back(mPrimRedSpecies.back());
                        }

                        // Erase the last element that I added to secondary species vector
                        mPrimRedSpecies.pop_back();

                        auxVector4.push_back(it->second);

                     }
                    else
                    {
                        auxVector8.push_back(it->second);
                    }

                }  

            }

		}   

    // Fix secondary non aqueous (they must have higher indices, it's simpler for "classic" components definition)
        
        if(auxVector6.size()!=0)
        {
            int maxIndex = this->mSpeciesIndices.size();
            for(int i=0; i!=auxVector6.size(); i++)
            {
                this->mSpeciesIndices.insert(pair<CSpecies*, int>(auxVector6.at(i),i + maxIndex));
            }
        }

        if(auxVector7.size()!=0)
        {
            int maxIndex = this->mSpeciesIndices.size();
            for(int j=0; j!=auxVector7.size(); j++)
            {
                this->mSpeciesIndices.insert(pair<CSpecies*, int>(auxVector7.at(j),j + maxIndex));
            }
        }

        if(auxVector8.size()!=0)
        {
            int maxIndex = this->mSpeciesIndices.size();
            for(int k=0; k!=auxVector8.size(); k++)
            {
                this->mSpeciesIndices.insert(pair<CSpecies*, int>(auxVector8.at(k),k + maxIndex));
            }
        }

				

	// Now fix the constant activity species
	if(auxVector2.size()!=0)
	{
		// erase primary surface
	}
	else
	{
		// Erase the last elements (Nc)
		for(int i=0; i!=constActSpcToLocChemSys.size()-casPrimary; i++)
		{
			auxVector1.pop_back();
			auxVector9.pop_back();

		}
	}

	int spcNum = mSpeciesIndices.size();
	for(int i=0; i!=constActSpcToLocChemSys.size(); i++)
	{
		// Continue filling the map of species indices
		this->mSpeciesIndices.insert(pair<CSpecies*, int>(constActSpcToLocChemSys[i],i + spcNum));

        // Update indices of constant activity species
        this->mRelConstActSpcIndices.push_back(i + spcNum);

        // Update primary species indices (reduced)
        this->mRelPrimSpcIndicesWithoutCA.pop_back();

	}  

        // Update secondary species indices
        this->mRelSecondarySpcIndices.clear();
        for(int i=0; i!=this->mEqReactions.size(); i++)
        {
            this->mRelSecondarySpcIndices.push_back(i+this->mRelPrimSpcIndicesWithoutCA.size());
        }


    this->mSpecies.insert(pair<eSpcToPhaseKind, vector<CSpecies*> >(eAqueous1,auxVector1));

    this->mSpecies.insert(pair<eSpcToPhaseKind, vector<CSpecies*> >(eAqueous1nc,auxVector9));

    this->mSpecies.insert(pair<eSpcToPhaseKind, vector<CSpecies*> >(eAqueous1c,auxVector10));

    this->mSpecies.insert(pair<eSpcToPhaseKind, vector<CSpecies*> >(eAqueous2,auxVector5));

    if(auxVector3.size() != 0)
    {
        this->mSpecies.insert(pair<eSpcToPhaseKind, vector<CSpecies*> >(eMineral1,auxVector3));
    }

    if(auxVector4.size() != 0)
    {
        this->mSpecies.insert(pair<eSpcToPhaseKind, vector<CSpecies*> >(eGas1,auxVector4));
    }

    if(auxVector6.size() != 0)
    {
        this->mSpecies.insert(pair<eSpcToPhaseKind, vector<CSpecies*> >(eSurface2,auxVector6));
    }
    if(auxVector7.size() != 0)
    {
        this->mSpecies.insert(pair<eSpcToPhaseKind, vector<CSpecies*> >(eMineral2,auxVector7));
    }
    if(auxVector8.size() != 0)
    {
        this->mSpecies.insert(pair<eSpcToPhaseKind, vector<CSpecies*> >(eGas2,auxVector8));
    }

    auxVector1.clear();
    auxVector2.clear();
    auxVector3.clear();
    auxVector4.clear();
    auxVector5.clear();
    auxVector6.clear();
    auxVector7.clear();
    auxVector8.clear();

    // Check if S2 is invertible
	this->mS2Inverse = MatrixXd::Zero(mEqReactions.size(), mEqReactions.size());

    // Compute S2
    this->ComputeS2(mS2Inverse);

    if(mS2Inverse.determinant()==0) 
    {
        // Then call a method to change primary and secondary species..to do!
    }

    // Add phases to LocalChemicalSystem 
    for(map<string,CPhase*>::iterator it = this->mGlobChemSysPointer->mAllPhases.begin() ; it != this->mGlobChemSysPointer->mAllPhases.end() ; it++)
    {
        if(it->second->mPhaseKind == eMineralPhase || it->second->mPhaseKind == eGasPhase)
        {
            set<string>::iterator myitFind2;
			myitFind2 = cas.find(it->first);
            if(myitFind2!=cas.end())
            {
                // Then add the phase to the local chemical system
                mPhases.insert(pair<ePhaseKind, CPhase* >(it->second->mPhaseKind,it->second));
            }
            else
            {
                mPhases.insert(pair<ePhaseKind, CPhase* >(it->second->mPhaseKind,it->second));
            }
        }
        else
        {
            mPhases.insert(pair<ePhaseKind, CPhase* >(it->second->mPhaseKind,it->second));
        }

    }

    // Assign attribute of constant activity species defined in the input
    this->mCAS = cas;

	// Initialize number of reduced primary and secondary species
    this->mPrimSpRedNum = this->mPrimSpeciesNum - this->mConstActSpeciesNum;
	this->mSecSpeciesNum = this->mEqReactions.size();

    // Declare and assign dimensions to dc2/dc1, the system matrix and the RHS 
    this->dc2_dc1 = MatrixXd::Zero(this->mSecSpeciesNum,this->mPrimSpRedNum);
    this->mSystemRHS = MatrixXd::Zero(this->mSecSpeciesNum, this->mPrimSpRedNum);
    this->mSystemMatrix = MatrixXd::Zero(this->mSecSpeciesNum, this->mSecSpeciesNum);

    // Declare jacbian and residual for N-R to evaluate c1
	this->jacobian = MatrixXd::Zero(this->mPrimSpRedNum, this->mPrimSpRedNum);
    this->residual = VectorXd::Zero(this->mPrimSpRedNum);

    // Compute the stoichiometric matrix S1* = -(1/S2)S1
    this->ComputeS1Mod();

    // Compute the vector k* = (1/S2)k
    this->ComputeLogkMod();

    // Compute component matrix
    this->ComputeComponentMatrix();

}

void
 CLocalChemicalSystemConstPrimary::ComputeS1Mod()
{
   // S1* initialization
    this->mS1Mod = MatrixXd::Zero(mEqReactions.size(), this->mPrimSpeciesNum - this->mConstActSpeciesNum);

    // Sec* initialization
    this->mSecMod = MatrixXd::Zero(mEqReactions.size(), this->mConstActSpeciesNum);

    if(this->mS2Inverse.size()==0) // If S2 hasn't been computed yet then do it now
    {
	    // Declare and allocate S2
	    this->mS2Inverse = MatrixXd::Zero(mEqReactions.size(), mEqReactions.size());

        // Compute S2
        this->ComputeS2(mS2Inverse);
    }

	// Check if S2 is -I
	if(this->mS2Inverse != MatrixXd::Identity(mEqReactions.size(), mEqReactions.size()) * -1.0)
	{
		// Then build mS1mod = -(1/S2)*S1

		// Declare and allocate S1 for aqueous species
		MatrixXd aS1ea = MatrixXd::Zero(mEqReactions.size(), this->mPrimSpeciesNum - this->mConstActSpeciesNum);

        // Loop over primary species indices (the species must be in order of index in S2)
        for(int i=0; i!=this->mPrimSpRedNum; i++)
        {
            // Get the species correspondent to the index
            for(map<CSpecies*, int>::iterator it= this->mSpeciesIndices.begin() ; it != this->mSpeciesIndices.end() ; it++)
            {
                if(mSpeciesIndices[it->first] == mRelPrimSpcIndicesWithoutCA[i])
                {
                    for ( int j = 0 ; j!=mEqReactions.size() ; j++ ) 
	                {
                        // Obtain pointer to current reaction
                        CReaction* aCurrReaction = mEqReactions.at(j);
    			        
                        // Obtain stoichiometric coefficient for the current species in the current reaction
			            double stoichCoeff = aCurrReaction->mRelSpeciesVsStoichCoeffs[it->first];

			            // Put it in the matrix
			            aS1ea(j,i) = stoichCoeff;

                    }
                 }
             }
         } 

		// Then build mSecmod = -(1/S2)*Sec

		// Declare and allocate S1 for aqueous species
        MatrixXd aSec = MatrixXd::Zero(mEqReactions.size(), this->mConstActSpeciesNum);

        if(this->mConstActSpeciesNum != 0)
		{
            // Loop over constant activity species indices (the species must be in order of index in S2)
            for(int i=0; i!=this->mConstActSpeciesNum; i++)
            {
                // Get the species correspondent to the index
                for(map<CSpecies*, int>::iterator it= this->mSpeciesIndices.begin() ; it != this->mSpeciesIndices.end() ; it++)
                {
                    if(mSpeciesIndices[it->first] == mRelConstActSpcIndices[i])
                    {
                        for ( int j = 0 ; j!=mEqReactions.size() ; j++ ) 
	                    {
                            // Obtain pointer to current reaction
                            CReaction* aCurrReaction = mEqReactions.at(j);
    			        
                            // Obtain stoichiometric coefficient for the current species in the current reaction
			                double stoichCoeff = aCurrReaction->mRelSpeciesVsStoichCoeffs[it->first];

			                // Put it in the matrix
			                aSec(j,i) = stoichCoeff;

                        }
                     }
                 }
             } 
        }


		// Evaluate the inverse of S2
		this->mS2Inverse = this->mS2Inverse.inverse();

		//Perform the multiplication to evaluate mS1Mod
		this->mS1Mod = - this->mS2Inverse * aS1ea;

        //Perform the multiplication to evaluate mSecMod
		if(this->mConstActSpeciesNum != 0) 
        {
            this->mSecMod = - this->mS2Inverse * aSec;
            aSec.resize(0,0);
        }
		aS1ea.resize(0,0);
	}
	else
	{
		// Build mS1mod = S1ec 

        // Loop over primary species indices (the species must be in order of index in S2)
        for(int i=0; i!=mRelPrimarySpcIndices.size() - this->mConstActSpeciesNum; i++)
        {
            // Get the species correspondent to the index
            for(map<CSpecies*, int>::iterator it= this->mSpeciesIndices.begin() ; it != this->mSpeciesIndices.end() ; it++)
            {
                if(mSpeciesIndices[it->first] == mRelPrimarySpcIndices[i])
                {
                    for ( int j = 0 ; j!=mEqReactions.size() ; j++ ) 
	                {
                        // Obtain pointer to current reaction
                        CReaction* aCurrReaction = mEqReactions.at(j);
    			        
                        // Obtain stoichiometric coefficient for the current species in the current reaction
			            double stoichCoeff = aCurrReaction->mRelSpeciesVsStoichCoeffs[it->first];

			            // Put it in the matrix
			            this->mS1Mod(j,i) = stoichCoeff;

                    }
                 }
             }
         }

		// Build mSecmod = Sec 

        // Loop over constant activity species indices (the species must be in order of index in S2)
        if(this->mConstActSpeciesNum != 0)
        {
            for(int i=0; i!=this->mConstActSpeciesNum; i++)
            {
                // Get the species correspondent to the index
                for(map<CSpecies*, int>::iterator it= this->mSpeciesIndices.begin() ; it != this->mSpeciesIndices.end() ; it++)
                {
                    if(mSpeciesIndices[it->first] == mRelPrimarySpcIndices[this->mPrimSpRedNum + i])
                    {
                        for ( int j = 0 ; j!=mEqReactions.size() ; j++ ) 
	                    {
                            // Obtain pointer to current reaction
                            CReaction* aCurrReaction = mEqReactions.at(j);
    			        
                            // Obtain stoichiometric coefficient for the current species in the current reaction
			                double stoichCoeff = aCurrReaction->mRelSpeciesVsStoichCoeffs[it->first];

			                // Put it in the matrix
			                this->mSecMod(j,i) = stoichCoeff;

                        }
                     }
                 }
             } 
        }

	}
}

void 
CLocalChemicalSystemConstPrimary::ComputeComponentMatrix()
{
    int aSpeciesNum = this->mPrimSpRedNum + mEqReactions.size();

    // Evaluate the component matrix: allocate and initialize 
    this->mComponentMatrix_c = MatrixXd::Zero(this->mConstActSpeciesNum, this->mConstActSpeciesNum + mEqReactions.size());
    this->mComponentMatrixInit =  MatrixXd::Zero(this->mPrimSpeciesNum, this->mSpecies[eAqueous1nc].size() + this->mSpecies[eAqueous1c].size() + this->mSpecies[eAqueous2].size());

    // Calculate the "classic" component matrix
    this->mComponentMatrixInit = this->ComputeComponentMatrixInit();

    // Assign values to component matrix
	this->mComponentMatrix = MatrixXd::Zero(this->mPrimSpRedNum, aSpeciesNum);
    for(int i=0; i!=this->mPrimSpRedNum; i++)
    {
        // Diagonal terms equal to 1
        mComponentMatrix(i,i) = 1.0;

        for(int j=0; j!=mEqReactions.size(); j++)
        {
            this->mComponentMatrix(i,j+this->mPrimSpRedNum) = this->mS1Mod(j,i);
        }
    }

    // Assign values to component matrix for constant activity species
    for(int i=0; i!=this->mConstActSpeciesNum; i++)
    {
        // Diagonal terms equal to 1
        mComponentMatrix_c(i,i) = 1.0;

        for(int j=0; j!=mEqReactions.size(); j++)
        {
            this->mComponentMatrix_c(i,j+this->mConstActSpeciesNum) = this->mSecMod(j,i);
        }
    }

}

MatrixXd 
CLocalChemicalSystemConstPrimary::Get_dc_dc1_ByPhase(QString aVarType)
{
    MatrixXd result;

    vector<CSpecies*> species1;
    vector<CSpecies*> species2;

    vector<int> indices;

    if(aVarType == "component_aq")
	{
		species1 = this->mSpecies[eAqueous1nc];
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

    result = MatrixXd::Zero(indices.size(), this->mPrimSpRedNum);
    result.topLeftCorner(this->mPrimSpRedNum, this->mPrimSpRedNum).setIdentity();

    // Fill part relative to secondary species
    for(int i=0; i!= species2.size(); i++)
    {
        it = this->mSpeciesIndices.find(species2.at(i));
        int sec_index = it->second;

        for(int j=0; j!=this->mPrimSpRedNum; j++)
        {
            result(sec_index,j) = this->dc2_dc1(sec_index - this->mPrimSpRedNum, j);
        }
    }


    return result;
}

MatrixXd
CLocalChemicalSystemConstPrimary::ComputeComponentMatrixInit()
{
	MatrixXd result;
    this->mComponentMatrix = MatrixXd::Zero(this->mPrimSpeciesNum,this->mPrimSpeciesNum + mEqReactions.size());
    MatrixXd stoichMatrix = MatrixXd::Zero(this->mEqReactions.size(),this->mSpeciesIndices.size());
    MatrixXd s2 = MatrixXd::Zero(this->mEqReactions.size(),this->mEqReactions.size());

    int index = 0;

    // Add: check if there is water then put its index to 0 and increase all the others of 1
    CSpecies* h2o;
    int index1 = 0;
    map<CSpecies*, int>::iterator it;
    map<string,CSpecies*>::iterator it1;
    map<CSpecies*,int> tempIndices;
    it1=this->mGlobChemSysPointer->mAllSpecies.find("h2o");

    if(it1!=this->mGlobChemSysPointer->mAllSpecies.end()) // If there's water in the system
    {
        h2o = it1->second;
        it=this->mSpeciesIndices.find(h2o);
        index1 = it->second;

        for(it = this->mSpeciesIndices.begin(); it != this->mSpeciesIndices.end(); it++)
        {
            if(it->second < index1) 
            {
                tempIndices.insert(pair<CSpecies*,int>(it->first,it->second+1));
            }
            else if(it->second==index1) 
            {
                tempIndices.insert(pair<CSpecies*,int>(it->first,0));
            }
            else
            {
                tempIndices.insert(pair<CSpecies*,int>(it->first,it->second));
            }
        }
  

		// Get the species correspondent to the index
		for(it = tempIndices.begin() ; it != tempIndices.end() ; it++)
		{
			index = it->second;
			for ( int j = 0 ; j!=mEqReactions.size() ; j++ ) 
			{
				// Obtain pointer to current reaction
				CReaction* aCurrReaction = mEqReactions.at(j);
    			        
				// Obtain stoichiometric coefficient for the current species in the current reaction
				double stoichCoeff = aCurrReaction->mRelSpeciesVsStoichCoeffs[it->first];

				// Put it in the matrix
				stoichMatrix(j,index) = stoichCoeff;
			}

		}

		s2 = stoichMatrix.block(0,this->mRelPrimarySpcIndices.size(),this->mEqReactions.size(),this->mEqReactions.size());
		s2 = - s2.inverse();

		stoichMatrix = s2 * stoichMatrix;

		mComponentMatrix.block(0,0,this->mRelPrimarySpcIndices.size(),this->mRelPrimarySpcIndices.size()) = MatrixXd::Identity(this->mRelPrimarySpcIndices.size(), this->mRelPrimarySpcIndices.size());
		mComponentMatrix.block(0,this->mRelPrimarySpcIndices.size(),this->mRelPrimarySpcIndices.size(),this->mEqReactions.size()) = stoichMatrix.topLeftCorner(this->mEqReactions.size(), this->mRelPrimarySpcIndices.size()).transpose();

		result = this->GetComponentInit(tempIndices);
	}
	else
	{
		// Get the species correspondent to the index
		for(it = this->mSpeciesIndices.begin() ; it != this->mSpeciesIndices.end() ; it++)
		{
			index = it->second;
			for ( int j = 0 ; j!=mEqReactions.size() ; j++ ) 
			{
				// Obtain pointer to current reaction
				CReaction* aCurrReaction = mEqReactions.at(j);
    			        
				// Obtain stoichiometric coefficient for the current species in the current reaction
				double stoichCoeff = aCurrReaction->mRelSpeciesVsStoichCoeffs[it->first];

				// Put it in the matrix
				stoichMatrix(j,index) = stoichCoeff;
			}

			// Add species to the "classic" components list
			if(index < this->mPrimSpeciesNum) this->mClassicComponents.insert(pair<string,int>(it->first->name().toStdString(),it->second));
        }

		s2 = stoichMatrix.block(0,this->mRelPrimarySpcIndices.size(),this->mEqReactions.size(),this->mEqReactions.size());
		s2 = - s2.inverse();

		stoichMatrix = s2 * stoichMatrix;

		mComponentMatrix.block(0,0,this->mRelPrimarySpcIndices.size(),this->mRelPrimarySpcIndices.size()) = MatrixXd::Identity(this->mRelPrimarySpcIndices.size(), this->mRelPrimarySpcIndices.size());
		mComponentMatrix.block(0,this->mRelPrimarySpcIndices.size(),this->mRelPrimarySpcIndices.size(),this->mEqReactions.size()) = stoichMatrix.topLeftCorner(this->mEqReactions.size(), this->mRelPrimarySpcIndices.size()).transpose();

		result = this->GetComponentInit(this->mSpeciesIndices);
		
	}

    tempIndices.clear();
	this->mComponentMatrix.resize(0,0);
    return result;
}

MatrixXd 
CLocalChemicalSystemConstPrimary::GetComponentInit(map<CSpecies*,int>& tempIndices)
{
    MatrixXd result;
       
    vector<CSpecies*> species1;
    vector<CSpecies*> species2;

    vector<int> indices;

	species1 = this->mSpecies[eAqueous1];
    species2 = this->mSpecies[eAqueous2];

    //Get useful iterator
    map<CSpecies*, int>::iterator it;

    result = MatrixXd::Zero(this->mComponentMatrix.rows(), species1.size()+species2.size()); 
                
    //for(int i=0; i!=result.rows(); i++)
    //{
        for(int j=0; j!= species1.size(); j++)
        {
			it = tempIndices.find(species1.at(j));
            result.col(it->second) = this->mComponentMatrix.col(it->second);
			// Add species to the "classic" components list
			this->mClassicComponents.insert(pair<string,int>(it->first->name().toStdString(),it->second));
        }

		for(int j=0; j!= species2.size(); j++)
		{
			it = tempIndices.find(species2.at(j));
			result.col(it->second) = this->mComponentMatrix.col(it->second);
			// Add species to the "classic" components list
			this->mClassicComponents.insert(pair<string,int>(it->first->name().toStdString(),it->second));
		}

    //}
 
    return result;
}

