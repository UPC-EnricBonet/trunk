#include "clocalchemicalsystemsaaltink.h"
#include "cglobalchemicalsystem.h"

CLocalChemicalSystemSaaltink::CLocalChemicalSystemSaaltink(void)
{
    this->mClassName = "CLocalChemicalSystemSaaltink";
}

CLocalChemicalSystemSaaltink::CLocalChemicalSystemSaaltink(QDomElement aNode)
: CLocalChemicalSystem(aNode)
{
    this->mClassName = "CLocalChemicalSystemSaaltink";
}

CLocalChemicalSystemSaaltink::~CLocalChemicalSystemSaaltink(void)
{
}

void
CLocalChemicalSystemSaaltink::SetLocal(vector<string>& primSpecies, set<string>& cas)
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
    for (int j=0; j!=primSpecies.size(); j++)
	{
		// Retrieve the species name
		const string aSpeciesName = primSpecies.at(j);

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
        this->mSpeciesIndices.insert(pair<CSpecies*, int>(mGlobChemSysPointer->mAllSpecies[aSpeciesName],j));

        // Update number of primary species
        this->mPrimSpeciesNum = j+1;

        // Fill vectors with primary species indices, with and without constant activity species
        this->mRelPrimarySpcIndices.push_back(j);

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

    int i = this->mPrimSpeciesNum;

    map<string,CSpecies*>::iterator it;
    for(it = this->mGlobChemSysPointer->mAllSpecies.begin() ; it != this->mGlobChemSysPointer->mAllSpecies.end() ; it++)
    {   
        // Check if the current species is already present in a reaction
        bool isNotPresent = true;
        map<string,CReaction*>::iterator it1;
        map<CSpecies*,int>::iterator it2;
        for(it1 = this->mGlobChemSysPointer->mAllReactions.begin(); it1!=this->mGlobChemSysPointer->mAllReactions.end(); it1++)
        {
            int m=0;
            CReaction* currReaction = it1->second;
            m = currReaction->mRelSpeciesVsStoichCoeffs.count(it->second);
            if(m!=0 && it1->first==it->first) isNotPresent = false;

        }
        // Check if I can find this species in the vector of primary species names
        myitFind = primSpeciesAux.find(it->first); 
        if(myitFind==primSpeciesAux.end() && isNotPresent) // Read the reaction if is not present
        {
            QString aReactionName = QString::fromStdString(it->first);
            this->mGlobChemSysPointer->CreateReaction(aReactionName);
            this->mEqReactions.push_back(this->mGlobChemSysPointer->mAllReactions[it->first]);
        }

        // Add the reaction to the corresponding vector
        //if(this->mGlobChemSysPointer->mAllReactions[it->first]->mIsKinetic)
        //{
            //this->mEqReactions.push_back(this->mGlobChemSysPointer->mAllReactions[it->first]);
        //}
        //else
        //{
        //    this->mEqReactions.push_back(this->mGlobChemSysPointer->mAllReactions[it->first]);
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

                this->mRelSecondarySpcIndices.push_back(i);

				// Continue filling the map of species indices
				this->mSpeciesIndices.insert(pair<CSpecies*, int>(it->second,i));

				i++;
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
					// Continue filling the map of species indices
					this->mSpeciesIndices.insert(pair<CSpecies*, int>(it->second,i));

                    this->mRelSecondarySpcIndices.push_back(i);

					i++;

                    auxVector7.push_back(it->second);
                }

                

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
					// Continue filling the map of species indices
					this->mSpeciesIndices.insert(pair<CSpecies*, int>(it->second,i));

				    this->mRelSecondarySpcIndices.push_back(i);

				    i++;

                    auxVector8.push_back(it->second);
                }

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
		for(int i=0; i!=constActSpcToLocChemSys.size(); i++)
		{
			auxVector1.pop_back();
			auxVector9.pop_back();

		}
	}

	// Number of species in mSpeciesIndices until now (useful to update the indices)
	int spcNum = mSpeciesIndices.size();
	for(int i=0; i!=constActSpcToLocChemSys.size(); i++)
	{
		// Continue filling the map of species indices
		this->mSpeciesIndices.insert(pair<CSpecies*, int>(constActSpcToLocChemSys[i],i + spcNum));

        // Update indices of constant activity species
        this->mRelConstActSpcIndices.push_back(i + spcNum);

        // Update primary species indices (reduced)
        this->mRelPrimSpcIndicesWithoutCA.pop_back();
		
        // Update secondary indices
        this->mRelSecondarySpcIndices.push_back(i + spcNum);

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
                // do nothing
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
 CLocalChemicalSystemSaaltink::ComputeS1Mod()
{
    // S1* initialization
    this->mS1Mod = MatrixXd::Zero(mEqReactions.size(), this->mPrimSpeciesNum);
    
    if(this->mS2Inverse.size()==0) // If S2 hasn't been computed yet then do it now
    {
	    // Declare and allocate S2
	    this->mS2Inverse = MatrixXd::Zero(mEqReactions.size(), mEqReactions.size());

        // Compute S2
        this->ComputeS2(mS2Inverse);
    }

	// Check if S2 is not -I
	if(this->mS2Inverse != MatrixXd::Identity(mEqReactions.size(), mEqReactions.size()) * -1.0)
	{
		// Then build mS1mod = -(1/S2)*S1

		// Declare and allocate S1 for aqueous species
		MatrixXd aS1ea = MatrixXd::Zero(mEqReactions.size(), this->mPrimSpeciesNum);

        // Loop over primary species indices (the species must be in order of index in S2)
        for(int i=0; i!=this->mPrimSpeciesNum; i++)
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
			            aS1ea(j,i) = stoichCoeff;

                    }
                 }
             }
         } 

		// Evaluate the inverse of S2
		this->mS2Inverse = this->mS2Inverse.inverse();

		//Perform the multiplication to evaluate mS1Mod
		this->mS1Mod = - this->mS2Inverse * aS1ea;

		aS1ea.resize(0,0);
	}
	else
	{
        // Loop over primary species indices (the species must be in order of index in S2)
        for(int i=0; i!=this->mPrimSpeciesNum; i++)
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
	}
}

void 
CLocalChemicalSystemSaaltink::ComputeComponentMatrix()
{
    int aSpeciesNum = this->mPrimSpeciesNum + mEqReactions.size();

    // Evaluate the component matrix: allocate and initialize 
    this->mComponentMatrixInit = MatrixXd::Zero(this->mPrimSpeciesNum, aSpeciesNum);
    this->mComponentMatrix = MatrixXd::Zero(this->mPrimSpeciesNum, aSpeciesNum);
    this->mComponentMatrix_c = MatrixXd::Zero(this->mConstActSpeciesNum, aSpeciesNum);

    // Assign values, first in all phases
    for(int i=0; i!= this->mPrimSpeciesNum; i++)
    {
        // Diagonal terms equal to 1
        this->mComponentMatrix(i,i) = 1.0;

        for(int j=0; j!=mEqReactions.size(); j++)
        {
            this->mComponentMatrix(i, j + this->mPrimSpeciesNum) = this->mS1Mod(j,i);
        }
    }

    // Then take only the aqueous part for initialization
    this->mComponentMatrixInit = this->GetComponentMatrixByPhase("component_aq");

	// Fill the map with "classic" components
	map<CSpecies*, int>::iterator it;
	for(it=this->mSpeciesIndices.begin(); it!=this->mSpeciesIndices.end(); it++)
	{
		this->mClassicComponents.insert(pair<string,int>(it->first->name().toStdString(),it->second));
	}


    // Evaluate the component matrix for constant activity species
    this->mComponentMatrix_c = this->mComponentMatrix.bottomRows(this->mConstActSpeciesNum);

    // Evaluate the aqueous component matrix **to be verified**
    //this->mComponentMatrix_aq = this->mComponentMatrixInit.topLeftCorner(this->mPrimSpeciesNum,aSpeciesNum-this->mConstActSpeciesNum);

	// Evaluate the elimination matrix
	this->ComputeEliminationMatrix();

	// Multiply the component matrix and the elimination matrix: U = E*U
	MatrixXd aEU = MatrixXd::Zero(this->mEliminationMatrix.rows(),this->mComponentMatrix.cols());
	aEU = this->mEliminationMatrix * this->mComponentMatrix;

	// Put EU in the component matrix
	this->mComponentMatrix= MatrixXd::Zero(aEU.rows(),aEU.cols());
	this->mComponentMatrix = aEU;


    // I need to update indices of species and mS1mod by now I do it here (to change when structure of the code will change)
    if(this->mConstActSpeciesNum!=0) this->UpdateLocal();

	aEU.resize(0,0);


}

void 
CLocalChemicalSystemSaaltink::ComputeEliminationMatrix()
{
	// First take the part of U relative to constant activity species 
    MatrixXd aUc = this->mComponentMatrix.rightCols(this->mConstActSpeciesNum);

	// Declare and evaluate Up1 and Up2 (parts of U relative to constant activity species)
	MatrixXd aUc1 = aUc.topRows(this->mPrimSpeciesNum - this->mConstActSpeciesNum);
	MatrixXd aUc2 = aUc.bottomRows(this->mConstActSpeciesNum);

	// Evaluate the inverse of Up2
	aUc2 = aUc2.inverse();

	// Multiply -Up1*(1/Up2) and put it in U*
	MatrixXd aUstar = - aUc1 * aUc2;

	// Build the elimination matrix 
	this->mEliminationMatrix = MatrixXd::Zero(this->mPrimSpRedNum, this->mPrimSpRedNum + this->mConstActSpeciesNum);

    // Assign values
    for(int i=0; i!=this->mPrimSpRedNum; i++)
    {
        // Diagonal terms equal to 1
        this->mEliminationMatrix(i,i) = 1.0;

        for(int j=0; j!=aUstar.cols(); j++)
        {
            this->mEliminationMatrix(i,this->mPrimSpRedNum + j) = aUstar(i,j);
        }
    }
	
	aUc1.resize(0,0);
	aUc2.resize(0,0);
	aUc.resize(0,0);
    aUstar.resize(0,0);
}

void 
CLocalChemicalSystemSaaltink::UpdateLocal()
{
    // Update secondary species indices:

    vector<int>::iterator it;
    int index = this->mRelSecondarySpcIndices.at(0) - 1;
    for(int i=0; i!=this->mConstActSpeciesNum; i++)
    {
        // First delete indices relative to constant activity species
        this->mRelSecondarySpcIndices.pop_back();

        // Then add indices of primary species that become secondary
        it = this->mRelSecondarySpcIndices.begin();
        this->mRelSecondarySpcIndices.insert(it, index);
        index--;
    }

    // Update also mS1mod (I get rid of the colums relative to constant activity species

    this->mS1Mod.conservativeResize(this->mS1Mod.rows(),this->mS1Mod.cols()-this->mConstActSpeciesNum);
    this->mS1Mod = this->mComponentMatrix.block(0,this->mPrimSpRedNum,this->mComponentMatrix.rows(),this->mSecSpeciesNum).transpose();

    this->mS2Inverse.resize(0,0);

    this->ComputeLogkMod();
    
}