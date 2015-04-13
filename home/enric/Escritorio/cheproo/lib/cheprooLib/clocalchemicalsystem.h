#ifndef CLOCALCHEMICALSYSTEM_H
#define CLOCALCHEMICALSYSTEM_H

#include "ccheproobase.h"

#include "cchemicalcomposition.h"
#include "creaction.h"
#include "cphase.h"

#include <QString>
#include <QFile>
#include <QCoreApplication>
#include <QTextStream>
#include <string>
#include <vector>
#include <map>
#include <math.h>

#define H2O_mol_mass 0.01801534
#define DEFAULT_ERROR 0.00000001

#include <Eigen/Dense>
using namespace Eigen;

using namespace std;

static const char* ERROR_SPCTOPHASE_KIND_NAME = QT_TRANSLATE_NOOP("CLocalChemicalSystem", "Phase kind name not allowed");
static const char* ERROR_SPCTOPHASE_KIND_NUMBER = QT_TRANSLATE_NOOP("CLocalChemicalSystem", "Phase kind number not allowed");
static const char* ERROR_CLASSICCOMPONENT_NOT_FOUND = QT_TRANSLATE_NOOP("CLocalChemicalSystem", "Error! Component not found: ");
static const char* ERROR_LOGDATA_NOT_ALLOWED = QT_TRANSLATE_NOOP("CLocalChemicalSystem", "Error! Value in log scale not accepted for condition: ");
#define ERROR_FILEOPEN "Error: error opening file"
#define ERROR_RISA_INITIALIZATION "Error in defining chemical composition for RISA: no redundant infos defined"

class CGlobalChemicalSystem;

/// \brief Enumeration containing a list of species kind (primary or secondary) and phases they may belong to
enum eSpcToPhaseKind
{
    eAqueous1 = 0,                ///< 0: Aqueous Primary
    eAqueous1nc = 1,              ///< 1: Aqueous Primary without constant activity species (typically, H2O)
    eAqueous1c = 2,               ///< 2: Aqueous Primary constant activity species (typically, H2O)
    eAqueous2 = 3,                ///< 3: Aqueous Secondary
    eMineral1 = 4,                ///< 4: Mineral Primary
    eMineral2 = 5,                ///< 5: Mineral Secondary
    eGas1 = 6,                    ///< 6: Gas Primary
    eGas2 = 7,                    ///< 7: Gas Secondary
    eSurface1 = 8,                ///< 8: Surface
    eSurface2 = 9,                ///< 9: Surface
    eLastSpcToPhaseKind = 10      ///< 10: Exception
};


class CLocalChemicalSystem : public CCheprooBase
{
public:

    int mPrimSpeciesNum; // Number of primary species associated to the LocalChemicalSystem (with constant activity species)

	int mSecSpeciesNum;  // Number of secondary species associated to the LocalChemicalSystem

    CGlobalChemicalSystem* mGlobChemSysPointer;  // Pointer to the GlobalChemicalSystem class

    vector<CReaction*> mEqReactions; // Vector containing the equilibrium reactions applying to the current LocalChemicalSystems

    vector<CReaction*> mKinReactions; // Vector containing the kinetic reactions applying to the current LocalChemicalSystems

    vector<CSpecies*> mPrimRedSpecies; // Vector containing the primary species applying to the current LocalChemicalSystems

    set<string> mCAS; // Vector containing the CAS names to check if a water belongs to this local or not

    map<eSpcToPhaseKind, vector<CSpecies*> > mSpecies; // Vector of pointers to the species of the LocalChemicalSystem, ordered by type ("primary" or "secondary")
                                                       // and phase they belong to (eSpcToPhaseKind)

    map<CSpecies*, int> mSpeciesIndices; // Map containing the indices of the species that belong to the LocalChemicalSystem

    multimap<ePhaseKind, CPhase*> mPhases; // Map containing the phases applying to the current LocalChemicalSystems

    vector<CChemicalComposition*> mRelChemicalCompositions;  // Vector containing the chemical compositions associated to the current LocalChemicalSystem

    map<string,CChemicalComposition*> mChemicalCompositions; 

	map<string,int> mClassicComponents;  // Map with "classic" components names and indices

protected:

	int mPrimSpRedNum;                        // Number of reduced primary species (primary species without constant activity species)

	int mConstActSpeciesNum;                  // Number of constant activity species

	bool mSpeciateWithNR;                     // If true Newton-Raphson is used to speciate c2, if false then Picard is used

    double mMaxRelativeError;                 // Convergence parameter for the speciation (defined by user)

    double mMaxResidual;                      // Convergence parameter for the speciation (defined by user)

    int mMaxIterNum;                          // Maximum number of iteration for Newton-Raphson
 
    vector<int> mRelPrimarySpcIndices;        // Vector of indices of the primary species contained in the LocalChemicalSystem, including constant activity species

    vector<int> mRelPrimSpcIndicesWithoutCA;  // Vector of indices of the primary species contained in the LocalChemicalSystem, without constant activity species

    vector<int> mRelSecondarySpcIndices;      // Vector of indices of the secondary species contained in the LocalChemicalSystem

    vector<int> mRelConstActSpcIndices;       // Vector of indices of the constant activity species contained in the LocalChemicalSystem

    ///	\brief Eigen library dynamic matrix implementation \n
	/// See http://eigen.tuxfamily.org for more 

    MatrixXd mS1Mod;                          // Matrix containing -(1/S2)*S1,nc

    MatrixXd mSecMod;                         // Matrix containing -(1/S2)*S1,c

    MatrixXd mSk;                             // Stoichiometric matrix for kinetic reactions

    VectorXd mLogkMod;                        // Vector containing (1/S2)*logk

	MatrixXd mS2Inverse;                      // Matrix containing -(1/S2), built only if is not -I

    MatrixXd mComponentMatrix;                // Component matrix for all phases

    MatrixXd mComponentMatrix_c;              // Part of the component matrix to evaluate constant activity species concentrations

    MatrixXd mComponentMatrixInit;            // Component matrix for the initialization

	MatrixXd dc2_dc1;                         // Matrix containing dc2/dc1

	MatrixXd jacobian;                        // Matrix containing the jacobian to evaluate c1

	VectorXd residual;						  // Matrix containing the residual to evaluate c1

    MatrixXd mSystemMatrix;                   // System matrix to evaluate c2 and dc2/dc1                             

    MatrixXd mSystemRHS;                      // Vector to the RHS to evaluate dc2/dc1



public:

    void Read(		
        const QDomElement aNode, ///< XML node with attributes
        vector<string>& primSpecies, 
        set<string>& cas
                      ); 

    /// \brief Set the attributes of the local 
    virtual void SetLocal(vector<string>& primSpecies, set<string>& cas){}

    /// \brief Set the convergence parameters from another local
    void SetConvParams(CLocalChemicalSystem* aLocal);

    /// \brief Speciate all the ChemicalCompositions of mRelChemicalCompActualTime
    void Speciate();

    /// \brief Speciate initial water 
    void SpeciateInitialWater(CChemicalComposition* aCurrChemComp);

    /// \brief Speciate from total aqueous concentrations
    void Speciate_uaq(CChemicalComposition* aCurrChemComp);

    /// \brief Speciate from total concentrations in all phases eliminating minerals
    void Speciate_utot(CChemicalComposition* aCurrChemComp);

	/// \brief Speciate from total concentrations in all phases without eliminating minerals, with "classic" definition of components
	void Speciate_utot_classic(CChemicalComposition* aCurrChemComp);

    ///  \brief Calculate u_tot for SIA by means of the classic definition
    VectorXd Calculate_utot_classic_SIA(CChemicalComposition* aCurrChemComp,
									    vector<double> &u_aq);

	///  \brief Calculate constant activity species concentration, typically minerals
	void CalculateCASConcentration(CChemicalComposition* aCurrChemComp);

	///  \brief Calculate u_tot for SIA by eliminating minerals, from u_aq "classic"
	VectorXd Calculate_ured_from_utot_classic(CChemicalComposition* aCurrChemComp,
									          VectorXd u_tot_classic);

    /// \brief Speciate with redundant conditions
    void RISA(CChemicalComposition* aCurrChemComp,
              bool writeHeader = true);

	/// \brief Mix chemical compositions
	void MixChemicalCompositions(map<CChemicalComposition*,double> aInitialWaters, ///< Map of initial waters and corresponding mixing ratios
				                 CChemicalComposition* aMixedWater,                ///< Mixed water 
				                 bool checkSI = false);                           ///< True if mixed water must be speciated

    /// \brief Mix chemical compositions
	void MixWaters(map<CChemicalComposition*,double> &aInitialWaters, ///< Map of initial waters and corresponding mixing ratios
				   CChemicalComposition* aMixedWater,                 ///< Mixed water 
				   bool checkSI = false);                             ///< True if mixed water must be speciated

    /// \brief Apply constraints to speciation to evaluate c1 (charge balance, activity of species, equilibrium with mineral or gases)
    void ApplyConstraints(CChemicalComposition* aCurrChemComp, VectorXd &c1, bool onPrimary, bool ispHAssigned = false, bool isPGasAssigned = false);

    /// \brief Get vector of primary species global indices (method useful for UnitTest)
    vector<int>& GetPrimarySpcIndices();

    /// \brief Get vector of secondary species global indices (method useful for UnitTest)
    vector<int>& GetSecondarySpcIndices();

    /// \brief Get matrix S1Mod (method useful for UnitTest)
    MatrixXd& GetS1Mod();

    /// \brief Get component matrix
    MatrixXd GetComponentMatrix();

    /// \brief Get the term USkrkdt
    vector<double> GetUSkrk(CChemicalComposition* aChemComp);

    /// \brief Get the term dUSkrkdt
    MatrixXd GetdUSkrkdt(const double dt, CChemicalComposition* aChemComp);

    /// \brief Sum ChemicalCompositions of the vector mRelChemicalCompActualTime and put them in another, that by now will be put in the last position of the vector mRelChemicalCompActualTime
    void SumChemicalComp();

    void ComputeAlkalinity(){}

    void ComputeMinArea(){}

    /// \brief Constructor from a XML file.
	CLocalChemicalSystem(
		   QDomElement aNode ///< XML node with attributes
		   );

    /// \brief Converts a string into a time target 
	static eSpcToPhaseKind str2PhaseKind(QString s   ///< String that has to be converted into a time target
                                        );
    /// \brief Converts a time target into a string
	static QString phaseKind2str(eSpcToPhaseKind e ///< phases that the species belongs to that has to be converted into a string
                                   );

    /// \brief Compute secondary concentrations (for Cheproo coupled with Proost)
	void ComputeSecondaries (CChemicalComposition* aCurrChemComp,
							 vector<double> &c1);

    void Compute_dc2_dc1(CChemicalComposition* aCurrChemComp,
                         VectorXd& c1);

    /// \brief Write chemical composition infos on an output file
    void WriteChemCompInfos(CChemicalComposition* aCurrChemComp, 
                            bool chemCompOnSeparateFiles = false   ///< True if the chemical compositions are printed on different output files
                            );

    /// \brief Write RISA infos on an output file
    void WriteRISAInfos(CChemicalComposition* aChemComp, int& iter, double& objFun, VectorXd& grad, double& norm, double& maxRelError, QString name);

    /// \brief Write concentrations of ChemicalCompositions with the same name on an output file
    void WriteChemCompInfos_conc(CChemicalComposition* aCurrChemComp, 
                                 bool isHeader,                         ///< True if the header needs to be printed
                                 bool chemCompOnSeparateFiles = false   ///< True if the chemical compositions are printed on different output files
                                 );

    /// \brief Write concentrations of ChemicalCompositions with the same name on an output file
    void WriteChemCompInfos_gamma(CChemicalComposition* aCurrChemComp, 
                                 bool isHeader,                         ///< True if the header needs to be printed
                                 bool chemCompOnSeparateFiles = false   ///< True if the chemical compositions are printed on different output files
                                 );

    /// \brief Method that returns u_aq as molarity (mol/m3_phase) 
    vector<double> Get_u_aq(CChemicalComposition* aChemComp); 

    /// \brief Method that returns u_gas as molarity (mol/m3_phase) 
    vector<double> Get_u_gas(CChemicalComposition* aChemComp); 

    /// \brief Method that returns u_min as molarity (mol/m3_phase) 
    vector<double> Get_u_min(CChemicalComposition* aChemComp); 

    /// \brief Method that returns u_surf as molarity (mol/m3_phase) 
    vector<double> Get_u_surf(CChemicalComposition* aChemComp);

    /// \brief Method that returns u_aq as molarity (mol/m3_phase) in "classic" definition
    vector<double> Get_u_aq_classic(CChemicalComposition* aChemComp); 

    /// \brief Method that returns u_gas as molarity (mol/m3_phase) in "classic" definition 
    vector<double> Get_u_gas_classic(CChemicalComposition* aChemComp); 

    /// \brief Method that returns u_min as molarity (mol/m3_phase) in "classic" definition
    vector<double> Get_u_min_classic(CChemicalComposition* aChemComp); 

    /// \brief Method that returns u_surf as molarity (mol/m3_phase) in "classic" definition
    vector<double> Get_u_surf_classic(CChemicalComposition* aChemComp);

    CLocalChemicalSystem();
	~CLocalChemicalSystem();


    MatrixXd GetComponentMatrixByPhase(QString aVarType);
    virtual MatrixXd Compute_dc_dc1_ByPhase(QString aVarType, CChemicalComposition* aCurrChemComp);
    virtual MatrixXd Get_dc_dc1_ByPhase(QString aVarType);

protected:

    void ComputeSecondaries(CChemicalComposition* aCurrentChemicalComposition,  // Chemical composition for which secondary concentrations are evaluated
                            VectorXd &c1,                                       // Vector of primary species
                            bool mustCalcDerivs,                                // If true dc2/dc1 is evaluated
                            int iterc1 = 0                                      // Iteration number of iterative process to evaluate c1
                            );

    /// \brief Compute the component matrix for all species except constant activity species
    virtual void ComputeComponentMatrix(){}

    /// \brief Compute the stoichiometric matrix S1* and Sec*.
    virtual void ComputeS1Mod(){}

    /// \brief Compute S2 matrix
    void ComputeS2(MatrixXd s2);

    /// \brief Compute the vector of equilibrium constants k*.
    void ComputeLogkMod();

    // Compute the residual for the N-R to evaluate primary species concentrations
    void ComputeResidualOfMAL(VectorXd &residual,
                              VectorXd &c1,
                              VectorXd &c2,
							  CChemicalComposition* aCurrChemComp);

    // Compute the residual for the N-R to evaluate primary species concentrations
    void ComputeResidualOfComponents(CChemicalComposition* aCurrChemComp,
                                     VectorXd &primaryConc,
                                     MatrixXd &compMatrixPrim,
                                     bool isInizialization = false);

    // Compute the residual for the N-R to evaluate primary species concentrations
    void ComputeResidualOfComponents_allPhases(CChemicalComposition* aCurrChemComp,
                                               VectorXd &primaryConc);

    // Compute the jacobian for the N-R to evaluate primary species concentrations
    void ComputeJacobian(MatrixXd &compMatrixPrim);

    // Compute the jacobian for the N-R to evaluate primary species concentrations with RISA
    void ComputeJacobian(CChemicalComposition* aCurrChemComp, VectorXd &c1, MatrixXd& aXi, MatrixXd& aLambda, MatrixXd& jacobian);

    // Compute the jacobian for the N-R to evaluate primary species concentrations
    void ComputeJacobian_allPhases(CChemicalComposition* aCurrChemComp);

    // Check convergence
    int CheckConvergence(VectorXd &cnew,
                         VectorXd &deltac,
                         VectorXd &residual);

    int CheckConvergence(VectorXd &cnew,
                         VectorXd &deltac,
                         VectorXd &residual,
                         double& relErrMax,
                         int index);

    // Method to check if a vector contains zeros
    bool IsZeroVector(VectorXd &vector);

    // Method to initialize the secondary concentration vector = S1Mod * ln(c1) + ln(k*)
    void InitializeSecondaries(VectorXd &c2,  
							   VectorXd &c1, 
							   CChemicalComposition* currChemComp);

	// Method to obtain the initial guess for a Newton-Raphson method
	void CalculateInitialGuess(CChemicalComposition* aCurrChemComp, 
						 VectorXd &c);

    // Methid to set the matrices for RISA algorithm
    void SetRISAInputs(CChemicalComposition* aChemComp, MatrixXd& aXi, MatrixXd& aLambda, VectorXd& aX, VectorXd& aV_diag);

    // Methid to compute fd for RISA
    void Computefd(CChemicalComposition* aCurrChemComp, VectorXd &c1, VectorXd  &x, MatrixXd& aXi, MatrixXd& aLambda, VectorXd &fd);

};

#endif
