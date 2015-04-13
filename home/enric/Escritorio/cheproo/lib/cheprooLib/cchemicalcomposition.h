#ifndef CCHEMICALCOMPOSITION_H
#define CCHEMICALCOMPOSITION_H

#include "ccheproobase.h"

#include <Eigen/Dense>
using namespace Eigen;

#include <string>
#include <QString>
#include <vector>
#include <map>
#include <set>
using namespace std;

static const char* ERROR_LOGSCALE_NOT_ALLOWED = QT_TRANSLATE_NOOP("CChemicalComposition", "Error! Value in log scale not accepted for charge balance condition.");

/// \brief Enumeration containing a list of phases kind 
enum eCheprooVariableType
{
    eAqComponent, 
    eMinComponent, 
    eSurfComponent, 
    eGasComponent,
    eLiqDensity, 
    eGasDensity,
    eLiqViscosity, 
    eGasViscosity, 
    eKineticSinkSource,
    ePrimaryConc,
    eSecondaryConc 
};

class CLocalChemicalSystem;

class CChemicalComposition : public CCheprooBase
{
public:

    CLocalChemicalSystem* mLocChemSysPointer; // Pointer to the LocalChemicalSystem that the ChemicalComposition belongs to
    
    vector<double> mConcentration;  // Vector containing the concentrations of the species ***convert to VectorXd sooner or later, better sooner****

    vector<double> mActivityCoeff;  // Vector containing the activity coefficients

    map<string, double> mVolumetricFractions; // Map containing the volumetric fractions for all phases present in the chemical composition (key=name of the phase)

    map<string, double> mSatIndices;  // Map containing the saturation indices of the minerals ***still to be evaluated***

    double mNodesVolume;

    double mTemp;

    double mGasDensity;

    double mGasViscosity;

    double mLiquidDensity;

    double mLiquidViscosity;

    double mLiquidPressure;

    // Think about those attributes, especially constrValue, if they need to be public or if they can be deleted at some point
    
    vector<string> iCon;
    VectorXd constrValue;
    vector<string> constraints;
	vector<string> refSpecies;

    vector<double> errors;
    vector<string> iCon_uncertain;
    vector<bool> isLogData;
    vector<bool> isRelError;
    VectorXd constrValue_uncertain;
    vector<string> constraints_uncertain;
	vector<string> refSpecies_uncertain;

    double measured_value; // Variable for RISA verification, delete later!

private:

    map<string,double> cGuess;

    
public:
    /// \brief default Constructor.
    CChemicalComposition();

    /// \brief default Constructor.
    CChemicalComposition(const CChemicalComposition& original )
    {
        throw("Error: not implemented");
    };

    /// \brief Constructor from a XML file.
	CChemicalComposition(
		   QDomElement aNode ///< XML node with attributes
           );

	/// \brief Constructor for RISA simulations
    CChemicalComposition( 
           CChemicalComposition* aChemComp,
           double& p, 
		   double& q, 
		   double& r, 
		   double& s,
           double& t,
           bool isPerturbedValue = false
           );


	virtual ~CChemicalComposition();

    /// \brief Read from a XML file.
    void Read(		
        const QDomElement aNode,    ///< XML node with attributes
        vector<string> &primSpecies, ///< primary species names
        set<string> &cas     ///< mineral names
                      );

    /// \brief Set the chemical composition from a local chemical system attributes
    void SetChemComp();

    /// \brief Get the concentration values of a sub-set of species, given their indices
    void GetSetOfConcVector(VectorXd &concVector, vector<int> &speciesIndices);

    /// \brief Get the concentration values of a sub-set of species or components, given the type of variable that are necessary (useful for coupling with CProost)
    vector<double> GetVector(eCheprooVariableType aVar, bool useClassicComponents = true);

    /// \brief Get the concentration values of a sub-set of species or components, given the type of variable that are necessary (useful for coupling with CProost)
    double GetScalar(eCheprooVariableType aVar);

    // Get the matrix containing dc2/dc1 for this chemical composition (useful for coupling with CProost)
    MatrixXd Get_dc_dc1_byPhase(QString aVarType);

    // Get the component matrix
    MatrixXd GetUByPhase(QString aVarType);

	/// \brief Get the natural logarithm of concentration values of a sub-set of species, given their indices
    void GetSetOfLnConcVector(VectorXd &concVector, vector<int> &speciesIndices);

    /// \brief Get the activity coefficient values of a sub-set of species, given their indices
    void GetSetOfActivityCoeffVector(VectorXd &actCoeffVector, vector<int> &speciesIndices);

    /// \brief Update the concentration values of a sub-set of species, given their indices
    void UpdateConcentrations(VectorXd &concVector, vector<int> &speciesIndices);

    /// \brief Update the concentration values of a sub-set of species, given their indices
    void UpdateConcentration(double &concValue, int &speciesIndices);

    /// \brief Update the concentration values of a sub-set of species, given their indices and their natural logarithm values
    void UpdateConcentrationsFromLn(VectorXd &concVector, vector<int> &speciesIndices);

    /// \brief Update the activity coeffieients values of a sub-set of species, given their indices
    void UpdateActivityCoeffs(vector<double> &actCoeffsVector, vector<int> &speciesIndices);

    /// \brief Update the concentration values of a sub-set of species, given their indices
    void UpdateVolumFraction(double value, string phaseName);

	/// \brief Given a set of c1, calls the speciation function of the CLocalChemicalSystem to evaluate c2
	void Set_c1(vector<double> &c1);

    /// \brief Copy a ChemicalComposition into this instance
    void Copy(CChemicalComposition* aChemComp);

    /// \brief Copy a ChemicalComposition into this instance for Reactive Transport
    void CopyForRT(CChemicalComposition* aChemComp);

    /// \brief Perform the chemical step for SIA (i.e.,  speciate with u_aq that is received from Proost)
    void ChemicalStepSIA(vector<double> &u_aq, double &f_aq, double &f_s, double &f_tot, bool useClassicComponents = false);


};


static eCheprooVariableType Vartype(QString in)
{



    if (in == "GasComponents")
    {
        return eGasComponent;
    }
    else if (in == "SurfaceComponents")
    {
        return eSurfComponent;
    }
    else if (in == "MineralComponents")
    {
        return eMinComponent;
    }
    else if (in == "AqueousComponents")
    {
        return eAqComponent;
    }
    else if (in == "LiquidDensity")
    {
        return eLiqDensity;
    }
    else if (in == "LiquidViscosity")
    {
        return eLiqViscosity;
    }
     else if (in == "GasdDensity")
    {
        return eGasDensity;
    }
    else if (in == "GasViscosity")
    {
        return eGasViscosity;
    }        
    else if (in == "PrimaryConcentrations")
    {
        return ePrimaryConc;
    }
    else if (in == "SecondaryConcentrations")
    {
        return eSecondaryConc;
    }
    else if (in == "KineticSinkSources")
    {
        return eKineticSinkSource;
    }

    // otherwise we're doomed
    GeneralException e("error: input string did not match any Cheproo variable");
    throw(e);

}
    static bool  IsVectorial(eCheprooVariableType variable)
    {
        switch (variable)
        {
            case  eLiqDensity:
            case eGasDensity:
            case eLiqViscosity:
            case eGasViscosity:
                return false;
                break;
            case eAqComponent:
            case eMinComponent:
            case eSurfComponent:
            case eGasComponent:
            case eKineticSinkSource:
            case ePrimaryConc:
            case eSecondaryConc:
                return true;
                break;
             default:
                GeneralException e("error: input unknown cheprooVariableType");
                throw(e);
          }
    }

#endif
