/* *********************************************************
 *
 *   Last Modified by:    $Author: VictorBez $
 *   Date:                $Date: 2015        $
 *   Revision:            $Revision: 1.0     $
 *
 * ***********************************************************/


#if !defined(KRATOS_FRACTIONAL_ITERATIVE_UNITARY_STRATEGY)
#define  KRATOS_FRACTIONAL_ITERATIVE_UNITARY_STRATEGY


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/geometry_utilities.h"
#include <string>
#include <stdio.h>
//Inheritance
#include "solving_strategies/strategies/solving_strategy.h"

//Configuration files of Strategy (to load objects of configuration attributes to deal with Matrix)
//#include "custom_strategies/solver_strategy_configuration.h"
//#include "custom_strategies/fractional_iterative_configuration.h"


//Application
#include "henry_application.h"

using namespace std;

namespace Kratos
{
    
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class FractionalIterativeUnitaryStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    
    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( FractionalIterativeUnitaryStrategy );

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType; 

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    
    
    /*
     * From incompresible_fluid_application/ strategies/ custom_strategies/ fractional_step_streategy.h 
     */
    
    //Doc-> constructor
    FractionalIterativeUnitaryStrategy(
    ModelPart& model_part,
    //SolverStrategyConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>& aSolverStrategyConfiguration,
    //bool aReformDofAtEachIteration = true,
    //double aPressureTol = 0.0001,
    //double aConcentrationTol = 0.0000001,
    //int aMaxIterPressure = 12,
    //int aMaxIterConcentration = 12,
    //unsigned int aTimeOrder = 1,
    unsigned int aDomainSize = 2
    //bool aPredictorCorrector = false,
    //bool aMoveMeshFlag = false, 
    //unsigned int aEchoLevel = 3
    )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, aDomainSize)//, mSolverStrategyConfiguration(aSolverStrategyConfiguration)
    {
        KRATOS_TRY

        //this->mPressureTol = aPressureTol;
        //this->mConcentrationTol = aConcentrationTol;
        //this->mMaxIterPressure = aMaxIterPressure;
        //this->mMaxIterConcentration = aMaxIterConcentration;
        //this->mPredictorOrder = aTimeOrder;
        //this->mTimeOrder = aTimeOrder;
        this->mDomainSize = aDomainSize;
        //this->mPredictorCorrector = aPredictorCorrector;
        //this->mReformDofAtEachIteration = aReformDofAtEachIteration;
        //this->mEchoLevel = aEchoLevel;  
        
        //performs checks to verify the quality of the input
        //this->Check();

        ProcessInfo& CurrentProcessInfo = model_part.GetProcessInfo();
        
        //initialize strategy
        //CurrentProcessInfo[FRACTIONAL_STEP] = 1;
        //this->mPressureStrategy = mSolverStrategyConfiguration.pGetStrategy(std::string("PressureStrategy"));
        //CurrentProcessInfo[FRACTIONAL_STEP] = 2;
        //this->mConcentrationStrategy = mSolverStrategyConfiguration.pGetStrategy(std::string("ConcentrationStrategy"));
        
        this->mStep = 1;

        KRATOS_CATCH("")
    }

    virtual ~FractionalIterativeUnitaryStrategy()
    {
        
    }

    double Solve()
    {
        KRATOS_TRY
        
        for(int elem=0; elem<static_cast<int>(BaseType::GetModelPart().Elements().size()); elem++)
        {
            ModelPart::ElementsContainerType::iterator itElem = BaseType::GetModelPart().ElementsBegin()+elem;
            itElem->CalculateUnitaryFunctions();//obtenemos valores calculados en cada una de las fucniones que hacemos testeo unitario
            //KRATOS_WATCH(elem);
        }
        //this->AssertUnitaryFunctions();

        KRATOS_CATCH("")
    }

   // void AssertUnitaryFunctions();
    
    /*void calculateDensityNodes()
    {
 
   
        const double pDensity_pConc = 700.0;
        const double referenceDensity = 1000.0;
        const double referenceConcentration = 0.0; 
        
        for(int elem=0; elem<static_cast<int>(BaseType::GetModelPart().Elements().size()); elem++)
        {
            ModelPart::ElementsContainerType::iterator itElem = BaseType::GetModelPart().ElementsBegin()+elem;
         
            double densityElement= 0.0;
            unsigned int numberOfNodes = itElem->GetGeometry().PointsNumber();
            double weight = 1.0 / ((double) numberOfNodes);
        
            for (unsigned int node = 0; node < numberOfNodes; node++)
            {    
               const double concNode = itElem->GetGeometry()[node].FastGetSolutionStepValue(CONCENTRATION);
               const double densNode =  (referenceDensity + pDensity_pConc * (concNode-referenceConcentration));

               itElem->GetGeometry()[node].FastGetSolutionStepValue(DENSITY) = densNode;
               densityElement += (weight *densNode);              
            }
            
            itElem->SetValue(DENSITY_ELEM,densityElement);
        }
        
    }*/
    
    
 
    


    void ReadFile(string file)
    {

         FILE *fichero;
         std::string str;
         str = file;
         int nelement;
         double density;
         const char * nombre = str.c_str();
         fichero = fopen( nombre, "r" );
         
         int initialSteps = 1;
         
         

         for(unsigned int elem=0; elem<static_cast<int>(BaseType::GetModelPart().Elements().size()); elem++)
         {
                 ModelPart::ElementsContainerType::iterator itElem = BaseType::GetModelPart().ElementsBegin()+elem;
                 Vector &varVectDens = itElem->GetValue(DENSITY_ELEM_TARGET);
                 varVectDens.resize(initialSteps,true);
                 fscanf(fichero,"%d",&nelement);
                 fscanf(fichero,"%lf",&density);
                 varVectDens[0] = density;

                 KRATOS_WATCH(varVectDens[0])
                 KRATOS_WATCH(itElem->GetValue(DENSITY_ELEM))
         }

                //////////////////////////////////////////////////
    }
    

    
protected:
    


    

private:    

    unsigned int mStep;
    unsigned int mDomainSize; 
       
    //SolverStrategyConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>& mSolverStrategyConfiguration;

    FractionalIterativeUnitaryStrategy(const FractionalIterativeUnitaryStrategy& Other);

}; /* Class FractionalIterativeStrategy */

} /* namespace Kratos.*/

#endif /* KRATOS_FRACTIONAL_ITERATIVE_UNITARY_STRATEGY  defined */

