#if !defined(KRATOS_FRACTIONAL_ITERATIVE_CONFIGURATION )
#define  KRATOS_FRACTIONAL_ITERATIVE_CONFIGURATION


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
//#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"

//#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_slip.h"

#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"


namespace Kratos
{

template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class FractionalIterativeConfiguration : public SolverStrategyConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
                         typedef typename BaseType::DofsArrayType DofsArrayType; 
    typedef Scheme< TSparseSpace, TDenseSpace > SchemeType;
    typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverTypePointer;
    
    /*
     * From incompresible_fluid_application/ strategies/ custom_strategies/ fractional_iterative_configuration.h 
     */
    
    FractionalIterativeConfiguration(ModelPart& model_part,
                                typename TLinearSolver::Pointer pNewPressureLinearSolver,
                                typename TLinearSolver::Pointer pNewConcentrationLinearSolver,
                                unsigned int mDomainSize
                               )
        : SolverStrategyConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, mDomainSize)
    {

        const bool CalculateReactions = false;
        const bool CalculateNormDxFlag = true;
        const bool ReformDofAtEachIteration = true;
        
        //FLOW SYSTEM
                
        this->mSchemePressure = typename SchemeType::Pointer
                                               (new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace > ());

        this->mPressureBuild = BuilderSolverTypePointer(
                                   new ResidualBasedBlockBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver > (pNewPressureLinearSolver) ); 
                                    //new ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace, TDenseSpace, TLinearSolver, Variable<double> >(pNewPressureLinearSolver, PRESSURE));
                                   //////
        
        this->mPressureStrategy = typename BaseType::Pointer( new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver > (model_part,mSchemePressure,pNewPressureLinearSolver,mPressureBuild,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
        this->mPressureStrategy->SetEchoLevel(2);
        
        //TRANSPORT SYSTEM
        
        this->mSchemeConcentration = typename SchemeType::Pointer
                                               (new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace > ());
        
        this->mConcentrationBuild = BuilderSolverTypePointer(
                                    new ResidualBasedBlockBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver > (pNewConcentrationLinearSolver) ); 
                                    //new ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace, TDenseSpace, TLinearSolver, Variable<double> >(pNewConcentrationLinearSolver, CONCENTRATION));
                                   //////
        
        this->mConcentrationStrategy = typename BaseType::Pointer( new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver > (model_part,mSchemeConcentration,pNewConcentrationLinearSolver,mConcentrationBuild,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
        this->mConcentrationStrategy->SetEchoLevel(2);    

    }


    typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer pGetStrategy(const std::string& strategy_name)
    {
        KRATOS_TRY

        if (strategy_name == std::string("PressureStrategy"))
            return mPressureStrategy;
        else if (strategy_name == std::string("ConcentrationStrategy"))
            return mConcentrationStrategy;
        else
            //KRATOS_ERROR(std::invalid_argument, "trying to get an inexisting strategy", "");

        //KRATOS_CATCH("")
    }

    void pSetUpDof(ProcessInfo& rCurrentProcessInfo)
    {
                KRATOS_TRY
        
        switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
	{
	   case 1:
		{
		    this->mPressureBuild->SetUpDofSet(mSchemePressure, mPressureStrategy->GetModelPart());
                    //DofsArrayType& rDofSetPressure =this->mPressureBuild->GetDofSet();
		    break;
		}
	    case 2:
		{
		    this->mConcentrationBuild->SetUpDofSet(mSchemeConcentration, mConcentrationStrategy->GetModelPart());
	            //DofsArrayType& rDofSetConcentration =this->mConcentrationBuild->GetDofSet();
                    break;
		}
	    default:
		{
			//KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
		}
        }

        //KRATOS_CATCH("")
                
    }
    /*
    typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer pGetBuildAndSolve(ProcessInfo& rCurrentProcessInfo)
    {
                        KRATOS_TRY
        
        switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
	{
	   case 1:
		{
                    //FlowSystem
		    return this->mPressureBuild;
		    break;
		}
	    case 2:
		{
		    return this->mConcentrationBuild;
	            break;
		}
	    default:
		{
			KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
		}
        }

        KRATOS_CATCH("")
        
    }
    */
    /*
    typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer pGetScheme(ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        
        switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
	{
	   case 1:
		{
                    //FlowSystem
		    return this->mSchemePressure;
		    break;
		}
	    case 2:
		{
		    return this->mSchemeConcentration;
	            break;
		}
	    default:
		{
			KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
		}
        }

        KRATOS_CATCH("")
    }
    */
protected:


private:

    typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer mPressureStrategy;
    typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer mConcentrationStrategy;

    BuilderSolverTypePointer mPressureBuild; 
    BuilderSolverTypePointer mConcentrationBuild; 
    
    typename SchemeType::Pointer mSchemePressure;
    typename SchemeType::Pointer mSchemeConcentration;

}; /* Class FractionalStepStrategy */

} /* namespace Kratos.*/

#endif /* KRATOS_FRACTIONAL_ITERATIVE_CONFIGURATION  defined */
