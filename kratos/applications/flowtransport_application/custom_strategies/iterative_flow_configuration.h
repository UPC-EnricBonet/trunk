#if !defined(KRATOS_ITERATIVE_FLOW_CONFIGURATION )
#define  KRATOS_ITERATIVE_FLOW_CONFIGURATION


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
class IterativeFluidConfiguration : public SolverStrategyConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
                         typedef typename BaseType::DofsArrayType DofsArrayType; 
    typedef Scheme< TSparseSpace, TDenseSpace > SchemeType;
    typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverTypePointer;
    
    /*
     * From incompresible_fluid_application/ strategies/ custom_strategies/ fractional_iterative_configuration.h 
     */
    
    IterativeFluidConfiguration(ModelPart& model_part,
                                typename TLinearSolver::Pointer pNewPressureLinearSolver,
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
        
        this->mPressureStrategy = typename BaseType::Pointer( new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver > (model_part,mSchemePressure,pNewPressureLinearSolver,mPressureBuild,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
        this->mPressureStrategy->SetEchoLevel(2);
        
       
    }


    typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer pGetStrategy(const std::string& strategy_name)
    {
        KRATOS_TRY

        if (strategy_name == std::string("PressureStrategy"))
            return mPressureStrategy;
        else
            KRATOS_ERROR(std::invalid_argument, "trying to get an inexisting strategy", "");

        KRATOS_CATCH("")
    }

    void pSetUpDof(ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        
		this->mPressureBuild->SetUpDofSet(mSchemePressure, mPressureStrategy->GetModelPart());


        KRATOS_CATCH("")
                
    }

protected:


private:

    typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer mPressureStrategy;

    BuilderSolverTypePointer mPressureBuild; 
    
    typename SchemeType::Pointer mSchemePressure;

}; /* Class FractionalStepStrategy */

} /* namespace Kratos.*/

#endif /* KRATOS_ITERATIVE_FLOW_CONFIGURATION  defined */
