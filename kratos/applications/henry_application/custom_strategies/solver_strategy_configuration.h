#if !defined(KRATOS_SOLVER_STRATEGY_CONFIGURATION )
#define  KRATOS_SOLVER_STRATEGY_CONFIGURATION


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"


namespace Kratos
{


template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class SolverStrategyConfiguration
{
public:

     //typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverTypePointer;
     typedef Scheme< TSparseSpace, TDenseSpace > SchemeType;
     
    /** Constructor.
    */
    SolverStrategyConfiguration(ModelPart& model_part, unsigned int domain_size)
        : mrModelPart(model_part), mDomainSize(domain_size)
    {
    }

    /** Destructor.
    */

    unsigned int GetDomainSize()
    {
        return this->mDomainSize;
    }
    
    ModelPart GetModelPart()
    {
        return this->mrModelPart;
    }

    virtual typename SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer pGetStrategy(const std::string& strategy_name )
    {
        //KRATOS_ERROR(std::logic_error,"accessing to the SolverStrategyConfiguration base class","");
    }
    
    virtual void pSetUpDof(ProcessInfo& rCurrentProcessInfo)
    {
        //KRATOS_ERROR(std::logic_error,"accessing to the SolverStrategyConfiguration base class","");
    }
    /* 
    virtual typename BuilderSolverTypePointer pGetBuildAndSolve(ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR(std::logic_error,"accessing to the SolverStrategyConfiguration base class","");
    }
    */
    /*
    virtual  SchemeType::Pointer pGetScheme(ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR(std::logic_error,"accessing to the SolverStrategyConfiguration base class","");
    }
    */  
protected:

    ModelPart& mrModelPart;
    unsigned int mDomainSize;

private:
    /** Copy constructor.
    */


}; /* Class FractionalStepStrategy */


}  /* namespace Kratos.*/

#endif /* KRATOS_SOLVER_STRATEGY_CONFIGURATION  defined */
