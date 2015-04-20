/*
==============================================================================
KratosFlowApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
//
//   Project Name:        Kratos
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:32 $
//   Revision:            $Revision: 1.2 $
//
//


// System includes


// External includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp>


// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"

//custom strategies & convergence criterias
#include "custom_strategies/iterative_flow_strategy.h"
#include "custom_strategies/iterative_transport_strategy.h"
//#include "custom_strategies/residualbased_flow_strategy.h"

//configuration files
#include "custom_strategies/solver_strategy_configuration.h"
#include "custom_strategies/iterative_flow_configuration.h"
#include "custom_strategies/iterative_transport_configuration.h"

//linear solvers
#include "linear_solvers/linear_solver.h"



namespace Kratos
{

namespace Python
{
using namespace boost::python;

void  AddCustomStrategiesToPython()
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    //typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

    //********************************************************************
    //********************************************************************
    //

    //Not used
    /*
    class_< ResidualBasedFlowStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
            bases< BaseSolvingStrategyType >,  boost::noncopyable >
            ("ResidualBasedFlowStrategy",
             init<	ModelPart&, LinearSolverType::Pointer,	bool, int, int	>() )
            .def("Clear",&ResidualBasedFlowStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear)
            ;
    */
    
    class_< SolverStrategyConfiguration<SparseSpaceType, LocalSpaceType, LinearSolverType >,
            boost::noncopyable >
            ("SolverStrategyConfiguration", init< ModelPart&, unsigned int>())
            ;

    class_< IterativeFluidConfiguration<SparseSpaceType, LocalSpaceType, LinearSolverType >,
            bases< SolverStrategyConfiguration<SparseSpaceType, LocalSpaceType, LinearSolverType > >,
            boost::noncopyable >
            ("IterativeFluidConfiguration", init< ModelPart&, LinearSolverType::Pointer,
             unsigned int >());
    class_< IterativeTransportConfiguration<SparseSpaceType, LocalSpaceType, LinearSolverType >,
            bases< SolverStrategyConfiguration<SparseSpaceType, LocalSpaceType, LinearSolverType > >,
            boost::noncopyable >
            ("IterativeTransportConfiguration", init< ModelPart&,LinearSolverType::Pointer,
             unsigned int >());
    
    class_< IterativeFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
            bases< BaseSolvingStrategyType >, boost::noncopyable >
            ("IterativeFluidStrategy",
             init < ModelPart&,
             SolverStrategyConfiguration<SparseSpaceType, LocalSpaceType, LinearSolverType >&,
             bool,
             double, 
             int,
             unsigned int, unsigned int,
             bool, bool, unsigned int
             >())
            .def("iterationPressure", &IterativeFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::iterationPressure)
            .def("StorageOldSPressIteration", &IterativeFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::StorageOldSPressIteration)
            .def("isConverged", &IterativeFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::isConverged)
            .def("CalculateFlowBalances", &IterativeFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::CalculateFlowBalances)
            .def("SetEchoLevel", &IterativeFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetEchoLevel)
            .def("Clear", &IterativeFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear)
            .def("AssignInitialStepValues", &IterativeFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::AssignInitialStepValues)
            .def("PredictSV", &IterativeFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::PredictSV)
            .def("InitializeFractionalStep", &IterativeFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::InitializeFractionalStep)
            
        ;
   
      class_< IterativeTransportStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
            bases< BaseSolvingStrategyType >, boost::noncopyable >
            ("IterativeTransportStrategy",
             init < ModelPart&,
             SolverStrategyConfiguration<SparseSpaceType, LocalSpaceType, LinearSolverType >&,
             bool,
             double,
             int,
             unsigned int, unsigned int,
             bool, bool, unsigned int
             >())
            .def("iterationConcentration", &IterativeTransportStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::iterationConcentration)
            .def("StorageOldSConcIteration", &IterativeTransportStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::StorageOldSConcIteration)
            .def("isConverged", &IterativeTransportStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::isConverged)
            .def("calculateDensityNodes", &IterativeTransportStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::calculateDensityNodes)
            .def("SetEchoLevel", &IterativeTransportStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetEchoLevel)
            .def("Clear", &IterativeTransportStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear)
            .def("AssignInitialStepValues", &IterativeTransportStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::AssignInitialStepValues)
            .def("PredictSV", &IterativeTransportStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::PredictSV)
            .def("InitializeFractionalStep", &IterativeTransportStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::InitializeFractionalStep)
            .def("CalculateDarcyFlowPressure", &IterativeTransportStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::CalculateDarcyFlowPressure)
        ;
        
}

}  // namespace Python.

} // Namespace Kratos

