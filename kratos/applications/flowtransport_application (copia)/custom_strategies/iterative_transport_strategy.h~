/* *********************************************************
 *
 *   Last Modified by:    $Author: VictorBez $
 *   Date:                $Date: 2015        $
 *   Revision:            $Revision: 1.0     $
 *
 * ***********************************************************/


#if !defined(KRATOS_ITERATIVE_TRANSPORT_STRATEGY)
#define  KRATOS_ITERATIVE_TRANSPORT_STRATEGY


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/geometry_utilities.h"

//Inheritance
#include "solving_strategies/strategies/solving_strategy.h"

//Configuration files of Strategy (to load objects of configuration attributes to deal with Matrix)
#include "custom_strategies/solver_strategy_configuration.h"
#include "custom_strategies/fractional_iterative_configuration.h"


//Application
#include "flowtransport_application.h"


namespace Kratos
{
    
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class IterativeTransportStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    
    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( IterativeTransportStrategy );

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
    IterativeTransportStrategy(
    ModelPart& model_part,
    SolverStrategyConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>& aSolverStrategyConfiguration,
    bool aReformDofAtEachIteration = true,
    double aConcentrationTol = 0.0000001,
    int aMaxIterConcentration = 12,
    unsigned int aTimeOrder = 1,
    unsigned int aDomainSize = 2,
    bool aPredictorCorrector = false,
    bool aMoveMeshFlag = false, 
    unsigned int aEchoLevel = 3
    )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, aMoveMeshFlag), mSolverStrategyConfiguration(aSolverStrategyConfiguration)
    {
        KRATOS_TRY

        this->mConcentrationTol = aConcentrationTol;
        this->mMaxIterConcentration = aMaxIterConcentration;
        this->mPredictorOrder = aTimeOrder;
        this->mTimeOrder = aTimeOrder;
        this->mDomainSize = aDomainSize;
        this->mPredictorCorrector = aPredictorCorrector;
        this->mReformDofAtEachIteration = aReformDofAtEachIteration;
        this->mEchoLevel = aEchoLevel;  
        
        //performs checks to verify the quality of the input
        this->Check();

        ProcessInfo& CurrentProcessInfo = model_part.GetProcessInfo();
        
        //initialize strategy
        this->mConcentrationStrategy = mSolverStrategyConfiguration.pGetStrategy(std::string("ConcentrationStrategy"));
        
        this->mStep = 1;

        KRATOS_CATCH("")
    }

    virtual ~IterativeTransportStrategy()
    {
        
    }

    double Solve()
    {
        KRATOS_TRY
        
        Timer::Start("solve");
        
        //Initial Guess, if not is by default method
        //Predict();
        
        //mCalculateNormDxFlag = false;
        //double Dp_norm = 0.0;
        
        //Doc -> assign the correct fractional step coefficients (BDF_COEFFICIENTS..)
        //by now, empty
        InitializeFractionalStep(this->mStep, this->mTimeOrder);
        
        //Set to zero Old iteration pressure
        PredictSV(this->mStep, this->mPredictorOrder);

        //Doc -> Save old iteration to compare in convergence with the new pressure iteration. Save always that repit iterations.
        //Assign Velocity To Fract Step Velocity and Node Area to Zero
        //AssignInitialStepValues();

        // Loop convergence
        double iteration = 0; 
        
        //solve for fractional step systems
        bool isGlobalConverged = false;
        while (!isGlobalConverged)
        {
            //perform one iteration over Pressure and concentration
            //if (rank == 0)
            ////////////////////////////////////////////////////////////////////
            std::cout << "New iteration. Iteration = " << iteration  << std::endl;
            std::ofstream myfile;
            myfile.open ("MatrixSystem.txt", std::ios::app);
            //myfile << "New iteration. Iteration = " << iteration  << std::endl;
            //////////////////////////////////////////////////////////////////////
            
            this->CalculateDarcyFlowPressure();

            this->calculateDensityNodes();
            
            if(this->mEchoLevel == 3)
                std::cout << "TRANSPORT EQUATION SOLVE: "   << std::endl;
            
            if(iteration != 0)
                this->StorageOldSConcIteration();
            
            double concentrationNormDx = iterationConcentration();

            //this->mConcentrationStrategy->Clear();
            
            std::cout << "Concentration norm = " << concentrationNormDx  << std::endl;
            
            iteration++;
            
            /////////////////////////////////////////////////
            std::cout << "End of iteration = " << iteration  << std::endl;
            myfile.open ("MatrixSystemTransport.txt", std::ios::app);
            myfile << "End of iteration = " << iteration  << std::endl;  
            myfile.close();
            /////////////////////////////////////////////////
            
            if( isConverged() )
            {
                isGlobalConverged = true;
                std::cout << "CONVERGENCE  ACHIEVED for time step number: " << this->mStep << std::endl;
                std::cout << "It reachs at " << iteration << " iteration" << std::endl;
                
                /////////////////////////////////////////////////
                myfile.open ("MatrixSystemTransport.txt", std::ios::app);
                myfile << "CONVERGENCE  ACHIEVED for time step number: " << this->mStep << std::endl;
                myfile << "It reachs at " << iteration << " iteration" << std::endl << std::endl;
                myfile.precision(10);
                for(int elem=0; elem<static_cast<int>(BaseType::GetModelPart().Elements().size()); elem++)
                {
                    ModelPart::ElementsContainerType::iterator itElem = BaseType::GetModelPart().ElementsBegin()+elem;
                    
                    const double densityElem = itElem->GetValue(DENSITY_ELEM);
                    myfile << "Elem: " <<  itElem->Id() << std::endl;
                    myfile << "DENSITY_ELEM   " << densityElem << std::endl;
                    const double darcyFlowX = itElem->GetValue(DARCY_FLOW_X);
                    const double darcyFlowY = itElem->GetValue(DARCY_FLOW_Y);
                    myfile << "Darcy_Flow_X  "<< darcyFlowX << "   Darcy_Flow_Y  "<< darcyFlowY << std::endl << std::endl;

                    for (unsigned int node = 0; node < itElem->GetGeometry().PointsNumber(); node++)
                    { 
                    
                    const double concNode = itElem->GetGeometry()[node].FastGetSolutionStepValue(CONCENTRATION);
                    const double densNode = itElem->GetGeometry()[node].FastGetSolutionStepValue(DENSITY);
                    
                    myfile << "Node: " <<  itElem->GetGeometry()[node].Id()  << "   CONCENTRATION   " << concNode << "   DENSITY_NODE   " << densNode << std::endl;
                    
                    }
                    
                    myfile << std::endl;
                }
                myfile << std::endl << std::endl;
           
                myfile.close();
                //////////////////////////////////////////////////
            }
            
            if(!isGlobalConverged && (iteration > this->mMaxIterConcentration))
            {
                std::cout << "ATTENTION: convergence NOT achieved iteration  " << iteration << " reached" << std::endl;
                
                //////////////////////////////////////////////////
                myfile.open ("MatrixSystem.txt", std::ios::app);
                myfile << "ATTENTION: convergence NOT achieved iteration  " << iteration << " reached" << std::endl;
                myfile.close();
                //////////////////////////////////////////////////
                        
                break;
            }
            
        }

        if (this->mReformDofAtEachIteration == true)
            this->Clear();

        this->mStep += 1;
        Timer::Stop("solve");
        
        return 0.0;
        
        KRATOS_CATCH("")
    }

    
    ///////////////////////////////////////
    
    void CalculateDarcyFlowPressure()
    {
        ProcessInfo& CurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
        
        unsigned int isBuoyancy = CurrentProcessInfo[IS_BUOYANCY];
        array_1d<double,2> gravity= CurrentProcessInfo[GRAVITY];
        
        for(int elem=0; elem<static_cast<int>(BaseType::GetModelPart().Elements().size()); elem++)
        {
            ModelPart::ElementsContainerType::iterator itElem = BaseType::GetModelPart().ElementsBegin()+elem;
            
            unsigned int numberOfNodes = itElem->GetGeometry().PointsNumber();
            
            const double permeability = itElem->GetProperties()[PERMEABILITY_WATER];

            const double density = itElem->GetValue(DENSITY_ELEM);

            //Area
            double area = 0.0;
            const int Ndim = 2;
            //GradN
            boost::numeric::ublas::bounded_matrix<double, 3, Ndim > DN_DX = ZeroMatrix(3,2);
            //N
            array_1d<double, 3 > N = ZeroVector(3);
            //Geometry staff
	    GeometryUtils::CalculateGeometryData(itElem->GetGeometry(), DN_DX, N, area);       
            
            //Compute gradient of head
            array_1d<double, Ndim > headGrad = ZeroVector(Ndim);
            
            for (unsigned int node = 0; node < numberOfNodes ; node++)
            {
               const double currentHeadLevel =  itElem->GetGeometry()[node].FastGetSolutionStepValue(PRESSURE);
               for (unsigned int dim = 0; dim < Ndim; dim++)
                {
                   headGrad[dim]+=currentHeadLevel*DN_DX.at_element(node,dim);
                }
               
            }
 
            //if we have a  buoyancy term, substract it from the gradient
            if ( isBuoyancy != 1)
            {
                for (unsigned int dim = 0; dim < Ndim; dim++)
                    headGrad[dim] -= gravity[dim] * (density) ;
            }

            //K*gradh
            boost::numeric::ublas::bounded_matrix<double,2,2> K = ZeroMatrix(2,2);
            K(0,0)= permeability;
            K(1,1)= permeability;
                 
            array_1d<double, Ndim > darcyFlow = ZeroVector(Ndim);
            
            darcyFlow = (-1.0)*prod(K,trans(headGrad));
            const double darcyFlowX = darcyFlow[0];
            const double darcyFlowY = darcyFlow[1];

            itElem->SetValue(DARCY_FLOW_X,darcyFlowX);
            itElem->SetValue(DARCY_FLOW_Y,darcyFlowY);
        } 
          
    }
    
    void calculateDensityNodes()
    {
        //rCurrentProcessInfo[FRACTIONAL_STEP] = 1;
        /*
        for(int node=0; node<static_cast<int>(BaseType::GetModelPart().Nodes().size()); node++)
        {
            ModelPart::NodesContainerType::iterator itNode = BaseType::GetModelPart().NodesBegin()+node;
            
            double densityElement= 0.0;
            double pDensity_pConc = 700.0;
            double referenceDensity = 1000.0;
            double referenceConcentration = 0.0;
        
            const double concNode = itNode->FastGetSolutionStepValue(CONCENTRATION);
            const double densNode =  (referenceDensity + pDensity_pConc * (concNode-referenceConcentration));
            
            itNode->FastGetSolutionStepValue(DENSITY) = densNode;
               
            std::cout << " conc: "<< concNode << " densNode:    " << densNode << std::endl;
            
        }
        */
        //// This part loop over Elements and after over nodes in order to compute
        //// density and densityElem. But it has more efficiency, loop over nodes in the
        //// this class, and later in the element class compute only densityElem (previuosly code)
   
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
        
    }
    
    
 
    

   double iterationConcentration()
    {
        KRATOS_TRY
        
        double normDx = this->mConcentrationStrategy->Solve();
        
        return normDx;

        KRATOS_CATCH("");
    }
    

   void StorageOldSConcIteration()
    {
        KRATOS_TRY

        for (ModelPart::NodeIterator node = BaseType::GetModelPart().NodesBegin();
                node != BaseType::GetModelPart().NodesEnd(); ++node)
        {
            const double c = (node)->FastGetSolutionStepValue(CONCENTRATION);
            (node)->FastGetSolutionStepValue(CONCENTRATION_OLD_IT) = c;
                
        }

        KRATOS_CATCH("")
    }
    
    virtual bool isConverged()
    {
        KRATOS_TRY;
        
        double maxUpdateConcentration = 0.0;

        for (ModelPart::NodeIterator node = BaseType::GetModelPart().NodesBegin();
                node != BaseType::GetModelPart().NodesEnd(); ++node)
        {
                
            const double oldIterConc = (node)->FastGetSolutionStepValue(CONCENTRATION_OLD_IT);
            const double newIterConc = (node)->FastGetSolutionStepValue(CONCENTRATION);
            const double currentDiffConcNode = std::abs (double (oldIterConc - newIterConc));
            
            if( currentDiffConcNode >= maxUpdateConcentration )
                maxUpdateConcentration = currentDiffConcNode;
                
        }
        
       
        if(maxUpdateConcentration < this->mConcentrationTol )
        {
            std::cout << "Convergence achieved, " << ", maxUpdateConcentration: "<< maxUpdateConcentration <<std::endl;
            return true;
        }
        
        std::cout << "Concentration max diffrence between iterations = " << maxUpdateConcentration << std::endl;
        
        return false;
        
         KRATOS_CATCH("")
        
    }
    
    virtual void SetEchoLevel(int Level)
    {
        mEchoLevel = Level;
        mConcentrationStrategy->SetEchoLevel(Level);
    }
    
    virtual void Clear()
    {
        KRATOS_WATCH("IterativeTransportStrategy Clear Function called");  
        this->mConcentrationStrategy->Clear();
    }
   
    virtual int Check()
    {
        KRATOS_TRY

        
        //veryfying that the model part has all the variables needed
        if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(CONCENTRATION) == false)
            KRATOS_ERROR(std::logic_error, "Add  ----CONCENTRATION---- variable!!!!!! ERROR", "");
        if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(CONCENTRATION_OLD_IT) == false)
            KRATOS_ERROR(std::logic_error, "Add  ----PRESSURE_OLD_IT---- variable!!!!!! ERROR", "");
        if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(DENSITY) == false)
            KRATOS_ERROR(std::logic_error, "Add  ----DENSITY---- variable!!!!!! ERROR", "");
        
        
        //check that the domain size is correctly prescribed
        if (this->mDomainSize != mSolverStrategyConfiguration.GetDomainSize())
            KRATOS_ERROR(std::logic_error, "domain size not coinciding", "")

            //verify buffer size
            if (BaseType::GetModelPart().GetBufferSize() < mTimeOrder + 1)
                KRATOS_ERROR(std::logic_error, "insufficient buffer size. Buffer size should be >= time_order+1", "");

        //check that, in the 2D case, the xy plane is used.
        if (this->mDomainSize == 2)
        {
            double zmin = BaseType::GetModelPart().NodesBegin()->Z();
            double zmax = zmin;
            for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
                    i != BaseType::GetModelPart().NodesEnd(); ++i)
            {
                if (i->Z() < zmin) zmin = i->Z();
                else if (i->Z() > zmax) zmax = i->Z();
            }
            if (fabs(zmax - zmin) > 1e-20)
                KRATOS_ERROR(std::logic_error, "2D model is not in the XY plane!", "")
        }

            const char ElementName[] = "Transport";
            Element const& ref_el = KratosComponents<Element>::Get(ElementName);

            for (ModelPart::ElementsContainerType::iterator it = BaseType::GetModelPart().ElementsBegin();
                    it != BaseType::GetModelPart().ElementsEnd(); ++it)
            {
                if (it->Id() < 1)
                    KRATOS_ERROR(std::logic_error, "Element Id can not be lesser than 1 (0 is not allowed as Id)", "");
                if (typeid (ref_el) != typeid (*it))
                {
                    std::cout << "wrong element found --> " << it->Id() << std::endl;
                    KRATOS_ERROR(std::logic_error, "Fractional step strategy requires Transport element for the 2D case", "");
                }
                it->Check(BaseType::GetModelPart().GetProcessInfo());
            }

        return 0;
        
        
        KRATOS_CATCH("")
    }
    
    void InitializeFractionalStep(const int aStep, const int aTimeOrder)
    {
        KRATOS_TRY;
        
        KRATOS_CATCH("");
    }
    
    void AssignInitialStepValues() 
    {
        KRATOS_TRY

        KRATOS_CATCH("");
    }
    
    void PredictSV(int step, int prediction_order)
    {
        KRATOS_TRY
 
        for (ModelPart::NodeIterator node = BaseType::GetModelPart().NodesBegin();
        node != BaseType::GetModelPart().NodesEnd(); ++node)
        {
            
            (node)->FastGetSolutionStepValue(CONCENTRATION_OLD_IT,1) = 0.0;
            (node)->FastGetSolutionStepValue(CONCENTRATION_OLD_IT) = 0.0;
            (node)->SetValue(CONCENTRATION_OLD_IT, 0.0);
        }
        
        KRATOS_CATCH("");
    }
    
        

    
protected:
    
    typename BaseType::Pointer mConcentrationStrategy; 

    double mConcentrationTol;
    int mMaxIterConcentration; 
    unsigned int mTimeOrder;
    unsigned int mPredictorOrder; 
    bool mPredictorCorrector; 
    bool mReformDofAtEachIteration;
    int mEchoLevel; 

private:    

    unsigned int mStep;
    unsigned int mDomainSize; 
       
    SolverStrategyConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>& mSolverStrategyConfiguration;

    IterativeTransportStrategy(const IterativeTransportStrategy& Other);

}; /* Class FractionalIterativeStrategy */

} /* namespace Kratos.*/

#endif /* KRATOS_ITERATIVE_TRANSPORT_STRATEGY  defined */

