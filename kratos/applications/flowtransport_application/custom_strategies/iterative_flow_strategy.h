/* *********************************************************
 *
 *   Last Modified by:    $Author: E.Bonet $
 *   Date:                $Date: 2015        $
 *   Revision:            $Revision: 1.0     $
 *
 * ***********************************************************/


#if !defined(KRATOS_ITERATIVE_FLOW_STRATEGY)
#define  KRATOS_ITERATIVE_FLOW_STRATEGY


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
class IterativeFluidStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    
    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( IterativeFluidStrategy );

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
    IterativeFluidStrategy(
    ModelPart& model_part,
    SolverStrategyConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>& aSolverStrategyConfiguration,
    bool aReformDofAtEachIteration = true,
    double aPressureTol = 0.0001,
    int aMaxIterPressure = 12,
    unsigned int aTimeOrder = 1,
    unsigned int aDomainSize = 2,
    bool aPredictorCorrector = false,
    bool aMoveMeshFlag = false, 
    unsigned int aEchoLevel = 3
    )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, aMoveMeshFlag), mSolverStrategyConfiguration(aSolverStrategyConfiguration)
    {
        KRATOS_TRY

        this->mPressureTol = aPressureTol;
        this->mMaxIterPressure = aMaxIterPressure;
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
        this->mPressureStrategy = mSolverStrategyConfiguration.pGetStrategy(std::string("PressureStrategy"));
        
        this->mStep = 1;

        KRATOS_CATCH("")
    }

    virtual ~IterativeFluidStrategy()
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
            //perform one iteration over Pressure
            //if (rank == 0)
            ////////////////////////////////////////////////////////////////////
            std::cout << "New iteration. Iteration = " << iteration  << std::endl;
            std::ofstream myfile;
            myfile.open ("MatrixSystem.txt", std::ios::app);
            //myfile << "New iteration. Iteration = " << iteration  << std::endl;
            //////////////////////////////////////////////////////////////////////
            
            if(this->mEchoLevel == 3)
                std::cout << "FLOW EQUATION SOLVE: "  << std::endl;
            
            
            if(iteration != 0)
                this->StorageOldSPressIteration();
            
            double pressureNormDx = iterationPressure();
            
            std::cout << "Pressure norm = " << pressureNormDx  << std::endl;
            
            //this->mPressureStrategy->Clear();
            
            //Calculate at element flowpressuretrans

            this->CalculateFlowBalances(); 
            
            
            iteration++;
            
            /////////////////////////////////////////////////
            std::cout << "End of iteration = " << iteration  << std::endl;
            myfile.open ("MatrixSystemFlow.txt", std::ios::app);
            myfile << "End of iteration = " << iteration  << std::endl;  
            myfile.close();
            /////////////////////////////////////////////////
            
            if( isConverged() )
            {
                isGlobalConverged = true;
                std::cout << "CONVERGENCE  ACHIEVED for time step number: " << this->mStep << std::endl;
                std::cout << "It reachs at " << iteration << " iteration" << std::endl;
                
                /////////////////////////////////////////////////
                myfile.open ("MatrixSystemFlow.txt", std::ios::app);
                myfile << "CONVERGENCE  ACHIEVED for time step number: " << this->mStep << std::endl;
                myfile << "It reachs at " << iteration << " iteration" << std::endl << std::endl;
                myfile.precision(10);
                for(int elem=0; elem<static_cast<int>(BaseType::GetModelPart().Elements().size()); elem++)
                {
                    ModelPart::ElementsContainerType::iterator itElem = BaseType::GetModelPart().ElementsBegin()+elem;
                    
                    const double densityElem = itElem->GetValue(DENSITY_ELEM);
                    myfile << "Elem: " <<  itElem->Id() << std::endl;
                    myfile << "DENSITY_ELEM   " << densityElem << std::endl;
                    
                    for (unsigned int node = 0; node < itElem->GetGeometry().PointsNumber(); node++)
                    { 
                    
                    const double pressNode = itElem->GetGeometry()[node].FastGetSolutionStepValue(PRESSURE);
                    
                    myfile << "Node: " <<  itElem->GetGeometry()[node].Id() << "   PRESSURE   " << pressNode  << std::endl;
                    
                    }
                    
                    myfile << std::endl;
                }
                myfile << std::endl;
                for (ModelPart::NodeIterator node = BaseType::GetModelPart().NodesBegin(); node != BaseType::GetModelPart().NodesEnd(); ++node)
                {
                    const double storageBal =  node->FastGetSolutionStepValue(STORAGE_BALANCE);
                    const double darcyFlowBal =  node->FastGetSolutionStepValue(DARCY_FLOW_BALANCE);
                    const double bal =  node->FastGetSolutionStepValue(SINKSOURCE_BALANCE);
                    
                    myfile << "BALANCE   " << std::endl;
                    myfile << "Node: " <<  node->Id() << "   ST   " << storageBal << "   DF   " << darcyFlowBal << "   SS   " << bal << std::endl;
                    
                }
                myfile << std::endl << std::endl;
           
                myfile.close();
                //////////////////////////////////////////////////
            }
            
            if(!isGlobalConverged && (iteration > this->mMaxIterPressure))
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
    
    
    
    
    /////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    // DarcyFlow is calculate at the flowpressuretrans elem (more efficient)
    

    
    
    void CalculateFlowBalances()
    {
        
        for (ModelPart::NodeIterator node = BaseType::GetModelPart().NodesBegin();
                node != BaseType::GetModelPart().NodesEnd(); ++node)
        {
                node->FastGetSolutionStepValue(STORAGE_BALANCE) = 0.0;
                node->FastGetSolutionStepValue(DARCY_FLOW_BALANCE) = 0.0;
                node->FastGetSolutionStepValue(SINKSOURCE_BALANCE) = 0.0;
        }   
        
        ProcessInfo& CurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
           
        unsigned int isFlowStationary = CurrentProcessInfo[IS_FLOW_STATIONARY];
        const double deltaTime = CurrentProcessInfo.GetValue(DELTA_TIME);
        
        array_1d<double,2> gravity= CurrentProcessInfo[GRAVITY];
        
        for(int elem=0; elem<static_cast<int>(BaseType::GetModelPart().Elements().size()); elem++)
        {
            
            ModelPart::ElementsContainerType::iterator itElem = BaseType::GetModelPart().ElementsBegin()+elem;

            unsigned int numberOfNodes = itElem->GetGeometry().PointsNumber();
                
            //area
            double area = 0.0;
            const int Ndim = 2;
            
            //GradGradN
            boost::numeric::ublas::bounded_matrix<double, 3, 3 > DN_DXX = ZeroMatrix(3,3);
            //GradN
            boost::numeric::ublas::bounded_matrix<double, 3, Ndim > DN_DX = ZeroMatrix(3,Ndim);
             //N
            array_1d<double, 3 > N = ZeroVector(3);
            //Geometry staff
            GeometryUtils::CalculateGeometryData(itElem->GetGeometry(),DN_DXX, DN_DX, N, area); 
            
            const double permeability = itElem->GetProperties()[PERMEABILITY_WATER];

            const double densityElem = itElem->GetValue(DENSITY_ELEM);
   
            //Matrix MatElemHeadLevel
            boost::numeric::ublas::bounded_matrix<double, 3, 3 > MatElemHeadLevel = ZeroMatrix(3,3);
            
            const int ntensor=3;
            //K
            array_1d<double, 3 > K = ZeroVector(3);
            K[0]= permeability;
            K[1]= 0.0;
            K[2]= permeability;

	    int k=0;    // Row of GradGrad matrix

	   for (int i = 0; i< numberOfNodes-1; i++)    // Loop of columns of local matrix
	   {
               	for (int j = i+1; j < numberOfNodes; j++)	// Loop of columns of local matrix
                {
		    for (int itensor = 0; itensor < ntensor; itensor++)	// Loop columns of GradGrad matrix (or of tensor values)
		    {
			double elem_ij = -DN_DXX.at_element(k, itensor)* K[itensor];
			MatElemHeadLevel(i,j) += elem_ij;
			MatElemHeadLevel(j,i) += elem_ij;
                        MatElemHeadLevel(i,i) -= elem_ij;
                        MatElemHeadLevel(j,j) -= elem_ij;

		    }
		    k++;
		}
            }

           array_1d<double, 3 > Pressurek_1 = ZeroVector(3);
           for (unsigned int node = 0; node < numberOfNodes; node++)
               Pressurek_1[node] =  itElem->GetGeometry()[node].FastGetSolutionStepValue(PRESSURE);

            array_1d<double, 3 > balanceHeadLevel = ZeroVector(3);
            noalias(balanceHeadLevel) = prod (MatElemHeadLevel,Pressurek_1);
   
            //Vector MatElemHeadLevel
             array_1d <double, 3 > MatElemK_Rho_g = ZeroVector(3);
            
            //loop over nodes of element
            for (int inode = 0; inode < numberOfNodes; inode++)
            {
                double result = 0;
                //double densityNode = itElem->GetGeometry()[inode].FastGetSolutionStepValue(DENSITY);
                //loop over spatial dimensions
                for (int idim = 0; idim < Ndim; idim++)
                {   
                    const int k_index =  idim*Ndim; //Only 2D           
                    result += K[k_index]* gravity[idim] * densityElem * DN_DX.at_element(inode,idim)*area;
                }

                //store result 
               MatElemK_Rho_g[inode] = result;
            
               
            }
            
             balanceHeadLevel += MatElemK_Rho_g;
             
                
            for (unsigned int node = 0; node < numberOfNodes; node++)
            {
                if( isFlowStationary != 1 )
                {
                    //Storage balance
                    const double specificStorage = itElem->GetProperties()[SPECIFIC_STORAGE];
                    double Pk = itElem->GetGeometry()[node].FastGetSolutionStepValue(PRESSURE, 1);
                    if(this->mStep == 1)
                          Pk = itElem->GetGeometry()[node].GetValue(PRESSURE);
                         
                    const double Pk_1 = itElem->GetGeometry()[node].FastGetSolutionStepValue(PRESSURE);
                    const double storageBalance = ( /*densityElem */ N[node] * area * specificStorage * (Pk_1 - Pk) ) / (deltaTime);

                    itElem->GetGeometry()[node].FastGetSolutionStepValue(STORAGE_BALANCE) += storageBalance;
                     KRATOS_WATCH(this->mStep);
                }
                
                const double darcyFlowBalance = balanceHeadLevel[node];
                itElem->GetGeometry()[node].FastGetSolutionStepValue(DARCY_FLOW_BALANCE) += darcyFlowBalance;

            }  
                    
        }
        
       for (ModelPart::NodeIterator node = BaseType::GetModelPart().NodesBegin();
                node != BaseType::GetModelPart().NodesEnd(); ++node)
        {
           const double storageBal =  node->FastGetSolutionStepValue(STORAGE_BALANCE);
           const double darcyFlowBal =  node->FastGetSolutionStepValue(DARCY_FLOW_BALANCE);
           const double bal =  storageBal - darcyFlowBal;
           node->FastGetSolutionStepValue(SINKSOURCE_BALANCE) = bal;
                                                                             
        } 
        
    }
    
    double iterationPressure()
    {
        KRATOS_TRY

        double normDx = this->mPressureStrategy->Solve();
        
        return normDx;

        KRATOS_CATCH("");
    }
    
    
    void StorageOldSPressIteration()
    {
        KRATOS_TRY

        for (ModelPart::NodeIterator node = BaseType::GetModelPart().NodesBegin();
                node != BaseType::GetModelPart().NodesEnd(); ++node)
        {
            //setting the old value of the pressure & concentration to the current one
            const double p1 = (node)->FastGetSolutionStepValue(PRESSURE);
            (node)->FastGetSolutionStepValue(PRESSURE_OLD_IT) = p1;
            
            const double p2 = (node)->FastGetSolutionStepValue(PRESSURE, 1);
            (node)->FastGetSolutionStepValue(PRESSURE_OLD_IT,1) = p2;
            
            const double p3 = (node)->GetValue(PRESSURE);
            (node)->SetValue(PRESSURE_OLD_IT,p3);
        }

        KRATOS_CATCH("")
    }
    
    
    virtual bool isConverged()
    {
        KRATOS_TRY;
        
        double maxUpdatePressure = 0.0;

        for (ModelPart::NodeIterator node = BaseType::GetModelPart().NodesBegin();
                node != BaseType::GetModelPart().NodesEnd(); ++node)
        {
            const double oldIterPressure = (node)->FastGetSolutionStepValue(PRESSURE_OLD_IT);
            const double newIterPressure = (node)->FastGetSolutionStepValue(PRESSURE);
            const double currentDiffPressNode = std::abs(double (oldIterPressure - newIterPressure));
            
            if( currentDiffPressNode >= maxUpdatePressure )
                maxUpdatePressure = currentDiffPressNode;
        
        }
        
       
        if(maxUpdatePressure < this->mPressureTol)
        {
            std::cout << "Convergence achieved, " <<"maxUpdatePressure: "<< maxUpdatePressure <<std::endl;
            return true;
        }
        
        std::cout << "Pressure max diffrence between iterations = " << maxUpdatePressure << std::endl;
        
        return false;
        
         KRATOS_CATCH("")
        
    }
    
    virtual void SetEchoLevel(int Level)
    {
        mEchoLevel = Level;
        mPressureStrategy->SetEchoLevel(Level);
    }
    
    virtual void Clear()
    {
        KRATOS_WATCH("FractionalIterativeStrategy Clear Function called");  
        this->mPressureStrategy->Clear();
    }
   
    virtual int Check()
    {
        KRATOS_TRY

        
        //veryfying that the model part has all the variables needed
        if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(PRESSURE) == false)
            KRATOS_ERROR(std::logic_error, "Add  ----PRESSURE---- variable!!!!!! ERROR", "");
        if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(PRESSURE_OLD_IT) == false)
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

            const char ElementName[] = "FlowDarcy";
            Element const& ref_el = KratosComponents<Element>::Get(ElementName);

            for (ModelPart::ElementsContainerType::iterator it = BaseType::GetModelPart().ElementsBegin();
                    it != BaseType::GetModelPart().ElementsEnd(); ++it)
            {
                if (it->Id() < 1)
                    KRATOS_ERROR(std::logic_error, "Element Id can not be lesser than 1 (0 is not allowed as Id)", "");
                if (typeid (ref_el) != typeid (*it))
                {
                    std::cout << "wrong element found --> " << it->Id() << std::endl;
                    KRATOS_ERROR(std::logic_error, "Fractional step strategy requires FlowDarcy2D element for the 2D case", "");
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
            (node)->FastGetSolutionStepValue(PRESSURE_OLD_IT,1) = 0.0;
            (node)->FastGetSolutionStepValue(PRESSURE_OLD_IT) = 0.0;
            (node)->SetValue(PRESSURE_OLD_IT, 0.0);
            
        }
        
        KRATOS_CATCH("");
    }
    
  
    
protected:
    
    typename BaseType::Pointer mPressureStrategy; 

    double mPressureTol; 
    int mMaxIterPressure; 
    unsigned int mTimeOrder;
    unsigned int mPredictorOrder; 
    bool mPredictorCorrector; 
    bool mReformDofAtEachIteration;
    int mEchoLevel; 

private:    

    unsigned int mStep;
    unsigned int mDomainSize; 
       
    SolverStrategyConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>& mSolverStrategyConfiguration;

    IterativeFluidStrategy(const IterativeFluidStrategy& Other);

}; /* Class FractionalIterativeStrategy */

} /* namespace Kratos.*/

#endif /* KRATOS_ITERATIVE_FLOW_STRATEGY  defined */

