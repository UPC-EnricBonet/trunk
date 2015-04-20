/* *********************************************************
 *
 *   Last Modified by:    $Author: rrossi $
 *   Date:                $Date: 2009-01-13 15:39:56 $
 *   Revision:            $Revision: 1.14 $
 *
 * ***********************************************************/


#if !defined(KRATOS_FRACTIONAL_ITERATIVE_STRATEGY)
#define  KRATOS_FRACTIONAL_ITERATIVE_STRATEGY


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
//#include "custom_strategies/solver_configuration.h"
#include "custom_strategies/fractional_iterative_configuration.h"
//#include "custom_strategies/fractionalstep_configuration.h"


namespace Kratos
{
    
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class FractionalIterativeStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    
    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( FractionalIterativeStrategy );

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
    FractionalIterativeStrategy(
    ModelPart& model_part,
    SolverStrategyConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>& aSolverStrategyConfiguration,
    bool aReformDofAtEachIteration = true,
    double aPressureTol = 0.00001,
    double aConcentrationTol = 0.0000001,
    int aMaxIterPressure = 12,
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

        this->mPressureTol = aPressureTol;
        this->mConcentrationTol = aConcentrationTol;
        this->mMaxIterPressure = aMaxIterPressure;
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
        
        //initialize strategy (Matriz management??)
        CurrentProcessInfo[FRACTIONAL_STEP] = 1;
        this->mPressureStrategy = mSolverStrategyConfiguration.pGetStrategy(std::string("PressureStrategy"));
        CurrentProcessInfo[FRACTIONAL_STEP] = 2;
        this->mConcentrationStrategy = mSolverStrategyConfiguration.pGetStrategy(std::string("ConcentrationStrategy"));
        
        //Doc -> fijar number step
        this->mStep = 1;

        KRATOS_CATCH("")
    }

    virtual ~FractionalIterativeStrategy()
    {
        
    }

    double Solve()
    {
        KRATOS_TRY
        
        Timer::Start("solve");
        
        //Initial Guess, if not is by default method
        //Predict();
        
        //
        
        //mCalculateNormDxFlag = false;
        double Dp_norm = 0.0;
        //Doc -> assign the correct fractional step coefficients (BDF_COEFFICIENTS..)
        //by now, empty
        InitializeFractionalStep(this->mStep, this->mTimeOrder);
        
        //Doc -> initial guess for velocity from interpolation from timestep before velocity
        //by now, empty
        PredictSV(this->mStep, this->mPredictorOrder);

        //Doc -> Save old iteration to compare in convergence with the new pressure iteration. Save always that repit iterations.
        //Assign Velocity To Fract Step Velocity and Node Area to Zero
        AssignInitialStepValues();

        // Loop convergence
        
        double iteration = 0; 
        bool isGlobalConverged = false;
        //solve for fractional step velocities
        while (!isGlobalConverged)
        {
            //perform one iteration over Pressure and concentration
            //if (rank == 0)
            std::cout << "New iteration. Iteration = " << iteration  << std::endl;
            
            if(this->mEchoLevel == 3)
                std::cout << "FLOW EQUATION SOLVE: "  << std::endl;
            
            this->StorageOldSPressIteration();
            
            this->calculateDensityNodes();
            
            double pressureNormDx = iterationPressure();
            
            this->CalculateDarcyFlowPressure();
            
            std::cout << "Pressure norm = " << pressureNormDx  << std::endl;
            
            if(this->mEchoLevel == 3)
                std::cout << "TRANSPORT EQUATION SOLVE: "   << std::endl;
            
            this->StorageOldSConcIteration(); 
            
            //this->CalculateFlowBalances();  
            
            double concentrationNormDx = iterationConcentration();

            
            std::cout << "Concentration norm = " << concentrationNormDx  << std::endl;
            
            iteration++;
            
            std::cout << "End of iteration = " << iteration  << std::endl;
            
            if( isConverged() )
            {
                isGlobalConverged = true;
                std::cout << "ATTENTION: convergence  ACHIEVED, it reachs at " << iteration << " iteration" << std::endl;
            }
            
            if(!isGlobalConverged && (iteration > this->mMaxIterPressure || iteration > this->mMaxIterConcentration))
            {
                std::cout << "ATTENTION: convergence NOT achieved iteration  " << iteration << " reached" << std::endl;
                break;
            }
            
        }

        if (this->mReformDofAtEachIteration == true)
            this->Clear();

        //Doc -> timeStep+1
        this->mStep += 1;
        Timer::Stop("solve");
        
        return 0.0;
        
        KRATOS_CATCH("")
    }
    
    void calculateDensityNodes()
    {
        //rCurrentProcessInfo[FRACTIONAL_STEP] = 1;
        
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
        
        //// This part loop over Elements and after over nodes in order to compute
        //// density and densityElem. But it has more efficiency, loop over nodes in the
        //// this class, and later in the element class compute only densityElem (previuosly code)
        /*
        for(int elem=0; elem<static_cast<int>(BaseType::GetModelPart().Elements().size()); elem++)
        {
            ModelPart::ElementsContainerType::iterator itElem = BaseType::GetModelPart().ElementsBegin()+elem;

            unsigned int numberOfNodes = itElem->GetGeometry().PointsNumber();
            
            double densityElement= 0.0;
            double pDensity_pConc = 700.0;
            double referenceDensity = 1000.0;
            double referenceConcentration = 0.0; 
            double weight = 1.0 / ((double) numberOfNodes);
        
            for (unsigned int node = 0; node < numberOfNodes; node++)
            {    
               const double concNode = itElem->GetGeometry()[node].FastGetSolutionStepValue(CONCENTRATION);
               const double densNode =  (referenceDensity + pDensity_pConc * (concNode-referenceConcentration));
               itElem->GetGeometry()[node].FastGetSolutionStepValue(DENSITY) = densNode;
               densityElement += (weight *densNode);
               
               std::cout << "node:  " << node << " conc: "<< concNode << " densNode:    " << densNode << std::endl;
               std::cout << std::endl;
            }
            
            itElem->GetValue(b/*DENSITY_ELEM) = densityElement;
            
            std::cout << "   densElement  " << densityElement << std::endl;
        }
        */
    }
  
    void CalculateDarcyFlowPressure()
    {
        
        ProcessInfo& CurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
        
        unsigned int isBuoyancy = 0;//CurrentProcessInfo[IS_BUOYANCY];
        array_1d<double,2> gravity= CurrentProcessInfo[GRAVITY];
        
        //#pragma omp parallel for firstprivate(rhs,lhs)
        for(int elem=0; elem<static_cast<int>(BaseType::GetModelPart().Elements().size()); elem++)
        {
            ModelPart::ElementsContainerType::iterator itElem = BaseType::GetModelPart().ElementsBegin()+elem;
            
            unsigned int numberOfNodes = itElem->GetGeometry().PointsNumber();
            
            const double permeability = itElem->GetProperties()[PERMEABILITY_WATER];

            const double density = itElem->GetValue(b/*DENSITY_ELEM*/);
            
            //Area
            double area;
            const int Ndim = 2;
            //GradN
            boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX;
            //N
            array_1d<double, 3 > N;
            //Geometry staff
	    GeometryUtils::CalculateGeometryData(itElem->GetGeometry(), DN_DX, N, area);       
            
            //Compute gradient of head
            array_1d<double, Ndim > headGrad;
            

            for (unsigned int node = 0; node < numberOfNodes ; node++)
            {
               double currentHeadLevel =  itElem->GetGeometry()[node].FastGetSolutionStepValue(PRESSURE);
               for (unsigned int dim = 0; dim < Ndim; dim++)
                {
                   headGrad[dim]+=currentHeadLevel*DN_DX.at_element(node,dim);
                }
               
            }
 
            //if we have a  buoyancy term, substract it from the gradient
            if ( isBuoyancy != 1)
            {
                array_1d<double, Ndim > grav;
                for (unsigned int dim = 0; dim < Ndim; dim++)
                {
                    if(dim == Ndim-1)
                         grav[dim]=9.8;
                     else
                         grav[dim]=0.0;

                    headGrad[dim] += grav[dim] * (density) ;
                }
            }
            

            array_1d<double, 2 > darcyFlow;
	    double darcyFlowX, darcyFlowY;
            
            //K*gradh
            boost::numeric::ublas::bounded_matrix<double,2,2> K = ZeroMatrix(2,2);
            K(0,0)= permeability;
            K(1,1)= permeability;
                 
            darcyFlow=(-1.0)*prod(K,trans(headGrad));
            darcyFlowX=darcyFlow[0];
            darcyFlowY=darcyFlow[1];
            
            itElem->GetValue(IMPOSED_VELOCITY_X_VALUE /*DARCY_FLOW_X*/) = darcyFlowX;
            itElem->GetValue(IMPOSED_VELOCITY_Y_VALUE /*DARCY_FLOW_Y*/) = darcyFlowY;
        } 
           
    }
    
    void CalculateFlowBalances()
    {
        /*
        ProcessInfo& CurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
            
        unsigned int isBuoyancy = 0;//CurrentProcessInfo[IS_BUOYANCY];
        unsigned int isFlowStationary = 1; //rCurrentProcessInfo[IS_FLOW_STATIONARY];
        const double deltaTime = CurrentProcessInfo.GetValue(DELTA_TIME);
        array_1d<double,2> gravity= CurrentProcessInfo[GRAVITY];
        
        //BULK_MODULUS

        for(int elem=0; elem<static_cast<int>(BaseType::GetModelPart().Elements().size()); elem++)
        {
            ModelPart::ElementsContainerType::iterator itElem = BaseType::GetModelPart().ElementsBegin()+elem;
            
            unsigned int numberOfNodes = itElem->GetGeometry().PointsNumber();
            
            const double permeability = itElem->GetProperties()[PERMEABILITY_WATER];
            
            double specificStorage;
            if(isFlowStationary!=1)
                specificStorage = itElem->GetProperties()[BULK_MODULUS/*SPECIFIC_STORAGE];
                
            const double density = itElem->GetValue(DENSITY_ELEM);
           
            //area
            double area;
            const int nDim = 2;
            
            //GradN
            boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX;
            //N
            array_1d<double, 3 > N;
            //Geometry staff
	    GeometryUtils::CalculateGeometryData(itElem->GetGeometry(), DN_DX, N, area);          
            
            //Compute gradient of head
            boost::numeric::ublas::bounded_matrix<double, nDim, 3> headGrad_i;

            for (unsigned int node = 0; node < numberOfNodes; node++)
            {
               double currentHeadLevel =  itElem->GetGeometry()[node].FastGetSolutionStepValue(PRESSURE);
               for (unsigned int dim = 0; dim < nDim; dim++)
                {
                   headGrad_i.at_element(node,dim)+=currentHeadLevel* DN_DX.at_element(node,dim);
                }  
            }

            
            //if we have a  buoyancy term, substract it from the gradient
            if ( isBuoyancy != 1)
            {
                array_1d<double, nDim > grav;
                for (unsigned int dim = 0; dim < nDim; dim++)
                {
                    if(dim == nDim-1)
                         grav[dim]=-9.8;
                     else
                         grav[dim]=0.0;
                    
                    for (unsigned int node = 0; node < numberOfNodes; node++)
                    {
                        const double density_i = itElem->GetGeometry()[node].FastGetSolutionStepValue(DENSITY);
                        headGrad_i.at_element(node,dim)+= grav[dim] * density_i;  
                    }
                }
            }
            
            //K*gradh
            boost::numeric::ublas::bounded_matrix<double,2,2> K = ZeroMatrix(2,2);
            K(0,0)= density * permeability;
            K(1,1)= density * permeability;
            
            boost::numeric::ublas::bounded_matrix<double, nDim, 3> darcyFlow_i;
	    double darcyFlowX_i, darcyFlowY_i;
            
            for (unsigned int node = 0; node < numberOfNodes; node++)
            {
                itElem->GetGeometry()[node].FastGetSolutionStepValue(STORAGE_BALANCE) = 0.0;
                itElem->GetGeometry()[node].FastGetSolutionStepValue(DARCY_FLOW_BALANCE) = 0.0;
                itElem->GetGeometry()[node].FastGetSolutionStepValue(SINKSOURCE_BALANCE) = 0.0;
                
                //Storage
                if( isFlowStationary != 1 )
                {
                    double Pk = itElem->GetGeometry()[node].FastGetSolutionStepValue(PRESSURE, 1);
                    double Pk_1 = itElem->GetGeometry()[node].FastGetSolutionStepValue(PRESSURE);
                    double storageBalance = ( specificStorage * (Pk_1 - Pk) ) / (deltaTime);
                    itElem->GetGeometry()[node].FastGetSolutionStepValue(STORAGE_BALANCE) += storageBalance;
                    itElem->GetGeometry()[node].FastGetSolutionStepValue(SINKSOURCE_BALANCE) += storageBalance;
                }
                
                //DarcyFlow
                darcyFlow_i[node]=-1*prod(K,trans(headGrad_i[node]));
                darcyFlowX_i=darcyFlow_i[node][0];
                darcyFlowY_i=darcyFlow_i[node][1];
                double darcyFlowBalance = darcyFlowX_i*DN_DX.at_element(node,0) +
                                        darcyFlowY_i*DN_DX.at_element(node,1);
                itElem->GetGeometry()[node].FastGetSolutionStepValue(DARCY_FLOW_BALANCE) += darcyFlowBalance;
                
                //Total balance
                itElem->GetGeometry()[node].FastGetSolutionStepValue(SINKSOURCE_BALANCE) -= darcyFlowBalance;
            }
            
        } 
       */
    }  
    
    double iterationPressure()
    {
        KRATOS_TRY

        ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();

        rCurrentProcessInfo[FRACTIONAL_STEP] = 1;

        //int solstep = pCurrentProcessInfo.GetCurrentSolutionStep();
        //DofsArrayType& rDofSet = mSolverStrategyConfiguration.pGetBuildAndSolve(rCurrentProcessInfo)->GetDofSet();
       
        //mSolverStrategyConfiguration.pSetUpDof(rCurrentProcessInfo);
        
        double normDx = this->mPressureStrategy->Solve();
        
        return normDx;

        KRATOS_CATCH("");
    }
    
   double iterationConcentration()
    {
        KRATOS_TRY

        ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();

        rCurrentProcessInfo[FRACTIONAL_STEP] = 2;
        
        //mSolverStrategyConfiguration.pSetUpDof(rCurrentProcessInfo);
        
        double normDx = this->mConcentrationStrategy->Solve();
        
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
            const double p = (node)->FastGetSolutionStepValue(PRESSURE);
            (node)->FastGetSolutionStepValue(PRESSURE_OLD_IT) = p;
     
        }

        KRATOS_CATCH("")
    }
    
        void StorageOldSConcIteration()
    {
        KRATOS_TRY

        for (ModelPart::NodeIterator node = BaseType::GetModelPart().NodesBegin();
                node != BaseType::GetModelPart().NodesEnd(); ++node)
        {
            const double c = (node)->FastGetSolutionStepValue(CONCENTRATION);
            (node)->FastGetSolutionStepValue(a/*CONCENTRATION_OLD_IT*/ ) = c;
                
        }

        KRATOS_CATCH("")
    }
    
    virtual bool isConverged()
    {
        KRATOS_TRY;
        
        double maxUpdatePressure = 0.0;
        double maxUpdateConcentration = 0.0;

        for (ModelPart::NodeIterator node = BaseType::GetModelPart().NodesBegin();
                node != BaseType::GetModelPart().NodesEnd(); ++node)
        {
            const double oldIterPressure = (node)->FastGetSolutionStepValue(PRESSURE_OLD_IT);
            const double newIterPressure = (node)->FastGetSolutionStepValue(PRESSURE);
            const double currentDiffPressNode = std::abs(double (oldIterPressure - newIterPressure));
            
            if( currentDiffPressNode >= maxUpdatePressure )
                maxUpdatePressure = currentDiffPressNode;
                
            const double oldIterConc = (node)->FastGetSolutionStepValue(a /*CONCENTRATION_OLD_IT*/);
            const double newIterConc = (node)->FastGetSolutionStepValue(CONCENTRATION);
            const double currentDiffConcNode = std::abs (double (oldIterConc - newIterConc));
            
            if( currentDiffConcNode >= maxUpdateConcentration )
                maxUpdateConcentration = currentDiffConcNode;
                
        }
        
       
        if(maxUpdatePressure < this->mPressureTol && maxUpdateConcentration < this->mConcentrationTol )
        {
            std::cout << "convergence achieved, " <<"maxUpdatePressure: "<< maxUpdatePressure << ", maxUpdateConcentration: "<< maxUpdateConcentration <<std::endl;
            return true;
        }
        
        std::cout << "Pressure max diffrence between iterations = " << maxUpdatePressure << std::endl;
        std::cout << "Concentration max diffrence between iterations = " << maxUpdateConcentration << std::endl;
        
        return false;
        
         KRATOS_CATCH("")
        
    }
    
    virtual void SetEchoLevel(int Level)
    {
        mEchoLevel = Level;
        mPressureStrategy->SetEchoLevel(Level);
        mConcentrationStrategy->SetEchoLevel(Level);
    }
    
    //Dubbio  
    virtual double GetStageResidualNorm(unsigned int step)
    {
        /*
        if (step <= 3)
            return this->mPressureStrategy->GetResidualNorm();
        if (step == 4)
            return this->mConcentrationStrategy->GetResidualNorm();
        else
         */
        
        return 0.0;
    }
    
    virtual void Clear()
    {
        KRATOS_WATCH("FractionalIterativeStrategy Clear Function called");  
        this->mPressureStrategy->Clear();
        this->mConcentrationStrategy->Clear();
    }
   
    virtual int Check()
    {
        KRATOS_TRY

        
        //veryfying that the model part has all the variables needed
        if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(PRESSURE) == false)
            KRATOS_ERROR(std::logic_error, "Add  ----PRESSURE---- variable!!!!!! ERROR", "");
        if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(CONCENTRATION) == false)
            KRATOS_ERROR(std::logic_error, "Add  ----CONCENTRATION---- variable!!!!!! ERROR", "");
        if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(PRESSURE_OLD_IT) == false)
            KRATOS_ERROR(std::logic_error, "Add  ----PRESSURE_OLD_IT---- variable!!!!!! ERROR", "");
        if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(a) == false)
            KRATOS_ERROR(std::logic_error, "Add  ----a---- variable!!!!!! ERROR", "");
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

            const char ElementName[] = "FlowPressureTrans2D";
            Element const& ref_el = KratosComponents<Element>::Get(ElementName);

            for (ModelPart::ElementsContainerType::iterator it = BaseType::GetModelPart().ElementsBegin();
                    it != BaseType::GetModelPart().ElementsEnd(); ++it)
            {
                if (it->Id() < 1)
                    KRATOS_ERROR(std::logic_error, "Element Id can not be lesser than 1 (0 is not allowed as Id)", "");
                if (typeid (ref_el) != typeid (*it))
                {
                    std::cout << "wrong element found --> " << it->Id() << std::endl;
                    KRATOS_ERROR(std::logic_error, "Fractional step strategy requires FlowPressureTrans2D element for the 2D case", "");
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
 
        KRATOS_CATCH("");
    }
    
        
    /*
    void calculateDensityDomain()
    {
        //#pragma omp parallel for firstprivate(rhs,lhs)
        for(int elem=0; elem<static_cast<int>(BaseType::GetModelPart().Elements().size()); elem++)
        {
            ModelPart::ElementsContainerType::iterator itElem = BaseType::GetModelPart().ElementsBegin()+elem;
            itElem->calculateDensity();
        }
    }
    
    void CalculateFlowBalancesDomain()
    {
        ProcessInfo& CurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
        
        //#pragma omp parallel for firstprivate(rhs,lhs)
        for(int elem=0; elem<static_cast<int>(BaseType::GetModelPart().Elements().size()); elem++)
        {
            ModelPart::ElementsContainerType::iterator itElem = BaseType::GetModelPart().ElementsBegin()+elem;
            itElem->CalculateFlowBalances(CurrentProcessInfo);
        }
    }
        
    void CalculateDarcyFlowPressureDomain()
    {
        ProcessInfo& CurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
        
        //#pragma omp parallel for firstprivate(rhs,lhs)
        for(int elem=0; elem<static_cast<int>(BaseType::GetModelPart().Elements().size()); elem++)
        {
            ModelPart::ElementsContainerType::iterator itElem = BaseType::GetModelPart().ElementsBegin()+elem;
            itElem->CalculateDarcyFlowPressure(CurrentProcessInfo);
        }
    }
 */   
    
    void calculateAtPointAndBalance()
    {
        /*
        ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
        
        rCurrentProcessInfo[FRACTIONAL_STEP] = 2;

        Vector rhs(this->mDomainSize+1);
        Matrix lhs(this->mDomainSize+1,this->mDomainSize+1);

            //#pragma omp parallel for firstprivate(rhs,lhs)
        for(int i=0; i<static_cast<int>(BaseType::GetModelPart().Elements().size()); i++)
        {
                ModelPart::ElementsContainerType::iterator it = BaseType::GetModelPart().ElementsBegin()+i;

                it->CalculateLocalSystem(lhs,rhs,rCurrentProcessInfo);
        }
        */
    }
    
protected:
    
    typename BaseType::Pointer mPressureStrategy; 
    typename BaseType::Pointer mConcentrationStrategy; 

    double mPressureTol; 
    double mConcentrationTol;
    int mMaxIterPressure; 
    int mMaxIterConcentration; 
    unsigned int mTimeOrder;
    unsigned int mPredictorOrder; 
    bool mPredictorCorrector; 
    bool mReformDofAtEachIteration;
    int mEchoLevel; 

private:    
    //Dubbio
    unsigned int mStep;
    unsigned int mDomainSize; 
    
    //Doc -> mSolverConfigAttribute    
    SolverStrategyConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>& mSolverStrategyConfiguration;

    FractionalIterativeStrategy(const FractionalIterativeStrategy& Other);

}; /* Class FractionalIterativeStrategy */

} /* namespace Kratos.*/

#endif /* KRATOS_FRACTIONAL_ITERATIVE_STRATEGY  defined */

