//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_FLOWTRANSPORT_APPLICATION_H_INCLUDED )
#define  KRATOS_FLOWTRANSPORT_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"


#include "includes/condition.h"  
#include "custom_conditions/pointflow.h" 
#include "custom_conditions/pointsource.h" 
#include "custom_conditions/pointheat.h"
#include "custom_conditions/sinksource.h"
#include "custom_conditions/sinksourcepressure.h"
#include "custom_conditions/leakage.h"
#include "custom_conditions/massflow.h"       
#include "includes/ublas_interface.h"

#include "custom_elements/poisson_2d.h" 
#include "custom_elements/diffconv_2d.h"
#include "custom_elements/diffusion_2d.h"
#include "custom_elements/diffconvst_2d.h"
#include "custom_elements/flow.h"
#include "custom_elements/darcy.h"
#include "custom_elements/flowDarcy.h"
#include "custom_elements/transport.h"


namespace Kratos
{

	///@name Kratos Globals
	///@{ 

	// Initial variables definition with the flow-transport problem
        KRATOS_DEFINE_VARIABLE(double, POINT_HEAT_SOURCE)
        KRATOS_DEFINE_VARIABLE(double, POINT_FLOW_SOURCE)

        //Variables of the iteration flow-transport problem
        KRATOS_DEFINE_VARIABLE(int, IS_FLOW_STATIONARY)
        KRATOS_DEFINE_VARIABLE(int, IS_TRANSPORT_STATIONARY)
        KRATOS_DEFINE_VARIABLE(int, IS_BUOYANCY)
        
        KRATOS_DEFINE_VARIABLE(double, HEAD_LEVEL)
	KRATOS_DEFINE_VARIABLE(double, SPECIFIC_STORAGE)
 
        KRATOS_DEFINE_VARIABLE(double, SINK_SOURCE)
        
        KRATOS_DEFINE_VARIABLE(double, LEAKAGE_COEFFICIENT)
        KRATOS_DEFINE_VARIABLE(double, LEVEL)
 
        KRATOS_DEFINE_VARIABLE(double, DENSITY_ELEM)
        KRATOS_DEFINE_VARIABLE(double, DARCY_FLOW_X)
        KRATOS_DEFINE_VARIABLE(double, DARCY_FLOW_Y)
        KRATOS_DEFINE_VARIABLE(double, DIFFUSION_COEFFICIENT)
        
        KRATOS_DEFINE_VARIABLE(double, PRESCRIBED_VALUE)

        KRATOS_DEFINE_VARIABLE(double, STORAGE_BALANCE)
        KRATOS_DEFINE_VARIABLE(double, DARCY_FLOW_BALANCE) 
        KRATOS_DEFINE_VARIABLE(double, SINKSOURCE_BALANCE)
        
        KRATOS_DEFINE_VARIABLE(double, CONCENTRATION_OLD_IT)
        
        // sinksource.cpp condition. Generalization of input (testing)
       // KRATOS_DEFINE_VARIABLE(string, SV)
       // KRATOS_DEFINE_VARIABLE(string, SV_X)
       // KRATOS_DEFINE_VARIABLE(string, SV_Y)
       // KRATOS_DEFINE_VARIABLE(string, SV_Z)


	///@} 
	///@name Type Definitions
	///@{ 

	///@} 
	///@name  Enum's
	///@{

	///@}
	///@name  Functions 
	///@{

	///@}
	///@name Kratos Classes
	///@{

	/// Short class definition.
	/** Detail class definition.
	*/
	class KratosFlowTransportApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{
		

		/// Pointer definition of KratosFlowTransportApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosFlowTransportApplication);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		KratosFlowTransportApplication();

		/// Destructor.
		virtual ~KratosFlowTransportApplication(){}


		///@}
		///@name Operators 
		///@{


		///@}
		///@name Operations
		///@{

		virtual void Register();



		///@}
		///@name Access
		///@{ 


		///@}
		///@name Inquiry
		///@{


		///@}      
		///@name Input and output
		///@{

		/// Turn back information as a string.
		virtual std::string Info() const
		{
			return "KratosFlowTransportApplication";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << Info();
			PrintData(rOStream);
		}

		///// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
      {
      	KRATOS_WATCH("in my application");
      	KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
		rOStream << "Variables:" << std::endl;
		KratosComponents<VariableData>().PrintData(rOStream);
		rOStream << std::endl;
		rOStream << "Elements:" << std::endl;
		KratosComponents<Element>().PrintData(rOStream);
		rOStream << std::endl;
		rOStream << "Conditions:" << std::endl;
		KratosComponents<Condition>().PrintData(rOStream);
      }


		///@}      
		///@name Friends
		///@{


		///@}

	protected:
		///@name Protected static Member Variables 
		///@{ 


		///@} 
		///@name Protected member Variables 
		///@{ 


		///@} 
		///@name Protected Operators
		///@{ 


		///@} 
		///@name Protected Operations
		///@{ 


		///@} 
		///@name Protected  Access 
		///@{ 


		///@}      
		///@name Protected Inquiry 
		///@{ 


		///@}    
		///@name Protected LifeCycle 
		///@{ 


		///@}

	private:

                //Initial variables of the flow-transport problem
                const Poisson2D  mPoisson2D; //and here is our element.
                const PointFlow mPointFlow; //and our condition
                const PointHeat mPointHeat; //and our condition
                const PointSource mPointSource;
                const Diffusion2D<2> mDiffusion2D;
                const Darcy mDarcy;
                const Flow< Darcy > mFlow;
                const DiffConv2D< Diffusion2D<2> > mDiffConv2D; // our template element 2D, terms convection+diffusion
                const DiffConvSt2D< DiffConv2D < Diffusion2D<2> > > mDiffConvSt2D; // our template element 2D, terms storage+convection+diffusion
                //const Flow< DiffConvSt2d < Diffusion2D<2> > > mFlow; // our template element 2D , Flow equation, terms storage+diffusion
 
                // Variables of the iteration flow-transport problem
                const FlowDarcy  mFlowDarcy; //and here is our element.
                const Transport  mTransport; //and here is our element.
               
                const SinkSource  mPointSinkSource; //and our condition
                const SinkSource  mLineSinkSource; 
                const SinkSource  mTriangleSinkSource;
                
                const SinkSourcePressure mPointSinkSourcePressure;
                const SinkSourcePressure mLineSinkSourcePressure;
                //const SinkSourcePressure mTriangleSinkSourcePressure;
               
                const MassFlow mMassFlow;
                
                const Leakage  mPointLeakage; 
                const Leakage  mLineLeakage; 
                const Leakage  mTriangleLeakage; 


		/// Assignment operator.
		KratosFlowTransportApplication& operator=(KratosFlowTransportApplication const& rOther);

		/// Copy constructor.
		KratosFlowTransportApplication(KratosFlowTransportApplication const& rOther);


		///@}    

	}; // Class KratosFlowTransportApplication 

	///@} 


	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 

	///@} 


}  // namespace Kratos.

#endif // KRATOS_FLOWTRANSPORT_APPLICATION_H_INCLUDED  defined 


