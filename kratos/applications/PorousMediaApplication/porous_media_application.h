//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  VictorBez  $
//   Date:                $Date:  2015         $
//   Revision:            $Revision: 1.0       $
//
//

 
#if !defined(KRATOS_POROUS_MEDIA_APPLICATION_H_INCLUDED )
#define  KRATOS_POROUS_MEDIA_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream> 

// External includes:none in this case 

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

#include "custom_elements/flow_2d.h" //including the file for the element
#include "custom_elements/flowtrans_2d.h" //including the file for the element
#include "custom_elements/flowpressuretrans_2d.h" //including the file for the element 
#include "custom_conditions/sinksource.h"
#include "custom_conditions/sinksourcepressure.h"
#include "custom_conditions/leakage.h"
#include "custom_conditions/massflow.h" 

#include "includes/variables.h"
#include "includes/condition.h"         //we'll also need conditions for the point heat loads

#include "includes/ublas_interface.h"

using namespace std;

namespace Kratos
{

 
	///@name Kratos Globals

	///@{ 


	// Variables definition
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
        KRATOS_DEFINE_VARIABLE(string, SV)
        KRATOS_DEFINE_VARIABLE(string, SV_X)
        KRATOS_DEFINE_VARIABLE(string, SV_Y)
        KRATOS_DEFINE_VARIABLE(string, SV_Z)
        
	class KratosPorousMediaApplication : public KratosApplication
	{
	public:

		/// Pointer definition of KratosPorousMediaApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosPorousMediaApplication);


		/// Default constructor.
		KratosPorousMediaApplication();


		/// Destructor.
		virtual ~KratosPorousMediaApplication(){} 


		virtual void Register();


		/// Turn back information as a string.
		virtual std::string Info() const
		{
			return "KratosPorousMediaApplication";
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
      			KRATOS_WATCH("in the custom porous media application");
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



	protected:


	private:

		const Flow2D  mFlow2D; //and here is our element.
                const FlowTrans2D  mFlowTrans2D; //and here is our element.
                const FlowPressureTrans2D  mFlowPressureTrans2D; //and here is our element.
               
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
		KratosPorousMediaApplication& operator=(KratosPorousMediaApplication  const& rOther);


		/// Copy constructor.
		KratosPorousMediaApplication(KratosPorousMediaApplication const& rOther);

	}; // Class KratosPorousMediaApplication 


}  // namespace Kratos.

#endif // KRATOS_POROUS_MEDIA_APPLICATION_H_INCLUDED  defined 


