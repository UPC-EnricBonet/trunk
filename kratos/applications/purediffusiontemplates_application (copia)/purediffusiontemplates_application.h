//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_PUREDIFFUSIONTEMPLATES_APPLICATION_H_INCLUDED )
#define  KRATOS_PUREDIFFUSIONTEMPLATES_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "custom_elements/poisson_2d.h" //including the file for the element
#include "includes/condition.h"         //we'll also need conditions for the point heat loads
#include "custom_conditions/pointsource.h"         
#include "includes/ublas_interface.h"
#include "custom_elements/keytemplate2d.h"
//#include "custom_elements/superkeytemplate2d.h"
#include "custom_elements/diffusion_convection_2d.h"
#include "custom_elements/diffusion_2d.h"
#include "custom_elements/DiffConvSt2d.h"
#include "custom_elements/DiffConv2d.h"



namespace Kratos
{


	// Variables definition 
        KRATOS_DEFINE_VARIABLE(double, POINT_HEAT_SOURCE)

	class KratosPureDiffusionTemplatesApplication : public KratosApplication
	{
	public:

		/// Pointer definition of KratosPureDiffusionTemplatesApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosPureDiffusionTemplatesApplication);

		/// Default constructor.
		KratosPureDiffusionTemplatesApplication();

		/// Destructor.
		virtual ~KratosPureDiffusionTemplatesApplication(){}


		virtual void Register();

		/// Turn back information as a string.
		virtual std::string Info() const
		{
			return "KratosPureDiffusionTemplatesApplication";
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


	protected:
	
	private:
                
                const Poisson2D  mPoisson2D; //and here is our element.
		const PointSource mPointSource; //and our condition
                const Diffusion2D<2> mDiffusion2D;
                const DiffusionConvection2D<2,3> mDiffusionConvection2D;
                const KeyTemplate2D< Diffusion2D<2> > mTemplatePureDiffusion2D; // our template element 2D only diffusion
                const KeyTemplate2D< DiffusionConvection2D <2> > mTemplateDiffusionConvection2D; // our template element 2D diffusion+convection
                //const SuperKeyTemplate2D< DiffusionConvection2D <2>, Diffusion2D<2> > mIntTempDiffConv2D; // our template element 2D diffusion+convection
                //const Convection2D<2> mDiffusionConvection2D;
                const DiffConv2d< Diffusion2D<2> > mTempDiffConv2d; // our template element 2D only diffusion
                const DiffConvSt2d< DiffConv2d < Diffusion2D<2> > > mTempDiffConvSt2d; // our template element 2D diffusion+convectionn*/ mTempDiffConvSt2d
               
		/// Assignment operator.
		KratosPureDiffusionTemplatesApplication& operator=(KratosPureDiffusionTemplatesApplication const& rOther);

		/// Copy constructor.
		KratosPureDiffusionTemplatesApplication(KratosPureDiffusionTemplatesApplication const& rOther);


		///@}    

	}; // Class KratosPureDiffusionTemplatesApplication 



}  // namespace Kratos.

#endif // KRATOS_PUREDIFFUSIONTEMPLATES_APPLICATION_H_INCLUDED  defined 


