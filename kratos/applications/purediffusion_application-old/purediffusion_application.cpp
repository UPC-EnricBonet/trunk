//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.3 $
//
// 



// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d.h"
#include "purediffusion_application.h"
#include "includes/variables.h"
#include "geometries/point_2d.h"

namespace Kratos
{
         KRATOS_CREATE_VARIABLE(double, POINT_HEAT_SOURCE) // the other variables  needed in this app dont need to be created since  
                //they're already included in the kernel  ( CONDUCTIVITY and TEMPERATURE)

	//Example
// 	KRATOS_CREATE_VARIABLE(double, AUX_MESH_VAR)
//	KRATOS_CREATE_VARIABLE(double, IS_INTERFACE);
//	KRATOS_CREATE_VARIABLE(double, NODAL_AREA);
//

 	KratosPureDiffusionApplication::KratosPureDiffusionApplication()://constructor  do not forget to add the ":" 
		mPoisson2D    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) ),
		mPointSource ( 0, Element::GeometryType::Pointer( new Point2D <Node<3>     >( Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) )

	{}
 	
 	void KratosPureDiffusionApplication::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();
 		std::cout << "Initializing KratosPureDiffusionApplication... " << std::endl;
 
// 		KRATOS_REGISTER_VARIABLE( AUX_MESH_VAR )
// 		KRATOS_REGISTER_VARIABLE(IS_INTERFACE);
// 		KRATOS_REGISTER_VARIABLE(NODAL_AREA);
                KRATOS_REGISTER_VARIABLE(POINT_HEAT_SOURCE)
                // Registering elements and conditions here
		KRATOS_REGISTER_ELEMENT("Poisson2D", mPoisson2D);  //and here is our element
		KRATOS_REGISTER_CONDITION("PointSource", mPointSource) //and our condition

 
 	}

}  // namespace Kratos.


