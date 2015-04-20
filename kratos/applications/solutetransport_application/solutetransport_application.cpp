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
#include "solutetransport_application.h"
#include "includes/variables.h"
#include "geometries/point_2d.h"


namespace Kratos
{
	//Example
// 	KRATOS_CREATE_VARIABLE(double, AUX_MESH_VAR)
//	KRATOS_CREATE_VARIABLE(double, IS_INTERFACE);
//	KRATOS_CREATE_VARIABLE(double, NODAL_AREA);
//
        KRATOS_CREATE_VARIABLE(double, SINK_SOURCE) 
 	KratosSoluteTransportApplication::KratosSoluteTransportApplication(): //constructor  do not forget to add the ":" 
		mAdvDiffEq2D    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) ),
		mSinkSource ( 0, Element::GeometryType::Pointer( new Point2D <Node<3>     >( Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) )
 	{}
 	
 	void KratosSoluteTransportApplication::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();
 		std::cout << "Initializing KratosSoluteTransportApplication... " << std::endl;

                KRATOS_REGISTER_VARIABLE(SINK_SOURCE)

		// Registering elements and conditions here
		KRATOS_REGISTER_ELEMENT("AdvDiffEq2D", mAdvDiffEq2D);  //and here is our element
		KRATOS_REGISTER_CONDITION("SinkSource", mSinkSource);  //and our condition
 
// 		KRATOS_REGISTER_VARIABLE( AUX_MESH_VAR )
// 		KRATOS_REGISTER_VARIABLE(IS_INTERFACE);
// 		KRATOS_REGISTER_VARIABLE(NODAL_AREA);

 
 	}

}  // namespace Kratos.


