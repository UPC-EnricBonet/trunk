//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.3 $
//
// 


//Prueba para ver como afecta al cambiar cosillas, CAMBIO GORDO PARA VER QUE TAL!!!!!!!!!!!!!!




// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d.h"
#include "flowtransport_application.h"
#include "includes/variables.h"
#include "geometries/point_2d.h"

namespace Kratos
{
        //Variables iniciales de la aplicación
        KRATOS_CREATE_VARIABLE(double, POINT_HEAT_SOURCE) 
        KRATOS_CREATE_VARIABLE(double, POINT_FLOW_SOURCE)

        //Variables de flow-transport con iteraciones
        KRATOS_CREATE_VARIABLE(int, IS_FLOW_STATIONARY)
        KRATOS_CREATE_VARIABLE(int, IS_TRANSPORT_STATIONARY)
        KRATOS_CREATE_VARIABLE(int, IS_BUOYANCY)
           
	KRATOS_CREATE_VARIABLE(double, HEAD_LEVEL) // PERMEABILITY_WATER, DENSITY & CONCENTRATION ALREADY EXIST AT KERNEL
        KRATOS_CREATE_VARIABLE(double, SPECIFIC_STORAGE)
        KRATOS_CREATE_VARIABLE(double, SINK_SOURCE)
        KRATOS_CREATE_VARIABLE(string, SV)
        KRATOS_CREATE_VARIABLE(string, SV_X)
        KRATOS_CREATE_VARIABLE(string, SV_Y)
        KRATOS_CREATE_VARIABLE(string, SV_Z)
        
        KRATOS_CREATE_VARIABLE(double, LEAKAGE_COEFFICIENT) 
        KRATOS_CREATE_VARIABLE(double, LEVEL)
        
        KRATOS_CREATE_VARIABLE(double, INPUT_FLOW)
        KRATOS_CREATE_VARIABLE(double, PRESCRIBED_VALUE)
        
        KRATOS_CREATE_VARIABLE(double, DENSITY_ELEM)
        KRATOS_CREATE_VARIABLE(double, DARCY_FLOW_X)
        KRATOS_CREATE_VARIABLE(double, DARCY_FLOW_Y)
        KRATOS_CREATE_VARIABLE(double, DIFFUSION_COEFFICIENT)
        
        KRATOS_CREATE_VARIABLE(double, CONCENTRATION_OLD_IT)
        
        KRATOS_CREATE_VARIABLE(double, STORAGE_BALANCE)
        KRATOS_CREATE_VARIABLE(double, DARCY_FLOW_BALANCE)
        KRATOS_CREATE_VARIABLE(double, SINKSOURCE_BALANCE)
//

 	KratosFlowTransportApplication::KratosFlowTransportApplication():
                mPointFlow(0, Element::GeometryType::Pointer( new Point2D <Node<3>     >( Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) ),
                mPointHeat(0, Element::GeometryType::Pointer( new Point2D <Node<3>     >( Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) ),
                mPointSource(0, Element::GeometryType::Pointer( new Point2D <Node<3>     >( Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) ),
                mPoisson2D(0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) ),
                mDiffusion2D(0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) ),
                mDiffConv2D(0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) ),
                mDiffConvSt2D(0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>     >( Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) ),
                mFlow(0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>     >( Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) ),
                mDarcy(0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>     >( Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) ),
                
                //Variables de flow-transport iteration
                mFlowDarcy(0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>     >( Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) ),
                mTransport(0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>     >( Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) ),
                
                mPointSinkSource ( 0, Element::GeometryType::Pointer( new Point2D <Node<3>     >( Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) ),
                mLineSinkSource ( 0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
                mTriangleSinkSource ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) ),
                
                mPointSinkSourcePressure( 0, Element::GeometryType::Pointer( new Point2D <Node<3>     >( Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) ),
                mLineSinkSourcePressure( 0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
                //mTriangleSinkSourcePressure( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) ),
                
                mPointLeakage ( 0, Element::GeometryType::Pointer( new Point2D <Node<3>     >( Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) ),
                mLineLeakage ( 0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
                mTriangleLeakage ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) ),
                
                mMassFlow ( 0, Element::GeometryType::Pointer( new Point2D <Node<3>     >( Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) )
        {}
 	
 	void KratosFlowTransportApplication::Register()
 	{
                //Initial variables of the flow-transport problem

 		// calling base class register to register Kratos components
 		KratosApplication::Register();
                //REGISTER VARIABLES
 		std::cout << "Initializing KratosFlowTransportApplication... " << std::endl;
                KRATOS_REGISTER_VARIABLE(  POINT_HEAT_SOURCE )
                KRATOS_REGISTER_VARIABLE(  POINT_FLOW_SOURCE )
                //REGISTER ELEMENTS AND CONDITIONS
		KRATOS_REGISTER_ELEMENT("Poisson2D", mPoisson2D);  //and here is our element
                KRATOS_REGISTER_ELEMENT("Darcy", mDarcy);  //and here is our element
                KRATOS_REGISTER_ELEMENT("Diffusion2D", mDiffusion2D);  //and here is our element
                KRATOS_REGISTER_ELEMENT("DiffConv2D", mDiffConv2D);  //and here is our elemen 
                KRATOS_REGISTER_ELEMENT("DiffConvSt2D", mDiffConvSt2D);  //and here is our elemen  
                KRATOS_REGISTER_ELEMENT("Flow", mFlow);  //and here is our elemen  
                KRATOS_REGISTER_CONDITION( "PointFlow", mPointFlow ); //and our condition
                KRATOS_REGISTER_CONDITION( "PointHeat", mPointHeat ); //and our condition
                KRATOS_REGISTER_CONDITION( "PointSource", mPointSource ); //and our condition


                //Variables de flow-transport iteration
                KRATOS_REGISTER_VARIABLE(  IS_FLOW_STATIONARY )
                KRATOS_REGISTER_VARIABLE(  IS_TRANSPORT_STATIONARY )
                KRATOS_REGISTER_VARIABLE(  IS_BUOYANCY )
                                        
                KRATOS_REGISTER_VARIABLE(  HEAD_LEVEL )
                KRATOS_REGISTER_VARIABLE(  SPECIFIC_STORAGE )
		
                //KRATOS_REGISTER_VARIABLE( SV )
                //KRATOS_REGISTER_VARIABLE( SV_X )
                //KRATOS_REGISTER_VARIABLE( SV_Y )
                //KRATOS_REGISTER_VARIABLE( SV_Z )
                     
                KRATOS_REGISTER_VARIABLE(  SINK_SOURCE )
                        
                KRATOS_REGISTER_VARIABLE( LEAKAGE_COEFFICIENT )
                KRATOS_REGISTER_VARIABLE( LEVEL )
                
                KRATOS_REGISTER_VARIABLE(  PRESCRIBED_VALUE )  
                
                KRATOS_REGISTER_VARIABLE( DENSITY_ELEM ) 
                KRATOS_REGISTER_VARIABLE( DARCY_FLOW_X )
                KRATOS_REGISTER_VARIABLE( DARCY_FLOW_Y )
                KRATOS_REGISTER_VARIABLE( DIFFUSION_COEFFICIENT )

                KRATOS_REGISTER_VARIABLE( CONCENTRATION_OLD_IT )

                KRATOS_REGISTER_VARIABLE( STORAGE_BALANCE )
                KRATOS_REGISTER_VARIABLE( DARCY_FLOW_BALANCE )
                KRATOS_REGISTER_VARIABLE( SINKSOURCE_BALANCE ) 
                        
		// Registering elements and conditions here
                KRATOS_REGISTER_ELEMENT("FlowDarcy", mFlowDarcy);  //and here is our element
                KRATOS_REGISTER_ELEMENT("Transport", mTransport);  //and here is our element
                
                KRATOS_REGISTER_CONDITION( "PointSinkSource", mPointSinkSource ); //condition 1
                KRATOS_REGISTER_CONDITION("LineSinkSource", mLineSinkSource);  //condition 2
                KRATOS_REGISTER_CONDITION( "TriangleSinkSource", mTriangleSinkSource ); //condition 3
                
                KRATOS_REGISTER_CONDITION("PointSinkSourcePressure", mPointSinkSourcePressure);  //condition 2
                KRATOS_REGISTER_CONDITION("LineSinkSourcePressure", mLineSinkSourcePressure);  //condition 2 
                //KRATOS_REGISTER_CONDITION("TriangleSinkSourcePressure", mTriangleSinkSourcePressure);  //condition 2 
                
                KRATOS_REGISTER_CONDITION( "PointLeakage", mPointLeakage ); //condition 4
                KRATOS_REGISTER_CONDITION("LineLeakage", mLineLeakage);  //condition 5
                KRATOS_REGISTER_CONDITION( "TriangleLeakage", mTriangleLeakage ); //condition 6
                
                KRATOS_REGISTER_CONDITION( "MassFlow", mMassFlow ); //condition 7

 
 	}

}  // namespace Kratos.


