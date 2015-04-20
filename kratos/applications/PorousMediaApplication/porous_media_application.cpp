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
#include "geometries/point_2d.h" 
#include "geometries/line_2d.h"
#include "porous_media_application.h"
#include "includes/variables.h"

using namespace std;

namespace Kratos
{
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
        
 	KratosPorousMediaApplication::KratosPorousMediaApplication(): //constructor  do not forget to add the ":" 
		
		mFlow2D    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) ),
                mFlowTrans2D    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) ),
                mFlowPressureTrans2D    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) ),
                
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
        
	void KratosPorousMediaApplication::Register()
	{
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing KratosPorousMediaApplication... " << std::endl;

                KRATOS_REGISTER_VARIABLE(  IS_FLOW_STATIONARY )
                KRATOS_REGISTER_VARIABLE(  IS_TRANSPORT_STATIONARY )
                KRATOS_REGISTER_VARIABLE(  IS_BUOYANCY )
                                        
                KRATOS_REGISTER_VARIABLE(  HEAD_LEVEL )
                KRATOS_REGISTER_VARIABLE(  SPECIFIC_STORAGE )
		
                KRATOS_REGISTER_VARIABLE( SV )
                KRATOS_REGISTER_VARIABLE( SV_X )
                KRATOS_REGISTER_VARIABLE( SV_Y )
                KRATOS_REGISTER_VARIABLE( SV_Z )
                     
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
		KRATOS_REGISTER_ELEMENT("Flow2D", mFlow2D);  //and here is our element
                KRATOS_REGISTER_ELEMENT("FlowTrans2D", mFlowTrans2D);  //and here is our element
                KRATOS_REGISTER_ELEMENT("FlowPressureTrans2D", mFlowPressureTrans2D);  //and here is our element
                
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

