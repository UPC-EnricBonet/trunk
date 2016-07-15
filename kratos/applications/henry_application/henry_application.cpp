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
#include "geometries/point_2d.h" 
#include "geometries/line_3d_2.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "henry_application.h"
#include "includes/variables.h"


namespace Kratos
{
	//Example
// 	KRATOS_CREATE_VARIABLE(double, AUX_MESH_VAR)
//	KRATOS_CREATE_VARIABLE(double, IS_INTERFACE);
//	KRATOS_CREATE_VARIABLE(double, NODAL_AREA);
//
        KRATOS_CREATE_VARIABLE(int, IS_FLOW_STATIONARY)
        KRATOS_CREATE_VARIABLE(int, IS_TRANSPORT_STATIONARY)
        KRATOS_CREATE_VARIABLE(int, IS_BUOYANCY)
        KRATOS_CREATE_VARIABLE(int, IS_TRANSIENT)
           
	KRATOS_CREATE_VARIABLE(double, HEAD_LEVEL) // PERMEABILITY_WATER, DENSITY & CONCENTRATION ALREADY EXIST AT KERNEL
        KRATOS_CREATE_VARIABLE(double, SPECIFIC_STORAGE)
        	

        KRATOS_CREATE_VARIABLE(double, DENSITY_ELEM)
        KRATOS_CREATE_VARIABLE(Vector, DENSITY_ELEM_TARGET)
        KRATOS_CREATE_VARIABLE(double, CONCENTRATION)
        KRATOS_CREATE_VARIABLE(double, PERMEABILITY)
        KRATOS_CREATE_VARIABLE(double, DIFFUSION_COEFFICIENT)

        KRATOS_CREATE_VARIABLE(double, PRESSURE_OLD_ITT)
        KRATOS_CREATE_VARIABLE(double, PRESCRIBED_VALUE)
        KRATOS_CREATE_VARIABLE(double, CONCENTRATION_OLD_IT)

        
        
        KRATOS_CREATE_VARIABLE(double, SINK_SOURCE)
        KRATOS_CREATE_VARIABLE(double, SINK_SOURCE_PRESS)


        
        KRATOS_CREATE_VARIABLE(double, LEAKAGE_COEFFICIENT) 
        KRATOS_CREATE_VARIABLE(double, LEVEL)
        

        

        KRATOS_CREATE_VARIABLE(double, DARCY_FLOW_X)
        KRATOS_CREATE_VARIABLE(Vector, DARCY_FLOW_X_TARGET) 
        KRATOS_CREATE_VARIABLE(double, DARCY_FLOW_Y)
        


        
        KRATOS_CREATE_VARIABLE(double, STORAGE_BALANCE)
        KRATOS_CREATE_VARIABLE(double, DARCY_FLOW_BALANCE)
        KRATOS_CREATE_VARIABLE(double, SINKSOURCE_BALANCE)
        



 	KratosHenryApplication::KratosHenryApplication()://hola

           mFlowPressureTrans2D( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3)))),
           mPointSinkSourcePressure( 0, Element::GeometryType::Pointer( new Point2D <Node<3>     >( Element::GeometryType::PointsArrayType (1)))),
           mLineSinkSourcePressure( 0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2)))),
           mMassFlow ( 0, Element::GeometryType::Pointer( new Point2D <Node<3>     >( Element::GeometryType::PointsArrayType (1))))
 	{}
 	
 	void KratosHenryApplication::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();
 		std::cout << "Initializing KratosHenryApplication... " << std::endl;
                

                KRATOS_REGISTER_VARIABLE(  IS_FLOW_STATIONARY )
                KRATOS_REGISTER_VARIABLE(  IS_TRANSPORT_STATIONARY )
                KRATOS_REGISTER_VARIABLE(  IS_BUOYANCY )
                KRATOS_REGISTER_VARIABLE( IS_TRANSIENT )  
                        
                KRATOS_REGISTER_VARIABLE(  HEAD_LEVEL )
                KRATOS_REGISTER_VARIABLE(  SPECIFIC_STORAGE )
		
 
                KRATOS_REGISTER_VARIABLE(PERMEABILITY)
                KRATOS_REGISTER_VARIABLE( DIFFUSION_COEFFICIENT )
                        

                     
                KRATOS_REGISTER_VARIABLE(  SINK_SOURCE )
                KRATOS_REGISTER_VARIABLE( SINK_SOURCE_PRESS )

    
                KRATOS_REGISTER_VARIABLE( LEAKAGE_COEFFICIENT )
                KRATOS_REGISTER_VARIABLE( LEVEL )
                

                
                KRATOS_REGISTER_VARIABLE( DENSITY_ELEM ) 
                KRATOS_REGISTER_VARIABLE( DENSITY_ELEM_TARGET )
                KRATOS_REGISTER_VARIABLE( DARCY_FLOW_X )
                KRATOS_REGISTER_VARIABLE( DARCY_FLOW_X_TARGET )
                KRATOS_REGISTER_VARIABLE( DARCY_FLOW_Y )
                KRATOS_REGISTER_VARIABLE( CONCENTRATION )
                        


                KRATOS_REGISTER_VARIABLE( STORAGE_BALANCE )
                KRATOS_REGISTER_VARIABLE( DARCY_FLOW_BALANCE )
                KRATOS_REGISTER_VARIABLE( SINKSOURCE_BALANCE ) 

 
                               
                KRATOS_REGISTER_VARIABLE (PRESSURE_OLD_ITT)
                KRATOS_REGISTER_VARIABLE (PRESCRIBED_VALUE)
                KRATOS_REGISTER_VARIABLE (CONCENTRATION_OLD_IT)
 
                        

                KRATOS_REGISTER_ELEMENT("FlowPressureTrans2D", mFlowPressureTrans2D);  //and here is our element
                

                

                
                KRATOS_REGISTER_CONDITION("PointSinkSourcePressure", mPointSinkSourcePressure);  //condition 2
                KRATOS_REGISTER_CONDITION("LineSinkSourcePressure", mLineSinkSourcePressure);  //condition 2 
                //KRATOS_REGISTER_CONDITION("TriangleSinkSourcePressure", mTriangleSinkSourcePressure);  //condition 2 
                
 
                KRATOS_REGISTER_CONDITION( "MassFlow", mMassFlow ); //condition 7



 
 	}

}  // namespace Kratos.


