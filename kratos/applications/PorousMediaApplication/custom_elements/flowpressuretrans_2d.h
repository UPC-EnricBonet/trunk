//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_FLOWPRESSURETRANS_2D_ELEM_H_INCLUDED)
#define  KRATOS_FLOWPRESSURETRANS_2D_ELEM_H_INCLUDED 

// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h" 

namespace Kratos
{

  class FlowPressureTrans2D
	  : public Element
   {
   public:
     
     /// Counted pointer of Flow2D
     KRATOS_CLASS_POINTER_DEFINITION(FlowPressureTrans2D);


    /// Default constructor.
     FlowPressureTrans2D(IndexType NewId, GeometryType::Pointer pGeometry);
     FlowPressureTrans2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

     /// Destructor.
     virtual ~ FlowPressureTrans2D();


     Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

     void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
     
     void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
     
     void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

     void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);

     //void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

     //void AddExplicitContribution(ProcessInfo& rCurrentProcessInfo);
      
   protected:
   
   private:

    //Fractional_steps IdVector
    void FlowPressureEquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo);
    void TransportEquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo);
    
    //Fractional_steps DofList
    void GetFlowPressureEquationDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);
    void GetTransportEquationDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);
    
    //Fractional_steps CalculateLocalSystem
    void CalculateFlowPressure(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
    void CalculateTransport(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
    
    //DensityElement (others commented in the strategy)
    //void CalculateDensityElement();
    //void CalculateDarcyFlowPressure(ProcessInfo& rCurrentProcessInfo);
    //void CalculateFlowBalances(ProcessInfo& rCurrentProcessInfo);
    
    friend class Serializer;

    FlowPressureTrans2D() : Element()
    {
    }
       
       
   }; // Class FlowTrans2D
}  // namespace Kratos.

#endif // KRATOS_FLOWPRESSURETRANS_2D_ELEM_H_INCLUDED  defined
