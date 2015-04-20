//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_FLOWTRANS_2D_ELEM_H_INCLUDED)
#define  KRATOS_FLOWTRANS_2D_ELEM_H_INCLUDED 

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

  class FlowTrans2D
	  : public Element
   {
   public:
     
     /// Counted pointer of Flow2D
     KRATOS_CLASS_POINTER_DEFINITION(FlowTrans2D);


    /// Default constructor.
     FlowTrans2D(IndexType NewId, GeometryType::Pointer pGeometry);
     FlowTrans2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

     /// Destructor.
     virtual ~ FlowTrans2D();


     Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

     void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
     
     void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
     
     void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

     void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);

     void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

     //void AddExplicitContribution(ProcessInfo& rCurrentProcessInfo);
      
   protected:
   
   private:
    
    unsigned int mNumber_of_nodes;
    unsigned int mStationary;
    
    double mDeltaTime;
    double mDensity;
    double mPorosity;
    double mPermeability;
    double mSpecificStorage;
    double mDiffusionCoefficient;
    
    void SetPropertiesElement();
    
    void CalculateFlow(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
    void CalculateTransport(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
    void CalculateDarcyFlow(ProcessInfo& rCurrentProcessInfo);
    
    void FlowEquationEquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo);
    void TransportEquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo);
    void GetFlowEquationDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);
    void GetTransportEquationDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);

    
    friend class Serializer;

    FlowTrans2D() : Element()
    {
    }
       
       
   }; // Class FlowTrans2D
}  // namespace Kratos.

#endif // KRATOS_FLOWTRANS_2D_ELEM_H_INCLUDED  defined
