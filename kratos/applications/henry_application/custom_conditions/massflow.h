#if !defined(KRATOS_MASSFLOW_CONDITION_H_INCLUDED )
#define  KRATOS_MASSFLOW_CONDITION_H_INCLUDED

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"


namespace Kratos
{
 class MassFlow : public Condition
    {
    public:
      ///@name Type Definitions
      ///@{
      
       /// Counted pointer of PointForce2D
       KRATOS_CLASS_POINTER_DEFINITION(MassFlow);
      
      /// Default constructor. 
      MassFlow(IndexType NewId, GeometryType::Pointer pGeometry);

      MassFlow(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~MassFlow();
      

      Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

      void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      
      void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

      //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);
      
      void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

      void GetDofList(DofsVectorType& rConditionalDofList,ProcessInfo& rCurrentProcessInfo);
      
   protected:
 
   private:
       
       void EquationIdVectorConcentration(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);
       
       void GetDofListConcentration(DofsVectorType& rConditionalDofList,ProcessInfo& rCurrentProcessInfo);
       
       void CalculateRightHandSideConcentration(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
       void CalculateLocalSystemConcentration(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
       
       void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, bool aCalculateStiffnessMatrixFlag, bool aCalculateResidualVectorFlag);
       
       friend class Serializer;

	// A private default constructor necessary for serialization  
	MassFlow() : Condition()
       {
       }
        

  }; // Class MassFlow 

} //namespace kratos 
#endif
