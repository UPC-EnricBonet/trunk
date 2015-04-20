#if !defined(KRATOS_SinkSourcePressure_CONDITION_H_INCLUDED )
#define  KRATOS_SinkSourcePressure_CONDITION_H_INCLUDED

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
 class SinkSourcePressure : public Condition
    {
    public:
      ///@name Type Definitions
      ///@{
      
       /// Counted pointer of PointForce2D
       KRATOS_CLASS_POINTER_DEFINITION(SinkSourcePressure);
      
      /// Default constructor. 
      SinkSourcePressure(IndexType NewId, GeometryType::Pointer pGeometry);

      SinkSourcePressure(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~SinkSourcePressure();
      

      Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

      void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      
      void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

      //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);
      
      void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

      void GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& rCurrentProcessInfo);
      
   protected:
 
   private:
       
       void GetDofListPressure(DofsVectorType& rConditionalDofList,ProcessInfo& rCurrentProcessInfo);

       void EquationIdVectorPressure(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);
       
       void CalculateRightHandSidePressure(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
       void CalculateLocalSystemPressure(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
       
       void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, bool aCalculateStiffnessMatrixFlag);
       void CalculateLocalSystemTriangle(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double aSinkSource);
       void CalculateLocalSystemLine(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double aSinkSource);
       void CalculateLocalSystemPoint(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double aSinkSource);
       
       friend class Serializer;

	// A private default constructor necessary for serialization  
        SinkSourcePressure() : Condition()
       {
       }
        

  }; // Class SinkSourcePressure 

} //namespace kratos 
#endif
