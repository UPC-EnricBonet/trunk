#if !defined(KRATOS_SinkSource_CONDITION_H_INCLUDED )
#define  KRATOS_SinkSource_CONDITION_H_INCLUDED

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
 class SinkSource : public Condition
    {
    public:
      ///@name Type Definitions
      ///@{
      
       /// Counted pointer of PointForce2D
       KRATOS_CLASS_POINTER_DEFINITION(SinkSource);
      
      /// Default constructor. 
      SinkSource(IndexType NewId, GeometryType::Pointer pGeometry);

      SinkSource(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~SinkSource();
      

      Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

      void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      
      void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

      //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);
      
      void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

      void GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo);
      
   protected:
 
   private:

       unsigned int mNumber_of_nodes; 
       
       std::vector<std::string> mSVs;
       
       bool mSinkSourceIsVectorial;
       
       double mSinkSource;

        void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, bool aCalculateStiffnessMatrixFlag);
        void CalculateLocalSystemTriangle(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
        void CalculateLocalSystemLine(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
        void CalculateLocalSystemPoint(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
        void CalculateNormalFace2D(array_1d<double,3>& aN);
        void CalculateNormalFace3D(array_1d<double,3>& aN);
       
	friend class Serializer;

	// A private default constructor necessary for serialization  
	SinkSource() : Condition()
       {
       }
        

  }; // Class SinkSource 

} //namespace kratos 
#endif
