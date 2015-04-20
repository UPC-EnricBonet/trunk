#if !defined(KRATOS_Leakage_CONDITION_H_INCLUDED )
#define  KRATOS_Leakage_CONDITION_H_INCLUDED

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
 class Leakage
  : public Condition
    {
    public:
      ///@name Type Definitions
      ///@{
      
       /// Counted pointer of PointForce2D
       KRATOS_CLASS_POINTER_DEFINITION(Leakage);
      
      /// Default constructor. 
      Leakage(IndexType NewId, GeometryType::Pointer pGeometry);

      Leakage(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~Leakage();
      

      Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

      void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      
      void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

      //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);
      
      void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

      void GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo);
      

   protected:
 
 
   private:
       
       //std::string mSV;
       
       unsigned int mNumber_of_nodes; 

       double mLeakage_coeff;
       
       double mLevel;
       
        void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
                        ProcessInfo& rCurrentProcessInfo,
                        bool CalculateStiffnessMatrixFlag,
                        bool CalculateResidualVectorFlag);
       
        void CalculateLocalSystemTriangle(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
                        ProcessInfo& rCurrentProcessInfo,
                        bool CalculateStiffnessMatrixFlag,
                        bool CalculateResidualVectorFlag);
        
        void CalculateLocalSystemLine(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
                        ProcessInfo& rCurrentProcessInfo,
                        bool CalculateStiffnessMatrixFlag,
                        bool CalculateResidualVectorFlag);
        
        void CalculateLocalSystemPoint(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
                        ProcessInfo& rCurrentProcessInfo,
                        bool CalculateStiffnessMatrixFlag,
                        bool CalculateResidualVectorFlag);


	friend class Serializer;

	// A private default constructor necessary for serialization  
	Leakage() : Condition()
       {
       }
        

  }; // Class Leakage 

} //namespace kratos 
#endif
