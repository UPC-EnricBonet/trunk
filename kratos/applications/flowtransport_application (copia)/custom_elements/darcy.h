///   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_DARCY_ELEM_H_INCLUDED)
#define  KRATOS_DARCY_ELEM_H_INCLUDED 

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

  class Darcy
	  : public Element
   {
   public:
     
     /// Counted pointer of Poisson2D
     KRATOS_CLASS_POINTER_DEFINITION(Darcy);


    /// Default constructor.
     Darcy(IndexType NewId=0):Element(NewId){};
     Darcy(IndexType NewId, const NodesArrayType& ThisNodes):Element(NewId,ThisNodes){};
     Darcy(IndexType NewId, GeometryType::Pointer pGeometry):Element(NewId,pGeometry){};
     Darcy(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties):Element(NewId,pGeometry,pProperties){};
     /*Darcy(IndexType NewId, GeometryType::Pointer pGeometry);
     Darcy(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);*/

     /// Destructor.
     virtual ~ Darcy(){}


     Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

     void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
     
     void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
     
     void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

     void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);

     void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);


      
   /*protected:
   
   private:
	friend class Serializer;

       Darcy() : Element()
       {
       }*/
       
       
   }; // Class Dacry
}  // namespace Kratos.

#endif // KRATOS_DARCY_ELEM_H_INCLUDED  defined
