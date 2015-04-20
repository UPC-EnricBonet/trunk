///   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_DIFFUSION_2D_ELEM_H_INCLUDED)
#define  KRATOS_DIFFUSION_2D_ELEM_H_INCLUDED 

// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h" 
#include "includes/serializer.h"

namespace Kratos
{

  template <unsigned int TDim, unsigned int TNumNodes = TDim + 1> class Diffusion2D: public Element
   {

   public:
     
     /// Counted pointer of Poisson2D
     KRATOS_CLASS_POINTER_DEFINITION(Diffusion2D);

    ///base type: an IndexedObject that automatically has a unique number
    typedef IndexedObject BaseType;

    ///definition of node type (default is: Node<3>)
    typedef Node < 3 > NodeType;

    /**
     * Properties are used to store any parameters
     * related to the constitutive law
     */
    typedef Properties PropertiesType;

    ///definition of the geometry type with given NodeType
    typedef Geometry<NodeType> GeometryType;

    ///definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;
    

   // typedef unsigned int TDim dimen; 

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;
    
    
   /// Type for shape function values container
    typedef Kratos::Vector ShapeFunctionsType;

    /// Type for a matrix containing the shape function gradients
    typedef Kratos::Matrix ShapeFunctionDerivativesType;

    /// Type for an array of shape function gradient matrices
    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;
    /// Default constructor.

     Diffusion2D(IndexType NewId=0):Element(NewId){};
     Diffusion2D(IndexType NewId, const NodesArrayType& ThisNodes):Element(NewId,ThisNodes){};
     Diffusion2D(IndexType NewId, GeometryType::Pointer pGeometry):Element(NewId,pGeometry){};
     Diffusion2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties):Element(NewId,pGeometry,pProperties){};

     /// Destructor.
     /*virtual*/ ~Diffusion2D(){};


     Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new Diffusion2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

     virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
     
     virtual void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
     
     virtual void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

     virtual void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);

     virtual void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);
     
     const unsigned int getDimension(ProcessInfo& rCurrentProcessInfo);

      
   protected:
   
   private:


    friend class Serializer;

   /* virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }*/
       
       
   }; // Class Poisson2D
}  // namespace Kratos.

#endif // KRATOS_POISSON_2D_ELEM_H_INCLUDED  defined
