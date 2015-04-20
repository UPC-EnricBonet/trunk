///   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_DIFFCONV2D_ELEM_H_INCLUDED)
#define  KRATOS_DIFFCONV2D_ELEM_H_INCLUDED 

// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h" 
#include "includes/serializer.h"
#include "diffusion_2d.h"

namespace Kratos
{

  template<class TypeBaseClass/*, int dimension*/> class DiffConv2d
	  : public TypeBaseClass
   {
   public:
     
     /// Counted pointer of Poisson2D
     KRATOS_CLASS_POINTER_DEFINITION(DiffConv2d);
     

     /// Node type (default is: Node<3>)
    typedef Node <3> NodeType;
    
    
    
    /*NodeType.Data();*/

    /// Geometry type (using with given NodeType)
    typedef Geometry<NodeType> GeometryType;

    /// Definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Properties PropertiesType;

    /// Vector type for local contributions to the linear system
    typedef Vector VectorType;

    /// Matrix type for local contributions to the linear system
    typedef Matrix MatrixType;
    
    //typedef  typename TypeBaseClass::dimen dimen;
    

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    typedef VectorMap<IndexType, DataValueContainer> SolutionStepsElementalDataContainerType;

        /// Type for shape function values container
    typedef Kratos::Vector ShapeFunctionsType;

    /// Type for a matrix containing the shape function gradients
    typedef Kratos::Matrix ShapeFunctionDerivativesType;

    /// Type for an array of shape function gradient matrices
    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

    /// Default constructor.
    DiffConv2d(IndexType NewId = 0):TypeBaseClass(NewId){}

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    
    DiffConv2d(IndexType NewId, const NodesArrayType& ThisNodes):TypeBaseClass(NewId, ThisNodes){}

    /// Constructor using a geometry object.
    DiffConv2d(IndexType NewId, GeometryType::Pointer pGeometry):TypeBaseClass(NewId, pGeometry){}
    

    /// Constuctor using geometry and properties.
    DiffConv2d(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):TypeBaseClass(NewId, pGeometry, pProperties){}


     /// Destructor.
    /*virtual*/ ~DiffConv2d(){};
    
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new DiffConv2d<TypeBaseClass/*, dimension*/>(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
	}
    

   
    void CalculateLocalSystem (MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {

        TypeBaseClass::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
    }


    /*void addTermConvectionDiffusion(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        TypeBaseClass::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
    }*/
    /*void addTerm(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
       {
            TypeBaseClass Troya;
            Troya.addTermDifussion(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
       }*/
      
   protected:
   
   private:
	friend class Serializer;

       virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TypeBaseClass);
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TypeBaseClass);
    }
    

       
       
   }; // Class Poisson2D
   
}  // namespace Kratos.

#endif // KRATOS_KEYTEMPLATE2D_ELEM_H_INCLUDED  defined
