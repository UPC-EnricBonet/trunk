///   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_FLOW_ELEM_H_INCLUDED)
#define  KRATOS_FLOW_ELEM_H_INCLUDED 

// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h" 
#include "includes/serializer.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{

  template<class TypeBaseClass/*, int dimension*/> class Flow
	  : public TypeBaseClass
   {
   public:
     
     /// Counted pointer of Poisson2D
     KRATOS_CLASS_POINTER_DEFINITION(Flow);

     /// Node type (default is: Node<3>)
    typedef Node <3> NodeType;

    /// Geometry type (using with given NodeType)
    typedef Geometry<NodeType> GeometryType;

    /// Definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Properties PropertiesType;

    /// Vector type for local contributions to the linear system
    typedef Vector VectorType;

    /// Matrix type for local contributions to the linear system
    typedef Matrix MatrixType;

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
    Flow(IndexType NewId = 0):TypeBaseClass(NewId){}

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    
    Flow(IndexType NewId, const NodesArrayType& ThisNodes):TypeBaseClass(NewId, ThisNodes){}

    /// Constructor using a geometry object.
    Flow(IndexType NewId, GeometryType::Pointer pGeometry):TypeBaseClass(NewId, pGeometry){}
    

    /// Constuctor using geometry and properties.
    Flow(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):TypeBaseClass(NewId, pGeometry, pProperties){}


     /// Destructor.
    /*virtual*/ ~Flow(){};
    
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new Flow<TypeBaseClass/*, dimension*/>(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
                //return Element::Pointer(new BinghamFluid<TBaseElement>(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
	}
    //Calculamos el término de almacenamiento y pasamos a la clase convección-difusión
    
    void CalculateLocalSystem (MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {

        // Go to Convection-Diffusion    
        TypeBaseClass::CalculateLocalSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
        KRATOS_WATCH("Ha pasado por Flow");
        //KRATOS_WATCH(rLeftHandSideMatrix);
        //KRATOS_WATCH(rRightHandSideVector);

    }

   
      
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

#endif // KRATOS_FLOW_ELEM_H_INCLUDED  defined
