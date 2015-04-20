///   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_KEYTEMPLATE2D_H_INCLUDED)
#define  KRATOS_KEYTEMPLATE2D_H_INCLUDED 

// System includes 





// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"

namespace Kratos
{

  template< class TypeBaseClass,
            class TShapeFunctionValues = typename TypeBaseClass::ShapeFunctionsType,
            class TShapeFunctionGradients = typename TypeBaseClass::ShapeFunctionDerivativesType> 
   class KeyTemplate2D : public TypeBaseClass
   {
   public:
     
     /// Counted pointer of Poisson2D
     KRATOS_CLASS_POINTER_DEFINITION(KeyTemplate2D);

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
    KeyTemplate2D(IndexType NewId = 0) : TypeBaseClass(NewId){}

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    
    KeyTemplate2D(IndexType NewId, const NodesArrayType& ThisNodes) : TypeBaseClass(NewId, ThisNodes){}

    /// Constructor using a geometry object.
    KeyTemplate2D(IndexType NewId, GeometryType::Pointer pGeometry) : TypeBaseClass(NewId, pGeometry){}
    

    /// Constuctor using geometry and properties.
    KeyTemplate2D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : TypeBaseClass(NewId, pGeometry, pProperties){}


     /// Destructor.
    virtual ~KeyTemplate2D(){}
    
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new KeyTemplate2D<TypeBaseClass>(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
                //return Element::Pointer(new BinghamFluid<TBaseElement>(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
	}

    /*void CalculateLocalSystem (MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        TypeBaseClass Troya;
        Troya.CalculateLocalSystem (rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
    }*/
      
   //protected:

    virtual int Check(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        int Error = 0;

        // Check that any required model parameters are defined


        // Call the underlying element's check routine
        Error = TypeBaseClass::Check(rCurrentProcessInfo);

        return Error;

        KRATOS_CATCH("");
    }


    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "KeyTemplate2D " ;
        buffer << TypeBaseClass::Info();
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "KeyTemplate2D ";
        TypeBaseClass::PrintInfo(rOStream);
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}
   
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

#endif // KRATOS_KEYTEMPLATE2D_H_INCLUDED  defined
