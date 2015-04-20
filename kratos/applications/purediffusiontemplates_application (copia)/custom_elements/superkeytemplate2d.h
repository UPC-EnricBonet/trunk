///   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_SUPERKEYTEMPLATE2D_ELEM_H_INCLUDED)
#define  KRATOS_SUPERKEYTEMPLATE2D_ELEM_H_INCLUDED 

// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
//#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h" 
#include "includes/serializer.h"
#include "custom_elements/diffusion_2d.h"
#include "custom_elements/diffusion_convection_2d.h"
namespace Kratos
{

  template<class TypeBaseClass1, class TypeBaseClass2 > class SuperKeyTemplate2D: public  TypeBaseClass1,  public TypeBaseClass2
   {
      public:
     
     /// Counted pointer of Poisson2D
     KRATOS_CLASS_POINTER_DEFINITION(SuperKeyTemplate2D);

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


    /// Default constructor.
    SuperKeyTemplate2D(IndexType NewId = 0):TypeBaseClass1(NewId),TypeBaseClass2(NewId){}

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    
    SuperKeyTemplate2D(IndexType NewId, const NodesArrayType& ThisNodes):TypeBaseClass1(NewId,ThisNodes),TypeBaseClass2(NewId,ThisNodes){}

    /// Constructor using a geometry object.
    SuperKeyTemplate2D(IndexType NewId, GeometryType::Pointer pGeometry):TypeBaseClass1(NewId,pGeometry),TypeBaseClass2(NewId,pGeometry){}
    

    /// Constuctor using geometry and properties.
    SuperKeyTemplate2D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):TypeBaseClass1(NewId,pGeometry,pProperties),TypeBaseClass2(NewId,pGeometry,pProperties){}


     /// Destructor.
    /*virtual*/ ~SuperKeyTemplate2D(){};
    
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new SuperKeyTemplate2D<TypeBaseClass1,TypeBaseClass2>(NewId, this->getGeometryData(ThisNodes), pProperties));
                //return Element::Pointer(new BinghamFluid<TBaseElement>(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
	}

    void CalculateLocalSystem (MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_WATCH("hola ConvectionDiffusionElement");
       /* TypeBaseClass1 Troya;
        Troya.CalculateLocalSystem (rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
        delete Troya;
        TypeBaseClass2 Troya2;
        Troya2.CalculateLocalSystem (rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
        delete Troya2;*/
   
    }

    GeometryType const& getGeometryData( NodesArrayType const& ThisNodes) const 
    {
        return this->getGeometryDataBase(ThisNodes);
    }
    
    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_WATCH("hola ConvectionDiffusionElement");
    }
     
     void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
     {
        KRATOS_WATCH("hola ConvectionDiffusionElement");
     }

     void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
     {
        KRATOS_WATCH("hola ConvectionDiffusionElement");
    }

     void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
     {
        KRATOS_WATCH("hola ConvectionDiffusionElement");
    }


      
 
   
   private:
	friend class Serializer;

    /*   virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TypeBaseClass);
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TypeBaseClass);
    }
       */
       
   }; // Class Poisson2D
}  // namespace Kratos.

#endif // KRATOS_KEYTEMPLATE2D_ELEM_H_INCLUDED  defined
