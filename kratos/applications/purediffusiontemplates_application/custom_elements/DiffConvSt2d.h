///   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_DIFFCONVST2D_ELEM_H_INCLUDED)
#define  KRATOS_DIFFUSIONCONVECTION2D_ELEM_H_INCLUDED 

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

  template<class TypeBaseClass/*, int dimension*/> class DiffConvSt2d
	  : public TypeBaseClass
   {
   public:
     
     /// Counted pointer of Poisson2D
     KRATOS_CLASS_POINTER_DEFINITION(DiffConvSt2d);

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
    DiffConvSt2d(IndexType NewId = 0):TypeBaseClass(NewId){}

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    
    DiffConvSt2d(IndexType NewId, const NodesArrayType& ThisNodes):TypeBaseClass(NewId, ThisNodes){}

    /// Constructor using a geometry object.
    DiffConvSt2d(IndexType NewId, GeometryType::Pointer pGeometry):TypeBaseClass(NewId, pGeometry){}
    

    /// Constuctor using geometry and properties.
    DiffConvSt2d(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):TypeBaseClass(NewId, pGeometry, pProperties){}


     /// Destructor.
    /*virtual*/ ~DiffConvSt2d(){};
    
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new DiffConvSt2d<TypeBaseClass/*, dimension*/>(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
                //return Element::Pointer(new BinghamFluid<TBaseElement>(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
	}
    //Calculamos el término de almacenamiento y pasamos a la clase convección-difusión
    
    void CalculateLocalSystem (MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        
                const unsigned int dimension=TypeBaseClass::getDimension(rCurrentProcessInfo);;
                const unsigned int number_of_points = TypeBaseClass::GetGeometry().size();


             
                boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX; 
                msDN_DX.resize(number_of_points,dimension,true);
                msDN_DX = ZeroMatrix(number_of_points,dimension);//initializing the matrix as zero
		boost::numeric::ublas::bounded_matrix<double,2,2> msD;
                msD = ZeroMatrix(dimension,dimension);//initializing the matrix as zero
  		array_1d<double,3> msN; //dimension = number of nodes
		array_1d<double,3> ms_temp; //dimension = number of nodes
                msN.resize(number_of_points,1);
                ms_temp.resize(number_of_points,1);

		

		if(rLeftHandSideMatrix.size1() != number_of_points)
			rLeftHandSideMatrix.resize(number_of_points,number_of_points,false);

		if(rRightHandSideVector.size() != number_of_points)
			rRightHandSideVector.resize(number_of_points,false);

		//getting data for the given geometry
		double Area;
                double Volume;
                Area=Volume;
                Volume=Area;
                GeometryType geom = TypeBaseClass::GetGeometry();

		GeometryUtils::CalculateGeometryData(geom, msDN_DX, msN, Area);

		//reading properties and conditions
		double permittivity = TypeBaseClass::GetProperties()[CONDUCTIVITY];
		for (int i=0; i<dimension; i++)
                {
                    msD(i,i)=permittivity+4.0;
                }
                KRATOS_WATCH(msD);
		// main loop	
		//const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

        noalias(rLeftHandSideMatrix) = prod(msDN_DX,Matrix(prod(msD,trans(msDN_DX))));// Bt D B
                
                for (int it = 0;  it < rLeftHandSideMatrix.size1(); it++)
                {
                      rRightHandSideVector[it];
                }
       /* KRATOS_WATCH(rLeftHandSideMatrix);
        for (int i=0; i<number_of_points;i++)
        {
            rLeftHandSideMatrix(i,i)= 1.0+Area*2.0;
        }
        KRATOS_WATCH(rLeftHandSideMatrix);*/
		rLeftHandSideMatrix *= Area;

		//subtracting the dirichlet term
		// RHS -= LHS*DUMMY_UNKNOWNs
		for(unsigned int iii = 0; iii<number_of_points; iii++)
			ms_temp[iii] = TypeBaseClass::GetGeometry()[iii].FastGetSolutionStepValue(TEMPERATURE);
		noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix,ms_temp);
                
                
		KRATOS_CATCH("");
                
               KRATOS_WATCH(rLeftHandSideMatrix);
               KRATOS_WATCH(rRightHandSideVector);
              
       
        // Go to Convection-Diffusion    
        TypeBaseClass::CalculateLocalSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
        KRATOS_WATCH("Ha pasado por Diffusion");
        KRATOS_WATCH(rLeftHandSideMatrix);
        KRATOS_WATCH(rRightHandSideVector);

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

#endif // KRATOS_KEYTEMPLATE2D_ELEM_H_INCLUDED  defined
