/*
==============================================================================
KratosALEApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: AMini $
//   Date:                $Date: Nov. 2013 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_STRUCTURAL_MESHMOVING_ELEM_2D_NONLINEAR_INCLUDED )
#define  KRATOS_STRUCTURAL_MESHMOVING_ELEM_2D_NONLINEAR_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "geometries/geometry.h"



namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{


/// This class implements a structural mesh-updating scheme in 2D using non-linear kinematics
/**
 *Implements a mesh-updating scheme in 2D treating the mesh as a structure using a linear elastic
 *material law. The kinemematics are implemented non-linearly. In Addition the solver
 *can be stabilized by an exponential law using an exponential law containing the
 *Jacobi determinant.
*/

class StructuralMeshMovingElem2DNonlin
    : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of StructuralMeshMovingElem2D
    KRATOS_CLASS_POINTER_DEFINITION(StructuralMeshMovingElem2DNonlin);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    StructuralMeshMovingElem2DNonlin(IndexType NewId, GeometryType::Pointer pGeometry);
    StructuralMeshMovingElem2DNonlin(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~StructuralMeshMovingElem2DNonlin();

    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
    //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);

    double CalculateMaximum(double , double , double );




    //void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo);

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
//      virtual String Info() const;

    /// Print information about this object.
//      virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
//      virtual void PrintData(std::ostream& rOStream) const;


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    /*		static boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX;
            static array_1d<double,3> msN; //dimension = number of nodes
            static array_1d<double,3> ms_temp_vec_np; //dimension = number of nodes*/

    ///@}
    ///@name Member Variables
    ///@{
        double mJold;
        double mJ0;
        double mxi;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    StructuralMeshMovingElem2DNonlin() {}

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    /// Assignment operator.
    //StructuralMeshMovingElem2D& operator=(const StructuralMeshMovingElem2D& rOther);

    /// Copy constructor.
    //StructuralMeshMovingElem2D(const StructuralMeshMovingElem2D& rOther);


    ///@}

}; // Class StructuralMeshMovingElem2DNonlin

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
                    StructuralMeshMovingElem2DNonlin& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
                    const StructuralMeshMovingElem2DNonlin& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_TRIANGULAR_STRUCTURAL_MESHMOVING_ELEM_2D_NONLINEAR_H_INCLUDED  defined

