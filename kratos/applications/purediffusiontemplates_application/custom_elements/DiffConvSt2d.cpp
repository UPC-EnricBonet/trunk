//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes 
#include "includes/define.h"
//#include "custom_elements/DiffConvSt2d.h"
#include "purediffusiontemplates_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 
#include <iostream>

using namespace std;

namespace Kratos
{
	
/**/
	//************************************************************************************
	//************************************************************************************
  /*  template <> void DiffConvSt2d<TypeBaseClass,2>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        
        //const unsigned int dim = TypeBaseClass::getDimension(rCurrentProcessInfo);
        TypeBaseClass::CalculateLocalSystem (rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
    }
   /* template <> void DiffConvSt2d<TypeBaseClass,3>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        
        //const unsigned int dim = TypeBaseClass::getDimension(rCurrentProcessInfo);
        TypeBaseClass::CalculateLocalSystem (rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
    }*/


} // Namespace Kratos
	