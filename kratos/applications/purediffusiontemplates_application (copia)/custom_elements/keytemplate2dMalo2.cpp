//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes 
#include "includes/define.h"
#include "custom_elements/keytemplate2d.h"
#include "purediffusiontemplates_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 
#include <iostream>

using namespace std;

namespace Kratos
{
	

        KeyTemplate2D::KeyTemplate2D(IndexType NewId)
		: Element(NewId)
	{		
		//DO NOT ADD DOFS HERE!!!
	}
	//************************************************************************************
	//************************************************************************************
	KeyTemplate2D::KeyTemplate2D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	KeyTemplate2D::KeyTemplate2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{

	}

	




} // Namespace Kratos
