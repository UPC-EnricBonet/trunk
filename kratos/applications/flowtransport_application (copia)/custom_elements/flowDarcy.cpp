//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes 
#include "includes/define.h"
#include "custom_elements/flowDarcy.h"
#include "flowtransport_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 
#include <iostream>

using namespace std;

namespace Kratos
{
	


	//************************************************************************************
	//************************************************************************************
	/*Darcy::Darcy(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	Darcy::Darcy(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{

	}*/

	Element::Pointer FlowDarcy::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new FlowDarcy(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	/*Darcy::~Darcy()
	{
	}*/

	//************************************************************************************
	//************************************************************************************
	void FlowDarcy::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		unsigned int number_of_nodes = GetGeometry().PointsNumber();
            const int NDim = 2;
            
            //Flow properties
            unsigned int isFlowStationary = rCurrentProcessInfo[IS_FLOW_STATIONARY];        
            unsigned int isBuoyancy = rCurrentProcessInfo[IS_BUOYANCY];
            
            double deltaTime = rCurrentProcessInfo.GetValue(DELTA_TIME);
            array_1d<double,2> gravity = rCurrentProcessInfo[GRAVITY];
            
            //Get properties -> all are constant for all Elements (if they belong to the same "group of elements")  
            double specificStorage = 0.0;
            if(isFlowStationary!=1)
                specificStorage = GetProperties()[SPECIFIC_STORAGE];
                //const double porosity = GetProperties()[POROSITY];
            const double permeability = GetProperties()[PERMEABILITY_WATER];
            
            // Get variable (non-constant, depend on each element- > ElementData)
            const double densityElem = GetValue(DENSITY_ELEM);
             
            //dimension of our local matrix (right & left)
            if (rLeftHandSideMatrix.size1() != number_of_nodes)
                rLeftHandSideMatrix.resize(number_of_nodes,number_of_nodes, false);
            if (rRightHandSideVector.size() != number_of_nodes)
                rRightHandSideVector.resize(number_of_nodes, false);
                
            ////Element Geometry (N,gradN & area)
            double area = 0.0;
            //GradN
            boost::numeric::ublas::bounded_matrix<double, 3, NDim > DN_DX = ZeroMatrix(3,NDim);
            //N
            array_1d<double, 3 > N = ZeroVector(3);
		
            GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, area);    

            ////LHS: Contribution from darcyFlow
            //K
            boost::numeric::ublas::bounded_matrix<double, NDim, NDim> K = ZeroMatrix(NDim,NDim); 
            K(0,0)= densityElem * permeability;
            K(1,1)= densityElem * permeability;
                
	    noalias(rLeftHandSideMatrix) = prod(DN_DX,Matrix(prod(K,trans(DN_DX))));
              
            ////LHS: Contribution from storage
            //MassFactor
            boost::numeric::ublas::bounded_matrix<double, 3, 3 > massFactors = (1.0/double(number_of_nodes))*IdentityMatrix(3, 3);
            //Inverse time
            double invDt = 1.0/ deltaTime;
            if(isFlowStationary != 1 )
                    noalias(rLeftHandSideMatrix) += invDt * specificStorage * densityElem * massFactors;
            
            rLeftHandSideMatrix *= area;
 
            //Residual
            // RHS -= LHS*DUMMY_UNKNOWNs
            //Pk+1
            array_1d<double,3> pressureStepK_1 = ZeroVector(3);
            for(unsigned int node = 0; node< number_of_nodes; node++)
		pressureStepK_1[node] = GetGeometry()[node].FastGetSolutionStepValue(PRESSURE);
            
            noalias(rRightHandSideVector) = (-1.0)*prod(rLeftHandSideMatrix,pressureStepK_1);
             
            
            //RHS: Contribution from bouyancy 
            if(isBuoyancy != 1)
            {
                array_1d <double, 2 > K_Rhog = prod (K,gravity);
                noalias(rRightHandSideVector) += densityElem * area * prod(DN_DX,K_Rhog);
            }
            
            ////Contribution from storage to RHS
            if(isFlowStationary != 1 )
            {
                //Pk
                array_1d<double, 3 > pressureStepK= ZeroVector(3);
                for (unsigned int node = 0; node < number_of_nodes; node++)
                    pressureStepK[node] =  GetGeometry()[node].FastGetSolutionStepValue(PRESSURE, 1);
                
                noalias(rRightHandSideVector) += invDt * specificStorage * densityElem * area * prod(massFactors, pressureStepK);
            }
            
	    KRATOS_CATCH("");
	}
                
	//************************************************************************************
	//************************************************************************************
	void FlowDarcy::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_ERROR(std::logic_error,  "method not implemented" , "");
	}	 

	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the 
	// fractional step procedure
	void FlowDarcy::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void FlowDarcy::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes,false);	

		for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
	}

	//************************************************************************************
	//************************************************************************************
	void FlowDarcy::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(ElementalDofList.size() != number_of_nodes)
			ElementalDofList.resize(number_of_nodes);	

		for (unsigned int i=0;i<number_of_nodes;i++)
			ElementalDofList[i] = GetGeometry()[i].pGetDof(PRESSURE);

	}



} // Namespace Kratos
