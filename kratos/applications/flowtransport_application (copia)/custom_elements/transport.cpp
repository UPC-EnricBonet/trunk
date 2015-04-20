//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes 
#include "includes/define.h"
#include "custom_elements/transport.h"
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

	Element::Pointer Transport::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new Transport(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	/*Darcy::~Darcy()
	{
	}*/

	//************************************************************************************
	//************************************************************************************
	void Transport::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
	    KRATOS_TRY

	    unsigned int number_of_nodes = GetGeometry().PointsNumber();
            const int Ndim = 2;
            
            //transport properties
            unsigned int isTransportStationary = rCurrentProcessInfo[IS_TRANSPORT_STATIONARY];
            double deltaTime = rCurrentProcessInfo.GetValue(DELTA_TIME);
            
            //Get properties -> all are constant for all Elements (if they belong to the same "group of elements") 
            const double porosity = GetProperties()[POROSITY];
            const double diffusionCoefficient = GetProperties()[DIFFUSION_COEFFICIENT];
            
            // Get variable (non-constant, depend on each element- > ElementData)
            const double densityElem = GetValue(DENSITY_ELEM);
            
            //Darcy flow from flow
            array_1d <double, Ndim > darcyFlow = ZeroVector(Ndim);
            darcyFlow[0]= GetValue(DARCY_FLOW_X);
            darcyFlow[1]= GetValue(DARCY_FLOW_Y);
            
            //dimension of our local matrix (right & left)
            if (rLeftHandSideMatrix.size1() != number_of_nodes)
                rLeftHandSideMatrix.resize(number_of_nodes, number_of_nodes, false);
            if (rRightHandSideVector.size() != number_of_nodes)
                rRightHandSideVector.resize(number_of_nodes, false);

            ////Geometry data
            double area = 0.0;
            //N
            array_1d<double, 3 > N = ZeroVector(3);
            //GradN
            boost::numeric::ublas::bounded_matrix<double, 3, Ndim > DN_DX = ZeroMatrix(3,Ndim);
            GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, area);       
                
            ////LHS: Advective contribution to the stiffness matrix 
            array_1d<double, 3 > q_DN = ZeroVector(3);
            noalias(q_DN) = prod(DN_DX, darcyFlow);
            noalias(rLeftHandSideMatrix) = densityElem * outer_prod(N, q_DN);
                
            ////LHS: Contribution from storage
            // MassFactors
            boost::numeric::ublas::bounded_matrix<double, 3, 3 > massFactors = (1.0 / double(number_of_nodes)) * IdentityMatrix(3, 3);
            //inverse time
            double invDt = 1.0/ deltaTime;
            if(isTransportStationary != 1)
                 noalias(rLeftHandSideMatrix) += invDt * densityElem * porosity * massFactors;
                
            ////LHS: Contribution from diffusion
            //D
            boost::numeric::ublas::bounded_matrix<double,Ndim,Ndim> D = ZeroMatrix(Ndim,Ndim);
            D(0,0)= diffusionCoefficient * porosity;
            D(1,1)= diffusionCoefficient * porosity;
            noalias(rLeftHandSideMatrix) += densityElem * prod(DN_DX,Matrix(prod(D,trans(DN_DX))));
                
            rLeftHandSideMatrix *= area;
                
            //RHS: residual
            // RHS -= LHS*DUMMY_UNKNOWNs
            //ck+1
            array_1d<double,3> concentrationStepK_1 = ZeroVector(3);
            for(unsigned int node = 0; node< number_of_nodes; node++)
		concentrationStepK_1[node] = GetGeometry()[node].FastGetSolutionStepValue(CONCENTRATION);
            noalias(rRightHandSideVector) = (-1.0)*prod(rLeftHandSideMatrix,concentrationStepK_1);
            
            //// RHS: Contribution from storage
            if(isTransportStationary != 1)
            {
                //ck
                array_1d<double, 3 > concentrationStepK = ZeroVector(3);
                for (unsigned int node = 0; node < number_of_nodes; node++)
                    concentrationStepK[node] =  GetGeometry()[node].FastGetSolutionStepValue(CONCENTRATION, 1);
                noalias(rRightHandSideVector) += invDt * densityElem * porosity * area * prod(massFactors, concentrationStepK);
            }
            
            
            KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void Transport::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_ERROR(std::logic_error,  "method not implemented" , "");
	}	 

	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the 
	// fractional step procedure
	void Transport::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void Transport::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes,false);	

		for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(CONCENTRATION).EquationId();
	}

	//************************************************************************************
	//************************************************************************************
	void Transport::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(ElementalDofList.size() != number_of_nodes)
			ElementalDofList.resize(number_of_nodes);	

		for (unsigned int i=0;i<number_of_nodes;i++)
			ElementalDofList[i] = GetGeometry()[i].pGetDof(CONCENTRATION);

	}



} // Namespace Kratos
