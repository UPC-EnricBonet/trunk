//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes 
#include "includes/define.h"
#include "custom_elements/darcy.h"
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

	Element::Pointer Darcy::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new Darcy(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	/*Darcy::~Darcy()
	{
	}*/

	//************************************************************************************
	//************************************************************************************
	void Darcy::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX;
		boost::numeric::ublas::bounded_matrix<double,2,2> msD = ZeroMatrix(2,2); //initializing the matrix as zero
  		array_1d<double,3> msN; //dimension = number of nodes
		array_1d<double,3> ms_temp; //dimension = number of nodes
                const unsigned int number_of_points = GetGeometry().size();

                // S Compresibility matrix
                boost::numeric::ublas::bounded_matrix<double, 3, 3 > S_matrix= ZeroMatrix(3,3);
                const int numberOfDofs=1;
                double delta_t = rCurrentProcessInfo[DELTA_TIME];
                const double porosity =0.05;
                const double one_third = 1.0/double(number_of_points); //in 3d we would need one_quarter
                unsigned int index=0;
		if(rLeftHandSideMatrix.size1() != number_of_points)
			rLeftHandSideMatrix.resize(number_of_points,number_of_points,false);
		rLeftHandSideMatrix=ZeroMatrix(number_of_points,number_of_points);

		if(rRightHandSideVector.size() != number_of_points)
			rRightHandSideVector.resize(number_of_points,false);
		rRightHandSideVector=ZeroVector(number_of_points);

		//getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);

		//reading properties and conditions
		double permittivity = GetProperties()[CONDUCTIVITY];
		msD(0,0)=permittivity;
		msD(1,1)=permittivity;
		
		// main loop	
		//const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

        noalias(rLeftHandSideMatrix) = prod(msDN_DX,Matrix(prod(msD,trans(msDN_DX))));// Bt D B
                
                for (int it = 0;  it < rLeftHandSideMatrix.size1(); it++)
                {
                      rRightHandSideVector[it];
                }

		rLeftHandSideMatrix *= Area;
 
                 //Calculate S contribution & ensemble
		index=0;
		for (unsigned int node=0; node<number_of_points ; node++)            
		{
                    S_matrix(index,index) = one_third*Area*porosity / delta_t; 
                    index+=numberOfDofs;          
		}
                //KRATOS_WATCH(S_matrix)
                //ADDING THE MASS MATRIX TO THE LHS        
                noalias(rLeftHandSideMatrix)+=S_matrix;

		for(unsigned int iii = 0; iii<number_of_points; iii++)
			ms_temp[iii] = GetGeometry()[iii].GetSolutionStepValue(PRESSURE);
		//ADDING THE MASS MATRIX MULTIPLIED BY THE LAST STEP TEMPERATURE TO THE RHS
		noalias(rRightHandSideVector) += prod(S_matrix,ms_temp);

		
		//subtracting the dirichlet term
		// RHS -= LHS*DUMMY_UNKNOWNs
		for(unsigned int iii = 0; iii<number_of_points; iii++)
			ms_temp[iii] = GetGeometry()[iii].FastGetSolutionStepValue(PRESSURE);
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp);
                KRATOS_WATCH("Ha pasado por Flux");
                KRATOS_WATCH(delta_t)
		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void Darcy::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_ERROR(std::logic_error,  "method not implemented" , "");
	}	 

	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the 
	// fractional step procedure
	void Darcy::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void Darcy::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes,false);	

		for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
	}

	//************************************************************************************
	//************************************************************************************
	void Darcy::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(ElementalDofList.size() != number_of_nodes)
			ElementalDofList.resize(number_of_nodes);	

		for (unsigned int i=0;i<number_of_nodes;i++)
			ElementalDofList[i] = GetGeometry()[i].pGetDof(PRESSURE);

	}



} // Namespace Kratos
