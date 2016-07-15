//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: Victor Bez $
//   Date:                $Date: 2015 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes 
#include "includes/define.h"
#include "custom_elements/flowpressuretrans_2d.h"
#include "henry_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

/*

#include "utilities/enrichment_utilities.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"

*/

namespace Kratos
{


	//************************************************************************************
	//************************************************************************************
	FlowPressureTrans2D::FlowPressureTrans2D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	FlowPressureTrans2D::FlowPressureTrans2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{

	}

	Element::Pointer FlowPressureTrans2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new FlowPressureTrans2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	FlowPressureTrans2D::~FlowPressureTrans2D()
	{
	}

	//************************************************************************************
	//************************************************************************************
        void FlowPressureTrans2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
		{
		case 1:
		{
			this->FlowPressureEquationIdVector(rResult,rCurrentProcessInfo);
			break;
		}
		case 2:
		{
			this->TransportEquationIdVector(rResult,rCurrentProcessInfo);
			break;
		}
		default:
		{
			//KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
		}
		}

		KRATOS_CATCH("");
	}
        void FlowPressureTrans2D::FlowPressureEquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
                
		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes,false);	

		for (unsigned int i=0;i< number_of_nodes ;i++)
                {	
                    rResult[i] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
                }   

	}
        void FlowPressureTrans2D::TransportEquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
                unsigned int number_of_nodes = GetGeometry().PointsNumber();
		
		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes,false);	

		for (unsigned int i=0;i<number_of_nodes;i++)
		{		
                    rResult[i] = GetGeometry()[i].GetDof(CONCENTRATION).EquationId();
                }  
	}
	//************************************************************************************
	//************************************************************************************
	void FlowPressureTrans2D::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
		{
		case 1:
		{
			this->GetFlowPressureEquationDofList(rElementalDofList,rCurrentProcessInfo);
			break;
		}
		case 2:
		{
			this->GetTransportEquationDofList(rElementalDofList,rCurrentProcessInfo);
			break;
		}
		default:
		{
			//KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index at the flowpressuretrans_2d element: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
		}
		}

		KRATOS_CATCH("");
	}
        void FlowPressureTrans2D::GetFlowPressureEquationDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
                 
		if(ElementalDofList.size() != number_of_nodes)
			ElementalDofList.resize(number_of_nodes);	
                
		for (unsigned int i=0;i<number_of_nodes;i++)
			ElementalDofList[i] = GetGeometry()[i].pGetDof(PRESSURE);
	}
        void FlowPressureTrans2D::GetTransportEquationDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
                unsigned int number_of_nodes = GetGeometry().PointsNumber();	
            
		if(ElementalDofList.size() != number_of_nodes)
			ElementalDofList.resize(number_of_nodes);	
                
		for (unsigned int i=0;i<number_of_nodes;i++)
			ElementalDofList[i] = GetGeometry()[i].pGetDof(CONCENTRATION);
	}
        //************************************************************************************
	//************************************************************************************
	void FlowPressureTrans2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		//KRATOS_ERROR(std::logic_error,  "method not implemented at the flowpressuretrans_2d element" , "");
	}	 
        //************************************************************************************
	//************************************************************************************
        void FlowPressureTrans2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
		{
                    case 1:
                    {
                            this->CalculateFlowPressure(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
                            break;
                    }
                    case 2:
                    {
                            this->CalculateTransport(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
                            break;
                    }
                    default:
                    {
                            //KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index at the flowpressuretrans_2d element: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
                    }
		}

		KRATOS_CATCH("");
	}
	void FlowPressureTrans2D::CalculateFlowPressure(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
            
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
            //if(isBuoypublic:ancy != 1)
            //{
            //    array_1d <double, 2 > K_Rhog = prod (K,gravity);
            //    noalias(rRightHandSideVector) += densityElem * area * prod(DN_DX,K_Rhog);
            //}
            
            ////Contribution from storage to RHS
            if(isFlowStationary != 1 )
            {
                //Pk
                array_1d<double, 3 > pressureStepK= ZeroVector(3);
                for (unsigned int node = 0; node < number_of_nodes; node++)
                    pressureStepK[node] =  GetGeometry()[node].FastGetSolutionStepValue(PRESSURE, 1);
                
                noalias(rRightHandSideVector) += invDt * specificStorage * densityElem * area * prod(massFactors, pressureStepK);
            }
            
            
            /*
            //Internals checkings
            std::ofstream myfile;myfile.precision(10);
            myfile.open ("MatrixSystem.txt", std::ios::app); 
            myfile << " Elem "  << Id()<< std::endl;
            const double pepe1 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE, 1);
            unsigned int pepe1d = GetGeometry()[0].Id();
            const double pepe2 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE, 1);
            unsigned int pepe2d = GetGeometry()[1].Id();
            const double pepe3 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE, 1);
            unsigned int pepe3d = GetGeometry()[2].Id();
            
            const double pepe1f = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE, 0);
            const double pepe2f = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE, 0);
            const double pepe3f = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE, 0);
            
            myfile << " Node " << pepe1d << " ,  "  << pepe1 << "   ,  " << pepe1f << std::endl;
            myfile << " Node " << pepe2d << " ,  "  << pepe2 << "   ,  " << pepe2f << std::endl;
            myfile << " Node "  << pepe3d << " ,  " << pepe3 << "   ,  " << pepe3f << std::endl << std::endl;
             
            myfile.close();
            */
	}

        void FlowPressureTrans2D::CalculateUnitaryFunctions()
        {
               this->CalculateDensityElement();
               KRATOS_WATCH("CalculateUnitaryFunctions Unit test");
        }

        void FlowPressureTrans2D::CalculateTransport(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
           
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
            
            /* 
            //Internals checkings
            std::ofstream myfile;
            myfile.open ("MatrixSystem.txt", std::ios::app);
            myfile << " Elem "  << Id()<< std::endl;
            const double pepe1 = GetGeometry()[0].FastGetSolutionStepValue(CONCENTRATION, 1);
            unsigned int pepe1d = GetGeometry()[0].Id();
            const double pepe2 = GetGeometry()[1].FastGetSolutionStepValue(CONCENTRATION, 1);
            unsigned int pepe2d = GetGeometry()[1].Id();
            const double pepe3 = GetGeometry()[2].FastGetSolutionStepValue(CONCENTRATION, 1);
            unsigned int pepe3d = GetGeometry()[2].Id();
            
            const double pepe1f = GetGeometry()[0].FastGetSolutionStepValue(CONCENTRATION, 0);
            const double pepe2f = GetGeometry()[1].FastGetSolutionStepValue(CONCENTRATION, 0);
            const double pepe3f = GetGeometry()[2].FastGetSolutionStepValue(CONCENTRATION, 0);
            
            myfile << " Node " << pepe1d << " ,  " << pepe1 << "   ,  " << pepe1f << std::endl;
            myfile << " Node " << pepe2d << " ,  " << pepe2 << "   ,  " << pepe2f << std::endl;
            myfile << " Node " << pepe3d << " ,  " << pepe3 << "   ,  " << pepe3f << std::endl << std::endl;  
            myfile.close();
            */
            
        }
        
        //************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the 
	// fractional step procedure
        /*
	void FlowPressureTrans2D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY;
                        
                switch ( CurrentProcessInfo[FRACTIONAL_STEP] )
		{
                    case 1:
                    {
                            //this->CalculateDensityElement();
                            break;
                    }
                    case 2:
                    {
                            //this->CalculateDarcyFlowPressure(CurrentProcessInfo);
                            //this->CalculateFlowBalances();
                            break;
                    }
                    default:
                    {
                            KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",CurrentProcessInfo[FRACTIONAL_STEP]);
                    }       
                }
                    
		KRATOS_CATCH("");
	} 
        */
        
        void FlowPressureTrans2D::CalculateDensityElement()
        {
            KRATOS_WATCH("Estoy en CalculateDensityElement")
            unsigned int numberOfNodes = GetGeometry().PointsNumber();
            
            double densityElement = 0.0;
            const double weight = 1.0 / ((double) numberOfNodes);
            for (unsigned int node = 0; node < numberOfNodes; node++)
            {   
                const double densNode =  this->GetGeometry()[node].FastGetSolutionStepValue(DENSITY);
                densityElement += (weight *densNode);
            }
            
            this->GetValue(DENSITY_ELEM) = densityElement;
             
            
            //// This part loop over Elements and after over nodes in order to compute
            //// density and densityElem. But it has more efficiency, loop over nodes out of
            //// element class, and later in the element class compute only densityElem (previuosly code)
            

            /*unsigned int numberOfNodes = GetGeometry().PointsNumber();
            
            double densityElement = 0.0;
            const double pDensity_pConc = 700.0;
            const double referenceDensity = 1000;
            const double referenceConcentration = 0; 
            const double weight = 1.0 / ((double) numberOfNodes);

            for (unsigned int node = 0; node < numberOfNodes; node++)
            {    
               const double concNode = this->GetGeometry()[node].FastGetSolutionStepValue(CONCENTRATION);
               const double densNode =  (referenceDensity + pDensity_pConc * (concNode-referenceConcentration));
               this->GetGeometry()[node].FastGetSolutionStepValue(DENSITY) = densNode;
               densityElement += (weight *densNode);
            }
            
            this->GetValue(DENSITY_ELEM) = densityElement;*/
           
        }
            
        
        ///Previous documentation 
        
        
        ////*******************************************************************//////////////////////////////////////
        ////*******************************************************************//////////////////////////////////////
        /*
        void FlowPressureTrans2D::CalculateDarcyFlowPressure(ProcessInfo& rCurrentProcessInfo)
        {
            unsigned int numberOfNodes = this->GetGeometry().PointsNumber();
            
            unsigned int isBuoyancy = rCurrentProcessInfo[IS_BUOYANCY];
            
            //array_1d<double,2> gravity = rCurrentProcessInfo[GRAVITY];
            
            const double permeability = this->GetProperties()[PERMEABILITY_WATER];
            
            const double density = this->GetValue(b/*DENSITY_ELEM);
            
            //Area
            double area;
            const int Ndim = 2;
            //GradN
            boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX;
            //N
            array_1d<double, 3 > N;
            //Geometry staff
	    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, area);       
            
            //Compute gradient of head
            array_1d<double, Ndim > headGrad;

            for (unsigned int node = 0; node < numberOfNodes ; node++)
            {
               double currentHeadLevel =  GetGeometry()[node].FastGetSolutionStepValue(PRESSURE);
               for (unsigned int dim = 0; dim < Ndim; dim++)
                {
                   headGrad[dim]+=currentHeadLevel*DN_DX.at_element(node,dim);
                }
               
            }
 
            //if we have a  buoyancy term, substract it from the gradient
            if ( isBuoyancy != 1)
            {
                array_1d<double, Ndim > grav;
                for (unsigned int dim = 0; dim < Ndim; dim++)
                {
                    if(dim == Ndim-1)
                         grav[dim]=9.8;
                     else
                         grav[dim]=0.0;

                    headGrad[dim] += grav[dim] * (density) ;
                }
            }
            
            //K*gradh
            boost::numeric::ublas::bounded_matrix<double,2,2> K = ZeroMatrix(2,2);
            K(0,0)= permeability;
            K(1,1)= permeability;
            
            array_1d<double, 2 > darcyFlow = (-1.0)*prod(K,trans(headGrad));
            const double darcyFlowX=darcyFlow[0];
            const double darcyFlowY=darcyFlow[1];
            
            this->GetValue(IMPOSED_VELOCITY_X_VALUE/*DARCY_FLOW_X) = darcyFlowX;
            this->GetValue(IMPOSED_VELOCITY_Y_VALUE/*DARCY_FLOW_Y) = darcyFlowY;
            
            std::ofstream myfile;
            myfile.open ("MatrixSystem.txt", std::ios::app);
            myfile << "[ADarcyFlow_i_X]  "<< darcyFlowX << "[ADarcyFlow_i_Y]  "<< darcyFlowY << std::endl; 
            
           
        }
        */
        /*
        void FlowPressureTrans2D::CalculateFlowBalances(ProcessInfo& rCurrentProcessInfo)
        {
            
            unsigned int numberOfNodes = this->GetGeometry().PointsNumber();
            
            unsigned int isBuoyancy = rCurrentProcessInfo[IS_BUOYANCY];
            unsigned int isFlowStationary = rCurrentProcessInfo[IS_FLOW_STATIONARY];
            //unsigned int isTransportStationary = rCurrentProcessInfo[IS_TRANSPORT_STATIONARY];
            
            const double deltaTime = rCurrentProcessInfo[DELTA_TIME];
            //array_1d<double,2> gravity = rCurrentProcessInfo[GRAVITY];
            
            const double permeability = this->GetProperties()[PERMEABILITY_WATER];
            const double specificStorage = this->GetProperties()[SPECIFIC_STORAGE];
            
            const double density = this->GetValue(DENSITY_ELEM);
            
            //area
            double area;
            const int nDim = 2;
            
            //GradN
            boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX;
            //N
            array_1d<double, 3 > N;
            //Geometry staff
	    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, area);       
            
            //Compute gradient of head
            boost::numeric::ublas::bounded_matrix<double, nDim, 3> headGradNode;

            for (unsigned int node = 0; node < numberOfNodes; node++)
            {
               const double currentHeadLevel =  this->GetGeometry()[node].FastGetSolutionStepValue(PRESSURE);
               
               for (unsigned int dim = 0; dim < nDim; dim++)
                {
                   headGradNode.at_element(node, dim)+= currentHeadLevel * DN_DX.at_element(node,dim);
                }  
            }

            
            //if we have a  buoyancy term, substract it from the gradient
            if ( isBuoyancy != 1)
            {
                array_1d<double, nDim > grav;
                for (unsigned int dim = 0; dim < nDim; dim++)
                {
                    if(dim == nDim-1)
                         grav[dim]=-9.8;
                     else
                         grav[dim]=0.0;
                    
                    for (unsigned int node = 0; node < numberOfNodes; node++)
                    {
                        const double densNode = this->GetGeometry()[node].FastGetSolutionStepValue(DENSITY);
                        headGradNode.at_element(node, dim)+= grav[dim] * densNode;  
                    }
                }
            }
            
            //K*gradh
            boost::numeric::ublas::bounded_matrix<double,2,2> K = ZeroMatrix(2,2);
            K(0,0)= density * permeability;
            K(1,1)= density * permeability;
            
            //boost::numeric::ublas::bounded_matrix<double,2,3> darcyFlowNode;
             array_1d<double, 2 > darcyFlowNode; //= (-1.0)*prod(K,trans(headGrad));
            for (unsigned int node = 0; node < numberOfNodes ; node++)
            {
                this->GetGeometry()[node].FastGetSolutionStepValue(STORAGE_BALANCE) = 0.0;
                this->GetGeometry()[node].FastGetSolutionStepValue(DARCY_FLOW_BALANCE) = 0.0;
                this->GetGeometry()[node].FastGetSolutionStepValue(SINKSOURCE_BALANCE) = 0.0;
                
                //Storage
                if( isFlowStationary != 1 )
                {
                    const double Pk = GetGeometry()[node].FastGetSolutionStepValue(PRESSURE, 1);
                    const double Pk_1 = GetGeometry()[node].FastGetSolutionStepValue(PRESSURE);
                    const double storageBalance = ( specificStorage * (Pk_1 - Pk) ) / (deltaTime);
                    this->GetGeometry()[node].FastGetSolutionStepValue(STORAGE_BALANCE) += storageBalance;
                    this->GetGeometry()[node].FastGetSolutionStepValue(SINKSOURCE_BALANCE) += storageBalance;
                }
                
                //DarcyFlow
                darcyFlowNode[node]=-1*prod(K,trans(headGradNode[node]));
                const double darcyFlowXNode = darcyFlowNode[node][0];
                const double darcyFlowYNode = darcyFlowNode[node][1];
                const double darcyFlowBalance = darcyFlowXNode*DN_DX.at_element(node,0) +
                                        darcyFlowYNode*DN_DX.at_element(node,1);
                this->GetGeometry()[node].FastGetSolutionStepValue(DARCY_FLOW_BALANCE) += darcyFlowBalance;
                
                //Total balance
                this->GetGeometry()[node].FastGetSolutionStepValue(SINKSOURCE_BALANCE) -= darcyFlowBalance;
            }
              
        }
        */
} // Namespace Kratos
